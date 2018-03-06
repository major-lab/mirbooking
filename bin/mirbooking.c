#include <errno.h>
#include <fcntl.h>
#include <gio/gio.h>
#include <glib.h>
#include <glib/gprintf.h>
#include <glib/gstdio.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mirbooking.h>

G_DEFINE_AUTOPTR_CLEANUP_FUNC (FILE, fclose)

static gchar   *mirnas_file      = NULL;
static gchar   *targets_file     = NULL;
static gchar   *cds_regions_file = NULL;
static gchar   *score_table_file = NULL;
static gchar   *quantities_file  = NULL;
static gchar   *output_file      = NULL;
static gdouble  threshold        = MIRBOOKING_BROKER_DEFAULT_THRESHOLD;
static gdouble  log_base         = MIRBOOKING_BROKER_DEFAULT_LOG_BASE;
static gsize    seed_offset      = MIRBOOKING_PRECOMPUTED_SCORE_TABLE_DEFAULT_SEED_OFFSET;
static gsize    seed_length      = MIRBOOKING_PRECOMPUTED_SCORE_TABLE_DEFAULT_SEED_LENGTH;
static gsize    prime5_footprint = MIRBOOKING_BROKER_DEFAULT_5PRIME_FOOTPRINT;
static gsize    prime3_footprint = MIRBOOKING_BROKER_DEFAULT_3PRIME_FOOTPRINT;

static GOptionEntry MIRBOOKING_OPTION_ENTRIES[] =
{
    {"mirnas",                0, 0, G_OPTION_ARG_FILENAME, &mirnas_file,      "MiRNAs sequences FASTA file",                                                                     "FILE"},
    {"targets",               0, 0, G_OPTION_ARG_FILENAME, &targets_file,     "Transcripts sequences FASTA file",                                                                "FILE"},
    {"cds-regions",           0, 0, G_OPTION_ARG_FILENAME, &cds_regions_file, "Coding regions as a two-column (accession, 1-based inclusive interval) TSV file",                 "FILE"},
    {"score-table",           0, 0, G_OPTION_ARG_FILENAME, &score_table_file, "Precomputed seed-MRE duplex table as a row-major big-endian float matrix file",                   "FILE"},
    {"quantities",            0, 0, G_OPTION_ARG_FILENAME, &quantities_file,  "MiRNA and targets quantities as a two-column (accession, quantity) TSV file (defaults to stdin)", "FILE"},
    {"output",                0, 0, G_OPTION_ARG_FILENAME, &output_file,      "Output destination file (defaults to stdout)",                                                    "FILE"},
    {"threshold",             0, 0, G_OPTION_ARG_DOUBLE,   &threshold,        "Probability threshold for site matching",                                                         G_STRINGIFY (MIRBOOKING_BROKER_DEFAULT_THRESHOLD)},
    {"log-base",              0, 0, G_OPTION_ARG_DOUBLE,   &log_base,         "Logarithm base for spreading quantites across sites",                                             G_STRINGIFY (MIRBOOKING_BROKER_DEFAULT_LOG_BASE)},
    {"seed-offset",           0, 0, G_OPTION_ARG_INT,      &seed_offset,      "MiRNA seed offset",                                                                               G_STRINGIFY (MIRBOOKING_PRECOMPUTED_SCORE_TABLE_DEFAULT_SEED_OFFSET)},
    {"seed-length",           0, 0, G_OPTION_ARG_INT,      &seed_length,      "MiRNA seed length",                                                                               G_STRINGIFY (MIRBOOKING_PRECOMPUTED_SCORE_TABLE_DEFAULT_SEED_LENGTH)},
    {"5prime-footprint",      0, 0, G_OPTION_ARG_INT,      &prime5_footprint, "Footprint in the MRE's 5' direction",                                                             G_STRINGIFY (MIRBOOKING_BROKER_DEFAULT_5PRIME_FOOTPRINT)},
    {"3prime-footprint",      0, 0, G_OPTION_ARG_INT,      &prime3_footprint, "Footprint in the MRE's 3' direction",                                                             G_STRINGIFY (MIRBOOKING_BROKER_DEFAULT_3PRIME_FOOTPRINT)},
    {NULL}
};

static gpointer
gpointer_from_two_guint16 (guint16 a, guint16 b)
{
    union {
        gpointer p;
        guint16 u[2];
    } ret;
    ret.u[0] = a;
    ret.u[1] = b;
    return ret.p;
}

static void
two_guint16_from_gpointer (gpointer ptr, guint16 out[2])
{
    union {
        gpointer p;
        guint16 u[2];
    } ret;
    ret.p = ptr;
    out[0] = ret.u[0];
    out[1] = ret.u[1];
}

static void
read_sequences_from_fasta (FILE        *file,
                           GMappedFile *mapped_file,
                           gboolean     accession_in_comment,
                           GHashTable  *sequences_hash)
{
    gchar *accession;
    gchar *name;
    gchar *seq;
    gchar line[1024];

    while (fgets (line, sizeof (line), file))
    {
        if (line[0] == '>')
        {
            accession = strtok (line + 1, " |");

            if (accession_in_comment)
            {
                name      = accession;
                accession = strtok (NULL, " ");
            }
            else if (strtok (NULL, "(") != NULL)
            {
                name = strtok (NULL, ")");
            }
            else
            {
                name = NULL;
            }

            seq = g_mapped_file_get_contents (mapped_file) + ftell (file);

            gsize remaining = g_mapped_file_get_length (mapped_file) - ftell (file);
            gchar *next_seq = memchr (seq, '>', remaining);
            gsize seq_len = next_seq == NULL ? remaining - 1 : next_seq - seq - 1;

            MirbookingSequence *sequence;

            if (g_str_has_prefix (accession, "MIMAT") || g_str_has_prefix (accession, "SYNTH"))
            {
                sequence = MIRBOOKING_SEQUENCE (mirbooking_mirna_new_with_name (accession, name));
            }
            else
            {
                sequence = MIRBOOKING_SEQUENCE (mirbooking_target_new_with_name (accession, name));
            }

            mirbooking_sequence_set_raw_sequence (sequence,
                                                  seq,
                                                  seq_len);

            g_hash_table_insert (sequences_hash,
                                 g_strdup (accession),
                                 sequence);
        }
    }
}

typedef enum _MirbookingRegion
{
    MIRBOOKING_REGION_5PRIME_UTR,
    MIRBOOKING_REGION_CDS,
    MIRBOOKING_REGION_3PRIME_UTR,
    MIRBOOKING_REGION_UNKNOWN
} MirbookingRegion;

MirbookingRegion
mirbooking_region_from_target_site (const MirbookingTargetSite *target_site,
                                    gsize                       a,
                                    gsize                       b)
{
    if (target_site->position < a)
    {
        return MIRBOOKING_REGION_5PRIME_UTR;
    }
    else if (target_site->position < b)
    {
        return MIRBOOKING_REGION_CDS;
    }
    else
    {
        return MIRBOOKING_REGION_3PRIME_UTR;
    }
}

static gchar *
mirbooking_region_to_string (MirbookingRegion region)
{
    switch (region)
    {
        case MIRBOOKING_REGION_5PRIME_UTR:
            return "5'UTR";
        case MIRBOOKING_REGION_CDS:
            return "CDS";
        case MIRBOOKING_REGION_3PRIME_UTR:
            return "3'UTR";
        case MIRBOOKING_REGION_UNKNOWN:
            return "N/A";
        default:
            g_assert_not_reached ();
    }
}

int
main (gint argc, gchar **argv)
{
    g_autoptr (GOptionContext) context = g_option_context_new ("[FILE]");

    g_option_context_add_main_entries (context,
                                       MIRBOOKING_OPTION_ENTRIES,
                                       NULL);

    GError *error = NULL;
    if (!g_option_context_parse (context,
                                 &argc,
                                 &argv,
                                 &error))
    {
        g_printerr ("%s (%s, %u).\n", error->message, g_quark_to_string (error->domain), error->code);
        return EXIT_FAILURE;
    }

    g_autoptr (MirbookingBroker) mirbooking = mirbooking_broker_new ();

    mirbooking_broker_set_threshold (mirbooking, threshold);
    mirbooking_broker_set_log_base (mirbooking, log_base);
    mirbooking_broker_set_5prime_footprint (mirbooking, prime5_footprint);
    mirbooking_broker_set_3prime_footprint (mirbooking, prime3_footprint);

    if (score_table_file == NULL)
    {
        g_printerr ("The '--score-table' argument is required.\n");
        return EXIT_FAILURE;
    }

    g_autoptr (GMappedFile) score_map = g_mapped_file_new (score_table_file, FALSE, &error);
    if (score_map == NULL)
    {
        g_printerr ("%s (%s, %u).\n", error->message, g_quark_to_string (error->domain), error->code);
        return EXIT_FAILURE;
    }

    if (g_mapped_file_get_length (score_map) != (1 << 2 * seed_length) * (1 << 2 * seed_length) * sizeof (gfloat))
    {
        g_printerr ("The specified '--seed-length' parameter would require a %luB score table, but %luB were provided.\n",
                    (1l << 2 * seed_length) * (1l << 2 * seed_length) * sizeof (gfloat),
                    g_mapped_file_get_length (score_map));
        return EXIT_FAILURE;
    }

    g_autoptr (MirbookingPrecomputedScoreTable) score_table = mirbooking_precomputed_score_table_new_from_bytes (g_mapped_file_get_bytes (score_map),
                                                                                                                 seed_offset,
                                                                                                                 seed_length);
    mirbooking_broker_set_score_table (mirbooking,
                                       g_object_ref (score_table));

    if (mirnas_file == NULL)
    {
        g_printerr ("The '--mirnas' argument is required.\n");
        return EXIT_FAILURE;
    }

    if (targets_file == NULL)
    {
        g_printerr ("The '--targets' argument is required.\n");
        return EXIT_FAILURE;
    }


    g_autoptr (FILE) mirnas_f = NULL;
    g_autoptr (FILE) targets_f = NULL;
    g_autoptr (FILE) quantities_f = NULL;
    g_autoptr (FILE) output_f = NULL;

    if ((mirnas_f = g_fopen (mirnas_file, "r")) == NULL)
    {
        g_printerr ("Could not open the miRNAs file '%s': %s.\n", mirnas_file, g_strerror (errno));
        return EXIT_FAILURE;
    }

    if ((targets_f = g_fopen (targets_file, "r")) == NULL)
    {
        g_printerr ("Could not open the targets file '%s': %s.\n", targets_file, g_strerror (errno));
        return EXIT_FAILURE;
    }

    if ((quantities_f = quantities_file == NULL ? stdin : g_fopen (quantities_file, "r")) == NULL)
    {
        g_printerr ("Could not open the quantities file '%s': %s.\n", quantities_file, g_strerror (errno));
        return EXIT_FAILURE;
    }

    if ((output_f = output_file == NULL ? stdout : g_fopen (output_file, "w")) == NULL)
    {
        g_printerr ("Could not open the output file '%s': %s.\n", output_file, g_strerror (errno));
        return EXIT_FAILURE;
    }

    // accession -> #MirbookingSequence
    g_autoptr (GHashTable) sequences_hash = g_hash_table_new_full (g_str_hash,
                                                                   g_str_equal,
                                                                   g_free,
                                                                   g_object_unref);

    // mapped sequences
    g_autoptr (GMappedFile) mirnas_map = g_mapped_file_new_from_fd (fileno (mirnas_f),
                                                                    FALSE,
                                                                    &error);

    if (mirnas_map == NULL)
    {
        g_printerr ("%s (%s, %u).\n", error->message, g_quark_to_string (error->domain), error->code);
        return EXIT_FAILURE;
    }

    g_autoptr (GMappedFile) targets_map = g_mapped_file_new_from_fd (fileno (targets_f),
                                                                     FALSE,
                                                                     &error);

    if (targets_map == NULL)
    {
        g_printerr ("%s (%s, %u).\n", error->message, g_quark_to_string (error->domain), error->code);
        return EXIT_FAILURE;
    }

    // precondition mirnas
    read_sequences_from_fasta (mirnas_f,
                               mirnas_map,
                               TRUE,
                               sequences_hash);

    // precondition targets
    read_sequences_from_fasta (targets_f,
                               targets_map,
                               FALSE,
                               sequences_hash);

    g_autoptr (GHashTable) cds_hash = g_hash_table_new_full (g_str_hash,
                                                             g_str_equal,
                                                             g_free,
                                                             NULL);

    if (cds_regions_file != NULL)
    {
        g_autoptr (FILE) cds_regions_f = g_fopen (cds_regions_file, "r");

        if (cds_regions_f == NULL)
        {
            g_printerr ("Could not open the coding regions file '%s': %s.\n", cds_regions_file, g_strerror (errno));
            return EXIT_FAILURE;
        }

        // extract cds
        gchar line[1024];
        guint lineno = 0;
        while (lineno++, fgets (line, sizeof (line), cds_regions_f))
        {
            if (lineno == 1 && g_str_has_prefix (line, "accession\tregion\n"))
            {
                continue;
            }

            gchar *accession = strtok (line, "\t");
            guint16 cds_start, cds_end;
            if (sscanf (strtok (NULL, "\t"), "%" G_GUINT16_FORMAT ".." "%" G_GUINT16_FORMAT, &cds_start, &cds_end) != 2)
            {
                g_printerr ("Malformed coding region interval for accession '%s' at line %u.\n", accession, lineno);
                return EXIT_FAILURE;
            }
            g_hash_table_insert (cds_hash,
                                 g_strdup (accession),
                                 gpointer_from_two_guint16 (cds_start - 1, cds_end)); /* convert inclusive index in proper slice */
        }
    }

    gchar line[1024];
    guint lineno = 0;
    while (lineno++, fgets (line, sizeof (line), quantities_f))
    {
        if (lineno == 1 && g_str_has_prefix (line, "accession\tquantity\n"))
        {
            continue;
        }

        gchar *accession = strtok (line, "\t");
        gchar *er = NULL;
        gdouble quantity = g_strtod (strtok (NULL, "\n"), &er);

        if (*er != '\0') // strtok replaces the '\n' by '\0'
        {
            g_printerr ("Malformed quantity for accession '%s' at line %u.\n", accession, lineno);
            return EXIT_FAILURE;
        }

        MirbookingSequence *sequence = g_hash_table_lookup (sequences_hash, accession);

        if (sequence == NULL)
        {
            g_printerr ("Unknown sequence with accession '%s'.\n", accession);
            return EXIT_FAILURE;
        }

        mirbooking_broker_set_sequence_quantity (mirbooking,
                                                 sequence,
                                                 quantity);
    }

    g_hash_table_unref (sequences_hash);

    if (!mirbooking_broker_run (mirbooking, &error))
    {
        g_printerr ("%s (%s, %u)\n", error->message, g_quark_to_string (error->domain), error->code);
        return EXIT_FAILURE;
    }

    g_fprintf (output_f, "target_accession\t"
                         "target_name\t"
                         "target_quantity\t"
                         "target_silencing\t"
                         "position\t"
                         "region\t"
                         "occupancy\t"
                         "mirna_accession\t"
                         "mirna_name\t"
                         "mirna_quantity\t"
                         "score\t"
                         "quantity\n");

    GArray *target_sites = mirbooking_broker_get_target_sites (mirbooking);

    gdouble target_silencing = 0;
    gfloat target_quantity = 0;

    const MirbookingTargetSite *target_site;
    MirbookingTarget *cur_target = NULL;
    for (target_site = &g_array_index (target_sites, MirbookingTargetSite, 0);
         target_site < &g_array_index (target_sites, MirbookingTargetSite, target_sites->len);
         target_site++)
    {
        // recompute each time the target changes
        if (cur_target != target_site->target)
        {
            cur_target = target_site->target;
            target_quantity = mirbooking_broker_get_sequence_quantity (mirbooking,
                                                                       MIRBOOKING_SEQUENCE (target_site->target));
            target_silencing = mirbooking_broker_get_target_silencing (mirbooking,
                                                                       target_site->target);
        }

        gfloat occupancy = (target_quantity - mirbooking_broker_get_target_site_vacancy (mirbooking, target_site)) / target_quantity;

        GSList *occupants;
        for (occupants = target_site->occupants; occupants != NULL; occupants = occupants->next)
        {
            gfloat score;
            MirbookingRegion region;
            gpointer cds_ptr;

            MirbookingOccupant *occupant = occupants->data;

            g_autoptr (GError) error = NULL;
            score = mirbooking_score_table_compute_score (MIRBOOKING_SCORE_TABLE (score_table),
                                                          occupant->mirna,
                                                          target_site->target,
                                                          target_site->position,
                                                          &error);
            if (error != NULL)
            {
                g_printerr ("%s (%s, %d)\n", error->message, g_quark_to_string (error->domain), error->code);
                return 1;
            }

            cds_ptr = g_hash_table_lookup (cds_hash,
                                           mirbooking_sequence_get_accession (MIRBOOKING_SEQUENCE (target_site->target)));
            if (cds_ptr != NULL)
            {
                guint16 cds[2];
                two_guint16_from_gpointer (cds_ptr, cds);
                region = mirbooking_region_from_target_site (target_site,
                                                             cds[0],
                                                             cds[1]);
            }
            else
            {
                region = MIRBOOKING_REGION_UNKNOWN;
            }

            #define COALESCE(x,d) (x == NULL ? (d) : (x))
            #define FLOAT_FORMAT "%6f"
            g_fprintf (output_f, "%s\t%s\t" FLOAT_FORMAT "\t" FLOAT_FORMAT "\t%lu\t%s\t" FLOAT_FORMAT "\t%s\t%s\t" FLOAT_FORMAT "\t" FLOAT_FORMAT "\t%u\n",
                       mirbooking_sequence_get_accession (MIRBOOKING_SEQUENCE (target_site->target)),
                       COALESCE (mirbooking_sequence_get_name (MIRBOOKING_SEQUENCE (target_site->target)), "N/A"),
                       target_quantity,
                       target_silencing,
                       target_site->position + 1, // 1-based
                       mirbooking_region_to_string (region),
                       occupancy,
                       mirbooking_sequence_get_accession (MIRBOOKING_SEQUENCE (occupant->mirna)),
                       COALESCE (mirbooking_sequence_get_name (MIRBOOKING_SEQUENCE (occupant->mirna)), "N/A"),
                       mirbooking_broker_get_sequence_quantity (mirbooking, MIRBOOKING_SEQUENCE (occupant->mirna)),
                       score,
                       occupant->quantity);
            #undef COALESCE
        }
    }

    return EXIT_SUCCESS;
}
