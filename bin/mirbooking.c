#include <glib.h>
#include <glib/gprintf.h>
#include <fcntl.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <glib/gstdio.h>
#include <gio/gio.h>
#include <math.h>

#include <mirbooking.h>

G_DEFINE_AUTOPTR_CLEANUP_FUNC (FILE, fclose)

static gchar   *mirnas_file      = NULL;
static gchar   *targets_file     = NULL;
static gchar   *cds_regions_file = NULL;
static gchar   *score_table_file = NULL;
static gchar   *quantities_file  = NULL;
static gchar   *output_file      = NULL;
static gdouble  threshold        = MIRBOOKING_DEFAULT_THRESHOLD;
static gdouble  log_base         = MIRBOOKING_DEFAULT_LOG_BASE;
static gsize    prime5_footprint = MIRBOOKING_DEFAULT_5PRIME_FOOTPRINT;
static gsize    prime3_footprint = MIRBOOKING_DEFAULT_3PRIME_FOOTPRINT;
static gdouble  utr5_multiplier  = 0.1;
static gdouble  cds_multiplier   = 0.1;
static gdouble  utr3_multiplier  = 1.0;

static GOptionEntry MIRBOOKING_OPTION_ENTRIES[] =
{
    {"mirnas",                0, 0, G_OPTION_ARG_FILENAME, &mirnas_file,      "MiRNAs sequences FASTA file",                                                                     "FILE"},
    {"targets",               0, 0, G_OPTION_ARG_FILENAME, &targets_file,     "Transcripts sequences FASTA file",                                                                "FILE"},
    {"cds-regions",           0, 0, G_OPTION_ARG_FILENAME, &cds_regions_file, "Coding regions as a two-column (accession, 1-based inclusive interval) TSV file",                 "FILE"},
    {"score-table",           0, 0, G_OPTION_ARG_FILENAME, &score_table_file, "Precomputed seed-MRE duplex table as a row-major big-endian float matrix file",                   "FILE"},
    {"quantities",            0, 0, G_OPTION_ARG_FILENAME, &quantities_file,  "MiRNA and targets quantities as a two-column (accession, quantity) TSV file (defaults to stdin)", "FILE"},
    {"output",                0, 0, G_OPTION_ARG_FILENAME, &output_file,      "Output destination file (defaults to stdout)",                                                    "FILE"},
    {"threshold",             0, 0, G_OPTION_ARG_DOUBLE,   &threshold,        "Probability threshold for site matching",                                                         G_STRINGIFY (MIRBOOKING_DEFAULT_THRESHOLD)},
    {"log-base",              0, 0, G_OPTION_ARG_DOUBLE,   &log_base,         "Logarithm base for spreading quantites across sites",                                             G_STRINGIFY (MIRBOOKING_DEFAULT_LOG_BASE)},
    {"5prime-footprint",      0, 0, G_OPTION_ARG_INT,      &prime5_footprint, "Footprint in the MRE's 5' direction",                                                             G_STRINGIFY (MIRBOOKING_DEFAULT_5PRIME_FOOTPRINT)},
    {"3prime-footprint",      0, 0, G_OPTION_ARG_INT,      &prime3_footprint, "Footprint in the MRE's 3' direction",                                                             G_STRINGIFY (MIRBOOKING_DEFAULT_3PRIME_FOOTPRINT)},
    {"5prime-utr-multiplier", 0, 0, G_OPTION_ARG_DOUBLE,   &utr5_multiplier,  "Silencing multiplier for the 3'UTR region",                                                       "0.1"},
    {"cds-multiplier",        0, 0, G_OPTION_ARG_DOUBLE,   &cds_multiplier,   "Silencing multiplier for the CDS region",                                                         "0.1"},
    {"3prime-utr-multiplier", 0, 0, G_OPTION_ARG_DOUBLE,   &utr3_multiplier,  "Silencing multiplier for the 5'UTR region",                                                       "1.0"},
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
    gchar *seq;
    gchar line[1024];

    while (fgets (line, sizeof (line), file))
    {
        if (line[0] == '>')
        {
            accession = strtok (line + 1, " ");

            if (accession_in_comment)
            {
                accession = strtok (NULL, " ");
            }

            seq = g_mapped_file_get_contents (mapped_file) + ftell (file);

            gsize remaining = g_mapped_file_get_length (mapped_file) - ftell (file);
            gchar *next_seq = memchr (seq, '>', remaining);
            gsize seq_len = next_seq == NULL ? remaining - 1 : next_seq - seq - 1;

            MirbookingSequence *sequence;

            if (g_str_has_prefix (accession, "MIMAT"))
            {
                sequence = MIRBOOKING_SEQUENCE (mirbooking_mirna_new (accession));
            }
            else
            {
                sequence = MIRBOOKING_SEQUENCE (mirbooking_target_new (accession));
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

static gfloat
mirbooking_region_get_multiplier (MirbookingRegion region)
{
    switch (region)
    {
        case MIRBOOKING_REGION_5PRIME_UTR:
            return utr5_multiplier;
        case MIRBOOKING_REGION_CDS:
            return cds_multiplier;
        case MIRBOOKING_REGION_3PRIME_UTR:
            return utr3_multiplier;
        case MIRBOOKING_REGION_UNKNOWN:
            return 1.0f;
        default:
            g_return_val_if_reached (1.0f);
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
        g_printerr ("%s (%s, %u)\n", error->message, g_quark_to_string (error->domain), error->code);
        return EXIT_FAILURE;
    }

    g_autoptr (Mirbooking) mirbooking = mirbooking_new ();

    mirbooking_set_threshold (mirbooking, threshold);
    mirbooking_set_log_base (mirbooking, log_base);
    mirbooking_set_5prime_footprint (mirbooking, prime5_footprint);
    mirbooking_set_3prime_footprint (mirbooking, prime3_footprint);

    if (score_table_file == NULL)
    {
        g_printerr ("The '--score-table' argument is required.\n");
        return EXIT_FAILURE;
    }

    g_autoptr (GMappedFile) score_map = g_mapped_file_new (score_table_file, FALSE, &error);
    if (score_map == NULL)
    {
        g_printerr ("%s (%s, %u)\n", error->message, g_quark_to_string (error->domain), error->code);
        return EXIT_FAILURE;
    }

    g_autoptr (MirbookingScoreTable) score_table = mirbooking_score_table_new_precomputed (g_mapped_file_get_bytes (score_map));
    mirbooking_set_score_table (mirbooking, score_table);

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

    g_autoptr (FILE) mirnas_f = g_fopen (mirnas_file, "r");
    g_autoptr (FILE) targets_f = g_fopen (targets_file, "r");
    g_autoptr (FILE) quantities_f = quantities_file == NULL ? stdin : g_fopen (quantities_file, "r");
    g_autoptr (FILE) output_f = output_file == NULL ? stdout : g_fopen (output_file, "w");

    g_return_val_if_fail (mirnas_f != NULL, EXIT_FAILURE);
    g_return_val_if_fail (targets_f != NULL, EXIT_FAILURE);
    g_return_val_if_fail (quantities_f != NULL, EXIT_FAILURE);
    g_return_val_if_fail (output_f != NULL, EXIT_FAILURE);

    // accession -> #MirbookingSequence
    g_autoptr (GHashTable) sequences_hash = g_hash_table_new_full (g_str_hash,
                                                                   g_str_equal,
                                                                   g_free,
                                                                   g_object_unref);

    // mapped sequences
    g_autoptr (GMappedFile) mirnas_map = g_mapped_file_new_from_fd (fileno (mirnas_f),
                                                                    FALSE,
                                                                    NULL);

    if (mirnas_map == NULL)
    {
        g_printerr ("Could not map the miRNAs file.\n");
        return EXIT_FAILURE;
    }

    g_autoptr (GMappedFile) targets_map = g_mapped_file_new_from_fd (fileno (targets_f),
                                                                     FALSE,
                                                                     NULL);

    if (targets_map == NULL)
    {
        g_printerr ("Could not map the targets file.\n");
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
            g_printerr ("Could not open the CDS regions file.\n");
            return EXIT_FAILURE;
        }

        // extract cds
        gchar line[1024];
        while (fgets (line, sizeof (line), cds_regions_f))
        {
            gchar *accession = strtok (line, "\t");
            guint16 cds_start, cds_end;
            sscanf (strtok (NULL, "\t"), "%" G_GUINT16_FORMAT ".." "%" G_GUINT16_FORMAT, &cds_start, &cds_end);
            g_hash_table_insert (cds_hash,
                                 g_strdup (accession),
                                 gpointer_from_two_guint16 (cds_start - 1, cds_end)); /* convert inclusive index in proper slice */
        }
    }

    gfloat total_mirna_quantity;
    gfloat total_target_quantity;

    gchar line[1024];
    while (fgets (line, sizeof (line), quantities_f))
    {
        gchar *accession = strtok (line, "\t");
        gdouble quantity = g_strtod (strtok (NULL, "\t"), NULL);

        MirbookingSequence *sequence = g_hash_table_lookup (sequences_hash, accession);

        if (sequence == NULL)
        {
            g_printerr ("Unknown sequence with accession '%s'.\n", accession);
            continue;
        }

        mirbooking_set_sequence_quantity (mirbooking,
                                          g_object_ref (sequence),
                                          quantity);

        if (MIRBOOKING_IS_MIRNA (sequence))
        {
            total_mirna_quantity += quantity;
        }

        if (MIRBOOKING_IS_TARGET (sequence))
        {
            total_target_quantity += quantity;
        }
    }

    g_hash_table_unref (sequences_hash);

    if (ABS (logf (total_target_quantity) - logf (total_mirna_quantity)) >= 1.0f)
    {
        g_warning ("The quantity of mirnas %f is not in the same scale as the quantity of target %f.",
                   total_mirna_quantity,
                   total_target_quantity);
    }

    if (!mirbooking_run (mirbooking, &error))
    {
        g_printerr ("%s (%s, %u)\n", error->message, g_quark_to_string (error->domain), error->code);
        return EXIT_FAILURE;
    }

    g_fprintf (output_f, "Target Accession\tMiRNA Accession\tPosition\tLocation\tProbability\tOccupancy\tSilencing\n");

    gsize target_sites_len;
    const MirbookingTargetSite *target_sites = mirbooking_get_target_sites (mirbooking,
                                                                            &target_sites_len);

    const MirbookingTargetSite *target_site = target_sites;
    while (target_site < target_sites + target_sites_len)
    {
        GSList *occupants;
        for (occupants = target_site->occupants; occupants != NULL; occupants = occupants->next)
        {
            gfloat probability;
            MirbookingRegion region;
            gpointer cds_ptr;

            MirbookingOccupant *occupant = occupants->data;

            probability = mirbooking_score_table_compute_score (score_table,
                                                                MIRBOOKING_SEQUENCE (occupant->mirna),
                                                                1,
                                                                MIRBOOKING_SEQUENCE (target_site->target),
                                                                target_site->position,
                                                                7);

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

            g_fprintf (output_f, "%s\t%s\t%ld\t%s\t%f\t%d\t%f\n",
                       mirbooking_sequence_get_accession (MIRBOOKING_SEQUENCE (target_site->target)),
                       mirbooking_sequence_get_accession (MIRBOOKING_SEQUENCE (occupant->mirna)),
                       target_site->position + 1, // 1-based
                       mirbooking_region_to_string (region),
                       probability,
                       occupant->quantity,
                       occupant->quantity * probability * mirbooking_region_get_multiplier (region));
        }
        ++target_site;
    }

    return EXIT_SUCCESS;
}
