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
#if HAVE_MPI
#include <mpi.h>
#endif
#include <stdarg.h>
#include <string.h>

#include <mirbooking.h>

G_DEFINE_AUTOPTR_CLEANUP_FUNC (FILE, fclose)

#define MIRBOOKING_DEFAULT_MAX_ITERATIONS 100
#define MIRBOOKING_DEFAULT_CUTOFF         100 // pM
#define MIRBOOKING_DEFAULT_REL_CUTOFF     0.0

typedef enum _MirbookingOutputFormat
{
    MIRBOOKING_OUTPUT_FORMAT_TSV,
    MIRBOOKING_OUTPUT_FORMAT_TSV_DETAILED,
    MIRBOOKING_OUTPUT_FORMAT_GFF3,
    MIRBOOKING_OUTPUT_FORMAT_WIG
} MirbookingOutputFormat;

typedef struct _MirbookingOutputFormatMeta
{
    MirbookingOutputFormat  output_format;
    const gchar            *nick;
    void (*write) (MirbookingBroker *broker, FILE *output);
} MirbookingOutputFormatMeta;

static void write_output_to_tsv (MirbookingBroker *broker, FILE *out);
static void write_output_to_tsv_detailed (MirbookingBroker *broker, FILE *out);
static void write_output_to_gff3 (MirbookingBroker *broker, FILE *out);
static void write_output_to_wiggle (MirbookingBroker *broker, FILE *out);

const MirbookingOutputFormatMeta MIRBOOKING_OUTPUT_FORMAT_META[] =
{
    {MIRBOOKING_OUTPUT_FORMAT_TSV,          "tsv",          write_output_to_tsv},
    {MIRBOOKING_OUTPUT_FORMAT_TSV_DETAILED, "tsv-detailed", write_output_to_tsv_detailed},
    {MIRBOOKING_OUTPUT_FORMAT_GFF3,         "gff3",         write_output_to_gff3},
    {MIRBOOKING_OUTPUT_FORMAT_WIG,          "wig",          write_output_to_wiggle}
};

static gchar                        **targets_files             = {NULL};
static gchar                        **mirnas_files              = {NULL};
static gchar                         *seed_scores_file          = NULL;
static gchar                         *supplementary_scores_file = NULL;
static gchar                         *accessibility_scores_file = NULL;
static gchar                         *input_file                = NULL;
static gchar                         *output_file               = NULL;
static MirbookingOutputFormat         output_format             = MIRBOOKING_OUTPUT_FORMAT_TSV;
static MirbookingBrokerSparseSolver   sparse_solver             = MIRBOOKING_BROKER_DEFAULT_SPARSE_SOLVER;
static guint64                        max_iterations            = MIRBOOKING_DEFAULT_MAX_ITERATIONS;
static gsize                          prime5_footprint          = MIRBOOKING_BROKER_DEFAULT_5PRIME_FOOTPRINT;
static gsize                          prime3_footprint          = MIRBOOKING_BROKER_DEFAULT_3PRIME_FOOTPRINT;
static gdouble                        cutoff                    = MIRBOOKING_DEFAULT_CUTOFF;
static gdouble                        rel_cutoff                = MIRBOOKING_DEFAULT_REL_CUTOFF;

static MirbookingDefaultScoreTableSupplementaryModel supplementary_model = MIRBOOKING_DEFAULT_SCORE_TABLE_DEFAULT_SUPPLEMENTARY_MODEL;

static gboolean
set_supplementary_model (const gchar   *key,
                         const gchar   *value,
                         gpointer       data,
                         GError      **error)
{
    if (g_strcmp0 (value, "wee-et-al-2012") == 0)
    {
        supplementary_model = MIRBOOKING_DEFAULT_SCORE_TABLE_SUPPLEMENTARY_MODEL_WEE_ET_AL_2012;
    }
    else if (g_strcmp0 (value, "yan-et-al-2018") == 0)
    {
        supplementary_model = MIRBOOKING_DEFAULT_SCORE_TABLE_SUPPLEMENTARY_MODEL_YAN_ET_AL_2018;
    }
    else
    {
        g_set_error (error,
                     G_OPTION_ERROR,
                     G_OPTION_ERROR_BAD_VALUE,
                     "No such supplementary model %s.", value);
        return FALSE;
    }

    return TRUE;
}

static gboolean
set_output_format (const gchar   *key,
                   const gchar   *value,
                   gpointer       data,
                   GError      **error)
{
    gint i;
    for (i = 0; i < sizeof (MIRBOOKING_OUTPUT_FORMAT_META) / sizeof (MIRBOOKING_OUTPUT_FORMAT_META[0]); i++)
    {
        if (g_strcmp0 (value, MIRBOOKING_OUTPUT_FORMAT_META[i].nick) == 0)
        {
            output_format = MIRBOOKING_OUTPUT_FORMAT_META[i].output_format;
            return TRUE;
        }
    }

    g_set_error (error,
                 G_OPTION_ERROR,
                 G_OPTION_ERROR_BAD_VALUE,
                 "No such output format %s.", value);
    return FALSE;
}

static gboolean
set_sparse_solver(const gchar   *key,
                  const gchar   *value,
                  gpointer       data,
                  GError      **error)
{
    g_autoptr (GEnumClass) sparse_solver_class = g_type_class_ref (MIRBOOKING_BROKER_SPARSE_SOLVER_ENUM);
    GEnumValue *eval = g_enum_get_value_by_nick (sparse_solver_class,
                                                 value);

    if (eval == NULL || !mirbooking_broker_sparse_solver_is_available (eval->value))
    {
        g_set_error (error,
                     G_OPTION_ERROR,
                     G_OPTION_ERROR_BAD_VALUE,
                     "No solver %s available.", value);
        return FALSE;
    }

    sparse_solver = eval->value;

    return TRUE;
}

static GOptionEntry MIRBOOKING_OPTION_ENTRIES[] =
{
    {"targets",              0, 0, G_OPTION_ARG_FILENAME_ARRAY, &targets_files,             "Targets FASTA files",                                                                             NULL},
    {"mirnas",               0, 0, G_OPTION_ARG_FILENAME_ARRAY, &mirnas_files,              "miRNA FASTA files",                                                                               NULL},
    {"seed-scores",          0, 0, G_OPTION_ARG_FILENAME,       &seed_scores_file,          "Precomputed seed::MRE binding free energy duplex table",                                          "FILE"},
    {"supplementary-model",  0, 0, G_OPTION_ARG_CALLBACK,       &set_supplementary_model,   "Supplementary bindings model to use",                                                             "yan-et-al-2018"},
    {"supplementary-scores", 0, 0, G_OPTION_ARG_FILENAME,       &supplementary_scores_file, "Precomputed supplementary::MRE binding free energy duplex table",                                 "FILE"},
    {"accessibility-scores", 0, 0, G_OPTION_ARG_FILENAME,       &accessibility_scores_file, "Accessibility scores as a variable columns (accession, positions...) TSV file",                   "FILE"},
    {"input",                0, 0, G_OPTION_ARG_FILENAME,       &input_file,                "MiRNA and targets quantities as a two-column (accession, quantity) TSV file (defaults to stdin)", "FILE"},
    {"output",               0, 0, G_OPTION_ARG_FILENAME,       &output_file,               "Output destination file (defaults to stdout)",                                                    "FILE"},
    {"output-format",        0, 0, G_OPTION_ARG_CALLBACK,       &set_output_format,         "Output format (i.e. 'tsv', 'gff3')",                                                              "tsv"},
    {"sparse-solver",        0, 0, G_OPTION_ARG_CALLBACK,       &set_sparse_solver,         "Sparse solver implementation to use",                                                             "superlu"},
    {"max-iterations",       0, 0, G_OPTION_ARG_INT,            &max_iterations,            "Maximum number of iterations",                                                                    G_STRINGIFY (MIRBOOKING_DEFAULT_MAX_ITERATIONS)},
    {"5prime-footprint",     0, 0, G_OPTION_ARG_INT,            &prime5_footprint,          "Footprint in the MRE's 5' direction",                                                             G_STRINGIFY (MIRBOOKING_BROKER_DEFAULT_5PRIME_FOOTPRINT)},
    {"3prime-footprint",     0, 0, G_OPTION_ARG_INT,            &prime3_footprint,          "Footprint in the MRE's 3' direction",                                                             G_STRINGIFY (MIRBOOKING_BROKER_DEFAULT_3PRIME_FOOTPRINT)},
    {"cutoff",               0, 0, G_OPTION_ARG_DOUBLE,         &cutoff,                    "Cutoff on the duplex concentration",                                                              G_STRINGIFY (MIRBOOKING_DEFAULT_CUTOFF)},
    {"relative-cutoff",      0, 0, G_OPTION_ARG_DOUBLE,         &rel_cutoff,                "Relative cutoff on the bound fraction",                                                           G_STRINGIFY (MIRBOOKING_DEFAULT_REL_CUTOFF)},
    {0}
};

typedef enum _FastaFormat
{
    FASTA_FORMAT_GENERIC,
    FASTA_FORMAT_NCBI,
    FASTA_FORMAT_GENCODE,
    FASTA_FORMAT_MIRBASE
} FastaFormat;

static FastaFormat
detect_fasta_format (const gchar *path)
{
    g_autofree const gchar *path_basename = g_path_get_basename (path);

    if (g_str_has_prefix (path_basename, "gencode."))
    {
        return FASTA_FORMAT_GENCODE;
    }
    else if (g_str_has_prefix (path_basename, "mature"))
    {
        return FASTA_FORMAT_MIRBASE;
    }
    else if (g_str_has_prefix (path_basename, "GCF_") || g_str_has_prefix (path_basename, "GCA_"))
    {
        return FASTA_FORMAT_NCBI;
    }
    else
    {
        return FASTA_FORMAT_GENERIC;
    }
}

static void
read_sequences_from_fasta (FILE        *file,
                           GMappedFile *mapped_file,
                           FastaFormat  fasta_format,
                           GHashTable  *sequences_hash,
                           gboolean     is_mirna)
{
    gchar *accession;
    gchar *name;
    gchar name_buffer[128]; /* in case the name does not appear literally */
    gchar *gene_accession = NULL;
    gchar *gene_name = NULL;
    gchar *seq;
    gchar line[1024];

    while (fgets (line, sizeof (line), file))
    {
        if (line[0] == '>')
        {
            if (fasta_format == FASTA_FORMAT_MIRBASE)
            {
                name      = strtok (line + 1, " ");
                accession = strtok (NULL, " ");
            }
            else if (fasta_format == FASTA_FORMAT_GENCODE)
            {
                // for GENCODE-style annotation, the name is in the sixth field
                accession = strtok (line + 1, "|");
                gene_accession = strtok (NULL, "|");

                gint i;
                for (i = 0; i < 3; i++)
                {
                    name = strtok (NULL, "|");
                }

                gene_name = strtok (NULL, "|");
            }
            else if (fasta_format == FASTA_FORMAT_NCBI)
            {
                accession = strtok (line + 1, " ");

                gchar *name_p = NULL;
                if (strtok (NULL, "(,") && (name_p = strtok (NULL, "),")))
                {
                    // if we parse the ")," following the gene name, name_p
                    // will be NULL and we should look for the transcript
                    // variant
                    if (name_p)
                    {
                        gene_name = name_p;
                    }
                }

                /* unfortunately, that's the best information we get from the FASTA */
                gene_accession = gene_name;

                /* construct the gene name with its variant number */
                guint variant = 0;
                sscanf (strtok (NULL, ","), " transcript variant %u", &variant);
                g_sprintf (name_buffer, "%s-%03u", gene_name, 200 + variant);
                name = name_buffer;
            }
            else
            {
                accession = strtok (line + 1, " ");
                name = strtok (NULL, "\n");
            }

            glong offset = ftell (file);
            g_return_if_fail (offset != -1);

            seq = g_mapped_file_get_contents (mapped_file) + offset;

            gsize remaining = g_mapped_file_get_length (mapped_file) - offset;
            gchar *next_seq = memchr (seq, '>', remaining);
            gsize seq_len = next_seq == NULL ? remaining - 1 : (gsize) (next_seq - seq) - 1;

            g_autoptr (MirbookingSequence) sequence = NULL;

            if (is_mirna)
            {
                sequence = MIRBOOKING_SEQUENCE (mirbooking_mirna_new_with_name (accession, name));
            }
            else
            {
                sequence = MIRBOOKING_SEQUENCE (mirbooking_target_new_with_name (accession, name));
            }

            // FIXME: add a destroy notify to clear this ref when the sequence is disposed
            g_mapped_file_ref (mapped_file);

            mirbooking_sequence_set_gene_accession (sequence, gene_accession);
            mirbooking_sequence_set_gene_name (sequence, gene_name);

            mirbooking_sequence_set_raw_sequence (sequence,
                                                  g_bytes_new_from_bytes (g_mapped_file_get_bytes (mapped_file), offset, seq_len));

            g_hash_table_insert (sequences_hash,
                                 g_strdup (accession),
                                 g_steal_pointer (&sequence));
        }
    }
}

gboolean
read_sequence_accessibility (GInputStream     *in,
                             GHashTable       *sequences_hash,
                             GError          **error)
{
    g_autoptr (GError) err = NULL;
    g_autoptr (GDataInputStream) data_in = g_data_input_stream_new (in);

    while (TRUE)
    {
        g_autofree gchar *accession = g_data_input_stream_read_upto (data_in,
                                                                     "\t",
                                                                     1,
                                                                     NULL,
                                                                     NULL,
                                                                     &err);

        if (accession == NULL && err)
        {
            g_propagate_error (error, err);
            return FALSE;
        }

        if (accession == NULL)
        {
            break; /* done */
        }

        MirbookingSequence *sequence = g_hash_table_lookup (sequences_hash,
                                                            accession);

        if (sequence != NULL)
        {
            gsize sequence_len = mirbooking_sequence_get_sequence_length (sequence);
            gsize position;
            for (position = 0; position < sequence_len; position++)
            {
                guchar ret = g_data_input_stream_read_byte (data_in, NULL, &err);
                if (!ret && err)
                {
                    g_propagate_error (error, err);
                    return FALSE;
                }

                g_return_val_if_fail (ret == '\t', FALSE);

                g_autofree gchar *score_str = g_data_input_stream_read_upto (data_in, "\t", 1, NULL, NULL, &err);

                if (score_str == NULL && err)
                {
                    g_propagate_error (error, err);
                    return FALSE;
                }

                g_return_val_if_fail (score_str != NULL, FALSE);

                gfloat score;
                if (*score_str == '\0')
                {
                    score = INFINITY;
                }
                else
                {
                    sscanf (score_str, "%f", &score);
                }

                mirbooking_target_set_accessibility_score (MIRBOOKING_TARGET (sequence),
                                                           position,
                                                           score);
            }
        }

        // consume the rest of the line
        g_autofree gchar *remaining = g_data_input_stream_read_line (data_in, NULL, NULL, &err);

        if (remaining == NULL && err)
        {
            g_propagate_error (error, err);
            return FALSE;
        }
    }

    return TRUE;
}

#define COALESCE(x,d) (x == NULL ? (d) : (x))

static void
write_output_to_tsv (MirbookingBroker *mirbooking,
                     FILE             *output_f)
{
    g_fprintf (output_f, "gene_accession\t"
                         "gene_name\t"
                         "target_accession\t"
                         "target_name\t"
                         "target_quantity\t"
                         "position\t"
                         "mirna_accession\t"
                         "mirna_name\t"
                         "mirna_quantity\t"
                         "score\t"
                         "quantity\n");

    const GArray *target_sites = mirbooking_broker_get_target_sites (mirbooking);

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
        }

        // report individual occupants
        GSList *occupants;
        for (occupants = target_site->occupants; occupants != NULL; occupants = occupants->next)
        {
            MirbookingOccupant *occupant = occupants->data;
            g_fprintf (output_f, "%s\t%s\t%s\t%s\t%e\t%lu\t%s\t%s\t%e\t%e\t%e\n",
                       COALESCE (mirbooking_sequence_get_gene_accession (MIRBOOKING_SEQUENCE (target_site->target)), "N/A"),
                       COALESCE (mirbooking_sequence_get_gene_name (MIRBOOKING_SEQUENCE (target_site->target)), "N/A"),
                       mirbooking_sequence_get_accession (MIRBOOKING_SEQUENCE (target_site->target)),
                       COALESCE (mirbooking_sequence_get_name (MIRBOOKING_SEQUENCE (target_site->target)), "N/A"),
                       target_quantity,
                       target_site->position + 1, // 1-based
                       mirbooking_sequence_get_accession (MIRBOOKING_SEQUENCE (occupant->mirna)),
                       COALESCE (mirbooking_sequence_get_name (MIRBOOKING_SEQUENCE (occupant->mirna)), "N/A"),
                       mirbooking_broker_get_sequence_quantity (mirbooking, MIRBOOKING_SEQUENCE (occupant->mirna)) + mirbooking_broker_get_bound_mirna_quantity (mirbooking, occupant->mirna),
                       MIRBOOKING_SCORE_KM (occupant->score) + (mirbooking_broker_get_target_site_kother (mirbooking, target_site) / occupant->score.kf),
                       mirbooking_broker_get_occupant_quantity (mirbooking, occupant));
        }
    }
}

static void
write_output_to_tsv_detailed (MirbookingBroker *mirbooking,
                     FILE             *output_f)
{
    g_fprintf (output_f, "gene_accession\t"
                         "gene_name\t"
                         "target_accession\t"
                         "target_name\t"
                         "target_quantity\t"
                         "position\t"
                         "mirna_accession\t"
                         "mirna_name\t"
                         "mirna_quantity\t"
                         "kf\t"
                         "kr\t"
                         "kcat\t"
                         "kother\t"
                         "kd\t"
                         "km\t"
                         "quantity\n");

    const GArray *target_sites = mirbooking_broker_get_target_sites (mirbooking);

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
        }

        // report individual occupants
        GSList *occupants;
        for (occupants = target_site->occupants; occupants != NULL; occupants = occupants->next)
        {
            MirbookingOccupant *occupant = occupants->data;
            g_fprintf (output_f, "%s\t%s\t%s\t%s\t%e\t%lu\t%s\t%s\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",
                       COALESCE (mirbooking_sequence_get_gene_accession (MIRBOOKING_SEQUENCE (target_site->target)), "N/A"),
                       COALESCE (mirbooking_sequence_get_gene_name (MIRBOOKING_SEQUENCE (target_site->target)), "N/A"),
                       mirbooking_sequence_get_accession (MIRBOOKING_SEQUENCE (target_site->target)),
                       COALESCE (mirbooking_sequence_get_name (MIRBOOKING_SEQUENCE (target_site->target)), "N/A"),
                       target_quantity,
                       target_site->position + 1, // 1-based
                       mirbooking_sequence_get_accession (MIRBOOKING_SEQUENCE (occupant->mirna)),
                       COALESCE (mirbooking_sequence_get_name (MIRBOOKING_SEQUENCE (occupant->mirna)), "N/A"),
                       mirbooking_broker_get_sequence_quantity (mirbooking, MIRBOOKING_SEQUENCE (occupant->mirna)) + mirbooking_broker_get_bound_mirna_quantity (mirbooking, occupant->mirna),
                       occupant->score.kf,
                       occupant->score.kr,
                       occupant->score.kcat,
                       mirbooking_broker_get_target_site_kother (mirbooking, target_site),
                       MIRBOOKING_SCORE_KD (occupant->score),
                       MIRBOOKING_SCORE_KM (occupant->score) + (mirbooking_broker_get_target_site_kother (mirbooking, target_site) / occupant->score.kf),
                       mirbooking_broker_get_occupant_quantity (mirbooking, occupant));
        }
    }
}

static void
write_output_to_gff3 (MirbookingBroker *mirbooking, FILE *output_f)
{
    g_fprintf (output_f, "##gff-version 3\n");

    const GArray *target_sites = mirbooking_broker_get_target_sites (mirbooking);

    gint i = 1;

    const MirbookingTargetSite *target_site;
    for (target_site = &g_array_index (target_sites, MirbookingTargetSite, 0);
         target_site < &g_array_index (target_sites, MirbookingTargetSite, target_sites->len);
         target_site++)
    {
        // report individual occupants
        GSList *occupants;
        for (occupants = target_site->occupants; occupants != NULL; occupants = occupants->next)
        {
            // the sequence ontology for 'miRNA_target_site' is 'SO:0000934'
            MirbookingOccupant *occupant = occupants->data;
            g_fprintf (output_f, "%s\tmiRBooking\tmiRNA_target_site\t%lu\t%lu\t%e\t.\t.\tID=%d;Name=%s;Alias=%s\n",
                       mirbooking_sequence_get_accession (MIRBOOKING_SEQUENCE (target_site->target)),
                       (gsize) MAX (1, (gssize) target_site->position + 1 - (gssize) prime5_footprint),
                       MIN (mirbooking_sequence_get_sequence_length (MIRBOOKING_SEQUENCE (target_site->target)), target_site->position + 1 + prime3_footprint),
                       mirbooking_broker_get_occupant_quantity (mirbooking, occupant),
                       i++,
                       mirbooking_sequence_get_name (MIRBOOKING_SEQUENCE (occupant->mirna)),
                       mirbooking_sequence_get_accession (MIRBOOKING_SEQUENCE (occupant->mirna)));
        }
    }
}

static void
write_output_to_wiggle (MirbookingBroker *broker, FILE *output_f)
{

    const GArray *target_sites = mirbooking_broker_get_target_sites (broker);

    g_fprintf (output_f, "track type=wiggle_0\n");

    MirbookingTarget *target = NULL;
    const MirbookingTargetSite *target_site;
    for (target_site = &g_array_index (target_sites, MirbookingTargetSite, 0);
         target_site < &g_array_index (target_sites, MirbookingTargetSite, target_sites->len);
         target_site++)
    {
        if (target != target_site->target)
        {
            g_fprintf (output_f, "variableStep chrom=%s\n",
                       mirbooking_sequence_get_accession (MIRBOOKING_SEQUENCE (target_site->target)));
            target = target_site->target;
        }

        gdouble St  = mirbooking_broker_get_sequence_quantity (broker, MIRBOOKING_SEQUENCE (target_site->target));
        gdouble Stp = mirbooking_broker_get_target_site_occupants_quantity (broker, target_site);

        // only report positions with activity
        if (Stp > 0)
        {
            g_fprintf (output_f, "%lu %f\n",
                       target_site->position + 1,
                       Stp / St);
        }
    }
}

static void
free_cutoff_filter_user_data (gpointer ud)
{
    MirbookingDefaultScoreTableCutoffFilterUserData *cutoff_filter_ud = ud;
    g_object_unref (cutoff_filter_ud->broker);
}

int
main (gint argc, gchar **argv)
{
    gint rank = 0;
#if HAVE_MPI
    gint provided;
    g_return_val_if_fail (MPI_Init_thread (&argc, &argv, MPI_THREAD_FUNNELED, &provided) == MPI_SUCCESS, EXIT_FAILURE);
    g_return_val_if_fail (provided == MPI_THREAD_FUNNELED, EXIT_FAILURE);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
#endif

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

    g_autoptr (MirbookingBroker) mirbooking = mirbooking_broker_new_with_rank (rank);

    mirbooking_broker_set_5prime_footprint (mirbooking, prime5_footprint);
    mirbooking_broker_set_3prime_footprint (mirbooking, prime3_footprint);
    mirbooking_broker_set_sparse_solver (mirbooking, sparse_solver);

    if (seed_scores_file == NULL)
    {
        g_printerr ("The '--seed-scores' argument is required.\n");
        return EXIT_FAILURE;
    }

    g_autoptr (GMappedFile) seed_scores_map = g_mapped_file_new (seed_scores_file, FALSE, &error);

    if (seed_scores_map == NULL)
    {
        g_printerr ("%s (%s, %u).\n", error->message, g_quark_to_string (error->domain), error->code);
        return EXIT_FAILURE;
    }

    g_autoptr (GBytes) seed_scores_map_bytes = g_mapped_file_get_bytes (seed_scores_map);

    g_return_val_if_fail (seed_scores_map_bytes != NULL, EXIT_FAILURE);

    gsize seed_scores_map_bytes_len;
    const gsize *seed_scores_data = g_bytes_get_data (seed_scores_map_bytes,
                                                      &seed_scores_map_bytes_len);

    // we require at least 'n' and 'nnz'
    if (seed_scores_map_bytes_len < 2 * sizeof (gsize) || (1l << 2 * 7) != *seed_scores_data)
    {
        g_printerr ("The specified seed scores file is invalid.\n");
        return EXIT_FAILURE;
    }

    g_autoptr (GMappedFile) supplementary_scores_file_map = NULL;
    g_autoptr (GBytes) supplementary_scores_map_bytes = NULL;
    if (supplementary_scores_file != NULL)
    {
        supplementary_scores_file_map = g_mapped_file_new (supplementary_scores_file, TRUE, &error);
        if (supplementary_scores_file_map == NULL)
        {
            g_printerr ("%s (%s, %u).\n", error->message, g_quark_to_string (error->domain), error->code);
            return EXIT_FAILURE;
        }

        supplementary_scores_map_bytes = g_mapped_file_get_bytes (supplementary_scores_file_map);
    }

    g_autoptr (MirbookingDefaultScoreTable) score_table = mirbooking_default_score_table_new (seed_scores_map_bytes,
                                                                                              supplementary_model,
                                                                                              supplementary_scores_map_bytes);

    MirbookingDefaultScoreTableCutoffFilterUserData user_data = {
        .broker          = g_object_ref (mirbooking),
        .cutoff          = cutoff,
        .relative_cutoff = rel_cutoff
    };

    mirbooking_default_score_table_set_filter (score_table,
                                               mirbooking_default_score_table_cutoff_filter,
                                               &user_data,
                                               free_cutoff_filter_user_data);

    mirbooking_broker_set_score_table (mirbooking,
                                       MIRBOOKING_SCORE_TABLE (score_table));

    g_autoptr (FILE) input_f = NULL;
    g_autoptr (FILE) output_f = NULL;

    if ((input_f = input_file == NULL ? stdin : g_fopen (input_file, "r")) == NULL)
    {
        g_printerr ("Could not open the quantities file '%s': %s.\n", input_file, g_strerror (errno));
        return EXIT_FAILURE;
    }

    if (rank == 0)
    {
        if ((output_f = output_file == NULL ? stdout : g_fopen (output_file, "w")) == NULL)
        {
            g_printerr ("Could not open the output file '%s': %s.\n", output_file, g_strerror (errno));
            return EXIT_FAILURE;
        }
    }

    // accession -> #MirbookingSequence
    g_autoptr (GHashTable) sequences_hash = g_hash_table_new_full (g_str_hash,
                                                                   g_str_equal,
                                                                   g_free,
                                                                   g_object_unref);

    // precondition all sequences
    gchar **cur_sequences_file = NULL;
    if (targets_files != NULL)
    {
        for (cur_sequences_file = targets_files; *cur_sequences_file != NULL; cur_sequences_file++)
        {
            g_autoptr (FILE) seq_f = g_fopen (*cur_sequences_file, "r");

            if (seq_f == NULL)
            {
                g_printerr ("Could not open the sequences file '%s': %s.\n", *cur_sequences_file, g_strerror (errno));
                return EXIT_FAILURE;
            }

            gint seq_fileno = fileno (seq_f);
            g_return_val_if_fail (seq_fileno != -1, EXIT_FAILURE);

            g_autoptr (GMappedFile) seq_map = g_mapped_file_new_from_fd (seq_fileno, FALSE, &error);

            if (seq_map == NULL)
            {
                g_printerr ("%s (%s, %u).\n", error->message, g_quark_to_string (error->domain), error->code);
                return EXIT_FAILURE;
            }

            read_sequences_from_fasta (seq_f,
                                       seq_map,
                                       detect_fasta_format (*cur_sequences_file),
                                       sequences_hash,
                                       FALSE);
        }
    }

    if (mirnas_files != NULL)
    {
        for (cur_sequences_file = mirnas_files; *cur_sequences_file != NULL; cur_sequences_file++)
        {
            g_autoptr (FILE) seq_f = g_fopen (*cur_sequences_file, "r");

            if (seq_f == NULL)
            {
                g_printerr ("Could not open the sequences file '%s': %s.\n", *cur_sequences_file, g_strerror (errno));
                return EXIT_FAILURE;
            }

            gint seq_fileno = fileno (seq_f);
            g_return_val_if_fail (seq_fileno != -1, EXIT_FAILURE);

            g_autoptr (GMappedFile) seq_map = g_mapped_file_new_from_fd (seq_fileno, FALSE, &error);

            if (seq_map == NULL)
            {
                g_printerr ("%s (%s, %u).\n", error->message, g_quark_to_string (error->domain), error->code);
                return EXIT_FAILURE;
            }

            read_sequences_from_fasta (seq_f,
                                       seq_map,
                                       detect_fasta_format (*cur_sequences_file),
                                       sequences_hash,
                                       TRUE);
        }
    }

    gchar line[1024];
    guint lineno = 0;
    while (lineno++, fgets (line, sizeof (line), input_f))
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

        if (quantity >= cutoff)
        {
            mirbooking_broker_set_sequence_quantity (mirbooking,
                                                     sequence,
                                                     quantity);
        }
        else
        {
            // clear unused entries immediatly for reducing the work
            // of read_sequence_accessibility
            // however, this does not clean sequences that were not quantified
            // in the first place
            g_hash_table_remove (sequences_hash, accession);
        }
    }

    // mark targets with  scores
    if (accessibility_scores_file != NULL)
    {
        g_autoptr (GFile) accessibility_scores_f = g_file_new_for_path (accessibility_scores_file);
        g_autoptr (GInputStream) accessibility_scores_in = G_INPUT_STREAM (g_file_read (accessibility_scores_f,
                                                                                        NULL,
                                                                                        &error));

        if (accessibility_scores_in == NULL)
        {
            g_printerr ("Could not open the accessibility scores file '%s': %s.\n", accessibility_scores_file, g_strerror (errno));
            return EXIT_FAILURE;
        }

        if (g_str_has_suffix (accessibility_scores_file, ".gz"))
        {
            g_autoptr (GZlibDecompressor) gzip_decompressor = g_zlib_decompressor_new (G_ZLIB_COMPRESSOR_FORMAT_GZIP);
            accessibility_scores_in = g_converter_input_stream_new (accessibility_scores_in,
                                                                    G_CONVERTER (gzip_decompressor));
        }

        guint64 accessibility_read_begin = g_get_monotonic_time ();

        // update sequence
        if (!read_sequence_accessibility (G_INPUT_STREAM (accessibility_scores_in),
                                          sequences_hash,
                                          &error))

        {
            if (error)
            {
                g_printerr ("Could not parse accessibility scores: %s (%s, %u).\n", error->message, g_quark_to_string (error->domain), error->code);
            }
            else
            {
                g_printerr ("Could not parse accessibility scores.\n");
            }
            return EXIT_FAILURE;
        }

        g_debug ("Done reading accessibility scores in %lums", 1000 * (g_get_monotonic_time () - accessibility_read_begin) / G_USEC_PER_SEC);
    }

    guint64 iteration = 0;
    guint64 iteration_begin, iteration_end, evaluate_begin, evaluate_end, step_begin, step_end;
    gdouble error_ratio;
    do
    {
        iteration_begin = g_get_monotonic_time ();

        iteration++;

        if (iteration > max_iterations)
        {
            g_print ("Did not converge after %lu iterations.", iteration);
            return EXIT_FAILURE;
        }

        evaluate_begin = g_get_monotonic_time ();

        if (!mirbooking_broker_evaluate (mirbooking,
                                         &error_ratio,
                                         &error))
        {
            g_printerr ("%s (%s, %u)\n", error->message, g_quark_to_string (error->domain), error->code);
            return EXIT_FAILURE;
        }

        evaluate_end = g_get_monotonic_time ();

        if (error_ratio <= 1)
        {
            break;
        }

        step_begin = g_get_monotonic_time ();

        if (!mirbooking_broker_step (mirbooking,
                                     MIRBOOKING_BROKER_STEP_MODE_SOLVE_STEADY_STATE,
                                     1,
                                     &error))
        {
            g_printerr ("%s (%s, %u)\n", error->message, g_quark_to_string (error->domain), error->code);
            return EXIT_FAILURE;
        }

        step_end = g_get_monotonic_time ();

        iteration_end = g_get_monotonic_time ();

        if (rank == 0)
        {
            g_debug ("iteration: %lu error-ratio: %.2e evaluate-time: %lums step-time: %lums total-time: %lums",
                     iteration,
                     error_ratio,
                     1000 * (evaluate_end - evaluate_begin) / G_USEC_PER_SEC,
                     1000 * (step_end - step_begin) / G_USEC_PER_SEC,
                     1000 * (iteration_end - iteration_begin) / G_USEC_PER_SEC);
        }
    }
    while (TRUE);

    if (rank == 0)
    {
        gint i;
        for (i = 0; i < sizeof (MIRBOOKING_OUTPUT_FORMAT_META); i++)
        {
            if (MIRBOOKING_OUTPUT_FORMAT_META[i].output_format == output_format)
            {
                MIRBOOKING_OUTPUT_FORMAT_META[i].write (mirbooking, output_f);
                break;
            }
        }
    }

#if HAVE_MPI
    MPI_Finalize ();
#endif

    return EXIT_SUCCESS;
}
