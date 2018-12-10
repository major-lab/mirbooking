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

#define MIRBOOKING_DEFAULT_TOLERANCE      1e-8
#define MIRBOOKING_DEFAULT_MAX_ITERATIONS 100

typedef enum _MirbookingOutputFormat
{
    MIRBOOKING_OUTPUT_FORMAT_TSV,
    MIRBOOKING_OUTPUT_FORMAT_GFF3
} MirbookingOutputFormat;

static gchar                        **targets_files             = {NULL};
static gchar                        **mirnas_files              = {NULL};
static gchar                         *seed_scores_file          = NULL;
static gchar                         *supplementary_scores_file = NULL;
static gchar                         *accessibility_scores_file = NULL;
static gchar                         *input_file                = NULL;
static gchar                         *output_file               = NULL;
static gboolean                       output_format             = MIRBOOKING_OUTPUT_FORMAT_TSV;
static MirbookingBrokerSparseSolver   sparse_solver             = MIRBOOKING_BROKER_DEFAULT_SPARSE_SOLVER;
static gdouble                        tolerance                 = MIRBOOKING_DEFAULT_TOLERANCE;
static guint64                        max_iterations            = MIRBOOKING_DEFAULT_MAX_ITERATIONS;
static gsize                          prime5_footprint          = MIRBOOKING_BROKER_DEFAULT_5PRIME_FOOTPRINT;
static gsize                          prime3_footprint          = MIRBOOKING_BROKER_DEFAULT_3PRIME_FOOTPRINT;
static gboolean                       verbose                   = FALSE;
static gboolean help = FALSE;

static gboolean
set_output_format (const gchar   *key,
                   const gchar   *value,
                   gpointer       data,
                   GError      **error)
{
    if (g_strcmp0 (value, "tsv") == 0)
    {
        output_format = MIRBOOKING_OUTPUT_FORMAT_TSV;
    }
    else if (g_strcmp0 (value, "gff3") == 0)
    {
        output_format = MIRBOOKING_OUTPUT_FORMAT_GFF3;
    }
    else
    {
        return FALSE;
    }

    return TRUE;
}

static gboolean
set_sparse_solver(const gchar   *key,
                  const gchar   *value,
                  gpointer       data,
                  GError      **error)
{
    GEnumClass *sparse_solver_class = g_type_class_ref (MIRBOOKING_BROKER_SPARSE_SOLVER_ENUM);
    GEnumValue *eval = g_enum_get_value_by_nick (sparse_solver_class,
                                                 value);

    if (eval == NULL)
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
    {"seed-scores",          0, 0, G_OPTION_ARG_FILENAME,       &seed_scores_file,          "Precomputed seed::MRE Gibbs free energy duplex table",                                            "FILE"},
    {"supplementary-scores", 0, 0, G_OPTION_ARG_FILENAME,       &supplementary_scores_file, "Precomputed supplementary::MRE Gibbs free energy duplex table",                                   "FILE"},
    {"accessibility-scores", 0, 0, G_OPTION_ARG_FILENAME,       &accessibility_scores_file, "Accessibility scores as a variable columns (accession, positions...) TSV file",                   "FILE"},
    {"input",                0, 0, G_OPTION_ARG_FILENAME,       &input_file,                "MiRNA and targets quantities as a two-column (accession, quantity) TSV file (defaults to stdin)", "FILE"},
    {"output",               0, 0, G_OPTION_ARG_FILENAME,       &output_file,               "Output destination file (defaults to stdout)",                                                    "FILE"},
    {"output-format",        0, 0, G_OPTION_ARG_CALLBACK,       &set_output_format,         "Output format (e.g. 'tsv', 'gff3')",                                                              "tsv"},
    {"sparse-solver",        0, 0, G_OPTION_ARG_CALLBACK,       &set_sparse_solver,         "Sparse solver implementation to use",                                                             G_STRINGIFY (MIRBOOKING_BROKER_DEFAULT_SPARSE_SOLVER)},
    {"tolerance",            0, 0, G_OPTION_ARG_DOUBLE,         &tolerance,                 "Absolute tolerance for the system norm to declare convergence",                                   G_STRINGIFY (MIRBOOKING_DEFAULT_TOLERANCE)},
    {"max-iterations",       0, 0, G_OPTION_ARG_INT,            &max_iterations,            "Maximum number of iterations",                                                                    G_STRINGIFY (MIRBOOKING_DEFAULT_MAX_ITERATIONS)},
    {"5prime-footprint",     0, 0, G_OPTION_ARG_INT,            &prime5_footprint,          "Footprint in the MRE's 5' direction",                                                             G_STRINGIFY (MIRBOOKING_BROKER_DEFAULT_5PRIME_FOOTPRINT)},
    {"3prime-footprint",     0, 0, G_OPTION_ARG_INT,            &prime3_footprint,          "Footprint in the MRE's 3' direction",                                                             G_STRINGIFY (MIRBOOKING_BROKER_DEFAULT_3PRIME_FOOTPRINT)},
    {"verbose",              0, 0, G_OPTION_ARG_NONE,           &verbose,                   "Turn on verbose output",                                                                          NULL},
    {"help",                 0, 0, G_OPTION_ARG_NONE,           &help,                      "Show this help message",                                                                          NULL},
    {NULL}
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
                gint i;
                for (i = 0; i < 5; i++)
                {
                    name = strtok (NULL, "|");
                }
            }
            else if (fasta_format == FASTA_FORMAT_NCBI)
            {
                accession = strtok (line + 1, " ");
                name = NULL;

                /* pick last opening brace */
                gchar *name_p;
                while (strtok (NULL, "(") && (name_p = strtok (NULL, ")")))
                {
                    name = name_p;
                }
            }
            else
            {
                accession = strtok (line + 1, " ");
                name = strtok (NULL, "\n");
            }

            seq = g_mapped_file_get_contents (mapped_file) + ftell (file);

            gsize remaining = g_mapped_file_get_length (mapped_file) - ftell (file);
            gchar *next_seq = memchr (seq, '>', remaining);
            gsize seq_len = next_seq == NULL ? remaining - 1 : next_seq - seq - 1;

            MirbookingSequence *sequence;

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

            mirbooking_sequence_set_raw_sequence (sequence,
                                                  seq,
                                                  seq_len);

            g_hash_table_insert (sequences_hash,
                                 g_strdup (accession),
                                 sequence);
        }
    }
}

void
read_sequence_accessibility (GInputStream     *in,
                             GHashTable       *sequences_hash,
                             MirbookingBroker *broker)
{
    g_autoptr (GDataInputStream) data_in = g_data_input_stream_new (in);

    while (TRUE)
    {
        g_autofree gchar *accession = g_data_input_stream_read_upto (data_in,
                                                                     "\t",
                                                                     1,
                                                                     NULL,
                                                                     NULL,
                                                                     NULL);

        if (accession == NULL)
        {
            return; /* done */
        }

        MirbookingSequence *sequence = g_hash_table_lookup (sequences_hash,
                                                            accession);

        if (sequence && mirbooking_broker_get_sequence_quantity (broker, sequence) > 0)
        {
            gsize sequence_len = mirbooking_sequence_get_sequence_length (sequence);
            gsize position;
            for (position = 0; position < sequence_len; position++)
            {
                g_return_if_fail (g_data_input_stream_read_byte (data_in, NULL, NULL) == '\t');

                g_autofree gchar *score_str = g_data_input_stream_read_upto (data_in, "\t", 1, NULL, NULL, NULL);

                gfloat score;
                if (*score_str == '\0')
                {
                    score = INFINITY;
                }
                else
                {
                    sscanf (score_str, "%f", &score);
                }

                // FIXME: have prefix instead of suffix positions
                mirbooking_target_set_accessibility_score (MIRBOOKING_TARGET (sequence),
                                                           (position - 7) % sequence_len,
                                                           score);
            }
        }

        // consume the rest of the line
        g_autofree gchar *remaining = g_data_input_stream_read_line (data_in, NULL, NULL, NULL);
    }
}

#define COALESCE(x,d) (x == NULL ? (d) : (x))

static void
write_output_to_tsv (MirbookingBroker *mirbooking,
                     FILE             *output_f)
{
    g_fprintf (output_f, "target_accession\t"
                         "target_name\t"
                         "target_quantity\t"
                         "position\t"
                         "region\t"
                         "occupancy\t"
                         "mirna_accession\t"
                         "mirna_name\t"
                         "mirna_quantity\t"
                         "score\t"
                         "quantity\n");

    GArray *target_sites = mirbooking_broker_get_target_sites (mirbooking);

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
            g_fprintf (output_f, "%s\t%s\t%e\t%lu\t%s\t%s\t%e\t%e\t%e\n",
                       mirbooking_sequence_get_accession (MIRBOOKING_SEQUENCE (target_site->target)),
                       COALESCE (mirbooking_sequence_get_name (MIRBOOKING_SEQUENCE (target_site->target)), "N/A"),
                       target_quantity,
                       target_site->position + 1, // 1-based
                       mirbooking_sequence_get_accession (MIRBOOKING_SEQUENCE (occupant->mirna)),
                       COALESCE (mirbooking_sequence_get_name (MIRBOOKING_SEQUENCE (occupant->mirna)), "N/A"),
                       mirbooking_broker_get_sequence_quantity (mirbooking, MIRBOOKING_SEQUENCE (occupant->mirna)),
                       occupant->enzymatic_score,
                       mirbooking_broker_get_occupant_quantity (mirbooking, occupant));
        }
    }
}

static void
write_output_to_gff3 (MirbookingBroker *mirbooking, FILE *output_f)
{
    g_fprintf (output_f, "##gff-version 3\n");

    GArray *target_sites = mirbooking_broker_get_target_sites (mirbooking);

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

int
main (gint argc, gchar **argv)
{
#if HAVE_MPI
    gint provided;
    MPI_Init_thread (&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    g_return_val_if_fail (provided != MPI_THREAD_FUNNELED, EXIT_FAILURE);
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

    g_autoptr (MirbookingBroker) mirbooking = mirbooking_broker_new ();

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
                                                                                              MIRBOOKING_DEFAULT_SCORE_TABLE_DEFAULT_SUPPLEMENTARY_MODEL,
                                                                                              supplementary_scores_map_bytes);
    mirbooking_broker_set_score_table (mirbooking,
                                       MIRBOOKING_SCORE_TABLE (g_object_ref (score_table)));

    g_autoptr (FILE) input_f = NULL;
    g_autoptr (FILE) output_f = NULL;

    if ((input_f = input_file == NULL ? stdin : g_fopen (input_file, "r")) == NULL)
    {
        g_printerr ("Could not open the quantities file '%s': %s.\n", input_file, g_strerror (errno));
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

    // precondition all sequences
    gchar **cur_sequences_file = NULL;
    for (cur_sequences_file = targets_files; *cur_sequences_file != NULL; cur_sequences_file++)
    {
        g_autoptr (FILE) seq_f = g_fopen (*cur_sequences_file, "r");

        if (seq_f == NULL)
        {
            g_printerr ("Could not open the sequences file '%s': %s.\n", *cur_sequences_file, g_strerror (errno));
            return EXIT_FAILURE;
        }

        g_autoptr (GMappedFile) seq_map = g_mapped_file_new_from_fd (fileno (seq_f), FALSE, &error);

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

    for (cur_sequences_file = mirnas_files; *cur_sequences_file != NULL; cur_sequences_file++)
    {
        g_autoptr (FILE) seq_f = g_fopen (*cur_sequences_file, "r");

        if (seq_f == NULL)
        {
            g_printerr ("Could not open the sequences file '%s': %s.\n", *cur_sequences_file, g_strerror (errno));
            return EXIT_FAILURE;
        }

        g_autoptr (GMappedFile) seq_map = g_mapped_file_new_from_fd (fileno (seq_f), FALSE, &error);

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
            continue;
            return EXIT_FAILURE;
        }

        mirbooking_broker_set_sequence_quantity (mirbooking,
                                                 sequence,
                                                 quantity);
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
        read_sequence_accessibility (G_INPUT_STREAM (accessibility_scores_in),
                                     sequences_hash,
                                     mirbooking);

        g_debug ("Done reading accessibility scores in %lums", 1000 * (g_get_monotonic_time () - accessibility_read_begin) / G_USEC_PER_SEC);
    }

    guint64 iteration = 0;
    guint64 iteration_begin = g_get_monotonic_time ();
    gdouble norm;
    do
    {
        iteration++;

        if (iteration > max_iterations)
        {
            g_print ("Did not converge after %lu iterations.", iteration);
            return EXIT_FAILURE;
        }

        iteration_begin = g_get_monotonic_time ();

        if (!mirbooking_broker_evaluate (mirbooking,
                                         &norm,
                                         &error))
        {
            g_printerr ("%s (%s, %u)\n", error->message, g_quark_to_string (error->domain), error->code);
            return EXIT_FAILURE;
        }

        g_debug ("iteration: %lu norm: %.2e time: %lums",
                 iteration,
                 norm,
                 1000 * (g_get_monotonic_time () - iteration_begin) / G_USEC_PER_SEC);

        if (!isfinite (norm))
        {
            g_printerr ("Failed to converge.\n");
            return EXIT_FAILURE;
        }

        if (norm <= tolerance)
        {
            break;
        }

        if (!mirbooking_broker_step (mirbooking,
                                     MIRBOOKING_BROKER_STEP_MODE_SOLVE_STEADY_STATE,
                                     1,
                                     &error))
        {
            g_printerr ("%s (%s, %u)\n", error->message, g_quark_to_string (error->domain), error->code);
            return EXIT_FAILURE;
        }
    }
    while (TRUE);

    switch (output_format)
    {
        case MIRBOOKING_OUTPUT_FORMAT_TSV:
            write_output_to_tsv (mirbooking,
                                 output_f);
            break;
        case MIRBOOKING_OUTPUT_FORMAT_GFF3:
            write_output_to_gff3 (mirbooking,
                                  output_f);
            break;
        default:
            return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
