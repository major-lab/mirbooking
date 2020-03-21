#include "config.h"

#include <errno.h>
#include <fcntl.h>
#include <gio/gio.h>
#include <gio/gunixoutputstream.h>
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

#if !GLIB_CHECK_VERSION(2,57,0)
G_DEFINE_AUTOPTR_CLEANUP_FUNC (GEnumClass, g_type_class_unref);
#endif

#define MIRBOOKING_DEFAULT_MAX_ITERATIONS          100
#define MIRBOOKING_DEFAULT_CUTOFF                  100 // pM
#define MIRBOOKING_DEFAULT_BOUND_FRACTION_CUTOFF   0.0
#define MIRBOOKING_DEFAULT_SPONGED_FRACTION_CUTOFF 0.0

static gchar                        **targets_files             = {NULL};
static gchar                        **mirnas_files              = {NULL};
static gchar                         *seed_scores_file          = MIRBOOKING_DEFAULT_SEED_SCORES_FILE;
static gchar                         *supplementary_scores_file = NULL;
static gchar                         *accessibility_scores_file = NULL;
static gchar                         *input_file                = NULL;
static gchar                         *output_file               = NULL;
static MirbookingBrokerOutputFormat   output_format             = MIRBOOKING_BROKER_OUTPUT_FORMAT_TSV;
static MirbookingBrokerSparseSolver   sparse_solver;
static guint64                        max_iterations            = MIRBOOKING_DEFAULT_MAX_ITERATIONS;
static gsize                          prime5_footprint          = MIRBOOKING_BROKER_DEFAULT_5PRIME_FOOTPRINT;
static gsize                          prime3_footprint          = MIRBOOKING_BROKER_DEFAULT_3PRIME_FOOTPRINT;
static gdouble                        cutoff                    = MIRBOOKING_DEFAULT_CUTOFF;
static gdouble                        bound_fraction_cutoff     = MIRBOOKING_DEFAULT_BOUND_FRACTION_CUTOFF;
static gdouble                        sponged_fraction_cutoff   = MIRBOOKING_DEFAULT_SPONGED_FRACTION_CUTOFF;
static gboolean                       version                   = FALSE;
static gchar                         *blacklist                 = NULL;

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
        if (supplementary_scores_file == NULL)
            supplementary_scores_file = MIRBOOKING_DEFAULT_WEE_ET_AL_2012_SUPPLEMENTARY_SCORES_FILE;
    }
    else if (g_strcmp0 (value, "yan-et-al-2018") == 0)
    {
        supplementary_model = MIRBOOKING_DEFAULT_SCORE_TABLE_SUPPLEMENTARY_MODEL_YAN_ET_AL_2018;
        if (supplementary_scores_file == NULL)
            supplementary_scores_file = MIRBOOKING_DEFAULT_YAN_ET_AL_2018_SUPPLEMENTARY_SCORES_FILE;
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
    g_autoptr (GEnumClass) output_format_class = g_type_class_ref (MIRBOOKING_BROKER_OUTPUT_FORMAT_ENUM);
    GEnumValue *eval = g_enum_get_value_by_nick (output_format_class,
                                                 value);

    if (eval == NULL)
    {
        g_set_error (error,
                     G_OPTION_ERROR,
                     G_OPTION_ERROR_BAD_VALUE,
                     "No such output format %s.", value);
        return FALSE;
    }

    output_format = eval->value;

    return TRUE;
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
    {"targets",                 0, 0, G_OPTION_ARG_FILENAME_ARRAY, &targets_files,             "Targets FASTA files",                                                                              NULL},
    {"mirnas",                  0, 0, G_OPTION_ARG_FILENAME_ARRAY, &mirnas_files,              "miRNA FASTA files",                                                                                NULL},
    {"seed-scores",             0, 0, G_OPTION_ARG_FILENAME,       &seed_scores_file,          "Precomputed seed::MRE binding free energy duplex table",                                           "FILE"},
    {"supplementary-model",     0, 0, G_OPTION_ARG_CALLBACK,       &set_supplementary_model,   "Supplementary bindings model to use",                                                              "none"},
    {"supplementary-scores",    0, 0, G_OPTION_ARG_FILENAME,       &supplementary_scores_file, "Precomputed supplementary::MRE binding free energy duplex table",                                  "FILE"},
    {"accessibility-scores",    0, 0, G_OPTION_ARG_FILENAME,       &accessibility_scores_file, "Accessibility scores as a variable columns (accession, positions...) TSV file",                    "FILE"},
    {"input",                   0, 0, G_OPTION_ARG_FILENAME,       &input_file,                "MiRNA and targets quantities as a two columns (accession, quantity) TSV file (defaults to stdin)", "FILE"},
    {"output",                  0, 0, G_OPTION_ARG_FILENAME,       &output_file,               "Output destination file (defaults to stdout)",                                                     "FILE"},
    {"output-format",           0, 0, G_OPTION_ARG_CALLBACK,       &set_output_format,         "Output format (i.e. 'tsv', 'tsv-detailed', 'gff3', 'wig')",                                        "tsv"},
    {"sparse-solver",           0, 0, G_OPTION_ARG_CALLBACK,       &set_sparse_solver,         "Sparse solver implementation to use",                                                              "best-available"},
    {"max-iterations",          0, 0, G_OPTION_ARG_INT,            &max_iterations,            "Maximum number of iterations",                                                                     G_STRINGIFY (MIRBOOKING_DEFAULT_MAX_ITERATIONS)},
    {"5prime-footprint",        0, 0, G_OPTION_ARG_INT,            &prime5_footprint,          "Footprint in the MRE's 5' direction",                                                              G_STRINGIFY (MIRBOOKING_BROKER_DEFAULT_5PRIME_FOOTPRINT)},
    {"3prime-footprint",        0, 0, G_OPTION_ARG_INT,            &prime3_footprint,          "Footprint in the MRE's 3' direction",                                                              G_STRINGIFY (MIRBOOKING_BROKER_DEFAULT_3PRIME_FOOTPRINT)},
    {"cutoff",                  0, 0, G_OPTION_ARG_DOUBLE,         &cutoff,                    "Cutoff on the duplex concentration",                                                               G_STRINGIFY (MIRBOOKING_DEFAULT_CUTOFF)},
    {"bound-fraction-cutoff",   0, 0, G_OPTION_ARG_DOUBLE,         &bound_fraction_cutoff,     "Cutoff on the bound fraction",                                                                     G_STRINGIFY (MIRBOOKING_DEFAULT_BOUND_FRACTION_CUTOFF)},
    {"sponged-fraction-cutoff", 0, 0, G_OPTION_ARG_DOUBLE,         &sponged_fraction_cutoff,   "Cutoff on the sponged fraction",                                                                   G_STRINGIFY (MIRBOOKING_DEFAULT_SPONGED_FRACTION_CUTOFF)},
    {"blacklist",               0, 0, G_OPTION_ARG_FILENAME,       &blacklist,                 "Interaction blacklist as a three columns (accession, position, quantity) TSV file",                NULL},
    {"version",                 0, 0, G_OPTION_ARG_NONE,           &version,                   "Show version and exit",                                                                            NULL},
    {0}
};

static MirbookingFastaFormat
detect_fasta_format (const gchar *path)
{
    g_autofree const gchar *path_basename = g_path_get_basename (path);

    if (g_str_has_prefix (path_basename, "gencode."))
    {
        return MIRBOOKING_FASTA_FORMAT_GENCODE;
    }
    else if (g_str_has_prefix (path_basename, "mature"))
    {
        return MIRBOOKING_FASTA_FORMAT_MIRBASE;
    }
    else if (g_str_has_prefix (path_basename, "GCF_") || g_str_has_prefix (path_basename, "GCA_"))
    {
        return MIRBOOKING_FASTA_FORMAT_NCBI;
    }
    else
    {
        return MIRBOOKING_FASTA_FORMAT_GENERIC;
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


gboolean
read_interaction_blacklist (GInputStream  *is,
                            GHashTable    *blacklist,
                            GHashTable    *sequences_hash,
                            GError       **error)
{
    g_autoptr (GDataInputStream) dis = g_data_input_stream_new (is);

    guint lineno = 0;
    do
    {
        g_autofree gchar *line = g_data_input_stream_read_line (dis,
                                                                NULL,
                                                                NULL,
                                                                error);

        if (!line)
        {
            return *error == NULL;
        }

        ++lineno;

        if (lineno == 1)
        {
            if (!g_str_equal (line, "target_accession\tposition\tmirna_accession"))
            {
                g_set_error (error,
                             G_FILE_ERROR,
                             G_FILE_ERROR_FAILED,
                             "Unexpected TSV header '%s' from blacklist at line %u.",
                             line,
                             lineno);
                return FALSE;
            }
            else
            {
                continue;
            }
        }

        gchar target_accession[64];
        gsize position = 0;
        gchar mirna_accession[64];
        MirbookingScore score = {0};
        if (sscanf (line, "%s\t%lu\t%s", target_accession, &position, mirna_accession) != 3)
        {
            g_set_error (error,
                         G_FILE_ERROR,
                         G_FILE_ERROR_FAILED,
                         "Could not parse row '%s' from blacklist at line %u.",
                         line,
                         lineno);
            return FALSE;
        }

        MirbookingTarget *target = g_hash_table_lookup (sequences_hash,
                                                        target_accession);
        MirbookingMirna *mirna = g_hash_table_lookup (sequences_hash,
                                                      mirna_accession);

        if (target == NULL)
        {
            g_set_error (error,
                         G_FILE_ERROR,
                         G_FILE_ERROR_FAILED,
                         "Unknown target '%s' in blacklist at line %u.",
                         target_accession,
                         lineno);
            return FALSE;
        }

        if (position >= mirbooking_sequence_get_sequence_length (MIRBOOKING_SEQUENCE(target)))
        {
            g_set_error (error,
                         G_FILE_ERROR,
                         G_FILE_ERROR_FAILED,
                         "Target '%s' does not have a site at %lu since it is only %lu nucleotides long in blacklist at line %u.",
                         target_accession,
                         position,
                         mirbooking_sequence_get_sequence_length (MIRBOOKING_SEQUENCE(target)),
                         lineno);
            return FALSE;
        }

        if (mirna == NULL)
        {
            g_set_error (error,
                         G_FILE_ERROR,
                         G_FILE_ERROR_FAILED,
                         "Unknown miRNA '%s' in blacklist at line %u.",
                         mirna_accession,
                         lineno);
            return FALSE;
        }

        MirbookingOccupant *occupant = g_new (MirbookingOccupant, 1);
        mirbooking_occupant_init (occupant,
                                  target,
                                  position - 1,
                                  mirna,
                                  score);

        g_hash_table_insert (blacklist,
                             occupant,
                             NULL);
    }
    while (TRUE);
}

typedef struct _MixedFilterUserData
{
    MirbookingDefaultScoreTableCutoffFilterUserData    cutoff_filter_ud;
    MirbookingDefaultScoreTableBlacklistFilterUserdata blacklist_filter_ud;
} MixedFilterUserData;

static gboolean
mixed_filter (MirbookingDefaultScoreTable *score_table,
              MirbookingMirna             *mirna,
              MirbookingTarget            *target,
              gssize                       position,
              gpointer                     user_data)
{
    MixedFilterUserData *mixed_filter_ud = user_data;
    gboolean filter_out = TRUE;

    filter_out &= mirbooking_default_score_table_blacklist_filter (score_table,
                                                                   mirna,
                                                                   target,
                                                                   position,
                                                                   &mixed_filter_ud->blacklist_filter_ud);

    filter_out &= mirbooking_default_score_table_cutoff_filter (score_table,
                                                                mirna,
                                                                target,
                                                                position,
                                                                &mixed_filter_ud->cutoff_filter_ud);

    return filter_out;
}

static void
free_mixed_filter_user_data (gpointer ud)
{
    MixedFilterUserData *mixed_filter_ud = ud;
    g_object_unref (mixed_filter_ud->cutoff_filter_ud.broker);
    g_hash_table_unref (mixed_filter_ud->blacklist_filter_ud.blacklist);
}

static void
mirbooking_occupant_free (MirbookingOccupant *occupant)
{
    mirbooking_occupant_clear (occupant);
    g_free (occupant);
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

    sparse_solver = MIRBOOKING_BROKER_DEFAULT_SPARSE_SOLVER;

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

    if (version)
    {
        g_print ("%s\n", MIRBOOKING_VERSION);
        return EXIT_SUCCESS;
    }

    if (!mirbooking_broker_sparse_solver_is_available (sparse_solver))
    {
        g_printerr ("No proper sparse solver could be determined; use the '--sparse-solver' argument explicitly.\n");
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
        if (supplementary_model == MIRBOOKING_DEFAULT_SCORE_TABLE_SUPPLEMENTARY_MODEL_NONE)
        {
            g_printerr ("The '--supplementary-model' argument must be set in order to provide a supplementary score file.\n");
            return EXIT_FAILURE;
        }

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

    MirbookingDefaultScoreTableCutoffFilterUserData cutoff_filter_ud =
    {
        .broker                  = g_object_ref (mirbooking),
        .cutoff                  = cutoff,
        .bound_fraction_cutoff   = bound_fraction_cutoff,
        .sponged_fraction_cutoff = sponged_fraction_cutoff
    };


    GHashTable *blacklist_map = g_hash_table_new_full ((GHashFunc)      mirbooking_occupant_hash,
                                                       (GEqualFunc)     mirbooking_occupant_equal,
                                                       (GDestroyNotify) mirbooking_occupant_free,
                                                       NULL);

    MirbookingDefaultScoreTableBlacklistFilterUserdata blacklist_filter_ud =
    {
        .blacklist = blacklist_map
    };

    MixedFilterUserData user_data =
    {
        .blacklist_filter_ud = blacklist_filter_ud,
        .cutoff_filter_ud    = cutoff_filter_ud
    };

    mirbooking_default_score_table_set_filter (score_table,
                                               mixed_filter,
                                               &user_data,
                                               free_mixed_filter_user_data);

    mirbooking_broker_set_score_table (mirbooking,
                                       MIRBOOKING_SCORE_TABLE (score_table));

    g_autoptr (FILE) input_f = NULL;

    if ((input_f = input_file == NULL ? stdin : g_fopen (input_file, "r")) == NULL)
    {
        g_printerr ("Could not open the quantities file '%s': %s.\n", input_file, g_strerror (errno));
        return EXIT_FAILURE;
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
            g_autoptr (GMappedFile) seq_map = g_mapped_file_new (*cur_sequences_file, FALSE, &error);

            if (seq_map == NULL)
            {
                g_printerr ("%s (%s, %u).\n", error->message, g_quark_to_string (error->domain), error->code);
                return EXIT_FAILURE;
            }

            g_autoptr (MirbookingSequenceIter) iter;
            mirbooking_read_sequences_from_mapped_file (seq_map,
                                                        detect_fasta_format (*cur_sequences_file),
                                                        MIRBOOKING_TYPE_TARGET,
                                                        &iter);

            while (mirbooking_sequence_iter_next (iter, &error))
            {
                MirbookingSequence *seq = mirbooking_sequence_iter_get_sequence (iter);
                g_hash_table_insert (sequences_hash,
                                     g_strdup (mirbooking_sequence_get_accession (seq)),
                                     g_object_ref (seq));
            }

            if (error != NULL)
            {
                g_printerr ("%s (%s, %u).\n", error->message, g_quark_to_string (error->domain), error->code);
                return EXIT_FAILURE;
            }
        }
    }

    if (mirnas_files != NULL)
    {
        for (cur_sequences_file = mirnas_files; *cur_sequences_file != NULL; cur_sequences_file++)
        {
            g_autoptr (GMappedFile) seq_map = g_mapped_file_new (*cur_sequences_file, FALSE, &error);

            if (seq_map == NULL)
            {
                g_printerr ("%s (%s, %u).\n", error->message, g_quark_to_string (error->domain), error->code);
                return EXIT_FAILURE;
            }

            g_autoptr (MirbookingSequenceIter) iter;
            mirbooking_read_sequences_from_mapped_file (seq_map,
                                                        detect_fasta_format (*cur_sequences_file),
                                                        MIRBOOKING_TYPE_MIRNA,
                                                        &iter);

            while (mirbooking_sequence_iter_next (iter, &error))
            {
                MirbookingSequence *seq = mirbooking_sequence_iter_get_sequence (iter);
                g_hash_table_insert (sequences_hash,
                                     g_strdup (mirbooking_sequence_get_accession (seq)),
                                     g_object_ref (seq));
            }

            if (error != NULL)
            {
                g_printerr ("%s (%s, %u).\n", error->message, g_quark_to_string (error->domain), error->code);
                return EXIT_FAILURE;
            }
        }
    }

    if (blacklist != NULL)
    {
        g_autoptr (GFile) bf = g_file_new_for_path (blacklist);
        if (!read_interaction_blacklist (G_INPUT_STREAM (g_file_read (bf, NULL, NULL)),
                                         blacklist_map,
                                         sequences_hash,
                                         &error))
        {
            g_printerr ("%s (%s, %u).\n", error->message, g_quark_to_string (error->domain), error->code);
            return EXIT_FAILURE;
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
        g_autoptr (GOutputStream) out = NULL;

        if (output_file == NULL)
        {
            out = g_unix_output_stream_new (fileno (stdout), FALSE);
        }
        else
        {
            g_autoptr (GFile) outputf = g_file_new_for_path (output_file);

            out = G_OUTPUT_STREAM (g_file_replace (outputf,
                                                   NULL,
                                                   FALSE,
                                                   G_FILE_CREATE_NONE,
                                                   NULL,
                                                   &error));
            if (out == NULL)
            {
                g_printerr ("%s (%s, %u)\n", error->message, g_quark_to_string (error->domain), error->code);
                return EXIT_FAILURE;
            }
        }

        if (!mirbooking_broker_write_output_to_stream (mirbooking,
                                                       out,
                                                       output_format,
                                                       &error))
        {
            g_printerr ("%s (%s, %u)\n", error->message, g_quark_to_string (error->domain), error->code);
            return EXIT_FAILURE;
        }
    }

#if HAVE_MPI
    MPI_Finalize ();
#endif

    return EXIT_SUCCESS;
}
