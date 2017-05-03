#include <glib.h>
#include <glib/gprintf.h>
#include <fcntl.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <glib/gstdio.h>
#include <gio/gio.h>

#include <mirbooking.h>

G_DEFINE_AUTOPTR_CLEANUP_FUNC (FILE, fclose)

static gchar *mirnas_file;
static gchar *targets_file;
static gchar *score_table_file;
static gchar *quantities_file;
static gchar *output_file;
static gdouble threshold;
static gdouble log_base;

static GOptionEntry MIRBOOKING_OPTION_ENTRIES[] =
{
    {"mirnas",      0, 0, G_OPTION_ARG_FILENAME, &mirnas_file,      "", "FILE"},
    {"targets",     0, 0, G_OPTION_ARG_FILENAME, &targets_file,     "", "FILE"},
    {"score-table", 0, 0, G_OPTION_ARG_FILENAME, &score_table_file, "", "FILE"},
    {"quantities",  0, 0, G_OPTION_ARG_FILENAME, &quantities_file,  "", "FILE"},
    {"output",      0, 0, G_OPTION_ARG_FILENAME, &output_file,      "", "FILE"},
    {"threshold",   0, 0, G_OPTION_ARG_DOUBLE,   &threshold,        "", G_STRINGIFY (MIRBOOKING_DEFAULT_THRESHOLD)},
    {"log-base",    0, 0, G_OPTION_ARG_DOUBLE,   &log_base,         "", G_STRINGIFY (MIRBOOKING_DEFAULT_LOG_BASE)},
    {NULL}
};

static void
read_sequences_from_fasta (FILE        *file,
                           GMappedFile *mapped_file,
                           guint        accession_column,
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

            gint i;
            for (i = 0; i < accession_column; i++)
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

            mirbooking_sequence_set_sequence (sequence,
                                              seq,
                                              seq_len);

            g_hash_table_insert (sequences_hash, g_strdup (accession), sequence);
        }
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

    if (score_table_file == NULL)
    {
        g_printerr ("The '--score-table' argument is required.\n");
        return EXIT_FAILURE;
    }

    g_autoptr (GMappedFile) score_map = g_mapped_file_new (score_table_file, FALSE, NULL);
    g_return_val_if_fail (score_map != NULL, EXIT_FAILURE);

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

    // precondition mirnas and targets
    read_sequences_from_fasta (mirnas_f, mirnas_map, 1, sequences_hash);
    read_sequences_from_fasta (targets_f, targets_map, 0, sequences_hash);

    gchar line[1024];
    while (fgets (line, sizeof (line), quantities_f))
    {
        gchar *accession = strtok (line, "\t");
        gdouble quantity = g_strtod (strtok (NULL, "\t"), NULL);

        MirbookingSequence *sequence = g_hash_table_lookup (sequences_hash, accession);

        if (sequence == NULL)
        {
            g_printerr ("Unknown sequence with accession '%s'.\n", accession);
        }

        if (MIRBOOKING_IS_MIRNA (sequence))
        {
            mirbooking_set_mirna_quantity (mirbooking,
                                           MIRBOOKING_MIRNA (sequence),
                                           quantity);
        }

        if (MIRBOOKING_IS_TARGET (sequence))
        {
            mirbooking_set_target_quantity (mirbooking,
                                            MIRBOOKING_TARGET (sequence),
                                            quantity);
        }
    }

    g_hash_table_unref (sequences_hash);

    if (!mirbooking_run (mirbooking, &error))
    {
        g_printerr ("%s (%s, %u)\n", error->message, g_quark_to_string (error->domain), error->code);
        return EXIT_FAILURE;
    }

    GSList *targets;
    GSList *target_sites;
    GSList *mirnas;

    targets = mirbooking_get_targets (mirbooking);

    g_fprintf (output_f, "Target Accession\tMirna Accession\tSite\tQuantity\n");

    for (; targets != NULL; targets = targets->next)
    {
        target_sites = mirbooking_get_target_sites (mirbooking, targets->data);

        for (; target_sites != NULL; target_sites = target_sites->next)
        {
            mirnas = mirbooking_target_site_get_mirnas (target_sites->data);

            for (; mirnas != NULL; mirnas = mirnas->next)
            {
                g_fprintf (output_f, "%s\t%s\t%ld\t%f\n", mirbooking_sequence_get_accession (targets->data),
                                                          mirbooking_sequence_get_accession (mirnas->data),
                                                          mirbooking_target_site_get_offset (target_sites->data),
                                                          mirbooking_target_site_get_mirna_quantity (target_sites->data, mirnas->data));
            }
        }
    }

    return EXIT_SUCCESS;
}
