#include <mirbooking.h>

static void
test_sequence_iter ()
{
    g_autoptr (GError) error = NULL;
    g_autoptr (GMappedFile) mf = g_mapped_file_new (g_test_get_filename (G_TEST_DIST, "data", "GCF_000001405.39_GRCh38.p13_rna.fna", NULL),
                                                    FALSE,
                                                    &error);
    g_autoptr (MirbookingSequenceIter) iter = NULL;
    mirbooking_read_sequences_from_mapped_file (mf,
                                                MIRBOOKING_FASTA_FORMAT_NCBI,
                                                MIRBOOKING_TYPE_TARGET,
                                                &iter);

    while (mirbooking_sequence_iter_next (iter, &error))
    {
        MirbookingSequence *sequence = mirbooking_sequence_iter_get_sequence (iter);
        g_assert_nonnull (sequence);
        g_assert_cmpstr (mirbooking_sequence_get_accession (sequence), ==, "NM_000014.5");
        g_assert_cmpstr (mirbooking_sequence_get_gene_name (sequence), ==, "A2M");
        g_assert_cmpstr (mirbooking_sequence_get_name (sequence), ==, "A2M-201");
    }

    g_assert_null (error);
}

gint
main (gint argc, gchar **argv)
{
    g_test_init (&argc, &argv, NULL);

    g_test_add_func ("/sequence-iter", test_sequence_iter);

    return g_test_run ();
}
