#include <glib.h>
#include <mirbooking.h>
#include <string.h>

static gfloat SCORE_TABLE[16384][16384];

static gchar *TARGET_SEQUENCE =
    "GCACACAGAGCAGCATAAAGCCCAGTTGCTTTGGGAAGTGTTTGGGACCAGATGGATTGT"
    "AGGGAGTAGGGTACAATACAGTCTGTTCTCCTCCAGCTCCTTCTTTCTGCAACATGGGGA"
    "AGAACAAACTCCTTCATCCAAGTCTGGTTCTTCTCCTCTTGGTCCTCCTGCCCACAGACG"
    "CCTCAGTCTCTGGAAAACCGCAGTATATGGTTCTGGTCCCCTCCCTGCTCCACACTGAGA"
    "CCACTGAGAAGGGCTGTGTCCTTCTGAGCTACCTGAATGAGACAGTGACTGTAAGTGCTT"
    "CCTTGGAGTCTGTCAGGGGAAACAGGAGCCTCTTCACTGACCTGGAGGCGGAGAATGACG"
    "TACTCCACTGTGTCGCCTTCGCTGTCCCAAAGTCTTCATCCAATGAGGAGGTAATGTTCC"
    "TCACTGTCCAAGTGAAAGGACCAACCCAAGAATTTAAGAAGCGGACCACAGTGATGGTTA"
    "AGAACGAGGACAGTCTGGTCTTTGTCCAGACAGACAAATCAATCTACAAACCAGGGCAGA"
    "CAGTGAAATTTCGTGTTGTCTCCATGGATGAAAACTTTCACCCCCTGAATGAGTTGATTC"
    "CACTAGTATACATTCAGGATCCCAAAGGAAATCGCATCGCACAATGGCAGAGTTTCCAGT"
    "TAGAGGGTGGCCTCAAGCAATTTTCTTTTCCCCTCTCATCAGAGCCCTTCCAGGGCTCCT"
    "ACAAGGTGGTGGTACAGAAGAAATCAGGTGGAAGGACAGAGCACCCTTTCACCGTGGAGG"
    "AATTTGTTCTTCCCAAGTTTGAAGTACAAGTAACAGTGCCAAAGATAATCACCATCTTGG"
    "AAGAAGAGATGAATGTATCAGTGTGTGGCCTATACACATATGGGAAGCCTGTCCCTGGAC"
    "ATGTGACTGTGAGCATTTGCAGAAAGTATAGTGACGCTTCCGACTGCCACGGTGAAGATT"
    "CACAGGCTTTCTGTGAGAAATTCAGTGGACAGCTAAACAGCCATGGCTGCTTCTATCAGC"
    "AAGTAAAAACCAAGGTCTTCCAGCTGAAGAGGAAGGAGTATGAAATGAAACTTCACACTG"
    "AGGCCCAGATCCAAGAAGAAGGAACAGTGGTGGAATTGACTGGAAGGCAGTCCAGTGAAA"
    "TCACAAGAACCATAACCAAACTCTCATTTGTGAAAGTGGACTCACACTTTCGACAGGGAA"
    "TTCCCTTCTTTGGGCAGGTGCGCCTAGTAGATGGGAAAGGCGTCCCTATACCAAATAAAG"
    "TCATATTCATCAGAGGAAATGAAGCAAACTATTACTCCAATGCTACCACGGATGAGCATG"
    "GCCTTGTACAGTTCTCTATCAACACCACCAATGTTATGGGTACCTCTCTTACTGTTAGGG"
    "TCAATTACAAGGATCGTAGTCCCTGTTACGGCTACCAGTGGGTGTCAGAAGAACACGAAG"
    "AGGCACATCACACTGCTTATCTTGTGTTCTCCCCAAGCAAGAGCTTTGTCCACCTTGAGC"
    "CCATGTCTCATGAACTACCCTGTGGCCATACTCAGACAGTCCAGGCACATTATATTCTGA"
    "ATGGAGGCACCCTGCTGGGGCTGAAGAAGCTCTCCTTCTATTATCTGATAATGGCAAAGG"
    "GAGGCATTGTCCGAACTGGGACTCATGGACTGCTTGTGAAGCAGGAAGACATGAAGGGCC"
    "ATTTTTCCATCTCAATCCCTGTGAAGTCAGACATTGCTCCTGTCGCTCGGTTGCTCATCT"
    "ATGCTGTTTTACCTACCGGGGACGTGATTGGGGATTCTGCAAAATATGATGTTGAAAATT"
    "GTCTGGCCAACAAGGTGGATTTGAGCTTCAGCCCATCACAAAGTCTCCCAGCCTCACACG"
    "CCCACCTGCGAGTCACAGCGGCTCCTCAGTCCGTCTGCGCCCTCCGTGCTGTGGACCAAA"
    "GCGTGCTGCTCATGAAGCCTGATGCTGAGCTCTCGGCGTCCTCGGTTTACAACCTGCTAC"
    "CAGAAAAGGACCTCACTGGCTTCCCTGGGCCTTTGAATGACCAGGACGATGAAGACTGCA"
    "TCAATCGTCATAATGTCTATATTAATGGAATCACATATACTCCAGTATCAAGTACAAATG"
    "AAAAGGATATGTACAGCTTCCTAGAGGACATGGGCTTAAAGGCATTCACCAACTCAAAGA"
    "TTCGTAAACCCAAAATGTGTCCACAGCTTCAACAGTATGAAATGCATGGACCTGAAGGTC"
    "TACGTGTAGGTTTTTATGAGTCAGATGTAATGGGAAGAGGCCATGCACGCCTGGTGCATG"
    "TTGAAGAGCCTCACACGGAGACCGTACGAAAGTACTTCCCTGAGACATGGATCTGGGATT"
    "TGGTGGTGGTAAACTCAGCAGGTGTGGCTGAGGTAGGAGTAACAGTCCCTGACACCATCA"
    "CCGAGTGGAAGGCAGGGGCCTTCTGCCTGTCTGAAGATGCTGGACTTGGTATCTCTTCCA"
    "CTGCCTCTCTCCGAGCCTTCCAGCCCTTCTTTGTGGAGCTCACAATGCCTTACTCTGTGA"
    "TTCGTGGAGAGGCCTTCACACTCAAGGCCACGGTCCTAAACTACCTTCCCAAATGCATCC"
    "GGGTCAGTGTGCAGCTGGAAGCCTCTCCCGCCTTCCTAGCTGTCCCAGTGGAGAAGGAAC"
    "AAGCGCCTCACTGCATCTGTGCAAACGGGCGGCAAACTGTGTCCTGGGCAGTAACCCCAA"
    "AGTCATTAGGAAATGTGAATTTCACTGTGAGCGCAGAGGCACTAGAGTCTCAAGAGCTGT"
    "GTGGGACTGAGGTGCCTTCAGTTCCTGAACACGGAAGGAAAGACACAGTCATCAAGCCTC"
    "TGTTGGTTGAACCTGAAGGACTAGAGAAGGAAACAACATTCAACTCCCTACTTTGTCCAT"
    "CAGGTGGTGAGGTTTCTGAAGAATTATCCCTGAAACTGCCACCAAATGTGGTAGAAGAAT"
    "CTGCCCGAGCTTCTGTCTCAGTTTTGGGAGACATATTAGGCTCTGCCATGCAAAACACAC"
    "AAAATCTTCTCCAGATGCCCTATGGCTGTGGAGAGCAGAATATGGTCCTCTTTGCTCCTA"
    "ACATCTATGTACTGGATTATCTAAATGAAACACAGCAGCTTACTCCAGAGATCAAGTCCA"
    "AGGCCATTGGCTATCTCAACACTGGTTACCAGAGACAGTTGAACTACAAACACTATGATG"
    "GCTCCTACAGCACCTTTGGGGAGCGATATGGCAGGAACCAGGGCAACACCTGGCTCACAG"
    "CCTTTGTTCTGAAGACTTTTGCCCAAGCTCGAGCCTACATCTTCATCGATGAAGCACACA"
    "TTACCCAAGCCCTCATATGGCTCTCCCAGAGGCAGAAGGACAATGGCTGTTTCAGGAGCT"
    "CTGGGTCACTGCTCAACAATGCCATAAAGGGAGGAGTAGAAGATGAAGTGACCCTCTCCG"
    "CCTATATCACCATCGCCCTTCTGGAGATTCCTCTCACAGTCACTCACCCTGTTGTCCGCA"
    "ATGCCCTGTTTTGCCTGGAGTCAGCCTGGAAGACAGCACAAGAAGGGGACCATGGCAGCC"
    "ATGTATATACCAAAGCACTGCTGGCCTATGCTTTTGCCCTGGCAGGTAACCAGGACAAGA"
    "GGAAGGAAGTACTCAAGTCACTTAATGAGGAAGCTGTGAAGAAAGACAACTCTGTCCATT"
    "GGGAGCGCCCTCAGAAACCCAAGGCACCAGTGGGGCATTTTTACGAACCCCAGGCTCCCT"
    "CTGCTGAGGTGGAGATGACATCCTATGTGCTCCTCGCTTATCTCACGGCCCAGCCAGCCC"
    "CAACCTCGGAGGACCTGACCTCTGCAACCAACATCGTGAAGTGGATCACGAAGCAGCAGA"
    "ATGCCCAGGGCGGTTTCTCCTCCACCCAGGACACAGTGGTGGCTCTCCATGCTCTGTCCA"
    "AATATGGAGCAGCCACATTTACCAGGACTGGGAAGGCTGCACAGGTGACTATCCAGTCTT"
    "CAGGGACATTTTCCAGCAAATTCCAAGTGGACAACAACAACCGCCTGTTACTGCAGCAGG"
    "TCTCATTGCCAGAGCTGCCTGGGGAATACAGCATGAAAGTGACAGGAGAAGGATGTGTCT"
    "ACCTCCAGACATCCTTGAAATACAATATTCTCCCAGAAAAGGAAGAGTTCCCCTTTGCTT"
    "TAGGAGTGCAGACTCTGCCTCAAACTTGTGATGAACCCAAAGCCCACACCAGCTTCCAAA"
    "TCTCCCTAAGTGTCAGTTACACAGGGAGCCGCTCTGCCTCCAACATGGCGATCGTTGATG"
    "TGAAGATGGTCTCTGGCTTCATTCCCCTGAAGCCAACAGTGAAAATGCTTGAAAGATCTA"
    "ACCATGTGAGCCGGACAGAAGTCAGCAGCAACCATGTCTTGATTTACCTTGATAAGGTGT"
    "CAAATCAGACACTGAGCTTGTTCTTCACGGTTCTGCAAGATGTCCCAGTAAGAGATCTGA"
    "AACCAGCCATAGTGAAAGTCTATGATTACTACGAGACGGATGAGTTTGCAATTGCTGAGT"
    "ACAATGCTCCTTGCAGCAAAGATCTTGGAAATGCTTGAAGACCACAAGGCTGAAAAGTGC"
    "TTTGCTGGAGTCCTGTTCTCAGAGCTCCACAGAAGACACGTGTTTTTGTATCTTTAAAGA"
    "CTTGATGAATAAACACTTTTTCTGGTCAATGTCAAAAAAAAAAAAAAAAAAAAAAAAA";

static gchar *MIRNA_SEQUENCE = "UGAGGUAGUAGGUUGUAUAGUU";

static void
test_mirbooking ()
{
    g_autoptr (MirbookingBroker) mirbooking = mirbooking_broker_new ();

    g_autoptr (GBytes) precomputed_table = g_bytes_new_static (SCORE_TABLE, sizeof (SCORE_TABLE));
    g_autoptr (MirbookingPrecomputedScoreTable) score_table = mirbooking_precomputed_score_table_new_from_bytes (precomputed_table, 1, 7);
    mirbooking_broker_set_score_table (mirbooking, g_object_ref (score_table));

    g_autoptr (MirbookingTarget) target = mirbooking_target_new ("NM_000014.4");
    g_autoptr (MirbookingMirna) mirna = mirbooking_mirna_new ("MIMAT0000001");

    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (target), TARGET_SEQUENCE, strlen (TARGET_SEQUENCE));
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (mirna), MIRNA_SEQUENCE, strlen (MIRNA_SEQUENCE));

    mirbooking_broker_set_sequence_quantity (mirbooking, MIRBOOKING_SEQUENCE (target), 10);
    mirbooking_broker_set_sequence_quantity (mirbooking, MIRBOOKING_SEQUENCE (mirna), 10);

    g_autoptr (GError) error = NULL;
    g_assert (mirbooking_broker_run (mirbooking, &error));
    g_assert_null (error);

    GArray *target_sites = mirbooking_broker_get_target_sites (mirbooking);

    guint total_occupancy = 0;
    gint i;
    for (i = 0; i < target_sites->len; i++)
    {
        guint site_occupancy = 0;
        GSList *occupants_list = NULL;
        for (occupants_list = g_array_index (target_sites, MirbookingTargetSite, i).occupants; occupants_list != NULL; occupants_list = occupants_list->next)
        {
            MirbookingOccupant *occupant = occupants_list->data;
            site_occupancy += occupant->quantity;
        }

        total_occupancy += site_occupancy;
    }

    g_assert_cmpint (total_occupancy, ==, 10);

    // only one execution is permitted
    g_assert_cmpfloat (mirbooking_broker_run (mirbooking, &error), ==, 0.0f);
    g_assert_nonnull (error);
}

static void
broker_callback (GObject      *broker,
                 GAsyncResult *result,
                 gpointer      user_data)
{
    GError *err = NULL;

    if (mirbooking_broker_run_finish (MIRBOOKING_BROKER (broker),
                                      result,
                                      &err))
    {
        g_print ("pass");
    }
    else
    {
        g_print ("fail");
    }

    g_main_loop_quit (user_data);
}

static void
test_mirbooking_run_async ()
{
    if (g_test_subprocess ())
    {
        g_autoptr (MirbookingBroker) broker = mirbooking_broker_new ();

        g_autoptr (GMainLoop) loop = g_main_loop_new (NULL, FALSE);

        g_autoptr (GBytes) precomputed_table = g_bytes_new_static (SCORE_TABLE, sizeof (SCORE_TABLE));
        g_autoptr (MirbookingPrecomputedScoreTable) score_table = mirbooking_precomputed_score_table_new_from_bytes (precomputed_table, 1, 7);
        mirbooking_broker_set_score_table (broker, g_object_ref (score_table));

        mirbooking_broker_run_async (broker,
                                     broker_callback,
                                     loop);

        g_main_loop_run (loop);

        return;
    }

    g_test_trap_subprocess (NULL, 0, 0);
    g_test_trap_assert_passed ();
    g_test_trap_assert_stdout ("pass");
}

static void
test_mirbooking_empty ()
{
    g_autoptr (MirbookingBroker) broker = mirbooking_broker_new ();

    g_autoptr (GBytes) precomputed_table = g_bytes_new_static (SCORE_TABLE, sizeof (SCORE_TABLE));
    g_autoptr (MirbookingPrecomputedScoreTable) score_table = mirbooking_precomputed_score_table_new_from_bytes (precomputed_table, 1, 7);
    mirbooking_broker_set_score_table (broker, g_object_ref (score_table));

    g_assert (mirbooking_broker_run (broker, NULL));
}

static void
test_mirbooking_bad_seed_range ()
{
    g_autoptr (MirbookingBroker) broker = mirbooking_broker_new ();

    g_autoptr (GBytes) precomputed_table = g_bytes_new_static (SCORE_TABLE, sizeof (SCORE_TABLE));
    g_autoptr (MirbookingPrecomputedScoreTable) score_table = mirbooking_precomputed_score_table_new_from_bytes (precomputed_table, 18, 7);
    mirbooking_broker_set_score_table (broker, g_object_ref (score_table));

    g_autoptr (MirbookingTarget) target = mirbooking_target_new ("NM_000014.4");
    g_autoptr (MirbookingMirna) mirna = mirbooking_mirna_new ("MIMAT0000001");

    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (target), TARGET_SEQUENCE, strlen (TARGET_SEQUENCE));
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (mirna), MIRNA_SEQUENCE, strlen (MIRNA_SEQUENCE));

    mirbooking_broker_set_sequence_quantity (broker, MIRBOOKING_SEQUENCE (target), 10);
    mirbooking_broker_set_sequence_quantity (broker, MIRBOOKING_SEQUENCE (mirna), 10);

    GError *error = NULL;
    g_assert (mirbooking_broker_run (broker, &error));

    // ensure that no MREs has been assigned
    GArray *target_sites = mirbooking_broker_get_target_sites (broker);
    gint i;
    for (i = 0; i < target_sites->len; i++)
    {
        g_assert_null (g_array_index (target_sites, MirbookingTargetSite, i).occupants);
    }
}

gint
main (gint argc, gchar **argv)
{
    g_test_init (&argc, &argv, NULL);

    g_test_add_func ("/mirbooking", test_mirbooking);
    g_test_add_func ("/mirbooking/run-async", test_mirbooking_run_async);
    g_test_add_func ("/mirbooking/empty", test_mirbooking_empty);
    g_test_add_func ("/mirbooking/bad-seed-range", test_mirbooking_bad_seed_range);

    return g_test_run ();
}

