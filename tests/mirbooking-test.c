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
    g_autoptr (Mirbooking) mirbooking = mirbooking_new ();

    mirbooking_set_threshold (mirbooking, 0.0f); // the score table is zeroed anyway

    g_autoptr (GBytes) precomputed_table = g_bytes_new_static (SCORE_TABLE, sizeof (SCORE_TABLE));
    g_autoptr (MirbookingScoreTable) score_table = mirbooking_score_table_new_precomputed (precomputed_table);
    mirbooking_set_score_table (mirbooking, score_table);

    g_autoptr (MirbookingTarget) target = mirbooking_target_new ("NM_000014.4");
    g_autoptr (MirbookingMirna) mirna = mirbooking_mirna_new ("MIMAT0000001");

    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (target), TARGET_SEQUENCE, strlen (TARGET_SEQUENCE));
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (mirna), MIRNA_SEQUENCE, strlen (MIRNA_SEQUENCE));

    mirbooking_set_sequence_quantity (mirbooking, MIRBOOKING_SEQUENCE (target), 10);
    mirbooking_set_sequence_quantity (mirbooking, MIRBOOKING_SEQUENCE (mirna), 10);

    GError *error = NULL;
    g_assert (mirbooking_run (mirbooking, &error));
    g_assert_null (error);

    gsize target_sites_len;
    MirbookingTargetSite *target_sites = mirbooking_get_target_sites (mirbooking,
                                                                      &target_sites_len);

    guint total_occupancy = 0;
    gint i;
    for (i = 0; i < target_sites_len; i++)
    {
        if (target_sites[i].occupancy)
            g_print ("%s %u\n", "a", target_sites[i].occupancy);

        total_occupancy += target_sites[i].occupancy;

        // g_print ("target %s position %d %lu\n", mirbooking_sequence_get_accession (MIRBOOKING_SEQUENCE (target_sites[i].target)), i, target_sites[i].occupancy);
        /*
        if (i == 4079 || i == 2461 || i == 4155 || i == 4234)
        {
            g_assert_nonnull (target_sites[i].occupancy);
            g_assert_cmpint (((MirbookingOccupant*) target_sites[i].occupants->data)->quantity, ==, 1);
        }
        else
        {
            g_print ("at site %i\n", i);
            g_assert_null (target_sites[i].occupants);
        }
        */
    }

    g_assert_cmpint (total_occupancy, ==, 10);
}

gint
main (gint argc, gchar **argv)
{
    g_test_init (&argc, &argv, NULL);

    g_test_add_func ("/mirbooking", test_mirbooking);

    return g_test_run ();
}
