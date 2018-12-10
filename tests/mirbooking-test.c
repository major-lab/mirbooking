#include <glib.h>
#include <mirbooking.h>
#include <string.h>
#include <math.h>

#define R 1.987203611e-3
#define T 310.15

typedef struct __attribute__ ((packed)) _SparseMatrixMem
{
    gsize  n;
    gsize  nnz;
    gsize  rowptr[16384 + 1];
    gsize  colind[4];
    gfloat data[4];
} SeedScoreLayout;

static SeedScoreLayout SEED_SCORES =
{
    16384,
    4,
    {[8882] = 0, [8883] = 4},
    {4370,    4371,   4372,   9284},
    {-7.0f, -9.0f, -8.0f, -9.0f},
};

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

    g_autoptr (GBytes) default_table = g_bytes_new_static (&SEED_SCORES, sizeof (SEED_SCORES));
    g_autoptr (MirbookingDefaultScoreTable) score_table = mirbooking_default_score_table_new (default_table, MIRBOOKING_DEFAULT_SCORE_TABLE_DEFAULT_SUPPLEMENTARY_MODEL, NULL);
    mirbooking_broker_set_score_table (mirbooking, MIRBOOKING_SCORE_TABLE (g_object_ref (score_table)));

    g_autoptr (MirbookingTarget) target = mirbooking_target_new ("NM_000014.4");
    g_autoptr (MirbookingMirna) mirna = mirbooking_mirna_new ("MIMAT0000001");

    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (target), "GCACACA", strlen ("GCACACA"));
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (mirna), MIRNA_SEQUENCE, strlen (MIRNA_SEQUENCE));

    gdouble E0 = 4e6;
    gdouble S0 = 6e5;

    mirbooking_broker_set_sequence_quantity (mirbooking, MIRBOOKING_SEQUENCE (target), S0);
    mirbooking_broker_set_sequence_quantity (mirbooking, MIRBOOKING_SEQUENCE (mirna), E0);

    mirbooking_broker_set_5prime_footprint (mirbooking, 0);
    mirbooking_broker_set_3prime_footprint (mirbooking, 0);

    gdouble norm, last_norm;
    last_norm = 1.0/0.0;
    g_autoptr (GError) error = NULL;

    do
    {
        g_assert (mirbooking_broker_evaluate (mirbooking, &norm, &error));
        g_assert (mirbooking_broker_step (mirbooking, MIRBOOKING_BROKER_STEP_MODE_SOLVE_STEADY_STATE, 1, &error));
        g_assert_null (error);
        g_assert_cmpfloat (norm, <=, last_norm);
        last_norm = norm;

        GArray *target_sites = mirbooking_broker_get_target_sites (mirbooking);
        MirbookingTargetSite target_site = g_array_index (target_sites, MirbookingTargetSite, 0);
        MirbookingOccupant *occupant = target_site.occupants->data;
        gdouble ES = mirbooking_broker_get_occupant_quantity (mirbooking, occupant);
        g_assert_cmpfloat (ES, <=, 6e5);
        g_assert_cmpfloat (ES, <=, 4e6);
    }
    while (norm > 1e-6);

    GArray *target_sites = mirbooking_broker_get_target_sites (mirbooking);

    MirbookingTargetSite target_site = g_array_index (target_sites, MirbookingTargetSite, 0);

    g_assert_nonnull (target_site.occupants);
    g_assert_null (target_site.occupants->next);

    MirbookingOccupant *occupant = target_site.occupants->data;

    g_assert_cmpfloat (occupant->score, ==, 1e12 * exp ((-9.0f - 5.43f) / (R*T)));

    /* analytical solution for a single reaction */
    gdouble kf = MIRBOOKING_SCORE_TABLE_DEFAULT_KF;
    gdouble kr = kf * occupant->score;
    gdouble kcat = MIRBOOKING_SCORE_TABLE_DEFAULT_KCAT;
    gdouble z = kf * (E0 + S0) + kr + kcat;

    gdouble ES = mirbooking_broker_get_occupant_quantity (mirbooking, occupant);
    gdouble ES_eq = (z - sqrt (pow (z, 2) - 4 * pow (kf, 2) * E0 * S0)) / (2 * kf);
    g_assert_cmpfloat (fabs (ES_eq - ES), <, 1e-8);
}

static void
test_mirbooking_empty ()
{
    g_autoptr (MirbookingBroker) broker = mirbooking_broker_new ();

    g_autoptr (GBytes) default_table = g_bytes_new_static (&SEED_SCORES, sizeof (SEED_SCORES));
    g_autoptr (MirbookingDefaultScoreTable) score_table = mirbooking_default_score_table_new (default_table, MIRBOOKING_DEFAULT_SCORE_TABLE_DEFAULT_SUPPLEMENTARY_MODEL, NULL);
    mirbooking_broker_set_score_table (broker, MIRBOOKING_SCORE_TABLE (g_object_ref (score_table)));

    gdouble norm;
    g_assert (mirbooking_broker_evaluate (broker, &norm, NULL));
    g_assert (mirbooking_broker_step (broker, MIRBOOKING_BROKER_STEP_MODE_SOLVE_STEADY_STATE, 1, NULL));
}

static void
test_mirbooking_bad_seed_range ()
{
    g_autoptr (MirbookingBroker) broker = mirbooking_broker_new ();

    g_autoptr (GBytes) default_table = g_bytes_new_static (&SEED_SCORES, sizeof (SEED_SCORES));
    g_autoptr (MirbookingDefaultScoreTable) score_table = mirbooking_default_score_table_new (default_table, MIRBOOKING_DEFAULT_SCORE_TABLE_DEFAULT_SUPPLEMENTARY_MODEL, NULL);
    mirbooking_broker_set_score_table (broker, MIRBOOKING_SCORE_TABLE (g_object_ref (score_table)));

    g_autoptr (MirbookingTarget) target = mirbooking_target_new ("NM_000014.4");
    g_autoptr (MirbookingMirna) mirna = mirbooking_mirna_new ("MIMAT0000001");

    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (target), TARGET_SEQUENCE, 6);
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (mirna), MIRNA_SEQUENCE, strlen (MIRNA_SEQUENCE));

    mirbooking_broker_set_sequence_quantity (broker, MIRBOOKING_SEQUENCE (target), 10);
    mirbooking_broker_set_sequence_quantity (broker, MIRBOOKING_SEQUENCE (mirna), 10);

    if (g_test_subprocess ())
    {
        gdouble norm;
        GError *error = NULL;
        mirbooking_broker_evaluate (broker, &norm, &error);
        return;
    }

    g_test_trap_subprocess (NULL, 0, 0);
    g_test_trap_assert_failed ();
    g_test_trap_assert_stderr ("*offset + len <= priv->sequence_len - priv->sequence_skips->len*");
}

static void
test_mirbooking_numerical_integration ()
{
    g_autoptr (MirbookingBroker) broker = mirbooking_broker_new ();

    g_autoptr (GBytes) default_table = g_bytes_new_static (&SEED_SCORES, sizeof (SEED_SCORES));
    g_autoptr (MirbookingDefaultScoreTable) score_table = mirbooking_default_score_table_new (default_table, MIRBOOKING_DEFAULT_SCORE_TABLE_DEFAULT_SUPPLEMENTARY_MODEL, NULL);
    mirbooking_broker_set_score_table (broker, MIRBOOKING_SCORE_TABLE (g_object_ref (score_table)));

    g_autoptr (MirbookingTarget) target = mirbooking_target_new ("NM_000014.4");
    g_autoptr (MirbookingMirna) mirna = mirbooking_mirna_new ("MIMAT0000001");

    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (target), TARGET_SEQUENCE, strlen (TARGET_SEQUENCE));
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (mirna), MIRNA_SEQUENCE, strlen (MIRNA_SEQUENCE));

    mirbooking_broker_set_sequence_quantity (broker, MIRBOOKING_SEQUENCE (target), 10);
    mirbooking_broker_set_sequence_quantity (broker, MIRBOOKING_SEQUENCE (mirna), 10);

    gdouble norm;
    GError *error = NULL;
    g_assert (mirbooking_broker_evaluate (broker, &norm, &error));
    g_assert (isfinite (norm));

    GArray *target_sites = mirbooking_broker_get_target_sites (broker);

    MirbookingTargetSite target_site = g_array_index (target_sites, MirbookingTargetSite, 0);

    g_assert_nonnull (target_site.occupants);
    g_assert_null (target_site.occupants->next);

    MirbookingOccupant *occupant = target_site.occupants->data;

    /* duplex concentration should steadily increase toward equilibrium */
    gint i;
    gfloat last_norm = 1.0/0.0;
    gdouble step_size = 1e-3;
    for (i = 0; i < 10; i++)
    {
        g_assert (mirbooking_broker_evaluate (broker, &norm, &error));
        g_assert_null (error);
        g_assert_cmpfloat (norm, <, last_norm);
        last_norm = norm;

        g_assert (mirbooking_broker_step (broker, MIRBOOKING_BROKER_STEP_MODE_INTEGRATE, step_size, &error));
        g_assert_null (error);
        gdouble time = mirbooking_broker_get_time (broker);
        gdouble expected_t = (step_size * (i + 1));
        g_assert_cmpfloat (time - expected_t, <=, 1e-6 * expected_t + 1e-12);

        gdouble ES = mirbooking_broker_get_occupant_quantity (broker, occupant);

        g_assert_cmpfloat (ES, <, 10);
    }
}

gint
main (gint argc, gchar **argv)
{
    g_test_init (&argc, &argv, NULL);

    g_test_add_func ("/mirbooking", test_mirbooking);
    g_test_add_func ("/mirbooking/empty", test_mirbooking_empty);
    g_test_add_func ("/mirbooking/bad-seed-range", test_mirbooking_bad_seed_range);
    g_test_add_func ("/mirbooking/numerical-integration", test_mirbooking_numerical_integration);

    return g_test_run ();
}
