#include <glib.h>
#include <mirbooking.h>
#include <string.h>
#include <math.h>

#define R 1.987203611e-3
#define T 310.15

#if !GLIB_CHECK_VERSION(2,58,0)
#define g_assert_cmpfloat_with_epsilon(n1,n2,epsilon) g_assert_true(fabs((n1) - (n2)) < (epsilon))
#endif

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

    g_assert_cmpint (mirbooking_broker_get_rank (mirbooking), ==, 0);

    g_autoptr (GBytes) default_table = g_bytes_new_static (&SEED_SCORES, sizeof (SEED_SCORES));
    g_autoptr (MirbookingDefaultScoreTable) score_table = mirbooking_default_score_table_new (default_table, MIRBOOKING_DEFAULT_SCORE_TABLE_DEFAULT_SUPPLEMENTARY_MODEL, NULL);
    mirbooking_broker_set_score_table (mirbooking, MIRBOOKING_SCORE_TABLE (g_object_ref (score_table)));

    g_autoptr (MirbookingTarget) target = mirbooking_target_new ("NM_000014.4");
    g_autoptr (MirbookingMirna) mirna = mirbooking_mirna_new ("MIMAT0000001");

    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (target), "GCACACA");
    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (mirna), MIRNA_SEQUENCE);

    gdouble E0 = 4e6;
    gdouble S0 = 6e5;

    mirbooking_broker_set_sequence_quantity (mirbooking, MIRBOOKING_SEQUENCE (target), S0);
    mirbooking_broker_set_sequence_quantity (mirbooking, MIRBOOKING_SEQUENCE (mirna), E0);

    g_assert_cmpfloat (mirbooking_broker_get_sequence_quantity (mirbooking, MIRBOOKING_SEQUENCE (target)), ==, S0);
    g_assert_cmpfloat (mirbooking_broker_get_sequence_quantity (mirbooking, MIRBOOKING_SEQUENCE (mirna)), ==, E0);

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

        const GArray *target_sites = mirbooking_broker_get_target_sites (mirbooking);
        MirbookingTargetSite target_site = g_array_index (target_sites, MirbookingTargetSite, 0);
        MirbookingOccupant *occupant = target_site.occupants->data;
        gdouble ES = mirbooking_broker_get_occupant_quantity (mirbooking, occupant);
        g_assert_cmpfloat (ES, <=, 6e5);
        g_assert_cmpfloat (ES, <=, 4e6);
    }
    while (norm > 1e-6);

    // at steady-state, the system is closed so we get equal incoming,
    // catalyzing and outgoing product
    gdouble ktr = mirbooking_broker_get_target_transcription_rate (mirbooking, target);
    gdouble kdeg = mirbooking_broker_get_product_degradation_rate (mirbooking, target);
    g_assert_cmpfloat (ktr, >, 0);
    g_assert_cmpfloat (kdeg, >, 0);
    g_assert_cmpfloat_with_epsilon (ktr, kdeg, 1e-8);
    g_assert_cmpfloat (mirbooking_broker_get_product_quantity (mirbooking, target), >, 0);

    const GArray *target_sites = mirbooking_broker_get_target_sites (mirbooking);

    MirbookingTargetSite *target_site = &g_array_index (target_sites, MirbookingTargetSite, 0);

    g_assert_nonnull (target_site->occupants);
    g_assert_null (target_site->occupants->next);

    MirbookingOccupant *occupant = target_site->occupants->data;

    g_assert_cmpfloat (MIRBOOKING_SCORE_KD (occupant->score), ==, 1e12 * exp ((-9.0f - 5.43f) / (R*T)));

    /* analytical solution for a single reaction */
    gdouble z = occupant->score.kf * (E0 + S0) + occupant->score.kr + occupant->score.kcat;

    gdouble ES = mirbooking_broker_get_occupant_quantity (mirbooking, occupant);
    gdouble ES_eq = (z - sqrt (pow (z, 2) - 4 * pow (occupant->score.kf, 2) * E0 * S0)) / (2 * occupant->score.kf);
    g_assert_cmpfloat (fabs (ES_eq - ES), <, 1e-8);

    gsize pmf_len;
    g_autofree gdouble *pmf = mirbooking_broker_get_target_occupants_pmf (mirbooking, target, &pmf_len);
    g_assert_cmpint (pmf_len, ==, 2);
    g_assert_cmpfloat (fabs (pmf[0] - (1 - (ES_eq / S0))), <=, 1e-12);
    g_assert_cmpfloat (fabs (pmf[1] - (ES_eq / S0)), <=, 1e-12);

    // substrate level is steady
    g_assert_cmpfloat (mirbooking_broker_get_sequence_quantity (mirbooking, MIRBOOKING_SEQUENCE (target)), ==, S0);

    // enzyme and position-wise concentration are conserved
    g_assert_cmpfloat_with_epsilon (mirbooking_broker_get_sequence_quantity (mirbooking, MIRBOOKING_SEQUENCE (mirna)), E0 - ES, 1e-1);
    g_assert_cmpfloat_with_epsilon (mirbooking_broker_get_target_site_quantity (mirbooking, target_site), S0 - ES, 1e-8);
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

    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (target), g_bytes_new_static (TARGET_SEQUENCE, 6));
    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (mirna), MIRNA_SEQUENCE);

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
    g_test_trap_assert_stderr ("*offset + len <= g_bytes_get_size (priv->sequence) - priv->sequence_skips->len*");
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

    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (target), TARGET_SEQUENCE);
    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (mirna), MIRNA_SEQUENCE);

    mirbooking_broker_set_sequence_quantity (broker, MIRBOOKING_SEQUENCE (target), 10);
    mirbooking_broker_set_sequence_quantity (broker, MIRBOOKING_SEQUENCE (mirna), 10);

    gdouble norm;
    GError *error = NULL;
    g_assert (mirbooking_broker_evaluate (broker, &norm, &error));
    g_assert (isfinite (norm));

    const GArray *target_sites = mirbooking_broker_get_target_sites (broker);

    MirbookingTargetSite target_site = g_array_index (target_sites, MirbookingTargetSite, 0);

    g_assert_nonnull (target_site.occupants);
    g_assert_null (target_site.occupants->next);

    MirbookingOccupant *occupant = target_site.occupants->data;

    /* duplex concentration should steadily increase toward equilibrium */
    gint i;
    gfloat last_norm = 1.0/0.0;
    gdouble step_size = 1;
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

static void
test_mirbooking_solve_and_integrate ()
{
    g_autoptr (MirbookingBroker) broker = mirbooking_broker_new ();

    g_autoptr (GBytes) default_table = g_bytes_new_static (&SEED_SCORES, sizeof (SEED_SCORES));
    g_autoptr (MirbookingDefaultScoreTable) score_table = mirbooking_default_score_table_new (default_table, MIRBOOKING_DEFAULT_SCORE_TABLE_DEFAULT_SUPPLEMENTARY_MODEL, NULL);
    mirbooking_broker_set_score_table (broker, MIRBOOKING_SCORE_TABLE (g_object_ref (score_table)));

    g_autoptr (MirbookingTarget) target = mirbooking_target_new ("NM_000014.4");
    g_autoptr (MirbookingMirna) mirna = mirbooking_mirna_new ("MIMAT0000001");

    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (target), TARGET_SEQUENCE);
    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (mirna), MIRNA_SEQUENCE);

    mirbooking_broker_set_sequence_quantity (broker, MIRBOOKING_SEQUENCE (target), 10);
    mirbooking_broker_set_sequence_quantity (broker, MIRBOOKING_SEQUENCE (mirna), 10);

    gdouble norm;
    GError *error = NULL;
    g_assert (mirbooking_broker_evaluate (broker, &norm, &error));
    g_assert (isfinite (norm));

    const GArray *target_sites = mirbooking_broker_get_target_sites (broker);

    MirbookingTargetSite *target_site = &g_array_index (target_sites, MirbookingTargetSite, 0);

    g_assert_nonnull (target_site->occupants);
    g_assert_null (target_site->occupants->next);

    MirbookingOccupant *occupant = target_site->occupants->data;

    do
    {
        mirbooking_broker_evaluate (broker, &norm, &error);
        g_assert_null (error);
        mirbooking_broker_step (broker, MIRBOOKING_BROKER_STEP_MODE_SOLVE_STEADY_STATE, 1.0, &error);
        g_assert_null (error);
    }
    while (norm > 1e-16);

    // steady-state assumption holds
    g_assert_cmpfloat (mirbooking_broker_get_sequence_quantity (broker, MIRBOOKING_SEQUENCE (target)), ==, mirbooking_broker_get_sequence_quantity (broker, MIRBOOKING_SEQUENCE (target)));

    gdouble ES = mirbooking_broker_get_occupant_quantity (broker, occupant);

    gdouble E = 10;
    const GArray *occupants = mirbooking_broker_get_occupants (broker);
    guint i;
    for (i = 0; i < occupants->len; i++)
    {
        E -= mirbooking_broker_get_occupant_quantity (broker, &g_array_index (occupants, MirbookingOccupant, i));
    }

    // enzyme is conserved
    g_assert_cmpfloat_with_epsilon (mirbooking_broker_get_sequence_quantity (broker, MIRBOOKING_SEQUENCE (mirna)), E, 1e-6);

    // complex is in equilibrium
    gdouble kother = mirbooking_broker_get_target_site_kother (broker, target_site);
    g_assert_cmpfloat (kother, >, 0);
    gdouble Stp = mirbooking_broker_get_target_site_quantity (broker, target_site);
    g_assert_cmpfloat_with_epsilon ((E * Stp) / ES, MIRBOOKING_SCORE_KM (occupant->score) + kother / occupant->score.kf, 1e-12);

    // steady-state is maintained
    gdouble ktr = mirbooking_broker_get_target_transcription_rate (broker, target);
    gdouble kdeg = mirbooking_broker_get_product_degradation_rate (broker, target);
    g_assert_cmpfloat (ktr, >, 0);
    g_assert_cmpfloat (kdeg, >, 0);
    g_assert_cmpfloat (fabs (ktr - kdeg), <=, 1e-12);

    // integrating at steady-state must preserve the steady-state
    g_assert (mirbooking_broker_step (broker, MIRBOOKING_BROKER_STEP_MODE_INTEGRATE, 120.0, &error));
    g_assert_cmpfloat_with_epsilon (mirbooking_broker_get_sequence_quantity (broker, MIRBOOKING_SEQUENCE (mirna)), E, 1e-6);
    g_assert_cmpfloat_with_epsilon (mirbooking_broker_get_occupant_quantity (broker, occupant), ES, 1e-6);
    mirbooking_broker_evaluate (broker, &norm, &error);
    g_assert_cmpfloat (norm, <=, 1e-12);

    // over-expression of MIMAT0000001 (10 -> 20)
    mirbooking_broker_set_sequence_quantity (broker, MIRBOOKING_SEQUENCE (mirna), 20);
    g_assert_cmpfloat (mirbooking_broker_get_sequence_quantity (broker, MIRBOOKING_SEQUENCE (mirna)), ==, 20);

    mirbooking_broker_evaluate (broker, &norm, &error);
    g_assert_cmpfloat (norm, >, 1e-6);

    gdouble step_size = 1;
    for (i = 120; i < 130; i++)
    {
        g_assert (mirbooking_broker_evaluate (broker, &norm, &error));
        g_assert_null (error);
        g_assert (isfinite (norm));

        g_assert (mirbooking_broker_step (broker, MIRBOOKING_BROKER_STEP_MODE_INTEGRATE, step_size, &error));
        g_assert_null (error);
        gdouble time = mirbooking_broker_get_time (broker);
        gdouble expected_t = (step_size * (i + 1));
        g_assert_cmpfloat (time - expected_t, <=, 1e-6 * expected_t + 1e-12);

        ES = mirbooking_broker_get_occupant_quantity (broker, occupant);

        g_assert_cmpfloat (ES, >=, 0);
        g_assert_cmpfloat (ES, <=, 10);
    }

    // some of the substrate will have degraded
    gdouble S = mirbooking_broker_get_sequence_quantity (broker, MIRBOOKING_SEQUENCE (target));
    gdouble P = mirbooking_broker_get_product_quantity (broker, target);
    g_assert_cmpfloat (P, >, 0);
    g_assert_cmpfloat (S, <, 10);
}

static void
test_mirbooking_set_occupant_quantity ()
{
    g_autoptr (MirbookingBroker) broker = mirbooking_broker_new ();

    g_autoptr (GBytes) default_table = g_bytes_new_static (&SEED_SCORES, sizeof (SEED_SCORES));
    g_autoptr (MirbookingDefaultScoreTable) score_table = mirbooking_default_score_table_new (default_table, MIRBOOKING_DEFAULT_SCORE_TABLE_DEFAULT_SUPPLEMENTARY_MODEL, NULL);
    mirbooking_broker_set_score_table (broker, MIRBOOKING_SCORE_TABLE (g_object_ref (score_table)));

    g_autoptr (MirbookingTarget) target = mirbooking_target_new ("NM_000014.4");
    g_autoptr (MirbookingMirna) mirna = mirbooking_mirna_new ("MIMAT0000001");

    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (target), "GCACACA");
    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (mirna), MIRNA_SEQUENCE);

    mirbooking_broker_set_sequence_quantity (broker, MIRBOOKING_SEQUENCE (target), 10);
    mirbooking_broker_set_sequence_quantity (broker, MIRBOOKING_SEQUENCE (mirna), 10);

    gdouble norm;
    GError *error = NULL;
    g_assert (mirbooking_broker_evaluate (broker, &norm, &error));
    g_assert (isfinite (norm));

    const GArray *occupants = mirbooking_broker_get_occupants (broker);
    MirbookingOccupant *occupant = &g_array_index (occupants, MirbookingOccupant, 0);

    {
        gdouble E = mirbooking_broker_get_sequence_quantity (broker, MIRBOOKING_SEQUENCE (occupant->mirna));
        gdouble ES = mirbooking_broker_get_occupant_quantity (broker, occupant);
        mirbooking_broker_set_sequence_quantity (broker, MIRBOOKING_SEQUENCE (occupant->mirna), E - 10);
        mirbooking_broker_set_occupant_quantity (broker, occupant, ES + 10);
    }

    g_assert_cmpfloat (mirbooking_broker_get_occupant_quantity (broker, occupant), ==, 10);

    do
    {
        mirbooking_broker_evaluate (broker, &norm, &error);
        g_assert_null (error);
        mirbooking_broker_step (broker, MIRBOOKING_BROKER_STEP_MODE_SOLVE_STEADY_STATE, 1.0, &error);
        g_assert_null (error);
    }
    while (norm > 1e-6);

    // ensure that we converge to the same equilibrium point
    gdouble ES = mirbooking_broker_get_occupant_quantity (broker, occupant);
    gdouble Z = 10 + 10 + MIRBOOKING_SCORE_KM (occupant->score);
    gdouble ES_eq = (Z - sqrt (pow (Z, 2) - 4 * 10 * 10)) / 2;

    g_assert_cmpfloat (fabs (ES - ES_eq), <=, 1e-6);
}

/**
 * The purpose of this test is to guarantee that it's possible to restore
 * broker from what is being publicly exposed.
 */
static void
test_mirbooking_restore_broker_state ()
{
    g_autoptr (MirbookingBroker) broker1 = mirbooking_broker_new ();
    g_autoptr (MirbookingBroker) broker2 = mirbooking_broker_new ();

    g_autoptr (GBytes) default_table = g_bytes_new_static (&SEED_SCORES, sizeof (SEED_SCORES));
    g_autoptr (MirbookingDefaultScoreTable) score_table = mirbooking_default_score_table_new (default_table, MIRBOOKING_DEFAULT_SCORE_TABLE_DEFAULT_SUPPLEMENTARY_MODEL, NULL);
    mirbooking_broker_set_score_table (broker1, MIRBOOKING_SCORE_TABLE (g_object_ref (score_table)));
    mirbooking_broker_set_score_table (broker2, MIRBOOKING_SCORE_TABLE (g_object_ref (score_table)));

    g_autoptr (MirbookingTarget) target = mirbooking_target_new ("NM_000014.4");
    g_autoptr (MirbookingMirna) mirna = mirbooking_mirna_new ("MIMAT0000001");

    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (target), "GCACACA");
    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (mirna), MIRNA_SEQUENCE);

    mirbooking_broker_set_sequence_quantity (broker1, MIRBOOKING_SEQUENCE (target), 10);
    mirbooking_broker_set_sequence_quantity (broker1, MIRBOOKING_SEQUENCE (mirna), 10);

    g_autoptr (GError) error = NULL;
    gdouble norm;
    mirbooking_broker_evaluate (broker1, &norm, &error);
    g_assert_null (error);
    mirbooking_broker_step (broker1, MIRBOOKING_BROKER_STEP_MODE_SOLVE_STEADY_STATE, 1.0, &error);
    g_assert_null (error);
    mirbooking_broker_evaluate (broker1, &norm, &error);

    // restore broker2 from broker1 state
    const GPtrArray *mirnas = mirbooking_broker_get_mirnas (broker1);
    const GPtrArray *targets = mirbooking_broker_get_targets (broker1);

    guint i;
    for (i = 0; i < mirnas->len; i++)
    {
        MirbookingMirna *mirna = g_ptr_array_index (mirnas, i);
        mirbooking_broker_set_sequence_quantity (broker2,
                                                 MIRBOOKING_SEQUENCE (mirna),
                                                 mirbooking_broker_get_sequence_quantity (broker1, MIRBOOKING_SEQUENCE (mirna)));
    }

    for (i = 0; i < targets->len; i++)
    {
        MirbookingTarget *target = g_ptr_array_index (targets, i);
        mirbooking_broker_set_sequence_quantity (broker2,
                                                 MIRBOOKING_SEQUENCE (target),
                                                 mirbooking_broker_get_sequence_quantity (broker1, MIRBOOKING_SEQUENCE (target)));
    }

    mirbooking_broker_evaluate (broker2, NULL, &error);
    g_assert_null (error);

    const GArray *broker1_occupants = mirbooking_broker_get_occupants (broker1);
    const GArray *broker2_occupants = mirbooking_broker_get_occupants (broker2);
    g_assert_cmpint (broker1_occupants->len, ==, broker2_occupants->len);

    for (i = 0; i < broker1_occupants->len; i++)
    {
        MirbookingOccupant *broker1_occupant = &g_array_index (broker1_occupants, MirbookingOccupant, i);
        MirbookingOccupant *broker2_occupant = &g_array_index (broker2_occupants, MirbookingOccupant, i);

        g_assert (broker1_occupant->mirna == broker2_occupant->mirna);
        g_assert (broker1_occupant->target == broker2_occupant->target);
        g_assert (broker1_occupant->position == broker2_occupant->position);

        mirbooking_broker_set_occupant_quantity (broker2, broker2_occupant, mirbooking_broker_get_occupant_quantity (broker1, broker1_occupant));
    }

    // trigger steady-state assumption to adjust transcription rates and
    // product quantity
    mirbooking_broker_step (broker2, MIRBOOKING_BROKER_STEP_MODE_SOLVE_STEADY_STATE, 0.0, &error);
    g_assert_null (error);

    for (i = 0; i < targets->len; i++)
    {
        MirbookingTarget *target = g_ptr_array_index (targets, i);
        g_assert_cmpfloat (mirbooking_broker_get_target_transcription_rate (broker1, target), ==,
                           mirbooking_broker_get_target_transcription_rate (broker2, target));
        g_assert_cmpfloat (mirbooking_broker_get_product_quantity (broker1, target), ==,
                           mirbooking_broker_get_product_quantity (broker2, target));
    }

    for (i = 0; i < broker1_occupants->len; i++)
    {
        MirbookingOccupant *broker1_occupant = &g_array_index (broker1_occupants, MirbookingOccupant, i);
        MirbookingOccupant *broker2_occupant = &g_array_index (broker2_occupants, MirbookingOccupant, i);

        g_assert_cmpfloat (mirbooking_broker_get_occupant_quantity (broker1, broker1_occupant), ==, mirbooking_broker_get_occupant_quantity (broker2, broker2_occupant));
    }

    gdouble broker2_norm;
    mirbooking_broker_evaluate (broker2, &broker2_norm, &error);
    g_assert_null (error);

    g_assert_cmpfloat_with_epsilon (norm, broker2_norm, 1e-6);
}

static void
test_mirbooking_target_knock_out ()
{
    g_autoptr (MirbookingBroker) broker = mirbooking_broker_new ();

    g_autoptr (GBytes) default_table = g_bytes_new_static (&SEED_SCORES, sizeof (SEED_SCORES));
    g_autoptr (MirbookingDefaultScoreTable) score_table = mirbooking_default_score_table_new (default_table, MIRBOOKING_DEFAULT_SCORE_TABLE_DEFAULT_SUPPLEMENTARY_MODEL, NULL);

    g_autoptr (MirbookingTarget) target = mirbooking_target_new ("NM_000014.4");
    g_autoptr (MirbookingMirna) mirna = mirbooking_mirna_new ("MIMAT0000001");

    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (target), "GCACACA");
    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (mirna), MIRNA_SEQUENCE);

    mirbooking_broker_set_sequence_quantity (broker, MIRBOOKING_SEQUENCE (target), 10);
    mirbooking_broker_set_sequence_quantity (broker, MIRBOOKING_SEQUENCE (mirna), 10);

    mirbooking_broker_set_score_table (broker, MIRBOOKING_SCORE_TABLE (score_table));

    gdouble norm;
    do
    {
        mirbooking_broker_step (broker, MIRBOOKING_BROKER_STEP_MODE_SOLVE_STEADY_STATE, 1.0, NULL);
        mirbooking_broker_evaluate (broker, &norm, NULL);
    }
    while (norm >= 1e-6);

    // stop target transcription
    mirbooking_broker_set_target_transcription_rate (broker, target, 0);

    mirbooking_broker_step (broker, MIRBOOKING_BROKER_STEP_MODE_INTEGRATE, 2 * 3600.0, NULL);
    mirbooking_broker_evaluate (broker, &norm, NULL);

    g_assert_cmpfloat_with_epsilon (mirbooking_broker_get_sequence_quantity (broker, MIRBOOKING_SEQUENCE (target)), 0, 1e-3);
}

gint
main (gint argc, gchar **argv)
{
    g_test_init (&argc, &argv, NULL);

    g_test_add_func ("/mirbooking", test_mirbooking);
    g_test_add_func ("/mirbooking/empty", test_mirbooking_empty);
    g_test_add_func ("/mirbooking/bad-seed-range", test_mirbooking_bad_seed_range);
    g_test_add_func ("/mirbooking/numerical-integration", test_mirbooking_numerical_integration);
    g_test_add_func ("/mirbooking/solve-and-integrate", test_mirbooking_solve_and_integrate);
    g_test_add_func ("/mirbooking/set-occupant-quantity", test_mirbooking_set_occupant_quantity);
    g_test_add_func ("/mirbooking/restore-broker-state", test_mirbooking_restore_broker_state);
    g_test_add_func ("/mirbooking/target-knock-out", test_mirbooking_target_knock_out);

    return g_test_run ();
}
