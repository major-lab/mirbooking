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
    gsize  colind[5];
    gfloat data[5];
} SeedScoreLayout;

static SeedScoreLayout SEED_SCORES =
{
    16384,
    5,
    {[8882] = 0, [8883] = 5},
    {4370,    4371,   4372,   7261,  9284},
    {-20.0f, -19.0f, -18.0f, -9.30f, -19.0f},
};

static void
test_score_table_compute_seed_score ()
{
    g_autoptr (GBytes) default_table = g_bytes_new_static (&SEED_SCORES, sizeof (SEED_SCORES));
    g_autoptr (MirbookingDefaultScoreTable) score_table = mirbooking_default_score_table_new_from_bytes (default_table, 1, 7);

    g_autoptr (MirbookingTarget) target = mirbooking_target_new ("NM");
    g_autoptr (MirbookingMirna) mirna = mirbooking_mirna_new ("MIMAT");

    gchar *target_seq = "GCACACAGAGCAGCATAAAGCCCAGTTGCTTTGGGAAGTGTTTGGGACCAGATGGATTGT";
    gchar *mirna_seq = "UGAGGUAGUAGGUUGUAUAGUU";

    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (target), target_seq, strlen (target_seq));
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (mirna), mirna_seq, strlen (mirna_seq));

    g_assert_cmpmem (mirbooking_sequence_get_subsequence (MIRBOOKING_SEQUENCE (target), 1, 7), 7, "CACACAG", 7);
    g_assert_cmpmem (mirbooking_sequence_get_subsequence (MIRBOOKING_SEQUENCE (mirna), 1, 7), 7, "GAGGUAG", 7);

    // test for a seed
    g_autoptr (GError) error = NULL;
    gdouble site_score = mirbooking_score_table_compute_score (MIRBOOKING_SCORE_TABLE (score_table),
                                                               mirna,     // GAGGUAG =>
                                                               target, 1, // CACACAG =>
                                                               &error);

    g_assert_cmpfloat (site_score, ==, 1e12 * exp ((-20.0f - 5.72) / (R * T)));
    g_assert_null (error);
}

static void
test_score_table_compute_seed_scores ()
{
    g_autoptr (GBytes) default_table = g_bytes_new_static (&SEED_SCORES, sizeof (SEED_SCORES));
    g_autoptr (MirbookingDefaultScoreTable) score_table = mirbooking_default_score_table_new_from_bytes (default_table, 1, 7);

    g_autoptr (MirbookingTarget) target = mirbooking_target_new ("NM");
    g_autoptr (MirbookingMirna) mirna = mirbooking_mirna_new ("MIMAT");

    gchar *target_seq = "GCACACAGAGCAGCATAAAGCCCAGTTGCTTTGGGAAGTGTTTGGGACCAGATGGATTGT";
    gchar *mirna_seq = "UGAGGUAGUAGGUUGUAUAGUU";

    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (target), target_seq, strlen (target_seq));
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (mirna), mirna_seq, strlen (mirna_seq));

    g_assert_cmpmem (mirbooking_sequence_get_subsequence (MIRBOOKING_SEQUENCE (target), 1, 7), 7, "CACACAG", 7);
    g_assert_cmpmem (mirbooking_sequence_get_subsequence (MIRBOOKING_SEQUENCE (mirna), 1, 7), 7, "GAGGUAG", 7);

    // test for a seed
    g_autoptr (GError) error = NULL;
    gsize *positions;
    gsize  positions_len;
    g_assert (mirbooking_score_table_compute_positions (MIRBOOKING_SCORE_TABLE (score_table),
                                                        mirna,     // GAGGUAG =>
                                                        target, &positions, // CACACAG =>
                                                        &positions_len,
                                                        &error));

    g_assert_nonnull (positions);
    g_assert_cmpint (positions_len, ==, 2);
    g_assert_cmpfloat (mirbooking_score_table_compute_score (score_table, mirna, target, positions[0], NULL), ==, 1e12 * exp ((-19.0f - 5.72) / (R * T)));
    g_assert_cmpfloat (mirbooking_score_table_compute_score (score_table, mirna, target, positions[1], NULL), ==, 1e12 * exp ((-20.0f - 5.72) / (R * T)));
    g_assert_null (error);

    g_free (positions);
}

static void
test_score_table_mcff ()
{
    g_autoptr (MirbookingMcffScoreTable) score_table = mirbooking_mcff_score_table_new ();
    g_assert_nonnull (score_table);

    g_autoptr (MirbookingTarget) target = mirbooking_target_new ("NM");
    g_autoptr (MirbookingMirna) mirna = mirbooking_mirna_new ("MIMAT");

    gchar *target_seq = "GCACACAGAGCAGCATAAAGCCCAGTTGCTTTGGGAAGTGTTTGGGACCAGATGGATTGT";
    gchar *mirna_seq = "UGAGGUAGUAGGUUGUAUAGUU";

    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (target), target_seq, strlen (target_seq));
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (mirna), mirna_seq, strlen (mirna_seq));

    // test for a seed
    g_autoptr (GError) error = NULL;
    gdouble site_score = mirbooking_score_table_compute_score (MIRBOOKING_SCORE_TABLE (score_table),
                                                               mirna,     // GAGGUAG =>
                                                               target, 1, // CACACAG =>
                                                               &error);

    g_assert_null (error);
    g_assert_cmpfloat (site_score, ==, 1e12 * exp (-7.09f / (R * T)));
}

static void
test_score_table_wee_et_al_2012 ()
{
    gdouble Kd, Km;
    // ensure that we reproduce experimental results from (Wee et al. 2012)

    g_autoptr (MirbookingMirna) mirna = mirbooking_mirna_new ("let-7b");
    g_autoptr (MirbookingTarget) target = mirbooking_target_new ("reporter");

    g_autoptr (GBytes) default_table = g_bytes_new_static (&SEED_SCORES, sizeof (SEED_SCORES));
    g_autoptr (MirbookingScoreTable) score_table = MIRBOOKING_SCORE_TABLE (mirbooking_default_score_table_new_from_bytes (default_table, 1, 7));

    // siRNA
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (mirna), "UGAGGUAGUAGGUUGUAUAUGU", strlen ("UGAGGUAGUAGGUUGUAUAUGU"));
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (target), "GAUACUAUACAACCUACUACCUCUAAAU", strlen ("GAUACUAUACAACCUACUACCUCUAAAU"));

    gsize *positions;
    gsize  positions_len;
    mirbooking_score_table_compute_positions (MIRBOOKING_SCORE_TABLE (score_table), mirna, target, &positions, &positions_len, NULL);
    g_assert_cmpint (positions_len, ==, 1);
    g_assert_cmpint (positions[0], ==, 16);

    Kd = mirbooking_score_table_compute_score (score_table, mirna, target, 16, NULL);
    Km = mirbooking_score_table_compute_enzymatic_score (score_table, mirna, target, 16, NULL);

    g_assert_cmpfloat (Kd, ==, 1e12 * exp ((-9.30f - 5.72) / (R * T)));
    g_assert_cmpfloat (Kd, >=, 20 - 10);
    g_assert_cmpfloat (Kd, <=, 20 + 10);

    g_assert_cmpfloat (Km, ==, Kd + (MIRBOOKING_SCORE_TABLE_DEFAULT_KCAT / MIRBOOKING_SCORE_TABLE_DEFAULT_KF));
    g_assert_cmpfloat (Km, >=, 100 - 60);
    g_assert_cmpfloat (Km, <=, 100 + 60);

    // seed-only
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (target), "GAAAAAAAAAAAAAAUCUACCUCUAAAU", strlen ("GAAAAAAAAAAAAAAUCUACCUCUAAAU"));

    Kd = mirbooking_score_table_compute_score (score_table, mirna, target, 16, NULL);
    Km = mirbooking_score_table_compute_enzymatic_score (score_table, mirna, target, 16, NULL);

    g_assert_cmpfloat (Kd, ==, 1e12 * exp ((-9.30f - 5.72) / (R * T)));
    g_assert_cmpfloat (Kd, >=, 26 - 2);
    g_assert_cmpfloat (Kd, <=, 26 + 2);

    g_assert_cmpfloat (Km, ==, Kd + (MIRBOOKING_SCORE_TABLE_DEFAULT_KCAT / MIRBOOKING_SCORE_TABLE_DEFAULT_KF));
    g_assert_cmpfloat (Km, >=, 100 - 60);
    g_assert_cmpfloat (Km, <=, 100 + 60);

    // seed and supplementary
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (target), "GAAAAAAAACAAAAAUCUACCUCUAAAU", strlen ("GAAAAAAAAAAAAAAUCUACCUCUAAAU"));

    // FIXME: add supplementary scores in the model
    SEED_SCORES.data[3] += -0.20f;

    Kd = mirbooking_score_table_compute_score (score_table, mirna, target, 16, NULL);
    Km = mirbooking_score_table_compute_enzymatic_score (score_table, mirna, target, 16, NULL);

    g_assert_cmpfloat (Kd, ==, 1e12 * exp ((-9.30f - 0.20f - 5.72) / (R * T)));
    // FIXME:
    // g_assert_cmpfloat (Kd, >=, 13 - 1);
    // g_assert_cmpfloat (Kd, <=, 13 + 1);

    g_assert_cmpfloat (Km, ==, Kd + (MIRBOOKING_SCORE_TABLE_DEFAULT_KCAT / MIRBOOKING_SCORE_TABLE_DEFAULT_KF));
    g_assert_cmpfloat (Km, >=, 100 - 60);
    g_assert_cmpfloat (Km, <=, 100 + 60);

    g_free (positions);
}

int main (int argc, gchar **argv)
{
    g_test_init (&argc, &argv, NULL);

    g_test_add_func ("/score-table/compute-seed-score", test_score_table_compute_seed_score);
    g_test_add_func ("/score-table/compute-seed-scores", test_score_table_compute_seed_scores);
    g_test_add_func ("/score-table/mcff", test_score_table_mcff);
    g_test_add_func ("/score-table/wee-et-al-2012", test_score_table_wee_et_al_2012);

    return g_test_run ();
}
