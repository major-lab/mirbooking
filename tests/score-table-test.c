#include <glib.h>

#include <mirbooking.h>
#include <string.h>
#include <math.h>

#include "../src/mirbooking-score-table-private.h"

#if !GLIB_CHECK_VERSION(2,58,0)
#define g_assert_cmpfloat_with_epsilon(n1,n2,epsilon) g_assert_true(abs((n1) - (n2)) < (epsilon))
#endif

typedef struct __attribute__ ((packed)) _SeedScoreLayout
{
    gsize  n;
    gsize  nnz;
    gsize  rowptr[16384 + 1];
    gsize  colind[11];
    gfloat data[11];
} SeedScoreLayout;

static SeedScoreLayout SEED_SCORES =
{
    16384,
    11,
    {[657] = 0, [658] = 2, [8882] = 2, [8883] = 9, [12152] = 9, [12153] = 10, [13391] = 10, [13392] = 11},
    {11871,  11927, 3228,   4370,  4371,   4372,   7261,   7421,   9284,   13441,  952},
    {-9.79f, -6.30, -2.12f, -0.77, -19.0f, -18.0f, -9.37f, -4.28f, -1.07f, -7.59f, -6.17f},
};

typedef struct __attribute__ ((packed)) _SupplementaryScore3Layout
{
    gsize  n;
    gsize  nnz;
    gsize  rowptr[64 + 1];
    gsize  colind[14];
    gfloat data[14];
} SupplementaryScore3Layout;

static void
test_score_table_compute_seed_score ()
{
    g_autoptr (GBytes) default_table = g_bytes_new_static (&SEED_SCORES, sizeof (SEED_SCORES));
    g_autoptr (MirbookingDefaultScoreTable) score_table = mirbooking_default_score_table_new (default_table,
                                                                                              MIRBOOKING_DEFAULT_SCORE_TABLE_DEFAULT_SUPPLEMENTARY_MODEL,
                                                                                              NULL);

    g_autoptr (MirbookingTarget) target = mirbooking_target_new ("NM");
    g_autoptr (MirbookingMirna) mirna = mirbooking_mirna_new ("MIMAT");

    gchar *target_seq = "GCACACAGAGCAGCATAAAGCCCAGTTGCTTTGGGAAGTGTTTGGGACCAGATGGATTGT";
    gchar *mirna_seq = "UGAGGUAGUAGGUUGUAUAGUU";

    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (target), target_seq);
    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (mirna), mirna_seq);

    g_assert_cmpmem (mirbooking_sequence_get_subsequence (MIRBOOKING_SEQUENCE (target), 1, 7), 7, "CACACAG", 7);
    g_assert_cmpmem (mirbooking_sequence_get_subsequence (MIRBOOKING_SEQUENCE (mirna), 1, 7), 7, "GAGGUAG", 7);

    // test for a seed
    g_autoptr (GError) error = NULL;
    MirbookingScore site_score;
    mirbooking_score_table_compute_score (MIRBOOKING_SCORE_TABLE (score_table),
                                          mirna,     // GAGGUAG =>
                                          target, 1, // CACACAG =>
                                          &site_score,
                                          &error);

    g_assert_cmpfloat_with_epsilon (MIRBOOKING_SCORE_KD (site_score), 1e12 * exp ((T1_ADENOSINE_SCORE + -0.77f + AGO2_SCORE) / (R * T)), 1e-6);
    g_assert_null (error);
}

static void
test_score_table_compute_seed_scores ()
{
    g_autoptr (GBytes) default_table = g_bytes_new_static (&SEED_SCORES, sizeof (SEED_SCORES));
    g_autoptr (GMappedFile) mapped_supplementary_scores = g_mapped_file_new (g_test_get_filename (G_TEST_DIST,  "..",  "data", "scores-3mer", NULL), FALSE, NULL);
    g_autoptr (GBytes) supplementary_scores = g_mapped_file_get_bytes (mapped_supplementary_scores);
    g_autoptr (MirbookingScoreTable) score_table = MIRBOOKING_SCORE_TABLE (mirbooking_default_score_table_new (default_table,
                                                                                                               MIRBOOKING_DEFAULT_SCORE_TABLE_SUPPLEMENTARY_MODEL_YAN_ET_AL_2018,
                                                                                                               supplementary_scores));

    g_autoptr (MirbookingTarget) target = mirbooking_target_new ("NM");
    g_autoptr (MirbookingMirna) mirna = mirbooking_mirna_new ("MIMAT");

    gchar *target_seq = "GCACACAGAGCAGCATAAAGCCCAGTTGCTTTGGGAAGTGTTTGGGACCAGATGGATTGT";
    gchar *mirna_seq = "UGAGGUAGUAGGUUGUAUAGUU";

    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (target), target_seq);
    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (mirna), mirna_seq);

    g_assert_cmpmem (mirbooking_sequence_get_subsequence (MIRBOOKING_SEQUENCE (target), 1, 7), 7, "CACACAG", 7);
    g_assert_cmpmem (mirbooking_sequence_get_subsequence (MIRBOOKING_SEQUENCE (mirna), 1, 7), 7, "GAGGUAG", 7);

    // test for a seed
    g_autoptr (GError) error = NULL;
    g_autofree gsize *positions = NULL;
    gsize  positions_len;
    g_assert (mirbooking_score_table_compute_positions (MIRBOOKING_SCORE_TABLE (score_table),
                                                        mirna,     // GAGGUAG =>
                                                        target, &positions, // CACACAG =>
                                                        &positions_len,
                                                        &error));

    g_assert_nonnull (positions);
    g_assert_cmpint (positions_len, ==, 2);
    MirbookingScore score;
    mirbooking_score_table_compute_score (score_table, mirna, target, positions[0], &score, NULL);
    g_assert_cmpfloat (MIRBOOKING_SCORE_KD (score), ==, 1e12 * exp ((-1.07f + AGO2_SCORE) / (R * T)));
    mirbooking_score_table_compute_score (score_table, mirna, target, positions[1], &score, NULL);
    g_assert_cmpfloat_with_epsilon (MIRBOOKING_SCORE_KD (score), 1e12 * exp ((T1_ADENOSINE_SCORE + -0.77f + AGO2_SCORE) / (R * T)), 1e-6);
    g_assert_null (error);
}

/**
 * Reference:
 * Liang Meng Wee et al., “Argonaute Divides Its RNA Guide into Domains with
 * Distinct Functions and RNA-Binding Properties,” Cell 151, no. 5 (November
 * 21, 2012): 1055–67, https://doi.org/10.1016/j.cell.2012.10.036.
 */
static void
test_score_table_wee_et_al_2012 ()
{
    MirbookingScore score;

    g_autoptr (MirbookingMirna) mirna = mirbooking_mirna_new ("let-7b");
    g_autoptr (MirbookingTarget) target = mirbooking_target_new ("reporter");

    g_autoptr (GBytes) default_table = g_bytes_new_static (&SEED_SCORES, sizeof (SEED_SCORES));
    g_autoptr (GMappedFile) mapped_supplementary_scores = g_mapped_file_new (g_test_get_filename (G_TEST_DIST,  "..",  "data", "scores-3mer", NULL), FALSE, NULL);
    g_autoptr (GBytes) supplementary_scores = g_mapped_file_get_bytes (mapped_supplementary_scores);
    g_autoptr (MirbookingScoreTable) score_table = MIRBOOKING_SCORE_TABLE (mirbooking_default_score_table_new (default_table,
                                                                                                               MIRBOOKING_DEFAULT_SCORE_TABLE_SUPPLEMENTARY_MODEL_YAN_ET_AL_2018,
                                                                                                               supplementary_scores));

    // siRNA
    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (mirna), "UGAGGUAGUAGGUUGUAUAUGU");
    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (target), "GAUACUAUACAACCUACUACCUCAACCU");

    g_autofree gsize *positions = NULL;
    gsize  positions_len;
    mirbooking_score_table_compute_positions (MIRBOOKING_SCORE_TABLE (score_table), mirna, target, &positions, &positions_len, NULL);
    g_assert_cmpint (positions_len, ==, 1);
    g_assert_cmpint (positions[0], ==, 16);

    mirbooking_score_table_compute_score (score_table, mirna, target, 16, &score, NULL);

    g_assert_cmpfloat_with_epsilon (MIRBOOKING_SCORE_KD (score), 1e12 * exp ((T1_ADENOSINE_SCORE - 9.37f + 0.231583f + AGO2_SCORE) / (R * T)), 1e-3);
    // FIXME: g_assert_cmpfloat_with_epsilon (MIRBOOKING_SCORE_KD (score), 20, 10);
    // FIXME: g_assert_cmpfloat_with_epsilon (score.kr, 7.7e-4, 2e-4);
    g_assert_cmpfloat_with_epsilon (MIRBOOKING_SCORE_KM (score), 100, 60);
    // FIXME: g_assert_cmpfloat_with_epsilon (score.kcat, 8.1e-4, 0.1e-4);

    // seed-only
    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (target), "GAAAAAAAAAAAAAAUCUACCUCUAAAU");

    mirbooking_score_table_compute_score (score_table, mirna, target, 16, &score, NULL);

    g_assert_cmpfloat_with_epsilon (MIRBOOKING_SCORE_KD (score), 1e12 * exp ((-9.37f + 0.037555f + AGO2_SCORE) / (R * T)), 1e-3);
    g_assert_cmpfloat (MIRBOOKING_SCORE_KD (score), >=, 26 - 2);
    // FIXME: g_assert_cmpfloat (MIRBOOKING_SCORE_KD (score), <=, 26 + 2);
    g_assert_cmpfloat (score.kcat, <=, 1e-5);
    g_assert_cmpfloat_with_epsilon (MIRBOOKING_SCORE_KM (score), MIRBOOKING_SCORE_KD (score) + (score.kcat / score.kf), 1e-12);

    // seed and supplementary
    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (target), "GAAAAAAAACAAAAAUCUACCUCUAAAU");

    mirbooking_score_table_compute_score (score_table, mirna, target, 16, &score, NULL);

    g_assert_cmpfloat_with_epsilon (MIRBOOKING_SCORE_KD (score), 1e12 * exp ((-9.37f + 0.045937f + AGO2_SCORE) / (R * T)), 1e-3);
    g_assert_cmpfloat (MIRBOOKING_SCORE_KD (score), >=, 13 - 1);
    // FIXME: g_assert_cmpfloat (MIRBOOKING_SCORE_KD (score), <=, 13 + 1);

    g_assert_cmpfloat_with_epsilon (MIRBOOKING_SCORE_KM (score), MIRBOOKING_SCORE_KD (score) + (score.kcat / score.kf), 1e-6);
    g_assert_cmpfloat_with_epsilon (MIRBOOKING_SCORE_KM (score), 100, 60);

    // g10g11 central internal loop 50±30
    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (target), "GAUACUAUACAACGAACUACCUCAACCU");

    mirbooking_score_table_compute_score (score_table, mirna, target, 16, &score, NULL);

    g_assert_cmpfloat_with_epsilon (MIRBOOKING_SCORE_KD (score), 1e12 * exp ((T1_ADENOSINE_SCORE -9.37f + 0.213990f + AGO2_SCORE) / (R * T)), 1e-3);
    // FIXME: g_assert_cmpfloat (MIRBOOKING_SCORE_KD (score), >=, 50 - 30);
    g_assert_cmpfloat (MIRBOOKING_SCORE_KD (score), <=, 50 + 30);

    // g15g16 mismatches in 3' supplementary region 30.83
    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (target), "GAUACUAUCGAACCUACUACCUCAACCU");

    mirbooking_score_table_compute_score (score_table, mirna, target, 16, &score, NULL);

    g_assert_cmpfloat_with_epsilon (MIRBOOKING_SCORE_KD (score), 1e12 * exp ((T1_ADENOSINE_SCORE -9.37f + 0.160218f + AGO2_SCORE) / (R * T)), 1e-3);
    g_assert_cmpfloat_with_epsilon (MIRBOOKING_SCORE_KD (score), 30, 20);

    // g1-g19 complementary 40.83
    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (target), "GAUUAUAUACAACCUACUACCUCAACCU");

    mirbooking_score_table_compute_score (score_table, mirna, target, 16, &score, NULL);

    g_assert_cmpfloat_with_epsilon (MIRBOOKING_SCORE_KD (score), 1e12 * exp ((T1_ADENOSINE_SCORE - 9.37f + 0.231753f + AGO2_SCORE) / (R * T)), 1e-3);
    // FIXME: g_assert_cmpfloat (MIRBOOKING_SCORE_KD (score), >=, 40 - 20);
    g_assert_cmpfloat (MIRBOOKING_SCORE_KD (score), <=, 40 + 20);

    // seed-only
    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (target), "GAAAAAAAAAAAAAAUCUACCUCUAAAU");

    mirbooking_score_table_compute_score (score_table, mirna, target, 16, &score, NULL);

    g_assert_cmpfloat_with_epsilon (MIRBOOKING_SCORE_KD (score), 1e12 * exp ((-9.37f + 0.037555f + AGO2_SCORE) / (R * T)), 1e-3);
    // FIXME: g_assert_cmpfloat_with_epsilon (MIRBOOKING_SCORE_KD (score), 30, 20);

    // seed plus g13-g16 3' supplementary 20±10
    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (target), "GAAAAAAAACAAAAAUCUACCUCUAAAU");

    mirbooking_score_table_compute_score (score_table, mirna, target, 16, &score, NULL);

    g_assert_cmpfloat_with_epsilon (MIRBOOKING_SCORE_KD (score), 1e12 * exp ((-9.37f + 0.045937f + AGO2_SCORE) / (R * T)), 1e-3);
    // FIXME: g_assert_cmpfloat_with_epsilon (MIRBOOKING_SCORE_KD (score), 20, 10);

    // g4g5 mismatches in seed 1e3±0.6e3
    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (target), "GAUACUAUACAACCUACUAUUUCAACCU");

    mirbooking_score_table_compute_score (score_table, mirna, target, 16, &score, NULL);

    g_assert_cmpfloat_with_epsilon (MIRBOOKING_SCORE_KD (score), 1e12 * exp ((T1_ADENOSINE_SCORE - 4.28f + 0.231583f + AGO2_SCORE) / (R * T)), 1e-1);
    g_assert_cmpfloat (MIRBOOKING_SCORE_KD (score), >=, 1e3 - 0.6e3);
    // FIXME: g_assert_cmpfloat (MIRBOOKING_SCORE_KD (score), <=, 1e3 + 0.6e3);

    // non-complementary luciferase target 2e3±1e3
    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (target), "GAUAAGGCAUUUCAUUAUAGCUAUACCU");

    mirbooking_score_table_compute_score (score_table, mirna, target, 16, &score, NULL);

    g_assert_cmpfloat_with_epsilon (MIRBOOKING_SCORE_KD (score), 1e12 * exp ((-2.12f + 0.052757f + AGO2_SCORE) / (R * T)), 10);
    g_assert_cmpfloat (MIRBOOKING_SCORE_KD (score), >=, 2e3 - 1e3);
    // FIXME: g_assert_cmpfloat (MIRBOOKING_SCORE_KD (score), <=, 2e3 + 1e3);
}

/**
 * Reference:
 * William E. Salomon et al., “Single-Molecule Imaging Reveals That Argonaute
 * Reshapes the Binding Properties of Its Nucleic Acid Guides,” Cell 162, no. 1
 * (July 2, 2015): 84–95, https://doi.org/10.1016/j.cell.2015.06.029.
 */
static void
test_score_table_salomon_et_al_2016 ()
{
    MirbookingScore score;

    g_autoptr (MirbookingMirna) mirna = mirbooking_mirna_new ("let-7b");
    g_autoptr (MirbookingTarget) target = mirbooking_target_new ("reporter");

    g_autoptr (GMappedFile) mapped_seed_scores = g_mapped_file_new (g_test_get_filename (G_TEST_DIST,  "..",  "data", "scores-7mer-3mismatch-ending", NULL), FALSE, NULL);
    g_autoptr (GBytes) default_table = g_mapped_file_get_bytes (mapped_seed_scores);

    g_autoptr (GMappedFile) mapped_supplementary_scores = g_mapped_file_new (g_test_get_filename (G_TEST_DIST,  "..",  "data", "scores-3mer", NULL), FALSE, NULL);
    g_autoptr (GBytes) supplementary_scores = g_mapped_file_get_bytes (mapped_supplementary_scores);
    g_autoptr (MirbookingScoreTable) score_table = MIRBOOKING_SCORE_TABLE (mirbooking_default_score_table_new (default_table,
                                                                                                               MIRBOOKING_DEFAULT_SCORE_TABLE_SUPPLEMENTARY_MODEL_YAN_ET_AL_2018,
                                                                                                               supplementary_scores));

    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (mirna),  "UGAGGUAGUAGGUUGUAUAGU");

    // Complete
    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (target), "ACUAUACAACCUACUACCUCA");
    mirbooking_score_table_compute_score (score_table, mirna, target, 13, &score, NULL);

    g_assert_cmpfloat_with_epsilon (MIRBOOKING_SCORE_KD (score), 1e12 * exp ((T1_ADENOSINE_SCORE - 9.15f + 0.236105f + AGO2_SCORE) / (R * T)), 1e-3);
    // dye control (Figure 3B)
    // FIXME: g_assert_cmpfloat_with_epsilon (MIRBOOKING_SCORE_KM (score), 1.7e3, 0.1e3);
    g_assert_cmpfloat_with_epsilon (score.kf, 3.9e-4, 0.5e-2);
    // FIXME: g_assert_cmpfloat_with_epsilon (score.kcat, 3.6e-2, 0.2e-2);

    // Seed plus 3'UTR
    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (target), "UGAUAACAAGGAUCUACCUCA");
    mirbooking_score_table_compute_score (score_table, mirna, target, 13, &score, NULL);
    g_assert_cmpfloat_with_epsilon (MIRBOOKING_SCORE_KD (score), 1e12 * exp ((T1_ADENOSINE_SCORE - 9.15f + 0.045951f + AGO2_SCORE) / (R * T)), 1e-3);
    // FIXME: g_assert_cmpfloat_with_epsilon (MIRBOOKING_SCORE_KD (score), 11, 2);
    g_assert_cmpfloat_with_epsilon (score.kf, 2.8e-4, 0.5e-2);
    // FIXME: g_assert_cmpfloat (score.kr, >=, 3e-3 - 0.4e-4);
    g_assert_cmpfloat (score.kr, <=, 3e-3 + 0.4e-3);

    // AAAAA...AAA
    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (target), "ACUAUACAACCAAAAAAAAAA");
    mirbooking_score_table_compute_score (score_table, mirna, target, 13, &score, NULL);
    // FIXME: g_assert_cmpfloat (score.kf, >=, 3.6e-5 - 0.2e-5);
    g_assert_cmpfloat (score.kf, <=, 3.6e-5 + 0.2e-5);
    g_assert_cmpfloat (score.kr, ==, INFINITY);

    // Seed only
    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (target), "UGAUAUGUUGGAUCUACCUCA");
    mirbooking_score_table_compute_score (score_table, mirna, target, 13, &score, NULL);
    // FIXME: g_assert_cmpfloat_with_epsilon (score.kf, 2.4e-4, 0.1e-4);
    // FIXME: g_assert_cmpfloat (score.kr, >=, 3.6e-3 - 0.3e-3);
    g_assert_cmpfloat (score.kr, <=, 3.6e-3 + 0.3e-3);
    g_assert_cmpfloat_with_epsilon (MIRBOOKING_SCORE_KD (score), 1e12 * exp ((T1_ADENOSINE_SCORE - 9.15f + 0.012648 + AGO2_SCORE) / (R * T)), 1e-3);
    // FIXME: g_assert_cmpfloat (MIRBOOKING_SCORE_KD (score), <=, 15 + 2);
    // FIXME: g_assert_cmpfloat (MIRBOOKING_SCORE_KD (score), >=, 15 - 2);

    // g2g3
    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (target), "UGAUAUGUUGGAUCUACCAGA");
    mirbooking_score_table_compute_score (score_table, mirna, target, 13, &score, NULL);
    // FIXME: g_assert_cmpfloat (score.kf, >=, 3.8e-5 - 0.1e-5);
    g_assert_cmpfloat (score.kf, <=, 3.8e-5 + 0.1e-5);
    g_assert_cmpfloat (score.kr, ==, INFINITY);

    // g3g4
    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (target), "UGAUAUGUUGGAUCUACGACA");
    mirbooking_score_table_compute_score (score_table, mirna, target, 13, &score, NULL);
    // FIXME: g_assert_cmpfloat (score.kf, >=, 3.4e-5 - 0.2e-5);
    g_assert_cmpfloat (score.kf, <=, 3.4e-5 + 0.2e-5);
    g_assert_cmpfloat (score.kr, ==, INFINITY);

    // g4g5
    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (target), "UGAUAUGUUGGAUCUAAAUCA");
    mirbooking_score_table_compute_score (score_table, mirna, target, 13, &score, NULL);
    // FIXME: g_assert_cmpfloat (score.kf, >=, 2.4e-5 - 0.1e-5);
    g_assert_cmpfloat (score.kf, <=, 2.4e-5 + 0.1e-5);
    g_assert_cmpfloat (score.kr, ==, INFINITY);

    // g5g6
    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (target), "UGAUAUGUUGGAUCUUGCUCA");
    mirbooking_score_table_compute_score (score_table, mirna, target, 13, &score, NULL);
    // FIXME: g_assert_cmpfloat (score.kf, >=, 3.9e-5 - 0.1e-5);
    g_assert_cmpfloat (score.kf, <=, 3.9e-5 + 0.1e-5);
    g_assert_cmpfloat (score.kr, ==, INFINITY);

    // g6g7
    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (target), "UGAUAUGUUGGAUCAUCCUCA");
    mirbooking_score_table_compute_score (score_table, mirna, target, 13, &score, NULL);
    // FIXME: g_assert_cmpfloat_with_epsilon (score.kf, 2.0e-4, 0.1e-4);
    // FIXME: g_assert_cmpfloat (score.kr, >=, 3.5 - 0.1);
    g_assert_cmpfloat (score.kr, <=, 3.5 + 0.1);
    g_assert_cmpfloat (MIRBOOKING_SCORE_KD (score), >=, 18e3 - 2e3);
    // FIXME: g_assert_cmpfloat (MIRBOOKING_SCORE_KD (score), <=, 18e3 + 2e3);

    // g7g8
    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (target), "UGAUAUGUUGGAUGAACCUCA");
    mirbooking_score_table_compute_score (score_table, mirna, target, 13, &score, NULL);
    // FIXME: g_assert_cmpfloat_with_epsilon (score.kf, 2.0e-4, 0.1e-4);
    // FIXME: g_assert_cmpfloat (score.kr, >=, 0.24 - 0.01);
    g_assert_cmpfloat (score.kr, <=, 0.24 + 0.01);
    // FIXME: g_assert_cmpfloat_with_epsilon (MIRBOOKING_SCORE_KD (score), 1.2e3, 0.2e3);

    // g8
    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (target), "UGAUAUGUUGGAUGUACCUCA");
    mirbooking_score_table_compute_score (score_table, mirna, target, 13, &score, NULL);
    // FIXME: g_assert_cmpfloat (score.kf, >=, 1.8e-4 - 0.1e-4);
    g_assert_cmpfloat (score.kf, <=, 1.8e-4 + 0.1e-4);
    // FIXME: g_assert_cmpfloat (score.kr, >=, 0.086 - 0.002);
    g_assert_cmpfloat (score.kr, <=, 0.086 + 0.002);
    // FIXME: g_assert_cmpfloat (MIRBOOKING_SCORE_KD (score), >=, 0.48e3 - 0.06e3);
    g_assert_cmpfloat (MIRBOOKING_SCORE_KD (score), <=, 0.48e3 + 0.06e3);
}

static void
test_score_table_schirle_et_al_2015 ()
{
    MirbookingScore score;
    g_autoptr (MirbookingTarget) target = mirbooking_target_new ("");

    g_autoptr (MirbookingMirna) mirna = mirbooking_mirna_new ("");
    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (mirna), "UUCACAUUGCCCAAGUCUCUU");

    g_autoptr (GMappedFile) mapped_seed_scores = g_mapped_file_new (g_test_get_filename (G_TEST_DIST,  "..",  "data", "scores-7mer-3mismatch-ending", NULL), FALSE, NULL);
    g_autoptr (GBytes) default_table = g_mapped_file_get_bytes (mapped_seed_scores);
    g_autoptr (GMappedFile) mapped_supplementary_scores = g_mapped_file_new (g_test_get_filename (G_TEST_DIST,  "..",  "data", "scores-3mer", NULL), FALSE, NULL);
    g_autoptr (GBytes) supplementary_scores = g_mapped_file_get_bytes (mapped_supplementary_scores);
    g_autoptr (MirbookingScoreTable) score_table = MIRBOOKING_SCORE_TABLE (mirbooking_default_score_table_new (default_table,
                                                                                                               MIRBOOKING_DEFAULT_SCORE_TABLE_SUPPLEMENTARY_MODEL_YAN_ET_AL_2018,
                                                                                                               supplementary_scores));

    // A
    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (target), "CAAUGUGAAAA");
    mirbooking_score_table_compute_score (score_table, mirna, target, 1, &score, NULL);
    g_assert_cmpfloat (MIRBOOKING_SCORE_KD (score), ==, 1e12 * exp ((T1_ADENOSINE_SCORE - 6.17f + AGO2_SCORE) / (R * T)));
    // FIXME: g_assert_cmpfloat (MIRBOOKING_SCORE_KD (score), <=, 0.75e3 + 0.04e3);
    g_assert_cmpfloat (MIRBOOKING_SCORE_KD (score), >=, 0.75e3 - 0.04e3);

    // U
    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (target), "CAAUGUGAUAA");
    mirbooking_score_table_compute_score (score_table, mirna, target, 1, &score, NULL);
    g_assert_cmpfloat (MIRBOOKING_SCORE_KD (score), ==, 1e12 * exp ((-6.17f + AGO2_SCORE) / (R * T)));
    // FIXME: g_assert_cmpfloat (MIRBOOKING_SCORE_KD (score), <=, 1.9e3 + 0.09e3);
    g_assert_cmpfloat (MIRBOOKING_SCORE_KD (score), >=, 1.9e3 - 0.09e3);

    // C
    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (target), "CAAUGUGACAA");
    mirbooking_score_table_compute_score (score_table, mirna, target, 1, &score, NULL);
    g_assert_cmpfloat (MIRBOOKING_SCORE_KD (score), ==, 1e12 * exp ((-6.17f + AGO2_SCORE) / (R * T)));
    // FIXME: g_assert_cmpfloat (MIRBOOKING_SCORE_KD (score), <=, 1.9e3 + 0.10e3);
    g_assert_cmpfloat (MIRBOOKING_SCORE_KD (score), >=, 1.9e3 - 0.10e3);

    // G
    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (target), "CAAUGUGAGAA");
    mirbooking_score_table_compute_score (score_table, mirna, target, 1, &score, NULL);
    g_assert_cmpfloat (MIRBOOKING_SCORE_KD (score), ==, 1e12 * exp ((-6.17f + AGO2_SCORE) / (R * T)));
    // FIXME: g_assert_cmpfloat (MIRBOOKING_SCORE_KD (score), <=, 1.8e3 + 0.12e3);
    g_assert_cmpfloat (MIRBOOKING_SCORE_KD (score), >=, 1.8e3 - 0.12e3);
}

static void
test_score_table_yan_et_al_2018 ()
{
    /*
     * The paper does not contains measurements, but we have to ensure that
     * we're capturing the correct boxes.
     */

    g_autoptr (MirbookingMirna) mirna = mirbooking_mirna_new ("miB");
    g_autoptr (MirbookingTarget) target = mirbooking_target_new ("FR-tat");

    g_autoptr (GMappedFile) mapped_seed_scores = g_mapped_file_new (g_test_get_filename (G_TEST_DIST,  "..",  "data", "scores-7mer-3mismatch-ending", NULL), FALSE, NULL);
    g_autoptr (GBytes) default_table = g_mapped_file_get_bytes (mapped_seed_scores);
    g_autoptr (GMappedFile) mapped_supplementary_scores = g_mapped_file_new (g_test_get_filename (G_TEST_DIST,  "..",  "data", "scores-3mer", NULL), FALSE, NULL);
    g_autoptr (GBytes) supplementary_scores = g_mapped_file_get_bytes (mapped_supplementary_scores);
    g_autoptr (MirbookingScoreTable) score_table = MIRBOOKING_SCORE_TABLE (mirbooking_default_score_table_new (default_table,
                                                                                                               MIRBOOKING_DEFAULT_SCORE_TABLE_SUPPLEMENTARY_MODEL_YAN_ET_AL_2018,
                                                                                                               supplementary_scores));

    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (target), "GACGAAGAGCUCAUCAGAACA");

    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (mirna),  "UGUUCUGAUGAGCUCUUCGUC");

    g_assert_cmpmem ("GUUCUGA", 7, mirbooking_sequence_get_subsequence (MIRBOOKING_SEQUENCE (mirna), 1, 7), 7);
    g_assert_cmpmem ("UGA", 3, mirbooking_sequence_get_subsequence     (MIRBOOKING_SEQUENCE (mirna), 8, 3), 3);
    g_assert_cmpmem ("GCU",3, mirbooking_sequence_get_subsequence      (MIRBOOKING_SEQUENCE (mirna), 11, 3), 3);
    g_assert_cmpmem ("CUU", 3, mirbooking_sequence_get_subsequence     (MIRBOOKING_SEQUENCE (mirna), 14, 3), 3);
    g_assert_cmpmem ("CGUC", 4, mirbooking_sequence_get_subsequence    (MIRBOOKING_SEQUENCE (mirna), 17, 4), 4);

    MirbookingScore score;

    mirbooking_score_table_compute_score (score_table,
                                               mirna,
                                               target,
                                               13,
                                               &score,
                                               NULL);

    g_assert_cmpfloat_with_epsilon (MIRBOOKING_SCORE_KD (score), 1e12 * exp ((T1_ADENOSINE_SCORE -7.59f - 0.602573f + AGO2_SCORE) / (R * T)), 1e-3);

    // A+C
    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (mirna),  "UGUUCUGAACUGCUGAACGUC");
    mirbooking_score_table_compute_score (score_table,
                                               mirna,
                                               target,
                                               13,
                                               &score,
                                               NULL);

    g_assert_cmpfloat_with_epsilon (MIRBOOKING_SCORE_KD (score), 1e12 * exp ((T1_ADENOSINE_SCORE -7.59f - 0.869875f + AGO2_SCORE) / (R * T)), 1e-3);

    // B+D
    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (mirna),  "UGUUCUGAUGACGACUUCGUC");
    mirbooking_score_table_compute_score (score_table,
                                               mirna,
                                               target,
                                               13,
                                               &score,
                                               NULL);

    g_assert_cmpfloat_with_epsilon (MIRBOOKING_SCORE_KD (score), 1e12 * exp ((T1_ADENOSINE_SCORE -7.59f + 0.022712f + AGO2_SCORE) / (R * T)), 1e-3);

    // A+D
    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (mirna),  "UGUUCUGAACUGCUCUUGCAG");
    mirbooking_score_table_compute_score (score_table,
                                               mirna,
                                               target,
                                               13,
                                               &score,
                                               NULL);

    g_assert_cmpfloat_with_epsilon (MIRBOOKING_SCORE_KD (score), 1e12 * exp ((T1_ADENOSINE_SCORE -7.59f - 0.765957f + AGO2_SCORE) / (R * T)), 1e-3);

    // A+B
    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (mirna),  "UGUUCUGAACUGCACUUCGUC");
    mirbooking_score_table_compute_score (score_table,
                                               mirna,
                                               target,
                                               13,
                                               &score,
                                               NULL);

    g_assert_cmpfloat_with_epsilon (MIRBOOKING_SCORE_KD (score), 1e12 * exp ((T1_ADENOSINE_SCORE -7.59f - 0.661824f + AGO2_SCORE) / (R * T)), 1e-3);

    // B+C
    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (mirna),  "UGUUCUGAUGACGAGAACGUC");
    mirbooking_score_table_compute_score (score_table,
                                               mirna,
                                               target,
                                               13,
                                               &score,
                                               NULL);

    g_assert_cmpfloat_with_epsilon (MIRBOOKING_SCORE_KD (score), 1e12 * exp ((T1_ADENOSINE_SCORE - 7.59f + 0.014053f + AGO2_SCORE) / (R * T)), 1e-3);

    // C+D
    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (mirna),  "UGUUCUGAUGAGCUGAAGCAG");
    mirbooking_score_table_compute_score (score_table,
                                               mirna,
                                               target,
                                               13,
                                               &score,
                                               NULL);

    g_assert_cmpfloat_with_epsilon (MIRBOOKING_SCORE_KD (score), 1e12 * exp ((T1_ADENOSINE_SCORE -7.59f - 0.862300f + AGO2_SCORE) / (R * T)), 1e-3);
}

/**
 * Reference:
 * Myung Hyun Jo et al., “Human Argonaute 2 Has Diverse Reaction Pathways on
 * Target RNAs,” Molecular Cell 59, no. 1 (July 2, 2015): 117–24,
 * https://doi.org/10.1016/j.molcel.2015.04.027.
 */
static void
test_score_table_jo_et_al_2015 ()
{
    g_autoptr (GMappedFile) mapped_seed_scores = g_mapped_file_new (g_test_get_filename (G_TEST_DIST,  "..",  "data", "scores-7mer-3mismatch-ending", NULL), FALSE, NULL);
    g_autoptr (GBytes) default_table = g_mapped_file_get_bytes (mapped_seed_scores);
    g_autoptr (GMappedFile) mapped_supplementary_scores = g_mapped_file_new (g_test_get_filename (G_TEST_DIST,  "..",  "data", "scores-3mer", NULL), FALSE, NULL);
    g_autoptr (GBytes) supplementary_scores = g_mapped_file_get_bytes (mapped_supplementary_scores);
    g_autoptr (MirbookingScoreTable) score_table = MIRBOOKING_SCORE_TABLE (mirbooking_default_score_table_new (default_table,
                                                                                                                          MIRBOOKING_DEFAULT_SCORE_TABLE_SUPPLEMENTARY_MODEL_YAN_ET_AL_2018,
                                                                                                                          supplementary_scores));

    g_autoptr (MirbookingMirna) mirna = mirbooking_mirna_new ("");
    g_autoptr (MirbookingTarget) target = mirbooking_target_new ("");

    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (mirna), "UGAGGUAGUAGGUUGUAUAGU");
    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (target), "ACUAUACAACCUACUACCUCG");

    MirbookingScore score;
    mirbooking_score_table_compute_score (score_table, mirna, target, 13, &score, NULL);

    // FIXME: g_assert_cmpfloat (score.kf, ==, 1.7e-5);
}

static void
test_score_table_chi_et_al_2012 ()
{
    g_autoptr (GBytes) default_table = g_bytes_new_static (&SEED_SCORES, sizeof (SEED_SCORES));
    g_autoptr (GMappedFile) mapped_supplementary_scores = g_mapped_file_new (g_test_get_filename (G_TEST_DIST,  "..",  "data", "scores-3mer", NULL), FALSE, NULL);
    g_autoptr (GBytes) supplementary_scores = g_mapped_file_get_bytes (mapped_supplementary_scores);
    g_autoptr (MirbookingScoreTable) score_table = MIRBOOKING_SCORE_TABLE (mirbooking_default_score_table_new (default_table,
                                                                                                                          MIRBOOKING_DEFAULT_SCORE_TABLE_SUPPLEMENTARY_MODEL_YAN_ET_AL_2018,
                                                                                                                          supplementary_scores));

    g_autoptr (MirbookingTarget) target = mirbooking_target_new ("");
    g_autoptr (MirbookingMirna) mirna = mirbooking_mirna_new ("");

    // G-buldged seed
    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (target), "GUGGCCUU");
    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (mirna),  "UAAGGCAC");

    MirbookingScore score;
    mirbooking_score_table_compute_score (score_table, mirna, target, 0, &score, NULL);

    gfloat G_seed = (exp (6.3f/(R*T)) * (-6.3f) + exp (8.59f/(R*T)) * (-8.59f)) / (exp (6.3f/(R*T)) + exp (8.59f/(R*T)));
    g_assert_cmpfloat_with_epsilon (MIRBOOKING_SCORE_KD (score), 1e12 * exp ((G_seed + AGO2_SCORE) / (R * T)), 1e-12);
}

int main (int argc, gchar **argv)
{
    g_test_init (&argc, &argv, NULL);

    g_test_add_func ("/score-table/compute-seed-score", test_score_table_compute_seed_score);
    g_test_add_func ("/score-table/compute-seed-scores", test_score_table_compute_seed_scores);
    g_test_add_func ("/score-table/wee-et-al-2012", test_score_table_wee_et_al_2012);
    g_test_add_func ("/score-table/salomon-et-al-2016", test_score_table_salomon_et_al_2016);
    g_test_add_func ("/score-table/schirle-et-al-2015", test_score_table_schirle_et_al_2015);
    g_test_add_func ("/score-table/yan-et-al-2018", test_score_table_yan_et_al_2018);
    g_test_add_func ("/score-table/jo-et-al-2015", test_score_table_jo_et_al_2015);
    g_test_add_func ("/score-table/chi-et-al-2012", test_score_table_chi_et_al_2012);

    return g_test_run ();
}
