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
    gsize  colind[7];
    gfloat data[7];
} SeedScoreLayout;

static SeedScoreLayout SEED_SCORES =
{
    16384,
    7,
    {[8882] = 0, [8883] = 7},
    {3228,  4370,   4371,   4372,   7261,   7421,   9284},
    {-1.6f, -20.0f, -19.0f, -18.0f, -9.30f, -4.20f, -19.0f},
};

typedef struct __attribute__ ((packed)) _SupplementaryScoreLayout
{
    gsize  n;
    gsize  nnz;
    gsize  rowptr[256 + 1];
    gsize  colind[5];
    gfloat data[5];
} SupplementaryScoreLayout;

static SupplementaryScoreLayout SUPPLEMENTARY_SCORES =
{
    256,
    5,
    {[251] = 0, [252] = 5},
    {0,    16,     63,   96,   239},
    {0.0f, -0.20f, 0.0f, 0.0f, 0.0f},
};

static void
test_score_table_compute_seed_score ()
{
    g_autoptr (GBytes) default_table = g_bytes_new_static (&SEED_SCORES, sizeof (SEED_SCORES));
    g_autoptr (MirbookingDefaultScoreTable) score_table = mirbooking_default_score_table_new_from_bytes (default_table, 1, 7, NULL);

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

    g_assert_cmpfloat (site_score, ==, 1e12 * exp ((-20.0f - 6.06f) / (R * T)));
    g_assert_null (error);
}

static void
test_score_table_compute_seed_scores ()
{
    g_autoptr (GBytes) default_table = g_bytes_new_static (&SEED_SCORES, sizeof (SEED_SCORES));
    g_autoptr (MirbookingDefaultScoreTable) score_table = mirbooking_default_score_table_new_from_bytes (default_table, 1, 7, NULL);

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
    g_assert_cmpfloat (mirbooking_score_table_compute_score (score_table, mirna, target, positions[0], NULL), ==, 1e12 * exp ((-19.0f - 6.06f) / (R * T)));
    g_assert_cmpfloat (mirbooking_score_table_compute_score (score_table, mirna, target, positions[1], NULL), ==, 1e12 * exp ((-20.0f - 6.06f) / (R * T)));
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
    g_assert_cmpfloat (site_score, ==, 1e12 * exp ((-7.09f + 3.02f) / (R * T)));
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
    gdouble Kd, Km;

    g_autoptr (MirbookingMirna) mirna = mirbooking_mirna_new ("let-7b");
    g_autoptr (MirbookingTarget) target = mirbooking_target_new ("reporter");

    g_autoptr (GBytes) default_table = g_bytes_new_static (&SEED_SCORES, sizeof (SEED_SCORES));
    g_autoptr (GBytes) supplementary_scores = g_bytes_new_static (&SUPPLEMENTARY_SCORES, sizeof (SUPPLEMENTARY_SCORES));
    g_autoptr (MirbookingScoreTable) score_table = MIRBOOKING_SCORE_TABLE (mirbooking_default_score_table_new_from_bytes (default_table, 1, 7, supplementary_scores));

    // siRNA
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (mirna), "UGAGGUAGUAGGUUGUAUAUGU", strlen ("UGAGGUAGUAGGUUGUAUAUGU"));
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (target), "GAUACUAUACAACCUACUACCUCAACCU", strlen ("GAUACUAUACAACCUACUACCUCAACCU"));

    gsize *positions;
    gsize  positions_len;
    mirbooking_score_table_compute_positions (MIRBOOKING_SCORE_TABLE (score_table), mirna, target, &positions, &positions_len, NULL);
    g_assert_cmpint (positions_len, ==, 1);
    g_assert_cmpint (positions[0], ==, 16);

    Kd = mirbooking_score_table_compute_score (score_table, mirna, target, 16, NULL);
    Km = mirbooking_score_table_compute_enzymatic_score (score_table, mirna, target, 16, NULL);

    g_assert_cmpfloat (Kd, ==, 1e12 * exp ((-9.30f - 0.20f - 6.06f) / (R * T)));
    g_assert_cmpfloat (Kd, >=, 20 - 10);
    g_assert_cmpfloat (Kd, <=, 20 + 10);

    g_assert_cmpfloat (Km, ==, Kd + (MIRBOOKING_SCORE_TABLE_DEFAULT_KCAT / MIRBOOKING_SCORE_TABLE_DEFAULT_KF));
    g_assert_cmpfloat (Km, >=, 100 - 60);
    // FIXME: g_assert_cmpfloat (Km, <=, 100 + 60);

    // seed-only
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (target), "GAAAAAAAAAAAAAAUCUACCUCUAAAU", strlen ("GAAAAAAAAAAAAAAUCUACCUCUAAAU"));

    Kd = mirbooking_score_table_compute_score (score_table, mirna, target, 16, NULL);
    Km = mirbooking_score_table_compute_enzymatic_score (score_table, mirna, target, 16, NULL);

    g_assert_cmpfloat (Kd, ==, 1e12 * exp ((-9.30f - 6.06f) / (R * T)));
    // FIXME: g_assert_cmpfloat (Kd, >=, 26 - 2);
    g_assert_cmpfloat (Kd, <=, 26 + 2);

    g_assert_cmpfloat (Km, ==, Kd + (MIRBOOKING_SCORE_TABLE_DEFAULT_KCAT / MIRBOOKING_SCORE_TABLE_DEFAULT_KF));
    g_assert_cmpfloat (Km, >=, 100 - 60);
    // FIXME: g_assert_cmpfloat (Km, <=, 100 + 60);

    // seed and supplementary
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (target), "GAAAAAAAACAAAAAUCUACCUCUAAAU", strlen ("GAAAAAAAAAAAAAAUCUACCUCUAAAU"));

    Kd = mirbooking_score_table_compute_score (score_table, mirna, target, 16, NULL);
    Km = mirbooking_score_table_compute_enzymatic_score (score_table, mirna, target, 16, NULL);

    g_assert_cmpfloat (Kd, ==, 1e12 * exp ((-9.30f - 0.20f - 6.06f) / (R * T)));
    // FIXME: g_assert_cmpfloat (Kd, >=, 13 - 1);
    g_assert_cmpfloat (Kd, <=, 13 + 1);

    g_assert_cmpfloat (Km, ==, Kd + (MIRBOOKING_SCORE_TABLE_DEFAULT_KCAT / MIRBOOKING_SCORE_TABLE_DEFAULT_KF));
    g_assert_cmpfloat (Km, >=, 100 - 60);
    // FIXME: g_assert_cmpfloat (Km, <=, 100 + 60);

    // g10g11 central internal loop 50±30
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (target), "GAUACUAUACAACGAACUACCUCAACCU", strlen ("GAUACUAUACAACGAACUACCUCAACCU"));

    Kd = mirbooking_score_table_compute_score (score_table, mirna, target, 16, NULL);
    Km = mirbooking_score_table_compute_enzymatic_score (score_table, mirna, target, 16, NULL);

    g_assert_cmpfloat (Kd, ==, 1e12 * exp ((-9.30f - 0.20f - 6.06f) / (R * T)));
    // FIXME: g_assert_cmpfloat (Kd, >=, 50 - 30);
    g_assert_cmpfloat (Kd, <=, 50 + 30);

    // g15g16 mismatches in 3' supplementary region 30±20
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (target), "GAUACUAUCGAACCUACUACCUCAACCU", strlen ("GAUACUAUCGAACCUACUACCUCAACCU"));

    Kd = mirbooking_score_table_compute_score (score_table, mirna, target, 16, NULL);
    Km = mirbooking_score_table_compute_enzymatic_score (score_table, mirna, target, 16, NULL);

    g_assert_cmpfloat (Kd, ==, 1e12 * exp ((-9.30f - 6.06f) / (R * T)));
    g_assert_cmpfloat (Kd, >=, 30 - 20);
    g_assert_cmpfloat (Kd, <=, 30 + 20);

    // g1-g19 complementary 40±20
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (target), "GAUUAUAUACAACCUACUACCUCAACCU", strlen ("GAUUAUAUACAACCUACUACCUCAACCU"));

    Kd = mirbooking_score_table_compute_score (score_table, mirna, target, 16, NULL);
    Km = mirbooking_score_table_compute_enzymatic_score (score_table, mirna, target, 16, NULL);

    g_assert_cmpfloat (Kd, ==, 1e12 * exp ((-9.30f - 0.20f - 6.06f) / (R * T)));
    // FIXME: g_assert_cmpfloat (Kd, >=, 40 - 20);
    g_assert_cmpfloat (Kd, <=, 40 + 20);

    // seed-only 30±20
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (target), "GAAAAAAAAAAAAAAUCUACCUCUAAAU", strlen ("GAAAAAAAAAAAAAAUCUACCUCUAAAU"));

    Kd = mirbooking_score_table_compute_score (score_table, mirna, target, 16, NULL);
    Km = mirbooking_score_table_compute_enzymatic_score (score_table, mirna, target, 16, NULL);

    g_assert_cmpfloat (Kd, ==, 1e12 * exp ((-9.30f - 6.06f) / (R * T)));
    g_assert_cmpfloat (Kd, >=, 30 - 20);
    g_assert_cmpfloat (Kd, <=, 30 + 20);

    // seed plus g13-g16 3' supplementary 20±10
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (target), "GAAAAAAAACAAAAAUCUACCUCUAAAU", strlen ("GAAAAAAAACAAAAAUCUACCUCUAAAU"));

    Kd = mirbooking_score_table_compute_score (score_table, mirna, target, 16, NULL);
    Km = mirbooking_score_table_compute_enzymatic_score (score_table, mirna, target, 16, NULL);

    g_assert_cmpfloat (Kd, ==, 1e12 * exp ((-9.30f - 0.20f - 6.06f) / (R * T)));
    g_assert_cmpfloat (Kd, >=, 20 - 10);
    g_assert_cmpfloat (Kd, <=, 20 + 10);

    // g4g5 mismatches in seed 1e3±0.6e3
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (target), "GAUACUAUACAACCUACUAUUUCAACCU", strlen ("GAUACUAUACAACCUACUAUUUCAACCU"));

    Kd = mirbooking_score_table_compute_score (score_table, mirna, target, 16, NULL);
    Km = mirbooking_score_table_compute_enzymatic_score (score_table, mirna, target, 16, NULL);

    g_assert_cmpfloat (Kd, ==, 1e12 * exp ((-4.20f - 0.20f - 6.06f) / (R * T)));
    g_assert_cmpfloat (Kd, >=, 1e3 - 0.6e3);
    // FIXME: g_assert_cmpfloat (Kd, <=, 1e3 + 0.6e3);

    // non-complementary luciferase target 2e3±1e3
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (target), "GAUAAGGCAUUUCAUUAUAGCUAUACCU", strlen ("GAUAAGGCAUUUCAUUAUAGCUAUACCU"));

    Kd = mirbooking_score_table_compute_score (score_table, mirna, target, 16, NULL);
    Km = mirbooking_score_table_compute_enzymatic_score (score_table, mirna, target, 16, NULL);

    g_assert_cmpfloat (Kd, ==, 1e12 * exp ((-1.6f - 6.06f) / (R * T)));
    g_assert_cmpfloat (Kd, >=, 2e3 - 1e3);
    // FIXME: g_assert_cmpfloat (Kd, <=, 2e3 + 1e3);

    g_free (positions);
}

/**
 * Reference:
 * William E. Salomon et al., “Single-Molecule Imaging Reveals That Argonaute
 * Reshapes the Binding Properties of Its Nucleic Acid Guides,” Cell 162, no. 1
 * (July 2, 2015): 84–95, https://doi.org/10.1016/j.cell.2015.06.029.
 */
static void
test_score_table_solomon_et_al_2016 ()
{
    gdouble Kd;

    g_autoptr (MirbookingMirna) mirna = mirbooking_mirna_new ("let-7b");
    g_autoptr (MirbookingTarget) target = mirbooking_target_new ("reporter");

    g_autoptr (GBytes) default_table = g_bytes_new_static (&SEED_SCORES, sizeof (SEED_SCORES));
    g_autoptr (GBytes) supplementary_scores = g_bytes_new_static (&SUPPLEMENTARY_SCORES, sizeof (SUPPLEMENTARY_SCORES));
    g_autoptr (MirbookingScoreTable) score_table = MIRBOOKING_SCORE_TABLE (mirbooking_default_score_table_new_from_bytes (default_table, 1, 7, supplementary_scores));

    // Complete
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (mirna),  "UGAGGUAGUAGGUUGUAUAGU", strlen ("UGAGGUAGUAGGUUGUAUAGU"));
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (target), "ACUAUACAACCUACUACCUCA", strlen ("ACUAUACAACCUACUACCUCA"));

    Kd = mirbooking_score_table_compute_score (score_table, mirna, target, 13, NULL);

    // Seed only
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (target), "UGAUAUGUUGGAUCUACCUCA", strlen ("ACUAUACAACCUACUACCUCA"));

    g_assert_cmpfloat (MIRBOOKING_SCORE_TABLE_DEFAULT_KF, >=, 2.4e-4 - 0.1e-4);
    g_assert_cmpfloat (MIRBOOKING_SCORE_TABLE_DEFAULT_KF, <=, 2.4e-4 + 0.1e-4);
    g_assert_cmpfloat (MIRBOOKING_SCORE_TABLE_DEFAULT_KCAT, >=, 3.6e-2 - 0.2e-2);
    g_assert_cmpfloat (MIRBOOKING_SCORE_TABLE_DEFAULT_KCAT, <=, 3.6e-2 + 0.2e-2);

    Kd = mirbooking_score_table_compute_score (score_table, mirna, target, 13, NULL);

    g_assert_cmpfloat (Kd, >=, 15 - 2);
    g_assert_cmpfloat (Kd, <=, 15 + 2);

    // Seed plus 3'UTR
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (target), "UGAUAACAAGGAUCUACCUCA", strlen ("ACUAUACAACCUACUACCUCA"));

    Kd = mirbooking_score_table_compute_score (score_table, mirna, target, 13, NULL);

    g_assert_cmpfloat (Kd, >=, 11 - 2);
    g_assert_cmpfloat (Kd, <=, 11 + 2);
}

int main (int argc, gchar **argv)
{
    g_test_init (&argc, &argv, NULL);

    g_test_add_func ("/score-table/compute-seed-score", test_score_table_compute_seed_score);
    g_test_add_func ("/score-table/compute-seed-scores", test_score_table_compute_seed_scores);
    g_test_add_func ("/score-table/mcff", test_score_table_mcff);
    g_test_add_func ("/score-table/wee-et-al-2012", test_score_table_wee_et_al_2012);
    g_test_add_func ("/score-table/solomon-et-al-2016", test_score_table_solomon_et_al_2016);

    return g_test_run ();
}
