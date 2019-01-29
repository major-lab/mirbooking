#include <glib.h>

#include <mirbooking.h>
#include <string.h>
#include <math.h>

#define R 1.987203611e-3
#define T 310.15

typedef struct __attribute__ ((packed)) _SeedScoreLayout
{
    gsize  n;
    gsize  nnz;
    gsize  rowptr[16384 + 1];
    gsize  colind[10];
    gfloat data[10];
} SeedScoreLayout;

static SeedScoreLayout SEED_SCORES =
{
    16384,
    10,
    {[657] = 0, [658] = 1, [8882] = 1, [8883] = 8, [12152] = 8, [12153] = 9, [13391] = 9, [13392] = 10},
    {11871,  3228,   4370,   4371,   4372,   7261,   7421,   9284,   13441,  952},
    {-9.79f, -2.12f, -20.0f, -19.0f, -18.0f, -9.37f, -4.28f, -19.0f, -7.79f, -6.34f}, // FIXME: 7.79f is incorrect
};

typedef struct __attribute__ ((packed)) _SupplementaryScore3Layout
{
    gsize  n;
    gsize  nnz;
    gsize  rowptr[64 + 1];
    gsize  colind[14];
    gfloat data[14];
} SupplementaryScore3Layout;

static SupplementaryScore3Layout SUPPLEMENTARY_SCORES =
{
    64,
    14,
    {
        [0] = 0, [1] = 1,
        [7] = 1, [8] = 2,
        [24] = 2, [25] = 3,
        [27] = 3, [28] = 4,
        [31] = 4,
        [32] = 5, [33] = 7,
        [36] = 7, [37] = 8,
        [39] = 9, [40] = 10,
        [44] = 10, [45] = 11,
        [47] = 11, [48] = 12,
        [50] = 12, [51] = 13,
        [56] = 13, [57] = 14
    },
    {63,   52,  9,   6,   2,     2,   33,  6,     9,     9,      49,  1,   28,  52},
    {2.05, 0.0, 0.0, 0.0, 1.16f, 0.0, 0.0, 1.20f, 1.96f, -1.11f, 0.0, 0.0, 0.0, -0.04}
};

typedef struct __attribute__ ((packed)) _SupplementaryScore4Layout
{
    gsize  n;
    gsize  nnz;
    gsize  rowptr[256 + 1];
    gsize  colind[5];
    gfloat data[5];
} SupplementaryScore4Layout;

static SupplementaryScore4Layout SUPPLEMENTARY_SCORES_4 =
{
    256,
    5,
    {[251] = 0, [252] = 5},
    {0,     16,     63,    96,    239},
    {1.58f, -0.83f, 2.63f, 0.29f, 0.0f},
};

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

    g_assert_cmpfloat (site_score, ==, 1e12 * exp ((-20.0f - 5.43f - 0.56f) / (R * T)));
    g_assert_null (error);
}

static void
test_score_table_compute_seed_scores ()
{
    g_autoptr (GBytes) default_table = g_bytes_new_static (&SEED_SCORES, sizeof (SEED_SCORES));
    g_autoptr (GBytes) supplementary_scores = g_bytes_new_static (&SUPPLEMENTARY_SCORES_4, sizeof (SUPPLEMENTARY_SCORES_4));
    g_autoptr (MirbookingScoreTable) score_table = MIRBOOKING_SCORE_TABLE (mirbooking_default_score_table_new (default_table,
                                                                                                               MIRBOOKING_DEFAULT_SCORE_TABLE_SUPPLEMENTARY_MODEL_ZAMORE_ET_AL_2012,
                                                                                                               supplementary_scores));

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
    g_autofree gsize *positions = NULL;
    gsize  positions_len;
    g_assert (mirbooking_score_table_compute_positions (MIRBOOKING_SCORE_TABLE (score_table),
                                                        mirna,     // GAGGUAG =>
                                                        target, &positions, // CACACAG =>
                                                        &positions_len,
                                                        &error));

    g_assert_nonnull (positions);
    g_assert_cmpint (positions_len, ==, 2);
    g_assert_cmpfloat (mirbooking_score_table_compute_score (score_table, mirna, target, positions[0], NULL), ==, 1e12 * exp ((-19.0f - 5.43f) / (R * T)));
    g_assert_cmpfloat (mirbooking_score_table_compute_score (score_table, mirna, target, positions[1], NULL), ==, 1e12 * exp ((-20.0f - 5.43f - 0.56f) / (R * T)));
    g_assert_null (error);
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
    g_assert_cmpfloat (site_score, ==, 1e12 * exp ((-7.691f + 4.40f) / (R * T)));
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
    g_autoptr (GBytes) supplementary_scores = g_bytes_new_static (&SUPPLEMENTARY_SCORES_4, sizeof (SUPPLEMENTARY_SCORES_4));
    g_autoptr (MirbookingScoreTable) score_table = MIRBOOKING_SCORE_TABLE (mirbooking_default_score_table_new (default_table,
                                                                                                               MIRBOOKING_DEFAULT_SCORE_TABLE_SUPPLEMENTARY_MODEL_ZAMORE_ET_AL_2012,
                                                                                                               supplementary_scores));

    // siRNA
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (mirna), "UGAGGUAGUAGGUUGUAUAUGU", strlen ("UGAGGUAGUAGGUUGUAUAUGU"));
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (target), "GAUACUAUACAACCUACUACCUCAACCU", strlen ("GAUACUAUACAACCUACUACCUCAACCU"));

    g_autofree gsize *positions = NULL;
    gsize  positions_len;
    mirbooking_score_table_compute_positions (MIRBOOKING_SCORE_TABLE (score_table), mirna, target, &positions, &positions_len, NULL);
    g_assert_cmpint (positions_len, ==, 1);
    g_assert_cmpint (positions[0], ==, 16);

    Kd = mirbooking_score_table_compute_score (score_table, mirna, target, 16, NULL);
    Km = mirbooking_score_table_compute_enzymatic_score (score_table, mirna, target, 16, NULL);

    g_assert_cmpfloat (Kd, ==, 1e12 * exp ((-0.56f - 9.37f - 0.83f - 5.43f) / (R * T)));
    // FIXME: g_assert_cmpfloat (Kd, >=, 20 - 10);
    g_assert_cmpfloat (Kd, <=, 20 + 10);

    g_assert_cmpfloat (Km, ==, Kd + (MIRBOOKING_SCORE_TABLE_DEFAULT_KCAT / MIRBOOKING_SCORE_TABLE_DEFAULT_KF));
    g_assert_cmpfloat (Km, >=, 100 - 60);
    // FIXME: g_assert_cmpfloat (Km, <=, 100 + 60);

    // seed-only
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (target), "GAAAAAAAAAAAAAAUCUACCUCUAAAU", strlen ("GAAAAAAAAAAAAAAUCUACCUCUAAAU"));

    Kd = mirbooking_score_table_compute_score (score_table, mirna, target, 16, NULL);
    Km = mirbooking_score_table_compute_enzymatic_score (score_table, mirna, target, 16, NULL);

    g_assert_cmpfloat (Kd, ==, 1e12 * exp ((-9.37f - 5.43f) / (R * T)));
    g_assert_cmpfloat (Kd, >=, 26 - 2);
    // FIXME: g_assert_cmpfloat (Kd, <=, 26 + 2);

    g_assert_cmpfloat (Km, ==, Kd + (MIRBOOKING_SCORE_TABLE_DEFAULT_KCAT / MIRBOOKING_SCORE_TABLE_DEFAULT_KF));
    g_assert_cmpfloat (Km, >=, 100 - 60);
    // FIXME: g_assert_cmpfloat (Km, <=, 100 + 60);

    // seed and supplementary
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (target), "GAAAAAAAACAAAAAUCUACCUCUAAAU", strlen ("GAAAAAAAAAAAAAAUCUACCUCUAAAU"));

    Kd = mirbooking_score_table_compute_score (score_table, mirna, target, 16, NULL);
    Km = mirbooking_score_table_compute_enzymatic_score (score_table, mirna, target, 16, NULL);

    g_assert_cmpfloat (Kd, ==, 1e12 * exp ((-9.37f - 0.83f - 5.43f) / (R * T)));
    // FIXME: g_assert_cmpfloat (Kd, >=, 13 - 1);
    g_assert_cmpfloat (Kd, <=, 13 + 1);

    g_assert_cmpfloat (Km, ==, Kd + (MIRBOOKING_SCORE_TABLE_DEFAULT_KCAT / MIRBOOKING_SCORE_TABLE_DEFAULT_KF));
    g_assert_cmpfloat (Km, >=, 100 - 60);
    // FIXME: g_assert_cmpfloat (Km, <=, 100 + 60);

    // g10g11 central internal loop 50±30
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (target), "GAUACUAUACAACGAACUACCUCAACCU", strlen ("GAUACUAUACAACGAACUACCUCAACCU"));

    Kd = mirbooking_score_table_compute_score (score_table, mirna, target, 16, NULL);
    Km = mirbooking_score_table_compute_enzymatic_score (score_table, mirna, target, 16, NULL);

    g_assert_cmpfloat (Kd, ==, 1e12 * exp ((-0.56f -9.37f - 0.83f - 5.43f) / (R * T)));
    // FIXME: g_assert_cmpfloat (Kd, >=, 50 - 30);
    g_assert_cmpfloat (Kd, <=, 50 + 30);

    // g15g16 mismatches in 3' supplementary region 30.83
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (target), "GAUACUAUCGAACCUACUACCUCAACCU", strlen ("GAUACUAUCGAACCUACUACCUCAACCU"));

    Kd = mirbooking_score_table_compute_score (score_table, mirna, target, 16, NULL);
    Km = mirbooking_score_table_compute_enzymatic_score (score_table, mirna, target, 16, NULL);

    // g_assert_cmpfloat (Kd, ==, 1e12 * exp ((-9.37f - 5.43f - 0.56f) / (R * T)));
    g_assert_cmpfloat (Kd, >=, 30 - 20);
    g_assert_cmpfloat (Kd, <=, 30 + 20);

    // g1-g19 complementary 40.83
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (target), "GAUUAUAUACAACCUACUACCUCAACCU", strlen ("GAUUAUAUACAACCUACUACCUCAACCU"));

    Kd = mirbooking_score_table_compute_score (score_table, mirna, target, 16, NULL);
    Km = mirbooking_score_table_compute_enzymatic_score (score_table, mirna, target, 16, NULL);

    g_assert_cmpfloat (Kd, ==, 1e12 * exp ((-0.56f - 9.37f - 0.83f - 5.43f) / (R * T)));
    // FIXME: g_assert_cmpfloat (Kd, >=, 40 - 20);
    g_assert_cmpfloat (Kd, <=, 40 + 20);

    // seed-only 30.83
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (target), "GAAAAAAAAAAAAAAUCUACCUCUAAAU", strlen ("GAAAAAAAAAAAAAAUCUACCUCUAAAU"));

    Kd = mirbooking_score_table_compute_score (score_table, mirna, target, 16, NULL);
    Km = mirbooking_score_table_compute_enzymatic_score (score_table, mirna, target, 16, NULL);

    g_assert_cmpfloat (Kd, ==, 1e12 * exp ((-9.37f - 5.43f) / (R * T)));
    g_assert_cmpfloat (Kd, >=, 30 - 20);
    g_assert_cmpfloat (Kd, <=, 30 + 20);

    // seed plus g13-g16 3' supplementary 20±10
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (target), "GAAAAAAAACAAAAAUCUACCUCUAAAU", strlen ("GAAAAAAAACAAAAAUCUACCUCUAAAU"));

    Kd = mirbooking_score_table_compute_score (score_table, mirna, target, 16, NULL);
    Km = mirbooking_score_table_compute_enzymatic_score (score_table, mirna, target, 16, NULL);

    g_assert_cmpfloat (Kd, ==, 1e12 * exp ((-9.37f - 0.83f - 5.43f) / (R * T)));
    // FIXME: g_assert_cmpfloat (Kd, >=, 20 - 10);
    g_assert_cmpfloat (Kd, <=, 20 + 10);

    // g4g5 mismatches in seed 1e3±0.6e3
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (target), "GAUACUAUACAACCUACUAUUUCAACCU", strlen ("GAUACUAUACAACCUACUAUUUCAACCU"));

    Kd = mirbooking_score_table_compute_score (score_table, mirna, target, 16, NULL);
    Km = mirbooking_score_table_compute_enzymatic_score (score_table, mirna, target, 16, NULL);

    g_assert_cmpfloat (Kd, ==, 1e12 * exp ((-0.56f - 4.28f - 0.83f - 5.43f) / (R * T)));
    g_assert_cmpfloat (Kd, >=, 1e3 - 0.6e3);
    // FIXME: g_assert_cmpfloat (Kd, <=, 1e3 + 0.6e3);

    // non-complementary luciferase target 2e3±1e3
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (target), "GAUAAGGCAUUUCAUUAUAGCUAUACCU", strlen ("GAUAAGGCAUUUCAUUAUAGCUAUACCU"));

    Kd = mirbooking_score_table_compute_score (score_table, mirna, target, 16, NULL);
    Km = mirbooking_score_table_compute_enzymatic_score (score_table, mirna, target, 16, NULL);

    g_assert_cmpfloat (Kd, ==, 1e12 * exp ((-2.12f - 5.43f) / (R * T)));
    g_assert_cmpfloat (Kd, >=, 2e3 - 1e3);
    // FIXME: g_assert_cmpfloat (Kd, <=, 2e3 + 1e3);
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
    g_autoptr (GBytes) supplementary_scores = g_bytes_new_static (&SUPPLEMENTARY_SCORES_4, sizeof (SUPPLEMENTARY_SCORES_4));
    g_autoptr (MirbookingScoreTable) score_table = MIRBOOKING_SCORE_TABLE (mirbooking_default_score_table_new (default_table,
                                                                                                               MIRBOOKING_DEFAULT_SCORE_TABLE_SUPPLEMENTARY_MODEL_ZAMORE_ET_AL_2012,
                                                                                                               supplementary_scores));

    // Complete
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (mirna),  "UGAGGUAGUAGGUUGUAUAGU", strlen ("UGAGGUAGUAGGUUGUAUAGU"));
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (target), "ACUAUACAACCUACUACCUCA", strlen ("ACUAUACAACCUACUACCUCA"));

    Kd = mirbooking_score_table_compute_score (score_table, mirna, target, 13, NULL);

    g_assert_cmpfloat (Kd, ==, 1e12 * exp ((-0.56f - 9.37f - 0.83f - 5.43f) / (R * T)));

    // Seed only
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (target), "UGAUAUGUUGGAUCUACCUCA", strlen ("ACUAUACAACCUACUACCUCA"));

    g_assert_cmpfloat (MIRBOOKING_SCORE_TABLE_DEFAULT_KF, >=, 2.4e-4 - 0.1e-4);
    g_assert_cmpfloat (MIRBOOKING_SCORE_TABLE_DEFAULT_KF, <=, 2.4e-4 + 0.1e-4);
    g_assert_cmpfloat (MIRBOOKING_SCORE_TABLE_DEFAULT_KCAT, >=, 3.6e-2 - 0.2e-2);
    g_assert_cmpfloat (MIRBOOKING_SCORE_TABLE_DEFAULT_KCAT, <=, 3.6e-2 + 0.2e-2);

    Kd = mirbooking_score_table_compute_score (score_table, mirna, target, 13, NULL);

    g_assert_cmpfloat (Kd, ==, 1e12 * exp ((-0.56f - 9.37f - 5.43f) / (R * T)));
    g_assert_cmpfloat (Kd, >=, 15 - 2);
    g_assert_cmpfloat (Kd, <=, 15 + 2);

    // Seed plus 3'UTR
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (target), "UGAUAACAAGGAUCUACCUCA", strlen ("ACUAUACAACCUACUACCUCA"));

    Kd = mirbooking_score_table_compute_score (score_table, mirna, target, 13, NULL);

    g_assert_cmpfloat (Kd, ==, 1e12 * exp ((-0.56f - 9.37f - 0.83f - 5.43f) / (R * T)));
    // FIXME: g_assert_cmpfloat (Kd, >=, 11 - 2);
    g_assert_cmpfloat (Kd, <=, 11 + 2);
}

static void
test_score_table_schirle_et_al_2015 ()
{
    gdouble Kd;
    g_autoptr (MirbookingTarget) target = mirbooking_target_new ("");

    g_autoptr (MirbookingMirna) mirna = mirbooking_mirna_new ("");
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (mirna), "UUCACAUUGCCCAAGUCUCUU", strlen ("UUCACAUUGCCCAAGUCUCUU"));

    g_autoptr (GBytes) default_table = g_bytes_new_static (&SEED_SCORES, sizeof (SEED_SCORES));
    g_autoptr (GBytes) supplementary_scores = g_bytes_new_static (&SUPPLEMENTARY_SCORES_4, sizeof (SUPPLEMENTARY_SCORES_4));
    g_autoptr (MirbookingScoreTable) score_table = MIRBOOKING_SCORE_TABLE (mirbooking_default_score_table_new (default_table,
                                                                                                               MIRBOOKING_DEFAULT_SCORE_TABLE_SUPPLEMENTARY_MODEL_ZAMORE_ET_AL_2012,
                                                                                                               supplementary_scores));

    // A
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (target), "CAAUGUGAAAA", strlen ("CAAUGUGAAAA"));
    Kd = mirbooking_score_table_compute_score (score_table, mirna, target, 1, NULL);
    g_assert_cmpfloat (Kd, ==, 1e12 * exp ((-0.56f - 6.34f - 5.43f) / (R * T)));
    // FIXME: g_assert_cmpfloat (Kd, <=, 0.75e3 + 0.04e3);
    g_assert_cmpfloat (Kd, >=, 0.75e3 - 0.04e3);

    // U
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (target), "CAAUGUGAUAA", strlen ("CAAUGUGAUAA"));
    Kd = mirbooking_score_table_compute_score (score_table, mirna, target, 1, NULL);
    g_assert_cmpfloat (Kd, ==, 1e12 * exp ((-6.34f - 5.43f) / (R * T)));
    // FIXME: g_assert_cmpfloat (Kd, <=, 1.9e3 + 0.09e3);
    g_assert_cmpfloat (Kd, >=, 1.9e3 - 0.09e3);

    // C
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (target), "CAAUGUGACAA", strlen ("CAAUGUGACAA"));
    Kd = mirbooking_score_table_compute_score (score_table, mirna, target, 1, NULL);
    g_assert_cmpfloat (Kd, ==, 1e12 * exp ((-6.34f - 5.43f) / (R * T)));
    // FIXME: g_assert_cmpfloat (Kd, <=, 1.9e3 + 0.10e3);
    g_assert_cmpfloat (Kd, >=, 1.9e3 - 0.10e3);

    // G
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (target), "CAAUGUGAGAA", strlen ("CAAUGUGAGAA"));
    Kd = mirbooking_score_table_compute_score (score_table, mirna, target, 1, NULL);
    g_assert_cmpfloat (Kd, ==, 1e12 * exp ((-6.34f - 5.43f) / (R * T)));
    // FIXME: g_assert_cmpfloat (Kd, <=, 1.8e3 + 0.12e3);
    g_assert_cmpfloat (Kd, >=, 1.8e3 - 0.12e3);
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

    g_autoptr (GBytes) default_table = g_bytes_new_static (&SEED_SCORES, sizeof (SEED_SCORES));
    g_autoptr (GBytes) supplementary_scores = g_bytes_new_static (&SUPPLEMENTARY_SCORES, sizeof (SUPPLEMENTARY_SCORES));
    g_autoptr (MirbookingScoreTable) score_table = MIRBOOKING_SCORE_TABLE (mirbooking_default_score_table_new (default_table,
                                                                                                               MIRBOOKING_DEFAULT_SCORE_TABLE_SUPPLEMENTARY_MODEL_YAN_ET_AL_2018,
                                                                                                               supplementary_scores));

    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (target), "GACGAAGAGCUCAUCAGAACA", strlen ("GACGAAGAGCUCAUCAGAACA"));

    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (mirna),  "UGUUCUGAUGAGCUCUUCGUC", strlen ("UGUUCUGAUGAGCUCUUCGUC"));

    g_assert_cmpmem ("GUUCUGA", 7, mirbooking_sequence_get_subsequence (MIRBOOKING_SEQUENCE (mirna), 1, 7), 7);
    g_assert_cmpmem ("UGA", 3, mirbooking_sequence_get_subsequence     (MIRBOOKING_SEQUENCE (mirna), 8, 3), 3);
    g_assert_cmpmem ("GCU",3, mirbooking_sequence_get_subsequence      (MIRBOOKING_SEQUENCE (mirna), 11, 3), 3);
    g_assert_cmpmem ("CUU", 3, mirbooking_sequence_get_subsequence     (MIRBOOKING_SEQUENCE (mirna), 14, 3), 3);
    g_assert_cmpmem ("CGUC", 4, mirbooking_sequence_get_subsequence    (MIRBOOKING_SEQUENCE (mirna), 17, 4), 4);

    gdouble Kd;

    Kd = mirbooking_score_table_compute_score (score_table,
                                               mirna,
                                               target,
                                               13,
                                               NULL);

    g_assert_cmpfloat (Kd, ==, 1e12 * exp ((-0.56f -7.79f - 1.11f - 5.43f) / (R * T)));

    // A+C
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (mirna),  "UGUUCUGAACUGCUGAACGUC", strlen ("UGUUCUGAUGAGCUCUUCGUC"));
    Kd = mirbooking_score_table_compute_score (score_table,
                                               mirna,
                                               target,
                                               13,
                                               NULL);

    g_assert_cmpfloat (Kd, ==, 1e12 * exp ((-0.56f -7.79f - 1.11f - 5.43f) / (R * T)));

    // B+D
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (mirna),  "UGUUCUGAUGACGACUUCGUC", strlen ("UGUUCUGAUGAGCUCUUCGUC"));
    Kd = mirbooking_score_table_compute_score (score_table,
                                               mirna,
                                               target,
                                               13,
                                               NULL);

    g_assert_cmpfloat (Kd, ==, 1e12 * exp ((-0.56f -7.79f - 5.43f) / (R * T)));

    // A+D
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (mirna),  "UGUUCUGAACUGCUCUUGCAG", strlen ("UGUUCUGAUGAGCUCUUCGUC"));
    Kd = mirbooking_score_table_compute_score (score_table,
                                               mirna,
                                               target,
                                               13,
                                               NULL);

    g_assert_cmpfloat (Kd, ==, 1e12 * exp ((-0.56f -7.79f - 1.11f - 5.43f) / (R * T)));

    // A+B
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (mirna),  "UGUUCUGAACUGCACUUCGUC", strlen ("UGUUCUGAUGAGCUCUUCGUC"));
    Kd = mirbooking_score_table_compute_score (score_table,
                                               mirna,
                                               target,
                                               13,
                                               NULL);

    g_assert_cmpfloat (Kd, ==, 1e12 * exp ((-0.56f -7.79f - 5.43f) / (R * T)));

    // B+C
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (mirna),  "UGUUCUGAUGACGAGAACGUC", strlen ("UGUUCUGAUGAGCUCUUCGUC"));
    Kd = mirbooking_score_table_compute_score (score_table,
                                               mirna,
                                               target,
                                               13,
                                               NULL);

    g_assert_cmpfloat (Kd, ==, 1e12 * exp ((-0.56f - 7.79f - 5.43f) / (R * T)));

    // C+D
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (mirna),  "UGUUCUGAUGAGCUGAAGCAG", strlen ("UGUUCUGAUGAGCUCUUCGUC"));
    Kd = mirbooking_score_table_compute_score (score_table,
                                               mirna,
                                               target,
                                               13,
                                               NULL);

    g_assert_cmpfloat (Kd, ==, 1e12 * exp ((-0.56f -7.79f - 1.11f - 5.43f) / (R * T)));
}

static void
test_score_table_jo_et_al_2015 ()
{
    g_autoptr (GBytes) default_table = g_bytes_new_static (&SEED_SCORES, sizeof (SEED_SCORES));
    g_autoptr (GBytes) supplementary_scores = g_bytes_new_static (&SUPPLEMENTARY_SCORES, sizeof (SUPPLEMENTARY_SCORES));
    g_autoptr (MirbookingScoreTable) score_table = MIRBOOKING_SCORE_TABLE (mirbooking_default_score_table_new (default_table,
                                                                                                                          MIRBOOKING_DEFAULT_SCORE_TABLE_SUPPLEMENTARY_MODEL_YAN_ET_AL_2018,
                                                                                                                          supplementary_scores));

    g_autoptr (MirbookingMirna) mirna = mirbooking_mirna_new ("");
    g_autoptr (MirbookingTarget) target = mirbooking_target_new ("");

    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (mirna), "UGAGGUAGUAGGUUGUAUAGU", 21);
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (target), "ACUAUACAACCUACUACCUCG", 21);

    // gdouble Kd;
    // Kd = mirbooking_score_table_compute_score (score_table, mirna, target, 13, NULL);

    // g_assert_cmpfloat (MIRBOOKING_SCORE_TABLE_DEFAULT_KF, ==, 1.7e-5);
    // g_assert_cmpfloat (Kd, ==, 17);
}

static void
test_score_table_chi_et_al_2012 ()
{
    g_autoptr (GBytes) default_table = g_bytes_new_static (&SEED_SCORES, sizeof (SEED_SCORES));
    g_autoptr (GBytes) supplementary_scores = g_bytes_new_static (&SUPPLEMENTARY_SCORES_4, sizeof (SUPPLEMENTARY_SCORES_4));
    g_autoptr (MirbookingScoreTable) score_table = MIRBOOKING_SCORE_TABLE (mirbooking_default_score_table_new (default_table,
                                                                                                                          MIRBOOKING_DEFAULT_SCORE_TABLE_SUPPLEMENTARY_MODEL_ZAMORE_ET_AL_2012,
                                                                                                                          supplementary_scores));

    g_autoptr (MirbookingTarget) target = mirbooking_target_new ("");
    g_autoptr (MirbookingMirna) mirna = mirbooking_mirna_new ("");

    // G-buldged seed
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (target), "GUGGCCUU", 8);
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (mirna),  "UAAGGCAC", 8);

    gdouble Kd = mirbooking_score_table_compute_score (score_table, mirna, target, 0, NULL);

    g_assert_cmpfloat (Kd, ==, 1e12 * exp ((-9.79f + 1.2f - 5.43f) / (R * T)));
}

int main (int argc, gchar **argv)
{
    g_test_init (&argc, &argv, NULL);

    g_test_add_func ("/score-table/compute-seed-score", test_score_table_compute_seed_score);
    g_test_add_func ("/score-table/compute-seed-scores", test_score_table_compute_seed_scores);
    g_test_add_func ("/score-table/mcff", test_score_table_mcff);
    g_test_add_func ("/score-table/wee-et-al-2012", test_score_table_wee_et_al_2012);
    g_test_add_func ("/score-table/solomon-et-al-2016", test_score_table_solomon_et_al_2016);
    g_test_add_func ("/score-table/schirle-et-al-2015", test_score_table_schirle_et_al_2015);
    g_test_add_func ("/score-table/yan-et-al-2018", test_score_table_yan_et_al_2018);
    g_test_add_func ("/score-table/jo-et-al-2015", test_score_table_jo_et_al_2015);
    g_test_add_func ("/score-table/chi-et-al-2012", test_score_table_chi_et_al_2012);

    return g_test_run ();
}
