#include <glib.h>

#include <mirbooking.h>
#include <string.h>
#include <math.h>

static gfloat SCORE_TABLE[16384][16384];

static gfloat
gfloat_to_be (gfloat f)
{
    union {
        gint32 i;
        gfloat f;
    } val;
    val.f = f;
    val.i = GINT32_TO_BE (val.i);
    return val.f;
}

static void
test_score_table_compute_seed_score ()
{
    g_autoptr (GBytes) precomputed_table = g_bytes_new_static (SCORE_TABLE, sizeof (SCORE_TABLE));
    g_autoptr (MirbookingPrecomputedScoreTable) score_table = mirbooking_precomputed_score_table_new_from_bytes (precomputed_table, 1, 7);

    g_autoptr (MirbookingTarget) target = mirbooking_target_new ("NM");
    g_autoptr (MirbookingMirna) mirna = mirbooking_mirna_new ("MIMAT");

    gchar *target_seq = "GCACACAGAGCAGCATAAAGCCCAGTTGCTTTGGGAAGTGTTTGGGACCAGATGGATTGT";
    gchar *mirna_seq = "UGAGGUAGUAGGUUGUAUAGUU";

    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (target), target_seq, strlen (target_seq));
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (mirna), mirna_seq, strlen (mirna_seq));

    g_assert_cmpmem (mirbooking_sequence_get_subsequence (MIRBOOKING_SEQUENCE (target), 1, 7), 7, "CACACAG", 7);
    g_assert_cmpmem (mirbooking_sequence_get_subsequence (MIRBOOKING_SEQUENCE (mirna), 1, 7), 7, "GAGGUAG", 7);

    // 3-state Botlzmann
    SCORE_TABLE[8882][4370] = gfloat_to_be (-20.0);
    SCORE_TABLE[8882][4371] = gfloat_to_be (-19.0);
    SCORE_TABLE[8882][4372] = gfloat_to_be (-18.0);

    // test for a seed
    g_autoptr (GError) error = NULL;
    gfloat site_score = mirbooking_score_table_compute_score (MIRBOOKING_SCORE_TABLE (score_table),
                                                              mirna,     // GAGGUAG =>
                                                              target, 1, // CACACAG =>
                                                              &error);

    g_assert_cmpfloat (site_score, ==, -20.0);
    g_assert_null (error);
}

static void
test_score_table_compute_seed_scores ()
{

}

int main (int argc, gchar **argv)
{
    g_test_init (&argc, &argv, NULL);

    g_test_add_func ("/score-table/compute-seed-score", test_score_table_compute_seed_score);
    g_test_add_func ("/score-table/compute-seed-scores", test_score_table_compute_seed_scores);

    return g_test_run ();
}
