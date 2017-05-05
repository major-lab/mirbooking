#include <glib.h>

#include <mirbooking.h>
#include <string.h>

static gfloat SCORE_TABLE[16384][16384];

static void
test_score_table_compute_seed_score ()
{
    g_autoptr (GBytes) precomputed_table = g_bytes_new_static (SCORE_TABLE, sizeof (SCORE_TABLE));
    g_autoptr (MirbookingScoreTable) score_table = mirbooking_score_table_new_precomputed (precomputed_table);

    g_autoptr (MirbookingTarget) target = mirbooking_target_new ("NM");
    g_autoptr (MirbookingMirna) mirna = mirbooking_mirna_new ("MIMAT");

    gchar *target_seq = "GCACACAGAGCAGCATAAAGCCCAGTTGCTTTGGGAAGTGTTTGGGACCAGATGGATTGT";
    gchar *mirna_seq = "UGAGGUAGUAGGUUGUAUAGUU";

    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (target), target_seq, strlen (target_seq));
    mirbooking_sequence_set_sequence (MIRBOOKING_SEQUENCE (mirna), mirna_seq, strlen (mirna_seq));

    g_assert_cmpstr (g_strndup (mirbooking_sequence_get_subsequence (MIRBOOKING_SEQUENCE (target), 1, 7), 7), ==, "CACACAG");
    g_assert_cmpstr (g_strndup (mirbooking_sequence_get_subsequence (MIRBOOKING_SEQUENCE (mirna), 1, 7), 7), ==, "GAGGUAG");

    union {
        gint32 i;
        gfloat f;
    } val;
    val.f = 5.6f;
    val.i = GINT32_TO_BE (val.i);
    SCORE_TABLE[4370][8882] = val.f;

    // test for a seed
    gfloat site_score = mirbooking_score_table_compute_score (score_table,
                                                              MIRBOOKING_SEQUENCE (target), 1, // CACACAG =>
                                                              MIRBOOKING_SEQUENCE (mirna), 1,  // GAGGUAG =>
                                                              7);

    g_assert_cmpfloat (site_score, ==, 5.6f);
}

int main (int argc, gchar **argv)
{
    g_test_init (&argc, &argv, NULL);

    g_test_add_func ("/score-table/compute-seed-score", test_score_table_compute_seed_score);

    return g_test_run ();
}
