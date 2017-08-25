#include <glib.h>
#include <mirbooking.h>
#include <string.h>

guint16 SCORE_INDEX[16384 * 16384 * 2];

static void
test_score_index ()
{
    g_autoptr (MirbookingPrecomputedScoreIndex) index = mirbooking_precomputed_score_index_new (SCORE_INDEX, 1, 7);

    g_autoptr (MirbookingScoreIndexIter) iter = mirbooking_score_index_iterator (MIRBOOKING_SCORE_INDEX (index));

    g_autoptr (MirbookingTarget) target = mirbooking_target_new ("NM");
    g_autoptr (MirbookingMirna) mirna = mirbooking_mirna_new ("MIMAT");

    gchar *target_seq = "GCACACAGCACACAGAGCAGCATAAAGCCCAGTTGCTTTGGGAAGTGTTTGGGACCAGATGGATTGT";
    gchar *mirna_seq = "UGAGGUAGUAGGUUGUAUAGUU";

    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (target), target_seq, strlen (target_seq));
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (mirna), mirna_seq, strlen (mirna_seq));

    SCORE_INDEX[0] = GUINT16_TO_BE ((guint16) mirbooking_sequence_get_subsequence_index (MIRBOOKING_SEQUENCE (mirna), 1, 7));
    SCORE_INDEX[1] = GUINT16_TO_BE ((guint16) mirbooking_sequence_get_subsequence_index (MIRBOOKING_SEQUENCE (target), 0, 7));

    mirbooking_score_index_set_sequence_quantity (MIRBOOKING_SCORE_INDEX (index), MIRBOOKING_SEQUENCE (target), 5.0f);
    mirbooking_score_index_set_sequence_quantity (MIRBOOKING_SCORE_INDEX (index), MIRBOOKING_SEQUENCE (mirna), 5.0f);

    MirbookingMirna *a;
    MirbookingTarget *b;
    gsize position;
    g_assert (mirbooking_score_index_iter_next (iter, &a, &b, &position));
    g_assert (mirna == a);
    g_assert (target == b);
    g_assert_cmpint (position, ==, 7);

    g_assert (mirbooking_score_index_iter_next (iter, &a, &b, &position));
    g_assert (mirna == a);
    g_assert (target == b);
    g_assert_cmpint (position, ==, 0);
}

int
main (gint argc, gchar **argv)
{
    g_test_init (&argc, &argv, NULL);

    g_test_add_func ("/score-index", test_score_index);

    return g_test_run ();
}
