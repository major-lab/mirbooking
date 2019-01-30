#include <glib.h>

#include <mirbooking.h>
#include <string.h>

static void
test_sequence_trim_line_feed ()
{
    gsize seq_len;
    const guint8 *seq;
    const guint8 *trimmed_seq;

    g_autoptr (MirbookingSequence) sequence = MIRBOOKING_SEQUENCE (mirbooking_target_new ("some-unique-accession"));

    mirbooking_sequence_set_sequence (sequence, "ACTG\nACTG");
    seq = mirbooking_sequence_get_raw_sequence (sequence, &seq_len);
    g_assert (seq_len == 9);
    g_assert_cmpmem (seq, 9, "ACTG\nACTG", 9);
    trimmed_seq = mirbooking_sequence_get_subsequence (sequence, 0, 8);
    g_assert_cmpmem (trimmed_seq, 8, "ACTGACTG", 8);

    // leading
    mirbooking_sequence_set_sequence (sequence, "\nACTGACTG");
    seq = mirbooking_sequence_get_raw_sequence (sequence, &seq_len);
    g_assert (seq_len == 9);
    g_assert_cmpmem ("\nACTGACTG", 9, seq, 9);
    trimmed_seq = mirbooking_sequence_get_subsequence (sequence, 0, 8);
    g_assert (trimmed_seq == seq + 1); // zero-copy

    // trailing
    mirbooking_sequence_set_sequence (sequence, "ACTGACTG\n");
    seq = mirbooking_sequence_get_raw_sequence (sequence, &seq_len);
    g_assert (seq_len == 9);
    g_assert_cmpmem ("ACTGACTG\n", 9, seq, 9);
    trimmed_seq = mirbooking_sequence_get_subsequence (sequence, 0, 8);
    g_assert (trimmed_seq == seq); // zero-copy

    // worst-case scenario
    mirbooking_sequence_set_sequence (sequence, "\nAC\nTG\nAC\nTG\n");
    seq = mirbooking_sequence_get_raw_sequence (sequence, &seq_len);
    g_assert (seq_len == 13);
    g_assert_cmpmem ("\nAC\nTG\nAC\nTG\n", 13, seq, 13);
    trimmed_seq = mirbooking_sequence_get_subsequence (sequence, 0, 8);
    g_assert_cmpmem (trimmed_seq, 8, "ACTGACTG", 8);

    // worst, worst-case scenario
    mirbooking_sequence_set_sequence (sequence, "\nAC\n\nTG\n");
    seq = mirbooking_sequence_get_raw_sequence (sequence, &seq_len);
    g_assert (seq_len == 8);
    g_assert_cmpmem ("\nAC\n\nTG\n", 8, seq, 8);
    trimmed_seq = mirbooking_sequence_get_subsequence (sequence, 0, 4);
    g_assert_cmpmem (trimmed_seq, 4, "ACTG", 4);

    // worst, worst-case scenario
    mirbooking_sequence_set_sequence (sequence, "\nAC\n\nTG\n\n\n\nATAT");
    seq = mirbooking_sequence_get_raw_sequence (sequence, &seq_len);
    g_assert (seq_len == 15);
    g_assert_cmpmem ("\nAC\n\nTG\n\n\n\nATAT", 15, seq, 15);
    trimmed_seq = mirbooking_sequence_get_subsequence (sequence, 2, 2);
    g_assert_cmpmem (trimmed_seq, 2, "TG", 2);
    g_assert (trimmed_seq == seq + 5); // zero-copy
}

static void
test_sequence_unreliable_subsequence_index (void)
{
    g_autoptr (MirbookingSequence) sequence = MIRBOOKING_SEQUENCE (mirbooking_target_new ("NM_1234"));

    mirbooking_sequence_set_sequence (sequence, "ACTGY");

    g_assert_cmpint (mirbooking_sequence_get_subsequence_index (sequence, 0, 4), ==, 30);
    g_assert_cmpint (mirbooking_sequence_get_subsequence_index (sequence, 0, 5), ==, -1);
}

int main (int argc, gchar **argv)
{
    g_test_init (&argc, &argv, NULL);

    g_test_add_func ("/sequence/trim-line-feed", test_sequence_trim_line_feed);
    g_test_add_func ("/sequence/unreliable-subsequence-index", test_sequence_unreliable_subsequence_index);

    return g_test_run ();
}
