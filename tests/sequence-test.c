#include <glib.h>

#include <mirbooking.h>
#include <string.h>

static void
test_sequence_trim_line_feed ()
{
    gsize seq_len;
    const gchar *seq;
    const gchar *trimmed_seq;

    g_autoptr (MirbookingSequence) sequence = MIRBOOKING_SEQUENCE (mirbooking_target_new ("some-unique-accession"));

    mirbooking_sequence_set_raw_sequence (sequence, "ACTG\nACTG", 9);
    seq = mirbooking_sequence_get_raw_sequence (sequence, &seq_len);
    g_assert (seq_len == 9);
    g_assert_cmpstr (seq, ==, "ACTG\nACTG");
    trimmed_seq = mirbooking_sequence_get_subsequence (sequence, 0, 8);
    g_assert_cmpmem (trimmed_seq, 8, "ACTGACTG", 8);

    // leading
    mirbooking_sequence_set_raw_sequence (sequence, "\nACTGACTG", 9);
    seq = mirbooking_sequence_get_raw_sequence (sequence, &seq_len);
    g_assert (seq_len == 9);
    g_assert (strncmp ("\nACTGACTG", seq, 9) == 0);
    trimmed_seq = mirbooking_sequence_get_subsequence (sequence, 0, 8);
    g_assert (trimmed_seq == seq + 1); // zero-copy

    // trailing
    mirbooking_sequence_set_raw_sequence (sequence, "ACTGACTG\n", 9);
    seq = mirbooking_sequence_get_raw_sequence (sequence, &seq_len);
    g_assert (seq_len == 9);
    g_assert (strncmp ("ACTGACTG\n", seq, 9) == 0);
    trimmed_seq = mirbooking_sequence_get_subsequence (sequence, 0, 8);
    g_assert (trimmed_seq == seq); // zero-copy

    // worst-case scenario
    mirbooking_sequence_set_raw_sequence (sequence, "\nAC\nTG\nAC\nTG\n", 13);
    seq = mirbooking_sequence_get_raw_sequence (sequence, &seq_len);
    g_assert (seq_len == 13);
    g_assert (strncmp ("\nAC\nTG\nAC\nTG\n", seq, 13) == 0);
    trimmed_seq = mirbooking_sequence_get_subsequence (sequence, 0, 8);
    g_assert_cmpmem (trimmed_seq, 8, "ACTGACTG", 8);

    // worst, worst-case scenario
    mirbooking_sequence_set_raw_sequence (sequence, "\nAC\n\nTG\n", 8);
    seq = mirbooking_sequence_get_raw_sequence (sequence, &seq_len);
    g_assert (seq_len == 8);
    g_assert (strncmp ("\nAC\n\nTG\n", seq, 8) == 0);
    trimmed_seq = mirbooking_sequence_get_subsequence (sequence, 0, 4);
    g_assert_cmpmem (trimmed_seq, 4, "ACTG", 4);

    // worst, worst-case scenario
    mirbooking_sequence_set_raw_sequence (sequence, "\nAC\n\nTG\n\n\n\nATAT", 15);
    seq = mirbooking_sequence_get_raw_sequence (sequence, &seq_len);
    g_assert (seq_len == 15);
    g_assert (strncmp ("\nAC\n\nTG\n\n\n\nATAT", seq, 15) == 0);
    trimmed_seq = mirbooking_sequence_get_subsequence (sequence, 2, 2);
    g_assert_cmpmem (trimmed_seq, 2, "TG", 2);
    g_assert (trimmed_seq == seq + 5); // zero-copy
}

static void
test_sequence_unreliable_subsequence_index (void)
{
    g_autoptr (MirbookingSequence) sequence = MIRBOOKING_SEQUENCE (mirbooking_target_new ("NM_1234"));

    mirbooking_sequence_set_raw_sequence (sequence, "ACTGY", 5);

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
