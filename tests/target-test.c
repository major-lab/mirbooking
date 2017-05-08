#include <glib.h>
#include <mirbooking.h>
#include <string.h>

void
test_target_cds ()
{
    g_autoptr (MirbookingTarget) target = mirbooking_target_new ("NM123123");

    gchar *seq = "ACTGACTG";
    mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (target), seq, strlen (seq));

    mirbooking_target_set_cds (target, 1, 3);
    mirbooking_target_set_cds (target, 1, 7);

    gsize cds_len;
    gsize cds = mirbooking_target_get_cds (target, &cds_len);
    g_assert_cmpint (cds, ==, 1);
    g_assert_cmpint (cds_len, ==, 7);
}

int main (int argc, gchar **argv)
{
    g_test_init (&argc, &argv, NULL);

    g_test_add_func ("/target/cds", test_target_cds);

    return g_test_run ();
}
