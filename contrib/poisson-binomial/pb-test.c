#include <pb.h>
#include <glib.h>

int main(void)
{
    PoissonBinomial pb;
    double p[] = {0.9, 0.3, 0.1};
    pb_init (&pb, p, 3);

    g_assert_cmpfloat ((gfloat)pb_pmf (&pb, 0), ==, 0.063f);
    g_assert_cmpfloat ((gfloat)pb_pmf (&pb, 1), ==, 0.601f);
    g_assert_cmpfloat ((gfloat)pb_pmf (&pb, 2), ==, 0.309f);
    g_assert_cmpfloat ((gfloat)pb_pmf (&pb, 3), ==, 0.027f);

    pb_destroy (&pb);

    return 0;
}
