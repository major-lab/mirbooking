#include <glib.h>
#include <gio/gio.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/mman.h>

G_DEFINE_AUTOPTR_CLEANUP_FUNC (FILE, fclose)

static gchar  *mcff           = "mcff";
static gsize   seed_length    = 7;
static gint    max_mismatches = 1;
static gchar  *output         = NULL;

static gchar *nt = "ACGT";

static const GOptionEntry MIRBOOKING_GENERATE_SCORE_TABLE_OPTIONS[] =
{
    {"mcff",           0, 0, G_OPTION_ARG_FILENAME,     &mcff,           NULL, "mcff"},
    {"seed-length",    0, 0, G_OPTION_ARG_INT,          &seed_length,    NULL, "7"},
    {"max-mismatches", 0, 0, G_OPTION_ARG_INT,          &max_mismatches, NULL, "1"},
    {"output",         0, 0, G_OPTION_ARG_FILENAME,     &output,         NULL, "FILE"},
    {NULL}
};

static gchar
rc (gchar c)
{
    switch (c)
    {
        case 'A':
            return 'T';
        case 'C':
            return 'G';
        case 'G':
            return 'C';
        case 'T':
            return 'A';
        default:
            g_assert_not_reached ();
    }
}

static gfloat
gfloat_to_be (gfloat f)
{
    union
    {
        gfloat f;
        gint i;
    } ret;

    ret.f = f;
    ret.i = GINT_TO_BE (ret.i);

    return ret.f;
}

gint
main (gint argc, gchar **argv)
{
    g_autoptr (GOptionContext) parser = g_option_context_new ("");

    g_option_context_add_main_entries (parser,
                                       MIRBOOKING_GENERATE_SCORE_TABLE_OPTIONS,
                                       NULL);

    if (!g_option_context_parse (parser,
                                 &argc,
                                 &argv,
                                 NULL))
    {
        return EXIT_FAILURE;
    }

    if (seed_length < 1 || seed_length >= 16)
    {
        g_printerr ("The '--seed-length' argument must be between 1 and 15.\n");
        return EXIT_FAILURE;
    }

    if (max_mismatches < 0 || max_mismatches > seed_length)
    {
        g_printerr ("The '--max-mismatches' argument must be between 0 and the seed length.\n");
        return EXIT_FAILURE;
    }

    if (output == NULL)
    {
        g_printerr ("The '--output' argument is required.\n");
        return EXIT_FAILURE;
    }

    gsize l = 1l << 2 * seed_length;

    g_autoptr (FILE) f = fopen (output, "w+");

    ftruncate (fileno (f), l * l * sizeof (gfloat));

    // big-endian floats
    gfloat *table = mmap (NULL,
                          l * l * sizeof (gfloat),
                          PROT_WRITE,
                          MAP_SHARED,
                          fileno (f),
                          0);

    if (table == MAP_FAILED)
    {
        g_printerr ("Failed to map %luB of memory.\n", l*l*sizeof(gfloat));
        return EXIT_FAILURE;
    }

    gsize i,j;

    gsize completed = 0;

    g_return_val_if_fail (seed_length < 16, EXIT_FAILURE);

    #pragma omp parallel for collapse(2) schedule(static)
    for (i = 0; i < l; i++)
    {
        for (j = 0; j < l; j++)
        {
            gsize k = i * l + j;

            gchar mir_seed[16];
            gchar mre_seed[16];

            gint z;
            for (z = 0; z < seed_length; z++)
            {
                mir_seed[seed_length - z - 1] = nt[(i & (3 << 2 * z)) >> (2 * z)];
                mre_seed[seed_length - z - 1] = nt[(j & (3 << 2 * z)) >> (2 * z)];
            }

            mir_seed[seed_length] = '\0';
            mre_seed[seed_length] = '\0';

            // reverse-complement hamming distance
            gint distance = 0;
            for (z = 0; z < seed_length; z++)
            {
                distance += mir_seed[z] != rc(mre_seed[seed_length - z - 1]);
            }

            gfloat mfe;

            if (distance <= max_mismatches)
            {
                gchar *mcff_argv[] = {mcff, "-seq", mre_seed, "-zzd", mir_seed, NULL};
                gchar *standard_output;
                gchar *standard_error;
                gint exit_status;
                GError *err = NULL;
                if (!g_spawn_sync (NULL,
                                   mcff_argv,
                                   NULL,
                                   G_SPAWN_DEFAULT | G_SPAWN_SEARCH_PATH,
                                   NULL,
                                   NULL,
                                   &standard_output,
                                   &standard_error,
                                   &exit_status,
                                   &err))
                {
                    g_printerr ("%s (%s, %u).\n", err->message, g_quark_to_string (err->domain), err->code);
                    exit (EXIT_FAILURE);
                }

                if (g_spawn_check_exit_status (exit_status,
                                               &err))
                {
                    sscanf (standard_output, "%f", &mfe);
                }
                else
                {
                    g_printerr ("%s (%s, %u).\n", err->message, g_quark_to_string (err->domain), err->code);
                    g_printerr ("%s", standard_error);
                    exit (EXIT_FAILURE);
                }
            }
            else
            {
                mfe = INFINITY;
            }

            table[k] = gfloat_to_be (mfe);

            #pragma omp atomic
            ++completed;

            #pragma omp critical
            if (completed % l == 0)
            {
                g_print ("Completed %f%% (%lu/%lu)\n", 100.0 * (gdouble) completed / (gdouble) (l * l), completed, l * l);
            }
        }
    }

    munmap (table, l * l * sizeof (gfloat));

    return EXIT_SUCCESS;
}
