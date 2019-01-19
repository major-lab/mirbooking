#include <glib.h>
#include <gio/gio.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <string.h>

#include <sparse.h>

G_DEFINE_AUTOPTR_CLEANUP_FUNC (FILE, fclose)

static gchar  *mcff   = "mcff";
static gchar  *mask   = NULL;
static gchar  *output = NULL;

static gchar *nt = "ACGT";

static const GOptionEntry MIRBOOKING_GENERATE_SCORE_TABLE_OPTIONS[] =
{
    {"mcff",   0, 0, G_OPTION_ARG_FILENAME, &mcff,   NULL, "mcff"},
    {"mask",   0, 0, G_OPTION_ARG_STRING,   &mask,   NULL, "....xxx"},
    {"output", 0, 0, G_OPTION_ARG_FILENAME, &output, NULL, "FILE"},
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

/**
 * Compute the fourth integer power efficiently using bit-shift.
 */
static gsize
pow4 (gsize k)
{
    return 1l << 2 * k;
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

    gsize seed_length = strlen (mask);

    if (seed_length < 2 || seed_length >= 16)
    {
        g_printerr ("The mask must be formed of at least 2 nucleotides.\n");
        return EXIT_FAILURE;
    }

    if (output == NULL)
    {
        g_printerr ("The '--output' argument is required.\n");
        return EXIT_FAILURE;
    }

    gsize n   = pow4 (seed_length);
    gsize nnz = 1;
    gint z;
    for (z = 0; z < seed_length; z++)
    {
        switch (mask[z])
        {
            case '.':
                nnz *= 4;
                break;
            case 'x':
                nnz *= 16;
                break;
            default:
                g_printerr ("Unknown symbol '%c' in mask at position %d.", mask[z], z);
                return EXIT_FAILURE;
        }
    }

    g_debug ("n: %lu nnz %lu\n", n, nnz);

    gsize file_len = sizeof (gsize) + sizeof (gsize) + (n + 1) * sizeof (gsize) + + nnz * sizeof (gsize) + nnz * sizeof (gfloat);

    g_autoptr (FILE) f = fopen (output, "w+");

    ftruncate (fileno (f), file_len);

    gsize *table = mmap (NULL,
                         file_len,
                         PROT_WRITE,
                         MAP_SHARED,
                         fileno (f),
                         0);

    if (table == MAP_FAILED)
    {
        g_printerr ("Failed to map %luB of memory.\n", file_len);
        return EXIT_FAILURE;
    }

    g_return_val_if_fail (seed_length < 16, EXIT_FAILURE);

    SparseMatrix sm;
    sm.storage        = SPARSE_MATRIX_STORAGE_CSR;
    sm.type           = SPARSE_MATRIX_TYPE_FLOAT;
    sm.shape[0]       = n;
    sm.shape[1]       = n;
    sm.s.csr.nnz      = nnz;
    sm.s.csr.rowptr   = table + 2;
    sm.s.csr.colind   = table + 2 + n + 1;
    sm.default_data.f = INFINITY;
    sm.data           = (gfloat*) (table + 2 + n + 1 + nnz);

    table[0] = n;
    table[1] = nnz;

    gsize i;

    // all rows are equally-sized
    #pragma omp parallel for
    for (i = 0; i <= n; i++)
    {
        sm.s.csr.rowptr[i] = i * (nnz / n);
    }

    gsize completed = 0;
    #pragma omp parallel for
    for (i = 0; i < n; i++)
    {
        gsize j;
        gsize k = 0;
        for (j = 0; j < n; j++)
        {
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

            // reverse-complement hamming distance for the seed for paired
            // positions
            gint distance = 0;
            for (z = 0; z < seed_length; z++)
            {
                if (mask[z] == '.')
                {
                    distance += mir_seed[z] != rc(mre_seed[seed_length - z - 1]);
                }
            }

            if (distance == 0)
            {
                gfloat mfe;
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

                sm.s.csr.colind[sm.s.csr.rowptr[i] + k]     = j;
                ((gfloat*)sm.data) [sm.s.csr.rowptr[i] + k] = mfe;

                ++k;
            }

            #pragma omp atomic
            ++completed;

            if (completed % n == 0)
            {
                #pragma omp critical
                g_print ("Completed %f%% (%lu/%lu)\n", 100.0 * (gdouble) completed / pow (n, 2), completed, n * n);
            }
        }
    }

    munmap (table, file_len);

    return EXIT_SUCCESS;
}
