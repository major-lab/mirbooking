#include <glib.h>
#include <gio/gio.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <string.h>
#if HAVE_OPENMP
#include <omp.h>
#else
static int omp_get_thread_num (void) { return 0; }
#endif

#include <sparse.h>

G_DEFINE_AUTOPTR_CLEANUP_FUNC (FILE, fclose)

static gchar  *method      = "RNAcofold";
static gdouble temperature = 310.5;
static gchar  *mask        = "||||...";
static gchar  *hard_mask   = NULL; /* defaults to mask */
static gchar  *output      = NULL;

static gchar *nt = "ACGT";

static gfloat (*fold_duplex) (gchar *a, gchar *b, gchar *a_mask, gchar *b_mask, GError **err);

static const GOptionEntry MIRBOOKING_GENERATE_SCORE_TABLE_OPTIONS[] =
{
    {"method",      0, 0, G_OPTION_ARG_STRING,   &method,      NULL, "RNAcofold"},
    {"temperature", 0, 0, G_OPTION_ARG_DOUBLE,   &temperature, NULL, "310.5"},
    {"mask",        0, 0, G_OPTION_ARG_STRING,   &mask,        NULL, "||||..."},
    {"hard-mask",   0, 0, G_OPTION_ARG_STRING,   &hard_mask,   NULL, "||||..."},
    {"output",      0, 0, G_OPTION_ARG_FILENAME, &output,      NULL, "FILE"},
    {0}
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

static gfloat
fold_duplex_RNAcofold (gchar *a, gchar *b, gchar *a_mask, gchar *b_mask, GError **err)
{
    gfloat binding_energy;
    g_autofree gchar *standard_input = g_strdup_printf ("%s&%s\n%s&%s", a, b, a_mask, b_mask);
    g_autofree gchar *standard_output;
    g_autofree gchar *standard_error;
    g_autofree gchar *temperature_str = g_strdup_printf ("%f", temperature - 272.15); // Kelvin -> Celsius

    g_autoptr (GRegex) vienna_binding_energy_regex = g_regex_new ("delta G binding=\\s?(.+)", 0, 0, NULL);

    g_autoptr (GSubprocess) proc = g_subprocess_new (G_SUBPROCESS_FLAGS_STDIN_PIPE | G_SUBPROCESS_FLAGS_STDOUT_PIPE | G_SUBPROCESS_FLAGS_STDERR_PIPE,
                                                     err,
                                                     "RNAcofold", "--noPS", "-p", "-C", "-T", temperature_str, NULL);

    if (!proc)
    {
        return NAN;
    }

    if (!g_subprocess_communicate_utf8 (proc,
                                        standard_input,
                                        NULL,
                                        &standard_output,
                                        &standard_error,
                                        err))
    {
        return NAN;
    }

    if (!g_subprocess_wait_check (proc,
                                  NULL,
                                  err))
    {
        return NAN;
    }

    g_autoptr (GMatchInfo) match_info;
    if (!g_regex_match (vienna_binding_energy_regex, standard_output, 0, &match_info))
    {
        return NAN;
    }

    g_autofree gchar *binding_energy_str = g_match_info_fetch (match_info, 1);

    sscanf (binding_energy_str, "%f", &binding_energy);

    return binding_energy;
}

static gfloat
fold_duplex_mcff (gchar *a, gchar *b, gchar *a_mask, gchar *b_mask, GError **err)
{
    g_autofree gchar *standard_output;
    gfloat mfe;
    g_autofree gchar *mask = g_strdup_printf ("%sxx%s", a_mask, b_mask);

    gint i;
    for (i = 0; i < strlen (mask); i++)
    {
        switch (mask[i])
        {
            case 'x':
                mask[i] = '.';
                break;
            case '.':
                mask[i] = 'x';
                break;
        }
    }

    g_autoptr (GSubprocess) proc = g_subprocess_new (G_SUBPROCESS_FLAGS_STDOUT_PIPE,
                                                     err,
                                                     "mcff", "-seq", a, "-sd", b, "-mask", mask, NULL);

    if (!proc)
    {
        return NAN;
    }

    if (!g_subprocess_communicate_utf8 (proc,
                                        NULL,
                                        NULL,
                                        &standard_output,
                                        NULL,
                                        err))
    {
        return NAN;
    }

    if (!g_subprocess_wait_check (proc,
                                  NULL,
                                  err))
    {
        return NAN;
    }

    sscanf (standard_output, "%f", &mfe);

    return mfe;
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

    if (g_strcmp0 (method, "RNAcofold") == 0)
    {
        fold_duplex = fold_duplex_RNAcofold;

    }
    else if (g_strcmp0 (method, "mcff") == 0)
    {
        fold_duplex = fold_duplex_mcff;
    }
    else
    {
        g_printerr ("Unknown folding method '%s'.\n", method);
        return EXIT_FAILURE;
    }

    gsize seed_length = strlen (mask);

    if (seed_length < 2 || seed_length >= 16)
    {
        g_printerr ("The mask must be formed of at least 2 nucleotides.\n");
        return EXIT_FAILURE;
    }

    if (!hard_mask)
    {
        hard_mask = mask;
    }

    if (strlen (mask) != strlen (hard_mask))
    {
        g_printerr ("The hard mask must be the same size as the mask.");
        return EXIT_FAILURE;
    }

    if (output == NULL)
    {
        g_printerr ("The '--output' argument is required.\n");
        return EXIT_FAILURE;
    }

    // convert the mask into a symmetric mask for folding
    g_autofree gchar *mre_mask = g_new0 (gchar, seed_length + 1);
    g_autofree gchar *mir_mask = g_new0 (gchar, seed_length + 1);
    {
        gint i;
        for (i = 0; i < seed_length; i++)
        {
            mir_mask[i] = (mask[i] == '|') ? ')' : mask[i];
            mre_mask[seed_length - i - 1] = (mask[i] == '|') ? '(' : mask[i];
        }
    }

    gsize n   = pow4 (seed_length);
    gsize nnz = 1;
    gsize z;
    for (z = 0; z < seed_length; z++)
    {
        switch (hard_mask[z])
        {
            case '|':
                nnz *= 4; /* canonical match */
                break;
            case 'x':
                nnz *= 12; /* canonical mismatch */
                break;
            case '.':
                nnz *= 16; /* no constraint */
                break;
            default:
                g_printerr ("Unknown symbol '%c' in mask at position %lu.", hard_mask[z], z);
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

    // shorthand to avoid aliasing
    gfloat* data = (gfloat*) (table + 2 + n + 1 + nnz);

    SparseMatrix sm;
    sm.storage        = SPARSE_MATRIX_STORAGE_CSR;
    sm.type           = SPARSE_MATRIX_TYPE_FLOAT;
    sm.shape[0]       = n;
    sm.shape[1]       = n;
    sm.s.csr.nnz      = nnz;
    sm.s.csr.rowptr   = table + 2;
    sm.s.csr.colind   = table + 2 + n + 1;
    sm.default_data.f = INFINITY;
    sm.data           = data;

    table[0] = n;
    table[1] = nnz;

    gsize i;

    // all rows are equally-sized
    #pragma omp parallel for
    for (i = 0; i <= n; i++)
    {
        sm.s.csr.rowptr[i] = i * (nnz / n);
    }

    guint64 begin = g_get_monotonic_time ();

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

            gsize z;
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
                switch (hard_mask[z])
                {
                    case '|':
                        distance += mir_seed[z] != rc(mre_seed[seed_length - z - 1]);
                        break;
                    case 'x':
                        distance += mir_seed[z] == rc(mre_seed[seed_length - z - 1]);
                        break;
                }
            }

            if (distance <= 0)
            {
                g_autoptr (GError) err = NULL;
                gfloat binding_energy = fold_duplex (mre_seed,
                                                     mir_seed,
                                                     mre_mask,
                                                     mir_mask,
                                                     &err);

                if (binding_energy == NAN)
                {
                    g_printerr ("%s (%s, %u).\n", err->message, g_quark_to_string (err->domain), err->code);
                    exit (EXIT_FAILURE);
                }

                sm.s.csr.colind[sm.s.csr.rowptr[i] + k] = j;
                data[sm.s.csr.rowptr[i] + k]            = binding_energy;

                ++k;

                #pragma omp atomic
                ++completed;
            }
        }

        if (omp_get_thread_num () == 0)
        {
            g_print ("\r%.2f%% %lu/%lu [%.2fit/sec]",
                     100.0 * (gdouble) completed / (gdouble) nnz,
                     completed,
                     nnz,
                     completed / ((gdouble) (g_get_monotonic_time () - begin) / (gdouble) G_USEC_PER_SEC));
        }
    }

    munmap (table, file_len);

    return EXIT_SUCCESS;
}
