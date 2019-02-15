#include "mirbooking-mcff-score-table.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#define R 1.987203611e-3
#define T 310.15

#define KF   2.4e-4 // pM^-1s^-1
#define KCAT 3.6e-2 // s^-1

/*
 * See @MirbookingDefaultScoreTable for the detail of the computation. Here,
 * mcff returned -19.765 kcal/mol.
 */
#define AGO2_SCORE (4.40f)

struct _MirbookingMcffScoreTable
{
    MirbookingScoreTable parent_instance;
};

G_DEFINE_TYPE (MirbookingMcffScoreTable, mirbooking_mcff_score_table, MIRBOOKING_TYPE_SCORE_TABLE)

static gboolean
compute_score (MirbookingScoreTable  *score_table,
               MirbookingMirna       *mirna,
               MirbookingTarget      *target,
               gsize                  position,
               MirbookingScore       *score,
               GError               **error)
{
    gchar mirna_seq[9]  = {0};
    gchar target_seq[9] = {0};

    memcpy (mirna_seq + 1, mirbooking_sequence_get_subsequence (MIRBOOKING_SEQUENCE (mirna), 1, 7), 7);
    memcpy (target_seq, mirbooking_sequence_get_subsequence (MIRBOOKING_SEQUENCE (target), position, position + 7), 7);

    // dangling end
    target_seq[7] = 'A';
    mirna_seq[0]  = 'A';

    gchar *argv[] = {"mcff", "-seq", target_seq, "-zzd", mirna_seq, NULL};

    g_autofree gchar *standard_output = NULL, *standard_error = NULL;
    gint exit_status;

    if (!g_spawn_sync (NULL,
                       argv,
                       NULL,
                       G_SPAWN_SEARCH_PATH,
                       NULL,
                       NULL,
                       &standard_output,
                       &standard_error,
                       &exit_status,
                       error))
    {
        return FALSE;
    }

    if (!g_spawn_check_exit_status (exit_status,
                                    error))
    {
        return FALSE;
    }

    gfloat mfe;
    sscanf (standard_output, "%f", &mfe);

    MirbookingScore ret = {KF, KF * (1e12 * exp ((mfe + AGO2_SCORE) / (R * T))), KCAT};
    *score = ret;

    return TRUE;
}

static void
mirbooking_mcff_score_table_init (MirbookingMcffScoreTable *self)
{
}

static void
mirbooking_mcff_score_table_class_init (MirbookingMcffScoreTableClass *klass)
{
    klass->parent_class.compute_score  = compute_score;
}

MirbookingMcffScoreTable *
mirbooking_mcff_score_table_new (void)
{
    return g_object_new (MIRBOOKING_TYPE_MCFF_SCORE_TABLE, NULL);
}
