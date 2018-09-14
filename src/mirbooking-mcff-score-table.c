#include "mirbooking-mcff-score-table.h"

#include <math.h>
#include <stdio.h>

#define R 1.987203611e-3
#define T 310.15

struct _MirbookingMcffScoreTable
{
    MirbookingScoreTable parent_instance;
};

G_DEFINE_TYPE (MirbookingMcffScoreTable, mirbooking_mcff_score_table, MIRBOOKING_TYPE_SCORE_TABLE);

static gdouble
compute_score (MirbookingScoreTable *score_table,
               MirbookingMirna      *mirna,
               MirbookingTarget     *target,
               gsize                 position,
               GError              **error)
{
    gchar mirna_seq[8] = {0};
    gchar target_seq[8] = {0};

    strncpy (mirna_seq, mirbooking_sequence_get_subsequence (MIRBOOKING_SEQUENCE (mirna), 1, 7), 7);
    strncpy (target_seq, mirbooking_sequence_get_subsequence (MIRBOOKING_SEQUENCE (target), position, position + 7), 7);

    gchar *argv[] = {MCFF, "-seq", target_seq, "-zzd", mirna_seq, NULL};

    gchar *standard_output, *standard_error;
    gint exit_status;

    if (!g_spawn_sync (NULL,
                       argv,
                       NULL,
                       G_SPAWN_DEFAULT,
                       NULL,
                       NULL,
                       &standard_output,
                       &standard_error,
                       &exit_status,
                       error))
    {
        return INFINITY;
    }

    gfloat mfe;
    sscanf (standard_output, "%f", &mfe);

    return 1e9 * exp (mfe / (R * T));
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
