#include "mirbooking-score-table.h"

#include <math.h>

G_DEFINE_TYPE (MirbookingScoreTable, mirbooking_score_table, G_TYPE_OBJECT)

void
mirbooking_score_table_init (MirbookingScoreTable *self)
{

}

static gboolean
default_compute_positions (MirbookingScoreTable *self,
                           MirbookingMirna       *mirna,
                           MirbookingTarget      *target,
                           gsize                **positions,
                           gsize                 *positions_len,
                           GError               **error)
{
    g_autofree gsize *_positions = NULL;

    gsize i;
    gsize j = 0;
    for (i = 0; i < mirbooking_sequence_get_sequence_length (MIRBOOKING_SEQUENCE (target)); i++)
    {
        MirbookingScore score;
        if (!mirbooking_score_table_compute_score (self, mirna, target, i, &score, error))
        {
            return FALSE;
        }

        if (MIRBOOKING_SCORE_IS_FINITE (score))
        {
            _positions = g_realloc (_positions, (j + 1) * sizeof (gsize));
            _positions[j] = i;
            ++j;
        }
    }

    *positions = g_steal_pointer (&_positions);
    *positions_len = j;

    return TRUE;
}

void
mirbooking_score_table_class_init (MirbookingScoreTableClass *klass)
{
    klass->compute_positions = default_compute_positions;
}

/**
 * mirbooking_score_table_compute_positions:
 * @mirna:
 * @target:
 * @positions: (array length=positions_len) (out): A vector if positions where
 * a score exists
 * @positions_len: (out): Length of @positions
 *
 * Compute scores for all positions where the #mirna might be encoutered on the
 * #target. This is much more efficient than traversing all the positions of a
 * target.
 */
gboolean
mirbooking_score_table_compute_positions (MirbookingScoreTable  *self,
                                          MirbookingMirna       *mirna,
                                          MirbookingTarget      *target,
                                          gsize                **positions,
                                          gsize                 *positions_len,
                                          GError               **error)
{
    g_return_val_if_fail (self != NULL, FALSE);
    g_return_val_if_fail (mirna != NULL, FALSE);
    g_return_val_if_fail (target != NULL, FALSE);

    MirbookingScoreTableClass *klass = MIRBOOKING_SCORE_TABLE_GET_CLASS (self);

    return klass->compute_positions (self,
                                     mirna,
                                     target,
                                     positions,
                                     positions_len,
                                     error);
}

/**
 * mirbooking_score_table_compute_score:
 * @mirna: (transfer none): A #MirbookingMirna from which a
 * reverse-complement is used against the second to compute a hybridization
 * score
 * @target: (transfer none): A #MirbookingTarget
 * @position: Identifies a MRE site in @target
 * @score: (out):
 *
 * Compute the score for a pair of #MirbookingSequence objects, extracting a
 * substring of a fixed length and evaluating the score table at its
 * corresponding index.
 *
 * Returns: %TRUE on success otherwise %FALSE and @error is set.
 */
gboolean
mirbooking_score_table_compute_score (MirbookingScoreTable  *self,
                                      MirbookingMirna       *mirna,
                                      MirbookingTarget      *target,
                                      gsize                  position,
                                      MirbookingScore       *score,
                                      GError               **error)
{
    g_return_val_if_fail (self != NULL, FALSE);
    g_return_val_if_fail (mirna != NULL, FALSE);
    g_return_val_if_fail (target != NULL, FALSE);
    g_return_val_if_fail (position < mirbooking_sequence_get_sequence_length (MIRBOOKING_SEQUENCE (target)), FALSE);
    g_return_val_if_fail (score != NULL, FALSE);

    MirbookingScoreTableClass *klass = MIRBOOKING_SCORE_TABLE_GET_CLASS (self);

    return klass->compute_score (self,
                                 mirna,
                                 target,
                                 position,
                                 score,
                                 error);
}
