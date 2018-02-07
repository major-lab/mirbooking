#include "mirbooking-score-table.h"

#include <math.h>

G_DEFINE_TYPE (MirbookingScoreTable, mirbooking_score_table, G_TYPE_OBJECT)

void
mirbooking_score_table_init (MirbookingScoreTable *self)
{

}

void
mirbooking_score_table_class_init (MirbookingScoreTableClass *klass)
{

}

/**
 * mirbooking_score_table_compute_score:
 * @mirna: (transfer none): A #MirbookingMirna from which a
 * reverse-complement is used against the second to compute a hybridization
 * score
 * @target: (transfer none): A #MirbookingTarget
 * @position: Identifies a MRE site in @target
 *
 * Compute the score for a pair of #MirbookingSequence objects, extracting a
 * substring of a fixed length and evaluating the score table at its
 * corresponding index.
 *
 * Returns: The corresponding score as a #gfloat or %0.0f and %error will be
 * set, unless the score was really zero.
 */
gfloat
mirbooking_score_table_compute_score (MirbookingScoreTable *self,
                                      MirbookingMirna      *mirna,
                                      MirbookingTarget     *target,
                                      gsize                 position,
                                      GError              **error)
{
    g_return_val_if_fail (self != NULL, 0.0f);
    g_return_val_if_fail (mirna != NULL, 0.0f);
    g_return_val_if_fail (target != NULL, 0.0f);

    MirbookingScoreTableClass *klass = MIRBOOKING_SCORE_TABLE_GET_CLASS (self);

    return klass->compute_score (self,
                                 mirna,
                                 target,
                                 position,
                                 error);
}

/**
 * mirbooking_score_table_compute_scores:
 * @param: (unowned) positions
 * @param: (unowned) positions_len
 *
 * Compute scores for all positions where the #mirna might be encoutered on the
 * #target. This is much more efficient than traversing all the positions of a
 * target.
 *
 * Returns: A #gfloat vector of #positions_len entries where corresponding
 * position on the target is given by the #positions out array.
 */
gfloat *
mirbooking_score_table_compute_scores (MirbookingScoreTable  *self,
                                       MirbookingMirna       *mirna,
                                       MirbookingTarget      *target,
                                       gsize                **positions,
                                       gsize                 *positions_len,
                                       GError               **error)
{
    g_return_val_if_fail (self != NULL, NULL);
    g_return_val_if_fail (mirna != NULL, NULL);
    g_return_val_if_fail (target != NULL, NULL);

    MirbookingScoreTableClass *klass = MIRBOOKING_SCORE_TABLE_GET_CLASS (self);

    return klass->compute_scores (self,
                                  mirna,
                                  target,
                                  positions,
                                  positions_len,
                                  error);
}
