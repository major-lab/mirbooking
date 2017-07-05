#include "mirbooking-score-table.h"

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
 * @a: (transfer none): The first #MirbookingSequence from which a
 * reverse-complement is used against the second to compute a hybridization
 * score; it is normally the miRNA seed
 * @a_offset: Offset in #a, which is normally %1 for the seed
 * @b: (transfer none): The second #MirbookingSequence which is normally the
 * target MRE site
 * @b_offset: Offset in #b which identify a specific target site
 * @len: Length of the subsequence considered for computing the hybridization
 * score
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
                                      MirbookingSequence   *a,
                                      gsize                 a_offset,
                                      MirbookingSequence   *b,
                                      gsize                 b_offset,
                                      gsize                 len,
                                      GError              **error)
{
    g_return_val_if_fail (self != NULL, 0.0f);
    g_return_val_if_fail (a != NULL, 0.0f);
    g_return_val_if_fail (b != NULL, 0.0f);

    MirbookingScoreTableClass *klass = MIRBOOKING_SCORE_TABLE_GET_CLASS (self);

    return klass->compute_score (self,
                                 a,
                                 a_offset,
                                 b,
                                 b_offset,
                                 len,
                                 error);
}
