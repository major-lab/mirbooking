#include "mirbooking-score-table.h"

#include <math.h>

G_DEFINE_TYPE (MirbookingScoreTable, mirbooking_score_table, G_TYPE_OBJECT)

void
mirbooking_score_table_init (MirbookingScoreTable *self)
{

}

static gdouble *
default_compute_scores (MirbookingScoreTable *self,
                        MirbookingMirna       *mirna,
                        MirbookingTarget      *target,
                        gsize                **positions,
                        gsize                 *positions_len,
                        GError               **error)
{
    gdouble *_scores   = g_new (gdouble, mirbooking_sequence_get_sequence_length (MIRBOOKING_SEQUENCE (target)));
    gsize *_positions = g_new (gsize, mirbooking_sequence_get_sequence_length (MIRBOOKING_SEQUENCE (target)));

    gint i;
    gint j = 0;
    for (i = 0; i < mirbooking_sequence_get_sequence_length (MIRBOOKING_SEQUENCE (target)); i++)
    {
        gdouble score = mirbooking_score_table_compute_score (self, mirna, target, i, error);

        g_return_val_if_fail (error != NULL, NULL);

        if (score < INFINITY)
        {
            _scores[j]    = score;
            _positions[j] = i;
            ++j;
        }
    }

    *positions = _positions;
    *positions_len = j;

    return _scores;
}

static gdouble
default_compute_enzymatic_score (MirbookingScoreTable *score_table,
                                 MirbookingMirna      *mirna,
                                 MirbookingTarget     *target,
                                 gsize                 position,
                                 GError              **error)
{
    // Reference: https://www.uniprot.org/uniprot/Q9UKV8
    // return 1.1;
    return mirbooking_score_table_compute_score (score_table, mirna, target, position, error);
}

void
mirbooking_score_table_class_init (MirbookingScoreTableClass *klass)
{
    klass->compute_scores          = default_compute_scores;
    klass->compute_enzymatic_score = default_compute_enzymatic_score;
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
 * The score correspond to the dissociation constant (Kd) expressed in
 * nanomolars units (nM).
 *
 * Returns: The corresponding score as a #gdouble or %0.0f and %error will be
 * set, unless the score was really zero.
 */
gdouble
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
 * @mirna:
 * @target:
 * @positions: (array length=positions_len) (out): A vector if positions where
 * a score exists
 * @positions_len: (out): Length of @positions
 *
 * Compute scores for all positions where the #mirna might be encoutered on the
 * #target. This is much more efficient than traversing all the positions of a
 * target.
 *
 * Returns: (array length=positions_len) (transfer full): A #gfloat vector of
 * #positions_len entries where corresponding position on the target is given
 * by the #positions out array.
 */
gdouble *
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

/**
 * mirbooking_score_table_compute_enzymatic_score:
 *
 * Compute an enzymatic score, which corresponds to the catalytic efficiency
 * (Km) expressed in nanomolar unit (nM).
 */
gdouble
mirbooking_score_table_compute_enzymatic_score (MirbookingScoreTable *self,
                                                MirbookingMirna      *mirna,
                                                MirbookingTarget     *target,
                                                gsize                 position,
                                                GError              **error)
{
    g_return_val_if_fail (self != NULL, 0.0f);
    g_return_val_if_fail (mirna != NULL, 0.0f);
    g_return_val_if_fail (target != NULL, 0.0f);

    MirbookingScoreTableClass *klass = MIRBOOKING_SCORE_TABLE_GET_CLASS (self);

    return klass->compute_enzymatic_score (self,
                                           mirna,
                                           target,
                                           position,
                                           error);
}
