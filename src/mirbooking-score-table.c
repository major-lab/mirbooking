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
        gdouble score = mirbooking_score_table_compute_score (self, mirna, target, i, error);

        // g_return_val_if_fail (error == NULL, FALSE);

        if (score < INFINITY)
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

static gdouble
default_compute_enzymatic_score (MirbookingScoreTable *score_table,
                                 MirbookingMirna      *mirna,
                                 MirbookingTarget     *target,
                                 gsize                 position,
                                 GError              **error)
{
    /*
     * We use a fixed catalytic constant of 8.1e-4 s^-1 for clevage the
     * reaction.
     *
     * Reference:
     * Liang Meng Wee et al., “Argonaute Divides Its RNA Guide into Domains
     * with Distinct Functions and RNA-Binding Properties,” Cell 151, no. 5
     * (November 21, 2012): 1055–67,
     * https://doi.org/10.1016/j.cell.2012.10.036.
     */
    return mirbooking_score_table_compute_score (score_table, mirna, target, position, error) + (MIRBOOKING_SCORE_TABLE_DEFAULT_KCAT / MIRBOOKING_SCORE_TABLE_DEFAULT_KF);
}

void
mirbooking_score_table_class_init (MirbookingScoreTableClass *klass)
{
    klass->compute_positions       = default_compute_positions;
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
    g_return_val_if_fail (self != NULL, INFINITY);
    g_return_val_if_fail (mirna != NULL, INFINITY);
    g_return_val_if_fail (target != NULL, INFINITY);
    g_return_val_if_fail (position < mirbooking_sequence_get_sequence_length (MIRBOOKING_SEQUENCE (target)), INFINITY);

    MirbookingScoreTableClass *klass = MIRBOOKING_SCORE_TABLE_GET_CLASS (self);

    return klass->compute_score (self,
                                 mirna,
                                 target,
                                 position,
                                 error);
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
