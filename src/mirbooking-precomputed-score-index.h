#ifndef __MIRBOOKING_PRECOMPUTED_SCORE_INDEX_H__
#define __MIRBOOKING_PRECOMPUTED_SCORE_INDEX_H__

#include "mirbooking-score-index.h"

G_BEGIN_DECLS

#define MIRBOOKING_TYPE_PRECOMPUTED_SCORE_INDEX mirbooking_precomputed_score_index_get_type ()
G_DECLARE_FINAL_TYPE (MirbookingPrecomputedScoreIndex, mirbooking_precomputed_score_index, MIRBOOKING, PRECOMPUTED_SCORE_INDEX, MirbookingScoreIndex)

#define MIRBOOKING_TYPE_PRECOMPUTED_SCORE_INDEX_ITER mirbooking_precomputed_score_index_iter_get_type ()
G_DECLARE_FINAL_TYPE (MirbookingPrecomputedScoreIndexIter, mirbooking_precomputed_score_index_iter, MIRBOOKING, PRECOMPUTED_SCORE_INDEX_ITER, MirbookingScoreIndexIter)

MirbookingPrecomputedScoreIndex * mirbooking_precomputed_score_index_new            (guint16 *data);
MirbookingPrecomputedScoreIndex * mirbooking_precomputed_score_index_new_from_bytes (GBytes *data);

G_END_DECLS

#endif /* __MIRBOOKING_PRECOMPUTED_SCORE_INDEX_H__ */
