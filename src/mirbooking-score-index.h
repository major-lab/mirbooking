#ifndef __MIRBOKING_SCORE_INDEX_H__
#define __MIRBOKING_SCORE_INDEX_H__

#include <glib.h>

#include "mirbooking-mirna.h"
#include "mirbooking-target.h"

G_BEGIN_DECLS

#define MIRBOOKING_TYPE_SCORE_INDEX mirbooking_score_index_get_type ()
G_DECLARE_DERIVABLE_TYPE (MirbookingScoreIndex, mirbooking_score_index, MIRBOOKING, SCORE_INDEX, GObject)

#define MIRBOOKING_TYPE_SCORE_INDEX_ITER mirbooking_score_index_iter_get_type ()
G_DECLARE_DERIVABLE_TYPE (MirbookingScoreIndexIter, mirbooking_score_index_iter, MIRBOOKING, SCORE_INDEX_ITER, GObject)

struct _MirbookingScoreIndexClass
{
    GObjectClass parent_class;
    void                       (*set_sequence_quantity) (MirbookingScoreIndex *self,
                                                         MirbookingSequence   *sequence,
                                                         gfloat                quantity);
    MirbookingScoreIndexIter * (*iterator)     (MirbookingScoreIndex *self);
};

struct _MirbookingScoreIndexIterClass
{
    GObjectClass parent_class;
    gboolean (*next) (MirbookingScoreIndexIter  *self,
                      MirbookingMirna          **mirna,
                      MirbookingTarget         **target,
                      gsize                     *position);
};

void                       mirbooking_score_index_set_sequence_quantity (MirbookingScoreIndex *self,
                                                                         MirbookingSequence   *sequence,
                                                                         gfloat                quantity);
MirbookingScoreIndexIter * mirbooking_score_index_iterator              (MirbookingScoreIndex *self);
MirbookingScoreIndex *     mirbooking_score_index_iter_get_score_index  (MirbookingScoreIndexIter *self);
gboolean                   mirbooking_score_index_iter_next             (MirbookingScoreIndexIter *self,
                                                                         MirbookingMirna         **a,
                                                                         MirbookingTarget        **b,
                                                                         gsize                    *position);


G_END_DECLS

#endif /* __MIRBOKING_SCORE_INDEX_H__ */
