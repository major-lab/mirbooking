#ifndef __MIRBOOKING_SCORE_TABLE_H__
#define __MIRBOOKING_SCORE_TABLE_H__

#include "mirbooking-error.h"
#include "mirbooking-target.h"
#include "mirbooking-mirna.h"

#include <glib-object.h>

G_BEGIN_DECLS

#define MIRBOOKING_TYPE_SCORE_TABLE mirbooking_score_table_get_type ()
G_DECLARE_DERIVABLE_TYPE (MirbookingScoreTable, mirbooking_score_table, MIRBOOKING, SCORE_TABLE, GObject)

struct _MirbookingScoreTableClass
{
    GObjectClass parent_class;

    gdouble (*compute_score)    (MirbookingScoreTable *self,
                                 MirbookingMirna      *mirna,
                                 MirbookingTarget     *target,
                                 gsize                 position,
                                 GError              **error);

    gdouble * (*compute_scores) (MirbookingScoreTable  *self,
                                 MirbookingMirna       *mirna,
                                 MirbookingTarget      *target,
                                 gsize                **positions,
                                 gsize                 *positions_len,
                                 GError               **error);

    gdouble (*compute_enzymatic_score) (MirbookingScoreTable *self,
                                        MirbookingMirna      *mirna,
                                        MirbookingTarget     *target,
                                        gsize                 position,
                                        GError              **error);

};

gdouble   mirbooking_score_table_compute_score  (MirbookingScoreTable *self,
                                                MirbookingMirna      *mirna,
                                                MirbookingTarget     *target,
                                                gsize                 position,
                                                GError              **error);

gdouble * mirbooking_score_table_compute_scores (MirbookingScoreTable  *self,
                                                MirbookingMirna       *mirna,
                                                MirbookingTarget      *target,
                                                gsize                **positions,
                                                gsize                 *positions_len,
                                                GError               **error);

gdouble mirbooking_score_table_compute_enzymatic_score (MirbookingScoreTable *self,
                                                         MirbookingMirna      *mirna,
                                                         MirbookingTarget     *target,
                                                         gsize                 position,
                                                         GError              **error);

G_END_DECLS

#endif /* __MIRBOOKING_SCORE_TABLE_H__ */
