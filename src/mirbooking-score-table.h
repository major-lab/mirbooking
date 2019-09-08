#ifndef __MIRBOOKING_SCORE_TABLE_H__
#define __MIRBOOKING_SCORE_TABLE_H__

#include "mirbooking-error.h"
#include "mirbooking-target.h"
#include "mirbooking-mirna.h"

#include <glib-object.h>

G_BEGIN_DECLS

#define MIRBOOKING_TYPE_SCORE_TABLE mirbooking_score_table_get_type ()
G_DECLARE_DERIVABLE_TYPE (MirbookingScoreTable, mirbooking_score_table, MIRBOOKING, SCORE_TABLE, GObject)

typedef struct _MirbookingScore
{
    gdouble kf;
    gdouble kr;
    gdouble kcleave;
    gdouble krelease;
    gdouble kcat;
} MirbookingScore;

#define MIRBOOKING_SCORE_KD(score) (score.kr / score.kf)
#define MIRBOOKING_SCORE_KM(score) ((score.kr + score.kcat) / score.kf)

struct _MirbookingScoreTableClass
{
    GObjectClass parent_class;

    gboolean (*compute_positions) (MirbookingScoreTable  *self,
                                   MirbookingMirna       *mirna,
                                   MirbookingTarget      *target,
                                   gsize                **positions,
                                   gsize                 *positions_len,
                                   GError               **error);

    gboolean (*compute_score)     (MirbookingScoreTable  *self,
                                   MirbookingMirna       *mirna,
                                   MirbookingTarget      *target,
                                   gsize                  position,
                                   MirbookingScore       *score,
                                   GError               **error);
};

gboolean mirbooking_score_table_compute_positions  (MirbookingScoreTable  *self,
                                                    MirbookingMirna       *mirna,
                                                    MirbookingTarget      *target,
                                                    gsize                **positions,
                                                    gsize                 *positions_len,
                                                    GError               **error);

gboolean mirbooking_score_table_compute_score (MirbookingScoreTable  *self,
                                               MirbookingMirna       *mirna,
                                               MirbookingTarget      *target,
                                               gsize                  position,
                                               MirbookingScore       *score,
                                               GError               **error);

G_END_DECLS

#endif /* __MIRBOOKING_SCORE_TABLE_H__ */
