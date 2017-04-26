#ifndef __MIRBOOKING_SCORE_TABLE_H__
#define __MIRBOOKING_SCORE_TABLE_H__

#include "mirbooking-sequence.h"

#include <glib-object.h>

G_BEGIN_DECLS

#define MIRBOOKING_TYPE_SCORE_TABLE mirbooking_score_table_get_type ()
G_DECLARE_FINAL_TYPE (MirbookingScoreTable, mirbooking_score_table, MIRBOOKING, SCORE_TABLE, GObject)

struct _MirbookingScoreTableClass
{
    GObject parent_class;
};

/* to decide later... */
MirbookingScoreTable * mirbooking_score_table_new             (void);
MirbookingScoreTable * mirbooking_score_table_new_precomputed (GBytes *bytes);
gfloat                 mirbooking_score_table_compute_score   (MirbookingScoreTable *self,
                                                               MirbookingSequence   *a,
                                                               gsize                 a_offset,
                                                               MirbookingSequence   *b,
                                                               gsize                 b_offset,
                                                               gsize                 len);

G_END_DECLS

#endif /* __MIRBOOKING_SCORE_TABLE_H__ */
