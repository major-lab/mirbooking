#ifndef __MIRBOOKING_PRECOMPUTED_SCORE_TABLE_H__
#define __MIRBOOKING_PRECOMPUTED_SCORE_TABLE_H__

#include "mirbooking-score-table.h"

#include <gio/gio.h>

G_BEGIN_DECLS

#define MIRBOOKING_TYPE_PRECOMPUTED_SCORE_TABLE mirbooking_precomputed_score_table_get_type ()
G_DECLARE_FINAL_TYPE (MirbookingPrecomputedScoreTable, mirbooking_precomputed_score_table, MIRBOOKING, PRECOMPUTED_SCORE_TABLE, MirbookingScoreTable)

struct _MirbookingPrecomputedScoreTableClass
{
    MirbookingScoreTableClass parent_class;
};

MirbookingPrecomputedScoreTable * mirbooking_precomputed_score_table_new             (gfloat *data);
MirbookingPrecomputedScoreTable * mirbooking_precomputed_score_table_new_from_bytes  (GBytes *bytes);
MirbookingPrecomputedScoreTable * mirbooking_precomputed_score_table_new_from_stream (GInputStream  *stream,
                                                                                      GError       **error);

G_END_DECLS

#endif /* __MIRBOOKING_PRECOMPUTED_SCORE_TABLE_H__ */
