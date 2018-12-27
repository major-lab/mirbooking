#ifndef __MIRBOOKING_DEFAULT_SCORE_TABLE_H__
#define __MIRBOOKING_DEFAULT_SCORE_TABLE_H__

#include "mirbooking-score-table.h"

#include <gio/gio.h>

G_BEGIN_DECLS

#define MIRBOOKING_TYPE_DEFAULT_SCORE_TABLE mirbooking_default_score_table_get_type ()
G_DECLARE_FINAL_TYPE (MirbookingDefaultScoreTable, mirbooking_default_score_table, MIRBOOKING, DEFAULT_SCORE_TABLE, MirbookingScoreTable)

struct _MirbookingDefaultScoreTableClass
{
    MirbookingScoreTableClass parent_class;
};

MirbookingDefaultScoreTable * mirbooking_default_score_table_new (GBytes *seed_bytes,
                                                                  GBytes *supp_bytes);

G_END_DECLS

#endif /* __MIRBOOKING_DEFAULT_SCORE_TABLE_H__ */
