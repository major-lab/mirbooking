#ifndef __MIRBOOKING_MCFF_SCORE_TABLE_H__
#define __MIRBOOKING_MCFF_SCORE_TABLE_H__

#include "mirbooking-score-table.h"

#include <gio/gio.h>

G_BEGIN_DECLS
#define MIRBOOKING_TYPE_MCFF_SCORE_TABLE mirbooking_mcff_score_table_get_type ()
G_DECLARE_FINAL_TYPE (MirbookingMcffScoreTable, mirbooking_mcff_score_table, MIRBOOKING, MCFF_SCORE_TABLE, MirbookingScoreTable)

struct _MirbookingMcffScoreTableClass
{
    MirbookingScoreTableClass parent_class;
};

MirbookingMcffScoreTable * mirbooking_mcff_score_table_new (void);

G_END_DECLS

#endif /* __MIRBOOKING_MCFF_SCORE_TABLE_H__ */
