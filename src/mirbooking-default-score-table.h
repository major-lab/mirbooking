#ifndef __MIRBOOKING_DEFAULT_SCORE_TABLE_H__
#define __MIRBOOKING_DEFAULT_SCORE_TABLE_H__

#include "mirbooking-score-table.h"

#include <gio/gio.h>

G_BEGIN_DECLS

#define MIRBOOKING_DEFAULT_SCORE_TABLE_DEFAULT_SEED_OFFSET 1
#define MIRBOOKING_DEFAULT_SCORE_TABLE_DEFAULT_SEED_LENGTH 7

#define MIRBOOKING_TYPE_DEFAULT_SCORE_TABLE mirbooking_default_score_table_get_type ()
G_DECLARE_FINAL_TYPE (MirbookingDefaultScoreTable, mirbooking_default_score_table, MIRBOOKING, DEFAULT_SCORE_TABLE, MirbookingScoreTable)

struct _MirbookingDefaultScoreTableClass
{
    MirbookingScoreTableClass parent_class;
};

MirbookingDefaultScoreTable * mirbooking_default_score_table_new             (gfloat *data,
                                                                              gsize   seed_offset,
                                                                              gsize   seed_len);
MirbookingDefaultScoreTable * mirbooking_default_score_table_new_from_bytes  (GBytes *bytes,
                                                                              gsize   seed_offset,
                                                                              gsize   seed_len,
                                                                              GBytes *supp_bytes);
MirbookingDefaultScoreTable * mirbooking_default_score_table_new_from_stream (GInputStream  *stream,
                                                                              gsize          seed_offset,
                                                                              gsize          seed_len,
                                                                              GError       **error);

G_END_DECLS

#endif /* __MIRBOOKING_DEFAULT_SCORE_TABLE_H__ */
