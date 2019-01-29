#ifndef __MIRBOOKING_DEFAULT_SCORE_TABLE_H__
#define __MIRBOOKING_DEFAULT_SCORE_TABLE_H__

#include "mirbooking-score-table.h"

#include <gio/gio.h>

G_BEGIN_DECLS

#define MIRBOOKING_DEFAULT_SCORE_TABLE_DEFAULT_SUPPLEMENTARY_MODEL MIRBOOKING_DEFAULT_SCORE_TABLE_SUPPLEMENTARY_MODEL_YAN_ET_AL_2018

typedef enum _MirbookingDefaultScoreTableSupplementaryModel
{
    MIRBOOKING_DEFAULT_SCORE_TABLE_SUPPLEMENTARY_MODEL_ZAMORE_ET_AL_2012,
    MIRBOOKING_DEFAULT_SCORE_TABLE_SUPPLEMENTARY_MODEL_YAN_ET_AL_2018
} MirbookingDefaultScoreTableSupplementaryModel;

#define MIRBOOKING_TYPE_DEFAULT_SCORE_TABLE mirbooking_default_score_table_get_type ()
G_DECLARE_FINAL_TYPE (MirbookingDefaultScoreTable, mirbooking_default_score_table, MIRBOOKING, DEFAULT_SCORE_TABLE, MirbookingScoreTable)

struct _MirbookingDefaultScoreTableClass
{
    MirbookingScoreTableClass parent_class;
};

typedef gboolean (*MirbookingDefaultScoreTableFilter) (MirbookingDefaultScoreTable *score_table,
                                                       MirbookingMirna             *mirna,
                                                       MirbookingTarget            *target,
                                                       gssize                       position,
                                                       gpointer                     user_data);


MirbookingDefaultScoreTable * mirbooking_default_score_table_new        (GBytes *seed_bytes,
                                                                         MirbookingDefaultScoreTableSupplementaryModel supplementary_model,
                                                                         GBytes *supp_bytes);
void                          mirbooking_default_score_table_set_filter (MirbookingDefaultScoreTable       *self,
                                                                         MirbookingDefaultScoreTableFilter  filter,
                                                                         gpointer                           user_data,
                                                                         GDestroyNotify                     destroy_func);

G_END_DECLS

#endif /* __MIRBOOKING_DEFAULT_SCORE_TABLE_H__ */
