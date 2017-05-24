#ifndef __MIRBOOKING_H__
#define __MIRBOOKING_H__

#include "mirbooking-score-table.h"
#include "mirbooking-sequence.h"
#include "mirbooking-target.h"
#include "mirbooking-mirna.h"

#include <glib-object.h>

G_BEGIN_DECLS

#define MIRBOOKING_DEFAULT_THRESHOLD 0.00125f
#define MIRBOOKING_DEFAULT_LOG_BASE  512.0f

#define MIRBOOKING_TYPE mirbooking_get_type ()
G_DECLARE_FINAL_TYPE (Mirbooking, mirbooking, MIRBOOKING, MIRBOOKING, GObject)

struct _MirbookingClass
{
    GObjectClass parent_class;
};

typedef struct _MirbookingTargetSite
{
    MirbookingTarget *target;
    gsize             site_offset;
    GSList           *mirna_quantities; // #GSList of #MirbookingMirnaQuantity
} MirbookingTargetSite;

typedef struct _MirbookingMirnaQuantity
{
    MirbookingMirna *mirna;
    gfloat           quantity;
} MirbookingMirnaQuantity;

Mirbooking * mirbooking_new                 (void);
void         mirbooking_set_threshold       (Mirbooking *self,
                                             gfloat      threshold);
void         mirbooking_set_log_base        (Mirbooking *self,
                                             gfloat      log_base);
void         mirbooking_set_cds_multiplier  (Mirbooking *self,
                                             gfloat      cds_multiplier);
void         mirbooking_set_score_table     (Mirbooking           *self,
                                             MirbookingScoreTable *score_table);
void         mirbooking_set_target_quantity (Mirbooking       *self,
                                             MirbookingTarget *target,
                                             gfloat            quantity);
void         mirbooking_set_mirna_quantity  (Mirbooking      *self,
                                             MirbookingMirna *mirna,
                                             gfloat           quantity);

gboolean     mirbooking_run                 (Mirbooking *self, GError **error);

MirbookingTargetSite * mirbooking_get_target_sites (Mirbooking *self,
                                                    gsize      *target_sites_len);

G_END_DECLS

#endif /* __MIRBOOKING_H__ */
