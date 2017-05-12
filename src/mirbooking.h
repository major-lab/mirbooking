#ifndef __MIRBOOKING_H__
#define __MIRBOOKING_H__

#include "mirbooking-score-table.h"
#include "mirbooking-sequence.h"
#include "mirbooking-target.h"
#include "mirbooking-mirna.h"

#include <glib-object.h>

G_BEGIN_DECLS

#define MIRBOOKING_DEFAULT_THRESHOLD 0.00125
#define MIRBOOKING_DEFAULT_LOG_BASE  2.0

#define MIRBOOKING_TYPE mirbooking_get_type ()
G_DECLARE_FINAL_TYPE (Mirbooking, mirbooking, MIRBOOKING, MIRBOOKING, GObject)

struct _MirbookingClass
{
    GObjectClass parent_class;
};

typedef struct _MirbookingTargetSite MirbookingTargetSite;

Mirbooking * mirbooking_new                 (void);
void         mirbooking_set_threshold       (Mirbooking *self,
                                             gdouble     threshold);
void         mirbooking_set_log_base        (Mirbooking *self,
                                             gdouble     log_base);
void         mirbooking_set_score_table     (Mirbooking           *self,
                                             MirbookingScoreTable *score_table);
void         mirbooking_set_target_quantity (Mirbooking       *self,
                                             MirbookingTarget *target,
                                             gfloat            quantity);
void         mirbooking_set_mirna_quantity  (Mirbooking      *self,
                                             MirbookingMirna *mirna,
                                             gfloat           quantity);

GSList *     mirbooking_get_targets         (Mirbooking *self);
GSList *     mirbooking_get_mirnas          (Mirbooking *self);

gboolean     mirbooking_run                 (Mirbooking *self, GError **error);

GSList *     mirbooking_get_target_sites    (Mirbooking       *self,
                                             MirbookingTarget *target);

gsize        mirbooking_target_site_get_site_offset    (MirbookingTargetSite *target_site);
GSList *     mirbooking_target_site_get_mirnas         (MirbookingTargetSite *target_site);
gfloat       mirbooking_target_site_get_mirna_quantity (MirbookingTargetSite *target_site,
                                                        MirbookingMirna      *mirna);

G_END_DECLS

#endif /* __MIRBOOKING_H__ */
