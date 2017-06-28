#ifndef __MIRBOOKING_H__
#define __MIRBOOKING_H__

#include "mirbooking-score-table.h"
#include "mirbooking-sequence.h"
#include "mirbooking-target.h"
#include "mirbooking-mirna.h"

#include <glib-object.h>

G_BEGIN_DECLS

#define MIRBOOKING_DEFAULT_THRESHOLD        0.00125f
#define MIRBOOKING_DEFAULT_LOG_BASE         512.0f
#define MIRBOOKING_DEFAULT_5PRIME_FOOTPRINT 26
#define MIRBOOKING_DEFAULT_3PRIME_FOOTPRINT 19

#define MIRBOOKING_TYPE mirbooking_get_type ()
G_DECLARE_FINAL_TYPE (Mirbooking, mirbooking, MIRBOOKING, MIRBOOKING, GObject)

struct _MirbookingClass
{
    GObjectClass parent_class;
};

typedef struct _MirbookingTargetSite
{
    MirbookingTarget *target;
    gsize             position;
    GSList           *occupants; // #GSList of #MirbookingOccupant
    guint             occupancy; // sum of all occupants' quantities
} MirbookingTargetSite;

typedef struct _MirbookingOccupant
{
    MirbookingMirna *mirna;
    guint            quantity;
} MirbookingOccupant;

Mirbooking * mirbooking_new                   (void);
void         mirbooking_set_threshold         (Mirbooking *self,
                                               gfloat      threshold);
void         mirbooking_set_log_base          (Mirbooking *self,
                                               gfloat      log_base);
void         mirbooking_set_5prime_footprint  (Mirbooking *self,
                                               gsize       footprint);
void         mirbooking_set_3prime_footprint  (Mirbooking *self,
                                               gsize       footprint);
void         mirbooking_set_score_table       (Mirbooking           *self,
                                               MirbookingScoreTable *score_table);
void         mirbooking_set_sequence_quantity (Mirbooking         *self,
                                               MirbookingSequence *sequence,
                                               gfloat              quantity);
gboolean     mirbooking_run                   (Mirbooking *self, GError **error);

/* results */
const MirbookingTargetSite * mirbooking_get_target_sites (Mirbooking *self,
                                                          gsize      *target_sites_len);

G_END_DECLS

#endif /* __MIRBOOKING_H__ */
