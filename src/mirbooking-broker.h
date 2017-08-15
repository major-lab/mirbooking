#ifndef __MIRBOOKING_BROKER_H__
#define __MIRBOOKING_BROKER_H__

#include "mirbooking-error.h"
#include "mirbooking-mirna.h"
#include "mirbooking-occupant.h"
#include "mirbooking-score-table.h"
#include "mirbooking-sequence.h"
#include "mirbooking-target-site.h"
#include "mirbooking-target.h"

#include <glib-object.h>
#include <gio/gio.h>

G_BEGIN_DECLS

#define MIRBOOKING_BROKER_DEFAULT_THRESHOLD        0.0179f
#define MIRBOOKING_BROKER_DEFAULT_LOG_BASE         512.0f
#define MIRBOOKING_BROKER_DEFAULT_5PRIME_FOOTPRINT 26
#define MIRBOOKING_BROKER_DEFAULT_3PRIME_FOOTPRINT 19

#define MIRBOOKING_BROKER_TYPE mirbooking_broker_get_type ()
G_DECLARE_FINAL_TYPE (MirbookingBroker, mirbooking_broker, MIRBOOKING, BROKER, GObject)

struct _MirbookingBrokerClass
{
    GObjectClass parent_class;
};

MirbookingBroker *     mirbooking_broker_new                     (void);
void                   mirbooking_broker_set_threshold           (MirbookingBroker *self,
                                                                  gfloat            threshold);
void                   mirbooking_broker_set_log_base            (MirbookingBroker *self,
                                                                  gfloat            log_base);
void                   mirbooking_broker_set_5prime_footprint    (MirbookingBroker *self,
                                                                  gsize             footprint);
void                   mirbooking_broker_set_3prime_footprint    (MirbookingBroker *self,
                                                                  gsize             footprint);
MirbookingScoreTable * mirbooking_broker_get_score_table         (MirbookingBroker *self);
void                   mirbooking_broker_set_score_table         (MirbookingBroker     *self,
                                                                  MirbookingScoreTable *score_table);
void                   mirbooking_broker_set_sequence_quantity   (MirbookingBroker   *self,
                                                                  MirbookingSequence *sequence,
                                                                  gfloat              quantity);
gboolean               mirbooking_broker_run                     (MirbookingBroker  *self,
                                                                  GError           **error);
void                   mirbooking_broker_run_async               (MirbookingBroker    *self,
                                                                  GAsyncReadyCallback  callback,
                                                                  gpointer             callback_data);
gboolean               mirbooking_broker_run_finish              (MirbookingBroker  *self,
                                                                  GAsyncResult      *result,
                                                                  GError           **error);
GArray *               mirbooking_broker_get_target_sites        (MirbookingBroker *self);
guint                  mirbooking_broker_get_target_site_vacancy (MirbookingBroker           *self,
                                                                  const MirbookingTargetSite *target_site);

G_END_DECLS

#endif /* __MIRBOOKING_BROKER_H__ */
