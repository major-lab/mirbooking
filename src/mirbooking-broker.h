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

#define MIRBOOKING_BROKER_DEFAULT_KAPPA            6.4135955621023e-5
#define MIRBOOKING_BROKER_DEFAULT_CUTOFF           500
#define MIRBOOKING_BROKER_DEFAULT_STEP_SIZE        1e-6 // 7e-6
#define MIRBOOKING_BROKER_DEFAULT_TOLERANCE        1e-3
#define MIRBOOKING_BROKER_DEFAULT_MAX_ITERATIONS   1e8
#define MIRBOOKING_BROKER_DEFAULT_5PRIME_FOOTPRINT 26
#define MIRBOOKING_BROKER_DEFAULT_3PRIME_FOOTPRINT 19

#define MIRBOOKING_BROKER_TYPE mirbooking_broker_get_type ()
G_DECLARE_FINAL_TYPE (MirbookingBroker, mirbooking_broker, MIRBOOKING, BROKER, GObject)

struct _MirbookingBrokerClass
{
    GObjectClass parent_class;
};

MirbookingBroker *     mirbooking_broker_new                     (void);
void                   mirbooking_broker_set_kappa               (MirbookingBroker *self,
                                                                  gdouble           kappa);
void                   mirbooking_broker_set_cutoff              (MirbookingBroker *self,
                                                                  gfloat            cutoff);
void                   mirbooking_broker_set_step_size           (MirbookingBroker *self,
                                                                  gdouble           step_size);
void                   mirbooking_broker_set_tolerance           (MirbookingBroker *self,
                                                                  gdouble           tolerance);
void                   mirbooking_broker_set_max_iterations      (MirbookingBroker *self,
                                                                  guint64           max_iterations);
void                   mirbooking_broker_set_5prime_footprint    (MirbookingBroker *self,
                                                                  gsize             footprint);
void                   mirbooking_broker_set_3prime_footprint    (MirbookingBroker *self,
                                                                  gsize             footprint);
MirbookingScoreTable * mirbooking_broker_get_score_table         (MirbookingBroker *self);
void                   mirbooking_broker_set_score_table         (MirbookingBroker     *self,
                                                                  MirbookingScoreTable *score_table);
gfloat                 mirbooking_broker_get_sequence_quantity   (MirbookingBroker   *self,
                                                                  MirbookingSequence *sequence);
void                   mirbooking_broker_set_sequence_quantity   (MirbookingBroker   *self,
                                                                  MirbookingSequence *sequence,
                                                                  gfloat              quantity);
gboolean               mirbooking_broker_iter                    (MirbookingBroker  *self,
                                                                  gdouble            step_size,
                                                                  GError           **error);
gboolean               mirbooking_broker_run                     (MirbookingBroker  *self,
                                                                  GError           **error);
void                   mirbooking_broker_run_async               (MirbookingBroker    *self,
                                                                  GAsyncReadyCallback  callback,
                                                                  gpointer             callback_data);
gboolean               mirbooking_broker_run_finish              (MirbookingBroker  *self,
                                                                  GAsyncResult      *result,
                                                                  GError           **error);
GArray *               mirbooking_broker_get_target_sites        (MirbookingBroker *self);
gdouble                mirbooking_broker_get_target_site_vacancy (MirbookingBroker           *self,
                                                                  const MirbookingTargetSite *target_site);
gdouble                mirbooking_broker_get_target_silencing    (MirbookingBroker *self,
                                                                  MirbookingTarget *target);

G_END_DECLS

#endif /* __MIRBOOKING_BROKER_H__ */
