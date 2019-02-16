#ifndef __MIRBOOKING_BROKER_H__
#define __MIRBOOKING_BROKER_H__

#include "mirbooking-broker-sparse-solver.h"
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

#define MIRBOOKING_BROKER_DEFAULT_SPARSE_SOLVER    MIRBOOKING_BROKER_SPARSE_SOLVER_SUPERLU

/*
 * The full effect seems to be occurring with at least 17 nucleotides distance
 * between consecutive seed starts. We attribute 7 nucleotides toward the 3'
 * segment to account for the seed bindings and attribute the rest to the tail.
 *
 * Note that the binding position (i.e. start of the seed) is always accounted
 * for.
 *
 * Reference: Pål Sætrom et al., “Distance Constraints between MicroRNA Target
 * Sites Dictate Efficacy and Cooperativity,” Nucleic Acids Research 35, no. 7
 * (April 2007): 2333–42, https://doi.org/10.1093/nar/gkm133.
 */
#define MIRBOOKING_BROKER_DEFAULT_5PRIME_FOOTPRINT 9
#define MIRBOOKING_BROKER_DEFAULT_3PRIME_FOOTPRINT 7

#define MIRBOOKING_BROKER_TYPE mirbooking_broker_get_type ()
G_DECLARE_FINAL_TYPE (MirbookingBroker, mirbooking_broker, MIRBOOKING, BROKER, GObject)

struct _MirbookingBrokerClass
{
    GObjectClass parent_class;
};

/**
 * MirbookingBrokerStepMode:
 * @MIRBOOKING_BROKER_STEP_MODE_INTEGRATE:
 * @MIRBOOKING_BROKER_STEP_MODE_SOLVE_STEADY_STATE:
 */
typedef enum _MirbookingBrokerStepMode
{
    MIRBOOKING_BROKER_STEP_MODE_INTEGRATE,
    MIRBOOKING_BROKER_STEP_MODE_SOLVE_STEADY_STATE
} MirbookingBrokerStepMode;

MirbookingBroker *     mirbooking_broker_new                     (void);
MirbookingBroker *     mirbooking_broker_new_with_rank           (gint rank);
gint                   mirbooking_broker_get_rank                (MirbookingBroker *self);
void                   mirbooking_broker_set_sparse_solver       (MirbookingBroker             *self,
                                                                  MirbookingBrokerSparseSolver  sparse_solver);
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
gdouble                mirbooking_broker_get_time                (MirbookingBroker *self);
void                   mirbooking_broker_set_time                (MirbookingBroker *self,
                                                                  gdouble           time);
gboolean               mirbooking_broker_evaluate                (MirbookingBroker          *self,
                                                                  gdouble                   *norm,
                                                                  GError                   **error);
gboolean               mirbooking_broker_step                    (MirbookingBroker          *self,
                                                                  MirbookingBrokerStepMode   step_mode,
                                                                  gdouble                    step_size,
                                                                  GError                   **error);
const GArray *         mirbooking_broker_get_target_sites        (MirbookingBroker *self);
const GArray *         mirbooking_broker_get_occupants           (MirbookingBroker *self);
gdouble                mirbooking_broker_get_occupant_quantity   (MirbookingBroker         *self,
                                                                  const MirbookingOccupant *occupant);
void                   mirbooking_broker_set_occupant_quantity   (MirbookingBroker         *self,
                                                                  const MirbookingOccupant *occupant,
                                                                  gdouble                   quantity);

gdouble *              mirbooking_broker_get_target_occupants_pmf (MirbookingBroker *self,
                                                                   MirbookingTarget  *target,
                                                                   gsize             *pmf_len);

G_END_DECLS

#endif /* __MIRBOOKING_BROKER_H__ */
