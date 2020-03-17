#ifndef __MIRBOOKING_BROKER_H__
#define __MIRBOOKING_BROKER_H__

#include "mirbooking-broker-sparse-solver.h"
#include "mirbooking-broker-output-format.h"
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

/*
 *
 */
#define MIRBOOKING_BROKER_DEFAULT_SPARSE_SOLVER mirbooking_broker_sparse_solver_get_default ()

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
gdouble                mirbooking_broker_get_sequence_quantity   (MirbookingBroker   *self,
                                                                  MirbookingSequence *sequence);
void                   mirbooking_broker_set_sequence_quantity   (MirbookingBroker   *self,
                                                                  MirbookingSequence *sequence,
                                                                  gdouble             quantity);
gdouble                mirbooking_broker_get_time                (MirbookingBroker *self);
void                   mirbooking_broker_set_time                (MirbookingBroker *self,
                                                                  gdouble           time);
gboolean               mirbooking_broker_evaluate                (MirbookingBroker          *self,
                                                                  gdouble                   *error_ratio,
                                                                  GError                   **error);
gboolean               mirbooking_broker_step                    (MirbookingBroker          *self,
                                                                  MirbookingBrokerStepMode   step_mode,
                                                                  gdouble                    step_size,
                                                                  GError                   **error);
const GArray *         mirbooking_broker_get_target_sites        (MirbookingBroker *self);
gdouble                mirbooking_broker_get_target_site_quantity (MirbookingBroker *self,
                                                                   const MirbookingTargetSite *target_site);
gdouble                mirbooking_broker_get_target_site_occupants_quantity (MirbookingBroker *self,
                                                                             const MirbookingTargetSite *target_site);
gdouble                mirbooking_broker_get_target_site_kother  (MirbookingBroker *self,
                                                                  const MirbookingTargetSite *target_site);

const GPtrArray *      mirbooking_broker_get_mirnas              (MirbookingBroker *self);
const GPtrArray *      mirbooking_broker_get_targets             (MirbookingBroker *self);
const GArray *         mirbooking_broker_get_occupants           (MirbookingBroker *self);
gdouble                mirbooking_broker_get_occupant_quantity   (MirbookingBroker         *self,
                                                                  const MirbookingOccupant *occupant);
void                   mirbooking_broker_set_occupant_quantity   (MirbookingBroker         *self,
                                                                  const MirbookingOccupant *occupant,
                                                                  gdouble                   quantity);

gdouble                mirbooking_broker_get_bound_mirna_quantity (MirbookingBroker *self,
                                                                   MirbookingMirna *mirna);

gdouble *              mirbooking_broker_get_target_occupants_pmf (MirbookingBroker *self,
                                                                   MirbookingTarget  *target,
                                                                   gsize             *pmf_len);
gdouble                mirbooking_broker_get_target_expressed_fraction (MirbookingBroker *self,
                                                                        MirbookingTarget *target);
gdouble                mirbooking_broker_get_product_quantity          (MirbookingBroker *self,
                                                                        MirbookingTarget *target);
gdouble                mirbooking_broker_get_mirna_transcription_rate  (MirbookingBroker *self,
                                                                        MirbookingMirna  *mirna);
void                   mirbooking_broker_set_mirna_transcription_rate  (MirbookingBroker *self,
                                                                        MirbookingMirna  *mirna,
                                                                        gdouble           transcription_rate);
gdouble                mirbooking_broker_get_mirna_degradation_rate    (MirbookingBroker *self,
                                                                        MirbookingMirna  *mirna);
gdouble                mirbooking_broker_get_target_transcription_rate (MirbookingBroker *self,
                                                                        MirbookingTarget *target);
void                   mirbooking_broker_set_target_transcription_rate (MirbookingBroker *self,
                                                                        MirbookingTarget *target,
                                                                        gdouble           transcription_rate);
gdouble                mirbooking_broker_get_target_degradation_rate   (MirbookingBroker *self,
                                                                        MirbookingTarget *target);
gdouble                mirbooking_broker_get_product_degradation_rate  (MirbookingBroker *self,
                                                                        MirbookingTarget *target);

gboolean               mirbooking_broker_write_output_to_stream (MirbookingBroker              *broker,
                                                                 GOutputStream                 *out,
                                                                 MirbookingBrokerOutputFormat   output_format,
                                                                 GError                       **error);
gboolean               mirbooking_broker_write_output_to_file   (MirbookingBroker              *self,
                                                                 GFile                         *output_file,
                                                                 MirbookingBrokerOutputFormat   output_format,
                                                                 GError                        **error);

G_END_DECLS

#endif /* __MIRBOOKING_BROKER_H__ */
