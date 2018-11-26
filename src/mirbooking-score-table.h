#ifndef __MIRBOOKING_SCORE_TABLE_H__
#define __MIRBOOKING_SCORE_TABLE_H__

#include "mirbooking-error.h"
#include "mirbooking-target.h"
#include "mirbooking-mirna.h"

#include <glib-object.h>

G_BEGIN_DECLS

/**
 * Forward rate constant for the [E] + [S] -> [ES] reaction, which is
 * indistinguishable from seed types or supplementary bindings (Wee et al. 2012).
 *
 * Reference:
 * Liang Meng Wee et al., “Argonaute Divides Its RNA Guide into Domains
 * with Distinct Functions and RNA-Binding Properties,” Cell 151, no. 5
 * (November 21, 2012): 1055–67,
 * https://doi.org/10.1016/j.cell.2012.10.036.
 */
#define MIRBOOKING_SCORE_TABLE_DEFAULT_LAMBDA 3.6e-5 // pM^-1s^-1

#define MIRBOOKING_TYPE_SCORE_TABLE mirbooking_score_table_get_type ()
G_DECLARE_DERIVABLE_TYPE (MirbookingScoreTable, mirbooking_score_table, MIRBOOKING, SCORE_TABLE, GObject)

struct _MirbookingScoreTableClass
{
    GObjectClass parent_class;

    gboolean (*compute_positions) (MirbookingScoreTable  *self,
                                   MirbookingMirna       *mirna,
                                   MirbookingTarget      *target,
                                   gsize                **positions,
                                   gsize                 *positions_len,
                                   GError               **error);

    gdouble (*compute_score)    (MirbookingScoreTable *self,
                                 MirbookingMirna      *mirna,
                                 MirbookingTarget     *target,
                                 gsize                 position,
                                 GError              **error);


    gdouble (*compute_enzymatic_score) (MirbookingScoreTable *self,
                                        MirbookingMirna      *mirna,
                                        MirbookingTarget     *target,
                                        gsize                 position,
                                        GError              **error);

};

gboolean mirbooking_score_table_compute_positions  (MirbookingScoreTable  *self,
                                                    MirbookingMirna       *mirna,
                                                    MirbookingTarget      *target,
                                                    gsize                **positions,
                                                    gsize                 *positions_len,
                                                    GError               **error);

gdouble   mirbooking_score_table_compute_score  (MirbookingScoreTable *self,
                                                 MirbookingMirna      *mirna,
                                                 MirbookingTarget     *target,
                                                 gsize                 position,
                                                 GError              **error);

gdouble mirbooking_score_table_compute_enzymatic_score (MirbookingScoreTable *self,
                                                         MirbookingMirna      *mirna,
                                                         MirbookingTarget     *target,
                                                         gsize                 position,
                                                         GError              **error);

G_END_DECLS

#endif /* __MIRBOOKING_SCORE_TABLE_H__ */
