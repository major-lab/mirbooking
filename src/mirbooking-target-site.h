#ifndef __MIRBOOKING_TARGET_SITE_H__
#define __MIRBOOKING_TARGET_SITE_H__

#include "mirbooking-mirna.h"
#include "mirbooking-target.h"

#include <glib-object.h>

G_BEGIN_DECLS

#define MIRBOOKING_TYPE_TARGET_SITE mirbooking_target_site_get_type ()
G_DECLARE_FINAL_TYPE (MirbookingTargetSite, mirbooking_target_site, MIRBOOKING, TARGET_SITE, GObject)

struct _MirbookingTargetSiteClass
{
    GObjectClass parent_class;
};

MirbookingTargetSite * mirbooking_target_site_new                (MirbookingTarget *target,
                                                                  gsize             offset);
MirbookingTarget *     mirbooking_target_site_get_target         (MirbookingTargetSite *target_site);
gsize                  mirbooking_target_site_get_offset         (MirbookingTargetSite *target_site);
gfloat                 mirbooking_target_site_get_mirna_quantity (MirbookingTargetSite *target_site,
                                                                  MirbookingMirna      *mirna);
void                   mirbooking_target_site_set_mirna_quantity (MirbookingTargetSite *target_site,
                                                                  MirbookingMirna      *mirna,
                                                                  gfloat                quantity);
GSList *               mirbooking_target_site_get_mirnas         (MirbookingTargetSite *target_site);

G_END_DECLS

#endif /* __MIRBOOKING_TARGET_SITE_H__ */
