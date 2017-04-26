#ifndef __MIRBOOKING_TARGET_H__
#define __MIRBOOKING_TARGET_H__

#include "mirbooking-sequence.h"

#include <glib-object.h>

G_BEGIN_DECLS

#define MIRBOOKING_TYPE_TARGET mirbooking_target_get_type ()
G_DECLARE_FINAL_TYPE (MirbookingTarget, mirbooking_target, MIRBOOKING, TARGET, MirbookingSequence)

struct _MirbookingTargetClass
{
    MirbookingSequenceClass parent_class;
};

MirbookingTarget * mirbooking_target_new (const gchar *accession);

G_END_DECLS

#endif /* __MIRBOOKING_TARGET_H__ */
