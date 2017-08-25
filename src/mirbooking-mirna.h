#ifndef __MIRBOOKING_MIRNA_H__
#define __MIRBOOKING_MIRNA_H__

#include "mirbooking-sequence.h"

#include <glib-object.h>

G_BEGIN_DECLS

#define MIRBOOKING_TYPE_MIRNA mirbooking_mirna_get_type ()
G_DECLARE_FINAL_TYPE (MirbookingMirna, mirbooking_mirna, MIRBOOKING, MIRNA, MirbookingSequence)

struct _MirbookingMirnaClass
{
    MirbookingSequenceClass parent_class;
};

MirbookingMirna * mirbooking_mirna_new           (const gchar *accession);
MirbookingMirna * mirbooking_mirna_new_with_name (const gchar *accession,
                                                  const gchar *name);

G_END_DECLS

#endif /* __MIRBOOKING_MIRNA_H__ */
