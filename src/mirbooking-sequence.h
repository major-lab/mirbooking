#ifndef __MIRBOOKING_SEQUENCE_H__
#define __MIRBOOKING_SEQUENCE_H__

#include <glib-object.h>

G_BEGIN_DECLS

#define MIRBOOKING_TYPE_SEQUENCE mirbooking_sequence_get_type ()
G_DECLARE_DERIVABLE_TYPE (MirbookingSequence, mirbooking_sequence, MIRBOOKING, SEQUENCE, GObject)

struct _MirbookingSequenceClass
{
    GObjectClass parent_class;
};

const gchar * mirbooking_sequence_get_accession (MirbookingSequence *self);
const gchar * mirbooking_sequence_get_sequence  (MirbookingSequence *self, gsize *sequence_len);
void          mirbooking_sequence_set_sequence  (MirbookingSequence *self,
                                                 const gchar        *sequence,
                                                 gsize               sequence_len);
gsize         mirbooking_sequence_get_sequence_length (MirbookingSequence *self);
const gchar * mirbooking_sequence_get_subsequence (MirbookingSequence *self,
                                                   gsize               subsequence_offset,
                                                   gsize               subsequence_len);

guint    mirbooking_sequence_hash  (const MirbookingSequence *a);
gboolean mirbooking_sequence_equal (const MirbookingSequence *a,
                                    const MirbookingSequence *b);

G_END_DECLS

#endif /* __MIRBOOKING_SEQUENCE_H__ */
