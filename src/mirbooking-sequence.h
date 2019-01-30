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

const gchar * mirbooking_sequence_get_accession         (MirbookingSequence *self);
const gchar * mirbooking_sequence_get_name              (MirbookingSequence *self);
gchar       * mirbooking_sequence_get_sequence          (MirbookingSequence *self);
void          mirbooking_sequence_set_sequence          (MirbookingSequence *self,
                                                         const gchar *sequence);
const guint8 * mirbooking_sequence_get_raw_sequence     (MirbookingSequence *self, gsize *sequence_len);
void           mirbooking_sequence_set_raw_sequence     (MirbookingSequence *self,
                                                         const guint8       *sequence,
                                                         gssize              sequence_len);
gsize          mirbooking_sequence_get_sequence_length  (MirbookingSequence *self);
const guint8 * mirbooking_sequence_get_subsequence      (MirbookingSequence *self,
                                                         gsize               subsequence_offset,
                                                         gsize               subsequence_len);

gssize        mirbooking_sequence_get_subsequence_index (MirbookingSequence *self,
                                                         gsize               offset,
                                                         gsize               len);

guint         mirbooking_sequence_hash                  (const MirbookingSequence *a);
gboolean      mirbooking_sequence_equal                 (const MirbookingSequence *a,
                                                         const MirbookingSequence *b);

G_END_DECLS

#endif /* __MIRBOOKING_SEQUENCE_H__ */
