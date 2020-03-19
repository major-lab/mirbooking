#ifndef __MIRBOOKING_UTILS_H__
#define __MIRBOOKING_UTILS_H__

#include "mirbooking-sequence.h"

#include <glib.h>
#include <gio/gio.h>

G_BEGIN_DECLS

typedef enum _MirbookingFastaFormat
{
    MIRBOOKING_FASTA_FORMAT_GENERIC,
    MIRBOOKING_FASTA_FORMAT_NCBI,
    MIRBOOKING_FASTA_FORMAT_GENCODE,
    MIRBOOKING_FASTA_FORMAT_MIRBASE
} MirbookingFastaFormat;

typedef struct _MirbookingSequenceIter MirbookingSequenceIter;

gboolean             mirbooking_sequence_iter_next         (MirbookingSequenceIter *iter, GError **error);
MirbookingSequence * mirbooking_sequence_iter_get_sequence (MirbookingSequenceIter *iter);
void                 mirbooking_sequence_iter_free         (MirbookingSequenceIter *iter);

gboolean mirbooking_read_sequences_from_mapped_file (GMappedFile             *mapped_file,
                                                     MirbookingFastaFormat    fasta_format,
                                                     GType                    sequence_type,
                                                     MirbookingSequenceIter **iter);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(MirbookingSequenceIter, mirbooking_sequence_iter_free)

G_END_DECLS

#endif /* __MIRBOOKING_UTILS_H__ */
