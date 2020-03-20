#include "mirbooking-utils.h"

#include <glib/gprintf.h>

G_DEFINE_POINTER_TYPE (MirbookingSequenceIter, mirbooking_sequence_iter)

typedef struct _MirbookingSequenceIter
{
    GMappedFile        *mapped_file;
    GDataInputStream   *dis;
    MirbookingFastaFormat fasta_format;
    GType               sequence_type;
    MirbookingSequence *cur_seq;
} MirbookingSequenceIter;

/**
 * mirbooking_read_sequences_from_mapped_file:
 * @iter: (out): A #MirbookingSequenceIter iterator capable of traversing the
 * sequence records in @mapped_file
 *
 * Read sequence records from a memory-mapped file. Each parsed sequence
 * records indirectly hold a reference to @mapped_file by referencing a slice
 * of the underlying memory via #GBytes.
 */
gboolean
mirbooking_read_sequences_from_mapped_file (GMappedFile             *mapped_file,
                                            MirbookingFastaFormat    fasta_format,
                                            GType                    sequence_type,
                                            MirbookingSequenceIter **iter)
{
    g_autoptr (MirbookingSequenceIter) ret = g_new0 (MirbookingSequenceIter, 1);

    g_return_val_if_fail (g_type_is_a (sequence_type, MIRBOOKING_TYPE_SEQUENCE), FALSE);

    const gchar *contents = g_mapped_file_get_contents (mapped_file);
    gsize contents_len = g_mapped_file_get_length (mapped_file);
    g_autoptr (GInputStream) in = g_memory_input_stream_new_from_data (contents, contents_len, NULL);
    g_autoptr (GDataInputStream) dis = g_data_input_stream_new (in);

    ret->mapped_file   = g_mapped_file_ref (mapped_file);
    ret->dis           = g_steal_pointer (&dis);
    ret->fasta_format  = fasta_format;
    ret->sequence_type = sequence_type;

    *iter = g_steal_pointer (&ret);

    return TRUE;
}

/**
 * mirbooking_sequence_iter_next:
 *
 * Move to the next sequence record in the iterator, parsing the record header
 * and obtaining the position of the stream where the sequence begins.
 *
 * Returns: %TRUE if a new record is available, %FALSE otherwise and @error
 * might be set.
 */
gboolean
mirbooking_sequence_iter_next (MirbookingSequenceIter  *iter,
                               GError                 **error)
{
    gchar *accession;
    gchar *name;
    gchar name_buffer[128]; /* in case the name does not appear literally */
    gchar *gene_accession = NULL;
    gchar *gene_name = NULL;
    gchar *seq;

    while (TRUE)
    {
        gsize line_len;
        g_autofree gchar *line = g_data_input_stream_read_line (iter->dis,
                                                                &line_len,
                                                                NULL,
                                                                error);

        if (line == NULL)
        {
            return FALSE;
        }

        g_return_val_if_fail (line_len > 0 && line[0] == '>', FALSE);

        if (iter->fasta_format == MIRBOOKING_FASTA_FORMAT_MIRBASE)
        {
            name      = strtok (line + 1, " ");
            accession = strtok (NULL, " ");
        }
        else if (iter->fasta_format == MIRBOOKING_FASTA_FORMAT_GENCODE)
        {
            // for GENCODE-style annotation, the name is in the sixth field
            accession = strtok (line + 1, "|");
            gene_accession = strtok (NULL, "|");

            gint i;
            for (i = 0; i < 3; i++)
            {
                name = strtok (NULL, "|");
            }

            gene_name = strtok (NULL, "|");
        }
        else if (iter->fasta_format == MIRBOOKING_FASTA_FORMAT_NCBI)
        {
            accession = strtok (line + 1, " ");

            gchar *name_p = NULL;
            if (strtok (NULL, "(,") && (name_p = strtok (NULL, "),")))
            {
                // if we parse the ")," following the gene name, name_p
                // will be NULL and we should look for the transcript
                // variant
                if (name_p)
                {
                    gene_name = name_p;
                }
            }

            /* unfortunately, that's the best information we get from the FASTA */
            gene_accession = gene_name;

            /* construct the gene name with its variant number */
            guint variant = 0;
            sscanf (strtok (NULL, ","), " transcript variant %u", &variant);
            g_sprintf (name_buffer, "%s-%03u", gene_name, 200 + variant);
            name = name_buffer;
        }
        else
        {
            accession = strtok (line + 1, " ");
            name = strtok (NULL, "\n");
        }

        goffset offset = g_seekable_tell (G_SEEKABLE (iter->dis));
        seq = g_mapped_file_get_contents (iter->mapped_file) + offset;

        gsize remaining = g_mapped_file_get_length (iter->mapped_file) - offset;
        gchar *next_seq = memchr (seq, '>', remaining);
        gsize seq_len = next_seq == NULL ? remaining - 1 : (gsize) (next_seq - seq) - 1;

        // skip to the next sequence
        if (!g_seekable_seek (G_SEEKABLE (iter->dis),
                              seq_len + 1,
                              G_SEEK_CUR,
                              NULL,
                              error))
        {
            return FALSE;
        }

        g_autoptr (MirbookingSequence) sequence = g_object_new (iter->sequence_type,
                                                                "accession",      accession,
                                                                "name",           name,
                                                                "gene-accession", gene_accession,
                                                                "gene-name",      gene_name,
                                                                NULL);

        mirbooking_sequence_set_raw_sequence (sequence,
                                              g_bytes_new_from_bytes (g_mapped_file_get_bytes (iter->mapped_file), offset, seq_len));

        if (iter->cur_seq)
        {
            g_object_unref (iter->cur_seq);
        }

        iter->cur_seq = g_steal_pointer (&sequence);

        return TRUE;
    }

    return FALSE;
}

/**
 * mirbooking_sequence_iter_get_sequence:
 * Returns: (transfer none)
 */
MirbookingSequence *
mirbooking_sequence_iter_get_sequence (MirbookingSequenceIter *iter)
{
    return iter->cur_seq;
}

void
mirbooking_sequence_iter_free (MirbookingSequenceIter *iter)
{
    g_mapped_file_unref (iter->mapped_file);
    g_object_unref (iter->dis);
    if (iter->cur_seq)
    {
        g_object_unref (iter->cur_seq);
    }
    g_free (iter);
}
