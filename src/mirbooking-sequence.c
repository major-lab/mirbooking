#include "mirbooking-sequence.h"

#include <ctype.h>
#include <string.h>
#include <glib-object.h>

typedef struct
{
    gchar       *accession;
    gchar       *name;
    gchar       *gene_accession;
    gchar       *gene_name;
    GBytes      *sequence;
    GPtrArray   *sequence_skips; /* all the pointers in 'sequence' to skip (i.e. line feeds) */
} MirbookingSequencePrivate;

enum
{
    PROP_ACCESSION = 1,
    PROP_NAME,
    PROP_GENE_ACCESSION,
    PROP_GENE_NAME,
    PROP_SEQUENCE
};

G_DEFINE_TYPE_WITH_PRIVATE (MirbookingSequence, mirbooking_sequence, G_TYPE_OBJECT)

static void
mirbooking_sequence_init (MirbookingSequence *self)
{
}

static void
mirbooking_sequence_set_property (GObject *object, guint property_id, const GValue *value, GParamSpec *pspec)
{
    MirbookingSequencePrivate *priv = mirbooking_sequence_get_instance_private (MIRBOOKING_SEQUENCE (object));

    const gchar *seq;

    switch (property_id)
    {
        case PROP_ACCESSION:
            priv->accession = g_value_dup_string ((value));
            break;
        case PROP_NAME:
            priv->name = g_value_dup_string (value);
            break;
        case PROP_GENE_ACCESSION:
            priv->gene_accession = g_value_dup_string (value);
            break;
        case PROP_GENE_NAME:
            priv->gene_name = g_value_dup_string (value);
            break;
        case PROP_SEQUENCE:
            seq = g_value_get_string (value);
            mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (object),
                                                  g_bytes_new (seq, strlen (seq)));
            break;
        default:
            G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
    }
}

static void
mirbooking_sequence_get_property (GObject *object, guint property_id, GValue *value, GParamSpec *pspec)
{
    MirbookingSequencePrivate *priv = mirbooking_sequence_get_instance_private (MIRBOOKING_SEQUENCE (object));

    switch (property_id)
    {
        case PROP_ACCESSION:
            g_value_set_string (value, priv->accession);
            break;
        case PROP_NAME:
            g_value_set_string (value, priv->name);
            break;
        case PROP_GENE_ACCESSION:
            g_value_set_string (value, priv->gene_accession);
            break;
        case PROP_GENE_NAME:
            g_value_set_string (value, priv->gene_name);
            break;
        case PROP_SEQUENCE:
            g_value_take_string (value, mirbooking_sequence_get_sequence (MIRBOOKING_SEQUENCE (object)));
            break;
        default:
            G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
    }
}

static void
mirbooking_sequence_finalize (GObject *object)
{
    MirbookingSequencePrivate *priv = mirbooking_sequence_get_instance_private (MIRBOOKING_SEQUENCE (object));

    g_free (priv->accession);
    g_ptr_array_unref (priv->sequence_skips);
    if (priv->sequence)
        g_bytes_unref (priv->sequence);

    G_OBJECT_CLASS (mirbooking_sequence_parent_class)->finalize (object);
}

static void
mirbooking_sequence_class_init (MirbookingSequenceClass *klass)
{
    GObjectClass *object_class = G_OBJECT_CLASS (klass);

    object_class->set_property = mirbooking_sequence_set_property;
    object_class->get_property = mirbooking_sequence_get_property;
    object_class->finalize     = mirbooking_sequence_finalize;

    g_object_class_install_property (object_class,
                                     PROP_ACCESSION,
                                     g_param_spec_string ("accession", "Accession", "", NULL, G_PARAM_CONSTRUCT_ONLY | G_PARAM_READWRITE));
    g_object_class_install_property (object_class,
                                     PROP_NAME,
                                     g_param_spec_string ("name", "name", "", NULL, G_PARAM_CONSTRUCT_ONLY | G_PARAM_READWRITE));
    g_object_class_install_property (object_class,
                                     PROP_SEQUENCE,
                                     g_param_spec_string ("sequence", "Sequence", "", NULL, G_PARAM_READWRITE));
}

const gchar *
mirbooking_sequence_get_accession (MirbookingSequence *self)
{
    MirbookingSequencePrivate *priv = mirbooking_sequence_get_instance_private (self);

    return priv->accession;
}

const gchar *
mirbooking_sequence_get_name (MirbookingSequence *self)
{
    MirbookingSequencePrivate *priv = mirbooking_sequence_get_instance_private (self);

    return priv->name;
}


const gchar *
mirbooking_sequence_get_gene_accession (MirbookingSequence *self)
{
    MirbookingSequencePrivate *priv = mirbooking_sequence_get_instance_private (self);

    return priv->gene_accession;
}

const gchar *
mirbooking_sequence_get_gene_name (MirbookingSequence *self)
{
    MirbookingSequencePrivate *priv = mirbooking_sequence_get_instance_private (self);

    return priv->gene_name;
}

void
mirbooking_sequence_set_gene_accession (MirbookingSequence *self, const gchar *gene_accession)
{
    MirbookingSequencePrivate *priv = mirbooking_sequence_get_instance_private (self);

    priv->gene_accession = g_strdup (gene_accession);
}

void
mirbooking_sequence_set_gene_name (MirbookingSequence *self, const gchar *gene_name)
{
    MirbookingSequencePrivate *priv = mirbooking_sequence_get_instance_private (self);

    priv->gene_name = g_strdup (gene_name);
}

/**
 * mirbooking_sequence_get_raw_sequence:
 * @self: A #MirbookingSequence
 *
 * Obtain the internal representation of the sequence, which corresponds to
 * what has been set previsouly via #mirbooking_sequence_set_raw_sequence.
 *
 * Returns: A view over the internal sequence.
 */
GBytes *
mirbooking_sequence_get_raw_sequence (MirbookingSequence *self)
{
    MirbookingSequencePrivate *priv = mirbooking_sequence_get_instance_private (self);

    return g_bytes_ref (priv->sequence);
}

/**
 * mirbooking_sequence_set_raw_sequence:
 * @self: A #MirbookingSequence
 * @sequence: A sequence of nucleotides
 *
 * Note that the passed sequence is not copied internally and thus, it must
 * live for as long as this object does.
 */
void
mirbooking_sequence_set_raw_sequence (MirbookingSequence *self, GBytes *sequence)
{
    MirbookingSequencePrivate *priv = mirbooking_sequence_get_instance_private (self);

    if (sequence != priv->sequence)
    {
        if (priv->sequence)
            g_bytes_unref (priv->sequence);

        priv->sequence = g_bytes_ref (sequence);

        if (priv->sequence_skips)
        {
            g_ptr_array_unref (priv->sequence_skips);
        }

        gsize sequence_len;
        const guint8* seq = g_bytes_get_data (sequence, &sequence_len);

        priv->sequence_skips = g_ptr_array_sized_new (sequence_len / 80); // line feed every 80 characters

        const guint8* ptr = seq;
        while ((ptr = memchr (ptr, '\n', sequence_len - (ptr - seq))))
        {
            g_ptr_array_add (priv->sequence_skips, (gpointer) ptr);
            ptr++; // jump right after the line feed
        }

        g_object_notify (G_OBJECT (self), "sequence");
    }
}

gchar *
mirbooking_sequence_get_sequence (MirbookingSequence *self)
{
    MirbookingSequencePrivate *priv = mirbooking_sequence_get_instance_private (self);

    gsize i, j = 0;
    gchar *ret;

    gsize seq_len;
    const gchar *seq = g_bytes_get_data (priv->sequence, &seq_len);

    ret = g_new0 (gchar, mirbooking_sequence_get_sequence_length (self) + 1);

    for (i = 0; i < seq_len; i++)
    {
        if (seq[i] != '\n')
        {
            ret[j++] = seq[i];
        }
    }

    return ret;
}

void
mirbooking_sequence_set_sequence (MirbookingSequence *self, const gchar *sequence)
{
    mirbooking_sequence_set_raw_sequence (self,
                                          g_bytes_new_static (sequence, strlen (sequence)));
}

/**
 * mirbooking_sequence_get_sequence_length:
 * @self: A #MirbookingSequence
 *
 * Obtain the actual sequence length, excluding any linefeed in the internal
 * representation.
 */
gsize
mirbooking_sequence_get_sequence_length (MirbookingSequence *self)
{
    MirbookingSequencePrivate *priv = mirbooking_sequence_get_instance_private (self);
    return g_bytes_get_size (priv->sequence) - priv->sequence_skips->len;
}

/**
 * mirbooking_sequence_get_subsequence:
 * @self: A #MirbookingSequence
 * @subsequence_offset: Offset in the sequence
 * @subsequence_len: Length of the queried subsequence
 *
 * Obtain a subsequence at given offset and length.
 *
 * Note that this function is thread-safe, but not always reentrant.
 *
 * Returns: (array zero-terminated=1) (transfer none): A view of the sequence with line feed characters
 * trimmed
 */
const guint8*
mirbooking_sequence_get_subsequence (MirbookingSequence *self, gsize subsequence_offset, gsize subsequence_len)
{
    MirbookingSequencePrivate *priv = mirbooking_sequence_get_instance_private (self);

    const guint8 *subsequence = (guint8*) g_bytes_get_data (priv->sequence, NULL) + subsequence_offset;

    static __thread guint8 subsequence_buffer[64];
    gsize subsequence_buffer_offset = 0;

    g_return_val_if_fail (subsequence_offset + subsequence_len <= g_bytes_get_size (priv->sequence) - priv->sequence_skips->len, NULL);
    g_return_val_if_fail (subsequence_len <= sizeof (subsequence_buffer), NULL);

    guint i;
    for (i = 0; i < priv->sequence_skips->len; i++)
    {
        const guint8 *linefeed = g_ptr_array_index (priv->sequence_skips, i);
        if (linefeed <= subsequence)
        {
            subsequence++; // move the subsequence right
        }
        else if (G_UNLIKELY (linefeed < subsequence + subsequence_len))
        {
            // length until the next linefeed
            gsize len_to_copy = linefeed - subsequence;

            memcpy (subsequence_buffer + subsequence_buffer_offset, subsequence, len_to_copy);

            subsequence_buffer_offset += len_to_copy;
            subsequence += len_to_copy + 1; // jump right after the linefeed
        }
        else
        {
            break; // linefeed supersedes the subsequence
        }
    }

    if (G_UNLIKELY (subsequence_buffer_offset > 0))
    {
        if (subsequence_buffer_offset < subsequence_len)
        {
            memcpy (subsequence_buffer + subsequence_buffer_offset, subsequence, subsequence_len - subsequence_buffer_offset);
        }

        return subsequence_buffer;
    }
    else
    {
        return subsequence;
    }
}

static gchar
rc (gchar c)
{
    switch (toupper (c))
    {
        case 'A':
            return 'T';
        case 'C':
            return 'G';
        case 'G':
            return 'C';
        case 'T':
        case 'U':
            return 'A';
        default:
            g_assert_not_reached ();
    }
}

const guint8*
mirbooking_sequence_get_subsequence_rc (MirbookingSequence *self, gsize subsequence_offset, gsize subsequence_len)
{
    static __thread guint8 subsequence_rc_buffer[64];

    const guint8 *subsequence_buffer = mirbooking_sequence_get_subsequence (self,
                                                                            subsequence_offset,
                                                                            subsequence_len);

    gsize i;
    for (i = 0; i < subsequence_len; i++)
    {
        subsequence_rc_buffer[i] = rc (subsequence_buffer[subsequence_len - i - 1]);
    }

    return subsequence_rc_buffer;
}

static gssize
sequence_index (const guint8 *seq, gsize seq_len)
{
    gsize i;
    gsize index = 0;

    // note: 4**i = 2**2i = 1 << 2*i
    for (i = 0; i < seq_len; i++)
    {
        gsize base = 1l << (2 * i);
        switch (seq[seq_len - i - 1])
        {
            case 'A':
            case 'a':
                index += 0 * base;
                break;
            case 'C':
            case 'c':
                index += 1 * base;
                break;
            case 'G':
            case 'g':
                index += 2 * base;
                break;
            case 'T':
            case 't':
            case 'U':
            case 'u':
                index += 3 * base;
                break;
            default:
                return -1;
        }
    }

    return index;
}

/**
 * mirbooking_sequence_get_subsequence_index:
 *
 * Compute an ordinal index for the given subsequence assuming a base 4 (i.e.
 * 'A', 'C', 'G', 'T' | 'U').
 *
 * Returns: The corresponding index or %-1 if it could not be determined
 * reliably.
 */
gssize
mirbooking_sequence_get_subsequence_index (MirbookingSequence *self, gsize offset, gsize len)
{
    MirbookingSequencePrivate *priv = mirbooking_sequence_get_instance_private (self);

    g_return_val_if_fail (offset + len <= g_bytes_get_size (priv->sequence) - priv->sequence_skips->len, -1);

    return sequence_index (mirbooking_sequence_get_subsequence (self, offset, len), len);
}

guint
mirbooking_sequence_hash (const MirbookingSequence *self)
{
    MirbookingSequencePrivate *priv = mirbooking_sequence_get_instance_private ((MirbookingSequence*) self);

    return g_str_hash (priv->accession);
}

gboolean
mirbooking_sequence_equal (const MirbookingSequence *a, const MirbookingSequence *b)
{
    MirbookingSequencePrivate *priva = mirbooking_sequence_get_instance_private ((MirbookingSequence*) a);
    MirbookingSequencePrivate *privb = mirbooking_sequence_get_instance_private ((MirbookingSequence*) b);

    return g_str_equal (priva->accession, privb->accession);
}
