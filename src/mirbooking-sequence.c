#include "mirbooking-sequence.h"

#include <string.h>
#include <glib-object.h>

typedef struct
{
    gchar       *accession;
    const gchar *sequence;
    gsize        sequence_len;
    GPtrArray   *sequence_skips; /* all the pointers in 'sequence' to skip (i.e. line feeds) */
} MirbookingSequencePrivate;

enum
{
    PROP_ACCESSION = 1,
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

    switch (property_id)
    {
        case PROP_ACCESSION:
            priv->accession = g_value_dup_string ((value));
            break;
        case PROP_SEQUENCE:
            mirbooking_sequence_set_raw_sequence (MIRBOOKING_SEQUENCE (object),
                                                  g_value_dup_string (value),
                                                  strlen (g_value_get_string (value)));
            break;
        default:
            g_assert_not_reached ();
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
        case PROP_SEQUENCE:
            g_value_set_string (value, priv->sequence);
            break;
        default:
            g_assert_not_reached ();
    }
}

static void
mirbooking_sequence_finalize (GObject *object)
{
    MirbookingSequencePrivate *priv = mirbooking_sequence_get_instance_private (MIRBOOKING_SEQUENCE (object));

    g_free (priv->accession);
    g_ptr_array_unref (priv->sequence_skips);

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
                                     PROP_SEQUENCE,
                                     g_param_spec_string ("sequence", "Sequence", "", NULL, G_PARAM_READWRITE));
}

const gchar *
mirbooking_sequence_get_accession (MirbookingSequence *self)
{
    MirbookingSequencePrivate *priv = mirbooking_sequence_get_instance_private (self);

    return priv->accession;
}

/**
 * mirbooking_sequence_get_raw_sequence:
 * @self: A #MirbookingSequence
 * @sequence_len: (out): (optional): The length of the returned raw sequence
 *
 * Obtain the internal representation of the sequence, which corresponds to
 * what has been set previsouly via #mirbooking_sequence_set_raw_sequence.
 */
const gchar *
mirbooking_sequence_get_raw_sequence (MirbookingSequence *self, gsize *sequence_len)
{
MirbookingSequencePrivate *priv = mirbooking_sequence_get_instance_private (self);

    if (sequence_len)
    {
        *sequence_len = priv->sequence_len;
    }

    return priv->sequence;
}

/**
 * mirbooking_sequence_set_raw_sequence:
 * @self: A #MirbookingSequence
 * @sequence: (array length=sequence_len): A sequence of nucleotides
 * @sequence_len: The length of the passed sequence or -1 for a null-terminated
 * string
 *
 * Note that the passed sequence is not copied internally and thus, it must
 * live for as long as this object does.
 */
void
mirbooking_sequence_set_raw_sequence (MirbookingSequence *self, const gchar *sequence, gssize sequence_len)
{
    MirbookingSequencePrivate *priv = mirbooking_sequence_get_instance_private (self);

    priv->sequence     = sequence;
    priv->sequence_len = sequence_len == -1 ? strlen (sequence) : sequence_len;

    if (priv->sequence_skips)
    {
        g_ptr_array_unref (priv->sequence_skips);
    }

    priv->sequence_skips = g_ptr_array_sized_new (sequence_len / 80); // line feed every 80 characters

    const gchar* seq = sequence;
    while ((seq = memchr (seq, '\n', sequence_len - (seq - sequence))))
    {
        g_ptr_array_add (priv->sequence_skips, (gpointer) seq);
        seq++; // jump right after the line feed
    }
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
    return priv->sequence_len - priv->sequence_skips->len;
}

/**
 * mirbooking_sequence_get_subsequence:
 *
 * Returns: (transfer none): A view of the sequence with line feed characters
 * trimmed
 */
const gchar*
mirbooking_sequence_get_subsequence (MirbookingSequence *self, gsize subsequence_offset, gsize subsequence_len)
{
    MirbookingSequencePrivate *priv = mirbooking_sequence_get_instance_private (self);
    gboolean subsequence_used = FALSE;

    static gchar subsequence[64];

    gsize subsequence_so_far = 0;

    g_return_val_if_fail (subsequence_offset + subsequence_len <= priv->sequence_len - priv->sequence_skips->len, NULL);
    g_return_val_if_fail (subsequence_len <= sizeof (subsequence), NULL);

    gint i;
    for (i = 0; i < priv->sequence_skips->len; i++)
    {
        gchar *linefeed = priv->sequence_skips->pdata[i];
        gsize  linefeed_offset = linefeed - priv->sequence;
        if (linefeed_offset <= subsequence_offset)
        {
            subsequence_offset++; // move the subsequence right
        }
        else if (linefeed_offset < subsequence_offset + subsequence_len)
        {
            // length until the next linefeed
            gsize len_to_copy = linefeed_offset - subsequence_offset;

            memcpy (subsequence + subsequence_so_far, priv->sequence + subsequence_offset, len_to_copy);

            subsequence_so_far += len_to_copy;
            subsequence_offset = subsequence_offset + len_to_copy + 1; // jump right after the linefeed
            subsequence_used = TRUE;
        }
        else
        {
            break; // linefeed supersedes the subsequence
        }
    }

    if (subsequence_used)
    {
        if (subsequence_so_far < subsequence_len)
        {
            memcpy (subsequence + subsequence_so_far, priv->sequence + subsequence_offset, subsequence_len - subsequence_so_far);
        }

        return subsequence;
    }
    else
    {
        return priv->sequence + subsequence_offset;
    }
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
