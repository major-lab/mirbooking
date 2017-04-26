#include "mirbooking-sequence.h"

#include <string.h>
#include <glib-object.h>

/*
 * The 'sequence' field is not owned by the instance and it's corresponding
 * memory is meant to be managed externally.
 *
 * This support *raw* sequences that may contain line feeds like those found in
 * a FASTA document.
 */

typedef struct
{
    gchar       *accession;
    const gchar *sequence;
    gsize        sequence_len;
    GPtrArray   *sequence_skips; /* all the pointers in 'sequence' to skip (i.e. line feeds) */
} MirbookingSequencePrivate;

enum
{
    PROP_ACCESSION = 1
};

G_DEFINE_TYPE_WITH_PRIVATE (MirbookingSequence, mirbooking_sequence, G_TYPE_OBJECT)

static void
mirbooking_sequence_init (MirbookingSequence *self)
{
    MirbookingSequencePrivate *priv = mirbooking_sequence_get_instance_private (MIRBOOKING_SEQUENCE (self));

    priv->sequence_skips = g_ptr_array_new ();
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
}

const gchar *
mirbooking_sequence_get_accession (MirbookingSequence *self)
{
    MirbookingSequencePrivate *priv = mirbooking_sequence_get_instance_private (self);

    return priv->accession;
}

const gchar *
mirbooking_sequence_get_sequence (MirbookingSequence *self, gsize *sequence_len)
{
    MirbookingSequencePrivate *priv = mirbooking_sequence_get_instance_private (self);

    *sequence_len = priv->sequence_len;
    return priv->sequence;
}

/**
 * mirbooking_sequence_set_sequence:
 *
 * Note that the passed sequence is not copied internally, unless it contains
 * line feeds.
 * #sequence: (transfer none)
 */
void
mirbooking_sequence_set_sequence (MirbookingSequence *self, const gchar *sequence, gsize sequence_len)
{
    MirbookingSequencePrivate *priv = mirbooking_sequence_get_instance_private (self);

    priv->sequence     = sequence;
    priv->sequence_len = sequence_len;

    // g_ptr_array_clear (priv->sequence_skips);

    const gchar* seq = priv->sequence;
    while ((seq = memchr (seq, '\n', priv->sequence_len)))
    {
        g_ptr_array_add (priv->sequence_skips, (gpointer) seq);
    }
}

/**
 * mirbooking_sequence_get_subsequence:
 *
 * Return: (transfer-none): a view of the sequence with line feed characters
 * trimmed
 */
const gchar*
mirbooking_sequence_get_subsequence (MirbookingSequence *self, gsize offset, gsize len)
{
    MirbookingSequencePrivate *priv = mirbooking_sequence_get_instance_private (self);

    static gchar subsequence[64];

    return subsequence;

    /*

    g_return_val_if_fail (len <= sizeof (subsequence), NULL);
    g_return_val_if_fail (offset + len < priv->sequence_len, NULL);

    gchar *trimmed_sequence = g_malloc (priv->sequence_len);

    GSList *list;
    for (list = priv->sequence_skips; list != NULL; list = list->next)
    {
    }

    if (memchr (priv->sequence, '\n', priv->sequence_len))
    {
        g_debug ("The sequence with accession '%s' contains line feeds and will be copied internally.", priv->accession);

        const gchar *original_sequence = priv->sequence;
        gsize  original_sequence_len   = priv->sequence_len;

        gchar *line_feed;
        ptrdiff_t offset = 0;
        while ((line_feed = memchr (original_sequence, '\n', original_sequence_len)))
        {
            offset = line_feed - original_sequence;
            memcpy (trimmed_sequence + offset, original_sequence, original_sequence_len - offset);

            original_sequence      = line_feed + 1; // skip the line feed
            original_sequence_len -= line_feed - original_sequence - 1;

            // drop the line feed
            priv->sequence_len--;
        }

        // copy any remaining nucleotides
        if (original_sequence_len > 0)
        {
            memcpy (trimmed_sequence, original_sequence, original_sequence_len);
        }

        // mark it to be freed internally
        priv->sequence = trimmed_sequence;
    }
    else
    {
        return priv->sequence + offset;
    }
    */
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
