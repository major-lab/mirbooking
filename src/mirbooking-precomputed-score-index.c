#include "mirbooking-precomputed-score-index.h"
#include "mirbooking-target-site.h"

typedef struct _MirbookingPrecomputedScoreIndexPrivate
{
    // the hash tables are mapping seed_index to #GSList of #MirbookingSequence
    GBytes     *score_index_bytes;
    GHashTable *mirnas_by_seed;
    GHashTable *target_sites_by_heptamer;
} MirbookingPrecomputedScoreIndexPrivate;

struct _MirbookingPrecomputedScoreIndex
{
    MirbookingScoreIndex parent_instance;
    MirbookingPrecomputedScoreIndexPrivate *priv;
};

G_DEFINE_TYPE_WITH_PRIVATE (MirbookingPrecomputedScoreIndex, mirbooking_precomputed_score_index, MIRBOOKING_TYPE_SCORE_INDEX)

typedef struct _MirbookingPrecomputedScoreIndexIterPrivate
{
    const guint16 *score_index;
    gsize          score_index_len;
    const guint16 *current_duplex_index;
    GSList        *current_mirna;
    GSList        *first_target_site;
    GSList        *current_target_site;
} MirbookingPrecomputedScoreIndexIterPrivate;

struct _MirbookingPrecomputedScoreIndexIter
{
    MirbookingScoreIndexIter                    parent_instance;
    MirbookingPrecomputedScoreIndexIterPrivate *priv;
};

G_DEFINE_TYPE_WITH_PRIVATE (MirbookingPrecomputedScoreIndexIter, mirbooking_precomputed_score_index_iter, MIRBOOKING_TYPE_SCORE_INDEX_ITER)

MirbookingPrecomputedScoreIndexIter *
mirbooking_precomputed_score_index_iter_new (MirbookingPrecomputedScoreIndex *self)
{
    return g_object_new (MIRBOOKING_TYPE_PRECOMPUTED_SCORE_INDEX_ITER,
                         "score-index", self,
                         NULL);
}

MirbookingPrecomputedScoreIndex *
mirbooking_precomputed_score_index_new (guint16 *data)
{
    return g_object_new (MIRBOOKING_TYPE_PRECOMPUTED_SCORE_INDEX,
                         "score-index", g_bytes_new_take (data, 16384 * 16384 * 2 * sizeof (guint16)),
                         NULL);
}

MirbookingPrecomputedScoreIndex *
mirbooking_precomputed_score_index_new_from_bytes (GBytes *data)
{
    return g_object_new (MIRBOOKING_TYPE_PRECOMPUTED_SCORE_INDEX,
                         "score-index", data,
                         NULL);
}

static void
add_sequence (MirbookingScoreIndex *self,
              MirbookingSequence   *sequence,
              gfloat                quantity)
{
    MirbookingPrecomputedScoreIndexPrivate *priv = MIRBOOKING_PRECOMPUTED_SCORE_INDEX (self)->priv;

    if (MIRBOOKING_IS_MIRNA (sequence))
    {
        gsize seed_index = mirbooking_sequence_get_subsequence_index (sequence, 1, 7);

        GSList *mirnas_list = g_hash_table_lookup (priv->mirnas_by_seed,
                                                   GSIZE_TO_POINTER (seed_index));

        g_hash_table_insert (priv->mirnas_by_seed,
                             GSIZE_TO_POINTER (seed_index),
                             g_slist_prepend (mirnas_list, g_object_ref (sequence)));
    }
    else
    {
        gsize heptamer_index;
        gint i;
        for (i = 0; i < mirbooking_sequence_get_sequence_length (sequence) - 7; i++)
        {
            heptamer_index = mirbooking_sequence_get_subsequence_index (sequence, i, 7);
            GSList *target_sites_list = g_hash_table_lookup (priv->target_sites_by_heptamer,
                                                             GSIZE_TO_POINTER (heptamer_index));

            // create a #MirbookingTargetSite to keep the target and the position
            MirbookingTargetSite *target_site = g_new0 (MirbookingTargetSite, 1);
            target_site->target   = g_object_ref (sequence);
            target_site->position = i;

            g_hash_table_insert (priv->target_sites_by_heptamer,
                                 GSIZE_TO_POINTER (heptamer_index),
                                 g_slist_prepend (target_sites_list, target_site));
        }
    }
}

static MirbookingScoreIndexIter *
iterator (MirbookingScoreIndex *self)
{
    return MIRBOOKING_SCORE_INDEX_ITER (mirbooking_precomputed_score_index_iter_new (MIRBOOKING_PRECOMPUTED_SCORE_INDEX (self)));
}

enum
{
    PROP_SCORE_INDEX = 1
};

static void
get_property (GObject *object, guint property_id, GValue *value, GParamSpec *pspec)
{
    switch (property_id)
    {
        case PROP_SCORE_INDEX:
            g_value_set_boxed (value, MIRBOOKING_PRECOMPUTED_SCORE_INDEX (object)->priv->score_index_bytes);
            break;
        default:
            g_assert_not_reached ();
    }
}

static void
set_property (GObject *object, guint property_id, const GValue *value, GParamSpec *pspec)
{
    MirbookingPrecomputedScoreIndex *self = MIRBOOKING_PRECOMPUTED_SCORE_INDEX (object);

    GBytes *score_index_bytes;

    switch (property_id)
    {
        case PROP_SCORE_INDEX:
            score_index_bytes = g_value_get_boxed (value);
            g_return_if_fail (g_bytes_get_size (score_index_bytes) == 16384 * 16384 * 2 * sizeof (guint16));
            self->priv->score_index_bytes = g_bytes_ref (score_index_bytes);
            break;
        default:
            g_assert_not_reached ();
    }

}

void
mirbooking_precomputed_score_index_class_init (MirbookingPrecomputedScoreIndexClass *klass)
{
    klass->parent_class.add_sequence = add_sequence;
    klass->parent_class.iterator     = iterator;

    GObjectClass *object_class = G_OBJECT_CLASS (klass);

    object_class->get_property = get_property;
    object_class->set_property = set_property;

    g_object_class_install_property (object_class,
                                     PROP_SCORE_INDEX,
                                     g_param_spec_boxed ("score-index", "", "", G_TYPE_BYTES, G_PARAM_CONSTRUCT | G_PARAM_READWRITE));
}

void
mirbooking_precomputed_score_index_init (MirbookingPrecomputedScoreIndex *self)
{
    self->priv = g_new (MirbookingPrecomputedScoreIndexPrivate, 1);

    self->priv->target_sites_by_heptamer = g_hash_table_new (g_direct_hash,
                                                             g_direct_equal);

    self->priv->mirnas_by_seed = g_hash_table_new (g_direct_hash,
                                                   g_direct_equal);
}

void
mirbooking_precomputed_score_index_finalize (MirbookingPrecomputedScoreIndex *self)
{
    g_hash_table_unref (self->priv->mirnas_by_seed);
    g_hash_table_unref (self->priv->target_sites_by_heptamer);
    g_free (self->priv);
}

static void
constructed (GObject *object)
{
    MirbookingPrecomputedScoreIndexIter *self                    = MIRBOOKING_PRECOMPUTED_SCORE_INDEX_ITER (object);
    MirbookingPrecomputedScoreIndex     *precomputed_score_index = MIRBOOKING_PRECOMPUTED_SCORE_INDEX (mirbooking_score_index_iter_get_score_index (MIRBOOKING_SCORE_INDEX_ITER (self)));

    self->priv = g_new0 (MirbookingPrecomputedScoreIndexIterPrivate, 1);
    gsize len;
    self->priv->score_index = g_bytes_get_data (precomputed_score_index->priv->score_index_bytes,
                                                &len);
    self->priv->score_index_len = len / sizeof (guint16);
    self->priv->current_duplex_index = self->priv->score_index;
}

static gboolean
next (MirbookingScoreIndexIter  *score_index_iter,
      MirbookingMirna          **mirna,
      MirbookingTarget         **target,
      gsize                     *position)
{
    MirbookingPrecomputedScoreIndexIter *self                    = MIRBOOKING_PRECOMPUTED_SCORE_INDEX_ITER (score_index_iter);
    MirbookingPrecomputedScoreIndex     *precomputed_score_index = MIRBOOKING_PRECOMPUTED_SCORE_INDEX (mirbooking_score_index_iter_get_score_index (MIRBOOKING_SCORE_INDEX_ITER (self)));

    // find a duplex with candidates
    if (self->priv->current_mirna == NULL)
    {
        while (self->priv->current_duplex_index < self->priv->score_index + self->priv->score_index_len)
        {
            // fetch the next set of candidates
            gsize seed_index, heptamer_index;

            seed_index     = (gsize) GUINT16_FROM_BE (self->priv->current_duplex_index[0]);
            heptamer_index = (gsize) GUINT16_FROM_BE (self->priv->current_duplex_index[1]);

            self->priv->current_mirna = g_hash_table_lookup (precomputed_score_index->priv->mirnas_by_seed,
                                                            GSIZE_TO_POINTER (seed_index));

            self->priv->current_target_site = g_hash_table_lookup (precomputed_score_index->priv->target_sites_by_heptamer,
                                                                GSIZE_TO_POINTER (heptamer_index));
            self->priv->first_target_site = self->priv->current_target_site;

            self->priv->current_duplex_index += 2;

            if (self->priv->current_mirna != NULL && self->priv->current_target_site != NULL)
            {
                g_assert (self->priv->current_mirna->data != NULL);
                g_assert (self->priv->current_target_site->data != NULL);
                break;
            }
        }
    }

    // all duplexes have been tested
    if (self->priv->current_duplex_index >= self->priv->score_index + self->priv->score_index_len)
    {
        return FALSE;
    }

    MirbookingTargetSite *target_site = self->priv->current_target_site->data;

    *mirna    = self->priv->current_mirna->data;
    *target   = target_site->target;
    *position = target_site->position;

    self->priv->current_target_site = self->priv->current_target_site->next;

    if (self->priv->current_target_site == NULL)
    {
        self->priv->current_mirna       = self->priv->current_mirna->next;
        self->priv->current_target_site = self->priv->first_target_site;
    }


    return TRUE;
}

void
mirbooking_precomputed_score_index_iter_class_init (MirbookingPrecomputedScoreIndexIterClass *klass)
{
    GObjectClass *object_class = G_OBJECT_CLASS (klass);
    object_class->constructed = constructed;
    klass->parent_class.next = next;
}

void
mirbooking_precomputed_score_index_iter_init (MirbookingPrecomputedScoreIndexIter *score_index_iter)
{
}
