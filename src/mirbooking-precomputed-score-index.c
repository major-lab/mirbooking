#include "mirbooking-precomputed-score-index.h"
#include "mirbooking-precomputed-score-table.h"

typedef struct _MirbookingPrecomputedScoreIndexPrivate
{
    // the hash tables are mapping seed_index to #GSList of #MirbookingSequence
    gsize       seed_offset;
    gsize       seed_len;
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
    const guint16 *current_duplex;
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
mirbooking_precomputed_score_index_new (guint16 *data, gsize seed_offset, gsize seed_len)
{
    return g_object_new (MIRBOOKING_TYPE_PRECOMPUTED_SCORE_INDEX,
                         "score-index", g_bytes_new_take (data, 16384 * 16384 * 2 * sizeof (guint16)),
                         "seed-offset", seed_offset,
                         "seed-length", seed_len,
                         NULL);
}

MirbookingPrecomputedScoreIndex *
mirbooking_precomputed_score_index_new_from_bytes (GBytes *data, gsize seed_offset, gsize seed_len)
{
    return g_object_new (MIRBOOKING_TYPE_PRECOMPUTED_SCORE_INDEX,
                         "score-index", data,
                         "seed-offset", seed_offset,
                         "seed-length", seed_len,
                         NULL);
}

typedef struct _MirbookingQuantifiedSequence
{
    MirbookingSequence *sequence;
    gfloat              quantity;
} MirbookingQuantifiedSequence;

static gint
mirbooking_quantified_sequence_quantity_cmp_desc (const MirbookingQuantifiedSequence *a, const MirbookingQuantifiedSequence *b)
{
    return (a->quantity < b->quantity) - (a->quantity > b->quantity);
}

/*
 * Here we use a more quantified storage for the target site than
 * #MirbookingTargetSite
 */
typedef struct _MirbookingQuantifiedTargetSite
{
    MirbookingQuantifiedSequence sequence;
    gsize                        position;
} MirbookingQuantifiedTargetSite;

typedef struct _MirbookingQuantifiedMirna
{
    MirbookingQuantifiedSequence sequence;
} MirbookingQuantifiedMirna;

static void
set_sequence_quantity (MirbookingScoreIndex *self,
                       MirbookingSequence   *sequence,
                       gfloat                quantity)
{
    MirbookingPrecomputedScoreIndexPrivate *priv = MIRBOOKING_PRECOMPUTED_SCORE_INDEX (self)->priv;

    if (MIRBOOKING_IS_MIRNA (sequence))
    {
        gsize seed_index = mirbooking_sequence_get_subsequence_index (sequence, priv->seed_offset, priv->seed_len);

        GSList *mirnas_list = g_hash_table_lookup (priv->mirnas_by_seed,
                                                   GSIZE_TO_POINTER (seed_index));

        g_hash_table_steal (priv->mirnas_by_seed, GSIZE_TO_POINTER (seed_index));

        MirbookingQuantifiedMirna *mirna = g_new (MirbookingQuantifiedMirna, 1);
        mirna->sequence.sequence = g_object_ref (sequence);
        mirna->sequence.quantity = quantity;

        g_hash_table_insert (priv->mirnas_by_seed,
                             GSIZE_TO_POINTER (seed_index),
                             g_slist_insert_sorted (mirnas_list, mirna, (GCompareFunc) mirbooking_quantified_sequence_quantity_cmp_desc));
    }
    else
    {
        gsize heptamer_index;
        gint i;
        for (i = 0; i < mirbooking_sequence_get_sequence_length (sequence) - priv->seed_len; i++)
        {
            heptamer_index = mirbooking_sequence_get_subsequence_index (sequence, i, priv->seed_len);
            GSList *target_sites_list = g_hash_table_lookup (priv->target_sites_by_heptamer,
                                                             GSIZE_TO_POINTER (heptamer_index));

            // create a #MirbookingQuantifiedTargetSite to keep the target and the position
            MirbookingQuantifiedTargetSite *target_site = g_new0 (MirbookingQuantifiedTargetSite, 1);
            target_site->sequence.sequence = g_object_ref (sequence);
            target_site->position = i;
            target_site->sequence.quantity = quantity;

            g_hash_table_steal (priv->target_sites_by_heptamer, GSIZE_TO_POINTER (heptamer_index));

            g_hash_table_insert (priv->target_sites_by_heptamer,
                                 GSIZE_TO_POINTER (heptamer_index),
                                 g_slist_insert_sorted (target_sites_list, target_site, (GCompareFunc) mirbooking_quantified_sequence_quantity_cmp_desc));
        }
    }
}

static MirbookingScoreIndexIter *
iterator (MirbookingScoreIndex *self)
{
    return MIRBOOKING_SCORE_INDEX_ITER (mirbooking_precomputed_score_index_iter_new (g_object_ref (self)));
}

enum
{
    PROP_SCORE_INDEX = 1,
    PROP_SEED_OFFSET,
    PROP_SEED_LENGTH
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
        case PROP_SEED_OFFSET:
            self->priv->seed_offset = g_value_get_uint (value);
            break;
        case PROP_SEED_LENGTH:
            self->priv->seed_len = g_value_get_uint (value);
            break;
        default:
            g_assert_not_reached ();
    }

}

static void
finalize (GObject *object)
{
    MirbookingPrecomputedScoreIndex *self = MIRBOOKING_PRECOMPUTED_SCORE_INDEX (object);

    g_bytes_unref (self->priv->score_index_bytes);
    g_hash_table_unref (self->priv->mirnas_by_seed);
    g_hash_table_unref (self->priv->target_sites_by_heptamer);
    g_free (self->priv);

    G_OBJECT_CLASS (mirbooking_precomputed_score_index_parent_class)->finalize (object);
}

void
mirbooking_precomputed_score_index_class_init (MirbookingPrecomputedScoreIndexClass *klass)
{
    GObjectClass *object_class = G_OBJECT_CLASS (klass);

    object_class->get_property = get_property;
    object_class->set_property = set_property;
    object_class->finalize     = finalize;

    klass->parent_class.set_sequence_quantity = set_sequence_quantity;
    klass->parent_class.iterator              = iterator;

    g_object_class_install_property (object_class,
                                     PROP_SCORE_INDEX,
                                     g_param_spec_boxed ("score-index", "", "", G_TYPE_BYTES, G_PARAM_CONSTRUCT | G_PARAM_READWRITE));
    g_object_class_install_property (object_class,
                                     PROP_SEED_OFFSET,
                                     g_param_spec_uint ("seed-offset", "", "", 0, G_MAXUINT, MIRBOOKING_PRECOMPUTED_SCORE_TABLE_DEFAULT_SEED_OFFSET, G_PARAM_CONSTRUCT | G_PARAM_READWRITE));
    g_object_class_install_property (object_class,
                                     PROP_SEED_LENGTH,
                                     g_param_spec_uint ("seed-length", "", "", 1, G_MAXUINT, MIRBOOKING_PRECOMPUTED_SCORE_TABLE_DEFAULT_SEED_LENGTH, G_PARAM_CONSTRUCT | G_PARAM_READWRITE));
}

static void
free_quantified_sequence (MirbookingQuantifiedSequence *sequence)
{
    g_object_unref (sequence->sequence);
    g_free (sequence);
}
static void
free_slist_of_quantified_sequences (GSList *list)
{
    g_slist_free_full (list,
                       (GDestroyNotify) free_quantified_sequence);
}

void
mirbooking_precomputed_score_index_init (MirbookingPrecomputedScoreIndex *self)
{
    self->priv = g_new (MirbookingPrecomputedScoreIndexPrivate, 1);

    self->priv->target_sites_by_heptamer = g_hash_table_new_full (g_direct_hash,
                                                                  g_direct_equal,
                                                                  NULL,
                                                                  (GDestroyNotify) free_slist_of_quantified_sequences);

    self->priv->mirnas_by_seed = g_hash_table_new_full (g_direct_hash,
                                                        g_direct_equal,
                                                        NULL,
                                                        (GDestroyNotify) free_slist_of_quantified_sequences);
}

static void
iter_constructed (GObject *object)
{
    MirbookingPrecomputedScoreIndexIter *self                    = MIRBOOKING_PRECOMPUTED_SCORE_INDEX_ITER (object);
    MirbookingPrecomputedScoreIndex     *precomputed_score_index = MIRBOOKING_PRECOMPUTED_SCORE_INDEX (mirbooking_score_index_iter_get_score_index (MIRBOOKING_SCORE_INDEX_ITER (self)));

    self->priv = g_new0 (MirbookingPrecomputedScoreIndexIterPrivate, 1);
    gsize len;
    self->priv->score_index = g_bytes_get_data (precomputed_score_index->priv->score_index_bytes,
                                                &len);
    self->priv->score_index_len = len / sizeof (guint16);
    self->priv->current_duplex = self->priv->score_index;

    G_OBJECT_CLASS (mirbooking_precomputed_score_index_iter_parent_class)->constructed (object);
}

static void
iter_finalize (GObject *object)
{
    MirbookingPrecomputedScoreIndexIter *self = MIRBOOKING_PRECOMPUTED_SCORE_INDEX_ITER (object);

    g_free (self->priv);

    G_OBJECT_CLASS (mirbooking_precomputed_score_index_iter_parent_class)->finalize (object);
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
        while (self->priv->current_duplex < self->priv->score_index + self->priv->score_index_len)
        {
            // fetch the next set of candidates
            gsize seed_index, heptamer_index;

            seed_index     = (gsize) GUINT16_FROM_BE (self->priv->current_duplex[0]);
            heptamer_index = (gsize) GUINT16_FROM_BE (self->priv->current_duplex[1]);

            self->priv->current_mirna = g_hash_table_lookup (precomputed_score_index->priv->mirnas_by_seed,
                                                            GSIZE_TO_POINTER (seed_index));

            self->priv->current_target_site = g_hash_table_lookup (precomputed_score_index->priv->target_sites_by_heptamer,
                                                                GSIZE_TO_POINTER (heptamer_index));
            self->priv->first_target_site = self->priv->current_target_site;

            self->priv->current_duplex += 2;

            if (self->priv->current_mirna != NULL && self->priv->current_target_site != NULL)
            {
                g_assert (self->priv->current_mirna->data != NULL);
                g_assert (self->priv->current_target_site->data != NULL);
                break;
            }
        }
    }

    // all duplexes have been tested
    if (self->priv->current_duplex >= self->priv->score_index + self->priv->score_index_len)
    {
        return FALSE;
    }

    MirbookingQuantifiedMirna *quantified_mirna = self->priv->current_mirna->data;
    MirbookingQuantifiedTargetSite *target_site = self->priv->current_target_site->data;

    *mirna    = MIRBOOKING_MIRNA (quantified_mirna->sequence.sequence);
    *target   = MIRBOOKING_TARGET (target_site->sequence.sequence);
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

    object_class->constructed = iter_constructed;
    object_class->finalize    = iter_finalize;

    klass->parent_class.next = next;
}

void
mirbooking_precomputed_score_index_iter_init (MirbookingPrecomputedScoreIndexIter *score_index_iter)
{
}
