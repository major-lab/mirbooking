#include "mirbooking-score-index.h"

G_DEFINE_ABSTRACT_TYPE (MirbookingScoreIndex, mirbooking_score_index, G_TYPE_OBJECT)

void
mirbooking_score_index_class_init (MirbookingScoreIndexClass *klass)
{

}

void
mirbooking_score_index_init (MirbookingScoreIndex *klass)
{

}

/**
 * mirbooking_score_index_add_sequence:
 *
 * Add a #MirbookingSequence to this index such that it can be yielded by the
 * #MirbookingScoreIndexIter
 */
void
mirbooking_score_index_add_sequence (MirbookingScoreIndex *self,
                                     MirbookingSequence   *sequence,
                                     gfloat                quantity)
{
    MirbookingScoreIndexClass *klass = MIRBOOKING_SCORE_INDEX_GET_CLASS (self);

    return klass->add_sequence (self, sequence, quantity);
}

/**
 * mirbooking_score_index_iterator:
 */
MirbookingScoreIndexIter *
mirbooking_score_index_iterator (MirbookingScoreIndex *self)
{
    MirbookingScoreIndexClass *klass = MIRBOOKING_SCORE_INDEX_GET_CLASS (self);

    return klass->iterator (self);
}

typedef struct _MirbookingScoreIndexIterPrivate
{
    MirbookingScoreIndex *score_index;
} MirbookingScoreIndexIterPrivate;

G_DEFINE_TYPE_WITH_PRIVATE (MirbookingScoreIndexIter, mirbooking_score_index_iter, G_TYPE_OBJECT)

enum
{
    PROP_SCORE_INDEX = 1
};

void
get_property (GObject *object, guint property_id, GValue *value, GParamSpec *pspec)
{
    MirbookingScoreIndexIterPrivate *priv = mirbooking_score_index_iter_get_instance_private (MIRBOOKING_SCORE_INDEX_ITER (object));
    switch (property_id)
    {
        case PROP_SCORE_INDEX:
            g_value_set_object (value, priv->score_index);
            break;
        default:
            g_assert_not_reached ();
    }
}

void
set_property (GObject *object, guint property_id, const GValue *value, GParamSpec *pspec)
{
    MirbookingScoreIndexIterPrivate *priv = mirbooking_score_index_iter_get_instance_private (MIRBOOKING_SCORE_INDEX_ITER (object));
    switch (property_id)
    {
        case PROP_SCORE_INDEX:
            priv->score_index = g_value_get_object (value);
            break;
        default:
            g_assert_not_reached ();
    }
}

void
mirbooking_score_index_iter_class_init (MirbookingScoreIndexIterClass *klass)
{
    GObjectClass *object_class = G_OBJECT_CLASS (klass);
    object_class->get_property = get_property;
    object_class->set_property = set_property;

    g_object_class_install_property (G_OBJECT_CLASS (klass),
                                     PROP_SCORE_INDEX,
                                     g_param_spec_object ("score-index", "", "", MIRBOOKING_TYPE_SCORE_INDEX, G_PARAM_CONSTRUCT_ONLY | G_PARAM_READWRITE));
}

void
mirbooking_score_index_iter_init (MirbookingScoreIndexIter *self)
{

}

/**
 * Returns: (transfer none)
 */
MirbookingScoreIndex *
mirbooking_score_index_iter_get_score_index (MirbookingScoreIndexIter  *self)
{
    MirbookingScoreIndexIterPrivate *priv = mirbooking_score_index_iter_get_instance_private (MIRBOOKING_SCORE_INDEX_ITER (self));
    return priv->score_index;
}

gboolean mirbooking_score_index_iter_next (MirbookingScoreIndexIter  *self,
                                           MirbookingMirna          **mirna,
                                           MirbookingTarget         **target,
                                           gsize                     *position)
{
    g_return_val_if_fail (self != NULL, FALSE);
    g_return_val_if_fail (mirna != NULL, FALSE);
    g_return_val_if_fail (target != NULL, FALSE);
    g_return_val_if_fail (position != NULL, FALSE);

    MirbookingScoreIndexIterClass *klass = MIRBOOKING_SCORE_INDEX_ITER_GET_CLASS (self);
    return klass->next (self, mirna, target, position);
}
