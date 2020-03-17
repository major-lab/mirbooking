#include "mirbooking-broker-output-format.h"

GType
mirbooking_broker_output_format_get_type (void)
{
    static gsize type_id_once;

    if (g_once_init_enter (&type_id_once))
    {
        static const GEnumValue values[] =
        {
            {0, "TSV",          "tsv"},
            {1, "TSV_DETAILED", "tsv-detailed"},
            {2, "GFF3",         "gff3"},
            {3, "WIG",          "wig"},
            {0, NULL,           NULL}
        };

        GType type = g_enum_register_static ("MirbookingBrokerOutputFormat",
                                             values);
        g_once_init_leave (&type_id_once, type);
    }

    return type_id_once;
}
