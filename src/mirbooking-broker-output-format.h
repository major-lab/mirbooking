#ifndef __MIRBOOKING_BROKER_OUTPUT_FORMAT_H__
#define __MIRBOOKING_BROKER_OUTPUT_FORMAT_H__

#include <glib-object.h>

G_BEGIN_DECLS

GType mirbooking_broker_output_format_get_type (void) G_GNUC_CONST;
#define MIRBOOKING_BROKER_OUTPUT_FORMAT_ENUM mirbooking_broker_output_format_get_type ()
typedef enum _MirbookingBrokerOutputFormat
{
    MIRBOOKING_BROKER_OUTPUT_FORMAT_TSV,
    MIRBOOKING_BROKER_OUTPUT_FORMAT_TSV_DETAILED,
    MIRBOOKING_BROKER_OUTPUT_FORMAT_GFF3,
    MIRBOOKING_BROKER_OUTPUT_FORMAT_WIG
} MirbookingBrokerOutputFormat;

G_END_DECLS

#endif /* __MIRBOOKING_BROKER_OUTPUT_FORMAT_H__ */
