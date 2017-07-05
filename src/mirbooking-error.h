#ifndef __MIRBOOKING_ERROR_H__
#define __MIRBOOKING_ERROR_H__

#include <glib.h>

/**
 * MirbookingError:
 * @MIRBOOKING_ERROR_FAILED: Unknown or unspecified error
 */
typedef enum
{
    MIRBOOKING_ERROR_FAILED = 0
} MirbookingError;

GQuark mirbooking_error_quark (void);
#define MIRBOOKING_ERROR mirbooking_error_quark ()

#endif /* __MIRBOOKING_ERROR_H__ */
