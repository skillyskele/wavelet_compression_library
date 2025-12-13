#ifndef COMPRESSION_TYPES_H
#define COMPRESSION_TYPES_H


#include <stdint.h>


typedef float COEFFICIENT_TYPE;
typedef int16_t SAMPLE_TYPE; // will actually be uint16 with the iadc, but for testing, it's easier to have it float
#define SIGNAL_LENGTH 32
#define COEFFICIENTS_LENGTH 42




#endif // COMPRESSION_TYPES_H

