// this customized header file allows
// ajdusting data types for working enviroment/project
//
// brief idea of how many bits to use
//
// 64bit is out of our reach, do not use it
//
// 32bit unsigned integer goes upto 4,294,967,296
// that's sufficient for HEALPix Nsise 18,918
//
// 16bit unsigned integer goes upto 65,536
// sufficent for HEALPix Nside 73

#ifndef HAMMURABI_TYPE_H
#define HAMMURABI_TYPE_H

#include <cstdint>

typedef std::uint_fast32_t ham_uint; // unsigned int
typedef std::int_fast32_t ham_int;   // singed int
typedef std::int_fast8_t ham_sint;   // short int
typedef double ham_float;            // floating piont

#endif
