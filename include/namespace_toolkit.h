// auxiliary functions

#ifndef HAMMURABI_TOOLKIT_H
#define HAMMURABI_TOOLKIT_H

#include <cassert>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include <fftw3.h>
#include <hvec.h>
#include <tinyxml2.h>

#include <cgs_units_file.h>

namespace toolkit {
// perpendicular component of a vector wrt LOS direction
// 1st argument: vector in Cartesian frame
// 2nd argument: polar angle (in rad)
// 3rd argument: azimuthal angle (in rad)
double perp2los(const hvec<3, double> &, const double &, const double &);
// (signed) parallel component of a vector wrt LOS direction
// 1st argument: vector in Cartesian frame
// 2nd argument: polar angle (in rad)
// 3rd argument: azimuthal angle (in rad)
double par2los(const hvec<3, double> &, const double &, const double &);
// intrinsic polarization angle (in rad)
// 1st argument: magnetic field in Cartesian frame
// 2nd argument: polar angle (in rad) of LOS direction
// 3rd argument: azimuthal angle (in rad) of LOS direction
// use with caution, since vector can be parallel to LOS direction
double intr_pol_ang(const hvec<3, double> &, const double &, const double &);
// Carteisan unit vector of given LOS direction
// 1st argument: polar angle (in rad)
// 2nd argument: azimuthal angle (in rad)
inline hvec<3, double> los_versor(const double &the_los,
                                  const double &phi_los) {
  return hvec<3, double>{cos(phi_los) * sin(the_los),
                         sin(phi_los) * sin(the_los), cos(the_los)};
}
// convert coordinate from Cartesian to cylindrical frame
// 1st argument: coordinate in Cartesian frame
// 2nd argument: r in cylindrical frame
// 3rd argument: phi in cylindrical frame
// 4th argument: z in cylindrical frame
void cart_coord2cyl_coord(const hvec<3, double> &, double &, double &,
                          double &);
// convert coordinate from Cartesian to cylindrical frame
// 1st argument: coordinate in Cartesian frame
// 2nd argument: vec3_t{r,phi,z} in cylindrical frame
void cart_coord2cyl_coord(const hvec<3, double> &, hvec<3, double> &);
// find index of 3D grid
inline std::size_t index3d(const std::size_t & /* n1 */, const std::size_t &n2,
                           const std::size_t &n3, const std::size_t &i,
                           const std::size_t &j, const std::size_t &l) {
  assert(j <= n2 and l <= n3);
  return (i * n2 * n3 + j * n3 + l);
}
// find index for 4D grid
inline std::size_t index4d(const std::size_t & /* n1 */, const std::size_t &n2,
                           const std::size_t &n3, const std::size_t &n4,
                           const std::size_t &e, const std::size_t &i,
                           const std::size_t &j, const std::size_t &l) {
  assert(i <= n2 and j <= n3 and l <= n4);
  return (e * n2 * n3 * n4 + i * n3 * n4 + j * n4 + l);
}
//
double mean(const double *, const std::size_t &);
//
double mean(const std::vector<double> &);
//
double variance(const double *, const std::size_t &);
//
double variance(const std::vector<double> &);
//
double covariance(const double *, const double *, const std::size_t &);
//
double covariance(const std::vector<double> &, const std::vector<double> &);
// converting brightness temp into thermal temp with T_0 = 2.725K,
// Prog.Theor.Exp.Phys. (2014) 2014 (6): 06B109.
inline double temp_convert(const double &temp_br, const double &freq) {
  const double p{CGS_U_h_planck * freq / (CGS_U_kB * 2.725)};
  return temp_br * (exp(p) - 1.) * (exp(p) - 1.) / (p * p * exp(p));
  ;
}
// use given seed number or generate random seed according to thread and clock
std::size_t random_seed(const int &);
// substract real part of a fftw complex array
// 1st argument: input fftw complex array
// 2nd argument: output double array
// 3rd argument: (output) array size
void complex2real(const fftw_complex *, double *, const std::size_t &);
// substract imaginary part of a fftw complex array
// 1st argument: input fftw complex array
// 2nd argument: output double array
// 3rd argument: (output) array size
void complex2imag(const fftw_complex *, double *, const std::size_t &);
// substract real and imaginary part of a fftw complex array
// 1st argument: input fftw complex array
// 2nd argument: output double array (real part)
// 3nd argument: output double array (imaginary part)
// 4th argument: (output) array size
void complex2rni(const fftw_complex *, double *, double *, const std::size_t &);
// load tinyxml2::XML file
// 1st argument: tinyxml2::XML file name (with dir)
std::unique_ptr<tinyxml2::XMLDocument> loadxml(const std::string &);
// trace down a key inside tinyxml2::XML "document"
// 1st argument: pointer to tinyxml2::XMLDocument
// 2nd argument: a vector of string, with key chain for tracing
// the last string in 2nd argument is the target key
// <root> is automatically included
tinyxml2::XMLElement *tracexml(tinyxml2::XMLDocument *,
                               const std::vector<std::string> &);
// get string attribute
// 1st argument: ptr to tinyxml2::XMLElement (key)
// 2nd argument: attribute name
// 3rd argument: sub-key under 1st argument
std::string fetchstring(tinyxml2::XMLElement *, const std::string &,
                        const std::string &);
// get string attribute
// 1st argument: ptr to tinyxml2::XMLElement (key)
// 2nd argument: attribute name
std::string fetchstring(tinyxml2::XMLElement *, const std::string &);
// get integer attribute
// 1st argument: ptr to tinyxml2::XMLElement (key)
// 2nd argument: attribute name
// 3rd argument: sub-key under 1st argument
int fetchint(tinyxml2::XMLElement *, const std::string &, const std::string &,
             const int &dft = 0);
// get integer attribute
// 1st argument: ptr to tinyxml2::XMLElement (key)
// 2nd argument: attribute name
int fetchint(tinyxml2::XMLElement *, const std::string &, const int &dft = 0);
// get unsigned integer attribute
// 1st argument: ptr to tinyxml2::XMLElement (key)
// 2nd argument: attribute name
// 3rd argument: sub-key under 1st argument
unsigned int fetchunsigned(tinyxml2::XMLElement *, const std::string &,
                           const std::string &, const unsigned &dft = 0);
// get unsigned integer attribute
// 1st argument: ptr to tinyxml2::XMLElement (key)
// 2nd argument: attribute name
unsigned int fetchunsigned(tinyxml2::XMLElement *, const std::string &,
                           const unsigned &dft = 0);
// get bool attribute
// 1st argument: ptr to tinyxml2::XMLElement (key)
// 2nd argument: attribute name
// 3rd argument: sub-key under 1st argument
// 4th argument: default bool, use int instead of bool type
// since bool can take char* as input
bool fetchbool(tinyxml2::XMLElement *, const std::string &, const std::string &,
               const int &dft = 0);
// get bool attribute
// 1st argument: ptr to tinyxml2::XMLElement (key)
// 2nd argument: attribute name
bool fetchbool(tinyxml2::XMLElement *, const std::string &, const int &dft = 0);
// get double attribute
// 1st argument: ptr to tinyxml2::XMLElement (key)
// 2nd argument: attribute name
// 3rd argument: sub-key under 1st argument
double fetchdouble(tinyxml2::XMLElement *, const std::string &,
                   const std::string &, const double &dft = 0);
// get double attribute
// 1st argument: ptr to tinyxml2::XMLElement (key)
// 2nd argument: attribute name
double fetchdouble(tinyxml2::XMLElement *, const std::string &,
                   const double &dft = 0);
} // namespace toolkit

#endif
