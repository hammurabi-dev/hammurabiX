// auxiliary functions

#ifndef HAMMURABI_TOOLKIT_H
#define HAMMURABI_TOOLKIT_H

#include <cstddef> // for std::size_t
#include <cassert>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include <hamvec.h>
#include <tinyxml2.h>

namespace toolkit {
// find index of 3D grid
// by convention, dimensions are ordered as x-y-z
// 1st argument: number of vertices in 1st dimension (not used)
// 2nd argument: number of vertices in 2nd direction
// 3rd argument: number of vertices in 3rd direction
// 4th argument: index for 1st dimension
// 5th argument: index of 2nd dimension
// 6th argument: index for 3rd dimenstion
inline std::size_t index3d(const std::size_t & /* n1 */, const std::size_t &n2,
                           const std::size_t &n3, const std::size_t &i,
                           const std::size_t &j, const std::size_t &l) {
  assert(j <= n2 and l <= n3);
  return (i * n2 * n3 + j * n3 + l);
}
// find index for 4D grid
// following the same idea as index3d function
// but adding the 4th dimension
inline std::size_t index4d(const std::size_t & /* n1 */, const std::size_t &n2,
                           const std::size_t &n3, const std::size_t &n4,
                           const std::size_t &e, const std::size_t &i,
                           const std::size_t &j, const std::size_t &l) {
  assert(i <= n2 and j <= n3 and l <= n4);
  return (e * n2 * n3 * n4 + i * n3 * n4 + j * n4 + l);
}
// mean value estimator for double array
// 1st argument: pointer for a double array
// 2nd argument: size of array
double mean(const double *, const std::size_t &);
// mean value estimator for double vector
// 1st argument: ref for a double vector
double mean(const std::vector<double> &);
// variance estimator for double array
// 1st argument: pointer for a double array
// 2nd argument: size of array
double variance(const double *, const std::size_t &);
// variance estimator for double vector
// 1st argument: ref for a double vector
double variance(const std::vector<double> &);
// covariance estimator for double arrays
// 1st argument: pointer for 1st double array
// 2nd argument: pointer for 2nd double array
// 3rd argument: size of array
double covariance(const double *, const double *, const std::size_t &);
// covariance estimator double vectors
// 1st argument: ref for 1st double vector
// 2nd argument: ref for 2nd double vector
double covariance(const std::vector<double> &, const std::vector<double> &);
// use given seed number or generate random seed according to thread and clock
std::size_t random_seed(const int &);
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
int fetchint(tinyxml2::XMLElement *, const std::string &, const std::string &);
// get integer attribute
// 1st argument: ptr to tinyxml2::XMLElement (key)
// 2nd argument: attribute name
int fetchint(tinyxml2::XMLElement *, const std::string &);
// get unsigned integer attribute
// 1st argument: ptr to tinyxml2::XMLElement (key)
// 2nd argument: attribute name
// 3rd argument: sub-key under 1st argument
unsigned int fetchunsigned(tinyxml2::XMLElement *, const std::string &,
                           const std::string &);
// get unsigned integer attribute
// 1st argument: ptr to tinyxml2::XMLElement (key)
// 2nd argument: attribute name
unsigned int fetchunsigned(tinyxml2::XMLElement *, const std::string &);
// get bool attribute
// 1st argument: ptr to tinyxml2::XMLElement (key)
// 2nd argument: attribute name
// 3rd argument: sub-key under 1st argument
// 4th argument: default bool, use int instead of bool type
// since bool can take char* as input
bool fetchbool(tinyxml2::XMLElement *, const std::string &,
               const std::string &);
// get bool attribute
// 1st argument: ptr to tinyxml2::XMLElement (key)
// 2nd argument: attribute name
bool fetchbool(tinyxml2::XMLElement *, const std::string &);
// get double attribute
// 1st argument: ptr to tinyxml2::XMLElement (key)
// 2nd argument: attribute name
// 3rd argument: sub-key under 1st argument
double fetchdouble(tinyxml2::XMLElement *, const std::string &,
                   const std::string &);
// get double attribute
// 1st argument: ptr to tinyxml2::XMLElement (key)
// 2nd argument: attribute name
double fetchdouble(tinyxml2::XMLElement *, const std::string &);
} // namespace toolkit

#endif
