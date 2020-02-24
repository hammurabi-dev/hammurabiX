// auxiliary functions

#ifndef HAMMURABI_TOOLKIT_H
#define HAMMURABI_TOOLKIT_H

#include <cassert>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include <hamtype.h>
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
inline ham_uint index3d(const ham_uint & /* n1 */, const ham_uint &n2,
                        const ham_uint &n3, const ham_uint &i,
                        const ham_uint &j, const ham_uint &l) {
  assert(j <= n2 and l <= n3);
  return (i * n2 * n3 + j * n3 + l);
}
// find index for 4D grid
// following the same idea as index3d function
// but adding the 4th dimension
inline ham_uint index4d(const ham_uint & /* n1 */, const ham_uint &n2,
                        const ham_uint &n3, const ham_uint &n4,
                        const ham_uint &e, const ham_uint &i, const ham_uint &j,
                        const ham_uint &l) {
  assert(i <= n2 and j <= n3 and l <= n4);
  return (e * n2 * n3 * n4 + i * n3 * n4 + j * n4 + l);
}
// mean value estimator for ham_float array
// 1st argument: pointer for a ham_float array
// 2nd argument: size of array
ham_float mean(const ham_float *, const ham_uint &);
// mean value estimator for ham_float vector
// 1st argument: ref for a ham_float vector
ham_float mean(const std::vector<ham_float> &);
// variance estimator for ham_float array
// 1st argument: pointer for a ham_float array
// 2nd argument: size of array
ham_float variance(const ham_float *, const ham_uint &);
// variance estimator for ham_float vector
// 1st argument: ref for a ham_float vector
ham_float variance(const std::vector<ham_float> &);
// covariance estimator for ham_float arrays
// 1st argument: pointer for 1st ham_float array
// 2nd argument: pointer for 2nd ham_float array
// 3rd argument: size of array
ham_float covariance(const ham_float *, const ham_float *, const ham_uint &);
// covariance estimator ham_float vectors
// 1st argument: ref for 1st ham_float vector
// 2nd argument: ref for 2nd ham_float vector
ham_float covariance(const std::vector<ham_float> &,
                     const std::vector<ham_float> &);
// use given seed number or generate random seed according to thread and clock
ham_uint random_seed(const ham_int &);
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
ham_int fetchint(tinyxml2::XMLElement *, const std::string &,
                 const std::string &);
// get integer attribute
// 1st argument: ptr to tinyxml2::XMLElement (key)
// 2nd argument: attribute name
ham_int fetchint(tinyxml2::XMLElement *, const std::string &);
// get ham_uinteger attribute
// 1st argument: ptr to tinyxml2::XMLElement (key)
// 2nd argument: attribute name
// 3rd argument: sub-key under 1st argument
ham_uint fetchuint(tinyxml2::XMLElement *, const std::string &,
                   const std::string &);
// get ham_uinteger attribute
// 1st argument: ptr to tinyxml2::XMLElement (key)
// 2nd argument: attribute name
ham_uint fetchuint(tinyxml2::XMLElement *, const std::string &);
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
// get ham_float attribute
// 1st argument: ptr to tinyxml2::XMLElement (key)
// 2nd argument: attribute name
// 3rd argument: sub-key under 1st argument
ham_float fetchfloat(tinyxml2::XMLElement *, const std::string &,
                     const std::string &);
// get ham_float attribute
// 1st argument: ptr to tinyxml2::XMLElement (key)
// 2nd argument: attribute name
ham_float fetchfloat(tinyxml2::XMLElement *, const std::string &);
} // namespace toolkit

#endif
