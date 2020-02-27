// auxiliary functions

#ifndef HAMMURABI_TOOLKIT_H
#define HAMMURABI_TOOLKIT_H

#include <cassert>
#include <cmath>
#include <ctime>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <thread>
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
inline ham_float mean(const ham_float *arr, const ham_uint &size) {
  ham_float avg{0};
  for (ham_uint i = 0; i != size; ++i) {
    avg += arr[i];
  }
  avg /= size;
  return avg;
}
// mean value estimator for ham_float vector
// 1st argument: ref for a ham_float vector
inline ham_float mean(const std::vector<ham_float> &vect) {
  assert(!vect.empty());
  ham_float avg{0};
  for (auto &i : vect) {
    avg += i;
  }
  avg /= vect.size();
  return avg;
}
// variance estimator for ham_float array
// 1st argument: pointer for a ham_float array
// 2nd argument: size of array
inline ham_float variance(const ham_float *arr, const ham_uint &size) {
  const ham_float avg{mean(arr, size)};
  ham_float var{0.};
  for (ham_uint i = 0; i != size; ++i) {
    var += (arr[i] - avg) * (arr[i] - avg);
  }
  var /= size;
  return var;
}
// variance estimator for ham_float vector
// 1st argument: ref for a ham_float vector
inline ham_float variance(const std::vector<ham_float> &vect) {
  assert(!vect.empty());
  const ham_float avg{mean(vect)};
  ham_float var{0.};
  for (auto &i : vect) {
    var += (i - avg) * (i - avg);
  }
  var /= vect.size();
  return var;
}
// covariance estimator for ham_float arrays
// 1st argument: pointer for 1st ham_float array
// 2nd argument: pointer for 2nd ham_float array
// 3rd argument: size of array
inline ham_float covariance(const ham_float *arr1, const ham_float *arr2,
                            const ham_uint &size) {
  ham_float avg1{mean(arr1, size)};
  ham_float avg2{mean(arr2, size)};
  ham_float covar{0.};
  for (ham_uint m = 0; m != size; ++m) {
    covar += (arr1[m] - avg1) * (arr2[m] - avg2);
  }
  covar /= size;
  return covar;
}
// covariance estimator ham_float vectors
// 1st argument: ref for 1st ham_float vector
// 2nd argument: ref for 2nd ham_float vector
inline ham_float covariance(const std::vector<ham_float> &vect1,
                            const std::vector<ham_float> &vect2) {
  assert(vect1.size() == vect2.size());
  ham_float avg1{mean(vect1)};
  ham_float avg2{mean(vect2)};
  ham_float covar{0.};
  for (decltype(vect1.size()) i = 0; i != vect1.size(); ++i) {
    covar += (vect1[i] - avg1) * (vect2[i] - avg2);
  }
  covar /= vect1.size();
  return covar;
}
// use given seed number or generate random seed according to thread and clock
inline ham_uint random_seed(const ham_int &s) {
  assert(s >= 0);
  if (s == 0) {
    auto p = std::chrono::system_clock::now();
    // valid until 19 January, 2038 03:14:08 UTC
    time_t today_time = std::chrono::system_clock::to_time_t(p);
    // casting thread id into unsinged long
    std::stringstream ss;
    ss << std::this_thread::get_id();
    auto th_id = std::stoul(ss.str());
    // precision in (thread,second)
    return (th_id + today_time);
  }
  return s;
}
// load tinyxml2::XML file
// 1st argument: tinyxml2::XML file name (with dir)
inline std::unique_ptr<tinyxml2::XMLDocument>
loadxml(const std::string &filename) {
  auto doc = std::make_unique<tinyxml2::XMLDocument>();
  doc->LoadFile(filename.c_str());
  assert(!doc->Error());
  return doc;
}
// trace down a key inside tinyxml2::XML "document"
// 1st argument: pointer to tinyxml2::XMLDocument
// 2nd argument: a vector of string, with key chain for tracing
// the last string in 2nd argument is the target key
// <root> is automatically included
inline tinyxml2::XMLElement *
tracexml(tinyxml2::XMLDocument *doc, const std::vector<std::string> &keychain) {
  tinyxml2::XMLElement *el{doc->FirstChildElement("root")};
  if (!keychain.empty()) {
    for (auto key : keychain) {
#ifdef VERBOSE
      std::cout << "key: " << key << std::endl;
#endif
      el = el->FirstChildElement(key.c_str());
    }
  }
  return el;
}
// get string attribute
// 1st argument: ptr to tinyxml2::XMLElement (key)
// 2nd argument: attribute name
// 3rd argument: sub-key under 1st argument
inline std::string fetchstring(tinyxml2::XMLElement *el,
                               const std::string &att_type,
                               const std::string &key) {
#ifdef VERBOSE
  std::cout << "key: " << key << "attrib: " << att_type << std::endl;
#endif
  return el->FirstChildElement(key.c_str())->Attribute(att_type.c_str());
}
// get string attribute
// 1st argument: ptr to tinyxml2::XMLElement (key)
// 2nd argument: attribute name
inline std::string fetchstring(tinyxml2::XMLElement *el,
                               const std::string &att_type) {
#ifdef VERBOSE
  std::cout << "attrib: " << att_type << std::endl;
#endif
  return el->Attribute(att_type.c_str());
}
// get integer attribute
// 1st argument: ptr to tinyxml2::XMLElement (key)
// 2nd argument: attribute name
// 3rd argument: sub-key under 1st argument
inline ham_int fetchint(tinyxml2::XMLElement *el, const std::string &att_type,
                        const std::string &key) {
#ifdef VERBOSE
  std::cout << "key: " << key << "attrib: " << att_type << std::endl;
#endif
  return el->FirstChildElement(key.c_str())->IntAttribute(att_type.c_str());
}
// get integer attribute
// 1st argument: ptr to tinyxml2::XMLElement (key)
// 2nd argument: attribute name
inline ham_int fetchint(tinyxml2::XMLElement *el, const std::string &att_type) {
#ifdef VERBOSE
  std::cout << "attrib: " << att_type << std::endl;
#endif
  return el->IntAttribute(att_type.c_str());
}
// get ham_uinteger attribute
// 1st argument: ptr to tinyxml2::XMLElement (key)
// 2nd argument: attribute name
// 3rd argument: sub-key under 1st argument
inline ham_uint fetchuint(tinyxml2::XMLElement *el, const std::string &att_type,
                          const std::string &key) {
#ifdef VERBOSE
  std::cout << "key: " << key << "attrib: " << att_type << std::endl;
#endif
  return el->FirstChildElement(key.c_str())
      ->UnsignedAttribute(att_type.c_str());
}
// get ham_uinteger attribute
// 1st argument: ptr to tinyxml2::XMLElement (key)
// 2nd argument: attribute name
inline ham_uint fetchuint(tinyxml2::XMLElement *el,
                          const std::string &att_type) {
#ifdef VERBOSE
  std::cout << "attrib: " << att_type << std::endl;
#endif
  return el->UnsignedAttribute(att_type.c_str());
}
// get bool attribute
// 1st argument: ptr to tinyxml2::XMLElement (key)
// 2nd argument: attribute name
// 3rd argument: sub-key under 1st argument
// 4th argument: default bool, use int instead of bool type
// since bool can take char* as input
inline bool fetchbool(tinyxml2::XMLElement *el, const std::string &att_type,
                      const std::string &key) {
#ifdef VERBOSE
  std::cout << "key: " << key << "attrib: " << att_type << std::endl;
#endif
  return el->FirstChildElement(key.c_str())->BoolAttribute(att_type.c_str());
}
// get bool attribute
// 1st argument: ptr to tinyxml2::XMLElement (key)
// 2nd argument: attribute name
inline bool fetchbool(tinyxml2::XMLElement *el, const std::string &att_type) {
#ifdef VERBOSE
  std::cout << "attrib: " << att_type << std::endl;
#endif
  return el->BoolAttribute(att_type.c_str());
}
// get ham_float attribute
// 1st argument: ptr to tinyxml2::XMLElement (key)
// 2nd argument: attribute name
// 3rd argument: sub-key under 1st argument
inline ham_float fetchfloat(tinyxml2::XMLElement *el,
                            const std::string &att_type,
                            const std::string &key) {
#ifdef VERBOSE
  std::cout << "key: " << key << "attrib: " << att_type << std::endl;
#endif
  return el->FirstChildElement(key.c_str())->DoubleAttribute(att_type.c_str());
}
// get ham_float attribute
// 1st argument: ptr to tinyxml2::XMLElement (key)
// 2nd argument: attribute name
inline ham_float fetchfloat(tinyxml2::XMLElement *el,
                            const std::string &att_type) {
#ifdef VERBOSE
  std::cout << "attrib: " << att_type << std::endl;
#endif
  return el->DoubleAttribute(att_type.c_str());
}
} // namespace toolkit

#endif

/*
 class toolkit::ksum{
 ham_float acc, c; // acc, the accumulator; c, the compensator
 public:
 ksum() : acc(0), c(0) {}
 virtual ~ksum() = default;
 void add (const ham_float &val) {
 const volatile ham_float y = val - c;
 const volatile ham_float t = acc + y;
 c = (t - acc) - y; // the missing compensator
 acc = t;
 }
 ham_float result() const {return acc;}
 };
 */
