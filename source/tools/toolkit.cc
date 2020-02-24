#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <omp.h>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include <hamtype.h>
#include <tinyxml2.h>
#include <toolkit.h>

namespace toolkit {
/*
// Kahan summation algorithm
class ksum;
*/
// mean for ham_float array
ham_float mean(const ham_float *arr, const ham_uint &size) {
  ham_float avg{0};
  for (ham_uint i = 0; i != size; ++i) {
    avg += arr[i];
  }
  avg /= size;
  return avg;
}
// mean for ham_float vector
ham_float mean(const std::vector<ham_float> &vect) {
  assert(!vect.empty());
  ham_float avg{0};
  for (auto &i : vect) {
    avg += i;
  }
  avg /= vect.size();
  return avg;
}
// variance for ham_float array
ham_float variance(const ham_float *arr, const ham_uint &size) {
  const ham_float avg{mean(arr, size)};
  ham_float var{0.};
  for (ham_uint i = 0; i != size; ++i) {
    var += (arr[i] - avg) * (arr[i] - avg);
  }
  var /= size;
  return var;
}
// variance for ham_float vector
ham_float variance(const std::vector<ham_float> &vect) {
  assert(!vect.empty());
  const ham_float avg{mean(vect)};
  ham_float var{0.};
  for (auto &i : vect) {
    var += (i - avg) * (i - avg);
  }
  var /= vect.size();
  return var;
}
// covariance for ham_float array
ham_float covariance(const ham_float *arr1, const ham_float *arr2,
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
// covariance for ham_float vector
ham_float covariance(const std::vector<ham_float> &vect1,
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
// random seed generator
ham_uint random_seed(const ham_int &s) {
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
// load xml file
std::unique_ptr<tinyxml2::XMLDocument> loadxml(const std::string &filename) {
  auto doc = std::make_unique<tinyxml2::XMLDocument>();
  doc->LoadFile(filename.c_str());
  assert(!doc->Error());
  return std::move(doc);
}
// trace xml element
tinyxml2::XMLElement *tracexml(tinyxml2::XMLDocument *doc,
                               const std::vector<std::string> &keychain) {
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
// get string value of attribute in child level
std::string fetchstring(tinyxml2::XMLElement *el, const std::string &att_type,
                        const std::string &key) {
#ifdef VERBOSE
  std::cout << "key: " << key << "attrib: " << att_type << std::endl;
#endif
  return el->FirstChildElement(key.c_str())->Attribute(att_type.c_str());
}
// get string value of attribute in current level
std::string fetchstring(tinyxml2::XMLElement *el, const std::string &att_type) {
#ifdef VERBOSE
  std::cout << "attrib: " << att_type << std::endl;
#endif
  return el->Attribute(att_type.c_str());
}
// get integer value of attribute in child level, with default return
ham_int fetchint(tinyxml2::XMLElement *el, const std::string &att_type,
                 const std::string &key) {
#ifdef VERBOSE
  std::cout << "key: " << key << "attrib: " << att_type << std::endl;
#endif
  return el->FirstChildElement(key.c_str())->IntAttribute(att_type.c_str());
}
// get integer value of attribute in current level, with default return
ham_int fetchint(tinyxml2::XMLElement *el, const std::string &att_type) {
#ifdef VERBOSE
  std::cout << "attrib: " << att_type << std::endl;
#endif
  return el->IntAttribute(att_type.c_str());
}
// get unsigned value of attribute in child level, with default return
ham_uint fetchuint(tinyxml2::XMLElement *el, const std::string &att_type,
                   const std::string &key) {
#ifdef VERBOSE
  std::cout << "key: " << key << "attrib: " << att_type << std::endl;
#endif
  return el->FirstChildElement(key.c_str())
      ->UnsignedAttribute(att_type.c_str());
}
// get unsigned value of attribute in current level, with default return
ham_uint fetchuint(tinyxml2::XMLElement *el, const std::string &att_type) {
#ifdef VERBOSE
  std::cout << "attrib: " << att_type << std::endl;
#endif
  return el->UnsignedAttribute(att_type.c_str());
}
// get Boolean value of attribute in child level, with default return
bool fetchbool(tinyxml2::XMLElement *el, const std::string &att_type,
               const std::string &key) {
#ifdef VERBOSE
  std::cout << "key: " << key << "attrib: " << att_type << std::endl;
#endif
  return el->FirstChildElement(key.c_str())->BoolAttribute(att_type.c_str());
}
// get Boolean value of attribute in current level, with default return
bool fetchbool(tinyxml2::XMLElement *el, const std::string &att_type) {
#ifdef VERBOSE
  std::cout << "attrib: " << att_type << std::endl;
#endif
  return el->BoolAttribute(att_type.c_str());
}
// get ham_float value of attribute in child level, with default return
ham_float fetchfloat(tinyxml2::XMLElement *el, const std::string &att_type,
                     const std::string &key) {
#ifdef VERBOSE
  std::cout << "key: " << key << "attrib: " << att_type << std::endl;
#endif
  return el->FirstChildElement(key.c_str())->DoubleAttribute(att_type.c_str());
}
// get ham_float value of attribute in current level, with default return
ham_float fetchfloat(tinyxml2::XMLElement *el, const std::string &att_type) {
#ifdef VERBOSE
  std::cout << "attrib: " << att_type << std::endl;
#endif
  return el->DoubleAttribute(att_type.c_str());
}

} // namespace toolkit

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
