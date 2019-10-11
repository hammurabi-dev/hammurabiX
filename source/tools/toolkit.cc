#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <omp.h>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include <tinyxml2.h>
#include <toolkit.h>

namespace toolkit {
// mean for double array
double mean(const double *arr, const std::size_t &size) {
  double avg{0};
  for (std::size_t i = 0; i != size; ++i) {
    avg += arr[i];
  }
  avg /= size;
  return avg;
}
// mean for double vector
double mean(const std::vector<double> &vect) {
  assert(!vect.empty());
  double avg{0};
  for (auto &i : vect) {
    avg += i;
  }
  avg /= vect.size();
  return avg;
}
// variance for double array
double variance(const double *arr, const std::size_t &size) {
  const double avg{mean(arr, size)};
  double var{0.};
  for (std::size_t i = 0; i != size; ++i) {
    var += (arr[i] - avg) * (arr[i] - avg);
  }
  var /= size;
  return var;
}
// variance for double vector
double variance(const std::vector<double> &vect) {
  assert(!vect.empty());
  const double avg{mean(vect)};
  double var{0.};
  for (auto &i : vect) {
    var += (i - avg) * (i - avg);
  }
  var /= vect.size();
  return var;
}
// covariance for double array
double covariance(const double *arr1, const double *arr2,
                  const std::size_t &size) {
  double avg1{mean(arr1, size)};
  double avg2{mean(arr2, size)};
  double covar{0.};
  for (std::size_t m = 0; m != size; ++m) {
    covar += (arr1[m] - avg1) * (arr2[m] - avg2);
  }
  covar /= size;
  return covar;
}
// covariance for double vector
double covariance(const std::vector<double> &vect1,
                  const std::vector<double> &vect2) {
  assert(vect1.size() == vect2.size());
  double avg1{mean(vect1)};
  double avg2{mean(vect2)};
  double covar{0.};
  for (decltype(vect1.size()) i = 0; i != vect1.size(); ++i) {
    covar += (vect1[i] - avg1) * (vect2[i] - avg2);
  }
  covar /= vect1.size();
  return covar;
}
// random seed generator
std::size_t random_seed(const int &s) {
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
int fetchint(tinyxml2::XMLElement *el, const std::string &att_type,
             const std::string &key) {
#ifdef VERBOSE
  std::cout << "key: " << key << "attrib: " << att_type << std::endl;
#endif
  return el->FirstChildElement(key.c_str())->IntAttribute(att_type.c_str());
}
// get integer value of attribute in current level, with default return
int fetchint(tinyxml2::XMLElement *el, const std::string &att_type) {
#ifdef VERBOSE
  std::cout << "attrib: " << att_type << std::endl;
#endif
  return el->IntAttribute(att_type.c_str());
}
// get unsigned value of attribute in child level, with default return
unsigned int fetchunsigned(tinyxml2::XMLElement *el,
                           const std::string &att_type,
                           const std::string &key) {
#ifdef VERBOSE
  std::cout << "key: " << key << "attrib: " << att_type << std::endl;
#endif
  return el->FirstChildElement(key.c_str())
      ->UnsignedAttribute(att_type.c_str());
}
// get unsigned value of attribute in current level, with default return
unsigned int fetchunsigned(tinyxml2::XMLElement *el,
                           const std::string &att_type) {
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
// get double value of attribute in child level, with default return
double fetchdouble(tinyxml2::XMLElement *el, const std::string &att_type,
                   const std::string &key) {
#ifdef VERBOSE
  std::cout << "key: " << key << "attrib: " << att_type << std::endl;
#endif
  return el->FirstChildElement(key.c_str())->DoubleAttribute(att_type.c_str());
}
// get double value of attribute in current level, with default return
double fetchdouble(tinyxml2::XMLElement *el, const std::string &att_type) {
#ifdef VERBOSE
  std::cout << "attrib: " << att_type << std::endl;
#endif
  return el->DoubleAttribute(att_type.c_str());
}

} // namespace toolkit
