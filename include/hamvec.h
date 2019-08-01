// hammurabi multi-dimensional (1D-3D) vector class based on std::vector

#ifndef HAMMURABI_VECTOR_H
#define HAMMURABI_VECTOR_H

#include <cassert>
#include <cmath>
#include <iostream>
#include <type_traits>
#include <vector>

template <int dim, typename T> class hamvec {
protected:
  std::vector<T> ele;

public:
  // default ctor
  // initialize with 0
  hamvec() {
    switch (dim) {
    case 1:
      this->ele = {T(0)};
      break;
    case 2:
      this->ele = {T(0), T(0)};
      break;
    case 3:
      this->ele = {T(0), T(0), T(0)};
      break;
    default:
      std::cerr << "unsupported dimension";
      break;
    }
  }
  virtual ~hamvec() = default;
  // 1D vector
  hamvec<dim, T>(const T &x) {
    assert(dim == 1);
    this->ele.push_back(x);
  }
  // 2D vector
  hamvec<dim, T>(const T &x, const T &y) {
    assert(dim == 2);
    this->ele.push_back(x);
    this->ele.push_back(y);
  }
  // 3D vector
  hamvec<dim, T>(const T &x, const T &y, const T &z) {
    assert(dim == 3);
    this->ele.push_back(x);
    this->ele.push_back(y);
    this->ele.push_back(z);
  }
  // copy ctor
  hamvec<dim, T>(const hamvec<dim, T> &v) { this->ele = v.content(); }
  // move ctor
  hamvec<dim, T>(hamvec<dim, T> &&v) : ele(std::move(v.content())) {}
  // copy assign
  hamvec<dim, T> &operator=(const hamvec<dim, T> &v) noexcept {
    this->ele = std::move(v.content());
    return *this;
  }
  // move assign
  hamvec<dim, T> &operator=(hamvec<dim, T> &&v) noexcept {
    this->ele = std::move(v.content());
    return *this;
  }
  // constant operator []
  T operator[](const int &i) const { return this->ele[i]; }
  // operator []
  T &operator[](const int &i) { return this->ele[i]; }
  // get constant std::vector<T> ele
  const std::vector<T> content() const { return this->ele; }
  // get std::vector<T> ele
  std::vector<T> &content() { return this->ele; }
  // operator +
  // cast argument to the same template type
  template <typename R>
  hamvec<dim, T> operator+(const hamvec<dim, R> &v) const {
    hamvec<dim, T> tmp(*this);
    for (unsigned int i = 0; i < dim; ++i) {
      tmp[i] += static_cast<T>(v[i]);
    }
    return tmp;
  }
  // operator +=
  // cast argument to the same template type
  template <typename R> hamvec<dim, T> &operator+=(const hamvec<dim, R> &v) {
    for (unsigned int i = 0; i < dim; ++i) {
      this->ele[i] += static_cast<T>(v[i]);
    }
    return *this;
  }
  // operator -
  // cast argument to the same template type
  template <typename R>
  hamvec<dim, T> operator-(const hamvec<dim, R> &v) const {
    hamvec<dim, T> tmp(*this);
    for (unsigned int i = 0; i < dim; ++i) {
      tmp[i] -= static_cast<T>(v[i]);
    }
    return tmp;
  }
  // operator -=
  // cast argument to the same template type
  template <typename R> hamvec<dim, T> &operator-=(const hamvec<dim, R> &v) {
    for (unsigned int i = 0; i < dim; ++i) {
      this->ele[i] -= static_cast<T>(v[i]);
    }
    return *this;
  }
  // operator *
  // cast argument to the same template type
  template <typename R> hamvec<dim, T> operator*(const R &s) const {
    hamvec<dim, T> tmp(*this);
    for (unsigned int i = 0; i < dim; ++i) {
      tmp[i] *= static_cast<T>(s);
    }
    return tmp;
  }
  // operaotr *=
  // cast argument to the same template type
  template <typename R> hamvec<dim, T> &operator*=(const R &s) {
    for (unsigned int i = 0; i < dim; ++i) {
      this->ele[i] *= static_cast<T>(s);
    }
    return *this;
  }
  // operator /
  // cast argument to the same template type
  template <typename R> hamvec<dim, T> operator/(const R &s) const {
    hamvec<dim, T> tmp(*this);
    for (unsigned int i = 0; i < dim; ++i) {
      tmp[i] /= static_cast<T>(s);
    }
    return tmp;
  }
  // operator /=
  // cast argument to the same template type
  template <typename R> hamvec<dim, T> &operator/=(const R &s) {
    assert(s != 0);
    for (unsigned int i = 0; i < dim; ++i) {
      this->ele[i] /= static_cast<T>(s);
    }
    return *this;
  }
  // operator ==
  bool operator==(const hamvec<dim, T> &v) {
    for (unsigned int i = 0; i < dim; ++i) {
      if (this->ele[i] != v[i]) {
        return false;
      }
    }
    return true;
  }
  // operator !=
  bool operator!=(const hamvec<dim, T> &v) {
    for (unsigned int i = 0; i < dim; ++i) {
      if (this->ele[i] != v[i]) {
        return true;
      }
    }
    return false;
  }
  // vector length
  // cast into double type
  double length() const {
    double tmp{0};
    for (unsigned int i = 0; i < dim; ++i) {
      const double cache = static_cast<double>(this->ele[i]);
      tmp += cache * cache;
    }
    return std::sqrt(tmp);
  }
  // vector squared length
  // cast into double type
  double lengthsq() const {
    double tmp{0};
    for (unsigned int i = 0; i < dim; ++i) {
      const double cache = static_cast<double>(this->ele[i]);
      tmp += cache * cache;
    }
    return tmp;
  }
  // flip sign
  // should not be used to unsigned type
  void flip() {
    assert(std::is_signed<T>::value);
    for (unsigned int i = 0; i < dim; ++i) {
      this->ele[i] *= static_cast<T>(-1.0);
    }
  }
  // versor
  // cast into double type
  hamvec<dim, double> versor() const {
    hamvec<dim, double> tmp;
    for (unsigned int i = 0; i < dim; ++i)
      tmp[i] = static_cast<double>(this->ele[i]);
    const auto l2{tmp.lengthsq()};
    if (l2 != 0)
      tmp /= std::sqrt(l2);
    return tmp;
  }
  // inner product
  // cast into double type
  template <typename R> double dotprod(const hamvec<dim, R> &v) const {
    double tmp{0};
    for (unsigned int i = 0; i < dim; ++i) {
      tmp += static_cast<double>(this->ele[i]) * static_cast<double>(v[i]);
    }
    return tmp;
  }
  // cross product, works in 3D only
  // cast into double type
  template <typename R>
  hamvec<dim, double> crossprod(const hamvec<dim, R> &v) const {
    assert(dim == 3);
    hamvec<dim, double> tmp;
    for (unsigned int i = 0; i < dim; ++i) {
      tmp[i] = (static_cast<double>(this->ele[(i + 1) % 3]) *
                    static_cast<double>(v[(i + 2) % 3]) -
                static_cast<double>(this->ele[(i + 2) % 3]) *
                    static_cast<double>(v[(i + 1) % 3]));
    }
    return tmp;
  }
  // osteam function
  friend std::ostream &operator<<(std::ostream &os, const hamvec<dim, T> &v) {
    os << dim << "D vector: ";
    for (unsigned int i = 0; i < dim; ++i) {
      os << v[i] << "\t";
    }
    os << std::endl;
    return os;
  }
};

#endif
