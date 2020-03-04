// hammurabi (3D) vector class
//
// dotprod, versor, length, lengthsq
// return floating datatype independent of inputs

#ifndef HAMMURABI_VECTOR_H
#define HAMMURABI_VECTOR_H

#include <array>
#include <cassert>
#include <cmath>
#include <hamtype.h>
#include <ostream>
#include <stdexcept>
#include <type_traits>

template <ham_sint dim, typename T> class Hamvec {
protected:
  std::array<T, 3> Ele = {T(0), T(0), T(0)};

public:
  // default ctor
  // initialize with 0
  Hamvec() = default;
  virtual ~Hamvec() = default;
  // 1D vector
  Hamvec<dim, T>(const T &x) {
    assert(dim == 1);
    this->Ele[0] = x;
  }
  // 2D vector
  Hamvec<dim, T>(const T &x, const T &y) {
    assert(dim == 2);
    this->Ele[0] = x;
    this->Ele[1] = y;
  }
  // 3D vector
  Hamvec<dim, T>(const T &x, const T &y, const T &z) {
    assert(dim == 3);
    this->Ele = {x, y, z};
  }
  // copy ctor
  Hamvec<dim, T>(const Hamvec<dim, T> &v) { this->Ele = v.Ele; }
  // move ctor
  Hamvec<dim, T>(Hamvec<dim, T> &&v) { this->Ele = std::move(v.Ele); }
  // copy assign
  Hamvec<dim, T> &operator=(const Hamvec<dim, T> &v) noexcept {
    this->Ele = std::move(v.Ele);
    return *this;
  }
  // move assign
  Hamvec<dim, T> &operator=(Hamvec<dim, T> &&v) noexcept {
    this->Ele = std::move(v.Ele);
    return *this;
  }
  // constant operator []
  T operator[](const ham_sint &i) const { return this->Ele[i]; }
  // operator []
  T &operator[](const ham_sint &i) { return this->Ele[i]; }
  // operator +
  // cast argument to the same template type
  Hamvec<dim, T> operator+(const Hamvec<dim, T> &v) const {
    Hamvec<dim, T> tmp(*this);
    for (ham_sint i = 0; i < dim; ++i) {
      tmp[i] += v[i];
    }
    return tmp;
  }
  // operator +=
  // cast argument to the same template type
  Hamvec<dim, T> &operator+=(const Hamvec<dim, T> &v) {
    for (ham_sint i = 0; i < dim; ++i) {
      this->Ele[i] += v[i];
    }
    return *this;
  }
  // operator -
  // cast argument to the same template type
  Hamvec<dim, T> operator-(const Hamvec<dim, T> &v) const {
    Hamvec<dim, T> tmp(*this);
    for (ham_sint i = 0; i < dim; ++i) {
      tmp[i] -= v[i];
    }
    return tmp;
  }
  // operator -=
  // cast argument to the same template type
  Hamvec<dim, T> &operator-=(const Hamvec<dim, T> &v) {
    for (ham_sint i = 0; i < dim; ++i) {
      this->Ele[i] -= v[i];
    }
    return *this;
  }
  // operator *
  // cast argument to the same template type
  Hamvec<dim, T> operator*(const T &s) const {
    Hamvec<dim, T> tmp(*this);
    for (ham_sint i = 0; i < dim; ++i) {
      tmp[i] *= s;
    }
    return tmp;
  }
  // operaotr *=
  // cast argument to the same template type
  Hamvec<dim, T> &operator*=(const T &s) {
    for (ham_sint i = 0; i < dim; ++i) {
      this->Ele[i] *= s;
    }
    return *this;
  }
  // operator /
  // cast argument to the same template type
  Hamvec<dim, T> operator/(const T &s) const {
    assert(s != 0);
    Hamvec<dim, T> tmp(*this);
    for (ham_sint i = 0; i < dim; ++i) {
      tmp[i] /= s;
    }
    return tmp;
  }
  // operator /=
  // cast argument to the same template type
  Hamvec<dim, T> &operator/=(const T &s) {
    assert(s != 0);
    for (ham_sint i = 0; i < dim; ++i) {
      this->Ele[i] /= s;
    }
    return *this;
  }
  // operator ==
  bool operator==(const Hamvec<dim, T> &v) {
    for (ham_sint i = 0; i < dim; ++i) {
      if (this->Ele[i] != v[i]) {
        return false;
      }
    }
    return true;
  }
  // operator !=
  bool operator!=(const Hamvec<dim, T> &v) {
    for (ham_sint i = 0; i < dim; ++i) {
      if (this->Ele[i] != v[i]) {
        return true;
      }
    }
    return false;
  }
  // vector length
  // cast into ham_float type
  inline ham_float length() const {
    ham_float tmp{T(0)};
    for (ham_sint i = 0; i < dim; ++i) {
      tmp += static_cast<ham_float>(this->Ele[i] * this->Ele[i]);
    }
    return std::sqrt(tmp);
  }
  // vector squared length
  // cast into ham_float type
  inline ham_float lengthsq() const {
    T tmp{T(0)};
    for (ham_sint i = 0; i < dim; ++i) {
      tmp += this->Ele[i] * this->Ele[i];
    }
    return static_cast<ham_float>(tmp);
  }
  // flip sign
  // should not be used to unsigned type
  inline void flip() {
    assert(std::is_signed<T>::value);
    for (ham_sint i = 0; i < dim; ++i) {
      this->Ele[i] = -this->Ele[i];
    }
  }
  // versor
  // cast into ham_float type
  inline Hamvec<dim, ham_float> versor() const {
    const auto l2{this->lengthsq()};
    Hamvec<dim, ham_float> tmp;
    if (l2 != 0) {
      const ham_float denominator = 1. / std::sqrt(l2);
      for (ham_sint i = 0; i < dim; ++i) {
        tmp[i] = this->Ele[i] * denominator;
      }
    }
    return tmp;
  }
  // inner product
  // cast into ham_float type
  inline ham_float dotprod(const Hamvec<dim, T> &v) const {
    ham_float tmp{0};
    for (ham_sint i = 0; i < dim; ++i) {
      tmp += this->Ele[i] * v[i];
    }
    return tmp;
  }
  // cross product, works in 3D only
  // cast into ham_float type
  Hamvec<dim, T> crossprod(const Hamvec<dim, T> &v) const {
    assert(dim == 3);
    Hamvec<dim, ham_float> tmp;
    for (ham_sint i = 0; i < dim; ++i) {
      tmp[i] = (this->Ele[(i + 1) % 3] * v[(i + 2) % 3] -
                this->Ele[(i + 2) % 3] * v[(i + 1) % 3]);
    }
    return tmp;
  }
  // osteam function
  friend std::ostream &operator<<(std::ostream &os, const Hamvec<dim, T> &v) {
    os << dim << "D vector: ";
    for (ham_sint i = 0; i < dim; ++i) {
      os << v[i] << "\t";
    }
    os << std::endl;
    return os;
  }
};

#endif
