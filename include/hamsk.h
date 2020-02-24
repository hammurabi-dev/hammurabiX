#ifndef HAMMURABI_MSK_H
#define HAMMRABI_MSK_H

#include <iostream>
#include <map>
#include <memory>
#include <omp.h>
#include <stdexcept>

#include <hamdis.h>
#include <hamtype.h>

template <typename T> class Hamsk {
private:
public:
  Hamsk() = default;
  Hamsk(const Hamsk<T> &) = delete;
  Hamsk &operator=(const Hamsk<T> &) = delete;
  Hamsk &operator=(Hamsk<T> &&) = delete;
  Hamsk(Hamsk<T> &&) = delete;
  virtual ~Hamsk() = default;
  virtual void duplicate(const ham_uint &) {
    throw std::runtime_error("wrong inheritance");
  }
  virtual T data(const ham_uint &, const ham_uint &) const {
    throw std::runtime_error("wrong inheritance");
  }
};

// for Hampix masking
template <typename T> class Hampisk final : public Hamsk<T> {
protected:
  std::unique_ptr<std::map<ham_uint, Hampix<T>>> Maskmaps;
  ham_uint Pivot_Nside = 0;

public:
  Hampisk() {
    this->Maskmaps = std::make_unique<std::map<ham_uint, Hampix<T>>>();
  }
  Hampisk(const Hampix<T> &m) {
    this->Pivot_Nside = m.nside();
    this->Maskmaps = std::make_unique<std::map<ham_uint, Hampix<T>>>();
    this->Maskmaps->insert({this->Pivot_Nside, m});
  }
  Hampisk(const Hampisk<T> &) = delete;
  Hampisk &operator=(const Hampisk<T> &) = delete;
  Hampisk &operator=(Hampisk<T> &&) = delete;
  Hampisk(Hampisk<T> &&) = delete;
  virtual ~Hampisk() = default;
  // check the pivot nside
  ham_uint pivot() const { return this->Pivot_Nside; }
  // make a copy at the new resolution
  // from the pivot(input)-resolution mask
  // 1st argument: target Nside
  void duplicate(const ham_uint &nside) override {
    auto tmp = std::make_unique<Hampix<ham_float>>(nside);
    if (this->Maskmaps->find(nside) != this->Maskmaps->end())
      return; // avoid existed copy
    if (this->Pivot_Nside > nside) {
      // prepare tmp map with "downgrading"
      const ham_uint fact{this->Pivot_Nside / nside};
      const ham_uint npix{12 * nside * nside};
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (ham_uint m = 0; m < npix; ++m) {
        ham_int x, y, f;
        tmp->rpixxyf(m, x, y, f);
        bool hit{false};
        for (ham_uint j = fact * y; j < fact * (y + 1); ++j) {
          for (ham_uint i = fact * x; i < fact * (x + 1); ++i) {
            ham_uint opix =
                this->Maskmaps->at(this->Pivot_Nside).xyfrpix(i, j, f);
            hit = (this->Maskmaps->at(this->Pivot_Nside).data(opix) > 0);
            if (hit)
              break;
          }
          if (hit)
            break;
        }
        tmp->data(m, static_cast<T>(hit));
      }
    } else if (this->Pivot_Nside < nside) {
      // prepare tmp map with "upgrading"
      const ham_uint fact{nside / this->Pivot_Nside};
      const ham_uint npix{12 * this->Pivot_Nside * this->Pivot_Nside};
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (ham_uint m = 0; m < npix; ++m) {
        const T val{this->Maskmaps->at(this->Pivot_Nside).data(m)};
        ham_int x, y, f;
        this->Maskmaps->at(this->Pivot_Nside).rpixxyf(m, x, y, f);
        for (ham_uint j = fact * y; j < fact * (y + 1); ++j) {
          for (ham_uint i = fact * x; i < fact * (x + 1); ++i) {
            ham_uint nupix{tmp->xyfrpix(i, j, f)};
            tmp->data(nupix, val);
          }
        }
      }
    }
    this->Maskmaps->insert({nside, *tmp});
  }
  // return the mask info at given resolution and index
  // 1st argument: target mask HEALPix Nside
  // 2nd argument: target map index
  inline T data(const ham_uint &nside, const ham_uint &idx) const {
    return (this->Maskmaps->at(nside)).data(idx);
  }
};

#endif
