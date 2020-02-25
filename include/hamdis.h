#ifndef HAMMURABI_DIS_H
#define HAMMURABI_DIS_H

#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <omp.h>
#include <stdexcept>
#include <vector>

#include <hamp.h>
#include <hamtype.h>
#include <hamunits.h>

// data type T
// key structure of each sample point
template <typename T> class Node {
protected:
  // pointing
  Hamp Pointing{0.0, 0.0};
  // pixel data
  T Data{static_cast<T>(0)};
  // pixel index
  ham_uint Index{0};

public:
  // constructor
  Node() = default;
  // direct constr
  Node(const Hamp &ptr, const T &val, const ham_uint &idx = 0) {
    this->Pointing = ptr;
    this->Data = val;
    this->Index = idx;
  }
  // copy assign
  Node &operator=(const Node<T> &n) {
    this->Pointing = n.Pointing;
    this->Data = n.Data;
    this->Index = n.Index;
    return *this;
  }
  // move assign
  Node &operator=(Node<T> &&n) noexcept {
    this->Pointing = std::move(n.Pointing);
    this->Data = std::move(n.Data);
    this->Index = std::move(n.Index);
    return *this;
  }
  // copy constr
  Node(const Node<T> &n) {
    this->Pointing = n.Pointing;
    this->Data = n.Data;
    this->Index = n.Index;
  }
  // move constr
  Node(Node<T> &&n) noexcept {
    this->Pointing = std::move(n.Pointing);
    this->Data = std::move(n.Data);
    this->Index = std::move(n.Index);
  }
  // destr
  virtual ~Node() = default;
  // extract sky position
  virtual Hamp pointing() const { return this->Pointing; }
  // extract sky information
  virtual T data() const { return this->Data; }
  // extract sky index
  virtual ham_uint index() const { return this->Index; }
  // update sky position
  virtual void pointing(const Hamp &new_pointing) {
    this->Pointing = new_pointing;
  }
  // update sky position
  virtual void pointing(const ham_float &theta, const ham_float &phi) {
    this->Pointing.theta(theta);
    this->Pointing.phi(phi);
  }
  // update sky information
  virtual void data(const T &new_data) { this->Data = new_data; }
  // update ksy index
  virtual void index(const ham_uint &new_idx) { this->Index = new_idx; }
};

// data type T
// hosts a vector of Node objects
// the indices of Nodes are trivial
template <typename T> class Hamdis {
protected:
  // map holder, unique pointer of vector of Nodes
  std::unique_ptr<std::vector<Node<T>>> Map;

public:
  // HEALPix underfined value for masking
  const ham_float undef{-1.6375e30};
  // dft constr
  Hamdis() { this->Map = std::make_unique<std::vector<Node<T>>>(); }
  // initialize map with given sample number N
  // Node's Data assigned by the given value
  // Node's Index assigned from 0 to N-1
  Hamdis(const ham_uint &N, const T &v = static_cast<T>(0)) {
    this->Map = std::make_unique<std::vector<Node<T>>>(N);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (ham_uint i = 0; i < N; ++i) {
      this->Map->at(i).index(i);
      this->Map->at(i).data(v);
    }
  }
  // copy constr
  Hamdis(const Hamdis<T> &m) {
    this->Map.reset(new std::vector<Node<T>>(*(m.Map.get())));
  }
  // move constr
  Hamdis(Hamdis<T> &&m) noexcept { this->Map = std::move(m.Map); }
  // move assign
  Hamdis &operator=(Hamdis<T> &&m) noexcept {
    this->Map = std::move(m.Map);
    return *this;
  }
  // copy assignment
  Hamdis &operator=(const Hamdis<T> &m) {
    this->Map.reset(new std::vector<Node<T>>(*(m.Map.get())));
    return *this;
  }
  // dft destr
  virtual ~Hamdis() = default;
  // extract node data
  virtual T data(const ham_uint &idx) const {
    return this->Map->at(idx).data();
  }
  // set node data
  virtual void data(const ham_uint &idx, const T &v) {
    this->Map->at(idx).data(v);
  }
  // extract node index
  virtual ham_uint index(const ham_uint &idx) const {
    return this->Map->at(idx).index();
  }
  // set node index
  virtual void index(const ham_uint &idx, const ham_uint &new_idx) {
    this->Map->at(idx).index(new_idx);
  }
  // extract node pointing
  virtual Hamp pointing(const ham_uint &idx) const {
    return this->Map->at(idx).pointing();
  }
  // set node pointing
  virtual void pointing(const ham_uint &idx, const Hamp &new_point) {
    this->Map->at(idx).pointing(new_point);
  }
  // set node pointing
  virtual void pointing(const ham_uint &idx, const ham_float &theta,
                        const ham_float &phi) {
    this->Map->at(idx).pointing(theta, phi);
  }
  // extract map size
  virtual ham_uint npix() const { return this->Map->size(); }
  // print to content of each pix to screen
  virtual void print() const {
    std::cout << "... printing Hamdis map information ..." << std::endl;
    // no need for multi-threading here
    for (auto i = this->Map->begin(); i < this->Map->end(); ++i) {
      std::cout << "index: " << i->index() << "\t"
                << "data: " << i->data() << "\t"
                << "pointing: theta " << i->pointing().theta() << " phi "
                << i->pointing().phi() << std::endl;
    }
  }
  // reset with given size and clean up data
  virtual void reset(const ham_uint &n = 0) {
    // cleaning an used map with correct size
    if (n == 0 or this->Map->size() == n) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (ham_uint i = 0; i < this->Map->size(); ++i) {
        this->Map->at(i).data(0.0);
      }
    } else {
      // there 2 cases when Nside != n
      // 1, map to be recycled with wrong size
      // 2, empty map initialized by the default constr
      this->Map = std::make_unique<std::vector<Node<T>>>((n));
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (ham_uint i = 0; i < n; ++i) {
        this->Map->at(i).index(i);
        this->Map->at(i).pointing(Hamp(0.0, 0.0));
        this->Map->at(i).data(0.0);
      }
    }
  }
  // rescale
  virtual void rescale(const ham_float &v) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (ham_uint i = 0; i < this->Map->size(); ++i) {
      T current_data = this->Map->at(i).data();
      this->Map->at(i).data(current_data * v);
    }
  }
  // undefine a certain Node
  virtual void undefine(const ham_uint &idx) {
    this->Map->at(idx).data(this->undef);
  }
  // undefine a list of Nodes
  virtual void undefine(const std::vector<ham_uint> &list) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (ham_uint i = 0; i < list.size(); ++i) {
      this->Map->at(list[i]).data(this->undef);
    }
  }
};

// using HEALPix RING ordering, with HEALPix Nside as power of 2
// also having trivial indexing
template <typename T> class Hampix final : public Hamdis<T> {
protected:
  // HEALPix nside_
  ham_uint Nside = 0;
  // HEALPix order_
  int Order = -1;
  // HEALPix npix_
  ham_uint Npix = 0;
  // HEALPix npface_
  ham_uint Npface = 0;
  // HEALPix ncap_
  ham_uint Ncap = 0;
  // HEALPix fact1_, fact2_
  ham_float Fact1 = 0.0;
  ham_float Fact2 = 0.0;
  std::array<ham_int, 12> Jrll = {2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4};
  std::array<ham_int, 12> Jpll = {1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7};
  // initialize constants
  // copied from healpix_base.cc SetNside function
  void prepare() {
    // order_
    this->Order = this->ilog2(this->Nside);
    // check order_
    if (this->Order < 0)
      throw std::runtime_error("unsupported HEALPix Nside");
    // npface_
    this->Npface = this->Nside * this->Nside;
    // npix_
    this->Npix = 12 * this->Npface;
    // ncap_
    this->Ncap = (this->Npface - this->Nside) << 1;
    // fact2
    this->Fact2 = 4. / this->Npix;
    // fact1
    this->Fact1 = (this->Nside << 1) * this->Fact2;
  }
  // copied from HEALPix cxxutils.h ilog2 function
  // Returns the largest integer n that fulfills 2^n<=arg.
  inline ham_uint ilog2(const ham_uint &n) {
    ham_uint tmp{n};
    ham_uint res{0};
    while (tmp > 0x0000FFFF) {
      res += 16;
      tmp >>= 16;
    }
    if (tmp > 0x000000FF) {
      res |= 8;
      tmp >>= 8;
    }
    if (tmp > 0x0000000F) {
      res |= 4;
      tmp >>= 4;
    }
    if (tmp > 0x00000003) {
      res |= 2;
      tmp >>= 2;
    }
    if (tmp > 0x00000001) {
      res |= 1;
    }
    return res;
  }
  // copied from HEALPix cxxutils.h isqrt function
  // Returns the integer n, which fulfills n*n<=arg<(n+1)*(n+1).
  inline ham_uint isqrt(const ham_uint &n) const {
    return static_cast<ham_uint>(
        std::floor(std::sqrt(static_cast<ham_float>(n) + 0.5)));
  }
  // copied from HEALPix healpix_base.h special_div function
  // low-level hack to accelerate divisions with a result of [0;3]
  inline ham_int special_div(const ham_int &a, const ham_int &b) const {
    ham_int p{a};
    ham_int q{b};
    ham_int t{(p >= (q << 1))};
    p -= t * (q << 1);
    return (t << 1) + (p >= q);
  }
  // assigning correct pointing position
  // using HEALPix implementation, equivalent to HEALPix pix2ang function
  Hamp fillpoint(const ham_uint &idx) {
    ham_float z{0.0};
    ham_float phi{0.0};
    ham_float sth{0.0};
    bool have_sth{false};
    // healpix RING ordering
    // copied from HEALPix healpix_base.cc pix2loc function
    // counted from North pole
    // North Polar cap
    if (idx < this->Ncap) {
      const ham_uint iring{
          (1 + static_cast<ham_uint>(this->isqrt(1 + 2 * idx))) >> 1};
      const ham_uint iphi{(idx + 1) - 2 * iring * (iring - 1)};
      const ham_float tmp{(iring * iring) * this->Fact2};
      z = 1.0 - tmp;
      if (z > 0.99) {
        sth = std::sqrt(tmp * (2.0 - tmp));
        have_sth = true;
      }
      phi = (iphi - 0.5) * cgs::halfpi / iring;
    }
    // Equatorial region
    else if (idx < (this->Npix - this->Ncap)) {
      const ham_uint nl4{4 * this->Nside};
      const ham_uint ip{idx - this->Ncap};
      const ham_uint tmp{(this->Order >= 0) ? ip >> (this->Order + 2)
                                            : ip / nl4};
      const ham_uint iring{tmp + this->Nside};
      const ham_uint iphi{ip - nl4 * tmp + 1};
      // 1 if iring+nside is odd, 1/2 otherwise
      const ham_float fodd{((iring + this->Nside) & 1) ? 1 : 0.5};
      // nside can be smaller than iring
      z = (static_cast<ham_float>(2 * this->Nside) -
           static_cast<ham_float>(iring)) *
          this->Fact1;
      phi = (iphi - fodd) * cgs::halfpi * 1.5 * this->Fact1;
    }
    // South Polar cap
    else {
      const ham_uint ip{this->Npix - idx};
      // counted from South pole
      const ham_uint iring{
          (1 + static_cast<ham_uint>(this->isqrt(2 * ip - 1))) >> 1};
      const ham_uint iphi{4 * iring + 1 - (ip - 2 * iring * (iring - 1))};
      const ham_float tmp{(iring * iring) * this->Fact2};
      z = tmp - 1.0;
      if (z < -0.99) {
        sth = std::sqrt(tmp * (2.0 - tmp));
        have_sth = true;
      }
      phi = (iphi - 0.5) * cgs::halfpi / iring;
    }
    // copied from healpix_base.h pix2ang function
    return have_sth ? Hamp(std::atan2(sth, z), phi) : Hamp(std::acos(z), phi);
  }
  // copied from HEALPix ring_above function
  ham_uint rabove(const ham_float &z) const {
    ham_float az{std::fabs(z)};
    if (az <= cgs::twothirds) // equatorial region
      return static_cast<ham_uint>(this->Nside * (2 - 1.5 * z));
    const ham_uint iring{
        static_cast<ham_uint>(this->Nside * std::sqrt(3 * (1 - az)))};
    return (z > 0) ? iring : 4 * this->Nside - iring - 1;
  }
  // copied from HEALPix get_ring_info2 function
  void rinfo(const ham_uint &ring, ham_uint &startpix, ham_int &ringpix,
             ham_float &theta, bool &shifted) const {
    const ham_uint northring{(ring > 2 * this->Nside) ? 4 * this->Nside - ring
                                                      : ring};
    if (northring < this->Nside) { // northring < Nside
      const ham_float tmp{northring * northring * this->Fact2};
      const ham_float costheta{1 - tmp};
      const ham_float sintheta{std::sqrt(tmp * (2 - tmp))};
      theta = std::atan2(sintheta, costheta);
      ringpix = 4 * northring;
      shifted = true;
      startpix = 2 * northring * (northring - 1);
    } else { // northring >= Nside
      theta = std::acos((static_cast<ham_float>(2 * this->Nside) -
                         static_cast<ham_float>(northring)) *
                        this->Fact1);
      ringpix = 4 * this->Nside;
      shifted = ((northring - this->Nside) & 1) == 0;
      startpix = this->Ncap + (northring - this->Nside) * ringpix;
    }
    // southern hemisphere extra correction
    if (northring != ring) {
      theta = cgs::pi - theta;
      startpix = this->Npix - startpix - ringpix;
    }
  }

public:
  // dft constr
  Hampix() : Hamdis<T>() {}
  // initialize map with given HEALPix Nside
  // Node's Data assigned by the given value
  // Node's Index assigned from 0 to N-1
  Hampix(const ham_uint &n, const T &v = static_cast<T>(0)) : Hamdis<T>() {
    this->Nside = n;
    this->prepare();
    this->Map = std::make_unique<std::vector<Node<T>>>(
        this->Npix);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (ham_uint i = 0; i < this->Npix; ++i) {
      this->Map->at(i).index(i);
      this->Map->at(i).pointing(this->fillpoint(i));
      this->Map->at(i).data(v);
    }
  }
  Hampix(const ham_uint &n, const std::vector<T> &v) : Hamdis<T>() {
    this->Nside = n;
    this->prepare();
    this->Map = std::make_unique<std::vector<Node<T>>>(
        this->Npix);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (ham_uint i = 0; i < this->Npix; ++i) {
      this->Map->at(i).index(i);
      this->Map->at(i).pointing(this->fillpoint(i));
      this->Map->at(i).data(v[i]);
    }
  }
  // copy constr
  Hampix(const Hampix<T> &m) : Hamdis<T>(m) {
    this->Nside = m.Nside;
    this->Order = m.Order;
    this->Npix = m.Npix;
    this->Npface = m.Npface;
    this->Ncap = m.Ncap;
    this->Fact1 = m.Fact1;
    this->Fact2 = m.Fact2;
  }
  // move constr
  Hampix(Hampix<T> &&m) : Hamdis<T>(std::move(m)) {
    this->Nside = std::move(m.Nside);
    this->Order = std::move(m.Order);
    this->Npix = std::move(m.Npix);
    this->Npface = std::move(m.Npface);
    this->Ncap = std::move(m.Ncap);
    this->Fact1 = std::move(m.Fact1);
    this->Fact2 = std::move(m.Fact2);
  }
  // move assign
  Hampix &operator=(Hampix<T> &&m) noexcept {
    Hamdis<T>::operator=(std::move(m));
    this->Nside = std::move(m.Nside);
    this->Order = std::move(m.Order);
    this->Npix = std::move(m.Npix);
    this->Npface = std::move(m.Npface);
    this->Ncap = std::move(m.Ncap);
    this->Fact1 = std::move(m.Fact1);
    this->Fact2 = std::move(m.Fact2);
    return *this;
  }
  // copy assignment
  Hampix &operator=(const Hampix<T> &m) {
    Hamdis<T>::operator=(m);
    this->Nside = m.Nside;
    this->Order = m.Order;
    this->Npix = m.Npix;
    this->Npface = m.Npface;
    this->Ncap = m.Ncap;
    this->Fact1 = m.Fact1;
    this->Fact2 = m.Fact2;
    return *this;
  }
  // dft destr
  virtual ~Hampix() = default;
  // nside
  ham_uint nside() const { return this->Nside; }
  // npix
  ham_uint npix() const override { return this->Npix; }
  // reset with given nside and clean up data
  void reset(const ham_uint &n = 0) override {
    // cleaning an used map with correct size
    if (n == 0 or this->Nside == n) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (ham_uint i = 0; i < this->Npix; ++i) {
        this->Map->at(i).data(0.0);
      }
    } else {
      // there 2 cases when Nside != n
      // 1, map to be recycled with wrong size
      // 2, empty map initialized by the default constr
      this->Nside = n;
      this->prepare();
      this->Map = std::make_unique<std::vector<Node<T>>>(
          this->Npix);
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (ham_uint i = 0; i < this->Npix; ++i) {
        this->Map->at(i).index(i);
        this->Map->at(i).pointing(this->fillpoint(i));
        this->Map->at(i).data(0.0);
      }
    }
  }
  // rescale
  void rescale(const ham_float &v) override {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (ham_uint i = 0; i < this->Npix; ++i) {
      T current_data = this->Map->at(i).data();
      this->Map->at(i).data(current_data * v);
    }
  }
  // auxliliary function for up/downgrading
  // copied from HEALPix healpix_base.cc ring2xyf function
  void rpixxyf(const ham_uint pix, ham_int &ix, ham_int &iy,
               ham_int &face_num) const {
    ham_int iring, iphi, kshift, nr;
    ham_int Nside_int{
        static_cast<ham_int>(this->Nside)}; // for datatype convenience
    ham_int nl2{2 * Nside_int};
    // North Polar cap
    if (pix < this->Ncap) {
      iring = (1 + this->isqrt(1 + 2 * pix)) >> 1; // counted from North pole
      iphi = (pix + 1) - 2 * iring * (iring - 1);
      kshift = 0;
      nr = iring;
      face_num = this->special_div(iphi - 1, nr);
    }
    // Equatorial region
    else if (pix < (this->Npix - this->Ncap)) {
      ham_int ip{static_cast<ham_int>(pix - this->Ncap)};
      ham_int tmp =
          (this->Order >= 0) ? ip >> (this->Order + 2) : ip / (4 * Nside_int);
      iring = tmp + Nside_int;
      iphi = ip - tmp * 4 * Nside_int + 1;
      kshift = (iring + Nside_int) & 1;
      nr = Nside_int;
      ham_int ire{tmp + 1};
      ham_int irm{nl2 + 1 - tmp};
      ham_int ifm{iphi - (ire >> 1) + Nside_int - 1};
      ham_int ifp{iphi - (irm >> 1) + Nside_int - 1};
      if (this->Order >= 0) {
        ifm >>= this->Order;
        ifp >>= this->Order;
      } else {
        ifm /= Nside_int;
        ifp /= Nside_int;
      }
      face_num = (ifp == ifm) ? (ifp | 4) : ((ifp < ifm) ? ifp : (ifm + 8));
    }
    // South Polar cap
    else {
      ham_int ip{static_cast<ham_int>(this->Npix - pix)};
      iring = (1 + this->isqrt(2 * ip - 1)) >> 1; // counted from South pole
      iphi = 4 * iring + 1 - (ip - 2 * iring * (iring - 1));
      kshift = 0;
      nr = iring;
      iring = 2 * nl2 - iring;
      face_num = this->special_div(iphi - 1, nr) + 8;
    }
    ham_int irt{iring - ((2 + (face_num >> 2)) * Nside_int) + 1};
    ham_int ipt{2 * iphi - this->Jpll[face_num] * nr - kshift - 1};
    //  ham_uint ipt = 2*iphi- (((face_num&3)<<1)+1-((face_num>>2)&1))*nr -
    //  kshift -1;
    if (ipt >= nl2)
      ipt -= 8 * Nside_int;
    ix = (ipt - irt) >> 1;
    iy = (-ipt - irt) >> 1;
  }
  // auxliliary function for up/downgrading
  // copied from HEALPix healpix_base.cc xyf2ring function
  ham_uint xyfrpix(const ham_int &ix, const ham_int &iy,
                   const ham_int &face_num) const {
    ham_uint nl4{4 * this->Nside};
    ham_uint jr{(this->Jrll[face_num] * this->Nside) - ix - iy - 1};
    ham_uint nr, n_before;
    ham_int kshift;
    bool shifted;
    // cpied from get_ring_info_small
    if (jr < this->Nside) {
      shifted = true;
      nr = 4 * jr;
      n_before = 2 * jr * (jr - 1);
    } else if (jr < 3 * this->Nside) {
      shifted = (((jr - this->Nside) & 1) == 0);
      nr = 4 * this->Nside;
      n_before = this->Ncap + (jr - this->Nside) * nr;
    } else {
      shifted = true;
      ham_uint tmp = 4 * this->Nside - jr;
      nr = 4 * tmp;
      n_before = this->Npix - 2 * tmp * (tmp + 1);
    }
    // end of get_ring_info_small
    nr >>= 2;
    kshift = 1 - shifted;
    ham_int jp{(this->Jpll[face_num] * static_cast<ham_int>(nr) + ix - iy + 1 +
                kshift) /
               2};
    assert(jp <= 4 * static_cast<ham_int>(nr));
    if (jp < 1)
      jp += static_cast<ham_int>(
          nl4); // assumption: if this triggers, then nl4==4*nr
    return n_before + jp - 1;
  }
  // interpolate map at given pointing position (linear interpolation with
  // nearby 4 pixels) copied from HEALPix healpix_map.h interpolated_value
  // functions
  T interpolate(const Hamp &point) const {
    // calculate surrounding pixel indices and weights
    // copied from HEALPix healpix_base.cc get_interpol function
    std::array<ham_uint, 4> pix{};
    std::array<ham_float, 4> wght{};
    const ham_float z{std::cos(point.theta())};
    const ham_uint ring1{rabove(z)};
    const ham_uint ring2{ring1 + 1};
    ham_float theta1{0.0}, theta2{0.0};
    // ring above and ring below
    if (ring1 > 0) {
      ham_uint start;
      ham_int ringpix;
      bool shift;
      this->rinfo(ring1, start, ringpix, theta1, shift);
      const ham_float dphi{6.283185307179586476925286766559005768394 / ringpix};
      const ham_float tmp{(point.phi() / dphi - 0.5 * shift)};
      ham_int i1{static_cast<ham_int>(tmp) - static_cast<ham_int>(tmp < 0)};
      const ham_float weight{(point.phi() - (i1 + 0.5 * shift) * dphi) / dphi};
      ham_int i2{i1 + 1};
      i1 += ringpix * (i1 < 0);
      i2 -= ringpix * (i2 >= ringpix);
      pix[0] = start + i1;
      pix[1] = start + i2;
      wght[0] = 1.0 - weight;
      wght[1] = weight;
    }
    if (ring2 < (4 * this->Nside)) {
      ham_uint start;
      ham_int ringpix;
      bool shift;
      this->rinfo(ring2, start, ringpix, theta2, shift);
      const ham_float dphi{6.283185307179586476925286766559005768394 / ringpix};
      const ham_float tmp{(point.phi() / dphi - 0.5 * shift)};
      ham_int i1{static_cast<ham_int>(tmp) - static_cast<ham_int>(tmp < 0)};
      const ham_float weight{(point.phi() - (i1 + 0.5 * shift) * dphi) / dphi};
      ham_int i2{i1 + 1};
      i1 += ringpix * (i1 < 0);
      i2 -= ringpix * (i2 >= ringpix);
      pix[2] = start + i1;
      pix[3] = start + i2;
      wght[2] = 1.0 - weight;
      wght[3] = weight;
    }
    // special cases + post correction
    if (ring1 == 0) {
      const ham_float wtheta{point.theta() / theta2};
      wght[2] *= wtheta;
      wght[3] *= wtheta;
      const ham_float fac{(1.0 - wtheta) * 0.25};
      wght[0] = fac;
      wght[1] = fac;
      wght[2] += fac;
      wght[3] += fac;
      pix[0] = (pix[2] + 2) & 3;
      pix[1] = (pix[3] + 2) & 3;
    } else if (ring2 == (4 * this->Nside)) {
      ham_float wtheta{(point.theta() - theta1) / (cgs::pi - theta1)};
      wght[0] *= (1.0 - wtheta);
      wght[1] *= (1.0 - wtheta);
      ham_float fac{wtheta * 0.25};
      wght[0] += fac;
      wght[1] += fac;
      wght[2] = fac;
      wght[3] = fac;
      pix[2] = ((pix[0] + 2) & 3) + this->Npix - 4;
      pix[3] = ((pix[1] + 2) & 3) + this->Npix - 4;
    } else {
      ham_float wtheta{(point.theta() - theta1) / (theta2 - theta1)};
      wght[0] *= (1.0 - wtheta);
      wght[1] *= (1.0 - wtheta);
      wght[2] *= wtheta;
      wght[3] *= wtheta;
    }
    // calculate interpolated result
    // copied from HEALPix healpix_map.h interpolation function
    ham_float wtot{0.0};
    T res{static_cast<T>(0)};
    for (int i = 0; i < 4; ++i) {
      T val{this->Map->at(pix[i]).data()};
      if (val > -1.63749e30) { // larger than undef, exclude the masked
        res += val * wght[i];
        wtot += wght[i];
      }
    }
    return (wtot == 0.0) ? static_cast<T>(this->undef)
                         : static_cast<T>(res / wtot);
  }
  // add maps
  void accumulate(const Hampix<T> &m) {
    // in same Nside
    if (this->Npix == m.npix()) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (ham_uint i = 0; i < this->Npix; ++i) {
        const T cache = this->Map->at(i).data();
        this->Map->at(i).data(cache + m.data(i));
      }
    }
    // in different Nside
    // use interpolation
    else {
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (ham_uint i = 0; i < this->Npix; ++i) {
        const T cache = this->Map->at(i).data();
        this->Map->at(i).data(cache +
                              m.interpolate(this->Map->at(i).pointing()));
      }
    }
  }
};

#endif
