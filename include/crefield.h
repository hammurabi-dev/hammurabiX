// cosmic ray electron flux field generator/interface
//
// class design:
//
// CREfield
//    |-- CRE_xxx
//    |-- CRE_yyy
//
// CRE class either defines CRE flux analytically
// or read from numerical input
// CRE synchrotron emissivity is also defined here
// by function ``read_emissivity_t`` and ``read_emissivity_p``
// for reading out total and polarized emissivity respectively

#ifndef HAMMURABI_CRE_H
#define HAMMURABI_CRE_H

#include <grid.h>
#include <hamvec.h>
#include <param.h>

// base class, all functions are implemented in derived class
class CRE {
public:
  CRE() = default;
  CRE(const CRE &) = delete;
  CRE(CRE &&) = delete;
  CRE &operator=(const CRE &) = delete;
  CRE &operator=(CRE &&) = delete;
  virtual ~CRE() = default;
  // get CRE synchrotron total emissivity
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  // 3rd argument: CRE grid object
  // 4th argument: perpendicular component of magnetic field wrt LoS direction
  virtual double read_emissivity_t(const hamvec<3, double> &, const Param *,
                                   const Grid_cre *, const double &) const;
  // get CRE synchrotron polarized emissivity
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  // 3rd argument: CRE grid object
  // 4th argument: perpendicular component of magnetic field wrt LoS direction
  virtual double read_emissivity_p(const hamvec<3, double> &, const Param *,
                                   const Grid_cre *, const double &) const;
  // read CRE flux from grid at given position
  // (E_index, sylindrical_r, sylindrical_z) with {r,z} in cgs units,
  // actual value of E is calculated from {E_index,Ek_min,Ek_fact}
  // in read_emissivity automatically select bi/trilinear interpolation
  // according to 2+1/3+1 spatial-spectral CRE flux grid
  // 1st argument: index in energy
  // 2nd argument: galactic centric Cartesian frame position
  // 3rd argument: parameter class object
  // 4th argument: CRE grid class object
  virtual double read_grid(const std::size_t &, const hamvec<3, double> &,
                           const Param *, const Grid_cre *) const;
  // fill the grid with CRE flux distribution
  // 1st argument: parameter class object
  // 2nd argument: CRE grid class object
  virtual void write_grid(const Param *, Grid_cre *) const;
  // calculate CRE flux at given CRE energy,
  // input CRE energy at CGS units,
  // output in [GeV m^2 s sr]^-1 units
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  // 3rd argument: CRE energy in GeV
  virtual double flux(const hamvec<3, double> &, const Param *,
                      const double &) const;
};

// uniform CRE flux
class CRE_unif : public CRE {
public:
  CRE_unif() = default;
  CRE_unif(const CRE_unif &) = delete;
  CRE_unif(CRE_unif &&) = delete;
  CRE_unif &operator=(const CRE_unif &) = delete;
  CRE_unif &operator=(CRE_unif &&) = delete;
  virtual ~CRE_unif() = default;
  double read_emissivity_t(const hamvec<3, double> &, const Param *,
                           const Grid_cre *, const double &) const override;
  double read_emissivity_p(const hamvec<3, double> &, const Param *,
                           const Grid_cre *, const double &) const override;
  double flux(const hamvec<3, double> &, const Param *,
              const double &) const override;
#ifdef NDEBUG
protected:
#endif
  // flux normalization at given position
  double flux_norm(const hamvec<3, double> &, const Param *) const;
  // flux index at given position
  double flux_idx(const hamvec<3, double> &, const Param *) const;
  // spatial CRE flux reprofiling
  double spatial_profile(const hamvec<3, double> &, const Param *) const;
};

// analytic CRE flux
class CRE_ana : public CRE {
public:
  CRE_ana() = default;
  CRE_ana(const CRE_ana &) = delete;
  CRE_ana(CRE_ana &&) = delete;
  CRE_ana &operator=(const CRE_ana &) = delete;
  CRE_ana &operator=(CRE_ana &&) = delete;
  virtual ~CRE_ana() = default;
  double read_emissivity_t(const hamvec<3, double> &, const Param *,
                           const Grid_cre *, const double &) const override;
  double read_emissivity_p(const hamvec<3, double> &, const Param *,
                           const Grid_cre *, const double &) const override;
  double flux(const hamvec<3, double> &, const Param *,
              const double &) const override;
#ifdef NDEBUG
protected:
#endif
  // flux normalization at given position
  double flux_norm(const hamvec<3, double> &, const Param *) const;
  // flux index at given position
  double flux_idx(const hamvec<3, double> &, const Param *) const;
  // spatial CRE flux reprofiling
  double spatial_profile(const hamvec<3, double> &, const Param *) const;
};

// use numerical CRE flux
class CRE_num final : public CRE {
public:
  CRE_num() = default;
  CRE_num(const CRE_num &) = delete;
  CRE_num(CRE_num &&) = delete;
  CRE_num &operator=(const CRE_num &) = delete;
  CRE_num &operator=(CRE_num &&) = delete;
  virtual ~CRE_num() = default;
  double read_emissivity_t(const hamvec<3, double> &, const Param *,
                           const Grid_cre *, const double &) const override;
  double read_emissivity_p(const hamvec<3, double> &, const Param *,
                           const Grid_cre *, const double &) const override;
};

#endif
