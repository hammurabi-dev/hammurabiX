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
  // read CRE flux from grid at given spatial position and energy
  // actual value of E is calculated from {E_index,Ek_min,Ek_fact}
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: energy
  // 3rd argument: parameter class object
  // 4th argument: CRE grid class object
  virtual double read_grid(const hamvec<3, double> &, const double &,
                           const Param *, const Grid_cre *) const;
  // read CRE flux from grid at given energy index and spatial position
  // mind its difference to ``read_grid``
  // for reading CRE flux at arbitrary energy, use base class read_field
  // actual value of E is calculated from {E_index,Ek_min,Ek_fact}
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: index in energy
  // 3rd argument: parameter class object
  // 4th argument: CRE grid class object
  double read_grid_num(const hamvec<3, double> &, const std::size_t &,
                       const Param *, const Grid_cre *) const;
  // fill the grid with CRE flux distribution
  // 1st argument: parameter class object
  // 2nd argument: CRE grid class object
  virtual void write_grid(const Param *, Grid_cre *) const;
  // get CRE flux at given CRE energy and spatial position
  // input CRE energy at CGS units
  // output in [GeV m^2 s sr]^-1 units
  // notice that ``read_field`` function may not be used
  // in calculating the synchrotron emissivity
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: CRE energy in GeV
  // 3rd argument: parameter class object
  // 4th argument: CRE grid class object
  virtual double read_field(const hamvec<3, double> &, const double &,
                            const Param *, const Grid_cre *) const;
  // assemble CRE phase-space density at given position
  // 1st argument: galactic centric Cartesian spatial position
  // 2nd argument: CRE energy in GeV
  // 3rd argument: parameter class object
  virtual double write_field(const hamvec<3, double> &, const double &,
                             const Param *) const;
  // flux normalization at given position
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  virtual double flux_norm(const hamvec<3, double> &, const Param *) const;
  // flux index at given position
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  virtual double flux_idx(const hamvec<3, double> &, const Param *) const;
  // spatial CRE flux reprofiling
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  virtual double spatial_profile(const hamvec<3, double> &,
                                 const Param *) const;
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
  double write_field(const hamvec<3, double> &, const double &,
                     const Param *) const override;
  // flux normalization at given position
  double flux_norm(const hamvec<3, double> &, const Param *) const override;
  // flux index at given position
  double flux_idx(const hamvec<3, double> &, const Param *) const override;
  // spatial CRE flux reprofiling
  double spatial_profile(const hamvec<3, double> &,
                         const Param *) const override;
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
  double write_field(const hamvec<3, double> &, const double &,
                     const Param *) const override;
  // flux normalization at given position
  double flux_norm(const hamvec<3, double> &, const Param *) const override;
  // flux index at given position
  double flux_idx(const hamvec<3, double> &, const Param *) const override;
  // spatial CRE flux reprofiling
  double spatial_profile(const hamvec<3, double> &,
                         const Param *) const override;
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
  // overload the base class read_grid function
};

#endif
