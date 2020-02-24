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
#include <hamtype.h>
#include <hamvec.h>
#include <param.h>

// base class, all functions are implemented in derived class
class CREfield {
public:
  CREfield() = default;
  CREfield(const CREfield &) = delete;
  CREfield(CREfield &&) = delete;
  CREfield &operator=(const CREfield &) = delete;
  CREfield &operator=(CREfield &&) = delete;
  virtual ~CREfield() = default;
  // read CRE flux from grid at given spatial position and energy
  // actual value of E is calculated from {E_index,Ek_min,Ek_fact}
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: energy
  // 3rd argument: parameter class object
  // 4th argument: CRE grid class object
  virtual ham_float read_grid(const Hamvec<3, ham_float> &, const ham_float &,
                              const Param *, const Grid_cre *) const;
  // read CRE flux from grid at given energy index and spatial position
  // mind its difference to ``read_grid``
  // for reading CRE flux at arbitrary energy, use base class read_field
  // actual value of E is calculated from {E_index,Ek_min,Ek_fact}
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: index in energy
  // 3rd argument: parameter class object
  // 4th argument: CRE grid class object
  ham_float read_grid_num(const Hamvec<3, ham_float> &, const ham_uint &,
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
  virtual ham_float read_field(const Hamvec<3, ham_float> &, const ham_float &,
                               const Param *, const Grid_cre *) const;
  // assemble CRE phase-space density at given position
  // 1st argument: galactic centric Cartesian spatial position
  // 2nd argument: CRE energy in GeV
  // 3rd argument: parameter class object
  virtual ham_float write_field(const Hamvec<3, ham_float> &, const ham_float &,
                                const Param *) const;
  // flux normalization at given position
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  virtual ham_float flux_norm(const Hamvec<3, ham_float> &,
                              const Param *) const;
  // flux index at given position
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  virtual ham_float flux_idx(const Hamvec<3, ham_float> &, const Param *) const;
  // spatial CRE flux reprofiling
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  virtual ham_float spatial_profile(const Hamvec<3, ham_float> &,
                                    const Param *) const;
};

// uniform CRE flux
class CRE_unif : public CREfield {
public:
  CRE_unif() = default;
  CRE_unif(const CRE_unif &) = delete;
  CRE_unif(CRE_unif &&) = delete;
  CRE_unif &operator=(const CRE_unif &) = delete;
  CRE_unif &operator=(CRE_unif &&) = delete;
  virtual ~CRE_unif() = default;
  ham_float write_field(const Hamvec<3, ham_float> &, const ham_float &,
                        const Param *) const override;
  // flux normalization at given position
  ham_float flux_norm(const Hamvec<3, ham_float> &,
                      const Param *) const override;
  // flux index at given position
  ham_float flux_idx(const Hamvec<3, ham_float> &,
                     const Param *) const override;
  // spatial CRE flux reprofiling
  ham_float spatial_profile(const Hamvec<3, ham_float> &,
                            const Param *) const override;
};

// analytic CRE flux
class CRE_ana : public CREfield {
public:
  CRE_ana() = default;
  CRE_ana(const CRE_ana &) = delete;
  CRE_ana(CRE_ana &&) = delete;
  CRE_ana &operator=(const CRE_ana &) = delete;
  CRE_ana &operator=(CRE_ana &&) = delete;
  virtual ~CRE_ana() = default;
  ham_float write_field(const Hamvec<3, ham_float> &, const ham_float &,
                        const Param *) const override;
  // flux normalization at given position
  ham_float flux_norm(const Hamvec<3, ham_float> &,
                      const Param *) const override;
  // flux index at given position
  ham_float flux_idx(const Hamvec<3, ham_float> &,
                     const Param *) const override;
  // spatial CRE flux reprofiling
  ham_float spatial_profile(const Hamvec<3, ham_float> &,
                            const Param *) const override;
};

// use numerical CRE flux
class CRE_num final : public CREfield {
public:
  CRE_num() = default;
  CRE_num(const CRE_num &) = delete;
  CRE_num(CRE_num &&) = delete;
  CRE_num &operator=(const CRE_num &) = delete;
  CRE_num &operator=(CRE_num &&) = delete;
  virtual ~CRE_num() = default;
};

#endif
