// line-of-sight radiative transfer integrator

#ifndef HAMMURABI_INT_H
#define HAMMURABI_INT_H

#include <vector>

#include <bfield.h>
#include <crefield.h>
#include <grid.h>
#include <hamp.h>
#include <hamtype.h>
#include <param.h>
#include <tefield.h>

class Integrator {
public:
  Integrator() = default;
  Integrator(const Integrator &) = delete;
  Integrator(Integrator &&) = delete;
  Integrator &operator=(const Integrator &) = delete;
  Integrator &operator=(Integrator &&) = delete;
  virtual ~Integrator() = default;
  // assmebling pixels/shells into sky map
  void write_grid(const Breg *, const Brnd *, const TEreg *, const TErnd *,
                  const CREfield *, const Grid_breg *, const Grid_brnd *,
                  const Grid_tereg *, const Grid_ternd *, const Grid_cre *,
                  Grid_obs *, const Param *) const;
#ifdef NDEBUG
protected:
#endif
  // Carteisan unit vector of given LoS direction
  // 1st argument: polar angle (in rad)
  // 2nd argument: azimuthal angle (in rad)
  Hamvec<3, ham_float> los_versor(const ham_float &, const ham_float &) const;
  // perpendicular projection of a vector wrt LoS direction
  // 1st argument: vector in Cartesian frame
  // 2nd argument: los versor
  inline ham_float los_perproj(const Hamvec<3, ham_float> &input,
                               const Hamvec<3, ham_float> &ver) const {
    return (input.crossprod(ver)).length();
  }
  // (signed) parallel projection of a vector wrt LoS direction
  // 1st argument: vector in Cartesian frame
  // 2nd argument: los versor
  inline ham_float los_parproj(const Hamvec<3, ham_float> &input,
                               const Hamvec<3, ham_float> &ver) const {
    return input.dotprod(ver);
  }
  // synchrotron emission intrinsic polarization angle (in rad)
  // check Rybicki & Lightman Sec.6.5 'polarization of synchrotron radiation'
  // 1st argument: magnetic field in Cartesian frame
  // 2nd argument: polar angle (in rad) of LOS direction
  // 3rd argument: azimuthal angle (in rad) of LOS direction
  // use with caution, since vector can be parallel to LoS direction
  ham_float sync_ipa(const Hamvec<3, ham_float> &, const ham_float &,
                     const ham_float &) const;
  // synchrotron total emissivity calculator
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  // 3rd argument: CRE object
  // 4th argument: CRE grid object
  // 5th argument: perpendicular component of magnetic field wrt LoS direction
  ham_float sync_emissivity_t(const Hamvec<3, ham_float> &, const Param *,
                              const CREfield *, const Grid_cre *,
                              const ham_float &) const;
  // synchrotron polarized emissivity calculator
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  // 3rd argument: CRE object
  // 4th argument: CRE grid object
  // 5th argument: perpendicular component of magnetic field wrt LoS direction
  ham_float sync_emissivity_p(const Hamvec<3, ham_float> &, const Param *,
                              const CREfield *, const Grid_cre *,
                              const ham_float &) const;
  // converting brightness temp into thermal temp with T_0 = 2.725K,
  // Prog.Theor.Exp.Phys. (2014) 2014 (6): 06B109.
  // 1st argument: brightness temperature
  // 2nd argument: observational frequency
  ham_float temp_convert(const ham_float &, const ham_float &) const;
#ifdef NDEBUG
private:
#endif
  // to hold temporary information of observables
  struct struct_observables {
    ham_float is, qs, us;
    ham_float fd;
    ham_float dm;
    ham_float ff;
  };
  // to hold temporary information of spherical shells
  struct struct_shell {
    ham_uint shell_num;
    ham_float d_start;
    ham_float d_stop;
    ham_float delta_d;
    ham_uint step;
    std::vector<ham_float> dist;
  };
  // conduct LOS integration in one pixel at given shell
  void radial_integration(const struct_shell *, const Hamp &,
                          struct_observables *, const Breg *, const Brnd *,
                          const TEreg *, const TErnd *, const CREfield *,
                          const Grid_breg *, const Grid_brnd *,
                          const Grid_tereg *, const Grid_ternd *,
                          const Grid_cre *, const Param *) const;
  // general upper boundary check
  // return false if 1st argument is larger than 2nd
  inline bool check_simulation_upper_limit(const ham_float &value,
                                           const ham_float &limit) const {
    return (value > limit);
  }
  // general lower boundary check
  // return false if 1st argument is smaller than 2nd
  inline bool check_simulation_lower_limit(const ham_float &value,
                                           const ham_float &limit) const {
    return (value < limit);
  }
  // assembling ``struct_shell``
  // this part may introduce precision loss
  void assemble_shell_ref(struct_shell *, const Param *,
                          const ham_uint &) const;
};

#endif
