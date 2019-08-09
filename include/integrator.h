// line-of-sight radiative transfer integrator

#ifndef HAMMURABI_INT_H
#define HAMMURABI_INT_H

#include <vector>

#include <bfield.h>
#include <crefield.h>
#include <grid.h>
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
  void write_grid(Breg *, Brnd *, TEreg *, TErnd *, CRE *, Grid_breg *,
                  Grid_brnd *, Grid_tereg *, Grid_ternd *, Grid_cre *,
                  Grid_obs *, const Param *) const;
#ifdef NDEBUG
protected:
#endif
  // Carteisan unit vector of given LoS direction
  // 1st argument: polar angle (in rad)
  // 2nd argument: azimuthal angle (in rad)
  hamvec<3, double> los_versor(const double &, const double &) const;
  // perpendicular projection of a vector wrt LoS direction
  // 1st argument: vector in Cartesian frame
  // 2nd argument: polar angle (in rad)
  // 3rd argument: azimuthal angle (in rad)
  double los_perproj(const hamvec<3, double> &, const double &,
                     const double &) const;
  // (signed) parallel projection of a vector wrt LoS direction
  // 1st argument: vector in Cartesian frame
  // 2nd argument: polar angle (in rad)
  // 3rd argument: azimuthal angle (in rad)
  double los_parproj(const hamvec<3, double> &, const double &,
                     const double &) const;
  // synchrotron emission intrinsic polarization angle (in rad)
  // check Rybicki & Lightman Sec.6.5 'polarization of synchrotron radiation'
  // 1st argument: magnetic field in Cartesian frame
  // 2nd argument: polar angle (in rad) of LOS direction
  // 3rd argument: azimuthal angle (in rad) of LOS direction
  // use with caution, since vector can be parallel to LoS direction
  double sync_ipa(const hamvec<3, double> &, const double &,
                  const double &) const;
  // synchrotron total emissivity calculator
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  // 3rd argument: CRE object
  // 4th argument: CRE grid object
  // 5th argument: perpendicular component of magnetic field wrt LoS direction
  double sync_emissivity_t(const hamvec<3, double> &, const Param *,
                           const CRE *, const Grid_cre *, const double &) const;
  // synchrotron polarized emissivity calculator
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  // 3rd argument: CRE object
  // 4th argument: CRE grid object
  // 5th argument: perpendicular component of magnetic field wrt LoS direction
  double sync_emissivity_p(const hamvec<3, double> &, const Param *,
                           const CRE *, const Grid_cre *, const double &) const;
  // dust emission intrinsic polarization angle (in rad)
  // check Planck XX, A&A 576, A105 (2015)
  // 1st argument: magnetic field in Cartesian frame
  // 2nd argument: polar angle (in rad) of LOS direction
  // 3rd argument: azimuthal angle (in rad) of LOS direction
  // use with caution, since vector can be parallel to LoS direction
  double dust_ipa(const hamvec<3, double> &, const double &,
                  const double &) const;
  // converting brightness temp into thermal temp with T_0 = 2.725K,
  // Prog.Theor.Exp.Phys. (2014) 2014 (6): 06B109.
  // 1st argument: brightness temperature
  // 2nd argument: observational frequency
  double temp_convert(const double &, const double &) const;
#ifdef NDEBUG
private:
#endif
  // to hold temporary information of observables
  struct struct_observables {
    double is, qs, us;
    double fd;
    double dm;
    double ff;
  };
  // to hold temporary information of spherical shells
  struct struct_shell {
    std::size_t shell_num;
    double d_start;
    double d_stop;
    double delta_d;
    std::size_t step;
    std::vector<double> dist;
  };
  // conduct LOS integration in one pixel at given shell
  void radial_integration(struct_shell *, pointing &, struct_observables &,
                          Breg *, Brnd *, TEreg *, TErnd *, CRE *, Grid_breg *,
                          Grid_brnd *, Grid_tereg *, Grid_ternd *, Grid_cre *,
                          const Param *) const;
  // general upper boundary check
  // return false if 1st argument is larger than 2nd
  inline bool check_simulation_upper_limit(const double &value,
                                           const double &limit) const {
    return (value > limit);
  }
  // general lower boundary check
  // return false if 1st argument is smaller than 2nd
  inline bool check_simulation_lower_limit(const double &value,
                                           const double &limit) const {
    return (value < limit);
  }
  // assembling ``struct_shell``
  // this part may introduce precision loss
  void assemble_shell_ref(struct_shell *, const Param *,
                          const std::size_t &) const;
};

#endif
