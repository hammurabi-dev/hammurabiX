// (Gaussian) random GMF generators
// field parameterized modelling is implemented in breg function
// interactions with grid are carried out via read_grid and write_grid
// notice that read/write_grid acts with internal grid in memory
// interactions with disk data are carried out in Grid_brnd class
// get_breg is the universal interface for retrieving field vector
// at given 3D spatial position
// get_breg knows where to read the information
// either from grid or parameterzed modelling
//
// Brnd_global and Brnd_local are ghosts (of no technical usage)
// designed for distinguish global and local random GMF generators
//
// it is impossible to describe Brnd_es/jaffe/mhd in details here
// please check the wiki page
//
// customized field modelling can be added as derived class of Brnd

#ifndef HAMMURABI_BRND_H
#define HAMMURABI_BRND_H

#include <breg.h>
#include <grid.h>
#include <hvec.h>
#include <param.h>

// base class with read_grid implemented
// get_brnd is invoked when no specific derived class object is instantiated
class Brnd {
public:
  Brnd() = default;
  Brnd(const Brnd &) = delete;
  Brnd(Brnd &&) = delete;
  Brnd &operator=(const Brnd &) = delete;
  Brnd &operator=(Brnd &&) = delete;
  virtual ~Brnd() = default;
  // get random field vector
  // return zero field when read_grid is not invoked
  // 1st argument: Galactic centric Cartesian frame position
  // 2nd argument: random GMF grid object
  virtual hvec<3, double> get_brnd(const hvec<3, double> &, const Param *,
                                   const Grid_brnd *) const;
  // read field from grid with trilinear interpolation
  // user has to call write_grid ahead in main routine
  // 1st argument: Galactic centric Cartesian frame position
  // 2nd argument: random GMF grid object
  virtual hvec<3, double> read_grid(const hvec<3, double> &, const Param *,
                                    const Grid_brnd *) const;
  // write field to grid (model dependent)
  // user can export_grid to binary file with Grid_xxx::export_grid
  // for dynamic binding only, implemented in derived class
  // 1st argument: parameter class object
  // 2nd argument: regular GMF class object
  // 3rd argument: regular GMF grid class object
  // 4th argument: random GMF grid class object
  virtual void write_grid(const Param *, const Breg *, const Grid_breg *,
                          Grid_brnd *) const;
};

// global (an)isotropic turbulent GMF
// this class is treated as a covering class for specific methods
class Brnd_global : public Brnd {
public:
  Brnd_global() = default;
  Brnd_global(const Brnd_global &) = delete;
  Brnd_global(Brnd_global &&) = delete;
  Brnd_global &operator=(const Brnd_global &) = delete;
  Brnd_global &operator=(Brnd_global &&) = delete;
  virtual ~Brnd_global() = default;
};

// local (an)isotropic turbulent GMF
// this class is treated as a covering class for specific methods
class Brnd_local : public Brnd {
public:
  Brnd_local() = default;
  Brnd_local(const Brnd_local &) = delete;
  Brnd_local(Brnd_local &&) = delete;
  Brnd_local &operator=(const Brnd_local &) = delete;
  Brnd_local &operator=(Brnd_local &&) = delete;
  virtual ~Brnd_local() = default;
};

// Ensslin-Steininger method of global (an)isotropic turbulent GMF
class Brnd_es final : public Brnd_global {
public:
  Brnd_es() = default;
  Brnd_es(const Brnd_es &) = delete;
  Brnd_es(Brnd_es &&) = delete;
  Brnd_es &operator=(const Brnd_es &) = delete;
  Brnd_es &operator=(Brnd_es &&) = delete;
  virtual ~Brnd_es() = default;
  // use triple Fourier transform scheme
  // check technical report for details
  void write_grid(const Param *, const Breg *, const Grid_breg *,
                  Grid_brnd *) const override;
#ifndef NDEBUG
protected:
#endif
  // isotropic power-spectrum
  // 1st argument: isotropic wave-vector magnitude
  // 2nd argument: parameter class object
  virtual double spec(const double &, const Param *) const;
  // field energy density rescaling factor
  // 1st argument: Galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  virtual double rescal(const hvec<3, double> &, const Param *) const;
  // anisotropy factor
  // check technical report for details
  // 1st argument: Galactic centric Cartesian frame position
  // 2rd argument: parameter class object
  // 3th argument: regular GMF class object
  // 4th argument: regular GMF grid class object
  hvec<3, double> anisotropy_direction(const hvec<3, double> &, const Param *,
                                       const Breg *, const Grid_breg *) const;
  // anisotropy ratio
  // 1st argument: Galactic centric Cartesian frame position
  // 2rd argument: parameter class object
  // 3th argument: regular GMF class object
  // 4th argument: regular GMF grid class object
  double anisotropy_ratio(const hvec<3, double> &, const Param *, const Breg *,
                          const Grid_breg *) const;
  // Gram-Schmidt orthogonalization process
  // 1st argument: wave-vector
  // 2nd arugment: input GMF vector (in Fourier space)
  // remark: real and imagine parts of complex GMF vector in Fourier space
  // handled separately
  hvec<3, double> gramschmidt(const hvec<3, double> &,
                              const hvec<3, double> &) const;
};

// Jaffe method of global (an)isotropic turbulent GMF
// class Brnd_Jaffe final : public Brnd_global

// local anisotropic random GMF
// in compressive MHD plasma
// check technical report for details
class Brnd_mhd final : public Brnd_local {
public:
  Brnd_mhd() = default;
  Brnd_mhd(const Brnd_mhd &) = delete;
  Brnd_mhd(Brnd_mhd &&) = delete;
  Brnd_mhd &operator=(const Brnd_mhd &) = delete;
  Brnd_mhd &operator=(Brnd_mhd &&) = delete;
  virtual ~Brnd_mhd() = default;
  // use vector field decomposition scheme
  // use regular GMF at position of the Sun
  void write_grid(const Param *, const Breg *, const Grid_breg *,
                  Grid_brnd *) const override;
#ifndef NDEBUG
protected:
#endif
  // dynamo number
  // 1st argument: plasma beta
  // 2nd argument: cosine of k-B pitch angle
  inline double dynamo(const double &beta, const double &cosa) const {
    return (1 + 0.5 * beta) * (1 + 0.5 * beta) - 2. * beta * cosa * cosa;
  }
  // fast mode anisotropic tensor structure
  // 1st argument: plasma beta
  // 2nd argument: cosine of k-B pitch angle
  double hf(const double &, const double &) const;
  // slow mode anisotropic tensor structure
  // 1st argument: plasma beta
  // 2nd argument: cosine of k-B pitch angle
  double hs(const double &, const double &) const;
  // fast mode anisotropy power factor
  // 1st argument: plasma March number
  // 2nd argument: cosine of k-B pitch angle
  inline double fa(const double &ma, const double &cosa) const {
    return std::exp(-1. * std::pow(ma, -1.33333333) * cosa * cosa /
                    std::pow(1 - cosa * cosa, 0.66666667));
  }
  // slow mode anisotropy power factor
  // 1st argument: plasma March number
  // 2nd argument: cosine of k-B pitch angle
  inline double fs(const double &ma, const double &cosa) const {
    return std::exp(-1. * std::pow(ma, -1.33333333) * cosa * cosa /
                    std::pow(1 - cosa * cosa, 0.66666667));
  }
  // cosine of pitch angle between wavevector and regular GMF
  // 1st argument: field vector
  // 2nd argument: wave vector
  inline double cosa(const hvec<3, double> &b, const hvec<3, double> &k) const {
    // dotprod is function from vec3.h
    return (b.versor()).dotprod(k.versor());
  }
  // direction of Alfven mode
  // 1st argument: regualr GMF vector
  // 2nd argument: wave vector
  hvec<3, double> eplus(const hvec<3, double> &, const hvec<3, double> &) const;
  // direction of slow and fast modes
  // 1st argument: regualr GMF vector
  // 2nd argument: wave vector
  hvec<3, double> eminus(const hvec<3, double> &,
                         const hvec<3, double> &) const;
  // isotropic part of power spectrum of Alfvenic mode
  // 1st argument: wave vector magnitude
  // 2nd argument: parameter class object
  double speca(const double &, const Param *) const;
  // isotropic part of power spectrum of fast mode
  // 1st argument: wave vector magnitude
  // 2nd argument: parameter class object
  double specf(const double &, const Param *) const;
  // isotropic part of power spectrum of slow mode
  // 1st argument: wave vector magnitude
  // 2nd argument: parameter class object
  double specs(const double &, const Param *) const;
};

#endif

// END
