// magnetic field generators/interface
//
// class design:
//
// Bfield
//   |-- Breg
//   |     |-- Breg_xxx
//   |     |-- Breg_yyy
//   |
//   |-- Brnd
//         |-- Brnd_xxx
//         |-- Brnd_yyy
//
// field parameterized modelling is implemented in ``write_field`` function
// interactions with grid are carried out via ``read_grid`` and ``write_grid``
// notice that ``read/write_grid`` acts with the internal grid in memory
// interactions with disk data are carried out within Grid_bfield class
// ``read_field`` is the universal interface for retrieving field vector
// at given 3D spatial position
// ``read_field`` knows where to read the information
// either from grid or parameterzed modelling

#ifndef HAMMURABI_BFIELD_H
#define HAMMURABI_BFIELD_H

#include <grid.h>
#include <hamvec.h>
#include <param.h>

// magnetic field base class
class Bfield {
public:
  Bfield() = default;
  Bfield(const Bfield &) = delete;
  Bfield(Bfield &&) = delete;
  Bfield &operator=(const Bfield &) = delete;
  Bfield &operator=(Bfield &&) = delete;
  virtual ~Bfield() = default;
};

// regular magnetic field
class Breg : public Bfield {
public:
  Breg() = default;
  Breg(const Breg &) = delete;
  Breg(Breg &&) = delete;
  Breg &operator=(const Breg &) = delete;
  Breg &operator=(Breg &&) = delete;
  virtual ~Breg() = default;
  // get regular magnetic field
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  // 3rd argument: magnetic grid class object
  // read from grid of field type if permitted
  // otherwise use the field assembler
  virtual hamvec<3, double> read_field(const hamvec<3, double> &, const Param *,
                                       const Grid_breg *) const;
  // assemble analytic magnetic field, specified only in derived class
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  virtual hamvec<3, double> write_field(const hamvec<3, double> &,
                                        const Param *) const;
  // read from field grid with trilinear interpolation
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  // 3rd argument: magnetic field grid class object
  virtual hamvec<3, double> read_grid(const hamvec<3, double> &, const Param *,
                                      const Grid_breg *) const;
  // write to field grid
  // 1st argument: parameter class object
  // 2nd argument: magnetic field grid class object
  virtual void write_grid(const Param *, Grid_breg *) const;
};

// random magnetic field
class Brnd {
public:
  Brnd() = default;
  Brnd(const Brnd &) = delete;
  Brnd(Brnd &&) = delete;
  Brnd &operator=(const Brnd &) = delete;
  Brnd &operator=(Brnd &&) = delete;
  virtual ~Brnd() = default;
  // get random magnetic field
  virtual hamvec<3, double> read_field(const hamvec<3, double> &, const Param *,
                                       const Grid_brnd *) const;
  // read from field grid with trilinear interpolation
  virtual hamvec<3, double> read_grid(const hamvec<3, double> &, const Param *,
                                      const Grid_brnd *) const;
  // write random field to grid (model dependent)
  virtual void write_grid(const Param *, const Breg *, const Grid_breg *,
                          Grid_brnd *) const;
};

//------------------------------ Breg DERIVED --------------------------------//

// uniform regular field modeling
// with fixed field orientation and strength
class Breg_unif final : public Breg {
public:
  Breg_unif() = default;
  Breg_unif(const Breg_unif &) = delete;
  Breg_unif(Breg_unif &&) = delete;
  Breg_unif &operator=(const Breg_unif &) = delete;
  Breg_unif &operator=(Breg_unif &&) = delete;
  virtual ~Breg_unif() = default;
  hamvec<3, double> write_field(const hamvec<3, double> &,
                                const Param *) const override;
};

// WMAP LSA modeling
// http://iopscience.iop.org/article/10.1086/513699/meta
// with errata for GMF modeling
// https://lambda.gsfc.nasa.gov/product/map/dr2/pub_papers/threeyear/polarization/errata.cfm
class Breg_wmap final : public Breg {
public:
  Breg_wmap() = default;
  Breg_wmap(const Breg_wmap &) = delete;
  Breg_wmap(Breg_wmap &&) = delete;
  Breg_wmap &operator=(const Breg_wmap &) = delete;
  Breg_wmap &operator=(Breg_wmap &&) = delete;
  virtual ~Breg_wmap() = default;
  hamvec<3, double> write_field(const hamvec<3, double> &,
                                const Param *) const override;
};

// Jaffe modeling
// https://academic.oup.com/mnras/article/401/2/1013/1150693
// https://www.aanda.org/articles/aa/abs/2016/12/aa28033-15/aa28033-15.html
class Breg_jaffe final : public Breg {
public:
  Breg_jaffe() = default;
  Breg_jaffe(const Breg_jaffe &) = delete;
  Breg_jaffe(Breg_jaffe &&) = delete;
  Breg_jaffe &operator=(const Breg_jaffe &) = delete;
  Breg_jaffe &operator=(Breg_jaffe &&) = delete;
  virtual ~Breg_jaffe() = default;
  hamvec<3, double> write_field(const hamvec<3, double> &,
                                const Param *) const override;
#ifndef NDEBUG
protected:
#endif
  // field orientation
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  hamvec<3, double> orientation(const hamvec<3, double> &, const Param *) const;
  // field amplitude radial scaling
  // 1st argument: Galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  double radial_scaling(const hamvec<3, double> &, const Param *) const;
  // spiral arm height scaling
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  // remark: inlined function
  inline double arm_scaling(const hamvec<3, double> &pos,
                            const Param *par) const {
    return 1. / (cosh(pos[2] / par->breg_jaffe.arm_z0) *
                 cosh(pos[2] / par->breg_jaffe.arm_z0));
  }
  // disk height scaling
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  // remark: inlined function
  inline double disk_scaling(const hamvec<3, double> &pos,
                             const Param *par) const {
    return 1. / (cosh(pos[2] / par->breg_jaffe.disk_z0) *
                 cosh(pos[2] / par->breg_jaffe.disk_z0));
  }
  // halo height scaling
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  // remark: inlined function
  inline double halo_scaling(const hamvec<3, double> &pos,
                             const Param *par) const {
    return 1. / (cosh(pos[2] / par->breg_jaffe.halo_z0) *
                 cosh(pos[2] / par->breg_jaffe.halo_z0));
  }
  // spiral arm compression factor, for each arm
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  std::vector<double> arm_compress(const hamvec<3, double> &,
                                   const Param *) const;
  // spiral arm compression factor for dust, for each arm
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  std::vector<double> arm_compress_dust(const hamvec<3, double> &,
                                        const Param *) const;
  // distance to each spiral arm
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  std::vector<double> dist2arm(const hamvec<3, double> &, const Param *) const;
};

//------------------------------ Brnd DERIVED --------------------------------//

// global (an)isotropic random magnetic field
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

// local (an)isotropic random magnetic field
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

// Ensslin-Steininger method of global (an)isotropic random magnetic field
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
  virtual double spectrum(const double &, const Param *) const;
  // field energy density reprofiling factor
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  virtual double spatial_profile(const hamvec<3, double> &,
                                 const Param *) const;
  // anisotropy factor
  // check technical report for details
  // 1st argument: galactic centric Cartesian frame position
  // 2rd argument: parameter class object
  // 3th argument: regular magnetic field class object
  // 4th argument: regular magnetic field grid class object
  hamvec<3, double> anisotropy_direction(const hamvec<3, double> &,
                                         const Param *, const Breg *,
                                         const Grid_breg *) const;
  // anisotropy ratio
  // 1st argument: Galactic centric Cartesian frame position
  // 2rd argument: parameter class object
  // 3th argument: regular magnetic field class object
  // 4th argument: regular magnetic field grid class object
  double anisotropy_ratio(const hamvec<3, double> &, const Param *,
                          const Breg *, const Grid_breg *) const;
  // Gram-Schmidt orthogonalization process
  // 1st argument: wave-vector
  // 2nd arugment: input magnetic field vector (in Fourier space)
  // remark: real and imagine parts of complex magnetic field vector
  // in Fourier space handled separately
  hamvec<3, double> gramschmidt(const hamvec<3, double> &,
                                const hamvec<3, double> &) const;
};

// local anisotropic random magnetic field
// in compressive MHD plasma
class Brnd_mhd final : public Brnd_local {
public:
  Brnd_mhd() = default;
  Brnd_mhd(const Brnd_mhd &) = delete;
  Brnd_mhd(Brnd_mhd &&) = delete;
  Brnd_mhd &operator=(const Brnd_mhd &) = delete;
  Brnd_mhd &operator=(Brnd_mhd &&) = delete;
  virtual ~Brnd_mhd() = default;
  // use vector field decomposition scheme
  // use regular magnetic field at position of the Sun
  void write_grid(const Param *, const Breg *, const Grid_breg *,
                  Grid_brnd *) const override;
#ifndef NDEBUG
protected:
#endif
  // dynamo number
  // 1st argument: plasma beta
  // 2nd argument: cosine of k-B pitch angle
  inline double dynamo(const double &beta, const double &cosa) const {
    return 1 + 0.25 * beta * beta + beta * (1. - 2. * cosa * cosa);
  }
  // fast mode anisotropic tensor structure
  // 1st argument: plasma beta
  // 2nd argument: cosine of k-B pitch angle
  double h_f(const double &, const double &) const;
  // slow mode anisotropic tensor structure
  // 1st argument: plasma beta
  // 2nd argument: cosine of k-B pitch angle
  double h_s(const double &, const double &) const;
  // fast mode anisotropy power factor
  // 1st argument: plasma Mach number
  // 2nd argument: cosine of k-B pitch angle
  inline double F_a(const double &ma, const double &cosa) const {
    return std::exp(-1. * std::pow(ma, -1.33333333) * std::fabs(cosa) /
                    std::pow(1 - cosa * cosa, 0.33333333));
  }
  // slow mode anisotropy power factor
  // 1st argument: plasma Mach number
  // 2nd argument: cosine of k-B pitch angle
  inline double F_s(const double &ma, const double &cosa) const {
    return std::exp(-1. * std::pow(ma, -1.33333333) * std::fabs(cosa) /
                    std::pow(1 - cosa * cosa, 0.33333333));
  }
  // cosine of pitch angle between wavevector and regular GMF
  // 1st argument: field vector
  // 2nd argument: wave vector
  inline double cosine(const hamvec<3, double> &b,
                       const hamvec<3, double> &k) const {
    // dotprod is function from vec3.h
    return (b.versor()).dotprod(k.versor());
  }
  // direction of Alfven mode
  // 1st argument: regualr GMF vector
  // 2nd argument: wave vector
  hamvec<3, double> e_plus(const hamvec<3, double> &,
                           const hamvec<3, double> &) const;
  // direction of slow and fast modes
  // 1st argument: regualr GMF vector
  // 2nd argument: wave vector
  hamvec<3, double> e_minus(const hamvec<3, double> &,
                            const hamvec<3, double> &) const;
  // isotropic part of power spectrum of Alfvenic mode
  // 1st argument: wave vector magnitude
  // 2nd argument: parameter class object
  double spectrum_a(const double &, const Param *) const;
  // isotropic part of power spectrum of fast mode
  // 1st argument: wave vector magnitude
  // 2nd argument: parameter class object
  double spectrum_f(const double &, const Param *) const;
  // isotropic part of power spectrum of slow mode
  // 1st argument: wave vector magnitude
  // 2nd argument: parameter class object
  double spectrum_s(const double &, const Param *) const;
};

#endif
