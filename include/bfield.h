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
#include <hamtype.h>
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
  virtual Hamvec<3, ham_float> read_field(const Hamvec<3, ham_float> &,
                                          const Param *,
                                          const Grid_breg *) const;
  // assemble analytic magnetic field, specified only in derived class
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  virtual Hamvec<3, ham_float> write_field(const Hamvec<3, ham_float> &,
                                           const Param *) const;
  // read from field grid with linear interpolation
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  // 3rd argument: magnetic field grid class object
  virtual Hamvec<3, ham_float> read_grid(const Hamvec<3, ham_float> &,
                                         const Param *,
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
  virtual Hamvec<3, ham_float> read_field(const Hamvec<3, ham_float> &,
                                          const Param *,
                                          const Grid_brnd *) const;
  // read from field grid with linear interpolation
  virtual Hamvec<3, ham_float> read_grid(const Hamvec<3, ham_float> &,
                                         const Param *,
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
  Hamvec<3, ham_float> write_field(const Hamvec<3, ham_float> &,
                                   const Param *) const override;
};

// uniform Cartesian field
// with fixed field orientation and strength
class Breg_cart final : public Breg {
public:
  Breg_cart() = default;
  Breg_cart(const Breg_cart &) = delete;
  Breg_cart(Breg_cart &&) = delete;
  Breg_cart &operator=(const Breg_cart &) = delete;
  Breg_cart &operator=(Breg_cart &&) = delete;
  virtual ~Breg_cart() = default;
  Hamvec<3, ham_float> write_field(const Hamvec<3, ham_float> &,
                                   const Param *) const override;
};

// Helical field
class Breg_helix final : public Breg {
public:
  Breg_helix() = default;
  Breg_helix(const Breg_helix &) = delete;
  Breg_helix(Breg_helix &&) = delete;
  Breg_helix &operator=(const Breg_helix &) = delete;
  Breg_helix &operator=(Breg_helix &&) = delete;
  virtual ~Breg_helix() = default;
  Hamvec<3, ham_float> write_field(const Hamvec<3, ham_float> &,
                                   const Param *) const override;
};

// Jansson & Farrer 2012 model
/// http://iopscience.iop.org/article/10.1088/0004-637X/757/1/14/meta
class Breg_jf12 final : public Breg {
public:
  Breg_jf12() = default;
  Breg_jf12(const Breg_jf12 &) = delete;
  Breg_jf12(Breg_jf12 &&) = delete;
  Breg_jf12 &operator=(const Breg_jf12 &) = delete;
  Breg_jf12 &operator=(Breg_jf12 &&) = delete;
  virtual ~Breg_jf12() = default;
  Hamvec<3, ham_float> write_field(const Hamvec<3, ham_float> &,
                                   const Param *) const override;
};

// WMAP-3yr LSA modeling
// http://iopscience.iop.org/article/10.1086/513699/meta
// with errata for GMF modeling
// https://lambda.gsfc.nasa.gov/product/map/dr2/pub_papers/threeyear/polarization/errata.cfm
class Breg_lsa final : public Breg {
public:
  Breg_lsa() = default;
  Breg_lsa(const Breg_lsa &) = delete;
  Breg_lsa(Breg_lsa &&) = delete;
  Breg_lsa &operator=(const Breg_lsa &) = delete;
  Breg_lsa &operator=(Breg_lsa &&) = delete;
  virtual ~Breg_lsa() = default;
  Hamvec<3, ham_float> write_field(const Hamvec<3, ham_float> &,
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
  Hamvec<3, ham_float> write_field(const Hamvec<3, ham_float> &,
                                   const Param *) const override;
#ifndef NDEBUG
protected:
#endif
  // field orientation
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  Hamvec<3, ham_float> orientation(const Hamvec<3, ham_float> &,
                                   const Param *) const;
  // field amplitude radial scaling
  // 1st argument: Galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  ham_float radial_scaling(const Hamvec<3, ham_float> &, const Param *) const;
  // spiral arm height scaling
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  // remark: inlined function
  inline ham_float arm_scaling(const Hamvec<3, ham_float> &pos,
                               const Param *par) const {
    return 1. / (cosh(pos[2] / par->breg_jaffe.arm_z0) *
                 cosh(pos[2] / par->breg_jaffe.arm_z0));
  }
  // disk height scaling
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  // remark: inlined function
  inline ham_float disk_scaling(const Hamvec<3, ham_float> &pos,
                                const Param *par) const {
    return 1. / (cosh(pos[2] / par->breg_jaffe.disk_z0) *
                 cosh(pos[2] / par->breg_jaffe.disk_z0));
  }
  // halo height scaling
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  // remark: inlined function
  inline ham_float halo_scaling(const Hamvec<3, ham_float> &pos,
                                const Param *par) const {
    return 1. / (cosh(pos[2] / par->breg_jaffe.halo_z0) *
                 cosh(pos[2] / par->breg_jaffe.halo_z0));
  }
  // spiral arm compression factor, for each arm
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  std::vector<ham_float> arm_compress(const Hamvec<3, ham_float> &,
                                      const Param *) const;
  // spiral arm compression factor for dust, for each arm
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  std::vector<ham_float> arm_compress_dust(const Hamvec<3, ham_float> &,
                                           const Param *) const;
  // distance to each spiral arm
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  std::vector<ham_float> dist2arm(const Hamvec<3, ham_float> &,
                                  const Param *) const;
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
  virtual ham_float spectrum(const ham_float &, const Param *) const;
  // field energy density reprofiling factor
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  virtual ham_float spatial_profile(const Hamvec<3, ham_float> &,
                                    const Param *) const;
  // anisotropy factor
  // check technical report for details
  // 1st argument: galactic centric Cartesian frame position
  // 2rd argument: parameter class object
  // 3th argument: regular magnetic field class object
  // 4th argument: regular magnetic field grid class object
  Hamvec<3, ham_float> anisotropy_direction(const Hamvec<3, ham_float> &,
                                            const Param *, const Breg *,
                                            const Grid_breg *) const;
  // anisotropy ratio
  // 1st argument: Galactic centric Cartesian frame position
  // 2rd argument: parameter class object
  // 3th argument: regular magnetic field class object
  // 4th argument: regular magnetic field grid class object
  ham_float anisotropy_ratio(const Hamvec<3, ham_float> &, const Param *,
                             const Breg *, const Grid_breg *) const;
  // Gram-Schmidt orthogonalization process
  // 1st argument: wave-vector
  // 2nd arugment: input magnetic field vector (in Fourier space)
  // remark: real and imagine parts of complex magnetic field vector
  // in Fourier space handled separately
  Hamvec<3, ham_float> gramschmidt(const Hamvec<3, ham_float> &,
                                   const Hamvec<3, ham_float> &) const;
};

// Jansson-Farrar 2012 method of global (an)isotropic random magnetic field
class Brnd_jf12 final : public Brnd_global {
public:
  Brnd_jf12() = default;
  Brnd_jf12(const Brnd_jf12 &) = delete;
  Brnd_jf12(Brnd_jf12 &&) = delete;
  Brnd_jf12 &operator=(const Brnd_jf12 &) = delete;
  Brnd_jf12 &operator=(Brnd_jf12 &&) = delete;
  virtual ~Brnd_jf12() = default;
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
  virtual ham_float spectrum(const ham_float &, const Param *) const;
  // field energy density reprofiling factor
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  virtual ham_float spatial_profile(const Hamvec<3, ham_float> &,
                                    const Param *) const;
  // anisotropy factor
  // check technical report for details
  // 1st argument: galactic centric Cartesian frame position
  // 2rd argument: parameter class object
  // 3th argument: regular magnetic field class object
  // 4th argument: regular magnetic field grid class object
  Hamvec<3, ham_float> anisotropy_direction(const Hamvec<3, ham_float> &,
                                            const Param *, const Breg *,
                                            const Grid_breg *) const;
  // anisotropy ratio
  // 1st argument: Galactic centric Cartesian frame position
  // 2rd argument: parameter class object
  // 3th argument: regular magnetic field class object
  // 4th argument: regular magnetic field grid class object
  ham_float anisotropy_ratio(const Hamvec<3, ham_float> &, const Param *,
                             const Breg *, const Grid_breg *) const;
  // Gram-Schmidt orthogonalization process
  // 1st argument: wave-vector
  // 2nd arugment: input magnetic field vector (in Fourier space)
  // remark: real and imagine parts of complex magnetic field vector
  // in Fourier space handled separately
  Hamvec<3, ham_float> gramschmidt(const Hamvec<3, ham_float> &,
                                   const Hamvec<3, ham_float> &) const;
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
  inline ham_float dynamo(const ham_float &beta, const ham_float &cosa) const {
    return 1 + 0.25 * beta * beta + beta * (1. - 2. * cosa * cosa);
  }
  // fast mode anisotropic tensor structure
  // 1st argument: plasma beta
  // 2nd argument: cosine of k-B pitch angle
  ham_float h_f(const ham_float &, const ham_float &) const;
  // slow mode anisotropic tensor structure
  // 1st argument: plasma beta
  // 2nd argument: cosine of k-B pitch angle
  ham_float h_s(const ham_float &, const ham_float &) const;
  // fast mode anisotropy power factor
  // 1st argument: plasma Mach number
  // 2nd argument: cosine of k-B pitch angle
  inline ham_float F_a(const ham_float &ma, const ham_float &cosa) const {
    return std::exp(-1. * std::pow(ma, -1.33333333) * std::fabs(cosa) /
                    std::pow(1 - cosa * cosa, 0.33333333));
  }
  // slow mode anisotropy power factor
  // 1st argument: plasma Mach number
  // 2nd argument: cosine of k-B pitch angle
  inline ham_float F_s(const ham_float &ma, const ham_float &cosa) const {
    return std::exp(-1. * std::pow(ma, -1.33333333) * std::fabs(cosa) /
                    std::pow(1 - cosa * cosa, 0.33333333));
  }
  // cosine of pitch angle between wavevector and regular GMF
  // 1st argument: field vector
  // 2nd argument: wave vector
  inline ham_float cosine(const Hamvec<3, ham_float> &b,
                          const Hamvec<3, ham_float> &k) const {
    // dotprod is function from vec3.h
    return (b.versor()).dotprod(k.versor());
  }
  // direction of Alfven mode
  // 1st argument: regualr GMF vector
  // 2nd argument: wave vector
  Hamvec<3, ham_float> e_plus(const Hamvec<3, ham_float> &,
                              const Hamvec<3, ham_float> &) const;
  // direction of slow and fast modes
  // 1st argument: regualr GMF vector
  // 2nd argument: wave vector
  Hamvec<3, ham_float> e_minus(const Hamvec<3, ham_float> &,
                               const Hamvec<3, ham_float> &) const;
  // isotropic part of power spectrum of Alfvenic mode
  // 1st argument: wave vector magnitude
  // 2nd argument: parameter class object
  ham_float spectrum_a(const ham_float &, const Param *) const;
  // isotropic part of power spectrum of fast mode
  // 1st argument: wave vector magnitude
  // 2nd argument: parameter class object
  ham_float spectrum_f(const ham_float &, const Param *) const;
  // isotropic part of power spectrum of slow mode
  // 1st argument: wave vector magnitude
  // 2nd argument: parameter class object
  ham_float spectrum_s(const ham_float &, const Param *) const;
};

#endif
