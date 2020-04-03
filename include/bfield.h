// magnetic field generators/interface
//
// class design:
//
// Bfield
//
// method design:
//
//    (model collection)
//       |   |
//       |   | write_field
//       |   |
//       |   |--> (grid)
//       |          |
//       |          | read_field
//       |----------------------> (grid position/index)
//
//
// Bmodel
//   |-- Bmodel_xxx
//   |-- Bmodel_yyy
//
//
// method design:
//
//    (model) -----------> (grid)
//       |      write_grid    |
//       |                    |
//       |                    | read_grid
//       | read_model         |
//       ------------------------> (grid position/index)

#ifndef HAMMURABI_BFIELD_H
#define HAMMURABI_BFIELD_H

#include <cassert>
#include <grid.h>
#include <hamtype.h>
#include <hamunits.h>
#include <hamvec.h>
#include <param.h>
#include <stdexcept>
#include <string>

// magnetic field base class
class Bmodel {
public:
  Bmodel() = default;
  Bmodel(const Bmodel &) = delete;
  Bmodel(Bmodel &&) = delete;
  Bmodel &operator=(const Bmodel &) = delete;
  Bmodel &operator=(Bmodel &&) = delete;
  virtual ~Bmodel() = default;
  // assemble magnetic field
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter set
  virtual Hamvec<3, ham_float> read_model(const Hamvec<3, ham_float> &,
                                          const Param *) const;
  // read from field grid with given grid index
  // 1st argument: magnetic field grid index
  // 2nd argument: magnetic field grid class object
  virtual Hamvec<3, ham_float> read_grid(const ham_uint &,
                                         const Grid_b *) const;
  // read from field grid with linear interpolation
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter set
  // 3rd argument: magnetic field grid class object
  virtual Hamvec<3, ham_float> read_grid(const Hamvec<3, ham_float> &,
                                         const Param *, const Grid_b *) const;
  // write magnetic field to grid
  // 1st argument: parameter set
  // 2nd argument: magnetic field grid class object
  virtual void write_grid(const Param *, Grid_b *) const;
};

// uniform regular field modeling
// with fixed field orientation and strength
class Bmodel_unif final : public Bmodel {
public:
  Bmodel_unif() = default;
  Bmodel_unif(const Bmodel_unif &) = delete;
  Bmodel_unif(Bmodel_unif &&) = delete;
  Bmodel_unif &operator=(const Bmodel_unif &) = delete;
  Bmodel_unif &operator=(Bmodel_unif &&) = delete;
  virtual ~Bmodel_unif() = default;
  Hamvec<3, ham_float> read_model(const Hamvec<3, ham_float> &,
                                  const Param *) const override;
};

// WMAP-3yr LSA modeling
// http://iopscience.iop.org/article/10.1086/513699/meta
// with errata for GMF modeling
// https://lambda.gsfc.nasa.gov/product/map/dr2/pub_papers/threeyear/polarization/errata.cfm
class Bmodel_lsa final : public Bmodel {
public:
  Bmodel_lsa() = default;
  Bmodel_lsa(const Bmodel_lsa &) = delete;
  Bmodel_lsa(Bmodel_lsa &&) = delete;
  Bmodel_lsa &operator=(const Bmodel_lsa &) = delete;
  Bmodel_lsa &operator=(Bmodel_lsa &&) = delete;
  virtual ~Bmodel_lsa() = default;
  Hamvec<3, ham_float> read_model(const Hamvec<3, ham_float> &,
                                  const Param *) const override;
};

// Jaffe modeling
// https://academic.oup.com/mnras/article/401/2/1013/1150693
// https://www.aanda.org/articles/aa/abs/2016/12/aa28033-15/aa28033-15.html
class Bmodel_jaffe final : public Bmodel {
public:
  Bmodel_jaffe() = default;
  Bmodel_jaffe(const Bmodel_jaffe &) = delete;
  Bmodel_jaffe(Bmodel_jaffe &&) = delete;
  Bmodel_jaffe &operator=(const Bmodel_jaffe &) = delete;
  Bmodel_jaffe &operator=(Bmodel_jaffe &&) = delete;
  virtual ~Bmodel_jaffe() = default;
  Hamvec<3, ham_float> read_model(const Hamvec<3, ham_float> &,
                                  const Param *) const override;
#ifndef NDEBUG
protected:
#endif
  // field orientation
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter set
  Hamvec<3, ham_float> orientation(const Hamvec<3, ham_float> &,
                                   const Param *) const;
  // field amplitude radial scaling
  // 1st argument: Galactic centric Cartesian frame position
  // 2nd argument: parameter set
  ham_float radial_scaling(const Hamvec<3, ham_float> &, const Param *) const;
  // spiral arm height scaling
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter set
  // remark: inlined function
  inline ham_float arm_scaling(const Hamvec<3, ham_float> &pos,
                               const Param *par) const {
    return 1. / (cosh(pos[2] / par->bmodel_jaffe.arm_z0) *
                 cosh(pos[2] / par->bmodel_jaffe.arm_z0));
  }
  // disk height scaling
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter set
  // remark: inlined function
  inline ham_float disk_scaling(const Hamvec<3, ham_float> &pos,
                                const Param *par) const {
    return 1. / (cosh(pos[2] / par->bmodel_jaffe.disk_z0) *
                 cosh(pos[2] / par->bmodel_jaffe.disk_z0));
  }
  // halo height scaling
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter set
  // remark: inlined function
  inline ham_float halo_scaling(const Hamvec<3, ham_float> &pos,
                                const Param *par) const {
    return 1. / (cosh(pos[2] / par->bmodel_jaffe.halo_z0) *
                 cosh(pos[2] / par->bmodel_jaffe.halo_z0));
  }
  // spiral arm compression factor, for each arm
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter set
  std::vector<ham_float> arm_compress(const Hamvec<3, ham_float> &,
                                      const Param *) const;
  // spiral arm compression factor for dust, for each arm
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter set
  std::vector<ham_float> arm_compress_dust(const Hamvec<3, ham_float> &,
                                           const Param *) const;
  // distance to each spiral arm
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter set
  std::vector<ham_float> dist2arm(const Hamvec<3, ham_float> &,
                                  const Param *) const;
};

// Ensslin-Steininger method of global (an)isotropic random magnetic field
// https://iopscience.iop.org/article/10.3847/1538-4365/ab72a2
class Bmodel_es final : public Bmodel {
public:
  Bmodel_es() = default;
  Bmodel_es(const Bmodel_es &) = delete;
  Bmodel_es(Bmodel_es &&) = delete;
  Bmodel_es &operator=(const Bmodel_es &) = delete;
  Bmodel_es &operator=(Bmodel_es &&) = delete;
  virtual ~Bmodel_es() = default;
  void write_grid(const Param *, Grid_b *) const override;
#ifndef NDEBUG
protected:
#endif
  // isotropic power-spectrum
  // 1st argument: isotropic wave-vector magnitude
  // 2nd argument: parameter set
  ham_float spectrum(const ham_float &, const Param *) const;
  // field energy density reprofiling factor (analytic definition)
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter set
  ham_float spatial_profile(const Hamvec<3, ham_float> &, const Param *) const;
  // field energy density reprofiling factor (numeric definition)
  // 1st argument: magnetic field grid index
  // 2rd argument: parameter set
  // 3rd argument: magnetic field grid class object
  ham_float spatial_profile(const ham_uint &, const Param *,
                            const Grid_b *) const;
  // anisotropy factor (analytic definition)
  // check technical report for details
  // 1st argument: magnetic field grid index
  // 2rd argument: magnetic field grid class object
  Hamvec<3, ham_float> anisotropy_direction(const Hamvec<3, ham_float> &,
                                            const Param *) const;
  // anisotropy factor (numeric definition)
  // 1st argument: magnetic field grid index
  // 2rd argument: parameter set
  // 3rd argument: magnetic field grid class object
  Hamvec<3, ham_float> anisotropy_direction(const ham_uint &, const Param *,
                                            const Grid_b *) const;
  // anisotropy ratio (analytic definition)
  // 1st argument: Galactic centric Cartesian frame position
  // 2rd argument: parameter set
  // 3th argument: magnetic field class object
  // 4th argument: magnetic field grid class object
  ham_float anisotropy_ratio(const Hamvec<3, ham_float> &, const Param *) const;
  // anisotropy ratio (numeric definition)
  // 1st argument: magnetic field grid index
  // 2rd argument: parameter set
  // 3rd argument: magnetic field grid class object
  ham_float anisotropy_ratio(const ham_uint &, const Param *,
                             const Grid_b *) const;
  // Gram-Schmidt orthogonalization process
  // 1st argument: wave-vector
  // 2nd arugment: input magnetic field vector (in Fourier space)
  // remark: real and imagine parts of complex magnetic field vector
  // in Fourier space handled separately
  Hamvec<3, ham_float> gramschmidt(const Hamvec<3, ham_float> &,
                                   const Hamvec<3, ham_float> &) const;
};

// local anisotropic random magnetic field in compressible MHD plasma
// https://iopscience.iop.org/article/10.3847/1538-4365/ab72a2
class Bmodel_mhd final : public Bmodel {
public:
  Bmodel_mhd() = default;
  Bmodel_mhd(const Bmodel_mhd &) = delete;
  Bmodel_mhd(Bmodel_mhd &&) = delete;
  Bmodel_mhd &operator=(const Bmodel_mhd &) = delete;
  Bmodel_mhd &operator=(Bmodel_mhd &&) = delete;
  virtual ~Bmodel_mhd() = default;
  void write_grid(const Param *, Grid_b *) const override;
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
    return std::exp(-1. * std::pow(ma, -cgs::onethird) * std::fabs(cosa) / ma /
                    std::pow(1 - cosa * cosa, cgs::onethird));
  }
  // slow mode anisotropy power factor
  // 1st argument: plasma Mach number
  // 2nd argument: cosine of k-B pitch angle
  inline ham_float F_s(const ham_float &ma, const ham_float &cosa) const {
    return std::exp(-1. * std::pow(ma, -cgs::onethird) * std::fabs(cosa) / ma /
                    std::pow(1 - cosa * cosa, cgs::onethird));
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
  // 2nd argument: parameter set
  ham_float spectrum_a(const ham_float &, const Param *) const;
  // isotropic part of power spectrum of fast mode
  // 1st argument: wave vector magnitude
  // 2nd argument: parameter set
  ham_float spectrum_f(const ham_float &, const Param *) const;
  // isotropic part of power spectrum of slow mode
  // 1st argument: wave vector magnitude
  // 2nd argument: parameter set
  ham_float spectrum_s(const ham_float &, const Param *) const;
};

// magnetic field collection class
class Bfield {
public:
  Bfield() = default;
  Bfield(const Bfield &) = delete;
  Bfield(Bfield &&) = delete;
  Bfield &operator=(const Bfield &) = delete;
  Bfield &operator=(Bfield &&) = delete;
  virtual ~Bfield() = default;
  // initialize magnetic field model collection
  // model ordering arranged during parameter parsing
  Bfield(Param *par) {
    if (par->bmodel_list.empty()) {
      model_list.push_back(std::make_unique<Bmodel>());
    } else {
      for (auto &name : par->bmodel_list) {
        if (name == "unif") {
          model_list.push_back(std::make_unique<Bmodel_unif>());
        } else if (name == "lsa") {
          model_list.push_back(std::make_unique<Bmodel_lsa>());
        } else if (name == "jaffe") {
          model_list.push_back(std::make_unique<Bmodel_jaffe>());
        } else if (name == "es") {
          model_list.push_back(std::make_unique<Bmodel_es>());
        } else if (name == "mhd") {
          model_list.push_back(std::make_unique<Bmodel_mhd>());
        } else {
          throw std::runtime_error("unknown B model");
        }
      }
    }
  }
  // read magnetic field at given grid index
  // 1st argument: magnetic field grid index
  // 2nd argument: magnetic field grid
  Hamvec<3, ham_float> read_field(const ham_uint &idx,
                                  const Grid_b *grid) const {
    return Hamvec<3, ham_float>(grid->bx[idx], grid->by[idx], grid->bz[idx]);
  }
  // read magnetic field at given position (with interpolation)
  // 1st argument: Cartesian frame Galactic centric position
  // 2nd argument: parameter set
  // 34d argument: magnetic field grid
  Hamvec<3, ham_float> read_field(const Hamvec<3, ham_float> &pos,
                                  const Param *par, const Grid_b *grid) const {
    Hamvec<3, ham_float> tmp(0., 0., 0.);
    // get field content from grid
    // model ordering arranged during parameter parsing
    if (par->grid_b.build_permission or par->grid_b.read_permission) {
      // read_grid defined in base class
      tmp = model_list[0]->read_grid(pos, par, grid);
    }
    // get field content directly from model
    // model ordering arranged during parameter parsing
    else {
      for (auto &model : model_list) {
        tmp += model->read_model(pos, par);
      }
    }
    return tmp;
  }
  // write magnetic field to grid before reading
  // 1st argument: parameter set
  // 2nd argument: magnetic field grid
  void write_field(const Param *par, Grid_b *grid) {
    // write field content to grid
    // model ordering arranged during parameter parsing
    assert(par->grid_b.build_permission or par->grid_b.write_permission or
           par->grid_b.read_permission);
    for (auto &model : model_list) {
      model->write_grid(par, grid);
    }
  }

protected:
  // list of magnetic field model pointers
  std::vector<std::unique_ptr<Bmodel>> model_list;
};

#endif
