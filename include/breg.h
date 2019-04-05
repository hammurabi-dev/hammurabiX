// regular GMF generators
// field parameterized modelling is implemented in breg function
// interactions with grid are carried out via read_grid and write_grid
// notice that read/write_grid acts with internal grid in memory
// interactions with disk data are carried out in Grid_breg class
// get_breg is the universal interface for retrieving field vector
// at given 3D spatial position
// get_breg knows where to read the information
// either from grid or parameterzed modelling
//
// customized field modelling can be added as derived class of Breg

#ifndef HAMMURABI_BREG_H
#define HAMMURABI_BREG_H

#include <grid.h>
#include <hvec.h>
#include <param.h>
#include <vector>

// base class of GMF generator,
// read_grid and write_grid are implemented here
class Breg {
public:
  Breg() = default;
  Breg(const Breg &) = delete;
  Breg(Breg &&) = delete;
  Breg &operator=(const Breg &) = delete;
  Breg &operator=(Breg &&) = delete;
  virtual ~Breg() = default;
  // fetch regular magnetic field
  // 1st argument: Galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  // 3rd argument: regular GMF grid class object
  // inovke read_grid regardless of field type if permitted, otherwise invoke
  // breg
  virtual hvec<3, double> get_breg(const hvec<3, double> &, const Param *,
                                   const Grid_breg *) const;
  // GMF assembler, specified only in derived class
  // 1st argument: Galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  virtual hvec<3, double> breg(const hvec<3, double> &, const Param *) const;
  // read from field grid with trilinear interpolation
  // 1st argument: Galactic centric Cartesian frame position
  // 2nd argument: regular GMF grid class object
  virtual hvec<3, double> read_grid(const hvec<3, double> &, const Param *,
                                    const Grid_breg *) const;
  // write to field grid
  // 1st argument: parameter class object
  // 2nd argument: regular GMF grid class object
  virtual void write_grid(const Param *, Grid_breg *) const;
};

// uniform field modeling
// with fixed field orientation and strength
class Breg_unif final : public Breg {
public:
  Breg_unif() = default;
  Breg_unif(const Breg_unif &) = delete;
  Breg_unif(Breg_unif &&) = delete;
  Breg_unif &operator=(const Breg_unif &) = delete;
  Breg_unif &operator=(Breg_unif &&) = delete;
  virtual ~Breg_unif() = default;
  hvec<3, double> breg(const hvec<3, double> &, const Param *) const override;
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
  hvec<3, double> breg(const hvec<3, double> &, const Param *) const override;
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
  hvec<3, double> breg(const hvec<3, double> &, const Param *) const override;
#ifndef NDEBUG
protected:
#endif
  // field orientation
  // 1st argument: Galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  hvec<3, double> orientation(const hvec<3, double> &, const Param *) const;
  // field amplitude radial scaling
  // 1st argument: Galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  double radial_scaling(const hvec<3, double> &, const Param *) const;
  // spiral arm height scaling
  // 1st argument: Galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  // remark: inlined function
  inline double arm_scaling(const hvec<3, double> &pos,
                            const Param *par) const {
    return 1. / (cosh(pos[2] / par->breg_jaffe.arm_z0) *
                 cosh(pos[2] / par->breg_jaffe.arm_z0));
  }
  // disk height scaling
  // 1st argument: Galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  // remark: inlined function
  inline double disk_scaling(const hvec<3, double> &pos,
                             const Param *par) const {
    return 1. / (cosh(pos[2] / par->breg_jaffe.disk_z0) *
                 cosh(pos[2] / par->breg_jaffe.disk_z0));
  }
  // halo height scaling
  // 1st argument: Galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  // remark: inlined function
  inline double halo_scaling(const hvec<3, double> &pos,
                             const Param *par) const {
    return 1. / (cosh(pos[2] / par->breg_jaffe.halo_z0) *
                 cosh(pos[2] / par->breg_jaffe.halo_z0));
  }
  // spiral arm compression factor, for each arm
  // 1st argument: Galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  std::vector<double> arm_compress(const hvec<3, double> &,
                                   const Param *) const;
  // spiral arm compression factor for dust, for each arm
  // 1st argument: Galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  std::vector<double> arm_compress_dust(const hvec<3, double> &,
                                        const Param *) const;
  // distance to each spiral arm
  // 1st argument: Galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  std::vector<double> dist2arm(const hvec<3, double> &, const Param *) const;
};

#endif

// END
