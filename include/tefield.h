// thermal electron field generators/interface
//
// class design:
//
// TEfield
//    |-- TEreg
//    |     |-- TEreg_xxx
//    |     |-- TEreg_yyy
//    |
//    |-- TErnd
//          |-- TErnd_xxx
//          |-- TErnd_yyy
//
// TEfield has similar structure and member functions as Bfield
// the major difference is thermal electron density is scalar field

#ifndef HAMMURABI_TE_H
#define HAMMURABI_TE_H

#include <grid.h>
#include <hamvec.h>
#include <param.h>

class TEfield {
public:
  TEfield() = default;
  TEfield(const TEfield &) = delete;
  TEfield(TEfield &&) = delete;
  TEfield &operator=(const TEfield &) = delete;
  TEfield &operator=(TEfield &&) = delete;
  virtual ~TEfield() = default;
};

class TEreg : public TEfield {
public:
  TEreg() = default;
  TEreg(const TEreg &) = delete;
  TEreg(TEreg &&) = delete;
  TEreg &operator=(const TEreg &) = delete;
  TEreg &operator=(TEreg &&) = delete;
  virtual ~TEreg() = default;
  // get thermal electron density
  // read from grid if granted, otherwise
  // calculate directly from density function
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  // 3rd argument: electron grid class object
  virtual double read_field(const hamvec<3, double> &, const Param *,
                            const Grid_tereg *) const;
  // assemble thermal electron density at given position
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  virtual double write_field(const hamvec<3, double> &, const Param *) const;
  // read from grid with linear interpolation
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: electron field grid class object
  virtual double read_grid(const hamvec<3, double> &, const Param *,
                           const Grid_tereg *) const;
  // write to grid
  // 1st argument: parameter class object
  // 2nd argument: electron field grid class object
  virtual void write_grid(const Param *, Grid_tereg *) const;
};

// base class with read_grid implemented
class TErnd : public TEfield {
public:
  TErnd() = default;
  TErnd(const TErnd &) = delete;
  TErnd(TErnd &&) = delete;
  TErnd &operator=(const TErnd &) = delete;
  TErnd &operator=(TErnd &&) = delete;
  virtual ~TErnd() = default;
  virtual double read_field(const hamvec<3, double> &, const Param *,
                            const Grid_ternd *) const;
  virtual double read_grid(const hamvec<3, double> &, const Param *,
                           const Grid_ternd *) const;
  virtual void write_grid(const Param *, const TEreg *, const Grid_tereg *,
                          Grid_ternd *) const;
};

//--------------------------- TEreg DERIVED ----------------------------------//

// uniform field modeling
class TEreg_unif final : public TEreg {
public:
  TEreg_unif() = default;
  TEreg_unif(const TEreg_unif &) = delete;
  TEreg_unif(TEreg_unif &&) = delete;
  TEreg_unif &operator=(const TEreg_unif &) = delete;
  TEreg_unif &operator=(TEreg_unif &&) = delete;
  virtual ~TEreg_unif() = default;
  double write_field(const hamvec<3, double> &, const Param *) const override;
};

// YMW16 modeling (ignore Fermi Bubble due to lack of observation)
class TEreg_ymw16 final : public TEreg {
public:
  TEreg_ymw16() = default;
  TEreg_ymw16(const TEreg_ymw16 &) = delete;
  TEreg_ymw16(TEreg_ymw16 &&) = delete;
  TEreg_ymw16 &operator=(const TEreg_ymw16 &) = delete;
  TEreg_ymw16 &operator=(TEreg_ymw16 &&) = delete;
  virtual ~TEreg_ymw16() = default;
  double write_field(const hamvec<3, double> &, const Param *) const override;
#ifdef NDEBUG
private:
#endif
  // thick disk
  double thick(const double &, const double &, const Param *) const;
  // thin disk
  double thin(const double &, const double &, const Param *) const;
  // spiral arms
  double spiral(const double &, const double &, const double &, const double &,
                const Param *) const;
  // galactic center
  double galcen(const double &, const double &, const double &,
                const Param *) const;
  // gum nebula
  double gum(const double &, const double &, const double &,
             const Param *) const;
  // local bubble
  double localbubble(const double &, const double &, const double &,
                     const double &, const double &, const Param *) const;
  // northern polar spurs
  double nps(const double &, const double &, const double &,
             const Param *) const;
};

//--------------------------- TErnd DERIVED ----------------------------------//

// global random FE generator
// this class is treated as a covering class for specific methods
class TErnd_global : public TErnd {
public:
  TErnd_global() = default;
  TErnd_global(const TErnd_global &) = delete;
  TErnd_global(TErnd_global &&) = delete;
  TErnd_global &operator=(const TErnd_global &) = delete;
  TErnd_global &operator=(TErnd_global &&) = delete;
  virtual ~TErnd_global() = default;
};

// default method of global random FE
class TErnd_dft final : public TErnd_global {
public:
  TErnd_dft() = default;
  TErnd_dft(const TErnd_dft &) = delete;
  TErnd_dft(TErnd_dft &&) = delete;
  TErnd_dft &operator=(const TErnd_dft &) = delete;
  TErnd_dft &operator=(TErnd_dft &&) = delete;
  virtual ~TErnd_dft() = default;
  // trivial Fourier transform, with rescaling applied in spatial space
  void write_grid(const Param *, const TEreg *, const Grid_tereg *,
                  Grid_ternd *) const override;

protected:
  // isotropic turubulent power spectrum
  virtual double spectrum(const double &, const Param *) const;
  // density variance rescaling factor
  virtual double spatial_profile(const hamvec<3, double> &,
                                 const Param *) const;
};

#endif
