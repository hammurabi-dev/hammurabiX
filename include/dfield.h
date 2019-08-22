// dust density field generators/interface
//
// class design:
//
// Dfield
//    |-- Dreg
//    |     |-- Dreg_xxx
//    |     |-- Dreg_yyy
//    |
//    |-- Drnd
//          |-- Drnd_xxx
//          |-- Drnd_yyy
//
// Dfield has similar structure and member functions as Dfield

#ifndef HAMMURABI_D_H
#define HAMMURABI_D_H

#include <grid.h>
#include <hamvec.h>
#include <param.h>

class Dfield {
public:
  Dfield() = default;
  Dfield(const Dfield &) = delete;
  Dfield(Dfield &&) = delete;
  Dfield &operator=(const Dfield &) = delete;
  Dfield &operator=(Dfield &&) = delete;
  virtual ~Dfield() = default;
};

class Dreg : public Dfield {
public:
  Dreg() = default;
  Dreg(const Dreg &) = delete;
  Dreg(Dreg &&) = delete;
  Dreg &operator=(const Dreg &) = delete;
  Dreg &operator=(Dreg &&) = delete;
  virtual ~Dreg() = default;
  // get dust density
  // read from grid if granted, otherwise
  // calculate directly from density function
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  // 3rd argument: electron grid class object
  virtual double read_field(const hamvec<3, double> &, const Param *,
                            const Grid_dreg *) const;
  // assemble dust density at given position
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter class object
  virtual double write_field(const hamvec<3, double> &, const Param *) const;
  // read from grid with linear interpolation
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: electron field grid class object
  virtual double read_grid(const hamvec<3, double> &, const Param *,
                           const Grid_dreg *) const;
  // write to grid
  // 1st argument: parameter class object
  // 2nd argument: electron field grid class object
  virtual void write_grid(const Param *, Grid_dreg *) const;
};

// base class with read_grid implemented
class Drnd : public Dfield {
public:
  Drnd() = default;
  Drnd(const Drnd &) = delete;
  Drnd(Drnd &&) = delete;
  Drnd &operator=(const Drnd &) = delete;
  Drnd &operator=(Drnd &&) = delete;
  virtual ~Drnd() = default;
  virtual double read_field(const hamvec<3, double> &, const Param *,
                            const Grid_drnd *) const;
  virtual double read_grid(const hamvec<3, double> &, const Param *,
                           const Grid_drnd *) const;
  virtual void write_grid(const Param *, Grid_drnd *) const;
};

//--------------------------- Dreg DERIVED ----------------------------------//

// uniform field modeling
class Dreg_unif final : public Dreg {
public:
  Dreg_unif() = default;
  Dreg_unif(const Dreg_unif &) = delete;
  Dreg_unif(Dreg_unif &&) = delete;
  Dreg_unif &operator=(const Dreg_unif &) = delete;
  Dreg_unif &operator=(Dreg_unif &&) = delete;
  virtual ~Dreg_unif() = default;
  double write_field(const hamvec<3, double> &, const Param *) const override;
};

//--------------------------- Drnd DERIVED ----------------------------------//

// global random dust generator
// this class is treated as a covering class for specific methods
class Drnd_global : public Drnd {
public:
  Drnd_global() = default;
  Drnd_global(const Drnd_global &) = delete;
  Drnd_global(Drnd_global &&) = delete;
  Drnd_global &operator=(const Drnd_global &) = delete;
  Drnd_global &operator=(Drnd_global &&) = delete;
  virtual ~Drnd_global() = default;
};

// default method of global random dust
class Drnd_dft final : public Drnd_global {
public:
  Drnd_dft() = default;
  Drnd_dft(const Drnd_dft &) = delete;
  Drnd_dft(Drnd_dft &&) = delete;
  Drnd_dft &operator=(const Drnd_dft &) = delete;
  Drnd_dft &operator=(Drnd_dft &&) = delete;
  virtual ~Drnd_dft() = default;
  // trivial Fourier transform, with rescaling applied in spatial space
  void write_grid(const Param *, Grid_drnd *) const override;

protected:
  // isotropic turubulent power spectrum
  virtual double spectrum(const double &, const Param *) const;
  // density variance rescaling factor
  virtual double spatial_profile(const hamvec<3, double> &,
                                 const Param *) const;
};

#endif
