// (an)isotropic random free electron field generator

#ifndef HAMMURABI_FERND_H
#define HAMMURABI_FERND_H

#include <hvec.h>

#include <fereg.h>
#include <grid.h>
#include <param.h>

// base class with read_grid implemented
class FErnd {
public:
  FErnd() = default;
  FErnd(const FErnd &) = delete;
  FErnd(FErnd &&) = delete;
  FErnd &operator=(const FErnd &) = delete;
  FErnd &operator=(FErnd &&) = delete;
  virtual ~FErnd() = default;
  // return 0 if no derived class is instantiated
  virtual double get_fernd(const hvec<3, double> &, const Param *,
                           const Grid_fernd *) const;
  virtual double read_grid(const hvec<3, double> &, const Param *,
                           const Grid_fernd *) const;
  virtual void write_grid(const Param *, Grid_fernd *) const;
};

// global random FE generator
// this class is treated as a covering class for specific methods
class FErnd_global : public FErnd {
public:
  FErnd_global() = default;
  FErnd_global(const FErnd_global &) = delete;
  FErnd_global(FErnd_global &&) = delete;
  FErnd_global &operator=(const FErnd_global &) = delete;
  FErnd_global &operator=(FErnd_global &&) = delete;
  virtual ~FErnd_global() = default;
};

// default method of global random FE
class FErnd_dft final : public FErnd_global {
public:
  FErnd_dft() = default;
  FErnd_dft(const FErnd_dft &) = delete;
  FErnd_dft(FErnd_dft &&) = delete;
  FErnd_dft &operator=(const FErnd_dft &) = delete;
  FErnd_dft &operator=(FErnd_dft &&) = delete;
  virtual ~FErnd_dft() = default;
  // trivial Fourier transform, with rescaling applied in spatial space
  void write_grid(const Param *, Grid_fernd *) const override;

protected:
  // isotropic turubulent power spectrum
  virtual double spec(const double &, const Param *) const;
  // density variance rescaling factor
  virtual double rescal(const hvec<3, double> &, const Param *) const;
};

#endif

// END
