///
/// iso/aniso-tropic random free electron field generator
///
#ifndef HAMMURABI_FERND_H
#define HAMMURABI_FERND_H

#include <fftw3.h>
#include <vec3.h>
#include <vector>
#include <gsl/gsl_integration.h>
#include <param.h>
#include <grid.h>
#include <cgs_units_file.h>
#include <fereg.h>

///
/// base class with read_grid implemented
///
class FErnd{
public:
    FErnd(void) = default;
    virtual ~FErnd(void) = default;
    ///
    /// return 0 if no derived class is instantiated
    ///
    virtual double get_fernd(const vec3_t<double> &,Grid_fernd *);
    virtual double read_grid(const vec3_t<double> &,Grid_fernd *);
    virtual void write_grid(Param *,Grid_fernd *);
};

///
/// global turbulent field
///
class FErnd_global : public FErnd{
public:
    FErnd_global(void) = default;
    virtual ~FErnd_global(void) = default;
    double get_fernd(const vec3_t<double> &,Grid_fernd *) override;
    ///
    /// trivial Fourier transform, with rescaling applied in spatial space
    ///
    void write_grid(Param *,Grid_fernd *) override;
    
protected:
    ///
    /// isotropic turubulent power spectrum
    ///
    virtual double spec(const double &,Param *);
    ///
    /// density variance rescaling factor
    ///
    virtual double rescal(const vec3_t<double> &,Param *);
};

#endif

// END
