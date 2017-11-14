///
/// iso/aniso-tropic random free electron field generator
///
#ifndef HAMMURABI_FERND_H
#define HAMMURABI_FERND_H

#include <fftw3.h>
#include <vec3.h>
#include <vector>
#include <gsl/gsl_integration.h>
#include "pond.h"
#include "grid.h"
#include "cgs_units_file.h"
#include "fereg.h"

///
/// base class with read_grid implemented
///
class FErnd{
public:
    FErnd(void) = default;
    virtual ~FErnd(void) = default;
    ///
    /// return 0 if no specific field is instantiated
    ///
    virtual double get_fernd(const vec3_t<double> &,Grid_fernd *);
    virtual double read_grid(const vec3_t<double> &,Grid_fernd *);
    virtual void write_grid_iso(Pond *,Grid_fernd *);
    virtual void write_grid_ani(Pond *,FEreg *,Grid_fereg *,Grid_fernd *);
};

///
/// isotropic turbulent field
///
class FErnd_iso : public FErnd{
public:
    FErnd_iso(void) = default;
    virtual ~FErnd_iso(void) = default;
    double get_fernd(const vec3_t<double> &,Grid_fernd *) override;
    void write_grid_iso(Pond *,Grid_fernd *) override;
    ///
    /// isotropic turubulent power spectrum
    ///
    virtual double fe_spec(const double &,Pond *);
    ///
    /// calculate rescaling factor
    ///
    virtual double rescal_fact(const vec3_t<double> &,Pond *);
    
protected:
    void complex2real(const fftw_complex *,double *,const std::size_t &);
};

#endif

// END
