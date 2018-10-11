// (an)isotropic random free electron field generator

#ifndef HAMMURABI_FERND_H
#define HAMMURABI_FERND_H

#include <vec3.h>

#include <param.h>
#include <grid.h>
#include <fereg.h>

// base class with read_grid implemented
class FErnd{
    
public:
    
    FErnd () = default;
    
    virtual ~FErnd () = default;
    
    // return 0 if no derived class is instantiated
    virtual double get_fernd (const vec3_t<double> &,
                              const Grid_fernd *) const;
    
    virtual double read_grid (const vec3_t<double> &,
                              const Grid_fernd *) const;
    
    virtual void write_grid (const Param *,
                             Grid_fernd *) const;
};

// global random FE generator
// this class is treated as a covering class for specific methods
class FErnd_global : public FErnd{
    
public:
    
    FErnd_global () = default;
    
    virtual ~FErnd_global () = default;
};

// default method of global random FE
class FErnd_dft final : public FErnd_global{
    
public:
    
    FErnd_dft () = default;
    
    virtual ~FErnd_dft () = default;
    
    // trivial Fourier transform, with rescaling applied in spatial space
    void write_grid (const Param *,
                     Grid_fernd *) const override;
    
protected:
    
    // isotropic turubulent power spectrum
    virtual double spec (const double &,
                         const Param *) const;
    
    // density variance rescaling factor
    virtual double rescal (const vec3_t<double> &,
                           const Param *) const;
};

#endif

// END
