///
/// free electron density distribution
///
#ifndef HAMMURABI_FE_H
#define HAMMURABI_FE_H

#include <iostream>
#include <vec3.h>
#include "pond.h"
#include "grid.h"

class FEreg{
public:
    FEreg(void) = default;
    virtual ~FEreg(void) = default;
    ///
    /// get free electron density, read from grid if granted, otherwise, calculate directly from density function
    ///
    virtual double get_density(const vec3_t<double> &,Pond *,Grid_fereg *);
    ///
    /// read from grid with trilinear interpolation
    ///
    virtual double read_grid(const vec3_t<double> &,Grid_fereg *);
    ///
    /// write to grid
    ///
    virtual void write_grid(Pond *,Grid_fereg *);
    ///
    /// calculate free electron density at given position
    ///
    virtual double density(const vec3_t<double> &,Pond *);
    ///
    /// gaussian blur free electron density (grid elemental scale as FWHM)
    ///
    virtual double density_blur(const vec3_t<double> &,Pond *,Grid_fereg *);
    
};

///
/// designed for verification
///
class FEreg_verify final : public FEreg{
public:
    FEreg_verify(void) = default;
    virtual ~FEreg_verify(void) = default;
    double density(const vec3_t<double> &, Pond *) override;
};

///
/// YMW16 modeling (ignore Fermi Bubble due to lack of observation)
///
class FEreg_ymw16 final : public FEreg{
public:
    FEreg_ymw16(void) = default;
    virtual ~FEreg_ymw16(void) = default;
    double density(const vec3_t<double> &, Pond *) override;
    
private:
    double thick(const double &,const double &,Pond *);
    double thin(const double &,const double &,Pond *);
    double spiral(const double &,const double &,const double &,const double &,Pond *);
    double galcen(const double &,const double &,const double&,Pond *);
    double gum(const double &,const double &,const double &,Pond *);
    double localbubble(const double &,const double &,const double &,const double &,const double &,Pond *);
    double nps(const double &,const double &,const double &,Pond *);
    // auxiliary functions
    inline double max(const double &a, const double &b){
        return (a>b) ? a : b;
    }
    inline double min(const double &a, const double &b){
        return (a<b) ? a : b;
    }
};

#endif

// END
