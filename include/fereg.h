///
/// free electron density distribution
///
#ifndef HAMMURABI_FE_H
#define HAMMURABI_FE_H

#include <iostream>
#include <vec3.h>
#include <param.h>
#include <grid.h>

class FEreg{
public:
    FEreg(void) = default;
    virtual ~FEreg(void) = default;
    ///
    /// get free electron density,
    /// read from grid if granted, otherwise,
    /// calculate directly from density function
    ///
    virtual double get_density(const vec3_t<double> &,Param *,Grid_fereg *);
    ///
    /// read from grid with trilinear interpolation
    ///
    virtual double read_grid(const vec3_t<double> &,Grid_fereg *);
    ///
    /// write to grid
    ///
    virtual void write_grid(Param *,Grid_fereg *);
    ///
    /// assemble free electron density at given position
    ///
    virtual double density(const vec3_t<double> &,Param *);
    ///
    /// gaussian blur free electron density (grid elemental scale as FWHM),
    /// computing time consuming function based on \p density
    ///
    virtual double density_blur(const vec3_t<double> &,Param *,Grid_fereg *);
    
};

///
/// designed for verification
///
class FEreg_verify final : public FEreg{
public:
    FEreg_verify(void) = default;
    virtual ~FEreg_verify(void) = default;
    double density(const vec3_t<double> &, Param *) override;
};

///
/// YMW16 modeling (ignore Fermi Bubble due to lack of observation)
///
class FEreg_ymw16 final : public FEreg{
public:
    FEreg_ymw16(void) = default;
    virtual ~FEreg_ymw16(void) = default;
    double density(const vec3_t<double> &, Param *) override;
    
private:
    ///
    /// thick disk
    ///
    double thick(const double &,const double &,Param *);
    ///
    /// thin disk
    ///
    double thin(const double &,const double &,Param *);
    ///
    /// spiral arms
    ///
    double spiral(const double &,const double &,const double &,const double &,Param *);
    ///
    /// galactic center
    ///
    double galcen(const double &,const double &,const double&,Param *);
    ///
    /// gum nebula
    ///
    double gum(const double &,const double &,const double &,Param *);
    ///
    /// local bubble
    ///
    double localbubble(const double &,const double &,const double &,const double &,const double &,Param *);
    ///
    /// northern polar spurs
    ///
    double nps(const double &,const double &,const double &,Param *);
};

#endif

// END
