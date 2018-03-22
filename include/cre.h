///
/// synchrotron emissivity calculator
///
#ifndef HAMMURABI_CRE_H
#define HAMMURABI_CRE_H

#include <gsl/gsl_sf_synchrotron.h>
#include <param.h>
#include <grid.h>

///
/// base class, all functions are implemented in derived class
///
class CRE{
public:
    CRE(void) = default;
    virtual ~CRE(void) = default;
    virtual double get_emissivity_t(const vec3_t<double> &,Param *,Grid_cre *,const double &);
    virtual double get_emissivity_p(const vec3_t<double> &,Param *,Grid_cre *,const double &);
    ///
    /// read CRE flux from grid at given position
    /// (E_index, sylindrical_r, sylindrical_z) with {r,z} in cgs units,
    /// actual value of E is calculated from {E_index,Ek_min,Ek_fact}
    /// in get_emissivity automatically select bi/trilinear interpolation
    /// according to 2+1/3+1 spatial-spectral CRE flux grid
    ///
    virtual double read_grid(const std::size_t &, const vec3_t<double> &,Grid_cre *);
    virtual void write_grid(Param *,Grid_cre *);
    ///
    /// CRE flux at given CRE energy,
    /// input CRE energy at CGS units,
    /// output at [GeV m^2 s sr]^-1 units
    ///
    virtual double flux(const vec3_t<double> &,Param *,const double &);
};

///
/// designed for verification
///
class CRE_verify : public CRE{
public:
    CRE_verify(void) = default;
    virtual ~CRE_verify(void) = default;
    double get_emissivity_t(const vec3_t<double> &,Param *,Grid_cre *,const double &) override;
    double get_emissivity_p(const vec3_t<double> &,Param *,Grid_cre *,const double &) override;
    double flux(const vec3_t<double> &,Param *,const double &) override;
private:
    ///
    /// flux normalization at given position
    ///
    double flux_norm(const vec3_t<double> &,Param *);
    ///
    /// flux index at given position
    ///
    double flux_idx(const vec3_t<double> &,Param *);
    ///
    /// spatial CRE flux rescaling
    ///
    double rescal(const vec3_t<double> &,Param *);
};

///
/// Analytical CRE modeling
///
class CRE_ana : public CRE{
public:
    CRE_ana(void) = default;
    virtual ~CRE_ana(void) = default;
    double get_emissivity_t(const vec3_t<double> &,Param *,Grid_cre *,const double &) override;
    double get_emissivity_p(const vec3_t<double> &,Param *,Grid_cre *,const double &) override;
    double flux(const vec3_t<double> &,Param *,const double &) override;
private:
    ///
    /// flux normalization at given position
    ///
    double flux_norm(const vec3_t<double> &,Param *);
    ///
    /// flux index at given position
    ///
    double flux_idx(const vec3_t<double> &,Param *);
    ///
    /// spatial CRE flux rescaling
    ///
    double rescal(const vec3_t<double> &,Param *);
};

///
/// use Numerical CRE flux
///
class CRE_num final: public CRE{
public:
    CRE_num(void) = default;
    virtual ~CRE_num(void) = default;
    double get_emissivity_t(const vec3_t<double> &,Param *,Grid_cre *,const double &) override;
    double get_emissivity_p(const vec3_t<double> &,Param *,Grid_cre *,const double &) override;
};
#endif

// END
