/*
 *@file: cre.h
 *@brief: synchrotron emissivity calculator
 */
#ifndef HAMMURABI_CRE_H
#define HAMMURABI_CRE_H

#include <gsl/gsl_sf_synchrotron.h>
#include "pond.h"
#include "grid.h"

class CRE{
public:
    CRE(void) = default;
    virtual ~CRE(void) = default;
    /*@get_emissivity
     * get J_tot or J_pol at given position
     */
    virtual double get_emissivity_t(const vec3_t<double> &,Pond *,Grid_cre *,const double &);
    virtual double get_emissivity_p(const vec3_t<double> &,Pond *,Grid_cre *,const double &);
    virtual double read_grid(const unsigned int &, const vec3_t<double> &,Grid_cre *);
    virtual void write_grid(Pond *,Grid_cre *);
    
};

// verification
class CRE_verify : public CRE{
public:
    CRE_verify(void) = default;
    virtual ~CRE_verify(void) = default;
    double get_emissivity_t(const vec3_t<double> &,Pond *,Grid_cre *,const double &) override;
    double get_emissivity_p(const vec3_t<double> &,Pond *,Grid_cre *,const double &) override;
private:
    double flux(const vec3_t<double> &,Pond *,const double &);
    void flux_param(const vec3_t<double> &,Pond *,double &,double &);
    void write_grid(Pond *,Grid_cre *) override;
};

// Analytical CRE
class CRE_ana : public CRE{
public:
    CRE_ana(void) = default;
    virtual ~CRE_ana(void) = default;
    /*@CRE_ana(vector)
     * reassign parameters in pond
     */
    double get_emissivity_t(const vec3_t<double> &,Pond *,Grid_cre *,const double &) override;
    double get_emissivity_p(const vec3_t<double> &,Pond *,Grid_cre *,const double &) override;
private:
    /*@flux
     * return CRE flux at given CRE energy
     * input CRE energy at CGS units
     * output at [GeV m^2 s sr]^-1 units
     */
    double flux(const vec3_t<double> &,Pond *,const double &);
    /*@flux_param
     * flux index and normalization wrt gamma (Lorentz factor) at given position
     */
    void flux_param(const vec3_t<double> &,Pond *,double &,double &);
    /*@write_grid
     * write analytical CRE FLUX to discrete grid
     */
    void write_grid(Pond *,Grid_cre *) override;
};

// Numerical CRE flux stored in grid_cre
class CRE_num final: public CRE{
public:
    CRE_num(void) = default;
    virtual ~CRE_num(void) = default;
    double get_emissivity_t(const vec3_t<double> &,Pond *,Grid_cre *,const double &) override;
    double get_emissivity_p(const vec3_t<double> &,Pond *,Grid_cre *,const double &) override;
private:
    /*@read_grid
     * read CRE flux from grid at given position
     * (E_index, sylindrical_r, sylindrical_z) with {r,z} in cgs units
     * actual value of E is calculated from {E_index,Ek_min,Ek_fact} in get_emissivity
     * automatically select bi/trilinear interpolation according to
     * 2+1/3+1 spatial-spectral CRE flux grid
     */
    double read_grid(const unsigned int &,const vec3_t<double> &,Grid_cre *) override;
    
};
#endif

// END
