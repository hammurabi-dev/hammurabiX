/*
 *@file: class_cre.h
 *@author: Jiaxin Wang
 *@email: jiwang@sissa.it
 *@brief: synchrotron emissivity calculator
 */
#ifndef HAMMURABI_CRE_H
#define HAMMURABI_CRE_H

#include <gsl/gsl_sf_synchrotron.h>
#include "class_pond.h"
#include "class_grid.h"

class CRE{
    public:
    CRE(void) = default;
    virtual ~CRE(void) = default;
    /*@get_emissivity
     * get J_tot(cue=1)/J_pol(cue=0) at given position
     */
    virtual double get_emissivity(const vec3 &,Pond *,Grid_cre *,const double &,const bool &);
    virtual double read_grid(const unsigned int &, const vec3 &,Grid_cre *);
    virtual void write_grid(Pond *,Grid_cre *);
    
};

// Analytical CRE
class CRE_ana : public CRE{
    public:
    CRE_ana(void) = default;
    virtual ~CRE_ana(void) = default;
    /*@CRE_ana(vector)
     * reassign parameters in pond
     */
    //CRE_ana(const std::vector<double>&,Pond *);
    double get_emissivity(const vec3 &,Pond *,Grid_cre *,const double &,const bool &) override;
    private:
    /*@flux
     * return CRE flux at given CRE energy
     * input CRE energy at CGS units
     * output at [GeV m^2 s sr]^-1 units
     */
    double flux(const vec3 &,Pond *,const double &);
    /*@flux_param
     * flux index and normalization wrt gamma (Lorentz factor) at given position
     */
    void flux_param(const vec3 &,Pond *,double &,double &);
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
    double get_emissivity(const vec3 &,Pond *,Grid_cre *,const double &,const bool &) override;
    private:
    /*@read_grid
     * read CRE flux from grid at given position
     * (E_index, sylindrical_r, sylindrical_z) with {r,z} in cgs units
     * actual value of E is calculated from {E_index,Ek_min,Ek_fact} in get_emissivity
     * automatically select bi/trilinear interpolation according to
     * 2+1/3+1 spatial-spectral CRE flux grid
     */
    double read_grid(const unsigned int &,const vec3 &,Grid_cre *) override;
    
};
#endif

// END
