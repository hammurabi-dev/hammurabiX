/*
 *@file: class_brnd.h
 *@author: Jiaxin Wang
 *@email: jiwang@sissa.it
 *@brief: iso/aniso-tropic random magnetic field generator
 */
#ifndef HAMMURABI_BRND_H
#define HAMMURABI_BRND_H

#include <fftw3.h>
#include <vec3.h>
#include <vector>
#include <gsl/gsl_integration.h>

#include "class_pond.h"
#include "class_grid.h"
#include "cgs_units_file.h"
#include "class_breg.h"

/* base class */
class Brnd{
    public:
    Brnd(void) = default;
    virtual ~Brnd(void) = default;
    /*@get_brnd
     * base class version returns zero field
     * derived class version call read_grid
     */
    virtual vec3 get_brnd(const vec3 &,Pond *,Grid_brnd *);
    /*@read_grid
     * read field from grid with trilinear interpolation
     * user has to write_grid ahead
     */
    virtual vec3 read_grid(const vec3 &,Grid_brnd *,Pond *);
    /*@write_grid
     * write field to grid (model dependent)
     * user can export_grid to binary file with Grid_xxx::export_grid
     */
    virtual void write_grid(Pond *,Grid_brnd *);
    /*@write_grid_plus
     * for dynamic binding only
     */
    virtual void write_grid_plus(Pond *,Breg *,Grid_breg *,Grid_brnd *);
    /*@ b_spec
     * isotropic power-spectrum
     */
    virtual double b_spec(const double &,Pond *);
    /*@get_energy
     * get theoretical energy density in cgs units with given b_spec
     */
    virtual double get_energy(Pond *,Grid_brnd *);
    /*@get_missing
     * get missing energy density (in cgs units given b_spec)
     * due to k_max in discrete simulation
     */
    virtual double get_missing(Pond *,Grid_brnd *);
    /*@xxx_rslt
     * get_energy & get_missing only calculated once
     * and stored in xxx_rslt in constructor
     */
    double get_energy_rslt, get_missing_rslt;
    /*@rescal_fact
     * field energy density rescal factor
     */
    virtual double rescal_fact(const vec3 &,Pond *);
};


/* gaussian random field */
// not final derived class
class Bgrnd : public Brnd{
    public:
    Bgrnd(void) = default;
    Bgrnd(Pond *,Grid_brnd *);
    /*@Bgrnd(vector)
     * reassign parameters in pond
     * designed for passing parameters during MCMC
     */
    //Bgrnd(const std::vector<double>&,Pond *,Grid_brnd *);
    vec3 get_brnd(const vec3 &,Pond *,Grid_brnd *) override;
    void write_grid(Pond *,Grid_brnd *) override;
    
    protected:
    /*@hermiticity
     * fix hermiticity of b(k) grid
     */
    void hermiticity(fftw_complex *,const unsigned int &,const unsigned int &,const unsigned int &);
    /*@complex2real
     * fetch real b(x) values from in-plane DFT
     */
    void complex2real(const fftw_complex *,double *,const unsigned long int &);
};

/* gaussian random field + anisotropic random field */
class Bfrnd final : public Bgrnd{
    public:
    Bfrnd(void) = default;
    Bfrnd(Pond *,Grid_brnd *);
    //Bfrnd(const std::vector<double>&, Pond *,Grid_brnd *);
    // use Breg::breg function
    void write_grid_plus(Pond *, Breg *, Grid_breg *, Grid_brnd *) override;
};

#endif
// END
