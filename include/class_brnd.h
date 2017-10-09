/*
 *@file: class_brnd.h
 *@author: Jiaxin Wang
 *@email: jiwang@sissa.it
 *@brief: iso/aniso-tropic turbulent magnetic field generator
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

class Brnd{
    public:
    Brnd(void) = default;
    virtual ~Brnd(void) = default;
    /*@get_brnd
     * base class version returns zero field
     * derived class version call read_grid
     */
    virtual vec3 get_brnd(const vec3 &,Grid_brnd *);
    /*@read_grid
     * read field from grid with trilinear interpolation
     * user has to write_grid ahead
     */
    virtual vec3 read_grid(const vec3 &,Grid_brnd *);
    /*@write_grid_iso
     * write field to grid (model dependent)
     * user can export_grid to binary file with Grid_xxx::export_grid
     */
    virtual void write_grid_iso(Pond *,Grid_brnd *);
    /*@write_grid_anig
     * for dynamic binding only
     */
    virtual void write_grid_ani(Pond *,Breg *,Grid_breg *,Grid_brnd *);
    /*@ b_spec
     * isotropic power-spectrum
     */
    virtual double b_spec(const double &,Pond *);
    /*@rescal_fact
     * field energy density rescal factor
     */
    virtual double rescal_fact(const vec3 &,Pond *);
};


// isotropic turbulent field
class Brnd_iso : public Brnd{
    public:
    Brnd_iso(void) = default;
    virtual ~Brnd_iso(void) = default;
    vec3 get_brnd(const vec3 &,Grid_brnd *) override;
    void write_grid_iso(Pond *,Grid_brnd *) override;
    
    protected:
    /*@complex2real
     * fetch real b(x) values from in-plane DFT
     */
    void complex2real(const fftw_complex *,double *,const unsigned long int &);
    /*@gramschmidt
     * orthogonalization process
     */
    vec3 gramschmidt(const vec3 &,const vec3 &);
};


// global anisotropic random field
class Brnd_anig final : public Brnd_iso{
    public:
    Brnd_anig(void) = default;
    virtual ~Brnd_anig(void) = default;
    void write_grid_ani(Pond *,Breg *,Grid_breg *,Grid_brnd *) override;
    
    private:
    /*@anisotropy
     * calculate anisotropy at given point
     */
    double anisotropy(const vec3 &,vec3 &,Pond *,Breg *,Grid_breg *);
};

// local anisotropic random field

class Brnd_anil final : public Brnd_iso{
    public:
    Brnd_anil(void) = default;
    virtual ~Brnd_anil(void) = default;
    // use Breg::breg function
    void write_grid_ani(Pond *,Breg *,Grid_breg *,Grid_brnd *) override;
    
    private:
    /*@anisotropy
     * offer anisotropy tensor
     */
    void anisotropy(double *,Pond *,Breg *,Grid_breg *);
};

#endif

// END
