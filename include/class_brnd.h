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
    virtual void write_grid_iso(Pond *,Grid_brnd *);
    /*@write_grid_plus
     * for dynamic binding only
     */
    //virtual void write_grid_ani(Pond *,Breg *,Grid_breg *,Grid_brnd *);
    /*@ b_spec
     * isotropic power-spectrum
     */
    virtual double b_spec(const double &,Pond *);
    /*@rescal_fact
     * field energy density rescal factor
     */
    virtual double rescal_fact(const vec3 &,Pond *);
};


/* isotropic random field */
// not final derived class
class Brnd_iso : public Brnd{
    public:
    Brnd_iso(void) = default;
    
    vec3 get_brnd(const vec3 &,Pond *,Grid_brnd *) override;
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


/* anisotropic random field
class Brnd_ani final : public Brnd_iso{
    public:
    Brnd_ani(void) = default;
    Brnd_ani(Pond *,Grid_brnd *);
    // use Breg::breg function
    void write_grid_ani(Pond *, Breg *, Grid_breg *, Grid_brnd *) override;
};
*/
#endif
// END
