/*
 *@file: class_fernd.h
 *@author: Jiaxin Wang
 *@email: jiwang@sissa.it
 *@brief: iso/aniso-tropic random free electron field generator
 *@note: this is just a scalar version of Brnd class
 */
#ifndef HAMMURABI_FERND_H
#define HAMMURABI_FERND_H

#include <fftw3.h>
#include <vec3.h>
#include <vector>
#include <gsl/gsl_integration.h>

#include "class_pond.h"
#include "class_grid.h"
#include "cgs_units_file.h"
#include "class_fe.h"

/* base class */
class FErnd{
    public:
    FErnd(void) = default;
    virtual ~FErnd(void) = default;
    virtual double get_fernd(const vec3 &,Pond *,Grid_fernd *);
    virtual double read_grid(const vec3 &,Grid_fernd *,Pond *);
    virtual void write_grid(Pond *,Grid_fernd *);
    
    virtual double fe_spec(const double &k,Pond *);
    /*@get_variance
     * variance of density field <n^2>
     */
    virtual double get_variance(Pond *,Grid_fernd *);
    /*@get_missing
     * missing part of density variance, necessary for Emission Measure
     */
    virtual double get_missing(Pond *,Grid_fernd *);
    /*@rescal_fact
     * spatial rescaling factor for variance <n^2>
     */
    virtual double rescal_fact(const vec3 &,Pond *);
    
    double get_variance_rslt, get_missing_rslt;
};

/* gaussian random field */
class FEgrnd : public FErnd{
    public:
    FEgrnd(Pond *,Grid_fernd *);
    //FEgrnd(const std::vector<double>&,Pond *,Grid_fernd *);
    double get_fernd(const vec3 &,Pond *,Grid_fernd *) override;
    void write_grid(Pond *,Grid_fernd *) override;
    protected:
    void hermiticity(fftw_complex *,const unsigned int &,const unsigned int &,const unsigned int &);
    void complex2real(const fftw_complex *,double *,const unsigned long int &);
};


#endif
// END
