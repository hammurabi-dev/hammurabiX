/*
 *@file: integrator.h
 *@brief: line-of-sight integration calculator
 */
#ifndef HAMMURABI_INT_H
#define HAMMURABI_INT_H

#include <omp.h>
#include <string>
#include <healpix_map.h>
#include "breg.h"
#include "brnd.h"
#include "cre.h"
#include "fereg.h"
#include "fernd.h"
#include "pond.h"


class Integrator {
public:
    
    Integrator(void) = default;
    virtual ~Integrator(void) = default;
    /*@write_grid
     * assmebling pixels/shells into Healpix map
     */
    void write_grid(Breg *,Brnd *,FEreg *,FErnd *,CRE *,Grid_breg *,Grid_brnd *,Grid_fereg *,Grid_fernd *,Grid_cre *,Grid_int *,Pond *);
    
private:
    
    struct struct_observables{
        double Is, Qs, Us;
        double fd;
        double dm;
        double ff;
    };
    struct struct_shell{
        unsigned int shell_num;
        double d_start;
        double d_stop;
        double delta_d;
    };
    /*@radial_integration
     * commit LOS integration in one pixel one shell
     */
    void radial_integration(struct_shell &,pointing &,struct_observables &,Breg *,Brnd *,FEreg *,FErnd *,CRE *,Grid_breg *,Grid_brnd *,Grid_fereg *,Grid_fernd *,Grid_cre *,Grid_int *,Pond *);
    
    // auxiliary functions
    double get_max_shell_radius(const unsigned int &,const unsigned int &,const double &) const;
    double get_min_shell_radius(const unsigned int &,const unsigned int &,const double &) const;
    
    inline bool check_simulation_upper_limit(const double &value, const double &limit) const {
        if(value>limit) return true;
        else return false;
    }
    
    inline bool check_simulation_lower_limit(const double &value, const double &limit) const {
        if(value<limit) return true;
        else return false;
    }
    
};
#endif
