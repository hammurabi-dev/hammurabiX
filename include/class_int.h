/*
 *@file: class_int.h
 *@author: Jiaxin Wang
 *@email: jiwang@sissa.it
 *@brief: line-of-sight integration calculator
 *@origin: https://sourceforge.net/projects/hammurabicode/
 * we add some new features inside
 * and delete which we don't need
 */
#ifndef HAMMURABI_INT_H
#define HAMMURABI_INT_H

#include <omp.h>
#include <string>
#include <healpix_map.h>

#include "class_breg.h"
#include "class_brnd.h"
#include "class_cre.h"
#include "class_fe.h"
#include "class_fernd.h"
#include "class_pond.h"


class Integrator {
    public:
    
    Integrator(void) = default;
    virtual ~Integrator(void) = default;
    /*@write_grid
     * assmebling pixels/shells into Healpix map
     */
    void write_grid(Breg *,Brnd *,FE *,FErnd *,CRE *,Grid_breg *,Grid_brnd *,Grid_fe *,Grid_fernd *,Grid_cre *,Grid_int *,Pond *);
    
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
    void radial_integration(struct_shell &,pointing &,struct_observables &,Breg *,Brnd *,FE *,FErnd *,CRE *,Grid_breg *,Grid_brnd *,Grid_fe *,Grid_fernd *,Grid_cre *,Grid_int *,Pond *);
    
    // auxiliary functions
    unsigned int get_shell_nside(const unsigned int &,const unsigned int &,const unsigned int &) const;
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
