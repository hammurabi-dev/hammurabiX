/*
 *@file: namespace_dragonize.h
 *@author: Jiaxin Wang
 *@email: jiwang@sissa.it
 *@brief: numerical tools to satisfy Dragon
 */
#ifndef DRAGON_DRAGONIZE_H
#define DRAGON_DRAGONIZE_H

#include <vec3.h>

#include "class_breg.h"
#include "class_grid.h"
#include "class_pond.h"

namespace dragonize{
    // Axial symmetrize regular magnetic field for 2D grid
    // with given Breg object
    // r and z in cgs units
    vec3 get_breg_2D(const double &, const double &, Breg *, Pond *,Grid_breg *);
    // get versor with given field
    vec3 get_versor(const vec3 &);
    
    // get energy density with given field (in cgs units)
    inline double get_energy(const vec3 &b){
        return b.SquaredLength()/(8.*CGS_U_pi);
    }
    
    // get random magnetic field energy density
    // with given regular magnetic field
    // including isotropic and anisotropic random fields
    // scaling included
    double get_b_energy(const vec3 &pos, const vec3&B, Brnd *, Pond *)
};

#endif
