#include <iostream>
#include <vec3.h>

#include "class_breg.h"
#include "class_grid.h"
#include "class_pond.h"

using namespace std;

namespace dragonize{
    
    vec3 get_breg_2D(const double &r, const double &z, Breg *breg, Pond *par, Grid_breg *gbreg){
        const vec3 pos = vec3(-r,0,z);
        return breg->get_breg(pos,par,gbreg);
    }
    
    vec3 get_versor(const vec3 &b){
        if(b.Length()==0.) {return vec3(0.,0.,0.);}
        b.Normalize();
        return b;
    }
    
    double get_b_energy(const vec3 &pos, const vec3 &B, Brnd *brnd, Pond *par){
        return (brnd->get_energy_rslt + par->bfrnd[par->bfrnd.size()-1]*pow(B.Length(),2)/(8*CGS_U_pi))*brnd->rescal_fact(pos,par);
    }
}
// END
