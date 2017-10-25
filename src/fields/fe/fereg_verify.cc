#include <iostream>
#include <cmath>
#include "fereg.h"
#include "pond.h"
#include "cgs_units_file.h"

using namespace std;

// FEverify
double FEreg_verify::density(const vec3 &pos, Pond *par){
    if(fabs(pos.z) > par->fereg_verify.z0*CGS_U_kpc){
        return 0.;
    }
    else if(pos.x*pos.x+pos.y*pos.y > par->fereg_verify.r0*par->fereg_verify.r0*CGS_U_kpc*CGS_U_kpc){
        return 0.;
    }
    return par->fereg_verify.n0;
}
