#include <iostream>
#include <cmath>
#include "fereg.h"
#include "pond.h"
#include "cgs_units_file.h"

using namespace std;

// FEverify
double FEreg_verify::density(const vec3 &pos, Pond *par){
    if((pos-par->SunPosition).Length() > par->fereg_verify.r0*CGS_U_kpc){
        return 0.;
    }
    return par->fereg_verify.n0;
}
