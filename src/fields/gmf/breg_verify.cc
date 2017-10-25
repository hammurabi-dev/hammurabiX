#include <iostream>
#include <vec3.h>
#include <vector>
#include <array>
#include <cmath>

#include "pond.h"
#include "grid.h"
#include "breg.h"
#include "cgs_units_file.h"
#include "namespace_toolkit.h"

using namespace std;

vec3 Breg_verify::breg(const vec3 &,Pond *par){
    // units
    const double b0 {par->breg_verify.b0*CGS_U_muGauss};
    const double l0 {par->breg_verify.l0*CGS_U_pi/180.};
    const double r {par->breg_verify.r};
    
    return vec3 {b0*cos(l0),b0*sin(l0),b0*r};
}
