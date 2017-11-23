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

vec3_t<double> Breg_verify::breg(const vec3_t<double> &,Pond *par){
    // units
    const double b0 {par->breg_verify.b0};
    const double l0 {par->breg_verify.l0};
    const double r {par->breg_verify.r};
    return vec3_t<double> {b0*cos(l0),b0*sin(l0),b0*r};
}
