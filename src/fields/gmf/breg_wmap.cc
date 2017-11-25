#include <iostream>
#include <vec3.h>
#include <vector>
#include <array>
#include <cmath>

#include "param.h"
#include "grid.h"
#include "breg.h"
#include "cgs_units_file.h"
#include "namespace_toolkit.h"

using namespace std;

/* wmap-3yr */
vec3_t<double> Breg_wmap::breg(const vec3_t<double> &pos,Param *par){
    vec3_t<double> b_vec3 {0.,0.,0.};
    const double r {sqrt(pos.x*pos.x + pos.y*pos.y)};
    if (r>(20.*CGS_U_kpc) or r<(3.*CGS_U_kpc )) {
        return b_vec3;
    }
    // units
    const double b0 {par->breg_wmap.b0};
    const double psi0 {par->breg_wmap.psi0};
    const double psi1 {par->breg_wmap.psi1};
    const double chi0 {par->breg_wmap.chi0};
    const double phi {atan2(pos.y,pos.x)};
    const double psi {psi0 + psi1*log(r/(8.*CGS_U_kpc))};
    const double chi {chi0*tanh(pos.z/(1.*CGS_U_kpc))};
    const vec3_t<double> b_cyl {b0*sin(psi)*cos(chi),
        b0*cos(psi)*cos(chi),
        b0*sin(chi)};
    toolkit::Cyl2Cart(phi,b_cyl,b_vec3);
    return b_vec3;
}