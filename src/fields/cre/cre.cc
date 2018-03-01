#include <iostream>
#include <vec3.h>
#include <vector>
#include <array>
#include <cmath>
#include <gsl/gsl_sf_synchrotron.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_integration.h>
#include <cre.h>
#include <param.h>
#include <grid.h>
#include <cgs_units_file.h>
#include <namespace_toolkit.h>
#include <ap_err.h>
using namespace std;

double CRE::get_emissivity_t(const vec3_t<double> &,Param *,Grid_cre *,const double &){
    ap_err("dynamic binding fial");
    exit(1);
}

double CRE::get_emissivity_p(const vec3_t<double> &,Param *,Grid_cre *,const double &){
    ap_err("dynamic binding fial");
    exit(1);
}

double CRE::read_grid(const std::size_t &,const vec3_t<double> &,Grid_cre *){
    ap_err("dynamic binding fial");
    exit(1);
}

void CRE::write_grid(Param *,Grid_cre *){
    ap_err("dynamic binding fial");
    exit(1);
}

// END
