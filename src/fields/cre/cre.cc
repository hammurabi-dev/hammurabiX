#include <iostream>
#include <vec3.h>
#include <vector>
#include <array>
#include <cmath>
#include <gsl/gsl_sf_synchrotron.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_integration.h>
#include "cre.h"
#include "pond.h"
#include "grid.h"
#include "cgs_units_file.h"
#include "namespace_toolkit.h"

using namespace std;

double CRE::get_emissivity(const vec3 &,Pond *,Grid_cre *,const double &,const bool &){
    cerr<<"ERR:"<<__FILE__
    <<" : in function "<<__func__<<endl
    <<" at line "<<__LINE__<<endl
    <<"DYNAMIC BINDING FAILURE"<<endl;
    exit(1);
}

double CRE::read_grid(const unsigned int &,const vec3 &,Grid_cre *){
    cerr<<"ERR:"<<__FILE__
    <<" : in function "<<__func__<<endl
    <<" at line "<<__LINE__<<endl
    <<"DYNAMIC BINDING FAILURE"<<endl;
    exit(1);
}

void CRE::write_grid(Pond *,Grid_cre *){
    cerr<<"ERR:"<<__FILE__
    <<" : in function "<<__func__<<endl
    <<" at line "<<__LINE__<<endl
    <<"DYNAMIC BINDING FAILURE"<<endl;
    exit(1);
}

// END
