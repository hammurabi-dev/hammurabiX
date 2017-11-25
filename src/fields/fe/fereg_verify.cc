#include <iostream>
#include <cmath>
#include <vec3.h>
#include "fereg.h"
#include "param.h"
#include "cgs_units_file.h"

using namespace std;

// FEverify
double FEreg_verify::density(const vec3_t<double> &pos, Param *par){
    if((pos-par->SunPosition).Length() > par->fereg_verify.r0){
        return 0.;
    }
    return par->fereg_verify.n0;
}
