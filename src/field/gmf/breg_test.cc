#include <cmath>

#include <vec3.h>

#include <param.h>
#include <breg.h>

vec3_t<double> Breg_test::breg (const vec3_t<double> &,
                                const Param *par) const{
    return vec3_t<double> {par->breg_test.b0*std::cos(par->breg_test.l0),
        par->breg_test.b0*std::sin(par->breg_test.l0),
        par->breg_test.b0*par->breg_test.r};
}

// END
