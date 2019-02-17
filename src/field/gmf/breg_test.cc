#include <cmath>

#include <hvec.h>

#include <param.h>
#include <breg.h>

hvec<3,double> Breg_test::breg (const hvec<3,double> &,
                                const Param *par) const{
    return hvec<3,double> {par->breg_test.b0*std::cos(par->breg_test.l0),
        par->breg_test.b0*std::sin(par->breg_test.l0),
        par->breg_test.b0*par->breg_test.r};
}

// END
