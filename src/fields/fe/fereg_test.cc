#include <vec3.h>

#include <fereg.h>
#include <param.h>

double FEreg_test::density (const vec3_t<double> &pos,
                            const Param *par) const{
    if ((pos-par->SunPosition).Length() > par->fereg_test.r0){
        return 0.;
    }
    return par->fereg_test.n0;
}

// END
