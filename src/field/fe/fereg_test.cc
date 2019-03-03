#include <hvec.h>

#include <fereg.h>
#include <param.h>

double FEreg_test::density (const hvec<3,double> &pos,
                            const Param *par) const{
    if ((pos-par->observer).length() > par->fereg_test.r0){
        return 0.;
    }
    return par->fereg_test.n0;
}

// END
