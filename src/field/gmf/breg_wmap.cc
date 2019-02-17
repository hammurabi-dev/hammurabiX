#include <cmath>

#include <hvec.h>

#include <param.h>
#include <breg.h>
#include <cgs_units_file.h>

// wmap-3yr
hvec<3,double> Breg_wmap::breg (const hvec<3,double> &pos,
                                const Param *par) const{
    const double r {std::sqrt(pos[0]*pos[0] + pos[1]*pos[1])};
    if (r>(20.*CGS_U_kpc) or r<(3.*CGS_U_kpc)){
        return hvec<3,double> {0.,0.,0.};
    }
    const double b0 {par->breg_wmap.b0};
    const double psi0 {par->breg_wmap.psi0};
    const double psi1 {par->breg_wmap.psi1};
    const double chi0 {par->breg_wmap.chi0};
    const double cos_p {pos[0]/r};
    const double sin_p {pos[1]/r};
    const double psi {psi0 + psi1*std::log(r/(8.*CGS_U_kpc))};
    const double chi {chi0*std::tanh(pos[2]/(1.*CGS_U_kpc))};
    return hvec<3,double> {
        b0*(std::sin(psi)*std::cos(chi)*cos_p - std::cos(psi)*std::cos(chi)*sin_p),
        b0*(std::sin(psi)*std::cos(chi)*sin_p + std::cos(psi)*std::cos(chi)*cos_p),
        b0*std::sin(chi)};
}

// END
