///
/// storing parameters (fixed or free) for physical models
///
#ifndef GENERIC_POND_H
#define GENERIC_POND_H

#include <string>
#include <vector>
#include <tinyxml2.h>
#include <vec3.h>

using namespace tinyxml2;

class Pond {
public:
    Pond(std::string);
    virtual ~Pond(void) = default;
    // observer
    vec3_t<double> SunPosition;
    
    // magnetic field
    // wmap3yr
    struct param_breg_wmap{
        double b0;
        double psi0;
        double psi1;
        double chi0;
    }breg_wmap;
    // verify
    struct param_breg_verify{
        double b0;
        double l0;
        double r;
    }breg_verify;
    // random seed
    std::size_t brnd_seed;
    // isotropic
    struct param_brnd_iso{
        double rms;
        double k0;
        double a0;
    }brnd_iso;
    // global aniso
    struct param_brnd_anig{
        double rho;
    }brnd_anig;
    // local aniso
    struct param_brnd_anil{
        double rms;
        double rf;
        double rs;
        
        double aa0;
        double af0;
        double as0;
        double k0;
        
        double ma;
        double beta;
    }brnd_anil;
    // rescaling parameters for random b
    struct param_brnd_scal{
        double r0;
        double z0;
    }brnd_scal;
    
    // FE
    // ymw16
    struct param_fereg_ymw16{
        double R_warp, R0;
        double t0_Gamma_w;
        double t1_Ad, t1_Bd, t1_n1, t1_H1;
        double t2_A2, t2_B2, t2_n2, t2_K2;
        double t3_B2s, t3_Ka, t3_narm[5], t3_warm[5], t3_Aa, t3_ncn, t3_wcn, t3_thetacn, t3_nsg, t3_wsg, t3_thetasg;
        double t4_ngc, t4_Agc, t4_Hgc;
        double t5_Kgn, t5_ngn, t5_Wgn, t5_Agn;
        double t6_J_LB, t6_nlb1, t6_detlb1, t6_wlb1, t6_hlb1, t6_thetalb1, t6_nlb2, t6_detlb2, t6_wlb2, t6_hlb2, t6_thetalb2;
        double t7_nLI, t7_RLI, t7_WLI, t7_detthetaLI, t7_thetaLI;
    }fereg_ymw16;
    // verify
    struct param_fereg_verify{
        double n0;
        double r0;
    }fereg_verify;
    // random seed
    std::size_t fernd_seed;
    // isotropic
    struct param_fernd_iso{
        double rms;
        double k0;
        double a0;
    }fernd_iso;
    // rescaling parameters for turbulent FE
    struct param_fernd_scal{
        double r0;
        double z0;
    }fernd_scal;
    
    // CRE
    // analytical
    double sim_freq;
    //std::vector<double> creana;
    struct param_cre_ana{
        double alpha;
        double beta;
        double theta;
        double hr;
        double hz;
        double je;
    }cre_ana;
    struct param_cre_verify{
        double alpha;
        double r0;
        double je;
    }cre_verify;
    
private:
    void b_pond(XMLDocument *);
    void fe_pond(XMLDocument *);
    void cre_pond(XMLDocument *);
    
    // auxiliary functions
    std::string FetchString(XMLElement *,std::string);
    int FetchInt(XMLElement *,std::string);
    unsigned int FetchUnsigned(XMLElement *,std::string);
    bool FetchBool(XMLElement *,std::string);
    double FetchDouble(XMLElement *,std::string);
};
#endif
