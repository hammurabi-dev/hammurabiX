/*
 *@file: class_pond.h
 *@author: Jiaxin Wang
 *@email: jiwang@sissa.it
 *@brief: storing parameters (fixed or free) for various physical models
 */
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
    vec3 SunPosition;
    
    // magnetic field
    // wmap3yr
    std::vector<double> bwmap;
    // local
    std::vector<double> blocal;
    // isotropic
    std::vector<double> brnd_iso;
    // global aniso
    std::vector<double> brnd_anig;
    // local aniso
    std::vector<double> brnd_anil;
    // rescaling parameters for random b
    std::vector<double> brnd_scal;
    
    // free electron
    // ymw16
    double R_warp, R0;
    double t0_Gamma_w, t0_z_Sun;
    double t1_Ad, t1_Bd, t1_n1, t1_H1;
    double t2_A2, t2_B2, t2_n2, t2_K2;
    double t3_B2s, t3_Ka, t3_narm[5], t3_warm[5], t3_Aa, t3_ncn, t3_wcn, t3_thetacn, t3_nsg, t3_wsg, t3_thetasg;
    double t4_ngc, t4_Agc, t4_Hgc;
    double t5_Kgn, t5_ngn, t5_Wgn, t5_Agn;
    double t6_J_LB, t6_nlb1, t6_detlb1, t6_wlb1, t6_hlb1, t6_thetalb1, t6_nlb2, t6_detlb2, t6_wlb2, t6_hlb2, t6_thetalb2;
    double t7_nLI, t7_RLI, t7_WLI, t7_detthetaLI, t7_thetaLI;
    // gaussian random
    std::vector<double> fernd_iso;
    // rescaling parameters for random fe
    std::vector<double> fernd_scal;
    
    // cre
    // analytical
    double sim_freq;
    std::vector<double> creana;
    
    private:
    std::string FetchString(XMLElement *,std::string);
    int FetchInt(XMLElement *,std::string);
    unsigned int FetchUnsigned(XMLElement *,std::string);
    bool FetchBool(XMLElement *,std::string);
    double FetchDouble(XMLElement *,std::string);
    
    void b_pond(XMLDocument *);
    void fe_pond(XMLDocument *);
    void cre_pond(XMLDocument *);
};
#endif
