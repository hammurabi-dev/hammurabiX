#include <string>
#include <cmath>
#include <memory>

#include <vec3.h>
#include <tinyxml2.h>

#include <param.h>
#include <cgs_units_file.h>
#include <namespace_toolkit.h>

Param::Param (const std::string file_name){
    std::unique_ptr<tinyxml2::XMLDocument> doc = toolkit::loadxml(file_name);
    // gc sun position
    tinyxml2::XMLElement *ptr {toolkit::tracexml(doc.get(),{"Grid","SunPosition"})};
    SunPosition = vec3_t<double> {CGS_U_kpc*toolkit::fetchdouble(ptr,"value","x"),
        CGS_U_kpc*toolkit::fetchdouble(ptr,"value","y"),
        CGS_U_pc*toolkit::fetchdouble(ptr,"value","z")};
    
    b_param (doc.get());
    fe_param (doc.get());
    cre_param (doc.get());
}

// magnetic field
void Param::b_param (tinyxml2::XMLDocument *doc){
    tinyxml2::XMLElement *ptr {toolkit::tracexml(doc,{"MagneticField"})};
    std::string breg_type {toolkit::fetchstring(ptr,"type","Regular")};
    // bwmap
    if(breg_type=="WMAP"){
        tinyxml2::XMLElement *subptr {toolkit::tracexml(doc,{"MagneticField","Regular","WMAP"})};
        breg_wmap.b0 = toolkit::fetchdouble(subptr,"value","b0")*CGS_U_muGauss; //microGauss
        breg_wmap.psi0 = toolkit::fetchdouble(subptr,"value","psi0")*CGS_U_rad; //rad
        breg_wmap.psi1 = toolkit::fetchdouble(subptr,"value","psi1")*CGS_U_rad; //rad
        breg_wmap.chi0 = toolkit::fetchdouble(subptr,"value","chi0")*CGS_U_rad; //rad
    }
    // bjaffe
    else if(breg_type=="Jaffe"){
        tinyxml2::XMLElement *subptr {toolkit::tracexml(doc,{"MagneticField","Regular","Jaffe"})};
        breg_jaffe.quadruple = toolkit::fetchbool(subptr,"cue","quadruple");
        breg_jaffe.bss = toolkit::fetchbool(subptr,"cue","bss");
        breg_jaffe.disk_amp = toolkit::fetchdouble(subptr,"value","disk_amp")*CGS_U_muGauss; //microG
        breg_jaffe.disk_z0 = toolkit::fetchdouble(subptr,"value","disk_z0")*CGS_U_kpc; //kpc
        breg_jaffe.halo_amp = toolkit::fetchdouble(subptr,"value","halo_amp")*CGS_U_muGauss; //microG
        breg_jaffe.halo_z0 = toolkit::fetchdouble(subptr,"value","halo_z0")*CGS_U_kpc; //kpc
        breg_jaffe.r_inner = toolkit::fetchdouble(subptr,"value","r_inner")*CGS_U_kpc; //kpc
        breg_jaffe.r_scale = toolkit::fetchdouble(subptr,"value","r_scale")*CGS_U_kpc; //kpc
        breg_jaffe.r_peak = toolkit::fetchdouble(subptr,"value","r_peak")*CGS_U_kpc; //kpc
        breg_jaffe.ring = toolkit::fetchbool(subptr,"cue","ring");
        breg_jaffe.ring_amp = toolkit::fetchdouble(subptr,"value","ring_amp")*CGS_U_muGauss; //microG
        breg_jaffe.ring_r = toolkit::fetchdouble(subptr,"value","ring_r")*CGS_U_kpc; //kpc
        breg_jaffe.bar = toolkit::fetchbool(subptr,"cue","bar");
        breg_jaffe.bar_amp = toolkit::fetchdouble(subptr,"value","bar_amp")*CGS_U_muGauss; //microG
        breg_jaffe.bar_a = toolkit::fetchdouble(subptr,"value","bar_a")*CGS_U_kpc; //kpc
        breg_jaffe.bar_b = toolkit::fetchdouble(subptr,"value","bar_b")*CGS_U_kpc; //kpc
        breg_jaffe.bar_phi0 = toolkit::fetchdouble(subptr,"value","bar_phi0")*CGS_U_rad; //rad
        breg_jaffe.arm_num = toolkit::fetchunsigned(subptr,"value","arm_num");
        breg_jaffe.arm_r0 = toolkit::fetchdouble(subptr,"value","arm_r0")*CGS_U_kpc; //kpc
        breg_jaffe.arm_z0 = toolkit::fetchdouble(subptr,"value","arm_z0")*CGS_U_kpc; //kpc
        breg_jaffe.arm_phi0.push_back(toolkit::fetchdouble(subptr,"value","arm_phi1")*CGS_U_rad); //rad
        breg_jaffe.arm_phi0.push_back(toolkit::fetchdouble(subptr,"value","arm_phi2")*CGS_U_rad); //rad
        breg_jaffe.arm_phi0.push_back(toolkit::fetchdouble(subptr,"value","arm_phi3")*CGS_U_rad); //rad
        breg_jaffe.arm_phi0.push_back(toolkit::fetchdouble(subptr,"value","arm_phi4")*CGS_U_rad); //rad
        breg_jaffe.arm_amp.push_back(toolkit::fetchdouble(subptr,"value","arm_amp1")*CGS_U_muGauss); //microG
        breg_jaffe.arm_amp.push_back(toolkit::fetchdouble(subptr,"value","arm_amp2")*CGS_U_muGauss); //microG
        breg_jaffe.arm_amp.push_back(toolkit::fetchdouble(subptr,"value","arm_amp3")*CGS_U_muGauss); //microG
        breg_jaffe.arm_amp.push_back(toolkit::fetchdouble(subptr,"value","arm_amp4")*CGS_U_muGauss); //microG
        breg_jaffe.arm_pitch = toolkit::fetchdouble(subptr,"value","arm_pitch")*CGS_U_rad; //rad
        breg_jaffe.comp_c = toolkit::fetchdouble(subptr,"value","comp_c");
        breg_jaffe.comp_d = toolkit::fetchdouble(subptr,"value","comp_d")*CGS_U_kpc; //kpc
        breg_jaffe.comp_r = toolkit::fetchdouble(subptr,"value","comp_r")*CGS_U_kpc; //kpc
        breg_jaffe.comp_p = toolkit::fetchdouble(subptr,"value","comp_p");
    }
#ifndef NDEBUG
    // testing
    else if(breg_type=="Test"){
        tinyxml2::XMLElement *subptr {toolkit::tracexml(doc,{"MagneticField","Regular","Test"})};
        breg_test.b0 = toolkit::fetchdouble(subptr,"value","b0")*CGS_U_muGauss; //microGauss
        breg_test.l0 = toolkit::fetchdouble(subptr,"value","l0")*CGS_U_rad; //rad
        breg_test.r = toolkit::fetchdouble(subptr,"value","r");
    }
#endif
    
    if(toolkit::fetchbool(ptr,"cue","Random")){
        // random seed
        brnd_seed = toolkit::fetchunsigned(ptr,"seed","Random");
        std::string brnd_type {toolkit::fetchstring(ptr,"type","Random")};
        // brnd_global
        if(brnd_type=="Global"){
            tinyxml2::XMLElement *subptr {toolkit::tracexml(doc,{"MagneticField","Random","Global"})};
            std::string brnd_method {toolkit::fetchstring(subptr,"type")};
            if(brnd_method=="ES"){
                subptr = toolkit::tracexml(doc,{"MagneticField","Random","Global","ES"});
                brnd_es.rms = toolkit::fetchdouble(subptr,"value","rms")*CGS_U_muGauss;
                brnd_es.k0 = toolkit::fetchdouble(subptr,"value","k0");
                brnd_es.a0 = toolkit::fetchdouble(subptr,"value","a0");
                brnd_es.rho = toolkit::fetchdouble(subptr,"value","rho");
                brnd_es.r0 = toolkit::fetchdouble(subptr,"value","r0")*CGS_U_kpc;
                brnd_es.z0 = toolkit::fetchdouble(subptr,"value","z0")*CGS_U_kpc;
            }
            else if(brnd_method=="Jaffe"){
                subptr = toolkit::tracexml(doc,{"MagneticField","Random","Global","Jaffe"});
                // to be implemented
            }
        }
        // brnd_local
        else if(brnd_type=="Local"){
            tinyxml2::XMLElement *subptr {toolkit::tracexml(doc,{"MagneticField","Random","Local"})};
            std::string brnd_method {toolkit::fetchstring(subptr,"type")};
            if(brnd_method=="MHD"){
                subptr = toolkit::tracexml(doc,{"MagneticField","Random","Local","MHD"});
                brnd_mhd.pa0 = toolkit::fetchdouble(subptr,"value","pa0")*CGS_U_muGauss*CGS_U_muGauss;
                brnd_mhd.pf0 = toolkit::fetchdouble(subptr,"value","pf0")*CGS_U_muGauss*CGS_U_muGauss;
                brnd_mhd.ps0 = toolkit::fetchdouble(subptr,"value","ps0")*CGS_U_muGauss*CGS_U_muGauss;
                brnd_mhd.aa0 = toolkit::fetchdouble(subptr,"value","aa0");
                brnd_mhd.af0 = toolkit::fetchdouble(subptr,"value","af0");
                brnd_mhd.as0 = toolkit::fetchdouble(subptr,"value","as0");
                brnd_mhd.k0 = toolkit::fetchdouble(subptr,"value","k0");
                brnd_mhd.beta = toolkit::fetchdouble(subptr,"value","beta");
                brnd_mhd.ma = toolkit::fetchdouble(subptr,"value","ma");
            }
        }
    }
}

void Param::fe_param (tinyxml2::XMLDocument *doc){
    tinyxml2::XMLElement *ptr {toolkit::tracexml(doc,{"FreeElectron"})};
    std::string fereg_type {toolkit::fetchstring(ptr,"type","Regular")};
    // YMW16
    if(fereg_type=="YMW16"){
        // Warp_Sun
        tinyxml2::XMLElement *subptr {toolkit::tracexml(doc,{"FreeElectron","Regular","YMW16","Warp"})};
        fereg_ymw16.R_warp = toolkit::fetchdouble(subptr,"value","R_warp")*CGS_U_kpc; //kpc
        fereg_ymw16.R0 = toolkit::fetchdouble(subptr,"value","R0")*CGS_U_kpc; //kpc
        fereg_ymw16.t0_Gamma_w = toolkit::fetchdouble(subptr,"value","Gamma_w");
        // thick disk
        subptr = toolkit::tracexml(doc,{"FreeElectron","Regular","YMW16","ThickDisk"});
        fereg_ymw16.t1_Ad = toolkit::fetchdouble(subptr,"value","Ad")*CGS_U_pc;//pc
        fereg_ymw16.t1_Bd = toolkit::fetchdouble(subptr,"value","Bd")*CGS_U_pc;//pc
        fereg_ymw16.t1_n1 = toolkit::fetchdouble(subptr,"value","n1");//pccm
        fereg_ymw16.t1_H1 = toolkit::fetchdouble(subptr,"value","H1")*CGS_U_pc;//pc
        // thin disk
        subptr = toolkit::tracexml(doc,{"FreeElectron","Regular","YMW16","ThinDisk"});
        fereg_ymw16.t2_A2 = toolkit::fetchdouble(subptr,"value","A2")*CGS_U_pc;//pc
        fereg_ymw16.t2_B2 = toolkit::fetchdouble(subptr,"value","B2")*CGS_U_pc;//pc
        fereg_ymw16.t2_n2 = toolkit::fetchdouble(subptr,"value","n2");//pccm
        fereg_ymw16.t2_K2 = toolkit::fetchdouble(subptr,"value","K2");
        // spiral arm
        subptr = toolkit::tracexml(doc,{"FreeElectron","Regular","YMW16","SpiralArm"});
        fereg_ymw16.t3_B2s = toolkit::fetchdouble(subptr,"value","B2s")*CGS_U_pc; //pc
        fereg_ymw16.t3_narm[0] = toolkit::fetchdouble(subptr,"value","Ele_arm_0");//pccm
        fereg_ymw16.t3_narm[1] = toolkit::fetchdouble(subptr,"value","Ele_arm_1");
        fereg_ymw16.t3_narm[2] = toolkit::fetchdouble(subptr,"value","Ele_arm_2");
        fereg_ymw16.t3_narm[3] = toolkit::fetchdouble(subptr,"value","Ele_arm_3");
        fereg_ymw16.t3_narm[4] = toolkit::fetchdouble(subptr,"value","Ele_arm_4");
        fereg_ymw16.t3_warm[0] = toolkit::fetchdouble(subptr,"value","Wid_arm_0")*CGS_U_pc;//pc
        fereg_ymw16.t3_warm[1] = toolkit::fetchdouble(subptr,"value","Wid_arm_1")*CGS_U_pc;
        fereg_ymw16.t3_warm[2] = toolkit::fetchdouble(subptr,"value","Wid_arm_2")*CGS_U_pc;
        fereg_ymw16.t3_warm[3] = toolkit::fetchdouble(subptr,"value","Wid_arm_3")*CGS_U_pc;
        fereg_ymw16.t3_warm[4] = toolkit::fetchdouble(subptr,"value","Wid_arm_4")*CGS_U_pc;
        fereg_ymw16.t3_rmin[0] = toolkit::fetchdouble(subptr,"value","Rref_arm_0")*CGS_U_kpc;//kpc
        fereg_ymw16.t3_rmin[1] = toolkit::fetchdouble(subptr,"value","Rref_arm_1")*CGS_U_kpc;
        fereg_ymw16.t3_rmin[2] = toolkit::fetchdouble(subptr,"value","Rref_arm_2")*CGS_U_kpc;
        fereg_ymw16.t3_rmin[3] = toolkit::fetchdouble(subptr,"value","Rref_arm_3")*CGS_U_kpc;
        fereg_ymw16.t3_rmin[4] = toolkit::fetchdouble(subptr,"value","Rref_arm_4")*CGS_U_kpc;
        fereg_ymw16.t3_phimin[0] = toolkit::fetchdouble(subptr,"value","Phiref_arm_0")*CGS_U_rad;//rad
        fereg_ymw16.t3_phimin[1] = toolkit::fetchdouble(subptr,"value","Phiref_arm_1")*CGS_U_rad;//rad
        fereg_ymw16.t3_phimin[2] = toolkit::fetchdouble(subptr,"value","Phiref_arm_2")*CGS_U_rad;//rad
        fereg_ymw16.t3_phimin[3] = toolkit::fetchdouble(subptr,"value","Phiref_arm_3")*CGS_U_rad;//rad
        fereg_ymw16.t3_phimin[4] = toolkit::fetchdouble(subptr,"value","Phiref_arm_4")*CGS_U_rad;//rad
        fereg_ymw16.t3_tpitch[0] = tan(toolkit::fetchdouble(subptr,"value","pitch_arm_0")*CGS_U_rad);
        fereg_ymw16.t3_tpitch[1] = tan(toolkit::fetchdouble(subptr,"value","pitch_arm_1")*CGS_U_rad);
        fereg_ymw16.t3_tpitch[2] = tan(toolkit::fetchdouble(subptr,"value","pitch_arm_2")*CGS_U_rad);
        fereg_ymw16.t3_tpitch[3] = tan(toolkit::fetchdouble(subptr,"value","pitch_arm_3")*CGS_U_rad);
        fereg_ymw16.t3_tpitch[4] = tan(toolkit::fetchdouble(subptr,"value","pitch_arm_4")*CGS_U_rad);
        fereg_ymw16.t3_cpitch[0] = cos(toolkit::fetchdouble(subptr,"value","pitch_arm_0")*CGS_U_rad);
        fereg_ymw16.t3_cpitch[1] = cos(toolkit::fetchdouble(subptr,"value","pitch_arm_1")*CGS_U_rad);
        fereg_ymw16.t3_cpitch[2] = cos(toolkit::fetchdouble(subptr,"value","pitch_arm_2")*CGS_U_rad);
        fereg_ymw16.t3_cpitch[3] = cos(toolkit::fetchdouble(subptr,"value","pitch_arm_3")*CGS_U_rad);
        fereg_ymw16.t3_cpitch[4] = cos(toolkit::fetchdouble(subptr,"value","pitch_arm_4")*CGS_U_rad);
        fereg_ymw16.t3_Aa = toolkit::fetchdouble(subptr,"value","Aa")*CGS_U_pc;//pc
        fereg_ymw16.t3_Ka = toolkit::fetchdouble(subptr,"value","Ka");
        fereg_ymw16.t3_ncn = toolkit::fetchdouble(subptr,"value","ncn");
        fereg_ymw16.t3_thetacn = toolkit::fetchdouble(subptr,"value","thetacn");//deg
        fereg_ymw16.t3_wcn = toolkit::fetchdouble(subptr,"value","wcn");//deg
        fereg_ymw16.t3_nsg = toolkit::fetchdouble(subptr,"value","nsg");
        fereg_ymw16.t3_thetasg = toolkit::fetchdouble(subptr,"value","thetasg");//deg
        fereg_ymw16.t3_wsg = toolkit::fetchdouble(subptr,"value","wsg");//deg
        // gc
        subptr = toolkit::tracexml(doc,{"FreeElectron","Regular","YMW16","GalCenter"});
        fereg_ymw16.t4_ngc = toolkit::fetchdouble(subptr,"value","ngc");//pccm
        fereg_ymw16.t4_Agc = toolkit::fetchdouble(subptr,"value","Agc")*CGS_U_pc;//pc
        fereg_ymw16.t4_Hgc = toolkit::fetchdouble(subptr,"value","Hgc")*CGS_U_pc;//pc
        // Gum
        subptr = toolkit::tracexml(doc,{"FreeElectron","Regular","YMW16","GumNebula"});
        fereg_ymw16.t5_ngn = toolkit::fetchdouble(subptr,"value","ngn");//pccm
        fereg_ymw16.t5_Wgn = toolkit::fetchdouble(subptr,"value","Wgn")*CGS_U_pc;//pc
        fereg_ymw16.t5_Agn = toolkit::fetchdouble(subptr,"value","Agn")*CGS_U_pc;//pc
        fereg_ymw16.t5_Kgn = toolkit::fetchdouble(subptr,"value","Kgn");
        // Local Bubble
        subptr = toolkit::tracexml(doc,{"FreeElectron","Regular","YMW16","LocalBubble"});
        fereg_ymw16.t6_J_LB = toolkit::fetchdouble(subptr,"value","J_LB");
        fereg_ymw16.t6_nlb1 = toolkit::fetchdouble(subptr,"value","nlb1");//pccm
        fereg_ymw16.t6_thetalb1 = toolkit::fetchdouble(subptr,"value","thetalb1");//deg
        fereg_ymw16.t6_detlb1 = toolkit::fetchdouble(subptr,"value","detlb1");//deg
        fereg_ymw16.t6_wlb1 = toolkit::fetchdouble(subptr,"value","wlb1")*CGS_U_pc;//pc
        fereg_ymw16.t6_hlb1 = toolkit::fetchdouble(subptr,"value","hlb1")*CGS_U_pc;//pc
        fereg_ymw16.t6_nlb2 = toolkit::fetchdouble(subptr,"value","nlb2");//pccm
        fereg_ymw16.t6_thetalb2 = toolkit::fetchdouble(subptr,"value","thetalb2");//deg
        fereg_ymw16.t6_detlb2 = toolkit::fetchdouble(subptr,"value","detlb2");//deg
        fereg_ymw16.t6_wlb2 = toolkit::fetchdouble(subptr,"value","wlb2")*CGS_U_pc;//pc
        fereg_ymw16.t6_hlb2 = toolkit::fetchdouble(subptr,"value","hlb2")*CGS_U_pc;//pc
        // Loop I
        subptr = toolkit::tracexml(doc,{"FreeElectron","Regular","YMW16","LoopI"});
        fereg_ymw16.t7_nLI = toolkit::fetchdouble(subptr,"value","nLI");//pccm
        fereg_ymw16.t7_RLI = toolkit::fetchdouble(subptr,"value","RLI")*CGS_U_pc;//pc
        fereg_ymw16.t7_WLI = toolkit::fetchdouble(subptr,"value","WLI")*CGS_U_pc;//pc
        fereg_ymw16.t7_detthetaLI = toolkit::fetchdouble(subptr,"value","detthetaLI");//deg
        fereg_ymw16.t7_thetaLI = toolkit::fetchdouble(subptr,"value","thetaLI");//deg
    }
#ifndef NDEBUG
    // testing
    else if(fereg_type=="Test"){
        tinyxml2::XMLElement *subptr {toolkit::tracexml(doc,{"FreeElectron","Regular","Test"})};
        fereg_test.n0 = toolkit::fetchdouble(subptr,"value","n0");
        fereg_test.r0 = toolkit::fetchdouble(subptr,"value","r0")*CGS_U_kpc; //kpc
    }
#endif
    
    if(toolkit::fetchbool(ptr,"cue","Random")){
        // random seed
        fernd_seed = toolkit::fetchunsigned(ptr,"seed","Random");
        std::string fernd_type {toolkit::fetchstring(ptr,"type","Random")};
        // global turbulent
        if(fernd_type=="Global"){
            tinyxml2::XMLElement *subptr {toolkit::tracexml(doc,{"FreeElectron","Random","Global"})};
            std::string fernd_method {toolkit::fetchstring(subptr,"type")};
            if(fernd_method=="DFT"){
                subptr = toolkit::tracexml(doc,{"FreeElectron","Random","Global","DFT"});
                fernd_dft.rms = toolkit::fetchdouble(subptr,"value","rms");
                fernd_dft.k0 = toolkit::fetchdouble(subptr,"value","k0");
                fernd_dft.a0 = toolkit::fetchdouble(subptr,"value","a0");
                fernd_dft.r0 = toolkit::fetchdouble(subptr,"value","r0")*CGS_U_kpc;
                fernd_dft.z0 = toolkit::fetchdouble(subptr,"value","z0")*CGS_U_kpc;
            }
        }
    }
}

void Param::cre_param (tinyxml2::XMLDocument *doc){
    tinyxml2::XMLElement *ptr {toolkit::tracexml(doc,{"CRE"})};
    if (toolkit::tracexml(doc,{"Obsout","Sync"})!=nullptr){
        tinyxml2::XMLElement *subptr {toolkit::tracexml(doc,{"Obsout"})};
        sim_freq = toolkit::fetchdouble(subptr,"freq","Sync")*CGS_U_GHz;
    }
    else{
        sim_freq = 0.;
    }
    std::string cre_type {ptr->Attribute("type")};
    // analytical
    if(cre_type=="Analytic"){
        tinyxml2::XMLElement *subptr {toolkit::tracexml(doc,{"CRE","Analytic"})};
        cre_ana.alpha = toolkit::fetchdouble(subptr,"value","alpha");
        cre_ana.beta = toolkit::fetchdouble(subptr,"value","beta");
        cre_ana.theta = toolkit::fetchdouble(subptr,"value","theta");
        cre_ana.r0 = toolkit::fetchdouble(subptr,"value","r0")*CGS_U_kpc; //kpc
        cre_ana.z0 = toolkit::fetchdouble(subptr,"value","z0")*CGS_U_kpc; //kpc
        cre_ana.E0 = toolkit::fetchdouble(subptr,"value","E0")*CGS_U_GeV; //GeV
        cre_ana.j0 = toolkit::fetchdouble(subptr,"value","j0");
    }
#ifndef NDEBUG
    // testing
    else if(cre_type=="Test"){
        tinyxml2::XMLElement *subptr {toolkit::tracexml(doc,{"CRE","Test"})};
        cre_test.alpha = toolkit::fetchdouble(subptr,"value","alpha");
        cre_test.r0 = toolkit::fetchdouble(subptr,"value","r0")*CGS_U_kpc; //kpc
        cre_test.E0 = toolkit::fetchdouble(subptr,"value","E0")*CGS_U_GeV; //GeV
        cre_test.j0 = toolkit::fetchdouble(subptr,"value","j0");
    }
#endif
}

// END
