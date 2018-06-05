#include <string>
#include <cmath>
#include <memory>
#include <vec3.h>
#include <tinyxml2.h>
#include <param.h>
#include <cgs_units_file.h>
#include <namespace_toolkit.h>

using namespace tinyxml2;

Param::Param(std::string file_name){
    std::unique_ptr<XMLDocument> doc = toolkit::loadxml(file_name);
    // gc sun position
    XMLElement *ptr {toolkit::tracexml(doc.get(),{"Grid","SunPosition"})};
    SunPosition = vec3_t<double> {CGS_U_kpc*toolkit::FetchDouble(ptr,"value","x"),
        CGS_U_kpc*toolkit::FetchDouble(ptr,"value","y"),
        CGS_U_pc*toolkit::FetchDouble(ptr,"value","z")};
    
    b_param(doc.get());
    fe_param(doc.get());
    cre_param(doc.get());
}

// magnetic field
void Param::b_param(XMLDocument *doc){
    XMLElement *ptr {toolkit::tracexml(doc,{"MagneticField"})};
    std::string breg_type {toolkit::FetchString(ptr,"type","Regular")};
    // bwmap
    if(breg_type=="WMAP"){
        XMLElement *subptr {toolkit::tracexml(doc,{"MagneticField","Regular","WMAP"})};
        breg_wmap.b0 = toolkit::FetchDouble(subptr,"value","b0")*CGS_U_muGauss; //microGauss
        breg_wmap.psi0 = toolkit::FetchDouble(subptr,"value","psi0")*CGS_U_rad; //rad
        breg_wmap.psi1 = toolkit::FetchDouble(subptr,"value","psi1")*CGS_U_rad; //rad
        breg_wmap.chi0 = toolkit::FetchDouble(subptr,"value","chi0")*CGS_U_rad; //rad
    }
    // bjaffe
    else if(breg_type=="Jaffe"){
        XMLElement *subptr {toolkit::tracexml(doc,{"MagneticField","Regular","Jaffe"})};
        breg_jaffe.quadruple = toolkit::FetchBool(subptr,"cue","quadruple");
        breg_jaffe.bss = toolkit::FetchBool(subptr,"cue","bss");
        breg_jaffe.disk_amp = toolkit::FetchDouble(subptr,"value","disk_amp")*CGS_U_muGauss; //microG
        breg_jaffe.disk_z0 = toolkit::FetchDouble(subptr,"value","disk_z0")*CGS_U_kpc; //kpc
        breg_jaffe.halo_amp = toolkit::FetchDouble(subptr,"value","halo_amp")*CGS_U_muGauss; //microG
        breg_jaffe.halo_z0 = toolkit::FetchDouble(subptr,"value","halo_z0")*CGS_U_kpc; //kpc
        breg_jaffe.r_inner = toolkit::FetchDouble(subptr,"value","r_inner")*CGS_U_kpc; //kpc
        breg_jaffe.r_scale = toolkit::FetchDouble(subptr,"value","r_scale")*CGS_U_kpc; //kpc
        breg_jaffe.r_peak = toolkit::FetchDouble(subptr,"value","r_peak")*CGS_U_kpc; //kpc
        breg_jaffe.ring = toolkit::FetchBool(subptr,"cue","ring");
        breg_jaffe.ring_amp = toolkit::FetchDouble(subptr,"value","ring_amp")*CGS_U_muGauss; //microG
        breg_jaffe.ring_r = toolkit::FetchDouble(subptr,"value","ring_r")*CGS_U_kpc; //kpc
        breg_jaffe.bar = toolkit::FetchBool(subptr,"cue","bar");
        breg_jaffe.bar_amp = toolkit::FetchDouble(subptr,"value","bar_amp")*CGS_U_muGauss; //microG
        breg_jaffe.bar_a = toolkit::FetchDouble(subptr,"value","bar_a")*CGS_U_kpc; //kpc
        breg_jaffe.bar_b = toolkit::FetchDouble(subptr,"value","bar_b")*CGS_U_kpc; //kpc
        breg_jaffe.bar_phi0 = toolkit::FetchDouble(subptr,"value","bar_phi0")*CGS_U_rad; //rad
        breg_jaffe.arm_num = toolkit::FetchUnsigned(subptr,"value","arm_num");
        breg_jaffe.arm_r0 = toolkit::FetchDouble(subptr,"value","arm_r0")*CGS_U_kpc; //kpc
        breg_jaffe.arm_z0 = toolkit::FetchDouble(subptr,"value","arm_z0")*CGS_U_kpc; //kpc
        breg_jaffe.arm_phi0.push_back(toolkit::FetchDouble(subptr,"value","arm_phi1")*CGS_U_rad); //rad
        breg_jaffe.arm_phi0.push_back(toolkit::FetchDouble(subptr,"value","arm_phi2")*CGS_U_rad); //rad
        breg_jaffe.arm_phi0.push_back(toolkit::FetchDouble(subptr,"value","arm_phi3")*CGS_U_rad); //rad
        breg_jaffe.arm_phi0.push_back(toolkit::FetchDouble(subptr,"value","arm_phi4")*CGS_U_rad); //rad
        breg_jaffe.arm_amp.push_back(toolkit::FetchDouble(subptr,"value","arm_amp1")*CGS_U_muGauss); //microG
        breg_jaffe.arm_amp.push_back(toolkit::FetchDouble(subptr,"value","arm_amp2")*CGS_U_muGauss); //microG
        breg_jaffe.arm_amp.push_back(toolkit::FetchDouble(subptr,"value","arm_amp3")*CGS_U_muGauss); //microG
        breg_jaffe.arm_amp.push_back(toolkit::FetchDouble(subptr,"value","arm_amp4")*CGS_U_muGauss); //microG
        breg_jaffe.arm_pitch = toolkit::FetchDouble(subptr,"value","arm_pitch")*CGS_U_rad; //rad
        breg_jaffe.comp_c = toolkit::FetchDouble(subptr,"value","comp_c");
        breg_jaffe.comp_d = toolkit::FetchDouble(subptr,"value","comp_d")*CGS_U_kpc; //kpc
        breg_jaffe.comp_r = toolkit::FetchDouble(subptr,"value","comp_r")*CGS_U_kpc; //kpc
        breg_jaffe.comp_p = toolkit::FetchDouble(subptr,"value","comp_p");
    }
    // bverify
    else if(breg_type=="Verify"){
        XMLElement *subptr {toolkit::tracexml(doc,{"MagneticField","Regular","Verify"})};
        breg_verify.b0 = toolkit::FetchDouble(subptr,"value","b0")*CGS_U_muGauss; //microGauss
        breg_verify.l0 = toolkit::FetchDouble(subptr,"value","l0")*CGS_U_rad; //rad
        breg_verify.r = toolkit::FetchDouble(subptr,"value","r");
    }
    
    if(toolkit::FetchBool(ptr,"cue","Random")){
        // random seed
        brnd_seed = toolkit::FetchUnsigned(ptr,"seed","Random");
        std::string brnd_type {toolkit::FetchString(ptr,"type","Random")};
        // brnd_global
        if(brnd_type=="Global"){
            XMLElement *subptr {toolkit::tracexml(doc,{"MagneticField","Random","Global"})};
            std::string brnd_method {toolkit::FetchString(subptr,"type")};
            if(brnd_method=="ES"){
                subptr = toolkit::tracexml(doc,{"MagneticField","Random","Global","ES"});
                brnd_es.rms = toolkit::FetchDouble(subptr,"value","rms")*CGS_U_muGauss;
                brnd_es.k0 = toolkit::FetchDouble(subptr,"value","k0");
                brnd_es.a0 = toolkit::FetchDouble(subptr,"value","a0");
                brnd_es.rho = toolkit::FetchDouble(subptr,"value","rho");
                brnd_es.r0 = toolkit::FetchDouble(subptr,"value","r0")*CGS_U_kpc;
                brnd_es.z0 = toolkit::FetchDouble(subptr,"value","z0")*CGS_U_kpc;
            }
            else if(brnd_method=="Jaffe"){
                subptr = toolkit::tracexml(doc,{"MagneticField","Random","Global","Jaffe"});
                // to be implemented
            }
        }
        // brnd_local
        else if(brnd_type=="Local"){
            XMLElement *subptr {toolkit::tracexml(doc,{"MagneticField","Random","Local"})};
            std::string brnd_method {toolkit::FetchString(subptr,"type")};
            if(brnd_method=="MHD"){
                subptr = toolkit::tracexml(doc,{"MagneticField","Random","Local","MHD"});
                brnd_mhd.pa0 = toolkit::FetchDouble(subptr,"value","pa0")*CGS_U_muGauss*CGS_U_muGauss;
                brnd_mhd.pf0 = toolkit::FetchDouble(subptr,"value","pf0")*CGS_U_muGauss*CGS_U_muGauss;
                brnd_mhd.ps0 = toolkit::FetchDouble(subptr,"value","ps0")*CGS_U_muGauss*CGS_U_muGauss;
                brnd_mhd.aa0 = toolkit::FetchDouble(subptr,"value","aa0");
                brnd_mhd.af0 = toolkit::FetchDouble(subptr,"value","af0");
                brnd_mhd.as0 = toolkit::FetchDouble(subptr,"value","as0");
                brnd_mhd.k0 = toolkit::FetchDouble(subptr,"value","k0");
                brnd_mhd.beta = toolkit::FetchDouble(subptr,"value","beta");
                brnd_mhd.ma = toolkit::FetchDouble(subptr,"value","ma");
            }
        }
    }
}

void Param::fe_param(XMLDocument *doc){
    XMLElement *ptr {toolkit::tracexml(doc,{"FreeElectron"})};
    std::string fereg_type {toolkit::FetchString(ptr,"type","Regular")};
    // YMW16
    if(fereg_type=="YMW16"){
        // Warp_Sun
        XMLElement *subptr {toolkit::tracexml(doc,{"FreeElectron","Regular","YMW16","Warp"})};
        fereg_ymw16.R_warp = toolkit::FetchDouble(subptr,"value","R_warp")*CGS_U_kpc; //kpc
        fereg_ymw16.R0 = toolkit::FetchDouble(subptr,"value","R0")*CGS_U_kpc; //kpc
        fereg_ymw16.t0_Gamma_w = toolkit::FetchDouble(subptr,"value","Gamma_w");
        // thick disk
        subptr = toolkit::tracexml(doc,{"FreeElectron","Regular","YMW16","ThickDisk"});
        fereg_ymw16.t1_Ad = toolkit::FetchDouble(subptr,"value","Ad")*CGS_U_pc;//pc
        fereg_ymw16.t1_Bd = toolkit::FetchDouble(subptr,"value","Bd")*CGS_U_pc;//pc
        fereg_ymw16.t1_n1 = toolkit::FetchDouble(subptr,"value","n1");//pccm
        fereg_ymw16.t1_H1 = toolkit::FetchDouble(subptr,"value","H1")*CGS_U_pc;//pc
        // thin disk
        subptr = toolkit::tracexml(doc,{"FreeElectron","Regular","YMW16","ThinDisk"});
        fereg_ymw16.t2_A2 = toolkit::FetchDouble(subptr,"value","A2")*CGS_U_pc;//pc
        fereg_ymw16.t2_B2 = toolkit::FetchDouble(subptr,"value","B2")*CGS_U_pc;//pc
        fereg_ymw16.t2_n2 = toolkit::FetchDouble(subptr,"value","n2");//pccm
        fereg_ymw16.t2_K2 = toolkit::FetchDouble(subptr,"value","K2");
        // spiral arm
        subptr = toolkit::tracexml(doc,{"FreeElectron","Regular","YMW16","SpiralArm"});
        fereg_ymw16.t3_B2s = toolkit::FetchDouble(subptr,"value","B2s")*CGS_U_pc; //pc
        fereg_ymw16.t3_narm[0] = toolkit::FetchDouble(subptr,"value","Ele_arm_0");//pccm
        fereg_ymw16.t3_narm[1] = toolkit::FetchDouble(subptr,"value","Ele_arm_1");
        fereg_ymw16.t3_narm[2] = toolkit::FetchDouble(subptr,"value","Ele_arm_2");
        fereg_ymw16.t3_narm[3] = toolkit::FetchDouble(subptr,"value","Ele_arm_3");
        fereg_ymw16.t3_narm[4] = toolkit::FetchDouble(subptr,"value","Ele_arm_4");
        fereg_ymw16.t3_warm[0] = toolkit::FetchDouble(subptr,"value","Wid_arm_0")*CGS_U_pc;//pc
        fereg_ymw16.t3_warm[1] = toolkit::FetchDouble(subptr,"value","Wid_arm_1")*CGS_U_pc;
        fereg_ymw16.t3_warm[2] = toolkit::FetchDouble(subptr,"value","Wid_arm_2")*CGS_U_pc;
        fereg_ymw16.t3_warm[3] = toolkit::FetchDouble(subptr,"value","Wid_arm_3")*CGS_U_pc;
        fereg_ymw16.t3_warm[4] = toolkit::FetchDouble(subptr,"value","Wid_arm_4")*CGS_U_pc;
        fereg_ymw16.t3_rmin[0] = toolkit::FetchDouble(subptr,"value","Rref_arm_0")*CGS_U_kpc;//kpc
        fereg_ymw16.t3_rmin[1] = toolkit::FetchDouble(subptr,"value","Rref_arm_1")*CGS_U_kpc;
        fereg_ymw16.t3_rmin[2] = toolkit::FetchDouble(subptr,"value","Rref_arm_2")*CGS_U_kpc;
        fereg_ymw16.t3_rmin[3] = toolkit::FetchDouble(subptr,"value","Rref_arm_3")*CGS_U_kpc;
        fereg_ymw16.t3_rmin[4] = toolkit::FetchDouble(subptr,"value","Rref_arm_4")*CGS_U_kpc;
        fereg_ymw16.t3_phimin[0] = toolkit::FetchDouble(subptr,"value","Phiref_arm_0")*CGS_U_rad;//rad
        fereg_ymw16.t3_phimin[1] = toolkit::FetchDouble(subptr,"value","Phiref_arm_1")*CGS_U_rad;//rad
        fereg_ymw16.t3_phimin[2] = toolkit::FetchDouble(subptr,"value","Phiref_arm_2")*CGS_U_rad;//rad
        fereg_ymw16.t3_phimin[3] = toolkit::FetchDouble(subptr,"value","Phiref_arm_3")*CGS_U_rad;//rad
        fereg_ymw16.t3_phimin[4] = toolkit::FetchDouble(subptr,"value","Phiref_arm_4")*CGS_U_rad;//rad
        fereg_ymw16.t3_tpitch[0] = tan(toolkit::FetchDouble(subptr,"value","pitch_arm_0")*CGS_U_rad);
        fereg_ymw16.t3_tpitch[1] = tan(toolkit::FetchDouble(subptr,"value","pitch_arm_1")*CGS_U_rad);
        fereg_ymw16.t3_tpitch[2] = tan(toolkit::FetchDouble(subptr,"value","pitch_arm_2")*CGS_U_rad);
        fereg_ymw16.t3_tpitch[3] = tan(toolkit::FetchDouble(subptr,"value","pitch_arm_3")*CGS_U_rad);
        fereg_ymw16.t3_tpitch[4] = tan(toolkit::FetchDouble(subptr,"value","pitch_arm_4")*CGS_U_rad);
        fereg_ymw16.t3_cpitch[0] = cos(toolkit::FetchDouble(subptr,"value","pitch_arm_0")*CGS_U_rad);
        fereg_ymw16.t3_cpitch[1] = cos(toolkit::FetchDouble(subptr,"value","pitch_arm_1")*CGS_U_rad);
        fereg_ymw16.t3_cpitch[2] = cos(toolkit::FetchDouble(subptr,"value","pitch_arm_2")*CGS_U_rad);
        fereg_ymw16.t3_cpitch[3] = cos(toolkit::FetchDouble(subptr,"value","pitch_arm_3")*CGS_U_rad);
        fereg_ymw16.t3_cpitch[4] = cos(toolkit::FetchDouble(subptr,"value","pitch_arm_4")*CGS_U_rad);
        fereg_ymw16.t3_Aa = toolkit::FetchDouble(subptr,"value","Aa")*CGS_U_pc;//pc
        fereg_ymw16.t3_Ka = toolkit::FetchDouble(subptr,"value","Ka");
        fereg_ymw16.t3_ncn = toolkit::FetchDouble(subptr,"value","ncn");
        fereg_ymw16.t3_thetacn = toolkit::FetchDouble(subptr,"value","thetacn");//deg
        fereg_ymw16.t3_wcn = toolkit::FetchDouble(subptr,"value","wcn");//deg
        fereg_ymw16.t3_nsg = toolkit::FetchDouble(subptr,"value","nsg");
        fereg_ymw16.t3_thetasg = toolkit::FetchDouble(subptr,"value","thetasg");//deg
        fereg_ymw16.t3_wsg = toolkit::FetchDouble(subptr,"value","wsg");//deg
        // gc
        subptr = toolkit::tracexml(doc,{"FreeElectron","Regular","YMW16","GalCenter"});
        fereg_ymw16.t4_ngc = toolkit::FetchDouble(subptr,"value","ngc");//pccm
        fereg_ymw16.t4_Agc = toolkit::FetchDouble(subptr,"value","Agc")*CGS_U_pc;//pc
        fereg_ymw16.t4_Hgc = toolkit::FetchDouble(subptr,"value","Hgc")*CGS_U_pc;//pc
        // Gum
        subptr = toolkit::tracexml(doc,{"FreeElectron","Regular","YMW16","GumNebula"});
        fereg_ymw16.t5_ngn = toolkit::FetchDouble(subptr,"value","ngn");//pccm
        fereg_ymw16.t5_Wgn = toolkit::FetchDouble(subptr,"value","Wgn")*CGS_U_pc;//pc
        fereg_ymw16.t5_Agn = toolkit::FetchDouble(subptr,"value","Agn")*CGS_U_pc;//pc
        fereg_ymw16.t5_Kgn = toolkit::FetchDouble(subptr,"value","Kgn");
        // Local Bubble
        subptr = toolkit::tracexml(doc,{"FreeElectron","Regular","YMW16","LocalBubble"});
        fereg_ymw16.t6_J_LB = toolkit::FetchDouble(subptr,"value","J_LB");
        fereg_ymw16.t6_nlb1 = toolkit::FetchDouble(subptr,"value","nlb1");//pccm
        fereg_ymw16.t6_thetalb1 = toolkit::FetchDouble(subptr,"value","thetalb1");//deg
        fereg_ymw16.t6_detlb1 = toolkit::FetchDouble(subptr,"value","detlb1");//deg
        fereg_ymw16.t6_wlb1 = toolkit::FetchDouble(subptr,"value","wlb1")*CGS_U_pc;//pc
        fereg_ymw16.t6_hlb1 = toolkit::FetchDouble(subptr,"value","hlb1")*CGS_U_pc;//pc
        fereg_ymw16.t6_nlb2 = toolkit::FetchDouble(subptr,"value","nlb2");//pccm
        fereg_ymw16.t6_thetalb2 = toolkit::FetchDouble(subptr,"value","thetalb2");//deg
        fereg_ymw16.t6_detlb2 = toolkit::FetchDouble(subptr,"value","detlb2");//deg
        fereg_ymw16.t6_wlb2 = toolkit::FetchDouble(subptr,"value","wlb2")*CGS_U_pc;//pc
        fereg_ymw16.t6_hlb2 = toolkit::FetchDouble(subptr,"value","hlb2")*CGS_U_pc;//pc
        // Loop I
        subptr = toolkit::tracexml(doc,{"FreeElectron","Regular","YMW16","LoopI"});
        fereg_ymw16.t7_nLI = toolkit::FetchDouble(subptr,"value","nLI");//pccm
        fereg_ymw16.t7_RLI = toolkit::FetchDouble(subptr,"value","RLI")*CGS_U_pc;//pc
        fereg_ymw16.t7_WLI = toolkit::FetchDouble(subptr,"value","WLI")*CGS_U_pc;//pc
        fereg_ymw16.t7_detthetaLI = toolkit::FetchDouble(subptr,"value","detthetaLI");//deg
        fereg_ymw16.t7_thetaLI = toolkit::FetchDouble(subptr,"value","thetaLI");//deg
    }
    // verify
    else if(fereg_type=="Verify"){
        XMLElement *subptr {toolkit::tracexml(doc,{"FreeElectron","Regular","Verify"})};
        fereg_verify.n0 = toolkit::FetchDouble(subptr,"value","n0");
        fereg_verify.r0 = toolkit::FetchDouble(subptr,"value","r0")*CGS_U_kpc; //kpc
    }
    
    if(toolkit::FetchBool(ptr,"cue","Random")){
        // random seed
        fernd_seed = toolkit::FetchUnsigned(ptr,"seed","Random");
        std::string fernd_type {toolkit::FetchString(ptr,"type","Random")};
        // global turbulent
        if(fernd_type=="Global"){
            XMLElement *subptr {toolkit::tracexml(doc,{"FreeElectron","Random","Global"})};
            std::string fernd_method {toolkit::FetchString(subptr,"type")};
            if(fernd_method=="DFT"){
                subptr = toolkit::tracexml(doc,{"FreeElectron","Random","Global","DFT"});
                fernd_dft.rms = toolkit::FetchDouble(subptr,"value","rms");
                fernd_dft.k0 = toolkit::FetchDouble(subptr,"value","k0");
                fernd_dft.a0 = toolkit::FetchDouble(subptr,"value","a0");
                fernd_dft.r0 = toolkit::FetchDouble(subptr,"value","r0")*CGS_U_kpc;
                fernd_dft.z0 = toolkit::FetchDouble(subptr,"value","z0")*CGS_U_kpc;
            }
        }
    }
}

void Param::cre_param(XMLDocument *doc){
    XMLElement *ptr {toolkit::tracexml(doc,{"CRE"})};
    if (toolkit::tracexml(doc,{"Obsout","Sync"})!=nullptr){
        XMLElement *subptr {toolkit::tracexml(doc,{"Obsout"})};
        sim_freq = toolkit::FetchDouble(subptr,"freq","Sync")*CGS_U_GHz;
    }
    else{
        sim_freq = 0.;
    }
    std::string cre_type {ptr->Attribute("type")};
    // analytical
    if(cre_type=="Analytic"){
        XMLElement *subptr {toolkit::tracexml(doc,{"CRE","Analytic"})};
        cre_ana.alpha = toolkit::FetchDouble(subptr,"value","alpha");
        cre_ana.beta = toolkit::FetchDouble(subptr,"value","beta");
        cre_ana.theta = toolkit::FetchDouble(subptr,"value","theta");
        cre_ana.r0 = toolkit::FetchDouble(subptr,"value","r0")*CGS_U_kpc; //kpc
        cre_ana.z0 = toolkit::FetchDouble(subptr,"value","z0")*CGS_U_kpc; //kpc
        cre_ana.E0 = toolkit::FetchDouble(subptr,"value","E0")*CGS_U_GeV; //GeV
        cre_ana.j0 = toolkit::FetchDouble(subptr,"value","j0");
    }
    // verification
    else if(cre_type=="Verify"){
        XMLElement *subptr {toolkit::tracexml(doc,{"CRE","Verify"})};
        cre_verify.alpha = toolkit::FetchDouble(subptr,"value","alpha");
        cre_verify.r0 = toolkit::FetchDouble(subptr,"value","r0")*CGS_U_kpc; //kpc
        cre_verify.E0 = toolkit::FetchDouble(subptr,"value","E0")*CGS_U_GeV; //GeV
        cre_verify.j0 = toolkit::FetchDouble(subptr,"value","j0");
    }
}

// END
