#include <string>
#include <vec3.h>
#include <tinyxml2.h>
#include <param.h>
#include <cgs_units_file.h>
#include <namespace_toolkit.h>
#include <cmath>
using namespace tinyxml2;
using namespace std;

Param::Param(std::string file_name){
    XMLDocument *doc = new XMLDocument();
    doc->LoadFile(file_name.c_str());
    // gc sun position
    XMLElement *ptr {doc->FirstChildElement("root")->FirstChildElement("Grid")->FirstChildElement("SunPosition")};
    SunPosition = vec3_t<double> {CGS_U_kpc*toolkit::FetchDouble(ptr,"x"),
        CGS_U_kpc*toolkit::FetchDouble(ptr,"y"),
        CGS_U_pc*toolkit::FetchDouble(ptr,"z")};
    
    b_param(doc);
    fe_param(doc);
    cre_param(doc);
    delete doc;
}

// magnetic field
void Param::b_param(XMLDocument *doc){
    XMLElement *ptr {doc->FirstChildElement("root")->FirstChildElement("MagneticField")};
    string breg_type {ptr->FirstChildElement("Regular")->Attribute("type")};
    // bwmap
    if(breg_type=="WMAP"){
        XMLElement *subptr {ptr->FirstChildElement("Regular")->FirstChildElement("WMAP")};
        breg_wmap.b0 = toolkit::FetchDouble(subptr,"b0")*CGS_U_muGauss; //microGauss
        breg_wmap.psi0 = toolkit::FetchDouble(subptr,"psi0")*CGS_U_rad; //rad
        breg_wmap.psi1 = toolkit::FetchDouble(subptr,"psi1")*CGS_U_rad; //rad
        breg_wmap.chi0 = toolkit::FetchDouble(subptr,"chi0")*CGS_U_rad; //rad
    }
    // bjaffe
    else if(breg_type=="Jaffe"){
        XMLElement *subptr {ptr->FirstChildElement("Regular")->FirstChildElement("Jaffe")};
        breg_jaffe.quadruple = toolkit::FetchBool(subptr,"quadruple");
        breg_jaffe.bss = toolkit::FetchBool(subptr,"bss");
        breg_jaffe.disk_amp = toolkit::FetchDouble(subptr,"disk_amp")*CGS_U_muGauss; //microG
        breg_jaffe.disk_z0 = toolkit::FetchDouble(subptr,"disk_z0")*CGS_U_kpc; //kpc
        breg_jaffe.halo_amp = toolkit::FetchDouble(subptr,"halo_amp")*CGS_U_muGauss; //microG
        breg_jaffe.halo_z0 = toolkit::FetchDouble(subptr,"halo_z0")*CGS_U_kpc; //kpc
        breg_jaffe.r_inner = toolkit::FetchDouble(subptr,"r_inner")*CGS_U_kpc; //kpc
        breg_jaffe.r_scale = toolkit::FetchDouble(subptr,"r_scale")*CGS_U_kpc; //kpc
        breg_jaffe.r_peak = toolkit::FetchDouble(subptr,"r_peak")*CGS_U_kpc; //kpc
        breg_jaffe.ring = toolkit::FetchBool(subptr,"ring");
        breg_jaffe.ring_amp = toolkit::FetchDouble(subptr,"ring_amp")*CGS_U_muGauss; //microG
        breg_jaffe.ring_r = toolkit::FetchDouble(subptr,"ring_r")*CGS_U_kpc; //kpc
        breg_jaffe.bar = toolkit::FetchBool(subptr,"bar");
        breg_jaffe.bar_amp = toolkit::FetchDouble(subptr,"bar_amp")*CGS_U_muGauss; //microG
        breg_jaffe.bar_a = toolkit::FetchDouble(subptr,"bar_a")*CGS_U_kpc; //kpc
        breg_jaffe.bar_b = toolkit::FetchDouble(subptr,"bar_b")*CGS_U_kpc; //kpc
        breg_jaffe.bar_phi0 = toolkit::FetchDouble(subptr,"bar_phi0")*CGS_U_rad; //rad
        breg_jaffe.arm_num = toolkit::FetchUnsigned(subptr,"arm_num");
        breg_jaffe.arm_r0 = toolkit::FetchDouble(subptr,"arm_r0")*CGS_U_kpc; //kpc
        breg_jaffe.arm_z0 = toolkit::FetchDouble(subptr,"arm_z0")*CGS_U_kpc; //kpc
        breg_jaffe.arm_phi0.push_back(toolkit::FetchDouble(subptr,"arm_phi1")*CGS_U_rad); //rad
        breg_jaffe.arm_phi0.push_back(toolkit::FetchDouble(subptr,"arm_phi2")*CGS_U_rad); //rad
        breg_jaffe.arm_phi0.push_back(toolkit::FetchDouble(subptr,"arm_phi3")*CGS_U_rad); //rad
        breg_jaffe.arm_phi0.push_back(toolkit::FetchDouble(subptr,"arm_phi4")*CGS_U_rad); //rad
        breg_jaffe.arm_amp.push_back(toolkit::FetchDouble(subptr,"arm_amp1")*CGS_U_muGauss); //microG
        breg_jaffe.arm_amp.push_back(toolkit::FetchDouble(subptr,"arm_amp2")*CGS_U_muGauss); //microG
        breg_jaffe.arm_amp.push_back(toolkit::FetchDouble(subptr,"arm_amp3")*CGS_U_muGauss); //microG
        breg_jaffe.arm_amp.push_back(toolkit::FetchDouble(subptr,"arm_amp4")*CGS_U_muGauss); //microG
        breg_jaffe.arm_pitch = toolkit::FetchDouble(subptr,"arm_pitch")*CGS_U_rad; //rad
        breg_jaffe.comp_c = toolkit::FetchDouble(subptr,"comp_c");
        breg_jaffe.comp_d = toolkit::FetchDouble(subptr,"comp_d")*CGS_U_kpc; //kpc
        breg_jaffe.comp_r = toolkit::FetchDouble(subptr,"comp_r")*CGS_U_kpc; //kpc
        breg_jaffe.comp_p = toolkit::FetchDouble(subptr,"comp_p");
    }
    // bverify
    else if(breg_type=="Verify"){
        XMLElement *subptr {ptr->FirstChildElement("Regular")->FirstChildElement("Verify")};
        breg_verify.b0 = toolkit::FetchDouble(subptr,"b0")*CGS_U_muGauss; //microGauss
        breg_verify.l0 = toolkit::FetchDouble(subptr,"l0")*CGS_U_rad; //rad
        breg_verify.r = toolkit::FetchDouble(subptr,"r");
    }
    
    if(ptr->FirstChildElement("Random")->BoolAttribute("cue")){
        // random seed
        brnd_seed = ptr->FirstChildElement("Random")->UnsignedAttribute("seed");
        string brnd_type {ptr->FirstChildElement("Random")->Attribute("type")};
        // brnd_global
        if(brnd_type=="Global"){
            XMLElement *subptr {ptr->FirstChildElement("Random")->FirstChildElement("Global")};
            brnd_global.rms = toolkit::FetchDouble(subptr,"rms")*CGS_U_muGauss;
            brnd_global.k0 = toolkit::FetchDouble(subptr,"k0");
            brnd_global.a0 = toolkit::FetchDouble(subptr,"a0");
            brnd_global.rho = toolkit::FetchDouble(subptr,"rho");
            brnd_global.r0 = toolkit::FetchDouble(subptr,"r0")*CGS_U_kpc;
            brnd_global.z0 = toolkit::FetchDouble(subptr,"z0")*CGS_U_kpc;
        }
        // brnd_local
        else if(brnd_type=="Local"){
            XMLElement *subptr {ptr->FirstChildElement("Random")->FirstChildElement("Local")};
            brnd_local.pa0 = toolkit::FetchDouble(subptr,"pa0")*CGS_U_muGauss*CGS_U_muGauss;
            brnd_local.pf0 = toolkit::FetchDouble(subptr,"pf0")*CGS_U_muGauss*CGS_U_muGauss;
            brnd_local.ps0 = toolkit::FetchDouble(subptr,"ps0")*CGS_U_muGauss*CGS_U_muGauss;
            brnd_local.aa0 = toolkit::FetchDouble(subptr,"aa0");
            brnd_local.af0 = toolkit::FetchDouble(subptr,"af0");
            brnd_local.as0 = toolkit::FetchDouble(subptr,"as0");
            brnd_local.k0 = toolkit::FetchDouble(subptr,"k0");
            brnd_local.beta = toolkit::FetchDouble(subptr,"beta");
            brnd_local.ma = toolkit::FetchDouble(subptr,"ma");
        }
    }
}

void Param::fe_param(XMLDocument *doc){
    XMLElement *ptr {doc->FirstChildElement("root")->FirstChildElement("FreeElectron")};
    string fereg_type {ptr->FirstChildElement("Regular")->Attribute("type")};
    // YMW16
    if(fereg_type=="YMW16"){
        // Warp_Sun
        XMLElement *subptr {ptr->FirstChildElement("Regular")->FirstChildElement("YMW16")->FirstChildElement("Warp")};
        fereg_ymw16.R_warp = toolkit::FetchDouble(subptr,"R_warp")*CGS_U_kpc; //kpc
        fereg_ymw16.R0 = toolkit::FetchDouble(subptr,"R0")*CGS_U_kpc; //kpc
        fereg_ymw16.t0_Gamma_w = toolkit::FetchDouble(subptr,"Gamma_w");
        // thick disk
        subptr = ptr->FirstChildElement("Regular")->FirstChildElement("YMW16")->FirstChildElement("ThickDisk");
        fereg_ymw16.t1_Ad = toolkit::FetchDouble(subptr,"Ad")*CGS_U_pc;//pc
        fereg_ymw16.t1_Bd = toolkit::FetchDouble(subptr,"Bd")*CGS_U_pc;//pc
        fereg_ymw16.t1_n1 = toolkit::FetchDouble(subptr,"n1");//pccm
        fereg_ymw16.t1_H1 = toolkit::FetchDouble(subptr,"H1")*CGS_U_pc;//pc
        // thin disk
        subptr = ptr->FirstChildElement("Regular")->FirstChildElement("YMW16")->FirstChildElement("ThinDisk");
        fereg_ymw16.t2_A2 = toolkit::FetchDouble(subptr,"A2")*CGS_U_pc;//pc
        fereg_ymw16.t2_B2 = toolkit::FetchDouble(subptr,"B2")*CGS_U_pc;//pc
        fereg_ymw16.t2_n2 = toolkit::FetchDouble(subptr,"n2");//pccm
        fereg_ymw16.t2_K2 = toolkit::FetchDouble(subptr,"K2");
        // spiral arm
        subptr = ptr->FirstChildElement("Regular")->FirstChildElement("YMW16")->FirstChildElement("SpiralArm");
        fereg_ymw16.t3_B2s = toolkit::FetchDouble(subptr,"B2s")*CGS_U_pc; //pc
        fereg_ymw16.t3_narm[0] = toolkit::FetchDouble(subptr,"Ele_arm_0");//pccm
        fereg_ymw16.t3_narm[1] = toolkit::FetchDouble(subptr,"Ele_arm_1");
        fereg_ymw16.t3_narm[2] = toolkit::FetchDouble(subptr,"Ele_arm_2");
        fereg_ymw16.t3_narm[3] = toolkit::FetchDouble(subptr,"Ele_arm_3");
        fereg_ymw16.t3_narm[4] = toolkit::FetchDouble(subptr,"Ele_arm_4");
        fereg_ymw16.t3_warm[0] = toolkit::FetchDouble(subptr,"Wid_arm_0")*CGS_U_pc;//pc
        fereg_ymw16.t3_warm[1] = toolkit::FetchDouble(subptr,"Wid_arm_1")*CGS_U_pc;
        fereg_ymw16.t3_warm[2] = toolkit::FetchDouble(subptr,"Wid_arm_2")*CGS_U_pc;
        fereg_ymw16.t3_warm[3] = toolkit::FetchDouble(subptr,"Wid_arm_3")*CGS_U_pc;
        fereg_ymw16.t3_warm[4] = toolkit::FetchDouble(subptr,"Wid_arm_4")*CGS_U_pc;
        fereg_ymw16.t3_rmin[0] = toolkit::FetchDouble(subptr,"Rref_arm_0")*CGS_U_kpc;//kpc
        fereg_ymw16.t3_rmin[1] = toolkit::FetchDouble(subptr,"Rref_arm_1")*CGS_U_kpc;
        fereg_ymw16.t3_rmin[2] = toolkit::FetchDouble(subptr,"Rref_arm_2")*CGS_U_kpc;
        fereg_ymw16.t3_rmin[3] = toolkit::FetchDouble(subptr,"Rref_arm_3")*CGS_U_kpc;
        fereg_ymw16.t3_rmin[4] = toolkit::FetchDouble(subptr,"Rref_arm_4")*CGS_U_kpc;
        fereg_ymw16.t3_phimin[0] = toolkit::FetchDouble(subptr,"Phiref_arm_0")*CGS_U_rad;//rad
        fereg_ymw16.t3_phimin[1] = toolkit::FetchDouble(subptr,"Phiref_arm_1")*CGS_U_rad;//rad
        fereg_ymw16.t3_phimin[2] = toolkit::FetchDouble(subptr,"Phiref_arm_2")*CGS_U_rad;//rad
        fereg_ymw16.t3_phimin[3] = toolkit::FetchDouble(subptr,"Phiref_arm_3")*CGS_U_rad;//rad
        fereg_ymw16.t3_phimin[4] = toolkit::FetchDouble(subptr,"Phiref_arm_4")*CGS_U_rad;//rad
        fereg_ymw16.t3_tpitch[0] = tan(toolkit::FetchDouble(subptr,"pitch_arm_0")*CGS_U_rad);
        fereg_ymw16.t3_tpitch[1] = tan(toolkit::FetchDouble(subptr,"pitch_arm_1")*CGS_U_rad);
        fereg_ymw16.t3_tpitch[2] = tan(toolkit::FetchDouble(subptr,"pitch_arm_2")*CGS_U_rad);
        fereg_ymw16.t3_tpitch[3] = tan(toolkit::FetchDouble(subptr,"pitch_arm_3")*CGS_U_rad);
        fereg_ymw16.t3_tpitch[4] = tan(toolkit::FetchDouble(subptr,"pitch_arm_4")*CGS_U_rad);
        fereg_ymw16.t3_cpitch[0] = cos(toolkit::FetchDouble(subptr,"pitch_arm_0")*CGS_U_rad);
        fereg_ymw16.t3_cpitch[1] = cos(toolkit::FetchDouble(subptr,"pitch_arm_1")*CGS_U_rad);
        fereg_ymw16.t3_cpitch[2] = cos(toolkit::FetchDouble(subptr,"pitch_arm_2")*CGS_U_rad);
        fereg_ymw16.t3_cpitch[3] = cos(toolkit::FetchDouble(subptr,"pitch_arm_3")*CGS_U_rad);
        fereg_ymw16.t3_cpitch[4] = cos(toolkit::FetchDouble(subptr,"pitch_arm_4")*CGS_U_rad);
        fereg_ymw16.t3_Aa = toolkit::FetchDouble(subptr,"Aa")*CGS_U_pc;//pc
        fereg_ymw16.t3_Ka = toolkit::FetchDouble(subptr,"Ka");
        fereg_ymw16.t3_ncn = toolkit::FetchDouble(subptr,"ncn");
        fereg_ymw16.t3_thetacn = toolkit::FetchDouble(subptr,"thetacn");//deg
        fereg_ymw16.t3_wcn = toolkit::FetchDouble(subptr,"wcn");//deg
        fereg_ymw16.t3_nsg = toolkit::FetchDouble(subptr,"nsg");
        fereg_ymw16.t3_thetasg = toolkit::FetchDouble(subptr,"thetasg");//deg
        fereg_ymw16.t3_wsg = toolkit::FetchDouble(subptr,"wsg");//deg
        // gc
        subptr = ptr->FirstChildElement("Regular")->FirstChildElement("YMW16")->FirstChildElement("GalCenter");
        fereg_ymw16.t4_ngc = toolkit::FetchDouble(subptr,"ngc");//pccm
        fereg_ymw16.t4_Agc = toolkit::FetchDouble(subptr,"Agc")*CGS_U_pc;//pc
        fereg_ymw16.t4_Hgc = toolkit::FetchDouble(subptr,"Hgc")*CGS_U_pc;//pc
        // Gum
        subptr = ptr->FirstChildElement("Regular")->FirstChildElement("YMW16")->FirstChildElement("GumNebula");
        fereg_ymw16.t5_ngn = toolkit::FetchDouble(subptr,"ngn");//pccm
        fereg_ymw16.t5_Wgn = toolkit::FetchDouble(subptr,"Wgn")*CGS_U_pc;//pc
        fereg_ymw16.t5_Agn = toolkit::FetchDouble(subptr,"Agn")*CGS_U_pc;//pc
        fereg_ymw16.t5_Kgn = toolkit::FetchDouble(subptr,"Kgn");
        // Local Bubble
        subptr = ptr->FirstChildElement("Regular")->FirstChildElement("YMW16")->FirstChildElement("LocalBubble");
        fereg_ymw16.t6_J_LB = toolkit::FetchDouble(subptr,"J_LB");
        fereg_ymw16.t6_nlb1 = toolkit::FetchDouble(subptr,"nlb1");//pccm
        fereg_ymw16.t6_thetalb1 = toolkit::FetchDouble(subptr,"thetalb1");//deg
        fereg_ymw16.t6_detlb1 = toolkit::FetchDouble(subptr,"detlb1");//deg
        fereg_ymw16.t6_wlb1 = toolkit::FetchDouble(subptr,"wlb1")*CGS_U_pc;//pc
        fereg_ymw16.t6_hlb1 = toolkit::FetchDouble(subptr,"hlb1")*CGS_U_pc;//pc
        fereg_ymw16.t6_nlb2 = toolkit::FetchDouble(subptr,"nlb2");//pccm
        fereg_ymw16.t6_thetalb2 = toolkit::FetchDouble(subptr,"thetalb2");//deg
        fereg_ymw16.t6_detlb2 = toolkit::FetchDouble(subptr,"detlb2");//deg
        fereg_ymw16.t6_wlb2 = toolkit::FetchDouble(subptr,"wlb2")*CGS_U_pc;//pc
        fereg_ymw16.t6_hlb2 = toolkit::FetchDouble(subptr,"hlb2")*CGS_U_pc;//pc
        // Loop I
        subptr = ptr->FirstChildElement("Regular")->FirstChildElement("YMW16")->FirstChildElement("LoopI");
        fereg_ymw16.t7_nLI = toolkit::FetchDouble(subptr,"nLI");//pccm
        fereg_ymw16.t7_RLI = toolkit::FetchDouble(subptr,"RLI")*CGS_U_pc;//pc
        fereg_ymw16.t7_WLI = toolkit::FetchDouble(subptr,"WLI")*CGS_U_pc;//pc
        fereg_ymw16.t7_detthetaLI = toolkit::FetchDouble(subptr,"detthetaLI");//deg
        fereg_ymw16.t7_thetaLI = toolkit::FetchDouble(subptr,"thetaLI");//deg
    }
    // verify
    else if(fereg_type=="Verify"){
        XMLElement *subptr {ptr->FirstChildElement("Regular")->FirstChildElement("Verify")};
        fereg_verify.n0 = toolkit::FetchDouble(subptr,"n0");
        fereg_verify.r0 = toolkit::FetchDouble(subptr,"r0")*CGS_U_kpc; //kpc
    }
    
    if(ptr->FirstChildElement("Random")->BoolAttribute("cue")){
        // random seed
        fernd_seed = ptr->FirstChildElement("Random")->UnsignedAttribute("seed");
        string fernd_type {ptr->FirstChildElement("Random")->Attribute("type")};
        // global turbulent
        if(fernd_type=="Global"){
            XMLElement *subptr {ptr->FirstChildElement("Random")->FirstChildElement("Global")};
            fernd_global.rms = toolkit::FetchDouble(subptr,"rms");
            fernd_global.k0 = toolkit::FetchDouble(subptr,"k0");
            fernd_global.a0 = toolkit::FetchDouble(subptr,"a0");
            fernd_global.r0 = toolkit::FetchDouble(subptr,"r0")*CGS_U_kpc;
            fernd_global.z0 = toolkit::FetchDouble(subptr,"z0")*CGS_U_kpc;
        }
    }
}

void Param::cre_param(XMLDocument *doc){
    XMLElement *ptr {doc->FirstChildElement("root")->FirstChildElement("CRE")};
    if (doc->FirstChildElement("root")->FirstChildElement("Obsout")->FirstChildElement("Sync")!=nullptr){
        sim_freq = doc->FirstChildElement("root")->FirstChildElement("Obsout")->FirstChildElement("Sync")->DoubleAttribute("freq")*CGS_U_GHz;
    }
    else{
        sim_freq = 0.;
    }
    string cre_type {ptr->Attribute("type")};
    // analytical
    if(cre_type=="Analytic"){
        XMLElement *subptr {ptr->FirstChildElement("Analytic")};
        cre_ana.alpha = toolkit::FetchDouble(subptr,"alpha");
        cre_ana.beta = toolkit::FetchDouble(subptr,"beta");
        cre_ana.theta = toolkit::FetchDouble(subptr,"theta");
        cre_ana.r0 = toolkit::FetchDouble(subptr,"r0")*CGS_U_kpc; //kpc
        cre_ana.z0 = toolkit::FetchDouble(subptr,"z0")*CGS_U_kpc; //kpc
        cre_ana.E0 = toolkit::FetchDouble(subptr,"E0")*CGS_U_GeV; //GeV
        cre_ana.j0 = toolkit::FetchDouble(subptr,"j0");
    }
    // verification
    else if(cre_type=="Verify"){
        XMLElement *subptr {ptr->FirstChildElement("Verify")};
        cre_verify.alpha = toolkit::FetchDouble(subptr,"alpha");
        cre_verify.r0 = toolkit::FetchDouble(subptr,"r0")*CGS_U_kpc; //kpc
        cre_verify.E0 = toolkit::FetchDouble(subptr,"E0")*CGS_U_GeV; //GeV
        cre_verify.j0 = toolkit::FetchDouble(subptr,"j0");
    }
}

// END
