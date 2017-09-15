/* generic parameter hoster, using TinyXML2
 
 b field parameters
 analytical cre field parameters
 free electron parameters (of YMW16)
 
 */

#include <string>
#include <tinyxml2.h>
#include "class_pond.h"
#include "cgs_units_file.h"

using namespace tinyxml2;
using namespace std;

Pond::Pond(std::string file_name){
    XMLDocument *doc = new XMLDocument();
    doc->LoadFile(file_name.c_str());
    // gc sun position
    XMLElement *ptr = doc->FirstChildElement("root")->FirstChildElement("Grid")->FirstChildElement("SunPosition");
    SunPosition.x = CGS_U_kpc*FetchDouble(ptr,"x");
    SunPosition.y = CGS_U_kpc*FetchDouble(ptr,"y");
    SunPosition.z = CGS_U_pc*FetchDouble(ptr,"z");
    
    b_pond(doc);
    fe_pond(doc);
    cre_pond(doc);
}

// magnetic field
void Pond::b_pond(XMLDocument *doc){
    XMLElement *ptr = doc->FirstChildElement("root")->FirstChildElement("Galaxy")->FirstChildElement("MagneticField");
    // bwmap
    XMLElement *subptr = ptr->FirstChildElement("Regular")->FirstChildElement("WMAP");
    bwmap.push_back(FetchDouble(subptr,"b0")); //microGauss
    bwmap.push_back(FetchDouble(subptr,"psi0")); //deg
    bwmap.push_back(FetchDouble(subptr,"psi1")); //deg
    bwmap.push_back(FetchDouble(subptr,"chi0")); //deg
    // blocal
    subptr = ptr->FirstChildElement("Regular")->FirstChildElement("Local");
    blocal.push_back(FetchDouble(subptr,"bd")); //microGauss
    blocal.push_back(FetchDouble(subptr,"l0")); //deg
    blocal.push_back(FetchDouble(subptr,"z0")); //kpc
    blocal.push_back(FetchDouble(subptr,"bn")); //microGauss
    blocal.push_back(FetchDouble(subptr,"bs")); //microGauss
    // brnd_iso
    subptr = ptr->FirstChildElement("Random")->FirstChildElement("Iso");
    brnd_iso.push_back(FetchDouble(subptr,"rms")); //in microGauss^2;
    // turning point
    brnd_iso.push_back(FetchDouble(subptr,"k0"));
    // index before turning (large scales)
    brnd_iso.push_back(FetchDouble(subptr,"a0"));
    // brnd_aniso
    subptr = ptr->FirstChildElement("Random")->FirstChildElement("Anisoglob");
    brnd_ani.push_back(FetchDouble(subptr,"rho"));
    // rescaling parameters
    subptr = ptr->FirstChildElement("Random")->FirstChildElement("Rescal");
    brnd_scal.push_back(FetchDouble(subptr,"r0"));
    brnd_scal.push_back(FetchDouble(subptr,"z0"));
}

void Pond::fe_pond(XMLDocument *doc){
    XMLElement *ptr = doc->FirstChildElement("root")->FirstChildElement("Galaxy")->FirstChildElement("FreeElectron");
    // YMW16
    // Warp_Sun
    XMLElement *subptr = ptr->FirstChildElement("YMW16")->FirstChildElement("Warp");
    R_warp = FetchDouble(subptr,"R_warp")*CGS_U_kpc; //kpc
    R0 = FetchDouble(subptr,"R0")*CGS_U_kpc; //kpc
    t0_Gamma_w = FetchDouble(subptr,"Gamma_w");
    t0_z_Sun = SunPosition.z; //kpc
    // thick disk
    subptr = ptr->FirstChildElement("YMW16")->FirstChildElement("ThickDisk");
    t1_Ad = FetchDouble(subptr,"Ad")*CGS_U_pc;//pc
    t1_Bd = FetchDouble(subptr,"Bd")*CGS_U_pc;//pc
    t1_n1 = FetchDouble(subptr,"n1");//pccm
    t1_H1 = FetchDouble(subptr,"H1")*CGS_U_pc;//pc
    // thin disk
    subptr = ptr->FirstChildElement("YMW16")->FirstChildElement("ThinDisk");
    t2_A2 = FetchDouble(subptr,"A2")*CGS_U_pc;//pc
    t2_B2 = FetchDouble(subptr,"B2")*CGS_U_pc;//pc
    t2_n2 = FetchDouble(subptr,"n2");//pccm
    t2_K2 = FetchDouble(subptr,"K2");
    // spiral arm
    subptr = ptr->FirstChildElement("YMW16")->FirstChildElement("SpiralArm");
    t3_B2s = FetchDouble(subptr,"B2s")*CGS_U_pc; //pc
    t3_narm[0] = FetchDouble(subptr,"Ele_arm_0");//pccm
    t3_narm[1] = FetchDouble(subptr,"Ele_arm_1");
    t3_narm[2] = FetchDouble(subptr,"Ele_arm_2");
    t3_narm[3] = FetchDouble(subptr,"Ele_arm_3");
    t3_narm[4] = FetchDouble(subptr,"Ele_arm_4");
    t3_warm[0] = FetchDouble(subptr,"Wid_arm_0")*CGS_U_pc;//pc
    t3_warm[1] = FetchDouble(subptr,"Wid_arm_1")*CGS_U_pc;
    t3_warm[2] = FetchDouble(subptr,"Wid_arm_2")*CGS_U_pc;
    t3_warm[3] = FetchDouble(subptr,"Wid_arm_3")*CGS_U_pc;
    t3_warm[4] = FetchDouble(subptr,"Wid_arm_4")*CGS_U_pc;
    t3_Aa = FetchDouble(subptr,"Aa")*CGS_U_pc;//pc
    t3_Ka = FetchDouble(subptr,"Ka");
    t3_ncn = FetchDouble(subptr,"ncn");
    t3_thetacn = FetchDouble(subptr,"thetacn");//deg
    t3_wcn = FetchDouble(subptr,"wcn");//deg
    t3_nsg = FetchDouble(subptr,"nsg");
    t3_thetasg = FetchDouble(subptr,"thetasg");//deg
    t3_wsg = FetchDouble(subptr,"wsg");//deg
    // gc
    subptr = ptr->FirstChildElement("YMW16")->FirstChildElement("GalCenter");
    t4_ngc = FetchDouble(subptr,"ngc");//pccm
    t4_Agc = FetchDouble(subptr,"Agc")*CGS_U_pc;//pc
    t4_Hgc = FetchDouble(subptr,"Hgc")*CGS_U_pc;//pc
    // Gum
    subptr = ptr->FirstChildElement("YMW16")->FirstChildElement("GumNebula");
    t5_ngn = FetchDouble(subptr,"ngn");//pccm
    t5_Wgn = FetchDouble(subptr,"Wgn")*CGS_U_pc;//pc
    t5_Agn = FetchDouble(subptr,"Agn")*CGS_U_pc;//pc
    t5_Kgn = FetchDouble(subptr,"Kgn");
    // Local Bubble
    subptr = ptr->FirstChildElement("YMW16")->FirstChildElement("LocalBubble");
    t6_J_LB = FetchDouble(subptr,"J_LB");
    t6_nlb1 = FetchDouble(subptr,"nlb1");//pccm
    t6_thetalb1 = FetchDouble(subptr,"thetalb1");//deg
    t6_detlb1 = FetchDouble(subptr,"detlb1");//deg
    t6_wlb1 = FetchDouble(subptr,"wlb1")*CGS_U_pc;//pc
    t6_hlb1 = FetchDouble(subptr,"hlb1")*CGS_U_pc;//pc
    t6_nlb2 = FetchDouble(subptr,"nlb2");//pccm
    t6_thetalb2 = FetchDouble(subptr,"thetalb2");//deg
    t6_detlb2 = FetchDouble(subptr,"detlb2");//deg
    t6_wlb2 = FetchDouble(subptr,"wlb2")*CGS_U_pc;//pc
    t6_hlb2 = FetchDouble(subptr,"hlb2")*CGS_U_pc;//pc
    // Loop I
    subptr = ptr->FirstChildElement("YMW16")->FirstChildElement("LoopI");
    t7_nLI = FetchDouble(subptr,"nLI");//pccm
    t7_RLI = FetchDouble(subptr,"RLI")*CGS_U_pc;//pc
    t7_WLI = FetchDouble(subptr,"WLI")*CGS_U_pc;//pc
    t7_detthetaLI = FetchDouble(subptr,"detthetaLI");//deg
    t7_thetaLI = FetchDouble(subptr,"thetaLI");//deg
    
    // gaussian random
    subptr = ptr->FirstChildElement("Random")->FirstChildElement("Iso");
    // nomalization
    fernd_iso.push_back(FetchDouble(subptr,"rms"));
    // turning point
    fernd_iso.push_back(FetchDouble(subptr,"k0"));
    // index before turning (large scales)
    fernd_iso.push_back(FetchDouble(subptr,"a0"));
    
    // rescaling parameters for random fe
    subptr = ptr->FirstChildElement("Random")->FirstChildElement("Rescal");
    fernd_scal.push_back(FetchDouble(subptr,"r0"));
    fernd_scal.push_back(FetchDouble(subptr,"z0"));
    
}

void Pond::cre_pond(XMLDocument *doc){
    XMLElement *ptr = doc->FirstChildElement("root")->FirstChildElement("CRE");
    sim_freq = ptr->DoubleAttribute("freq")*CGS_U_GHz;
    // analytical
    XMLElement *subptr = ptr->FirstChildElement("Analytic");
    // the default parameter setting follows wmap3yr model
    creana.push_back(FetchDouble(subptr,"alpha"));
    creana.push_back(FetchDouble(subptr,"beta"));
    creana.push_back(FetchDouble(subptr,"theta"));
    creana.push_back(FetchDouble(subptr,"hr")); //kpc
    creana.push_back(FetchDouble(subptr,"hz")); //kpc
    // differential flux normalization factor at gamma 10
    creana.push_back(FetchDouble(subptr,"je"));
}

// auxiliary functions
std::string Pond::FetchString(XMLElement* el, string obj){
    XMLElement* el1 = el->FirstChildElement(obj.c_str());
    return el1->Attribute("value");
}

int Pond::FetchInt(XMLElement* el, string obj){
    XMLElement* el1 = el->FirstChildElement(obj.c_str());
    return el1->IntAttribute("value");
}

unsigned int Pond::FetchUnsigned(XMLElement* el, string obj){
    XMLElement* el1 = el->FirstChildElement(obj.c_str());
    return el1->UnsignedAttribute("value");
}

bool Pond::FetchBool(XMLElement* el, string obj){
    XMLElement* el1 = el->FirstChildElement(obj.c_str());
    return el1->BoolAttribute("value");
}

double Pond::FetchDouble(XMLElement* el, string obj){
    XMLElement* el1 = el->FirstChildElement(obj.c_str());
    return el1->DoubleAttribute("value");
}

// END
