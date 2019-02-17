#include <string>
#include <cmath>
#include <memory>

#include <hvec.h>
#include <tinyxml2.h>

#include <param.h>
#include <cgs_units_file.h>
#include <namespace_toolkit.h>

Param::Param (const std::string file_name){
    std::unique_ptr<tinyxml2::XMLDocument> doc = toolkit::loadxml(file_name);
    // gc sun position
    tinyxml2::XMLElement *ptr {toolkit::tracexml(doc.get(),{"Grid","SunPosition"})};
    SunPosition = hvec<3,double> {CGS_U_kpc*toolkit::fetchdouble(ptr,"value","x"),
        CGS_U_kpc*toolkit::fetchdouble(ptr,"value","y"),
        CGS_U_pc*toolkit::fetchdouble(ptr,"value","z")};
    //
    grid_param (doc.get());
    b_param (doc.get());
    fe_param (doc.get());
    cre_param (doc.get());
}

// grid
void Param::grid_param (tinyxml2::XMLDocument *doc){
    // breg box
    tinyxml2::XMLElement *ptr {toolkit::tracexml(doc,{"Grid","Box_GMF"})};
    grid_breg.nx = toolkit::fetchunsigned (ptr,"value","nx");
    grid_breg.ny = toolkit::fetchunsigned (ptr,"value","ny");
    grid_breg.nz = toolkit::fetchunsigned (ptr,"value","nz");
    grid_breg.full_size = grid_breg.nx*grid_breg.ny*grid_breg.nz;
    grid_breg.x_max = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","x_max");
    grid_breg.x_min = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","x_min");
    grid_breg.y_max = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","y_max");
    grid_breg.y_min = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","y_min");
    grid_breg.z_max = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","z_max");
    grid_breg.z_min = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","z_min");
    // brnd box
    ptr = toolkit::tracexml(doc,{"Grid","Box_GMF"});
    grid_brnd.nx = toolkit::fetchunsigned (ptr,"value","nx");
    grid_brnd.ny = toolkit::fetchunsigned (ptr,"value","ny");
    grid_brnd.nz = toolkit::fetchunsigned (ptr,"value","nz");
    grid_brnd.full_size = grid_brnd.nx*grid_brnd.ny*grid_brnd.nz;
    grid_brnd.x_max = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","x_max");
    grid_brnd.x_min = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","x_min");
    grid_brnd.y_max = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","y_max");
    grid_brnd.y_min = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","y_min");
    grid_brnd.z_max = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","z_max");
    grid_brnd.z_min = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","z_min");
    // gmf permission
    ptr = toolkit::tracexml (doc,{"MagneticField"});
    grid_breg.build_permission = toolkit::fetchbool (ptr,"cue","Regular");
    grid_brnd.build_permission = toolkit::fetchbool (ptr,"cue","Random");
    // fereg box
    ptr = toolkit::tracexml (doc,{"Grid","Box_FE"});
    grid_fereg.nx = toolkit::fetchunsigned (ptr,"value","nx");
    grid_fereg.ny = toolkit::fetchunsigned (ptr,"value","ny");
    grid_fereg.nz = toolkit::fetchunsigned (ptr,"value","nz");
    grid_fereg.full_size = grid_fereg.nx*grid_fereg.ny*grid_fereg.nz;
    grid_fereg.x_max = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","x_max");
    grid_fereg.x_min = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","x_min");
    grid_fereg.y_max = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","y_max");
    grid_fereg.y_min = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","y_min");
    grid_fereg.z_max = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","z_max");
    grid_fereg.z_min = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","z_min");
    // fernd box
    ptr = toolkit::tracexml (doc,{"Grid","Box_FE"});
    grid_fernd.nx = toolkit::fetchunsigned (ptr,"value","nx");
    grid_fernd.ny = toolkit::fetchunsigned (ptr,"value","ny");
    grid_fernd.nz = toolkit::fetchunsigned (ptr,"value","nz");
    grid_fernd.full_size = grid_fernd.nx*grid_fernd.ny*grid_fernd.nz;
    grid_fernd.x_max = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","x_max");
    grid_fernd.x_min = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","x_min");
    grid_fernd.y_max = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","y_max");
    grid_fernd.y_min = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","y_min");
    grid_fernd.z_max = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","z_max");
    grid_fernd.z_min = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","z_min");
    // fe permission
    ptr = toolkit::tracexml (doc,{"FreeElectron"});
    grid_fereg.build_permission = toolkit::fetchbool (ptr,"cue","Regular");
    grid_fernd.build_permission = toolkit::fetchbool (ptr,"cue","Random");
    // cre box
    ptr = toolkit::tracexml (doc,{"Grid","Box_CRE"});
    grid_cre.E_min = CGS_U_GeV*toolkit::fetchdouble (ptr,"value","E_min");
    grid_cre.E_max = CGS_U_GeV*toolkit::fetchdouble (ptr,"value","E_max");
    grid_cre.nE = toolkit::fetchunsigned (ptr,"value","nE");
    grid_cre.E_fact = std::log(grid_cre.E_max/grid_cre.E_min)/(grid_cre.nE-1);
    grid_cre.nz = toolkit::fetchunsigned (ptr,"value","nz");
    grid_cre.nx = toolkit::fetchunsigned (ptr,"value","nx");
    grid_cre.ny = toolkit::fetchunsigned (ptr,"value","ny");
    grid_cre.x_max = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","x_max");
    grid_cre.x_min = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","x_min");
    grid_cre.y_max = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","y_max");
    grid_cre.y_min = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","y_min");
    grid_cre.z_max = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","z_max");
    grid_cre.z_min = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","z_min");
    grid_cre.cre_size = grid_cre.nE*grid_cre.nx*grid_cre.ny*grid_cre.nz;
    // int box
    ptr = toolkit::tracexml (doc,{"Grid","Shell"});
	grid_int.nside_sim = toolkit::fetchunsigned (ptr,"value","nside_sim");
    grid_int.npix_sim = 12*grid_int.nside_sim*grid_int.nside_sim;
    grid_int.ec_r_max = toolkit::fetchdouble (ptr,"value","ec_r_max")*CGS_U_kpc;
    grid_int.gc_r_max = toolkit::fetchdouble (ptr,"value","gc_r_max")*CGS_U_kpc;
    grid_int.gc_z_max = toolkit::fetchdouble (ptr,"value","gc_z_max")*CGS_U_kpc;
    grid_int.radial_res = toolkit::fetchdouble (ptr,"value","ec_r_res")*CGS_U_kpc;
    grid_int.lat_lim = toolkit::fetchdouble (ptr,"value","lat_lim")*CGS_U_rad;
	if (toolkit::fetchstring (ptr,"type","layer")=="auto"){
        auto subptr = toolkit::tracexml (doc,{"Grid","Shell","layer","auto"});
        grid_int.total_shell = toolkit::fetchunsigned (subptr,"value","shell_num");
        for (std::size_t i=0;i!=grid_int.total_shell;++i){
            grid_int.nside_shell.push_back (pow(2,i)*toolkit::fetchunsigned (subptr,"value","nside_min"));
        }
        grid_int.radii_shell.push_back (0.);
        for (std::size_t i=0;i<grid_int.total_shell;++i) {
            grid_int.radii_shell.push_back (grid_int.ec_r_max*std::pow(0.5,grid_int.total_shell-i-1));
        }
    }
    else if (toolkit::fetchstring (ptr,"type","layer")=="manual"){
        auto subptr = toolkit::tracexml (doc,{"Grid","Shell","layer","manual"});
        grid_int.total_shell = 0;
        for (auto e = subptr->FirstChildElement("nside");e!=nullptr;e=e->NextSiblingElement("nside")){
            grid_int.total_shell++;
            grid_int.nside_shell.push_back (toolkit::fetchunsigned(e,"value"));
        }
        for (auto e = subptr->FirstChildElement("cut");e!=nullptr;e=e->NextSiblingElement("cut")){
            grid_int.cut_shell.push_back (toolkit::fetchdouble(e,"value"));
        }
		grid_int.cut_shell.push_back (1.);
        assert (grid_int.cut_shell.size() == grid_int.total_shell);
        assert (grid_int.nside_shell.size() == grid_int.total_shell);
		grid_int.radii_shell.push_back (0.);
        for (auto& i : grid_int.cut_shell){
            grid_int.radii_shell.push_back (grid_int.ec_r_max*i);
        }
    }
    else {
        assert (false);
    }
    // output file name
    ptr = toolkit::tracexml (doc,{"Obsout"});
    if (ptr->FirstChildElement("DM")!=nullptr){
        grid_int.do_dm = toolkit::fetchbool (ptr,"cue","DM");
        grid_int.sim_dm_name = toolkit::fetchstring (ptr,"filename","DM");
    }
    else {
        grid_int.do_dm = false;
    }
    if (ptr->FirstChildElement("Faraday")!=nullptr){
        grid_int.do_fd = toolkit::fetchbool (ptr,"cue","Faraday");
        grid_int.sim_fd_name = toolkit::fetchstring (ptr,"filename","Faraday");
    }
    else {
        grid_int.do_fd = false;
    }
    if (ptr->FirstChildElement("Sync")!=nullptr){
        auto subptr = toolkit::tracexml (doc,{"Obsout","Sync"});
        grid_int.do_sync.push_back (toolkit::fetchbool (subptr,"cue"));
		grid_int.sim_sync_freq.push_back (toolkit::fetchdouble (subptr,"freq")*CGS_U_GHz);
        grid_int.sim_sync_name.push_back (toolkit::fetchstring (subptr,"filename"));
        for (auto e = subptr->NextSiblingElement("Sync");e!=nullptr;e=e->NextSiblingElement("Sync")){
            grid_int.do_sync.push_back (toolkit::fetchbool (e,"cue"));
            grid_int.sim_sync_freq.push_back (toolkit::fetchdouble (e,"freq")*CGS_U_GHz);
            grid_int.sim_sync_name.push_back (toolkit::fetchstring (e,"filename"));
        }
    }
    else {
        grid_int.do_sync.push_back (false);
    }
    // output field params
    ptr = toolkit::tracexml (doc,{"Fieldout"});
    grid_breg.read_permission = toolkit::fetchbool (ptr,"read","breg_grid");
    grid_breg.write_permission = toolkit::fetchbool (ptr,"write","breg_grid");
    grid_breg.filename = toolkit::fetchstring(ptr,"filename","breg_grid");
    grid_brnd.read_permission = toolkit::fetchbool (ptr,"read","brnd_grid");
    grid_brnd.write_permission = toolkit::fetchbool (ptr,"write","brnd_grid");
    grid_brnd.filename = toolkit::fetchstring (ptr,"filename","brnd_grid");
    grid_fereg.read_permission = toolkit::fetchbool (ptr,"read","fereg_grid");
    grid_fereg.write_permission = toolkit::fetchbool (ptr,"write","fereg_grid");
    grid_fereg.filename = toolkit::fetchstring (ptr,"filename","fereg_grid");
    grid_fernd.read_permission = toolkit::fetchbool (ptr,"read","fernd_grid");
    grid_fernd.write_permission = toolkit::fetchbool (ptr,"write","fernd_grid");
    grid_fernd.filename = toolkit::fetchstring (ptr,"filename","fernd_grid");
    grid_cre.read_permission = toolkit::fetchbool (ptr,"read","cre_grid");
    grid_cre.write_permission = toolkit::fetchbool (ptr,"write","cre_grid");
    grid_cre.filename = toolkit::fetchstring (ptr,"filename","cre_grid");
}

// magnetic field
void Param::b_param (tinyxml2::XMLDocument *doc){
    tinyxml2::XMLElement *ptr {toolkit::tracexml(doc,{"MagneticField"})};
    breg_type = toolkit::fetchstring(ptr,"type","Regular");
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
        brnd_type = toolkit::fetchstring(ptr,"type","Random");
        // brnd_global
        if(brnd_type=="Global"){
            tinyxml2::XMLElement *subptr {toolkit::tracexml(doc,{"MagneticField","Random","Global"})};
            brnd_method = toolkit::fetchstring(subptr,"type");
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
            brnd_method = toolkit::fetchstring(subptr,"type");
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
    fereg_type = toolkit::fetchstring(ptr,"type","Regular");
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
        fernd_type = toolkit::fetchstring(ptr,"type","Random");
        // global turbulent
        if(fernd_type=="Global"){
            tinyxml2::XMLElement *subptr {toolkit::tracexml(doc,{"FreeElectron","Random","Global"})};
            fernd_method = toolkit::fetchstring(subptr,"type");
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
    cre_type = ptr->Attribute("type");
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
