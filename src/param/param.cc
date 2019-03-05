#include <string>
#include <cmath>
#include <memory>

#include <hvec.h>
#include <tinyxml2.h>

#include <param.h>
#include <cgs_units_file.h>
#include <namespace_toolkit.h>

Param::Param (const std::string file_name){
    // load xml file
    std::unique_ptr<tinyxml2::XMLDocument> doc {toolkit::loadxml(file_name)};
    // observer position
    tinyxml2::XMLElement* ptr {toolkit::tracexml(doc.get(),{"grid","observer"})};
    observer = hvec<3,double> {CGS_U_kpc*toolkit::fetchdouble(ptr,"value","x",-8.3),
        CGS_U_kpc*toolkit::fetchdouble(ptr,"value","y",0),
        CGS_U_pc*toolkit::fetchdouble(ptr,"value","z",6)};
    // collect parameters
    obs_param (doc.get());
    breg_param (doc.get());
    brnd_param (doc.get());
    fereg_param (doc.get());
    fernd_param (doc.get());
    cre_param (doc.get());
}

void Param::obs_param (tinyxml2::XMLDocument* doc){
    // observable base path
    tinyxml2::XMLElement* ptr {toolkit::tracexml (doc,{"observable"})};
    // controller for reading shell parameters
    grid_int.write_permission = false;
    // if dispersion measure is required
    if (ptr->FirstChildElement("dm")!=nullptr){
        grid_int.write_permission = true;
        grid_int.do_dm = toolkit::fetchbool (ptr,"cue","dm",0);
        grid_int.sim_dm_name = toolkit::fetchstring (ptr,"filename","dm");
        grid_int.nside_dm = toolkit::fetchunsigned (ptr,"nside","dm");
        grid_int.npix_dm = 12*grid_int.nside_dm*grid_int.nside_dm;
    }
    else {
        grid_int.do_dm = false;
    }
    // if faraday depth is required
    if (ptr->FirstChildElement("faraday")!=nullptr){
        grid_int.write_permission = true;
        grid_int.do_fd = toolkit::fetchbool (ptr,"cue","faraday",0);
        grid_int.sim_fd_name = toolkit::fetchstring (ptr,"filename","faraday");
        grid_int.nside_fd = toolkit::fetchunsigned (ptr,"nside","faraday");
        grid_int.npix_fd = 12*grid_int.nside_fd*grid_int.nside_fd;
    }
    else {
        grid_int.do_fd = false;
    }
    // if synchrotron emission is required
    if (ptr->FirstChildElement("sync")!=nullptr){
        grid_int.write_permission = true;
        tinyxml2::XMLElement* subptr {toolkit::tracexml (doc,{"observable","sync"})};
        grid_int.do_sync.push_back (toolkit::fetchbool (subptr,"cue",0));
        grid_int.sim_sync_freq.push_back (toolkit::fetchdouble (subptr,"freq")*CGS_U_GHz);
        grid_int.sim_sync_name.push_back (toolkit::fetchstring (subptr,"filename"));
        grid_int.nside_sync.push_back (toolkit::fetchunsigned (subptr,"nside"));
        grid_int.npix_sync.push_back (12*grid_int.nside_sync.back()*grid_int.nside_sync.back());
        for (auto e = subptr->NextSiblingElement("sync");e!=nullptr;e=e->NextSiblingElement("sync")){
            grid_int.do_sync.push_back (toolkit::fetchbool (e,"cue",0));
            grid_int.sim_sync_freq.push_back (toolkit::fetchdouble (e,"freq")*CGS_U_GHz);
            grid_int.sim_sync_name.push_back (toolkit::fetchstring (e,"filename"));
            grid_int.nside_sync.push_back (toolkit::fetchunsigned (e,"nside"));
            grid_int.npix_sync.push_back (12*grid_int.nside_sync.back()*grid_int.nside_sync.back());
        }
    }
    else {
        grid_int.do_sync.push_back (false);
    }
    // if any observable is requried
    if (grid_int.write_permission){
        ptr = toolkit::tracexml (doc,{"grid","shell"});
        grid_int.ec_r_min = toolkit::fetchdouble (ptr,"value","ec_r_min",0)*CGS_U_kpc;
        grid_int.ec_r_max = toolkit::fetchdouble (ptr,"value","ec_r_max",10)*CGS_U_kpc;
        grid_int.gc_r_min = toolkit::fetchdouble (ptr,"value","gc_r_min",0)*CGS_U_kpc;
        grid_int.gc_r_max = toolkit::fetchdouble (ptr,"value","gc_r_max",20)*CGS_U_kpc;
        grid_int.gc_z_min = toolkit::fetchdouble (ptr,"value","gc_z_min",0)*CGS_U_kpc;
        grid_int.gc_z_max = toolkit::fetchdouble (ptr,"value","gc_z_max",10)*CGS_U_kpc;
        grid_int.radial_res = toolkit::fetchdouble (ptr,"value","ec_r_res",0.01)*CGS_U_kpc;
        grid_int.lat_min = toolkit::fetchdouble (ptr,"value","lat_min",0)*CGS_U_rad;
        grid_int.lat_max = toolkit::fetchdouble (ptr,"value","lat_max",90)*CGS_U_rad;
        grid_int.lon_min = toolkit::fetchdouble (ptr,"value","lon_min",0)*CGS_U_rad;
        grid_int.lon_max = toolkit::fetchdouble (ptr,"value","lon_max",360)*CGS_U_rad;
        // auto shell cutting
        if (toolkit::fetchstring (ptr,"type","layer")=="auto"){
            tinyxml2::XMLElement* subptr {toolkit::tracexml (doc,{"grid","shell","layer","auto"})};
            grid_int.total_shell = toolkit::fetchunsigned (subptr,"value","shell_num",1);
            for (std::size_t i=0;i!=grid_int.total_shell;++i){
                grid_int.nside_shell.push_back (pow(2,i)*toolkit::fetchunsigned (subptr,"value","nside_min",32));
            }
            std::cout<<"auto shell nside size "<<grid_int.nside_shell.size()<<std::endl;
            grid_int.radii_shell.push_back (grid_int.ec_r_min);
            for (std::size_t i=0;i<grid_int.total_shell;++i) {
                grid_int.radii_shell.push_back ((grid_int.ec_r_max-grid_int.ec_r_min)*std::pow(0.5,grid_int.total_shell-i-1) + grid_int.ec_r_min);
            }
        }
        // manual shell cutting
        else if (toolkit::fetchstring (ptr,"type","layer")=="manual"){
            tinyxml2::XMLElement* subptr {toolkit::tracexml (doc,{"grid","shell","layer","manual"})};
            grid_int.total_shell = 0;
            for (auto e = subptr->FirstChildElement("nside");e!=nullptr;e=e->NextSiblingElement("nside")){
                grid_int.total_shell++;
                grid_int.nside_shell.push_back (toolkit::fetchunsigned(e,"value"));
            }
            std::cout<<"manual shell nside size "<<grid_int.nside_shell.size()<<std::endl;
            for (auto e = subptr->FirstChildElement("cut");e!=nullptr;e=e->NextSiblingElement("cut")){
                grid_int.cut_shell.push_back (toolkit::fetchdouble(e,"value"));
            }
            grid_int.cut_shell.push_back (1.);
            assert (grid_int.cut_shell.size() == grid_int.total_shell);
            assert (grid_int.nside_shell.size() == grid_int.total_shell);
            grid_int.radii_shell.push_back (grid_int.ec_r_min);
            for (auto& i : grid_int.cut_shell){
                grid_int.radii_shell.push_back ((grid_int.ec_r_max-grid_int.ec_r_min)*i+grid_int.ec_r_min);
            }
        }
        else {
            throw std::runtime_error("unsupported layer option");
        }
    }
}

void Param::breg_param (tinyxml2::XMLDocument* doc){
    // breg io
    tinyxml2::XMLElement* ptr {toolkit::tracexml (doc,{"fieldio"})};
    if (ptr->FirstChildElement("breg")!=nullptr){
        grid_breg.read_permission = toolkit::fetchbool (ptr,"read","breg");
        grid_breg.write_permission = toolkit::fetchbool (ptr,"write","breg");
        grid_breg.filename = toolkit::fetchstring(ptr,"filename","breg");
    }
    // breg internal
    ptr = toolkit::tracexml(doc,{"magneticfield"});
    grid_breg.build_permission = toolkit::fetchbool (ptr,"cue","regular",0);
    // if no external read and internal model is active
    if (grid_breg.build_permission and not grid_breg.read_permission){
        breg_type = toolkit::fetchstring(ptr,"type","regular");
        // bwmap
        if (breg_type=="wmap"){
            tinyxml2::XMLElement* subptr {toolkit::tracexml(doc,{"magneticfield","regular","wmap"})};
            breg_wmap.b0 = toolkit::fetchdouble(subptr,"value","b0")*CGS_U_muGauss; //microGauss
            breg_wmap.psi0 = toolkit::fetchdouble(subptr,"value","psi0")*CGS_U_rad; //rad
            breg_wmap.psi1 = toolkit::fetchdouble(subptr,"value","psi1")*CGS_U_rad; //rad
            breg_wmap.chi0 = toolkit::fetchdouble(subptr,"value","chi0")*CGS_U_rad; //rad
        }
        // bjaffe
        else if (breg_type=="jaffe"){
            tinyxml2::XMLElement* subptr {toolkit::tracexml(doc,{"magneticfield","regular","jaffe"})};
            breg_jaffe.quadruple = toolkit::fetchbool(subptr,"cue","quadruple",0);
            breg_jaffe.bss = toolkit::fetchbool(subptr,"cue","bss",0);
            breg_jaffe.disk_amp = toolkit::fetchdouble(subptr,"value","disk_amp")*CGS_U_muGauss; //microG
            breg_jaffe.disk_z0 = toolkit::fetchdouble(subptr,"value","disk_z0")*CGS_U_kpc; //kpc
            breg_jaffe.halo_amp = toolkit::fetchdouble(subptr,"value","halo_amp")*CGS_U_muGauss; //microG
            breg_jaffe.halo_z0 = toolkit::fetchdouble(subptr,"value","halo_z0")*CGS_U_kpc; //kpc
            breg_jaffe.r_inner = toolkit::fetchdouble(subptr,"value","r_inner")*CGS_U_kpc; //kpc
            breg_jaffe.r_scale = toolkit::fetchdouble(subptr,"value","r_scale")*CGS_U_kpc; //kpc
            breg_jaffe.r_peak = toolkit::fetchdouble(subptr,"value","r_peak")*CGS_U_kpc; //kpc
            breg_jaffe.ring = toolkit::fetchbool(subptr,"cue","ring",0);
            breg_jaffe.ring_amp = toolkit::fetchdouble(subptr,"value","ring_amp")*CGS_U_muGauss; //microG
            breg_jaffe.ring_r = toolkit::fetchdouble(subptr,"value","ring_r")*CGS_U_kpc; //kpc
            breg_jaffe.bar = toolkit::fetchbool(subptr,"cue","bar",0);
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
        else if (breg_type=="unif"){
            tinyxml2::XMLElement* subptr {toolkit::tracexml(doc,{"magneticfield","regular","unif"})};
            breg_unif.b0 = toolkit::fetchdouble(subptr,"value","b0")*CGS_U_muGauss; //microGauss
            breg_unif.l0 = toolkit::fetchdouble(subptr,"value","l0")*CGS_U_rad; //rad
            breg_unif.r = toolkit::fetchdouble(subptr,"value","r");
        }
        else{
            throw std::runtime_error("unsupported breg model");
        }
    }
    // breg io box
    if (grid_breg.read_permission or grid_breg.write_permission){
        // breg box
        tinyxml2::XMLElement* subptr {toolkit::tracexml(doc,{"grid","box_breg"})};
        grid_breg.nx = toolkit::fetchunsigned (subptr,"value","nx",800);
        grid_breg.ny = toolkit::fetchunsigned (subptr,"value","ny",800);
        grid_breg.nz = toolkit::fetchunsigned (subptr,"value","nz",160);
        grid_breg.full_size = grid_breg.nx*grid_breg.ny*grid_breg.nz;
        grid_breg.x_max = CGS_U_kpc*toolkit::fetchdouble (subptr,"value","x_max",20);
        grid_breg.x_min = CGS_U_kpc*toolkit::fetchdouble (subptr,"value","x_min",-20);
        grid_breg.y_max = CGS_U_kpc*toolkit::fetchdouble (subptr,"value","y_max",20);
        grid_breg.y_min = CGS_U_kpc*toolkit::fetchdouble (subptr,"value","y_min",-20);
        grid_breg.z_max = CGS_U_kpc*toolkit::fetchdouble (subptr,"value","z_max",4);
        grid_breg.z_min = CGS_U_kpc*toolkit::fetchdouble (subptr,"value","z_min",-4);
    }
}

void Param::brnd_param (tinyxml2::XMLDocument* doc){
    // brnd io
    tinyxml2::XMLElement* ptr {toolkit::tracexml (doc,{"fieldio"})};
    if (ptr->FirstChildElement("brnd")!=nullptr){
        grid_brnd.read_permission = toolkit::fetchbool (ptr,"read","brnd");
        grid_brnd.write_permission = toolkit::fetchbool (ptr,"write","brnd");
        grid_brnd.filename = toolkit::fetchstring (ptr,"filename","brnd");
    }
    // brnd internal
    ptr = toolkit::tracexml(doc,{"magneticfield"});
    grid_brnd.build_permission = toolkit::fetchbool (ptr,"cue","random",0);
    // if no external read and internal model is active
    if (grid_brnd.build_permission and not grid_brnd.read_permission){
        // random seed
        brnd_seed = toolkit::fetchunsigned(ptr,"seed","random");
        brnd_type = toolkit::fetchstring(ptr,"type","random");
        // brnd_global
        if (brnd_type=="global"){
            tinyxml2::XMLElement* subptr {toolkit::tracexml(doc,{"magneticfield","random","global"})};
            brnd_method = toolkit::fetchstring(subptr,"type");
            if (brnd_method=="es"){
                subptr = toolkit::tracexml(doc,{"magneticfield","random","global","es"});
                brnd_es.rms = toolkit::fetchdouble(subptr,"value","rms")*CGS_U_muGauss;
                brnd_es.k0 = toolkit::fetchdouble(subptr,"value","k0");
                brnd_es.a0 = toolkit::fetchdouble(subptr,"value","a0");
                brnd_es.rho = toolkit::fetchdouble(subptr,"value","rho");
                brnd_es.r0 = toolkit::fetchdouble(subptr,"value","r0")*CGS_U_kpc;
                brnd_es.z0 = toolkit::fetchdouble(subptr,"value","z0")*CGS_U_kpc;
            }
            else if (brnd_method=="jaffe"){
                subptr = toolkit::tracexml(doc,{"magneticfield","random","global","jaffe"});
                // to be implemented
            }
            else{
                throw std::runtime_error("unsupported brnd model");
            }
        }
        // brnd_local
        else if (brnd_type=="local"){
            tinyxml2::XMLElement* subptr {toolkit::tracexml(doc,{"magneticfield","random","local"})};
            brnd_method = toolkit::fetchstring(subptr,"type");
            if (brnd_method=="mhd"){
                subptr = toolkit::tracexml(doc,{"magneticfield","random","local","mhd"});
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
            else{
                throw std::runtime_error("unsupported brnd model");
            }
        }
        else{
            throw std::runtime_error("unsupported brnd type");
        }
    }
    // brnd io box
    if (grid_brnd.read_permission or grid_brnd.write_permission or grid_brnd.build_permission){
        // brnd box
        ptr = toolkit::tracexml(doc,{"grid","box_brnd"});
        grid_brnd.nx = toolkit::fetchunsigned (ptr,"value","nx",800);
        grid_brnd.ny = toolkit::fetchunsigned (ptr,"value","ny",800);
        grid_brnd.nz = toolkit::fetchunsigned (ptr,"value","nz",160);
        grid_brnd.full_size = grid_brnd.nx*grid_brnd.ny*grid_brnd.nz;
        grid_brnd.x_max = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","x_max",20);
        grid_brnd.x_min = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","x_min",-20);
        grid_brnd.y_max = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","y_max",20);
        grid_brnd.y_min = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","y_min",-20);
        grid_brnd.z_max = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","z_max",4);
        grid_brnd.z_min = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","z_min",-4);
    }
}

void Param::fereg_param (tinyxml2::XMLDocument* doc){
    // fereg io
    tinyxml2::XMLElement* ptr {toolkit::tracexml (doc,{"fieldio"})};
    if (ptr->FirstChildElement("fereg")!=nullptr){
        grid_fereg.read_permission = toolkit::fetchbool (ptr,"read","fereg");
        grid_fereg.write_permission = toolkit::fetchbool (ptr,"write","fereg");
        grid_fereg.filename = toolkit::fetchstring (ptr,"filename","fereg");
    }
    // fereg internal
    ptr = toolkit::tracexml(doc,{"freeelectron"});
    grid_fereg.build_permission = toolkit::fetchbool (ptr,"cue","regular",0);
    if (grid_fereg.build_permission){
        fereg_type = toolkit::fetchstring(ptr,"type","regular");
        // YMW16
        if (fereg_type=="ymw16"){
            // Warp_Sun
            tinyxml2::XMLElement* subptr {toolkit::tracexml(doc,{"freeelectron","regular","ymw16","warp"})};
            fereg_ymw16.r_warp = toolkit::fetchdouble(subptr,"value","r_warp",8.4)*CGS_U_kpc; //kpc
            fereg_ymw16.r0 = toolkit::fetchdouble(subptr,"value","r0",8.3)*CGS_U_kpc; //kpc
            fereg_ymw16.t0_gamma_w = toolkit::fetchdouble(subptr,"value","gamma_w",0.14);
            // thick disk
            subptr = toolkit::tracexml(doc,{"freeelectron","regular","ymw16","thickdisk"});
            fereg_ymw16.t1_ad = toolkit::fetchdouble(subptr,"value","ad",2500)*CGS_U_pc;//pc
            fereg_ymw16.t1_bd = toolkit::fetchdouble(subptr,"value","bd",15000)*CGS_U_pc;//pc
            fereg_ymw16.t1_n1 = toolkit::fetchdouble(subptr,"value","n1",0.01132);//pccm
            fereg_ymw16.t1_h1 = toolkit::fetchdouble(subptr,"value","h1",1673)*CGS_U_pc;//pc
            // thin disk
            subptr = toolkit::tracexml(doc,{"freeelectron","regular","ymw16","thindisk"});
            fereg_ymw16.t2_a2 = toolkit::fetchdouble(subptr,"value","a2",1200)*CGS_U_pc;//pc
            fereg_ymw16.t2_b2 = toolkit::fetchdouble(subptr,"value","b2",4000)*CGS_U_pc;//pc
            fereg_ymw16.t2_n2 = toolkit::fetchdouble(subptr,"value","n2",0.404);//pccm
            fereg_ymw16.t2_k2 = toolkit::fetchdouble(subptr,"value","k2",1.54);
            // spiral arm
            subptr = toolkit::tracexml(doc,{"freeelectron","regular","ymw16","spiralarm"});
            fereg_ymw16.t3_b2s = toolkit::fetchdouble(subptr,"value","b2s",4000)*CGS_U_pc; //pc
            fereg_ymw16.t3_narm[0] = toolkit::fetchdouble(subptr,"value","ele_arm_0",0.135000);//pccm
            fereg_ymw16.t3_narm[1] = toolkit::fetchdouble(subptr,"value","ele_arm_1",0.129000);
            fereg_ymw16.t3_narm[2] = toolkit::fetchdouble(subptr,"value","ele_arm_2",0.103000);
            fereg_ymw16.t3_narm[3] = toolkit::fetchdouble(subptr,"value","ele_arm_3",0.116000);
            fereg_ymw16.t3_narm[4] = toolkit::fetchdouble(subptr,"value","ele_arm_4",0.005700);
            fereg_ymw16.t3_warm[0] = toolkit::fetchdouble(subptr,"value","wid_arm_0",300)*CGS_U_pc;//pc
            fereg_ymw16.t3_warm[1] = toolkit::fetchdouble(subptr,"value","wid_arm_1",500)*CGS_U_pc;
            fereg_ymw16.t3_warm[2] = toolkit::fetchdouble(subptr,"value","wid_arm_2",300)*CGS_U_pc;
            fereg_ymw16.t3_warm[3] = toolkit::fetchdouble(subptr,"value","wid_arm_3",500)*CGS_U_pc;
            fereg_ymw16.t3_warm[4] = toolkit::fetchdouble(subptr,"value","wid_arm_4",300)*CGS_U_pc;
            fereg_ymw16.t3_rmin[0] = toolkit::fetchdouble(subptr,"value","rref_arm_0",3.35)*CGS_U_kpc;//kpc
            fereg_ymw16.t3_rmin[1] = toolkit::fetchdouble(subptr,"value","rref_arm_1",3.707)*CGS_U_kpc;
            fereg_ymw16.t3_rmin[2] = toolkit::fetchdouble(subptr,"value","rref_arm_2",3.56)*CGS_U_kpc;
            fereg_ymw16.t3_rmin[3] = toolkit::fetchdouble(subptr,"value","rref_arm_3",3.670)*CGS_U_kpc;
            fereg_ymw16.t3_rmin[4] = toolkit::fetchdouble(subptr,"value","rref_arm_4",8.21)*CGS_U_kpc;
            fereg_ymw16.t3_phimin[0] = toolkit::fetchdouble(subptr,"value","phiref_arm_0",44.4)*CGS_U_rad;//rad
            fereg_ymw16.t3_phimin[1] = toolkit::fetchdouble(subptr,"value","phiref_arm_1",120.0)*CGS_U_rad;//rad
            fereg_ymw16.t3_phimin[2] = toolkit::fetchdouble(subptr,"value","phiref_arm_2",218.6)*CGS_U_rad;//rad
            fereg_ymw16.t3_phimin[3] = toolkit::fetchdouble(subptr,"value","phiref_arm_3",330.3)*CGS_U_rad;//rad
            fereg_ymw16.t3_phimin[4] = toolkit::fetchdouble(subptr,"value","phiref_arm_4",55.1)*CGS_U_rad;//rad
            fereg_ymw16.t3_tpitch[0] = tan(toolkit::fetchdouble(subptr,"value","pitch_arm_0",11.43)*CGS_U_rad);
            fereg_ymw16.t3_tpitch[1] = tan(toolkit::fetchdouble(subptr,"value","pitch_arm_1",9.84)*CGS_U_rad);
            fereg_ymw16.t3_tpitch[2] = tan(toolkit::fetchdouble(subptr,"value","pitch_arm_2",10.38)*CGS_U_rad);
            fereg_ymw16.t3_tpitch[3] = tan(toolkit::fetchdouble(subptr,"value","pitch_arm_3",10.54)*CGS_U_rad);
            fereg_ymw16.t3_tpitch[4] = tan(toolkit::fetchdouble(subptr,"value","pitch_arm_4",2.77)*CGS_U_rad);
            fereg_ymw16.t3_cpitch[0] = cos(toolkit::fetchdouble(subptr,"value","pitch_arm_0",11.43)*CGS_U_rad);
            fereg_ymw16.t3_cpitch[1] = cos(toolkit::fetchdouble(subptr,"value","pitch_arm_1",9.84)*CGS_U_rad);
            fereg_ymw16.t3_cpitch[2] = cos(toolkit::fetchdouble(subptr,"value","pitch_arm_2",10.38)*CGS_U_rad);
            fereg_ymw16.t3_cpitch[3] = cos(toolkit::fetchdouble(subptr,"value","pitch_arm_3",10.54)*CGS_U_rad);
            fereg_ymw16.t3_cpitch[4] = cos(toolkit::fetchdouble(subptr,"value","pitch_arm_4",2.77)*CGS_U_rad);
            fereg_ymw16.t3_aa = toolkit::fetchdouble(subptr,"value","aa",11680)*CGS_U_pc;//pc
            fereg_ymw16.t3_ka = toolkit::fetchdouble(subptr,"value","ka",5.015);
            fereg_ymw16.t3_ncn = toolkit::fetchdouble(subptr,"value","ncn",2.4);
            fereg_ymw16.t3_thetacn = toolkit::fetchdouble(subptr,"value","thetacn",109);//deg
            fereg_ymw16.t3_wcn = toolkit::fetchdouble(subptr,"value","wcn",8.2);//deg
            fereg_ymw16.t3_nsg = toolkit::fetchdouble(subptr,"value","nsg",0.626);
            fereg_ymw16.t3_thetasg = toolkit::fetchdouble(subptr,"value","thetasg",75.8);//deg
            fereg_ymw16.t3_wsg = toolkit::fetchdouble(subptr,"value","wsg",20);//deg
            // gc
            subptr = toolkit::tracexml(doc,{"freeelectron","regular","ymw16","galcenter"});
            fereg_ymw16.t4_ngc = toolkit::fetchdouble(subptr,"value","ngc",6.2);//pccm
            fereg_ymw16.t4_agc = toolkit::fetchdouble(subptr,"value","agc",160)*CGS_U_pc;//pc
            fereg_ymw16.t4_hgc = toolkit::fetchdouble(subptr,"value","hgc",35)*CGS_U_pc;//pc
            // Gum
            subptr = toolkit::tracexml(doc,{"freeelectron","regular","ymw16","gumnebula"});
            fereg_ymw16.t5_ngn = toolkit::fetchdouble(subptr,"value","ngn",1.84);//pccm
            fereg_ymw16.t5_wgn = toolkit::fetchdouble(subptr,"value","wgn",15.1)*CGS_U_pc;//pc
            fereg_ymw16.t5_agn = toolkit::fetchdouble(subptr,"value","agn",125.8)*CGS_U_pc;//pc
            fereg_ymw16.t5_kgn = toolkit::fetchdouble(subptr,"value","kgn",1.4);
            // Local Bubble
            subptr = toolkit::tracexml(doc,{"freeelectron","regular","ymw16","localbubble"});
            fereg_ymw16.t6_j_lb = toolkit::fetchdouble(subptr,"value","j_lb",0.480);
            fereg_ymw16.t6_nlb1 = toolkit::fetchdouble(subptr,"value","nlb1",1.094);//pccm
            fereg_ymw16.t6_thetalb1 = toolkit::fetchdouble(subptr,"value","thetalb1",195.4);//deg
            fereg_ymw16.t6_detlb1 = toolkit::fetchdouble(subptr,"value","detlb1",28.4);//deg
            fereg_ymw16.t6_wlb1 = toolkit::fetchdouble(subptr,"value","wlb1",14.2)*CGS_U_pc;//pc
            fereg_ymw16.t6_hlb1 = toolkit::fetchdouble(subptr,"value","hlb1",112.9)*CGS_U_pc;//pc
            fereg_ymw16.t6_nlb2 = toolkit::fetchdouble(subptr,"value","nlb2",2.33);//pccm
            fereg_ymw16.t6_thetalb2 = toolkit::fetchdouble(subptr,"value","thetalb2",278.2);//deg
            fereg_ymw16.t6_detlb2 = toolkit::fetchdouble(subptr,"value","detlb2",14.7);//deg
            fereg_ymw16.t6_wlb2 = toolkit::fetchdouble(subptr,"value","wlb2",15.6)*CGS_U_pc;//pc
            fereg_ymw16.t6_hlb2 = toolkit::fetchdouble(subptr,"value","hlb2",43.6)*CGS_U_pc;//pc
            // Loop I
            subptr = toolkit::tracexml(doc,{"freeelectron","regular","ymw16","loopi"});
            fereg_ymw16.t7_nli = toolkit::fetchdouble(subptr,"value","nli",1.907);//pccm
            fereg_ymw16.t7_rli = toolkit::fetchdouble(subptr,"value","rli",80.0)*CGS_U_pc;//pc
            fereg_ymw16.t7_wli = toolkit::fetchdouble(subptr,"value","wli",15.0)*CGS_U_pc;//pc
            fereg_ymw16.t7_detthetali = toolkit::fetchdouble(subptr,"value","detthetali",30.0);//deg
            fereg_ymw16.t7_thetali = toolkit::fetchdouble(subptr,"value","thetali",40.0);//deg
        }
        // uniform
        else if (fereg_type=="unif"){
            tinyxml2::XMLElement* subptr {toolkit::tracexml(doc,{"freeelectron","regular","unif"})};
            fereg_unif.n0 = toolkit::fetchdouble(subptr,"value","n0");
            fereg_unif.r0 = toolkit::fetchdouble(subptr,"value","r0")*CGS_U_kpc; //kpc
        }
        else{
            throw std::runtime_error("unsupported fereg model");
        }
    }
    // breg io box
    if (grid_fereg.read_permission or grid_fereg.write_permission){
        ptr = toolkit::tracexml (doc,{"grid","box_fereg"});
        grid_fereg.nx = toolkit::fetchunsigned (ptr,"value","nx",800);
        grid_fereg.ny = toolkit::fetchunsigned (ptr,"value","ny",800);
        grid_fereg.nz = toolkit::fetchunsigned (ptr,"value","nz",800);
        grid_fereg.full_size = grid_fereg.nx*grid_fereg.ny*grid_fereg.nz;
        grid_fereg.x_max = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","x_max",20);
        grid_fereg.x_min = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","x_min",-20);
        grid_fereg.y_max = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","y_max",20);
        grid_fereg.y_min = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","y_min",-20);
        grid_fereg.z_max = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","z_max",4);
        grid_fereg.z_min = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","z_min",-4);
    }
}

void Param::fernd_param (tinyxml2::XMLDocument* doc){
    // fernd io
    tinyxml2::XMLElement* ptr {toolkit::tracexml (doc,{"fieldio"})};
    if (ptr->FirstChildElement("fernd")!=nullptr){
        grid_fernd.read_permission = toolkit::fetchbool (ptr,"read","fernd");
        grid_fernd.write_permission = toolkit::fetchbool (ptr,"write","fernd");
        grid_fernd.filename = toolkit::fetchstring (ptr,"filename","fernd");
    }
    // fernd internal
    ptr = toolkit::tracexml(doc,{"freeelectron"});
    grid_fernd.build_permission = toolkit::fetchbool (ptr,"cue","random",0);
    if (grid_fernd.build_permission){
        // random seed
        fernd_seed = toolkit::fetchunsigned(ptr,"seed","random");
        fernd_type = toolkit::fetchstring(ptr,"type","random");
        // global turbulent
        if (fernd_type=="global"){
            tinyxml2::XMLElement* subptr {toolkit::tracexml(doc,{"freeelectron","random","global"})};
            fernd_method = toolkit::fetchstring(subptr,"type");
            if (fernd_method=="dft"){
                subptr = toolkit::tracexml(doc,{"freeelectron","random","global","dft"});
                fernd_dft.rms = toolkit::fetchdouble(subptr,"value","rms");
                fernd_dft.k0 = toolkit::fetchdouble(subptr,"value","k0");
                fernd_dft.a0 = toolkit::fetchdouble(subptr,"value","a0");
                fernd_dft.r0 = toolkit::fetchdouble(subptr,"value","r0")*CGS_U_kpc;
                fernd_dft.z0 = toolkit::fetchdouble(subptr,"value","z0")*CGS_U_kpc;
            }
            else{
                throw std::runtime_error("unsupported fernd model");
            }
        }
        else{
            throw std::runtime_error("unsupported fernd type");
        }
    }
    // fernd io grid
    if (grid_fernd.read_permission or grid_fernd.write_permission or grid_fernd.build_permission){
        ptr = toolkit::tracexml (doc,{"grid","box_fernd"});
        grid_fernd.nx = toolkit::fetchunsigned (ptr,"value","nx",800);
        grid_fernd.ny = toolkit::fetchunsigned (ptr,"value","ny",800);
        grid_fernd.nz = toolkit::fetchunsigned (ptr,"value","nz",160);
        grid_fernd.full_size = grid_fernd.nx*grid_fernd.ny*grid_fernd.nz;
        grid_fernd.x_max = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","x_max",20);
        grid_fernd.x_min = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","x_min",-20);
        grid_fernd.y_max = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","y_max",20);
        grid_fernd.y_min = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","y_min",-20);
        grid_fernd.z_max = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","z_max",4);
        grid_fernd.z_min = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","z_min",-4);
    }
}

void Param::cre_param (tinyxml2::XMLDocument* doc){
    // cre io
    tinyxml2::XMLElement* ptr {toolkit::tracexml (doc,{"fieldio"})};
    if (ptr->FirstChildElement("cre")!=nullptr){
        grid_cre.read_permission = toolkit::fetchbool (ptr,"read","cre");
        grid_cre.write_permission = toolkit::fetchbool (ptr,"write","cre");
        grid_cre.filename = toolkit::fetchstring (ptr,"filename","cre");
    }
    // cre internal
    ptr = toolkit::tracexml(doc,{"cre"});
    grid_cre.build_permission = toolkit::fetchbool (ptr,"cue",0);
    if (grid_cre.build_permission){
        cre_type = ptr->Attribute("type");
        // analytical
        if (cre_type=="analytic"){
            tinyxml2::XMLElement* subptr {toolkit::tracexml(doc,{"cre","analytic"})};
            cre_ana.alpha = toolkit::fetchdouble(subptr,"value","alpha");
            cre_ana.beta = toolkit::fetchdouble(subptr,"value","beta");
            cre_ana.theta = toolkit::fetchdouble(subptr,"value","theta");
            cre_ana.r0 = toolkit::fetchdouble(subptr,"value","r0")*CGS_U_kpc; //kpc
            cre_ana.z0 = toolkit::fetchdouble(subptr,"value","z0")*CGS_U_kpc; //kpc
            cre_ana.E0 = toolkit::fetchdouble(subptr,"value","E0",20.6)*CGS_U_GeV; //GeV
            cre_ana.j0 = toolkit::fetchdouble(subptr,"value","j0",0.0217);
        }
        // uniform
        else if (cre_type=="unif"){
            tinyxml2::XMLElement* subptr {toolkit::tracexml(doc,{"cre","unif"})};
            cre_unif.alpha = toolkit::fetchdouble(subptr,"value","alpha");
            cre_unif.r0 = toolkit::fetchdouble(subptr,"value","r0")*CGS_U_kpc; //kpc
            cre_unif.E0 = toolkit::fetchdouble(subptr,"value","E0",20.6)*CGS_U_GeV; //GeV
            cre_unif.j0 = toolkit::fetchdouble(subptr,"value","j0",0.0217);
        }
        else{
            throw std::runtime_error("unsupported cre model");
        }
    }
    // cre io grid
    if (grid_cre.read_permission or grid_cre.write_permission){
        ptr = toolkit::tracexml (doc,{"grid","box_cre"});
        grid_cre.E_min = CGS_U_GeV*toolkit::fetchdouble (ptr,"value","E_min",0.1);
        grid_cre.E_max = CGS_U_GeV*toolkit::fetchdouble (ptr,"value","E_max",100.0);
        grid_cre.nE = toolkit::fetchunsigned (ptr,"value","nE",80);
        grid_cre.E_fact = std::log(grid_cre.E_max/grid_cre.E_min)/(grid_cre.nE-1);
        grid_cre.nz = toolkit::fetchunsigned (ptr,"value","nz",80);
        grid_cre.nx = toolkit::fetchunsigned (ptr,"value","nx",80);
        grid_cre.ny = toolkit::fetchunsigned (ptr,"value","ny",80);
        grid_cre.x_max = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","x_max",0);
        grid_cre.x_min = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","x_min",0);
        grid_cre.y_max = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","y_max",0);
        grid_cre.y_min = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","y_min",0);
        grid_cre.z_max = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","z_max",4);
        grid_cre.z_min = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","z_min",-4);
        grid_cre.cre_size = grid_cre.nE*grid_cre.nx*grid_cre.ny*grid_cre.nz;
    }
}

// END
