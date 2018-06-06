#include <iostream>
#include <memory>
#include <array>
#include <string>
#include <vector>
#include <tinyxml2.h>
#include <fitshandle.h>
#include <fstream>
#include <healpix_map_fitsio.h>
#include <fitsio.h>
#include <grid.h>
#include <cgs_units_file.h>
#include <namespace_toolkit.h>
#include <cassert>


/* line of sight integrator */
Grid_int::Grid_int(const std::string &file_name){
    std::unique_ptr<tinyxml2::XMLDocument> doc = toolkit::loadxml(file_name);
    build_grid(doc.get());
}

void Grid_int::build_grid(tinyxml2::XMLDocument *doc){
    tinyxml2::XMLElement *ptr {toolkit::tracexml(doc,{"Grid","Shell"})};
    // get shell Nside
    std::string shell_type {toolkit::FetchString(ptr,"type","layer")};
    if(shell_type=="auto"){
        tinyxml2::XMLElement *subptr {toolkit::tracexml(doc,{"Grid","Shell","layer","auto"})};
        total_shell = toolkit::FetchUnsigned(subptr,"value","shell_num");
        size_t nside_min {toolkit::FetchUnsigned(subptr,"value","nside_min")};
        for(std::size_t i=0;i!=total_shell;++i){
            nside_shell.push_back(pow(2,i)*nside_min);
        }
    }
    else if(shell_type=="manual"){
        tinyxml2::XMLElement *subptr {toolkit::tracexml(doc,{"Grid","Shell","layer","manual"})};
        total_shell = 0;
        for(auto e = subptr->FirstChildElement("nside");e!=nullptr;e=e->NextSiblingElement("nside")){
            total_shell++;
            nside_shell.push_back(toolkit::FetchUnsigned(e,"value"));
        }
        assert(nside_shell.size()==total_shell);
    }
    else{
        assert(false);
    }
    nside_sim = toolkit::FetchUnsigned(ptr,"value","nside_sim");
    npix_sim = 12*nside_sim*nside_sim;
    ec_r_max = CGS_U_kpc*toolkit::FetchDouble(ptr,"value","ec_r_max");
    gc_r_max = toolkit::FetchDouble(ptr,"value","gc_r_max")*CGS_U_kpc;
    gc_z_max = toolkit::FetchDouble(ptr,"value","gc_z_max")*CGS_U_kpc;
    radial_res = toolkit::FetchDouble(ptr,"value","ec_r_res")*CGS_U_kpc;
    lat_lim = toolkit::FetchDouble(ptr,"value","lat_lim")*CGS_U_rad;
    // output file name
    
    ptr = toolkit::tracexml(doc,{"Obsout"});
    if(ptr->FirstChildElement("DM")!=nullptr){
        do_dm = toolkit::FetchBool(ptr,"cue","DM");
        sim_dm_name = toolkit::FetchString(ptr,"filename","DM");
    }
    else{
        do_dm = false;
    }
    if(ptr->FirstChildElement("Faraday")!=nullptr){
        do_fd = toolkit::FetchBool(ptr,"cue","Faraday");
        sim_fd_name = toolkit::FetchString(ptr,"filename","Faraday");
    }
    else{
        do_fd = false;
    }
    if(ptr->FirstChildElement("Sync")!=nullptr){
        do_sync = toolkit::FetchBool(ptr,"cue","Sync");
        sim_sync_name = toolkit::FetchString(ptr,"filename","Sync");
    }
    else{
        do_sync = false;
    }
    assert(do_dm or do_fd or do_sync);
}

void Grid_int::export_grid(void){
    if(do_dm){
        // in units pc/cm^3, conventional units
        dm_map.Scale(CGS_U_ccm/CGS_U_pc);
        fitshandle out_dm;
        if(std::ifstream(sim_dm_name.c_str())){
            remove(sim_dm_name.c_str());
        }
        out_dm.create(sim_dm_name);
        write_Healpix_map_to_fits(out_dm,
                                  dm_map,
                                  PLANCK_FLOAT64);
        out_dm.close();
    }
    if(do_sync){
        fitshandle out_sc;
        if(std::ifstream(sim_sync_name.c_str())){
            remove(sim_sync_name.c_str());
        }
        out_sc.create(sim_sync_name);
        write_Healpix_map_to_fits(out_sc,
                                  Is_map,
                                  Qs_map,
                                  Us_map,
                                  PLANCK_FLOAT64);
        out_sc.close();
    }
    if(do_fd){
        // FD units is rad*m^(-2) in our calculation
        fd_map.Scale(CGS_U_m*CGS_U_m);
        fitshandle out_fd;
        if(std::ifstream(sim_fd_name.c_str())){
            remove(sim_fd_name.c_str());
        }
        out_fd.create(sim_fd_name);
        write_Healpix_map_to_fits(out_fd,
                                  fd_map,
                                  PLANCK_FLOAT64);
        out_fd.close();
    }
}
