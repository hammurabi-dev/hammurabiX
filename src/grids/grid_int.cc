#include <sstream>
#include <iostream>
#include <memory>
#include <fftw3.h>
#include <array>
#include <string>
#include <vector>
#include <tinyxml2.h>
#include <fitshandle.h>
#include <fstream>
#include <healpix_map_fitsio.h>
#include <fitsio.h>
#include "grid.h"
#include "cgs_units_file.h"
#include "namespace_toolkit.h"

using namespace tinyxml2;
using namespace std;

/* line of sight integrator */
Grid_int::Grid_int(string file_name){
    unique_ptr<XMLDocument> doc = unique_ptr<XMLDocument> (new XMLDocument());
    doc->LoadFile(file_name.c_str());
    build_grid(doc.get());
}

void Grid_int::build_grid(XMLDocument *doc){
    XMLElement *ptr {doc->FirstChildElement("root")->FirstChildElement("Grid")->FirstChildElement("Integration")};
    // get shell Nside
    string shell_type {ptr->FirstChildElement("shell")->Attribute("type")};
    if(shell_type=="auto"){
        XMLElement *subptr {ptr->FirstChildElement("shell")->FirstChildElement("auto")};
        total_shell = FetchUnsigned(subptr,"shell_num");
        unsigned int nside_min {FetchUnsigned(subptr,"nside_min")};
        for(unsigned int i=0;i!=total_shell;++i){
            nside_shell.push_back(pow(2,i)*nside_min);
        }
    }
    else if(shell_type=="manual"){
        XMLElement *subptr {ptr->FirstChildElement("shell")->FirstChildElement("manual")};
        total_shell = 0;
        for(auto e = subptr->FirstChildElement("nside");e!=NULL;e=e->NextSiblingElement("nside")){
            total_shell++;
            nside_shell.push_back(e->UnsignedAttribute("value"));
        }
#ifndef NDEBUG
        if(nside_shell.size()!=total_shell){
            cerr<<"ERR:"<<__FILE__
            <<" : in function "<<__func__<<endl
            <<" at line "<<__LINE__<<endl
            <<"WRONG SHELL NUM"<<endl;
            exit(1);
        }
#endif
    }
    else{
        cerr<<"ERR:"<<__FILE__
        <<" : in function "<<__func__<<endl
        <<" at line "<<__LINE__<<endl
        <<"INVALID SHELL TYPE"<<endl;
        exit(1);
    }
    nside_sim = FetchUnsigned(ptr,"nside_sim");
    npix_sim = 12*nside_sim*nside_sim;
    ec_r_max = CGS_U_kpc*FetchDouble(ptr,"ec_r_max");
    gc_r_max = FetchDouble(ptr,"gc_r_max")*CGS_U_kpc;
    gc_z_max = FetchDouble(ptr,"gc_z_max")*CGS_U_kpc;
    radial_res = FetchDouble(ptr,"ec_r_res")*CGS_U_kpc;
    lat_lim = FetchDouble(ptr,"lat_lim")*CGS_U_pi/180.;
    // output file name
    
    ptr = doc->FirstChildElement("root")->FirstChildElement("Output");
    cout<<"Starting"<<endl;
	cout<<ptr->FirstChildElement("DM")<<endl;
	
	if (ptr->FirstChildElement("DM") != 0) {
		do_dm = ptr->FirstChildElement("DM")->BoolAttribute("cue");
		sim_dm_name = ptr->FirstChildElement("DM")->Attribute("filename");
		cout<<"dm"<<endl;
	}
	else
	{
		do_dm = false;
		cout<<"no_dm"<<endl;
	}
	
	if (ptr->FirstChildElement("Faraday") != NULL) {
		do_fd = ptr->FirstChildElement("Faraday")->BoolAttribute("cue");
		sim_fd_name = ptr->FirstChildElement("Faraday")->Attribute("filename");
		cout<<"fd"<<endl;
	}
	else
	{
		do_fd = false;
		cout<<"no_fd"<<endl;		
	}
	
	if (ptr->FirstChildElement("Sync") != NULL) {
		do_sync = ptr->FirstChildElement("Sync")->BoolAttribute("cue");
		sim_sync_name = ptr->FirstChildElement("Sync")->Attribute("filename");
		cout<<"sync"<<endl;
	}
	else
	{
		do_sync = false;
		cout<<"no_sync"<<endl;		
	}
    
    if (not (do_dm or do_fd or do_sync)) {
        cerr<<"PEACE: NO OUTPUT REQUIRED"<<endl;
        exit(1);
    }
    
    

    
}

void Grid_int::export_grid(void){
    if(do_dm){
        // in units pc/cm^3, conventional units
        dm_map.Scale(CGS_U_ccm/CGS_U_pc);
        fitshandle out_dm;
        std::stringstream outfname_dm;
        // cleanup content of outfname
        outfname_dm.str("");
        // give content to outfname
        outfname_dm << sim_dm_name;
        std::ifstream in_dm(outfname_dm.str().data());
        if(in_dm.is_open()){
            out_dm.delete_file(outfname_dm.str());
        }
        out_dm.create(outfname_dm.str());
        write_Healpix_map_to_fits(out_dm, dm_map, PLANCK_FLOAT64);
        out_dm.close();
    }
    if(do_sync){
        fitshandle out_sc;
        std::stringstream outfname_sc;
        outfname_sc.str("");
        outfname_sc << sim_sync_name;
        std::ifstream in_sc(outfname_sc.str().data());
        if (in_sc.is_open()) {
            out_sc.delete_file(outfname_sc.str());
        }
        out_sc.create(outfname_sc.str());
        write_Healpix_map_to_fits(out_sc, Is_map, Qs_map, Us_map, PLANCK_FLOAT64);
        out_sc.close();
    }
    if(do_fd){
        // FD units is rad*m^(-2) in our calculation
        fd_map.Scale(CGS_U_m*CGS_U_m);
        fitshandle out_fd;
        std::stringstream outfname_fd;
        outfname_fd.str("");
        outfname_fd << sim_fd_name;
        std::ifstream in_fd(outfname_fd.str().data());
        if (in_fd.is_open()) {
            out_fd.delete_file(outfname_fd.str());
        }
        out_fd.create(outfname_fd.str());
        write_Healpix_map_to_fits(out_fd, fd_map, PLANCK_FLOAT64);
        out_fd.close();
    }
}
