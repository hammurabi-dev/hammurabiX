#include <iostream>
#include <string>
#include <vector>
#include <tinyxml2.h>
#include <cstdlib>
#include <memory>

#include "pond.h"
#include "grid.h"
#include "integrator.h"
#include "breg.h"
#include "brnd.h"
#include "cre.h"
#include "fereg.h"
#include "fernd.h"
#include "namespace_toolkit.h"


using namespace std;
using namespace tinyxml2;

int main(int , char **argv) {
#ifndef NDEBUG
    double time = toolkit::timestamp();
    cout<<"...STARTING HAMMURABI..."<<endl;
#endif
    
    string file_name {argv[1]};
    
    // BUILD SPECIFIED MODULES
    unique_ptr<Pond> par = unique_ptr<Pond> (new Pond(file_name));
    unique_ptr<Grid_fereg> grid_fereg = unique_ptr<Grid_fereg> (new Grid_fereg(file_name));
    unique_ptr<Grid_breg> grid_breg = unique_ptr<Grid_breg> (new Grid_breg(file_name));
    unique_ptr<Grid_brnd> grid_brnd = unique_ptr<Grid_brnd> (new Grid_brnd(file_name));
    unique_ptr<Grid_fernd> grid_fernd = unique_ptr<Grid_fernd> (new Grid_fernd(file_name));
    unique_ptr<Grid_cre> grid_cre = unique_ptr<Grid_cre> (new Grid_cre(file_name));
    unique_ptr<Grid_int> grid_int = unique_ptr<Grid_int> (new Grid_int(file_name));
#ifndef NDEBUG
    cout<<"INFO: ALL GRIDS BUILT"<<endl;
#endif
    unique_ptr<FEreg> fereg;
    unique_ptr<Breg> breg;
    unique_ptr<FErnd> fernd;
    unique_ptr<Brnd> brnd;
    unique_ptr<CRE> cre;
    
    unique_ptr<XMLDocument> doc = unique_ptr<XMLDocument> (new XMLDocument());
    doc->LoadFile(file_name.c_str());
    
    // regular FE field
    string fetype {doc->FirstChildElement("root")->FirstChildElement("Galaxy")->FirstChildElement("FreeElectron")->FirstChildElement("Regular")->Attribute("type")};
    // if import from file, no need to build specific fe class
    if(grid_fereg->read_permission){
        grid_fereg->import_grid();
        fereg = unique_ptr<FEreg> (new FEreg());
    }
    else if(fetype=="YMW16"){
#ifndef NDEBUG
        cout<<"INFO: USING YMW16 FE MODEL"<<endl;
#endif
        fereg = unique_ptr<FEreg> (new FEreg_ymw16());
    }
    else if(fetype=="Verify"){
#ifndef NDEBUG
        cout<<"INFO: USING VERIFICATION FE MODEL"<<endl;
#endif
        fereg = unique_ptr<FEreg> (new FEreg_verify());
    }
    else{return EXIT_FAILURE;}
    // if export to file
    if(grid_fereg->write_permission){
        // write out binary file and exit
        fereg->write_grid(par.get(),grid_fereg.get());
        grid_fereg->export_grid();
    }
    
    // regular B field, must before random B field
    string bregtype {doc->FirstChildElement("root")->FirstChildElement("Galaxy")->FirstChildElement("MagneticField")->FirstChildElement("Regular")->Attribute("type")};
    if(grid_breg->read_permission){
        grid_breg->import_grid();
        breg = unique_ptr<Breg> (new Breg());
    }
    else if(bregtype=="WMAP"){
#ifndef NDEBUG
        cout<<"INFO: USING WMAP3YR REGUALR B MODEL"<<endl;
#endif
        breg = unique_ptr<Breg> (new Breg_wmap());
    }
    else if(bregtype=="Verify"){
#ifndef NDEBUG
        cout<<"INFO: USING VERIFICATION REGUALR B MODEL"<<endl;
#endif
        breg = unique_ptr<Breg> (new Breg_verify());
    }
    else{return EXIT_FAILURE;}
    // if export to file
    if(grid_breg->write_permission){
        breg->write_grid(par.get(),grid_breg.get());
        grid_breg->export_grid();
    }
    
    // random FE field
    // if import from file, no need to build specific fe_rnd class
    if(grid_fernd->read_permission){
        grid_fernd->import_grid();
        fernd = unique_ptr<FErnd> (new FErnd());
    }
    else if(doc->FirstChildElement("root")->FirstChildElement("Galaxy")->FirstChildElement("FreeElectron")->FirstChildElement("Random")->BoolAttribute("cue")){
        string ferndtype {doc->FirstChildElement("root")->FirstChildElement("Galaxy")->FirstChildElement("FreeElectron")->FirstChildElement("Random")->Attribute("type")};
        if(ferndtype=="Iso"){
#ifndef NDEBUG
            cout<<"INFO: USING ISOTROPIC RANDOM FE MODEL"<<endl;
#endif
            // non default constructor
            fernd = unique_ptr<FErnd> (new FErnd_iso());
            // fill grid with random fields
            fernd->write_grid_iso(par.get(),grid_fernd.get());
        }
        else{return EXIT_FAILURE;}
    }
    else{
#ifndef NDEBUG
        cout<<"INFO: NO RANDOM FE FIELD"<<endl;
#endif
        fernd = unique_ptr<FErnd> (new FErnd());
    }
    // if export to file
    if(grid_fernd->write_permission){
        grid_fernd->export_grid();
    }
    
    // random B field
    if(grid_brnd->read_permission){
        grid_brnd->import_grid();
        brnd = unique_ptr<Brnd> (new Brnd());
    }
    else if(doc->FirstChildElement("root")->FirstChildElement("Galaxy")->FirstChildElement("MagneticField")->FirstChildElement("Random")->BoolAttribute("cue")){
        string brndtype {doc->FirstChildElement("root")->FirstChildElement("Galaxy")->FirstChildElement("MagneticField")->FirstChildElement("Random")->Attribute("type")};
        if(brndtype=="Iso"){
#ifndef NDEBUG
            cout<<"INFO: USING ISOTROPIC RANDOM B MODEL"<<endl;
#endif
            brnd = unique_ptr<Brnd> (new Brnd_iso());
            // fill grid with random fields
            brnd->write_grid_iso(par.get(),grid_brnd.get());
        }
        else if(brndtype=="Anisoglob"){
#ifndef NDEBUG
            cout<<"INFO: USING GLOBAL ANISOTROPIC RANDOM B MODEL"<<endl;
#endif
            brnd = unique_ptr<Brnd> (new Brnd_anig());
            // fill grid with random fields
            brnd->write_grid_ani(par.get(),breg.get(),grid_breg.get(),grid_brnd.get());
        }
        else if(brndtype=="Anisolocal"){
#ifndef NDEBUG
            cout<<"INFO: USING LOCAL ANISOTROPIC RANDOM B MODEL"<<endl;
#endif
            brnd = unique_ptr<Brnd> (new Brnd_anil());
            // fill grid with random fields
            brnd->write_grid_ani(par.get(),breg.get(),grid_breg.get(),grid_brnd.get());
        }
        else{return EXIT_FAILURE;}
        
    }
    else{
#ifndef NDEBUG
        cout<<"INFO: NO RANDOM B FIELD"<<endl;
#endif
        //without read permission, return zeros
        brnd = unique_ptr<Brnd> (new Brnd());
    }
    // if export to file
    if(grid_brnd->write_permission){
        grid_brnd->export_grid();
    }
    
    // cre
    string cretype {doc->FirstChildElement("root")->FirstChildElement("CRE")->Attribute("type")};
    if(cretype=="Analytic"){
#ifndef NDEBUG
        cout<<"INFO: USING ANALYTIC CRE"<<endl;
#endif
        cre = unique_ptr<CRE> (new CRE_ana());
    }
    else if(cretype=="Verify"){
#ifndef NDEBUG
        cout<<"INFO: USING VERIFICATION CRE"<<endl;
#endif
        cre = unique_ptr<CRE> (new CRE_verify());
    }
    else if(cretype=="Numeric"){
#ifndef NDEBUG
        cout<<"INFO: USING NUMERIC CRE"<<endl;
#endif
        grid_cre->import_grid();
        cre = unique_ptr<CRE> (new CRE_num());
    }
    else{return EXIT_FAILURE;}
    // if export to file
    if(grid_cre->write_permission){
        cre->write_grid(par.get(),grid_cre.get());
        grid_cre->export_grid();
    }
    
    // START INTEGRATION
    unique_ptr<Integrator> intobj = unique_ptr<Integrator> (new Integrator());
#ifndef NDEBUG
    cout<<"...ALL MODULES BUILT..."<<endl;
#endif
    
    intobj->write_grid(breg.get(),brnd.get(),fereg.get(),fernd.get(),cre.get(),grid_breg.get(),grid_brnd.get(),grid_fereg.get(),grid_fernd.get(),grid_cre.get(),grid_int.get(),par.get());
#ifndef NDEBUG
    cout<<"...PRODUCING MAPS..."<<endl;
#endif
    grid_int->export_grid();
    
    // INSERT MULTI_OUTPUT TAKEOVER LOOP
    if (doc->FirstChildElement("root")->FirstChildElement("Output")->FirstChildElement("Sync")!=nullptr){
        for(auto e = doc->FirstChildElement("root")->FirstChildElement("Output")->FirstChildElement("Sync")->NextSiblingElement("Sync");e!=NULL;e=e->NextSiblingElement("Sync")){
            // stop producing dm and fd
            grid_int->do_dm = false;
            grid_int->do_fd = false;
            grid_int->do_sync = e->BoolAttribute("cue");
            par->sim_freq = e->DoubleAttribute("freq")*CGS_U_GHz;
            grid_int->sim_sync_name = e->Attribute("filename");
#ifndef NDEBUG
            cout<<"...MULTI OUTPUT MODE..."<<endl;
#endif
            intobj->write_grid(breg.get(),brnd.get(),fereg.get(),fernd.get(),cre.get(),grid_breg.get(),grid_brnd.get(),grid_fereg.get(),grid_fernd.get(),grid_cre.get(),grid_int.get(),par.get());
#ifndef NDEBUG
            cout<<"...PRODUCING MAPS..."<<endl;
#endif
            grid_int->export_grid();
        }
    }
    
    // CLEANING
#ifndef NDEBUG
    cout<<"...ENDING HAMMURABI..."<<endl
    <<"INFO:TIME ELAPSE "<<(toolkit::timestamp()-time)<<"sec"<<endl;
#endif 
    return EXIT_SUCCESS;
}

// END
