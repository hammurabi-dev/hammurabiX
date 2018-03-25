#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <memory>

#include <param.h>
#include <grid.h>
#include <integrator.h>
#include <breg.h>
#include <brnd.h>
#include <cre.h>
#include <fereg.h>
#include <fernd.h>
#include <namespace_toolkit.h>
#include <tinyxml2.h>
#include <timer.h>

using namespace std;
using namespace tinyxml2;

int main(int /*argc*/,char **argv) {
#ifndef NDEBUG
    Timer tmr;
    tmr.start("main");
#endif
    
    string file_name {argv[1]};
    // BUILD SPECIFIED MODULES
    unique_ptr<Param> par = unique_ptr<Param> (new Param(file_name));
    unique_ptr<Grid_fereg> grid_fereg = unique_ptr<Grid_fereg> (new Grid_fereg(file_name));
    unique_ptr<Grid_breg> grid_breg = unique_ptr<Grid_breg> (new Grid_breg(file_name));
    unique_ptr<Grid_brnd> grid_brnd = unique_ptr<Grid_brnd> (new Grid_brnd(file_name));
    unique_ptr<Grid_fernd> grid_fernd = unique_ptr<Grid_fernd> (new Grid_fernd(file_name));
    unique_ptr<Grid_cre> grid_cre = unique_ptr<Grid_cre> (new Grid_cre(file_name));
    unique_ptr<Grid_int> grid_int = unique_ptr<Grid_int> (new Grid_int(file_name));
    
    unique_ptr<FEreg> fereg;
    unique_ptr<Breg> breg;
    unique_ptr<FErnd> fernd;
    unique_ptr<Brnd> brnd;
    unique_ptr<CRE> cre;
    
    unique_ptr<XMLDocument> doc = toolkit::loadxml(file_name);
    XMLElement *ptr;
    string field_type;
    
    // regular FE field
    ptr = toolkit::tracexml(doc.get(),{"FreeElectron"});
    field_type = toolkit::FetchString(ptr,"type","Regular");
    // if import from file, no need to build specific fe class
    if(grid_fereg->read_permission){
        grid_fereg->import_grid();
        fereg = unique_ptr<FEreg> (new FEreg());
    }
    else if(field_type=="YMW16"){
        fereg = unique_ptr<FEreg> (new FEreg_ymw16());
    }
    else if(field_type=="Verify"){
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
    ptr = toolkit::tracexml(doc.get(),{"MagneticField"});
    field_type = toolkit::FetchString(ptr,"type","Regular");
    if(grid_breg->read_permission){
        grid_breg->import_grid();
        breg = unique_ptr<Breg> (new Breg());
    }
    else if(field_type=="WMAP"){
        breg = unique_ptr<Breg> (new Breg_wmap());
    }
    else if(field_type=="Jaffe"){
        breg = unique_ptr<Breg> (new Breg_jaffe());
    }
    else if(field_type=="Verify"){
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
    else if(grid_fernd->build_permission){
        ptr = toolkit::tracexml(doc.get(),{"FreeElectron"});
        field_type = toolkit::FetchString(ptr,"type","Random");
        if(field_type=="Global"){
            // non default constructor
            fernd = unique_ptr<FErnd> (new FErnd_global());
            // fill grid with random fields
            fernd->write_grid(par.get(),grid_fernd.get());
        }
        else{return EXIT_FAILURE;}
    }
    else{
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
    else if(grid_brnd->build_permission){
        ptr = toolkit::tracexml(doc.get(),{"MagneticField"});
        field_type = toolkit::FetchString(ptr,"type","Random");
        if(field_type=="Global"){
            brnd = unique_ptr<Brnd> (new Brnd_global());
            // fill grid with random fields
            brnd->write_grid(par.get(),breg.get(),grid_breg.get(),grid_brnd.get());
        }
        else if(field_type=="Local"){
            brnd = unique_ptr<Brnd> (new Brnd_local());
            // fill grid with random fields
            brnd->write_grid(par.get(),breg.get(),grid_breg.get(),grid_brnd.get());
        }
        else{return EXIT_FAILURE;}
        
    }
    else{
        //without read permission, return zeros
        brnd = unique_ptr<Brnd> (new Brnd());
    }
    // if export to file
    if(grid_brnd->write_permission){
        grid_brnd->export_grid();
    }
    
    // cre
    ptr = toolkit::tracexml(doc.get(),{});
    field_type = toolkit::FetchString(ptr,"type","CRE");
    if(field_type=="Analytic"){
        cre = unique_ptr<CRE> (new CRE_ana());
    }
    else if(field_type=="Verify"){
        cre = unique_ptr<CRE> (new CRE_verify());
    }
    else if(field_type=="Numeric"){
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
    
    intobj->write_grid(breg.get(),brnd.get(),fereg.get(),fernd.get(),cre.get(),grid_breg.get(),grid_brnd.get(),grid_fereg.get(),grid_fernd.get(),grid_cre.get(),grid_int.get(),par.get());
    
    grid_int->export_grid();
    
    // INSERT MULTI_OUTPUT TAKEOVER LOOP
    ptr = toolkit::tracexml(doc.get(),{"Obsout","Sync"});
    if (ptr!=nullptr){
        for(auto e = ptr->NextSiblingElement("Sync");e!=nullptr;e=e->NextSiblingElement("Sync")){
            // stop producing dm and fd
            grid_int->do_dm = false;
            grid_int->do_fd = false;
            grid_int->do_sync = toolkit::FetchBool(e,"cue");
            par->sim_freq = toolkit::FetchDouble(e,"freq")*CGS_U_GHz;
            grid_int->sim_sync_name = toolkit::FetchString(e,"filename");
            
            intobj->write_grid(breg.get(),brnd.get(),fereg.get(),fernd.get(),cre.get(),grid_breg.get(),grid_brnd.get(),grid_fereg.get(),grid_fernd.get(),grid_cre.get(),grid_int.get(),par.get());
            
            grid_int->export_grid();
        }
    }
    
    // CLEANING
#ifndef NDEBUG
    tmr.stop("main");
    tmr.print();
#endif 
    return EXIT_SUCCESS;
}

// END
