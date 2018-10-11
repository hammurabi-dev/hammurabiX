#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <memory>
#include <cassert>
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

class Pipeline{
public:
    Pipeline (std::string);
    void assemble_grid();
    void assemble_fereg();
    void assemble_breg();
    void assemble_fernd();
    void assemble_brnd();
    void assemble_cre();
    void assemble_obs();
private:
    std::string file_name;
    std::unique_ptr<tinyxml2::XMLDocument> doc;
    std::unique_ptr<Param> par;
    std::unique_ptr<Grid_fereg> grid_fereg;
    std::unique_ptr<Grid_breg> grid_breg;
    std::unique_ptr<Grid_brnd> grid_brnd;
    std::unique_ptr<Grid_fernd> grid_fernd;
    std::unique_ptr<Grid_cre> grid_cre;
    std::unique_ptr<Grid_int> grid_int;
    std::unique_ptr<FEreg> fereg;
    std::unique_ptr<Breg> breg;
    std::unique_ptr<FErnd> fernd;
    std::unique_ptr<Brnd> brnd;
    std::unique_ptr<CRE> cre;
    std::unique_ptr<Integrator> intobj;
};

// constructor
Pipeline::Pipeline(std::string name) {
    file_name = name;
    doc = toolkit::loadxml(file_name);
}

void Pipeline::assemble_grid (){
    par = std::make_unique<Param> (file_name);
    grid_fereg = std::make_unique<Grid_fereg> (file_name);
    grid_breg = std::make_unique<Grid_breg> (file_name);
    grid_brnd = std::make_unique<Grid_brnd> (file_name);
    grid_fernd = std::make_unique<Grid_fernd> (file_name);
    grid_cre = std::make_unique<Grid_cre> (file_name);
    grid_int = std::make_unique<Grid_int> (file_name);
}

// regular FE field
void Pipeline::assemble_fereg (){
    if (!grid_fereg->build_permission){
        fereg = std::make_unique<FEreg>();
        return;
    }
    tinyxml2::XMLElement *ptr {toolkit::tracexml(doc.get(),{"FreeElectron"})};
    std::string field_type {toolkit::fetchstring(ptr,"type","Regular")};
    // if import from file, no need to build specific fe class
    if (grid_fereg->read_permission){
        grid_fereg->import_grid();
        fereg = std::make_unique<FEreg>();
    }
    else if (field_type=="YMW16"){
        fereg = std::make_unique<FEreg_ymw16>();
    }
#ifndef NDEBUG
    else if (field_type=="Test"){
        fereg = std::make_unique<FEreg_test>();
    }
#endif
    else assert (false);
    // if export to file
    if (grid_fereg->write_permission){
        // write out binary file and exit
        fereg->write_grid (par.get(),grid_fereg.get());
        grid_fereg->export_grid();
    }
}

// regular B field, must before random B field
void Pipeline::assemble_breg (){
    if (!grid_breg->build_permission){
        breg = std::make_unique<Breg>();
        return;
    }
    tinyxml2::XMLElement *ptr {toolkit::tracexml(doc.get(),{"MagneticField"})};
    std::string field_type {toolkit::fetchstring(ptr,"type","Regular")};
    if (grid_breg->read_permission){
        grid_breg->import_grid();
        breg = std::make_unique<Breg>();
    }
    else if (field_type=="WMAP"){
        breg = std::make_unique<Breg_wmap>();
    }
    else if (field_type=="Jaffe"){
        breg = std::make_unique<Breg_jaffe>();
    }
#ifndef NDEBUG
    else if (field_type=="Test"){
        breg = std::make_unique<Breg_test>();
    }
#endif
    else assert (false);
    // if export to file
    if (grid_breg->write_permission){
        breg->write_grid(par.get(),grid_breg.get());
        grid_breg->export_grid();
    }
}

// random FE field
void Pipeline::assemble_fernd (){
    // if import from file, no need to build specific fe_rnd class
    if (grid_fernd->read_permission){
        grid_fernd->import_grid();
        fernd = std::make_unique<FErnd>();
    }
    else if (grid_fernd->build_permission){
        tinyxml2::XMLElement *ptr {toolkit::tracexml(doc.get(),{"FreeElectron","Random"})};
        std::string field_type {toolkit::fetchstring(ptr,"type")};
        if (field_type=="Global"){
            std::string field_method {toolkit::fetchstring(ptr,"type","Global")};
            if (field_method=="DFT"){
                fernd = std::make_unique<FErnd_dft>();
            }
            // fill grid with random fields
            fernd->write_grid (par.get(),grid_fernd.get());
        }
        else assert (false);
    }
    else {
        fernd = std::make_unique<FErnd>();
    }
    // if export to file
    if (grid_fernd->write_permission){
        grid_fernd->export_grid();
    }
}

// random B field
void Pipeline::assemble_brnd (){
    if (grid_brnd->read_permission){
        grid_brnd->import_grid();
        brnd = std::make_unique<Brnd>();
    }
    else if (grid_brnd->build_permission){
        tinyxml2::XMLElement *ptr {toolkit::tracexml(doc.get(),{"MagneticField","Random"})};
        std::string field_type {toolkit::fetchstring(ptr,"type")};
        if (field_type=="Global"){
            std::string field_method {toolkit::fetchstring(ptr,"type","Global")};
            if (field_method=="ES"){
                brnd = std::make_unique<Brnd_es>();
            }
            //else if (field_method=="Jaffe"){
            // to be implemented
            //brnd = std::make_unique<Brnd_jaffe>();
            //}
            // fill grid with random fields
            brnd->write_grid (par.get(),
                              breg.get(),
                              grid_breg.get(),
                              grid_brnd.get());
        }
        else if (field_type=="Local"){
            std::string field_method {toolkit::fetchstring(ptr,"type","Local")};
            if (field_method=="MHD"){
                brnd = std::make_unique<Brnd_mhd>();
            }
            // fill grid with random fields
            brnd->write_grid (par.get(),
                              breg.get(),
                              grid_breg.get(),
                              grid_brnd.get());
        }
        else assert (false);
    }
    else {
        //without read permission, return zeros
        brnd = std::make_unique<Brnd>();
    }
    // if export to file
    if (grid_brnd->write_permission){
        grid_brnd->export_grid();
    }
}

// cre
void Pipeline::assemble_cre (){
    tinyxml2::XMLElement *ptr {toolkit::tracexml(doc.get(),{})};
    std::string field_type {toolkit::fetchstring(ptr,"type","CRE")};
    if (field_type=="Analytic"){
        cre = std::make_unique<CRE_ana>();
    }
#ifndef NDEBUG
    else if (field_type=="Test"){
        cre = std::make_unique<CRE_test>();
    }
#endif
    else if (field_type=="Numeric"){
        grid_cre->import_grid();
        cre = std::make_unique<CRE_num>();
    }
    else assert (false);
    // if export to file
    if (grid_cre->write_permission){
        cre->write_grid (par.get(),
                         grid_cre.get());
        grid_cre->export_grid();
    }
}
// LOS integration for observables
void Pipeline::assemble_obs (){
    intobj = std::make_unique<Integrator>();
    intobj->write_grid (breg.get(),
                        brnd.get(),
                        fereg.get(),
                        fernd.get(),
                        cre.get(),
                        grid_breg.get(),
                        grid_brnd.get(),
                        grid_fereg.get(),
                        grid_fernd.get(),
                        grid_cre.get(),
                        grid_int.get(),
                        par.get());
    grid_int->export_grid();
    // INSERT MULTI_OUTPUT TAKEOVER LOOP
    tinyxml2::XMLElement *ptr {toolkit::tracexml(doc.get(),{"Obsout","Sync"})};
    if (ptr!=nullptr){
        for (auto e = ptr->NextSiblingElement("Sync");e!=nullptr;e=e->NextSiblingElement("Sync")){
            // stop producing dm and fd
            grid_int->do_dm = false;
            grid_int->do_fd = false;
            grid_int->do_sync = toolkit::fetchbool(e,"cue");
            par->sim_freq = toolkit::fetchdouble(e,"freq")*CGS_U_GHz;
            grid_int->sim_sync_name = toolkit::fetchstring(e,"filename");
            intobj->write_grid (breg.get(),
                                brnd.get(),
                                fereg.get(),
                                fernd.get(),
                                cre.get(),
                                grid_breg.get(),
                                grid_brnd.get(),
                                grid_fereg.get(),
                                grid_fernd.get(),
                                grid_cre.get(),
                                grid_int.get(),
                                par.get());
            grid_int->export_grid();
        }
    }
}

int main (int /*argc*/,char **argv) {
#ifndef NTIMING
    Timer tmr;
    tmr.start ("main");
#endif
    std::string xmlfile {argv[1]};
    Pipeline run (xmlfile);
    run.assemble_grid();
    run.assemble_fereg();
    run.assemble_breg();
    run.assemble_fernd();
    run.assemble_brnd();
    run.assemble_cre();
    run.assemble_obs();
#ifndef NTIMING
    tmr.stop ("main");
    tmr.print();
#endif 
    return EXIT_SUCCESS;
}

// END
