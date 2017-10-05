/* hammurabi stand alone routine
 *@author: Jiaxin Wang
 *@email: jiwang@sissa.it
 */
#include <iostream>
#include <string>
#include <vector>
#include <tinyxml2.h>
#include <cstdlib>
#include <memory>

#include "class_pond.h"
#include "class_grid.h"
#include "class_int.h"
#include "class_breg.h"
#include "class_brnd.h"
#include "class_cre.h"
#include "class_fe.h"
#include "class_fernd.h"
#include "namespace_toolkit.h"

using namespace std;
using namespace tinyxml2;

int main(int , char **argv) {
    toolkit::timestamp();
    cout<<"...STARTING HAMMURABI..."<<endl;
    string file_name {argv[1]};
    
    // BUILD SPECIFIED MODULES
    Pond *par = new Pond(file_name);
    Grid_fe *grid_fe = new Grid_fe(file_name);
    Grid_breg *grid_breg = new Grid_breg(file_name);
    Grid_brnd *grid_brnd = new Grid_brnd(file_name);
    Grid_fernd *grid_fernd = new Grid_fernd(file_name);
    Grid_cre *grid_cre = new Grid_cre(file_name);
    Grid_int *grid_int = new Grid_int(file_name);
    cout<<"INFO: ALL GRIDS BUILT"<<endl;
    XMLDocument *doc = new XMLDocument();
    doc->LoadFile(file_name.c_str());
    FE *fe; Breg *breg; FErnd *fernd; Brnd *brnd; CRE *cre;
    // reminder for cleaning up allocated grid memory
    bool fe_reminder {false};
    bool breg_reminder {false};
    bool fernd_reminder {false};
    bool brnd_reminder {false};
    bool cre_reminder {false};
    
    // regular FE field
    string fetype {doc->FirstChildElement("root")->FirstChildElement("Galaxy")->FirstChildElement("FreeElectron")->Attribute("type")};
    // if import from file, no need to build specific fe class
    if(grid_fe->read_permission){
        grid_fe->import_grid();
        fe = new FE();
        fe_reminder = true;
    }
    else if(fetype=="YMW16"){
        cout<<"INFO: USING YMW16 FE MODEL"<<endl;
        fe = new YMW16();
    }
    else{return EXIT_FAILURE;}
    // if export to file
    if(grid_fe->write_permission){
        // write out binary file and exit
        fe->write_grid(par,grid_fe);
        grid_fe->export_grid();
    }
    
    // regular B field, must before random B field
    string bregtype {doc->FirstChildElement("root")->FirstChildElement("Galaxy")->FirstChildElement("MagneticField")->FirstChildElement("Regular")->Attribute("type")};
    if(grid_breg->read_permission){
        grid_breg->import_grid();
        breg = new Breg();
        breg_reminder = true;
    }
    else if(bregtype=="WMAP"){
        cout<<"INFO: USING WMAP3YR REGUALR B MODEL"<<endl;
        breg = new Bwmap();
    }
    else if(bregtype=="Local"){
        cout<<"INFO: USING LOCAL REGULAR B TESTING MODEL"<<endl;
        breg = new Blocal();
    }
    else{return EXIT_FAILURE;}
    // if export to file
    if(grid_breg->write_permission){
        breg->write_grid(par,grid_breg);
        grid_breg->export_grid();
    }
    
    // random FE field
    // if import from file, no need to build specific fe_rnd class
    if(grid_fernd->read_permission){
        grid_fernd->import_grid();
        fernd = new FErnd();
        fernd_reminder = true;
    }
    else if(doc->FirstChildElement("root")->FirstChildElement("Galaxy")->FirstChildElement("FreeElectron")->FirstChildElement("Random")->BoolAttribute("cue")){
        string ferndtype {doc->FirstChildElement("root")->FirstChildElement("Galaxy")->FirstChildElement("FreeElectron")->FirstChildElement("Random")->Attribute("type")};
        if(ferndtype=="Iso"){
            cout<<"INFO: USING ISOTROPIC RANDOM FE MODEL"<<endl;
            // non default constructor
            fernd = new FErnd_iso();
            // fill grid with random fields
            fernd->write_grid_iso(par,grid_fernd);
            fernd_reminder = true;
        }
        else{return EXIT_FAILURE;}
    }
    else{
        cout<<"INFO: NO RANDOM FE FIELD"<<endl;
        fernd = new FErnd();
    }
    // if export to file
    if(grid_fernd->write_permission){
        grid_fernd->export_grid();
    }
    
    // random B field
    if(grid_brnd->read_permission){
        grid_brnd->import_grid();
        brnd = new Brnd();
        brnd_reminder = true;
    }
    else if(doc->FirstChildElement("root")->FirstChildElement("Galaxy")->FirstChildElement("MagneticField")->FirstChildElement("Random")->BoolAttribute("cue")){
        string brndtype {doc->FirstChildElement("root")->FirstChildElement("Galaxy")->FirstChildElement("MagneticField")->FirstChildElement("Random")->Attribute("type")};
        if(brndtype=="Iso"){
            cout<<"INFO: USING ISOTROPIC RANDOM B MODEL"<<endl;
            brnd = new Brnd_iso();
            // fill grid with random fields
            brnd->write_grid_iso(par,grid_brnd);
            brnd_reminder = true;
        }
        else if(brndtype=="Anisoglob"){
            cout<<"INFO: USING GLOBAL ANISOTROPIC RANDOM B MODEL"<<endl;
            brnd = new Brnd_anig();
            // fill grid with random fields
            brnd->write_grid_ani(par,breg,grid_breg,grid_brnd);
            brnd_reminder = true;
        }
        else{return EXIT_FAILURE;}
        
    }
    else{
        cout<<"INFO: NO RANDOM B FIELD"<<endl;
        //without read permission, return zeros
        brnd = new Brnd();
    }
    // if export to file
    if(grid_brnd->write_permission){
        grid_brnd->export_grid();
    }
    
    // cre
    string cretype {doc->FirstChildElement("root")->FirstChildElement("CRE")->Attribute("type")};
    if(cretype=="Analytic"){
        cout<<"INFO: USING ANALYTIC CRE"<<endl;
        cre = new CRE_ana();
    }
    else if(cretype=="Numeric"){
        cout<<"INFO: USING NUMERIC CRE"<<endl;
        grid_cre->import_grid();
        cre = new CRE_num();
        cre_reminder = true;
    }
    else{return EXIT_FAILURE;}
    // if export to file
    if(grid_cre->write_permission){
        cre->write_grid(par,grid_cre);
        grid_cre->export_grid();
    }
    
    // START CALCULATION
    Integrator *intobj = new Integrator();
    cout<<"...ALL MODULES BUILT..."<<endl;
    
    intobj->write_grid(breg,brnd,fe,fernd,cre,grid_breg,grid_brnd,grid_fe,grid_fernd,grid_cre,grid_int,par);
    cout<<"...PRODUCING MAPS..."<<endl;
    grid_int->export_grid();
    
    
    // INSERT MULTI_OUTPUT TAKEOVER LOOP
    for(auto e = doc->FirstChildElement("root")->FirstChildElement("Output")->FirstChildElement("Sync")->NextSiblingElement("Sync");e!=NULL;e=e->NextSiblingElement("Sync")){
        // stop producing dm and fd
        grid_int->do_dm = false;
        grid_int->do_fd = false;
        par->sim_freq = e->DoubleAttribute("freq")*CGS_U_GHz;
        grid_int->sim_sync_name = e->Attribute("filename");
        
        cout<<"...MULTI OUTPUT MODE..."<<endl;
        intobj->write_grid(breg,brnd,fe,fernd,cre,grid_breg,grid_brnd,grid_fe,grid_fernd,grid_cre,grid_int,par);
        cout<<"...PRODUCING MAPS..."<<endl;
        grid_int->export_grid();
    }
    
    
    // CLEANING
    cout<<"...ENDING HAMMURABI..."<<endl;
    if(breg_reminder) grid_breg->clean_grid();
    if(brnd_reminder) grid_brnd->clean_grid();
    if(fe_reminder) grid_fe->clean_grid();
    if(fernd_reminder) grid_fernd->clean_grid();
    if(cre_reminder) grid_cre->clean_grid();
    delete par; delete grid_breg; delete grid_brnd;
    delete grid_fe;delete grid_fernd; delete grid_cre; delete grid_int;
    delete fe;delete fernd; delete breg; delete brnd;
    delete cre; delete doc; delete intobj;
    
    toolkit::timestamp();
    return EXIT_SUCCESS;
}

// END
