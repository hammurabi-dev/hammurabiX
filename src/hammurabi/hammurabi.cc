/* hammurabi stand alone routine
 *@author: Jiaxin Wang
 *@email: jiwang@sissa.it
 */
#include <iostream>
#include <string>
#include <vector>
#include <tinyxml2.h>
#include <cstdlib>

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

int main(int argc, char **argv) {
    toolkit::timestamp();
    cout<<"...STARTING HAMMURABI..."<<endl;
    string file_name(argv[1]);
    
    /* BUILD SPECIFIED MODULES */
    Pond *par = new Pond(file_name);
    Grid_breg *grid_breg = new Grid_breg(file_name);
    Grid_brnd *grid_brnd = new Grid_brnd(file_name);
    Grid_fe *grid_fe = new Grid_fe(file_name);
    Grid_fernd *grid_fernd = new Grid_fernd(file_name);
    Grid_cre *grid_cre = new Grid_cre(file_name);
    Grid_int *grid_int = new Grid_int(file_name);
    cout<<"INFO: ALL GRIDS BUILT"<<endl;
    
    XMLDocument *doc = new XMLDocument();
    doc->LoadFile(file_name.c_str());
    
    FE *fe = new FE();
    FErnd *fernd = new FErnd();
    Breg *breg = new Breg();
    Brnd *brnd = new Brnd();
    CRE *cre = new CRE();
    // fe
    string fetype = doc->FirstChildElement("root")->FirstChildElement("Galaxy")->FirstChildElement("FreeElectron")->Attribute("type");
    if(fetype=="YMW16"){
        cout<<"INFO: USING YMW16 FE MODEL"<<endl;
        fe = new YMW16();
        // write out binary file and exit
        //fe->write_grid(par,grid_fe);
        //grid_fe->export_grid();
    }
    else{return EXIT_FAILURE;}
    // random FE field
    if(doc->FirstChildElement("root")->FirstChildElement("Galaxy")->FirstChildElement("FreeElectron")->FirstChildElement("Random")->BoolAttribute("cue")){
        string ferndtype = doc->FirstChildElement("root")->FirstChildElement("Galaxy")->FirstChildElement("FreeElectron")->FirstChildElement("Random")->Attribute("type");
        if(ferndtype=="Gauss"){
            cout<<"INFO: USING GAUSSIAN RANDOM FE MODEL"<<endl;
            // non default constructor
            fernd = new FEgrnd(par,grid_fernd);
            // fill grid with random fields
            fernd->write_grid(par,grid_fernd);
        }
        else{return EXIT_FAILURE;}
    }
    else{
        cout<<"INFO: NO RANDOM FE FIELD"<<endl;
    }
    
    // regular B field, must before random B field
    string bregtype = doc->FirstChildElement("root")->FirstChildElement("Galaxy")->FirstChildElement("MagneticField")->FirstChildElement("Regular")->Attribute("type");
    if(bregtype=="WMAP"){
        cout<<"INFO: USING WMAP3YR REGUALR B MODEL"<<endl;
        breg = new Bwmap();
    }
    else if(bregtype=="Local"){
        cout<<"INFO: USING LOCAL REGULAR B TESTING MODEL"<<endl;
        breg = new Blocal();
    }
    else{return EXIT_FAILURE;}
    
    // random B field
    if(doc->FirstChildElement("root")->FirstChildElement("Galaxy")->FirstChildElement("MagneticField")->FirstChildElement("Random")->BoolAttribute("cue")){
        string brndtype = doc->FirstChildElement("root")->FirstChildElement("Galaxy")->FirstChildElement("MagneticField")->FirstChildElement("Random")->Attribute("type");
        if(brndtype=="Gauss"){
            cout<<"INFO: USING GAUSSIAN RANDOM B MODEL"<<endl;
            brnd = new Bgrnd(par,grid_brnd);
            // fill grid with random fields
            brnd->write_grid(par,grid_brnd);
        }
        else if(brndtype=="GaussPlus"){
            cout<<"INFO: USING GAUSSIAN PLUS ANISOTROPIC RANDOM B MODEL"<<endl;
            brnd = new Bfrnd(par,grid_brnd);
            // fill grid with random fields
            brnd->write_grid_plus(par,breg,grid_breg,grid_brnd);
        }
        else{return EXIT_FAILURE;}
        
    }
    else{
        cout<<"INFO: NO RANDOM B FIELD"<<endl;
    }
    
    // cre
    string cretype = doc->FirstChildElement("root")->FirstChildElement("CRE")->Attribute("type");
    if(cretype=="Analytic"){
        cout<<"INFO: USING ANALYTICAL CRE MODEL"<<endl;
        cre = new CRE_ana();
    }
    else if(cretype=="Dragon"){
        cout<<"INFO: USING CRE FROM DRAGON"<<endl;
        cre = new CRE_num();
    }
    else{return EXIT_FAILURE;}
    
    /* IMPORT/EXPORT BINARY FILES */
    if(doc->FirstChildElement("root")->FirstChildElement("Interface")->FirstChildElement("fe_grid")->BoolAttribute("cue")){
        //fe->write_grid(par,grid_fe);
        //grid_fe->export_grid();
        grid_fe->import_grid();
    }
    if(doc->FirstChildElement("root")->FirstChildElement("Interface")->FirstChildElement("breg_grid")->BoolAttribute("cue")){
        //breg->write_grid(par,grid_breg);
        //grid_breg->export_grid();
        grid_breg->import_grid();
    }
    if(doc->FirstChildElement("root")->FirstChildElement("Interface")->FirstChildElement("fernd_grid")->BoolAttribute("cue")){
        // DESIGNED FOR UNUSUAL NEED
        //grid_fernd->export_grid();
        //grid_fernd->export_grid();
    }
    if(doc->FirstChildElement("root")->FirstChildElement("Interface")->FirstChildElement("brnd_grid")->BoolAttribute("cue")){
        // DESIGNED FOR UNUSUAL NEED
        //grid_brnd->export_grid();
        //grid_brnd->import_grid();
    }
    /* START CALCULATION */
    
    Integrator *intobj = new Integrator();
    cout<<"...ALL MODULES BUILT..."<<endl;
    
    intobj->write_grid(breg,brnd,fe,fernd,cre,grid_breg,grid_brnd,grid_fe,grid_fernd,grid_cre,grid_int,par);
    cout<<"...PRODUCING MAPS..."<<endl;
    grid_int->export_grid();
    
    /* CLEANING */
    cout<<"...ENDING HAMMURABI..."<<endl;
    delete par; delete grid_breg; delete grid_brnd;
    delete grid_fe;delete grid_fernd; delete grid_cre; delete grid_int;
    delete fe;delete fernd; delete breg; delete brnd;
    delete cre; delete intobj; delete doc;
    
    toolkit::timestamp();
    return EXIT_SUCCESS;
}

// END
