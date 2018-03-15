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
#include <grid.h>
#include <cgs_units_file.h>
#include <namespace_toolkit.h>
#include <ap_err.h>
using namespace tinyxml2;
using namespace std;

Grid_cre::Grid_cre(string file_name){
    unique_ptr<XMLDocument> doc = unique_ptr<XMLDocument> (new XMLDocument());
    doc->LoadFile(file_name.c_str());
    XMLElement *ptr {doc->FirstChildElement("root")->FirstChildElement("Fieldout")->FirstChildElement("cre_grid")};
    read_permission = ptr->BoolAttribute("read");
    write_permission = ptr->BoolAttribute("write");
    // security check
    string type_check {doc->FirstChildElement("root")->FirstChildElement("CRE")->Attribute("type")};
    if(read_permission and type_check=="Analytic"){
        ap_err("analytic CRE does not read");
        exit(1);
    }
    if(read_permission or write_permission){
#ifdef DEBUG
        cout<<"IFNO: CRE I/O ACTIVE"<<endl;
#endif
        filename = ptr->Attribute("filename");
        build_grid(doc.get());
    }
}

void Grid_cre::build_grid(XMLDocument *doc){
    XMLElement *ptr {doc->FirstChildElement("root")->FirstChildElement("CRE")->FirstChildElement("Numeric")};
    E_min = CGS_U_GeV*toolkit::FetchDouble(ptr,"E_min");
    E_max = CGS_U_GeV*toolkit::FetchDouble(ptr,"E_max");
    E_fact = toolkit::FetchDouble(ptr,"E_fact");
    nE = ceil(log(E_max/E_min)/E_fact);
    // spatial 2D
    if(strcmp(ptr->Attribute("type"),"2D")==0){
        nr = toolkit::FetchUnsigned(ptr,"nr");
        nz = toolkit::FetchUnsigned(ptr,"nz");
        nx = 0; ny = 0;
        r_max = CGS_U_kpc*toolkit::FetchDouble(ptr,"r_max");
        z_max = CGS_U_kpc*toolkit::FetchDouble(ptr,"z_max");
        z_min = CGS_U_kpc*toolkit::FetchDouble(ptr,"z_min");
        x_max = 0.; x_min = 0.; y_max = 0.; y_min = 0.;
        cre_size = nE*nr*nz;
    }
    // spatial 3D
    else if(strcmp(ptr->Attribute("type"),"3D")==0){
        nr = 0;
        nz = toolkit::FetchUnsigned(ptr,"nz");
        nx = toolkit::FetchUnsigned(ptr,"nx");
        ny = toolkit::FetchUnsigned(ptr,"ny");
        x_max = CGS_U_kpc*toolkit::FetchDouble(ptr,"x_max");
        x_min = CGS_U_kpc*toolkit::FetchDouble(ptr,"x_min");
        y_max = CGS_U_kpc*toolkit::FetchDouble(ptr,"y_max");
        y_min = CGS_U_kpc*toolkit::FetchDouble(ptr,"y_min");
        z_max = CGS_U_kpc*toolkit::FetchDouble(ptr,"z_max");
        z_min = CGS_U_kpc*toolkit::FetchDouble(ptr,"z_min");
        r_max = 0.;
        cre_size = nE*nx*ny*nz;
    }
    else{
        ap_err("unexpected grid dimension");
        exit(1);
    }
    // memory check
#ifdef DEBUG
    const double bytes {cre_size*8.};
    cout<<"INFO: CRE REQUIRING "<<bytes/1.e9<<" GB MEMORY"<<endl;
#endif
    cre_flux = unique_ptr<double[]> (new double[cre_size]);
}

void Grid_cre::export_grid(void){
    if(filename.empty()){
        ap_err("invalid filename");
        exit(1);
    }
    ofstream output(filename.c_str(), std::ios::out|std::ios::binary);
    if (!output.is_open()){
        ap_err("open file fail");
        exit(1);
    }
    double tmp;
    for(decltype(cre_size) i=0;i!=cre_size;++i){
        if (output.eof()) {
            ap_err("unexpected");
            exit(1);
        }
        tmp = cre_flux[i];
        output.write(reinterpret_cast<char*>(&tmp),sizeof(double));
    }
    output.close();
    // exit program
#ifdef DEBUG
    if(nx!=0){
        cout<<"...COSMIC-RAY ELECTRON IN 4D GRID..."<<endl;
    }
    else{
        cout<<"...COSMIC-RAY ELECTRON IN 3D GRID..."<<endl;
    }
    cout<<"...FIELD EXPORTED AND CLEANED..."<<endl;
#endif
    exit(0);
}

void Grid_cre::import_grid(void){
    if(filename.empty()){
        ap_err("invalid filename");
        exit(1);
    }
    ifstream input(filename.c_str(), std::ios::in|std::ios::binary);
    if (!input.is_open()){
        ap_err("open file fail");
        exit(1);
    }
    double tmp;
    for(decltype(cre_size) i=0;i!=cre_size;++i){
        if (input.eof()) {
            ap_err("unexpected");
            exit(1);
        }
        input.read(reinterpret_cast<char *>(&tmp),sizeof(double));
        cre_flux[i] = tmp;
    }
    auto eof = input.tellg();
    input.seekg (0, input.end);
    if (eof != input.tellg()){
        ap_err("incorrect length");
        exit(1);
    }
    input.close();
}
