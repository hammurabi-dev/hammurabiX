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
        cerr<<"ERR:"<<__FILE__
        <<" : in function "<<__func__<<endl
        <<" at line "<<__LINE__<<endl
        <<"ANALYTIC CRE DO NOT READ"<<endl;
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
    E_min = CGS_U_GeV*FetchDouble(ptr,"E_min");
    E_max = CGS_U_GeV*FetchDouble(ptr,"E_max");
    E_fact = FetchDouble(ptr,"E_fact");
    nE = ceil(log(E_max/E_min)/E_fact);
    // spatial 2D
    if(strcmp(ptr->Attribute("type"),"2D")==0){
        nr = FetchUnsigned(ptr,"nr");
        nz = FetchUnsigned(ptr,"nz");
        nx = 0; ny = 0;
        r_max = CGS_U_kpc*FetchDouble(ptr,"r_max");
        z_max = CGS_U_kpc*FetchDouble(ptr,"z_max");
        z_min = CGS_U_kpc*FetchDouble(ptr,"z_min");
        x_max = 0.; x_min = 0.; y_max = 0.; y_min = 0.;
        cre_size = nE*nr*nz;
    }
    // spatial 3D
    else if(strcmp(ptr->Attribute("type"),"3D")==0){
        nr = 0;
        nz = FetchUnsigned(ptr,"nz");
        nx = FetchUnsigned(ptr,"nx");
        ny = FetchUnsigned(ptr,"ny");
        x_max = CGS_U_kpc*FetchDouble(ptr,"x_max");
        x_min = CGS_U_kpc*FetchDouble(ptr,"x_min");
        y_max = CGS_U_kpc*FetchDouble(ptr,"y_max");
        y_min = CGS_U_kpc*FetchDouble(ptr,"y_min");
        z_max = CGS_U_kpc*FetchDouble(ptr,"z_max");
        z_min = CGS_U_kpc*FetchDouble(ptr,"z_min");
        r_max = 0.;
        cre_size = nE*nx*ny*nz;
    }
    else{
        cerr<<"ERR:"<<__FILE__
        <<" : in function "<<__func__<<endl
        <<" at line "<<__LINE__<<endl
        <<"UNEXPECTED GRID DIMENSION"<<endl;
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
        cerr<<"ERR:"<<__FILE__
        <<" : in function "<<__func__<<endl
        <<" at line "<<__LINE__<<endl
        <<"NONEXIST FILE"<<endl;
        exit(1);
    }
    ofstream output(filename.c_str(), std::ios::out|std::ios::binary);
    if (!output.is_open()){
        cerr<<"ERR:"<<__FILE__
        <<" : in function "<<__func__<<endl
        <<" at line "<<__LINE__<<endl
        <<"COULD NOT OPEN: "<<filename<<endl;
        exit(1);
    }
    double tmp;
    for(decltype(cre_size) i=0;i!=cre_size;++i){
        if (output.eof()) {
            cerr<<"ERR:"<<__FILE__
            <<" : in function "<<__func__<<endl
            <<" at line "<<__LINE__<<endl
            <<"UNEXPECTED END AT: "<<i<<endl;
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
        cerr<<"ERR:"<<__FILE__
        <<" : in function "<<__func__<<endl
        <<" at line "<<__LINE__<<endl
        <<"NONEXIST FILE"<<endl;
        exit(1);
    }
    ifstream input(filename.c_str(), std::ios::in|std::ios::binary);
    if (!input.is_open()){
        cerr<<"ERR:"<<__FILE__
        <<" : in function "<<__func__<<endl
        <<" at line "<<__LINE__<<endl
        <<"COULD NOT OPEN: "<<filename<<endl;
        exit(1);
    }
    double tmp;
    for(decltype(cre_size) i=0;i!=cre_size;++i){
        if (input.eof()) {
            cerr<<"ERR:"<<__FILE__
            <<" : in function "<<__func__<<endl
            <<" at line "<<__LINE__<<endl
            <<"UNEXPECTED END AT: "<<i<<endl;
            exit(1);
        }
        input.read(reinterpret_cast<char *>(&tmp),sizeof(double));
        cre_flux[i] = tmp;
    }
    auto eof = input.tellg();
    input.seekg (0, input.end);
    if (eof != input.tellg()){
        cerr<<"ERR:"<<__FILE__
        <<" : in function "<<__func__<<endl
        <<" at line "<<__LINE__<<endl
        <<"INCORRECT LENGTH"<<endl;
        exit(1);
    }
    input.close();
}
