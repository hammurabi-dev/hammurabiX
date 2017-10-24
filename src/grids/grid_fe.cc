#include <sstream>
#include <iostream>
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

Grid_fe::Grid_fe(string file_name){
    XMLDocument *doc = new XMLDocument();
    doc->LoadFile(file_name.c_str());
    XMLElement *ptr {doc->FirstChildElement("root")->FirstChildElement("Interface")->FirstChildElement("fe_grid")};
    read_permission = ptr->BoolAttribute("read");
    write_permission = ptr->BoolAttribute("write");
    if(read_permission or write_permission){
        cout<<"IFNO: FE I/O ACTIVE"<<endl;
        filename = ptr->Attribute("filename");
        build_grid(doc);
    }
    delete doc;
    doc = nullptr;
}

void Grid_fe::build_grid(XMLDocument *doc){
    XMLElement *ptr {doc->FirstChildElement("root")->FirstChildElement("Grid")->FirstChildElement("Box")};
    // Cartesian grid
    nx = FetchUnsigned(ptr,"nx");
    ny = FetchUnsigned(ptr,"ny");
    nz = FetchUnsigned(ptr,"nz");
    full_size = nx*ny*nz;
    // box limit for filling field
    x_max = CGS_U_kpc*FetchDouble(ptr,"x_max");
    x_min = CGS_U_kpc*FetchDouble(ptr,"x_min");
    y_max = CGS_U_kpc*FetchDouble(ptr,"y_max");
    y_min = CGS_U_kpc*FetchDouble(ptr,"y_min");
    z_max = CGS_U_kpc*FetchDouble(ptr,"z_max");
    z_min = CGS_U_kpc*FetchDouble(ptr,"z_min");
    // memory check
    const double bytes {full_size*8.};
    cout<<"INFO: FE REQUIRING "<<bytes/1.e9<<" GB MEMORY"<<endl;
    fe = new double[full_size];
}

void Grid_fe::clean_grid(void){
    delete [] fe;
    fe = nullptr;
}

void Grid_fe::export_grid(void){
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
    for(decltype(full_size) i=0;i!=full_size;++i){
        if (output.eof()) {
            cerr<<"ERR:"<<__FILE__
            <<" : in function "<<__func__<<endl
            <<" at line "<<__LINE__<<endl
            <<"UNEXPECTED END AT: "<<i<<endl;
            exit(1);
        }
        tmp = fe[i];
        if(tmp<0) {
            cerr<<"ERR:"<<__FILE__
            <<" : in function "<<__func__<<endl
            <<" at line "<<__LINE__<<endl
            <<"WRONG VALUE"<<endl;
            exit(1);
        }
        output.write(reinterpret_cast<char*>(&tmp),sizeof(double));
    }
    output.close();
    // exit program
    clean_grid();
    cout<<"...FREE ELECTRON FIELD EXPORTED AND CLEANED..."<<endl;
    exit(0);
}

void Grid_fe::import_grid(void){
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
    for(decltype(full_size) i=0;i!=full_size;++i){
        if (input.eof()) {
            cerr<<"ERR:"<<__FILE__
            <<" : in function "<<__func__<<endl
            <<" at line "<<__LINE__<<endl
            <<"UNEXPECTED END AT: "<<i<<endl;
            exit(1);
        }
        input.read(reinterpret_cast<char *>(&tmp),sizeof(double));
        fe[i] = tmp;
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
