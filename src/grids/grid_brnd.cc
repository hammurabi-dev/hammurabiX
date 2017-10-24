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

Grid_brnd::Grid_brnd(string file_name){
    XMLDocument *doc = new XMLDocument();
    doc->LoadFile(file_name.c_str());
    XMLElement *ptr {doc->FirstChildElement("root")->FirstChildElement("Galaxy")->FirstChildElement("MagneticField")->FirstChildElement("Random")};
    // sometimes users don't want to write out random field
    // but generation of random field needs grid
    bool build_permission {ptr->BoolAttribute("cue")};
    ptr = doc->FirstChildElement("root")->FirstChildElement("Interface")->FirstChildElement("brnd_grid");
    read_permission = ptr->BoolAttribute("read");
    write_permission = ptr->BoolAttribute("write");
    if(build_permission or read_permission){
        build_grid(doc);
    }
    if(read_permission or write_permission){
        cout<<"IFNO: GRID_BRND I/O ACTIVE"<<endl;
        filename = ptr->Attribute("filename");
    }
    delete doc;
    doc = nullptr;
}

void Grid_brnd::build_grid(XMLDocument *doc){
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
    // memory check (double complex + double + double)
    const double bytes {full_size*(3.*16.+3.*8.)};
    cout<<"INFO: BRND REQUIRING "<<bytes/1.e9<<" GB MEMORY"<<endl;
    // real 3D random b field
    fftw_b_x = new double[full_size];
    fftw_b_y = new double[full_size];
    fftw_b_z = new double[full_size];
    // complex random b field in k-space
    fftw_b_kx = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex)*full_size));
    fftw_b_ky = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex)*full_size));
    fftw_b_kz = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex)*full_size));
    // DFT plans
    // backword plan
    fftw_px_bw = fftw_plan_dft_3d(nx,ny,nz,fftw_b_kx,fftw_b_kx,FFTW_BACKWARD,FFTW_ESTIMATE);
    fftw_py_bw = fftw_plan_dft_3d(nx,ny,nz,fftw_b_ky,fftw_b_ky,FFTW_BACKWARD,FFTW_ESTIMATE);
    fftw_pz_bw = fftw_plan_dft_3d(nx,ny,nz,fftw_b_kz,fftw_b_kz,FFTW_BACKWARD,FFTW_ESTIMATE);
    // forward plan
    fftw_px_fw = fftw_plan_dft_3d(nx,ny,nz,fftw_b_kx,fftw_b_kx,FFTW_FORWARD,FFTW_ESTIMATE);
    fftw_py_fw = fftw_plan_dft_3d(nx,ny,nz,fftw_b_ky,fftw_b_ky,FFTW_FORWARD,FFTW_ESTIMATE);
    fftw_pz_fw = fftw_plan_dft_3d(nx,ny,nz,fftw_b_kz,fftw_b_kz,FFTW_FORWARD,FFTW_ESTIMATE);
}

void Grid_brnd::clean_grid(void){
    fftw_destroy_plan(fftw_px_bw);
    fftw_destroy_plan(fftw_py_bw);
    fftw_destroy_plan(fftw_pz_bw);
    fftw_destroy_plan(fftw_px_fw);
    fftw_destroy_plan(fftw_py_fw);
    fftw_destroy_plan(fftw_pz_fw);
    fftw_free(fftw_b_kx);
    fftw_free(fftw_b_ky);
    fftw_free(fftw_b_kz);
    delete [] fftw_b_x;
    fftw_b_x = nullptr;
    delete [] fftw_b_y;
    fftw_b_y = nullptr;
    delete [] fftw_b_z;
    fftw_b_z = nullptr;
}

void Grid_brnd::export_grid(void){
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
        tmp = fftw_b_x[i];
        output.write(reinterpret_cast<char*>(&tmp),sizeof(double));
        tmp = fftw_b_y[i];
        output.write(reinterpret_cast<char*>(&tmp),sizeof(double));
        tmp = fftw_b_z[i];
        output.write(reinterpret_cast<char*>(&tmp),sizeof(double));
    }
    output.close();
    // exit program
    clean_grid();
    cout<<"...RANDOM MAGNETIC FIELD EXPORTED AND CLEANED..."<<endl;
    exit(0);
}

void Grid_brnd::import_grid(void){
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
        fftw_b_x[i] = tmp;
        input.read(reinterpret_cast<char *>(&tmp),sizeof(double));
        fftw_b_y[i] = tmp;
        input.read(reinterpret_cast<char *>(&tmp),sizeof(double));
        fftw_b_z[i] = tmp;
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
