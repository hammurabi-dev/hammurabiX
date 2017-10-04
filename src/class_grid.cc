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
#include "class_grid.h"
#include "cgs_units_file.h"
#include "namespace_toolkit.h"

using namespace tinyxml2;
using namespace std;

void Grid::build_grid(XMLDocument *){
    cerr<<"ERR:"<<__FILE__
    <<" : in function "<<__func__<<endl
    <<" at line "<<__LINE__<<endl
    <<"DYNAMIC BINDING FAILURE"<<endl;
    exit(1);
}

void Grid::clean_grid(void){
    cerr<<"ERR:"<<__FILE__
    <<" : in function "<<__func__<<endl
    <<" at line "<<__LINE__<<endl
    <<"DYNAMIC BINDING FAILURE"<<endl;
    exit(1);
}

void Grid::export_grid(void){
    cerr<<"ERR:"<<__FILE__
    <<" : in function "<<__func__<<endl
    <<" at line "<<__LINE__<<endl
    <<"DYNAMIC BINDING FAILURE"<<endl;
    exit(1);
}

void Grid::import_grid(void){
    cerr<<"ERR:"<<__FILE__
    <<" : in function "<<__func__<<endl
    <<" at line "<<__LINE__<<endl
    <<"DYNAMIC BINDING FAILURE"<<endl;
    exit(1);
}

// auxiliary functions
std::string Grid::FetchString(XMLElement* el, string obj){
    return el->FirstChildElement(obj.c_str())->Attribute("value");
}

int Grid::FetchInt(XMLElement* el, string obj){
    return el->FirstChildElement(obj.c_str())->IntAttribute("value");
}

unsigned int Grid::FetchUnsigned(XMLElement* el, string obj){
    return el->FirstChildElement(obj.c_str())->UnsignedAttribute("value");
}

bool Grid::FetchBool(XMLElement* el, string obj){
    return el->FirstChildElement(obj.c_str())->BoolAttribute("value");
}

double Grid::FetchDouble(XMLElement* el, string obj){
    return el->FirstChildElement(obj.c_str())->DoubleAttribute("value");
}

Grid_breg::Grid_breg(string file_name){
    XMLDocument *doc = new XMLDocument();
    doc->LoadFile(file_name.c_str());
    XMLElement *ptr {doc->FirstChildElement("root")->FirstChildElement("Interface")->FirstChildElement("breg_grid")};
    read_permission = ptr->BoolAttribute("read");
    write_permission = ptr->BoolAttribute("write");
    // build up grid when have read or write permission
    if(read_permission or write_permission){
        cout<<"IFNO: GRID_BREG I/O ACTIVE"<<endl;
        filename = ptr->Attribute("filename");
        build_grid(doc);
    }
    
    delete doc;
}

void Grid_breg::build_grid(XMLDocument *doc){
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
    const double bytes {full_size*(3.*8.)};
    cout<<"INFO: BREG REQUIRING "<<bytes/1.e9<<" GB MEMORY"<<endl;
    // real 3D regular b field
    reg_b_x = new double[full_size];
    reg_b_y = new double[full_size];
    reg_b_z = new double[full_size];
}

void Grid_breg::clean_grid(void){
    delete [] reg_b_x;
    delete [] reg_b_y;
    delete [] reg_b_z;
}

void Grid_breg::export_grid(void){
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
        tmp = reg_b_x[i];
        output.write(reinterpret_cast<char*>(&tmp),sizeof(double));
        tmp = reg_b_y[i];
        output.write(reinterpret_cast<char*>(&tmp),sizeof(double));
        tmp = reg_b_z[i];
        output.write(reinterpret_cast<char*>(&tmp),sizeof(double));
    }
    output.close();
    // exit program
    clean_grid();
    cout<<"...REGULAR MAGNETIC FIELD EXPORTED AND CLEANED..."<<endl;
    exit(0);
}

void Grid_breg::import_grid(void){
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
        reg_b_x[i] = tmp;
        input.read(reinterpret_cast<char *>(&tmp),sizeof(double));
        reg_b_y[i] = tmp;
        input.read(reinterpret_cast<char *>(&tmp),sizeof(double));
        reg_b_z[i] = tmp;
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
    delete [] fftw_b_y;
    delete [] fftw_b_z;
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

// turbulent free electron density field
Grid_fernd::Grid_fernd(string file_name){
    XMLDocument *doc = new XMLDocument();
    doc->LoadFile(file_name.c_str());
    XMLElement *ptr {doc->FirstChildElement("root")->FirstChildElement("Galaxy")->FirstChildElement("FreeElectron")->FirstChildElement("Random")};
    bool build_permission {ptr->BoolAttribute("cue")};
    // sometimes users don't want to write out random field
    // but generation of random field needs grid
    ptr = doc->FirstChildElement("root")->FirstChildElement("Interface")->FirstChildElement("fernd_grid");
    read_permission = ptr->BoolAttribute("read");
    write_permission = ptr->BoolAttribute("write");
    if(build_permission or read_permission){
        build_grid(doc);
    }
    if(read_permission or write_permission){
        cout<<"IFNO: FE_RND I/O ACTIVE"<<endl;
        filename = ptr->Attribute("filename");
    }
    delete doc;
}

void Grid_fernd::build_grid(XMLDocument *doc){
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
    const double bytes {full_size*(16.+ 8.)};
    cout<<"INFO: FERND REQUIRING "<<bytes/1.e9<<" GB MEMORY"<<endl;
    // real random fe field
    fftw_fe = new double[full_size];
    // complex random b field in k-space
    fftw_fe_k = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex)*full_size));
    // DFT plan
    fftw_p = fftw_plan_dft_3d(nx,ny,nz,fftw_fe_k,fftw_fe_k,FFTW_BACKWARD,FFTW_MEASURE);
}

void Grid_fernd::clean_grid(void){
    fftw_destroy_plan(fftw_p);
    fftw_free(fftw_fe_k);
    delete [] fftw_fe;
}

void Grid_fernd::export_grid(void){
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
        tmp = fftw_fe[i];
        output.write(reinterpret_cast<char*>(&tmp),sizeof(double));
    }
    output.close();
    // exit program
    clean_grid();
    cout<<"...RANDOM FREE ELECTRON FIELD EXPORTED AND CLEANED..."<<endl;
    exit(0);
}

void Grid_fernd::import_grid(void){
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
        fftw_fe[i] = tmp;
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


Grid_cre::Grid_cre(string file_name){
    XMLDocument *doc = new XMLDocument();
    doc->LoadFile(file_name.c_str());
    XMLElement *ptr {doc->FirstChildElement("root")->FirstChildElement("Interface")->FirstChildElement("cre_grid")};
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
        cout<<"IFNO: CRE I/O ACTIVE"<<endl;
        filename = ptr->Attribute("filename");
        build_grid(doc);
    }
    delete doc;
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
    const double bytes {cre_size*8.};
    cout<<"INFO: CRE REQUIRING "<<bytes/1.e9<<" GB MEMORY"<<endl;
    cre_flux = new double[cre_size];
}

void Grid_cre::clean_grid(void){
    delete [] cre_flux;
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
    clean_grid();
    if(nx!=0){
        cout<<"...COSMIC-RAY ELECTRON IN 4D GRID..."<<endl;
    }
    else{
        cout<<"...COSMIC-RAY ELECTRON IN 3D GRID..."<<endl;
    }
    cout<<"...FIELD EXPORTED AND CLEANED..."<<endl;
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

/* line of sight integrator */
Grid_int::Grid_int(string file_name){
    XMLDocument *doc = new XMLDocument();
    doc->LoadFile(file_name.c_str());
    build_grid(doc);
    
    delete doc;
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
    bin_num = floor((ec_r_max/(FetchDouble(ptr,"ec_r_res")*CGS_U_kpc))/pow(2.,total_shell-1));
    lat_lim = FetchDouble(ptr,"lat_lim")*CGS_U_pi/180.;
    // output file name
    
    ptr = doc->FirstChildElement("root")->FirstChildElement("Output");
    
    do_dm = ptr->FirstChildElement("DM")->BoolAttribute("cue");
    do_fd = ptr->FirstChildElement("Faraday")->BoolAttribute("cue");
    do_sync = ptr->FirstChildElement("Sync")->BoolAttribute("cue");
    
    if (not (do_dm or do_fd or do_sync)) {
        cerr<<"PEACE: NO OUTPUT REQUIRED"<<endl;
        exit(1);
    }
    sim_dm_name = ptr->FirstChildElement("DM")->Attribute("filename");
    sim_sync_name = ptr->FirstChildElement("Sync")->Attribute("filename");
    sim_fd_name = ptr->FirstChildElement("Faraday")->Attribute("filename");
    
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

//END
