#include <sstream>
#include <iostream>
#include <memory>
#include <fftw3.h>
#include <array>
#include <string>
#include <vector>
#include <omp.h>
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

// turbulent free electron density field
Grid_fernd::Grid_fernd(string file_name){
    unique_ptr<XMLDocument> doc = unique_ptr<XMLDocument> (new XMLDocument());
    doc->LoadFile(file_name.c_str());
    XMLElement *ptr {doc->FirstChildElement("root")->FirstChildElement("FreeElectron")->FirstChildElement("Random")};
    build_permission = ptr->BoolAttribute("cue");
    // sometimes users don't want to write out random field
    // but generation of random field needs grid
    ptr = doc->FirstChildElement("root")->FirstChildElement("Fieldout")->FirstChildElement("fernd_grid");
    read_permission = ptr->BoolAttribute("read");
    write_permission = ptr->BoolAttribute("write");
    if(build_permission or read_permission){
        build_grid(doc.get());
    }
    if(read_permission or write_permission){
#ifdef DEBUG
        cout<<"IFNO: FE_RND I/O ACTIVE"<<endl;
#endif
        filename = ptr->Attribute("filename");
    }
}

void Grid_fernd::build_grid(XMLDocument *doc){
    XMLElement *ptr {doc->FirstChildElement("root")->FirstChildElement("Grid")->FirstChildElement("Box")};
    // Cartesian grid
    nx = toolkit::FetchUnsigned(ptr,"nx");
    ny = toolkit::FetchUnsigned(ptr,"ny");
    nz = toolkit::FetchUnsigned(ptr,"nz");
    full_size = nx*ny*nz;
    // box limit for filling field
    x_max = CGS_U_kpc*toolkit::FetchDouble(ptr,"x_max");
    x_min = CGS_U_kpc*toolkit::FetchDouble(ptr,"x_min");
    y_max = CGS_U_kpc*toolkit::FetchDouble(ptr,"y_max");
    y_min = CGS_U_kpc*toolkit::FetchDouble(ptr,"y_min");
    z_max = CGS_U_kpc*toolkit::FetchDouble(ptr,"z_max");
    z_min = CGS_U_kpc*toolkit::FetchDouble(ptr,"z_min");
#ifdef DEBUG
    // memory check (double complex + double + double)
    const double bytes {full_size*(16.+ 8.)};
    cout<<"INFO: FERND REQUIRING "<<bytes/1.e9<<" GB MEMORY"<<endl;
#endif
    // real random fe field
    fftw_fe = unique_ptr<double[]> (new double[full_size]);
    // complex random b field in k-space
    fftw_fe_k = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex)*full_size));
    // DFT plan
#ifdef _OPENMP
    fftw_init_threads();
    fftw_plan_with_nthreads(omp_get_max_threads());
#endif
    fftw_p_bw = fftw_plan_dft_3d(nx,ny,nz,fftw_fe_k,fftw_fe_k,FFTW_BACKWARD,FFTW_ESTIMATE);
}

void Grid_fernd::export_grid(void){
    if(filename.empty()){
        ap_err("invalid filename");
        exit(1);
    }
    ofstream output(filename.c_str(), std::ios::out|std::ios::binary);
    if (!output.is_open()){
        ap_err("file open fail");
        exit(1);
    }
    double tmp;
    for(decltype(full_size) i=0;i!=full_size;++i){
        if (output.eof()) {
            ap_err("unexpected");
            exit(1);
        }
        tmp = fftw_fe[i];
        output.write(reinterpret_cast<char*>(&tmp),sizeof(double));
    }
    output.close();
    // exit program
#ifdef DEBUG
    cout<<"...RANDOM FREE ELECTRON FIELD EXPORTED AND CLEANED..."<<endl;
#endif
    exit(0);
}

void Grid_fernd::import_grid(void){
    if(filename.empty()){
        ap_err("invalid filename");
        exit(1);
    }
    ifstream input(filename.c_str(), std::ios::in|std::ios::binary);
    if (!input.is_open()){
        ap_err("file open fail");
        exit(1);
    }
    double tmp;
    for(decltype(full_size) i=0;i!=full_size;++i){
        if (input.eof()) {
            ap_err("unexpected");
            exit(1);
        }
        input.read(reinterpret_cast<char *>(&tmp),sizeof(double));
        fftw_fe[i] = tmp;
    }
    auto eof = input.tellg();
    input.seekg (0, input.end);
    if (eof != input.tellg()){
        ap_err("incorrect length");
        exit(1);
    }
    input.close();
}

// END
