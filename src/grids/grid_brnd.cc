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
#include <cassert>
using namespace tinyxml2;
using namespace std;

Grid_brnd::Grid_brnd(string file_name){
    unique_ptr<XMLDocument> doc = unique_ptr<XMLDocument> (new XMLDocument());
    doc->LoadFile(file_name.c_str());
    XMLElement *ptr {doc->FirstChildElement("root")->FirstChildElement("MagneticField")->FirstChildElement("Random")};
    // sometimes users don't want to write out random field
    // but generation of random field needs grid
    build_permission = ptr->BoolAttribute("cue");
    ptr = doc->FirstChildElement("root")->FirstChildElement("Fieldout")->FirstChildElement("brnd_grid");
    read_permission = ptr->BoolAttribute("read");
    write_permission = ptr->BoolAttribute("write");
    if(build_permission or read_permission){
        build_grid(doc.get());
    }
    if(read_permission or write_permission){
        filename = ptr->Attribute("filename");
    }
}

void Grid_brnd::build_grid(XMLDocument *doc){
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
    // real 3D random b field
    fftw_b_x = unique_ptr<double[]> (new double[full_size]);
    fftw_b_y = unique_ptr<double[]> (new double[full_size]);
    fftw_b_z = unique_ptr<double[]> (new double[full_size]);
    // complex random b field in k-space
    fftw_b_kx = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex)*full_size));
    fftw_b_ky = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex)*full_size));
    fftw_b_kz = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex)*full_size));
    // DFT plans
#ifdef _OPENMP
    fftw_init_threads();
    fftw_plan_with_nthreads(omp_get_max_threads());
#endif
    // backword plan
    fftw_px_bw = fftw_plan_dft_3d(nx,ny,nz,fftw_b_kx,fftw_b_kx,FFTW_BACKWARD,FFTW_ESTIMATE);
    fftw_py_bw = fftw_plan_dft_3d(nx,ny,nz,fftw_b_ky,fftw_b_ky,FFTW_BACKWARD,FFTW_ESTIMATE);
    fftw_pz_bw = fftw_plan_dft_3d(nx,ny,nz,fftw_b_kz,fftw_b_kz,FFTW_BACKWARD,FFTW_ESTIMATE);
    // forward plan
    fftw_px_fw = fftw_plan_dft_3d(nx,ny,nz,fftw_b_kx,fftw_b_kx,FFTW_FORWARD,FFTW_ESTIMATE);
    fftw_py_fw = fftw_plan_dft_3d(nx,ny,nz,fftw_b_ky,fftw_b_ky,FFTW_FORWARD,FFTW_ESTIMATE);
    fftw_pz_fw = fftw_plan_dft_3d(nx,ny,nz,fftw_b_kz,fftw_b_kz,FFTW_FORWARD,FFTW_ESTIMATE);
}

void Grid_brnd::export_grid(void){
    assert(!filename.empty());
    ofstream output(filename.c_str(), std::ios::out|std::ios::binary);
    assert(output.is_open());
    double tmp;
    for(decltype(full_size) i=0;i!=full_size;++i){
        assert(!output.eof());
        tmp = fftw_b_x[i];
        output.write(reinterpret_cast<char*>(&tmp),sizeof(double));
        tmp = fftw_b_y[i];
        output.write(reinterpret_cast<char*>(&tmp),sizeof(double));
        tmp = fftw_b_z[i];
        output.write(reinterpret_cast<char*>(&tmp),sizeof(double));
    }
    output.close();
    exit(0);
}

void Grid_brnd::import_grid(void){
    assert(!filename.empty());
    ifstream input(filename.c_str(), std::ios::in|std::ios::binary);
    assert(input.is_open());
    double tmp;
    for(decltype(full_size) i=0;i!=full_size;++i){
        assert(!input.eof());
        input.read(reinterpret_cast<char *>(&tmp),sizeof(double));
        fftw_b_x[i] = tmp;
        input.read(reinterpret_cast<char *>(&tmp),sizeof(double));
        fftw_b_y[i] = tmp;
        input.read(reinterpret_cast<char *>(&tmp),sizeof(double));
        fftw_b_z[i] = tmp;
    }
#ifndef NDEBUG
    auto eof = input.tellg();
    input.seekg (0, input.end);
#endif
    assert(eof==input.tellg());
    input.close();
}
