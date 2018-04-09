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
    unique_ptr<XMLDocument> doc = toolkit::loadxml(file_name);
    XMLElement *ptr {toolkit::tracexml(doc.get(),{"MagneticField"})};
    // sometimes users don't want to write out random field
    // but generation of random field needs grid
    build_permission = toolkit::FetchBool(ptr,"cue","Random");
    ptr = toolkit::tracexml(doc.get(),{"Fieldout"});
    read_permission = toolkit::FetchBool(ptr,"read","brnd_grid");
    write_permission = toolkit::FetchBool(ptr,"write","brnd_grid");
    if(build_permission or read_permission){
        build_grid(doc.get());
    }
    if(read_permission or write_permission){
        filename = toolkit::FetchString(ptr,"filename","brnd_grid");
    }
}

void Grid_brnd::build_grid(XMLDocument *doc){
    XMLElement *ptr {toolkit::tracexml(doc,{"Grid","Box"})};
    // Cartesian grid
    nx = toolkit::FetchUnsigned(ptr,"value","nx");
    ny = toolkit::FetchUnsigned(ptr,"value","ny");
    nz = toolkit::FetchUnsigned(ptr,"value","nz");
    full_size = nx*ny*nz;
    // box limit for filling field
    x_max = CGS_U_kpc*toolkit::FetchDouble(ptr,"value","x_max");
    x_min = CGS_U_kpc*toolkit::FetchDouble(ptr,"value","x_min");
    y_max = CGS_U_kpc*toolkit::FetchDouble(ptr,"value","y_max");
    y_min = CGS_U_kpc*toolkit::FetchDouble(ptr,"value","y_min");
    z_max = CGS_U_kpc*toolkit::FetchDouble(ptr,"value","z_max");
    z_min = CGS_U_kpc*toolkit::FetchDouble(ptr,"value","z_min");
    // real 3D random b field
    bx = unique_ptr<double[]> (new double[full_size]);
    by = unique_ptr<double[]> (new double[full_size]);
    bz = unique_ptr<double[]> (new double[full_size]);
    // complex a field in k-space
    c0 = fftw_alloc_complex(full_size);
    c0[0][0]=0;c0[0][1]=0; // 0th term should be zero
    c1 = fftw_alloc_complex(full_size);
    c1[0][0]=0;c1[0][1]=0; // 0th term should be zero
    // DFT plans
#ifdef _OPENMP
    fftw_init_threads();
    fftw_plan_with_nthreads(omp_get_max_threads());
#endif
    // backword plan
    plan_c0_bw = fftw_plan_dft_3d(nx,ny,nz,c0,c0,FFTW_BACKWARD,FFTW_ESTIMATE);
    plan_c1_bw = fftw_plan_dft_3d(nx,ny,nz,c1,c1,FFTW_BACKWARD,FFTW_ESTIMATE);
    // forward plan
    plan_c0_fw = fftw_plan_dft_3d(nx,ny,nz,c0,c0,FFTW_FORWARD,FFTW_ESTIMATE);
    plan_c1_fw = fftw_plan_dft_3d(nx,ny,nz,c1,c1,FFTW_FORWARD,FFTW_ESTIMATE);
}

void Grid_brnd::export_grid(void){
    assert(!filename.empty());
    ofstream output(filename.c_str(), std::ios::out|std::ios::binary);
    assert(output.is_open());
    double tmp;
    for(decltype(full_size) i=0;i!=full_size;++i){
        assert(!output.eof());
        tmp = bx[i];
        output.write(reinterpret_cast<char*>(&tmp),sizeof(double));
        tmp = by[i];
        output.write(reinterpret_cast<char*>(&tmp),sizeof(double));
        tmp = bz[i];
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
        bx[i] = tmp;
        input.read(reinterpret_cast<char *>(&tmp),sizeof(double));
        by[i] = tmp;
        input.read(reinterpret_cast<char *>(&tmp),sizeof(double));
        bz[i] = tmp;
    }
#ifndef NDEBUG
    auto eof = input.tellg();
    input.seekg (0, input.end);
#endif
    assert(eof==input.tellg());
    input.close();
}
