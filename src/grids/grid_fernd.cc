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

// turbulent free electron density field
Grid_fernd::Grid_fernd(string file_name){
    unique_ptr<XMLDocument> doc = toolkit::loadxml(file_name);
    XMLElement *ptr {toolkit::tracexml(doc.get(),{"FreeElectron"})};
    build_permission = toolkit::FetchBool(ptr,"cue","Random");
    // sometimes users don't want to write out random field
    // but generation of random field needs grid
    ptr = toolkit::tracexml(doc.get(),{"Fieldout"});
    read_permission = toolkit::FetchBool(ptr,"read","fernd_grid");
    write_permission = toolkit::FetchBool(ptr,"write","fernd_grid");
    if(build_permission or read_permission){
        build_grid(doc.get());
    }
    if(read_permission or write_permission){
        filename = toolkit::FetchString(ptr,"filename","fernd_grid");
    }
}

void Grid_fernd::build_grid(XMLDocument *doc){
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
    assert(!filename.empty());
    ofstream output(filename.c_str(), std::ios::out|std::ios::binary);
    assert(output.is_open());
    double tmp;
    for(decltype(full_size) i=0;i!=full_size;++i){
        assert(!output.eof());
        tmp = fftw_fe[i];
        output.write(reinterpret_cast<char*>(&tmp),sizeof(double));
    }
    output.close();
    exit(0);
}

void Grid_fernd::import_grid(void){
    assert(!filename.empty());
    ifstream input(filename.c_str(), std::ios::in|std::ios::binary);
    assert(input.is_open());
    double tmp;
    for(decltype(full_size) i=0;i!=full_size;++i){
        assert(!input.eof());
        input.read(reinterpret_cast<char *>(&tmp),sizeof(double));
        fftw_fe[i] = tmp;
    }
#ifndef NDEBUG
    auto eof = input.tellg();
    input.seekg (0, input.end);
#endif
    assert(eof==input.tellg());
    input.close();
}

// END
