#include <sstream>
#include <array>
#include <string>
#include <vector>
#include <fstream>
#include <memory>
#include <omp.h>
#include <cassert>

#include <fftw3.h>
#include <tinyxml2.h>

#include <grid.h>
#include <cgs_units_file.h>
#include <namespace_toolkit.h>

Grid_brnd::Grid_brnd (const std::string &file_name){
    std::unique_ptr<tinyxml2::XMLDocument> doc = toolkit::loadxml (file_name);
    tinyxml2::XMLElement *ptr {toolkit::tracexml (doc.get(),{"MagneticField"})};
    // sometimes users don't want to write out random field
    // but generation of random field needs grid
    build_permission = toolkit::fetchbool (ptr,"cue","Random");
    ptr = toolkit::tracexml (doc.get(),{"Fieldout"});
    read_permission = toolkit::fetchbool (ptr,"read","brnd_grid");
    write_permission = toolkit::fetchbool (ptr,"write","brnd_grid");
    if (build_permission or read_permission){
        build_grid(doc.get());
    }
    if (read_permission or write_permission){
        filename = toolkit::fetchstring (ptr,"filename","brnd_grid");
    }
}

void Grid_brnd::build_grid (tinyxml2::XMLDocument *doc){
    tinyxml2::XMLElement *ptr {toolkit::tracexml(doc,{"Grid","Box_GMF"})};
    // Cartesian grid
    nx = toolkit::fetchunsigned (ptr,"value","nx");
    ny = toolkit::fetchunsigned (ptr,"value","ny");
    nz = toolkit::fetchunsigned (ptr,"value","nz");
    full_size = nx*ny*nz;
    // box limit for filling field
    x_max = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","x_max");
    x_min = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","x_min");
    y_max = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","y_max");
    y_min = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","y_min");
    z_max = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","z_max");
    z_min = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","z_min");
    // real 3D random b field
    bx = std::make_unique<double[]> (full_size);
    by = std::make_unique<double[]> (full_size);
    bz = std::make_unique<double[]> (full_size);
    // complex a field in k-space
    c0 = fftw_alloc_complex (full_size);
    c0[0][0]=0;c0[0][1]=0; // 0th term should be zero
    c1 = fftw_alloc_complex (full_size);
    c1[0][0]=0;c1[0][1]=0; // 0th term should be zero
    // DFT plans
#ifdef _OPENMP
    fftw_init_threads();
    fftw_plan_with_nthreads (omp_get_max_threads());
#endif
    // backword plan
    plan_c0_bw = fftw_plan_dft_3d (nx,ny,nz,c0,c0,FFTW_BACKWARD,FFTW_ESTIMATE);
    plan_c1_bw = fftw_plan_dft_3d (nx,ny,nz,c1,c1,FFTW_BACKWARD,FFTW_ESTIMATE);
    // forward plan
    plan_c0_fw = fftw_plan_dft_3d (nx,ny,nz,c0,c0,FFTW_FORWARD,FFTW_ESTIMATE);
    plan_c1_fw = fftw_plan_dft_3d (nx,ny,nz,c1,c1,FFTW_FORWARD,FFTW_ESTIMATE);
}

void Grid_brnd::export_grid (){
    assert (!filename.empty());
    std::ofstream output (filename.c_str(),
                          std::ios::out|std::ios::binary);
    assert (output.is_open());
    double tmp;
    for (decltype(full_size) i=0;i!=full_size;++i){
        assert (!output.eof());
        tmp = bx[i];
        output.write (reinterpret_cast<char*>(&tmp),
                      sizeof(double));
        tmp = by[i];
        output.write (reinterpret_cast<char*>(&tmp),
                      sizeof(double));
        tmp = bz[i];
        output.write (reinterpret_cast<char*>(&tmp),
                      sizeof(double));
    }
    output.close();
    exit (0);
}

void Grid_brnd::import_grid (){
    assert (!filename.empty());
    std::ifstream input (filename.c_str(),
                         std::ios::in|std::ios::binary);
    assert (input.is_open());
    double tmp;
    for (decltype(full_size) i=0;i!=full_size;++i){
        assert (!input.eof());
        input.read (reinterpret_cast<char *>(&tmp),
                    sizeof(double));
        bx[i] = tmp;
        input.read (reinterpret_cast<char *>(&tmp),
                    sizeof(double));
        by[i] = tmp;
        input.read (reinterpret_cast<char *>(&tmp),
                    sizeof(double));
        bz[i] = tmp;
    }
#ifndef NDEBUG
    auto eof = input.tellg();
    input.seekg (0,input.end);
#endif
    assert (eof==input.tellg());
    input.close();
}
