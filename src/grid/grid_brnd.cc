#include <sstream>
#include <array>
#include <string>
#include <vector>
#include <fstream>
#include <memory>
#include <omp.h>
#include <cassert>

#include <fftw3.h>

#include <grid.h>
#include <cgs_units_file.h>
#include <namespace_toolkit.h>
#include <param.h>

Grid_brnd::Grid_brnd (const Param *par){
    if (par->grid_brnd.build_permission or par->grid_brnd.read_permission){
        build_grid (par);
        clean_switch = true;
    }
}

void Grid_brnd::build_grid (const Param *par){
    // real 3D random b field
    bx = std::make_unique<double[]> (par->grid_brnd.full_size);
    by = std::make_unique<double[]> (par->grid_brnd.full_size);
    bz = std::make_unique<double[]> (par->grid_brnd.full_size);
    // complex a field in k-space
    c0 = fftw_alloc_complex (par->grid_brnd.full_size);
    c0[0][0]=0; c0[0][1]=0; // 0th term should be zero
    c1 = fftw_alloc_complex (par->grid_brnd.full_size);
    c1[0][0]=0; c1[0][1]=0; // 0th term should be zero
    // DFT plans
#ifdef _OPENMP
    fftw_init_threads();
    fftw_plan_with_nthreads (omp_get_max_threads());
#endif
    // backword plan
    plan_c0_bw = fftw_plan_dft_3d (par->grid_brnd.nx,
                                   par->grid_brnd.ny,
                                   par->grid_brnd.nz,
                                   c0,
                                   c0,
                                   FFTW_BACKWARD,
                                   FFTW_ESTIMATE);
    plan_c1_bw = fftw_plan_dft_3d (par->grid_brnd.nx,
                                   par->grid_brnd.ny,
                                   par->grid_brnd.nz,
                                   c1,
                                   c1,
                                   FFTW_BACKWARD,
                                   FFTW_ESTIMATE);
    // forward plan
    plan_c0_fw = fftw_plan_dft_3d (par->grid_brnd.nx,
                                   par->grid_brnd.ny,
                                   par->grid_brnd.nz,
                                   c0,
                                   c0,
                                   FFTW_FORWARD,
                                   FFTW_ESTIMATE);
    plan_c1_fw = fftw_plan_dft_3d (par->grid_brnd.nx,
                                   par->grid_brnd.ny,
                                   par->grid_brnd.nz,
                                   c1,
                                   c1,
                                   FFTW_FORWARD,
                                   FFTW_ESTIMATE);
}

void Grid_brnd::export_grid (const Param *par){
    assert (!par->grid_brnd.filename.empty());
    std::ofstream output (par->grid_brnd.filename.c_str(),
                          std::ios::out|std::ios::binary);
    assert (output.is_open());
    double tmp;
    for (decltype(par->grid_brnd.full_size) i=0;i!=par->grid_brnd.full_size;++i){
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
}

void Grid_brnd::import_grid (const Param *par){
    assert (!par->grid_brnd.filename.empty());
    std::ifstream input (par->grid_brnd.filename.c_str(),
                         std::ios::in|std::ios::binary);
    assert (input.is_open());
    double tmp;
    for (decltype(par->grid_brnd.full_size) i=0;i!=par->grid_brnd.full_size;++i){
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
