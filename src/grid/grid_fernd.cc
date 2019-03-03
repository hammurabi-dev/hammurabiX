#include <sstream>
#include <fstream>
#include <memory>
#include <array>
#include <string>
#include <vector>
#include <cassert>

#include <fftw3.h>
#include <omp.h>

#include <grid.h>
#include <cgs_units_file.h>
#include <namespace_toolkit.h>
#include <param.h>

// turbulent free electron density field
Grid_fernd::Grid_fernd (const Param *par){
    if (par->grid_fernd.build_permission or par->grid_fernd.read_permission){
        build_grid (par);
        clean_switch = true;
    }
}

void Grid_fernd::build_grid (const Param *par){
    // real random fe field
    fe = std::make_unique<double[]> (par->grid_fernd.full_size);
    // complex random b field in k-space
    fe_k = fftw_alloc_complex (par->grid_fernd.full_size);
    fe_k[0][0]=0;fe_k[0][1]=0; // 0th term should be zero
    // DFT plan
#ifdef _OPENMP
    fftw_init_threads();
    fftw_plan_with_nthreads (omp_get_max_threads());
#endif
    plan_fe_bw = fftw_plan_dft_3d (par->grid_fernd.nx,
                                   par->grid_fernd.ny,
                                   par->grid_fernd.nz,
                                   fe_k,
                                   fe_k,
                                   FFTW_BACKWARD,
                                   FFTW_ESTIMATE);
}

void Grid_fernd::export_grid (const Param *par){
    assert (!par->grid_fernd.filename.empty());
    std::ofstream output (par->grid_fernd.filename.c_str(),
                          std::ios::out|std::ios::binary);
    assert (output.is_open());
    double tmp;
    for (decltype(par->grid_fernd.full_size) i=0;i!=par->grid_fernd.full_size;++i){
        assert (!output.eof());
        tmp = fe[i];
        output.write (reinterpret_cast<char*>(&tmp),
                      sizeof(double));
    }
    output.close();
}

void Grid_fernd::import_grid (const Param *par){
    assert (!par->grid_fernd.filename.empty());
    std::ifstream input (par->grid_fernd.filename.c_str(),
                         std::ios::in|std::ios::binary);
    assert (input.is_open());
    double tmp;
    for (decltype(par->grid_fernd.full_size) i=0;i!=par->grid_fernd.full_size;++i){
        assert (!input.eof());
        input.read (reinterpret_cast<char *>(&tmp),
                    sizeof(double));
        fe[i] = tmp;
    }
#ifndef NDEBUG
    auto eof = input.tellg();
    input.seekg (0,input.end);
#endif
    assert (eof==input.tellg());
    input.close();
}

// END
