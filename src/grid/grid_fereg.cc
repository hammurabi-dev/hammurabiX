#include <sstream>
#include <fstream>
#include <memory>
#include <array>
#include <string>
#include <vector>
#include <cassert>

#include <grid.h>
#include <cgs_units_file.h>
#include <namespace_toolkit.h>
#include <param.h>

Grid_fereg::Grid_fereg (const Param *par){
    if (par->grid_fereg.read_permission or par->grid_fereg.write_permission){
        build_grid (par);
    }
}

void Grid_fereg::build_grid (const Param *par){
    fe = std::make_unique<double[]> (par->grid_fereg.full_size);
}

void Grid_fereg::export_grid (const Param *par){
    assert (!par->grid_fereg.filename.empty());
    std::ofstream output (par->grid_fereg.filename.c_str(),
                          std::ios::out|std::ios::binary);
    assert (output.is_open());
    double tmp;
    for (decltype(par->grid_fereg.full_size) i=0;i!=par->grid_fereg.full_size;++i){
        assert (!output.eof());
        tmp = fe[i];
        assert (tmp>=0);
        output.write (reinterpret_cast<char*>(&tmp),
                      sizeof(double));
    }
    output.close();
}

void Grid_fereg::import_grid (const Param *par){
    assert (!par->grid_fereg.filename.empty());
    std::ifstream input (par->grid_fereg.filename.c_str(),
                         std::ios::in|std::ios::binary);
    assert (input.is_open());
    double tmp;
    for (decltype(par->grid_fereg.full_size) i=0;i!=par->grid_fereg.full_size;++i){
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
