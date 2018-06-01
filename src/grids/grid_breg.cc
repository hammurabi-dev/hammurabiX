#include <sstream>
#include <iostream>
#include <memory>
#include <fftw3.h>
#include <array>
#include <string>
#include <vector>
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

Grid_breg::Grid_breg(const std::string &file_name){
    std::unique_ptr<XMLDocument> doc = toolkit::loadxml(file_name);
    XMLElement *ptr {toolkit::tracexml(doc.get(),{"MagneticField"})};
    build_permission = toolkit::FetchBool(ptr,"cue","Regular");
    ptr = toolkit::tracexml(doc.get(),{"Fieldout"});
    read_permission = toolkit::FetchBool(ptr,"read","breg_grid");
    write_permission = toolkit::FetchBool(ptr,"write","breg_grid");
    // build up grid when have read or write permission
    if(read_permission or write_permission){
        filename = toolkit::FetchString(ptr,"filename","breg_grid");
        build_grid(doc.get());
    }
}

void Grid_breg::build_grid(XMLDocument *doc){
    XMLElement *ptr {toolkit::tracexml(doc,{"Grid","Box_GMF"})};
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
    // real 3D regular b field
    bx = std::make_unique<double[]>(full_size);
    by = std::make_unique<double[]>(full_size);
    bz = std::make_unique<double[]>(full_size);
}

void Grid_breg::export_grid(void){
    assert(!filename.empty());
    std::ofstream output(filename.c_str(),std::ios::out|std::ios::binary);
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

void Grid_breg::import_grid(void){
    assert(!filename.empty());
    std::ifstream input(filename.c_str(),std::ios::in|std::ios::binary);
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
    input.seekg(0,input.end);
#endif
    assert(eof==input.tellg());
    input.close();
}
