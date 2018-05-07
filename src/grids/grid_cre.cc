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
using namespace std;

Grid_cre::Grid_cre(string file_name){
    unique_ptr<XMLDocument> doc = toolkit::loadxml(file_name);
    XMLElement *ptr {toolkit::tracexml(doc.get(),{"Fieldout"})};
    read_permission = toolkit::FetchBool(ptr,"read","cre_grid");
    write_permission = toolkit::FetchBool(ptr,"write","cre_grid");
    // build up grid when have read or write permission
    if(read_permission or write_permission){
        filename = toolkit::FetchString(ptr,"filename","cre_grid");
        build_grid(doc.get());
    }
}

void Grid_cre::build_grid(XMLDocument *doc){
    XMLElement *ptr {toolkit::tracexml(doc,{"CRE","Numeric"})};
    E_min = CGS_U_GeV*toolkit::FetchDouble(ptr,"value","E_min");
    E_max = CGS_U_GeV*toolkit::FetchDouble(ptr,"value","E_max");
    nE = toolkit::FetchUnsigned(ptr,"value","nE");
    E_fact = log(E_max/E_min)/(nE-1);
    // spatial 3D
    nz = toolkit::FetchUnsigned(ptr,"value","nz");
    nx = toolkit::FetchUnsigned(ptr,"value","nx");
    ny = toolkit::FetchUnsigned(ptr,"value","ny");
    x_max = CGS_U_kpc*toolkit::FetchDouble(ptr,"value","x_max");
    x_min = CGS_U_kpc*toolkit::FetchDouble(ptr,"value","x_min");
    y_max = CGS_U_kpc*toolkit::FetchDouble(ptr,"value","y_max");
    y_min = CGS_U_kpc*toolkit::FetchDouble(ptr,"value","y_min");
    z_max = CGS_U_kpc*toolkit::FetchDouble(ptr,"value","z_max");
    z_min = CGS_U_kpc*toolkit::FetchDouble(ptr,"value","z_min");
    cre_size = nE*nx*ny*nz;
    cre_flux = unique_ptr<double[]> (new double[cre_size]);
}

void Grid_cre::export_grid(void){
    assert(!filename.empty());
    ofstream output(filename.c_str(), std::ios::out|std::ios::binary);
    assert(output.is_open());
    double tmp;
    for(decltype(cre_size) i=0;i!=cre_size;++i){
        assert(!output.eof());
        tmp = cre_flux[i];
        output.write(reinterpret_cast<char*>(&tmp),sizeof(double));
    }
    output.close();
    exit(0);
}

void Grid_cre::import_grid(void){
    assert(!filename.empty());
    ifstream input(filename.c_str(), std::ios::in|std::ios::binary);
    assert(input.is_open());
    double tmp;
    for(decltype(cre_size) i=0;i!=cre_size;++i){
        assert(!input.eof());
        input.read(reinterpret_cast<char *>(&tmp),sizeof(double));
        cre_flux[i] = tmp;
    }
#ifndef NDEBUG
    auto eof = input.tellg();
    input.seekg (0, input.end);
#endif
    assert(eof==input.tellg());
    input.close();
}
