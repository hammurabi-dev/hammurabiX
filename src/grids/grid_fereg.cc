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

Grid_fereg::Grid_fereg(string file_name){
    unique_ptr<XMLDocument> doc = unique_ptr<XMLDocument> (new XMLDocument());
    doc->LoadFile(file_name.c_str());
    XMLElement *ptr {doc->FirstChildElement("root")->FirstChildElement("Fieldout")->FirstChildElement("fereg_grid")};
    read_permission = ptr->BoolAttribute("read");
    write_permission = ptr->BoolAttribute("write");
    if(read_permission or write_permission){
        filename = ptr->Attribute("filename");
        build_grid(doc.get());
    }
}

void Grid_fereg::build_grid(XMLDocument *doc){
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
    fe = unique_ptr<double[]> (new double[full_size]);
}

void Grid_fereg::export_grid(void){
    assert(!filename.empty());
    ofstream output(filename.c_str(), std::ios::out|std::ios::binary);
    assert(output.is_open());
    double tmp;
    for(decltype(full_size) i=0;i!=full_size;++i){
        assert(!output.eof());
        tmp = fe[i];
        assert(tmp>=0);
        output.write(reinterpret_cast<char*>(&tmp),sizeof(double));
    }
    output.close();
    exit(0);
}

void Grid_fereg::import_grid(void){
    assert(!filename.empty());
    ifstream input(filename.c_str(), std::ios::in|std::ios::binary);
    assert(input.is_open());
    double tmp;
    for(decltype(full_size) i=0;i!=full_size;++i){
        assert(!input.eof());
        input.read(reinterpret_cast<char *>(&tmp),sizeof(double));
        fe[i] = tmp;
    }
#ifndef NDEBUG
    auto eof = input.tellg();
    input.seekg (0, input.end);
#endif
    assert(eof==input.tellg());
    input.close();
}
