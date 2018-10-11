#include <sstream>
#include <fstream>
#include <memory>
#include <array>
#include <string>
#include <vector>
#include <cassert>

#include <tinyxml2.h>

#include <grid.h>
#include <cgs_units_file.h>
#include <namespace_toolkit.h>

Grid_cre::Grid_cre (const std::string &file_name){
    std::unique_ptr<tinyxml2::XMLDocument> doc = toolkit::loadxml (file_name);
    tinyxml2::XMLElement *ptr {toolkit::tracexml (doc.get(),{"Fieldout"})};
    read_permission = toolkit::fetchbool (ptr,"read","cre_grid");
    write_permission = toolkit::fetchbool (ptr,"write","cre_grid");
    // build up grid when have read or write permission
    if (read_permission or write_permission){
        filename = toolkit::fetchstring (ptr,"filename","cre_grid");
        build_grid (doc.get());
    }
}

void Grid_cre::build_grid (tinyxml2::XMLDocument *doc){
    tinyxml2::XMLElement *ptr {toolkit::tracexml (doc,{"Grid","Box_CRE"})};
    E_min = CGS_U_GeV*toolkit::fetchdouble (ptr,"value","E_min");
    E_max = CGS_U_GeV*toolkit::fetchdouble (ptr,"value","E_max");
    nE = toolkit::fetchunsigned (ptr,"value","nE");
    E_fact = std::log(E_max/E_min)/(nE-1);
    // spatial 3D
    nz = toolkit::fetchunsigned (ptr,"value","nz");
    nx = toolkit::fetchunsigned (ptr,"value","nx");
    ny = toolkit::fetchunsigned (ptr,"value","ny");
    x_max = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","x_max");
    x_min = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","x_min");
    y_max = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","y_max");
    y_min = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","y_min");
    z_max = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","z_max");
    z_min = CGS_U_kpc*toolkit::fetchdouble (ptr,"value","z_min");
    cre_size = nE*nx*ny*nz;
    cre_flux = std::make_unique<double[]> (cre_size);
}

void Grid_cre::export_grid (){
    assert (!filename.empty());
    std::ofstream output (filename.c_str(),std::ios::out|std::ios::binary);
    assert (output.is_open());
    double tmp;
    for (decltype(cre_size) i=0;i!=cre_size;++i){
        assert (!output.eof());
        tmp = cre_flux[i];
        output.write (reinterpret_cast<char*>(&tmp),
                      sizeof(double));
    }
    output.close();
    exit (0);
}

void Grid_cre::import_grid (){
    assert (!filename.empty());
    std::ifstream input (filename.c_str(),std::ios::in|std::ios::binary);
    assert (input.is_open());
    double tmp;
    for (decltype(cre_size) i=0;i!=cre_size;++i){
        assert (!input.eof());
        input.read (reinterpret_cast<char *>(&tmp),
                    sizeof(double));
        cre_flux[i] = tmp;
    }
#ifndef NDEBUG
    auto eof = input.tellg();
    input.seekg (0,input.end);
#endif
    assert (eof==input.tellg());
    input.close();
}
