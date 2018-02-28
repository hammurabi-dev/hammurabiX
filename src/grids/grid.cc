#include <sstream>
#include <iostream>
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

using namespace tinyxml2;
using namespace std;

void Grid::build_grid(XMLDocument *){
    cerr<<"ERR:"<<__FILE__
    <<" : in function "<<__func__<<endl
    <<" at line "<<__LINE__<<endl
    <<"DYNAMIC BINDING FAILURE"<<endl;
    exit(1);
}

void Grid::export_grid(void){
    cerr<<"ERR:"<<__FILE__
    <<" : in function "<<__func__<<endl
    <<" at line "<<__LINE__<<endl
    <<"DYNAMIC BINDING FAILURE"<<endl;
    exit(1);
}

void Grid::import_grid(void){
    cerr<<"ERR:"<<__FILE__
    <<" : in function "<<__func__<<endl
    <<" at line "<<__LINE__<<endl
    <<"DYNAMIC BINDING FAILURE"<<endl;
    exit(1);
}

// auxiliary functions
std::string Grid::FetchString(XMLElement* el, string obj){
    return el->FirstChildElement(obj.c_str())->Attribute("value");
}

int Grid::FetchInt(XMLElement* el, string obj){
    return el->FirstChildElement(obj.c_str())->IntAttribute("value");
}

unsigned int Grid::FetchUnsigned(XMLElement* el, string obj){
    return el->FirstChildElement(obj.c_str())->UnsignedAttribute("value");
}

bool Grid::FetchBool(XMLElement* el, string obj){
    return el->FirstChildElement(obj.c_str())->BoolAttribute("cue");
}

double Grid::FetchDouble(XMLElement* el, string obj){
    return el->FirstChildElement(obj.c_str())->DoubleAttribute("value");
}

//END
