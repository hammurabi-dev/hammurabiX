#ifndef __AP_ERROR_H__
#define __AP_ERROR_H__

#include <iostream>
#include <sstream>
#include <string>
using namespace std;

#define ap_err(...) cerr<<"ERROR! "<<__VA_ARGS__<<endl \
    <<" in file "<<__FILE__<<endl \
    <<" at line "<<__LINE__<<endl \
    <<" from function "<<__PRETTY_FUNCTION__<<endl<<endl \

#endif
