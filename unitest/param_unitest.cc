/**
 * unit-test for class Param
 * feel free to add more rational testing blocks
 */
#include <iostream>
#include <cstdlib>
#include <param.h>
#include <cassert>
using namespace std;

int main(void){
    cout<<endl<<"param unitest"<<endl;
    
    Param par;
    
    par.brnd_seed = 23;
    assert(par.brnd_seed==23);
    
    
    // if all testing blocks pass
    cout<<"class Param ...... pass"<<endl<<endl;
    return EXIT_SUCCESS;
}
