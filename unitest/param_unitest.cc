///
/// unit-test for class Param
/// feel free to add more rational testing blocks
///
#include <iostream>
#include <cstdlib>
#include <param.h>
#include <ap_err.h>
using namespace std;

int main(void){
	cout<<"param unitest"<<endl;
	
	Param par;
	
	par.brnd_seed = 23;
	if(par.brnd_seed != 23){
        ap_err("assign value ...... fail");
		exit(1);
	}
    
    
    // if all testing blocks pass
	cout<<"class Param ...... pass"<<endl;
	return EXIT_SUCCESS;
}
