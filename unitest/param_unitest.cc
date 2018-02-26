///
/// unit-test for class Param
/// feel free to add more rational testing blocks
///
#include<iostream>
#include<cstdlib>
#include<param.h>

using namespace std;

int main(void){
	cout<<"param unitest"<<endl;
	
	Param par;
	cerr<<"default constr ...... pass"<<endl;
	par.brnd_seed = 23;
	if(par.brnd_seed != 23){
		cerr<<"assign value ...... fail"<<endl;
		exit(1);
	}
    
    
    // if all testing blocks pass
	cout<<"class Param ...... pass"<<endl;
	return EXIT_SUCCESS;
}
