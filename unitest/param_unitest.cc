#include<iostream>
#include<stdlib.h>
#include<param.h>

using namespace std;

int main(void){
	cout<<"param unitest"<<endl;
	
	Param par;
	cout<<"default constr ...... pass"<<endl;
	par.brnd_seed = 23;
	if(par.brnd_seed != 23){
		cout<<"assign value ...... fail"<<endl;
		exit(1);
	}
	cout<<"assign value ...... pass"<<endl;
		
	return EXIT_SUCCESS;
}
