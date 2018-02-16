#include<iostream>
#include<param.h>
#include<stdlib.h>
#include<memory>
using namespace std;

int main(void){
	cout<<"hello unitest"<<endl;
	unique_ptr<Param> par;
	par->brnd_seed = 23;
	cout<<"brnd_seed "<<par->brnd_seed<<endl;	
	return EXIT_SUCCESS;
}
