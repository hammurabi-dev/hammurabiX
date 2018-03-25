#include <timer.h>
#include <cassert>
#include <iostream>
#include <unistd.h>
#include <cstdlib>
#include <cmath>
using namespace std;

template<typename T>
bool compare(const T &,const T &,const double &);

template<typename T>
bool compare(const T &a,const T &b,const double &precision){
    return (fabs(a-b)<precision);
}

int main(void){
    cout<<endl<<"timer unitest"<<endl;
    
    Timer t;
    t.start("main");
    
    t.start("sub");
    usleep(3000000);
    t.stop("sub");
    
    usleep(3000000);
    t.stop("main");
    
    assert(compare(t.record.find("sub")->second.second,3000.,1.));
    assert(compare(t.record.find("main")->second.second,6000.,1.));
    
    // if all testing blocks pass
    cout<<"timer ...... pass"<<endl<<endl;
    return EXIT_SUCCESS;
}
