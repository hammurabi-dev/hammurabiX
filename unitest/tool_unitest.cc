///
/// unit-test for namespace toolkit
/// feel free to add more rational testing blocks
///
#include<iostream>
#include<cstdlib>
#include<cmath>
#include<namespace_toolkit.h>
#include<vec3.h>
#include<cgs_units_file.h>
using namespace std;
using namespace toolkit;

template<typename T>
bool compare(const T &,const T &,const double &);

template<typename T>
bool compare(const T &a,const T &b,const double &precision){
    return (abs(a-b)<precision);
}

int main(void){
    cout<<"toolkit unitest"<<endl;
    
    double theta[3] = {0.,90.*CGS_U_rad,90.*CGS_U_rad};
    double phi[3] = {0.,180.*CGS_U_rad,270.*CGS_U_rad};
    vec3_t<double> A {1.,0.,0.};
    //vec3_t<double> B {0.,0.,1.};
    bool test=true;
    
    // test get_LOS_unit_vec
    test &= (compare(get_LOS_unit_vec(theta[0],phi[0]).x,0.,1.e-10) and
             compare(get_LOS_unit_vec(theta[0],phi[0]).y,0.,1.e-10) and
             compare(get_LOS_unit_vec(theta[0],phi[0]).z,1.,1.e-10));
    test &= (compare(get_LOS_unit_vec(theta[1],phi[1]).x,-1.,1.e-10) and
             compare(get_LOS_unit_vec(theta[1],phi[1]).y,0.,1.e-10) and
             compare(get_LOS_unit_vec(theta[1],phi[1]).z,0.,1.e-10));
    test &= (compare(get_LOS_unit_vec(theta[2],phi[2]).x,0.,1.e-10) and
             compare(get_LOS_unit_vec(theta[2],phi[2]).y,-1.,1.e-10) and
             compare(get_LOS_unit_vec(theta[2],phi[2]).z,0.,1.e-10));
    if (!test){
        cerr<<"ERR:"<<__FILE__
        <<" : in function "<<__func__<<endl
        <<" at line "<<__LINE__<<endl
        <<"get_LOS_unit_vec ...... fail"<<endl;
        exit(1);
    }
    
    // test get_par2LOS
    
    test &= compare(get_par2LOS(A,theta[0],phi[0]),0.,1.e-10);
    test &= compare(get_par2LOS(A,theta[1],phi[1]),-1.,1.e-10);
    test &= compare(get_par2LOS(A,theta[2],phi[2]),0.,1.e-10);
    if (!test){
        cerr<<"ERR:"<<__FILE__
        <<" : in function "<<__func__<<endl
        <<" at line "<<__LINE__<<endl
        <<"get_par2LOS ...... fail"<<endl;
        exit(1);
    }
    
    // test get_per2LOS
    test &= compare(get_perp2LOS(A,theta[0],phi[0]),1.,1.e-10);
    test &= compare(get_perp2LOS(A,theta[1],phi[1]),0.,1.e-10);
    test &= compare(get_perp2LOS(A,theta[2],phi[2]),1.,1.e-10);
    if (!test){
        cerr<<"ERR:"<<__FILE__
        <<" : in function "<<__func__<<endl
        <<" at line "<<__LINE__<<endl
        <<"get_per2LOS ...... fail"<<endl;
        exit(1);
    }
    
    // test get_intr_pol_ang
    test &= compare(get_intr_pol_ang(A,theta[0],phi[0]),-90.*CGS_U_rad,1.e-10);
    test &= compare(get_intr_pol_ang(A,theta[2],phi[2]),180.*CGS_U_rad,1.e-10);
    if (!test){
        cerr<<"ERR:"<<__FILE__
        <<" : in function "<<__func__<<endl
        <<" at line "<<__LINE__<<endl
        <<"get_intr_pol_ang ...... fail"<<endl;
        exit(1);
    }
    
    
    // if all testing blocks pass
    cout<<"namespace toolkit ...... pass"<<endl;
    return EXIT_SUCCESS;
}
