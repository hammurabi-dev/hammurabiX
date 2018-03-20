/**
 * unit-test for namespace toolkit
 * feel free to add more rational testing blocks
 */
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <namespace_toolkit.h>
#include <vec3.h>
#include <cgs_units_file.h>
#include <cassert>

using namespace std;
using namespace toolkit;

template<typename T>
bool compare(const T &,const T &,const double &);

template<typename T>
bool compare(const T &a,const T &b,const double &precision){
    return (fabs(a-b)<precision);
}

int main(void){
    cout<<endl<<"toolkit unitest"<<endl;
    
    double theta[3] = {0.,90.*CGS_U_rad,90.*CGS_U_rad};
    double phi[3] = {0.,180.*CGS_U_rad,270.*CGS_U_rad};
    vec3_t<double> A {1.,0.,0.};
    vec3_t<double> B {0.,0.,1.};
    vec3_t<double> C {0.,1.,0.};
    vec3_t<double> tmp;
    
    // test get_LOS_unit_vec
    assert(compare(get_LOS_unit_vec(theta[0],phi[0]).x,0.,1e-10));
    assert(compare(get_LOS_unit_vec(theta[0],phi[0]).y,0.,1e-10));
    assert(compare(get_LOS_unit_vec(theta[0],phi[0]).z,1.,1e-10));
    assert(compare(get_LOS_unit_vec(theta[1],phi[1]).x,-1.,1e-10));
    assert(compare(get_LOS_unit_vec(theta[1],phi[1]).y,0.,1e-10));
    assert(compare(get_LOS_unit_vec(theta[1],phi[1]).z,0.,1e-10));
    assert(compare(get_LOS_unit_vec(theta[2],phi[2]).x,0.,1e-10));
    assert(compare(get_LOS_unit_vec(theta[2],phi[2]).y,-1.,1e-10));
    assert(compare(get_LOS_unit_vec(theta[2],phi[2]).z,0.,1e-10));
    
    // test get_par2LOS
    
    assert(compare(get_par2LOS(A,theta[0],phi[0]),0.,1e-10));
    assert(compare(get_par2LOS(A,theta[1],phi[1]),-1.,1e-10));
    assert(compare(get_par2LOS(A,theta[2],phi[2]),0.,1e-10));
    
    // test get_per2LOS
    assert(compare(get_perp2LOS(A,theta[0],phi[0]),1.,1e-10));
    assert(compare(get_perp2LOS(A,theta[1],phi[1]),0.,1e-10));
    assert(compare(get_perp2LOS(A,theta[2],phi[2]),1.,1e-10));
    
    // test get_intr_pol_ang
    assert(compare(get_intr_pol_ang(A,theta[0],phi[0]),-90.*CGS_U_rad,1e-10));
    assert(compare(get_intr_pol_ang(A,theta[2],phi[2]),180.*CGS_U_rad,1e-10));
    
    // test cart_coord2cyl_coord
    cart_coord2cyl_coord(A,tmp);
    assert(compare(tmp.x,1.,1e-10) and compare(tmp.y,0.,1e-10));
    assert(compare(tmp.z,0.,1.e-10));
    
    cart_coord2cyl_coord(B,tmp);
    assert(compare(tmp.x,0.,1e-10));
    assert(compare(tmp.z,1.,1e-10));
    
    cart_coord2cyl_coord(C,tmp);
    assert(compare(tmp.x,1.,1e-10));
    assert(compare(tmp.y,90.*CGS_U_rad,1e-10));
    assert(compare(tmp.z,0.,1e-10));
    
    // test versor
    tmp = vec3_t<double> {0.1,0.00002,0.0000003};
    assert(compare(crossprod(tmp,versor(tmp)).Length(),0.,1.e-10));
    
    
    // if all testing blocks pass
    cout<<"namespace toolkit ...... pass"<<endl<<endl;
    return EXIT_SUCCESS;
}
