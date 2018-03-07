///
/// unit-test for class Integration
/// feel free to add more rational testing blocks
///
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <integrator.h>
#include <grid.h>
#include <ap_err.h>
using namespace std;

int main(void){
    cout<<endl<<"integrator unitest"<<endl;
    
    Integrator test;
    
    // testing get_max/min_shell_radius
    size_t total_shell = 200;
    double R0 = 10.;
    for(size_t i=1;i<=total_shell;++i){
        double radius_max = R0*pow(2.,-int(total_shell-i));
        double radius_min = R0*pow(2.,-int(total_shell-i+1));
        if(i==1) radius_min = 0.;
        if(test.get_max_shell_radius(i,total_shell,R0)!=radius_max){
            ap_err("get_max_shell_radius ...... fail");
            exit(1);
        }
        if(test.get_min_shell_radius(i,total_shell,R0)!=radius_min){
            ap_err("get_min_shell_radius ...... fail");
            exit(1);
        }
    }
    
    // testing boundary check
    double R_lim = 1.0001*R0;
    if(test.check_simulation_upper_limit(R0,R_lim)){
        ap_err("check_simulation_upper_limit ...... fail");
        exit(1);
    }
    if(not test.check_simulation_lower_limit(R0,R_lim)){
        ap_err("check_simulation_lower_limit ...... fail");
        exit(1);
    }
    
    // testing shell_ref assembling
    // need Grid_int class
    unique_ptr<Grid_int> gint = unique_ptr<Grid_int>(new Grid_int());
    gint->total_shell = 1;
    gint->ec_r_max = 10.;
    gint->radial_res = 0.03;
    unique_ptr<Integrator::struct_shell> tmp = unique_ptr<Integrator::struct_shell>(new Integrator::struct_shell);
    test.assemble_shell_ref(tmp.get(),gint.get(),1);
    if(tmp->step!=667){
        ap_err("assemble_shell_ref ...... fail");
        exit(1);
    }
    gint->ec_r_max = 10.02;
    test.assemble_shell_ref(tmp.get(),gint.get(),1);
    if(tmp->step!=669){
        ap_err("assemble_shell_ref ...... fail");
        exit(1);
    }
    
    // LOS integration will be tested through integrated tests
    
    // if all testing blocks pass
    cout<<"class Integration ...... pass"<<endl<<endl;
    return EXIT_SUCCESS;
}
