/**
 * unit-test for class Integration
 * here we test auxiliary functions in Integrator class
 * LOS integration must be checked in integration test
 * feel free to add more rational testing blocks
 */
#include <gtest/gtest.h>

#include <memory>
#include <cmath>
#include <integrator.h>
#include <grid.h>


TEST(integrator, shell_radius){
    Integrator test;
    size_t total_shell = 200;
    double R0 = 10.;
    for(size_t i=1;i<=total_shell;++i){
        double radius_max = R0*pow(2.,-int(total_shell-i));
        double radius_min = R0*pow(2.,-int(total_shell-i+1));
        if(i==1) radius_min = 0.;
        EXPECT_LT(fabs(test.get_max_shell_radius(i,total_shell,R0) - radius_max),1e-10);
        EXPECT_LT(fabs(test.get_min_shell_radius(i,total_shell,R0) - radius_min),1e-10);
    }
}

TEST(integrator, boundary_check){
    Integrator test;
    double R0 = 10.;
    double R_lim = 1.0001*R0;
    EXPECT_FALSE(test.check_simulation_upper_limit(R0,R_lim));
    EXPECT_TRUE(test.check_simulation_lower_limit(R0,R_lim));
}

TEST(integrator, shell_ref_assembling){
    Integrator test;
    std::unique_ptr<Grid_int> gint = std::unique_ptr<Grid_int>(new Grid_int());
    gint->total_shell = 1;
    gint->ec_r_max = 10.;
    gint->radial_res = 0.03;
    std::unique_ptr<Integrator::struct_shell> tmp = std::unique_ptr<Integrator::struct_shell>(new Integrator::struct_shell);
    test.assemble_shell_ref(tmp.get(),gint.get(),1);
    unsigned int test_unsigned = 667;
    EXPECT_EQ(tmp->step,test_unsigned);
    gint->ec_r_max = 10.02;
    test.assemble_shell_ref(tmp.get(),gint.get(),1);
    test_unsigned += 2;
    EXPECT_EQ(tmp->step,test_unsigned);
}
