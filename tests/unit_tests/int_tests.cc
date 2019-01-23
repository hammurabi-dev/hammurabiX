#include <gtest/gtest.h>

#include <memory>
#include <cmath>
#include <integrator.h>
#include <grid.h>

TEST(integrator, boundary_check){
    Integrator test;
    double R0 = 10.;
    double R_lim = 1.0001*R0;
    EXPECT_FALSE(test.check_simulation_upper_limit(R0,R_lim));
    EXPECT_TRUE(test.check_simulation_lower_limit(R0,R_lim));
}

TEST(integrator, shell_ref_assembling){
    auto test = std::unique_ptr<Integrator> ();
    auto gint = std::unique_ptr<Grid_int> ();
    gint->init ("int_tests_01.xml");
    
    auto tmp = std::unique_ptr<Integrator::struct_shell> ();
    test->assemble_shell_ref(tmp.get(),gint.get(),1);
    unsigned int test_unsigned = 667;
    EXPECT_EQ(tmp->step,test_unsigned);
    gint->ec_r_max = 10.02;
    test->assemble_shell_ref(tmp.get(),gint.get(),1);
    test_unsigned += 2;
    EXPECT_EQ(tmp->step,test_unsigned);
    
    gint->init ("int_tests_02.xml");
    test
}
