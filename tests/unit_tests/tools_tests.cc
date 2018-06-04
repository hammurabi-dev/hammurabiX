/**
 * unit tests for namespace toolkit
 * feel free to add more rational testing blocks
 */

#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <namespace_toolkit.h>
#include <vec3.h>
#include <cgs_units_file.h>
#include <tinyxml2.h>

using namespace tinyxml2;

TEST(toolkit, los_versor){
    double theta[3] = {0.,90.*CGS_U_rad,90.*CGS_U_rad};
    double phi[3] = {0.,180.*CGS_U_rad,270.*CGS_U_rad};
    
    EXPECT_LT(fabs(toolkit::los_versor(theta[0],phi[0]).x - 0.),1e-10);
    EXPECT_LT(fabs(toolkit::los_versor(theta[0],phi[0]).y - 0.),1e-10);
    EXPECT_LT(fabs(toolkit::los_versor(theta[0],phi[0]).z - 1.),1e-10);
    EXPECT_LT(fabs(toolkit::los_versor(theta[1],phi[1]).x + 1.),1e-10);
    EXPECT_LT(fabs(toolkit::los_versor(theta[1],phi[1]).y - 0.),1e-10);
    EXPECT_LT(fabs(toolkit::los_versor(theta[1],phi[1]).z - 0.),1e-10);
    EXPECT_LT(fabs(toolkit::los_versor(theta[2],phi[2]).x - 0.),1e-10);
    EXPECT_LT(fabs(toolkit::los_versor(theta[2],phi[2]).y + 1.),1e-10);
    EXPECT_LT(fabs(toolkit::los_versor(theta[2],phi[2]).z - 0.),1e-10);
}

TEST(toolkit, par2los){
    double theta[3] = {0.,90.*CGS_U_rad,90.*CGS_U_rad};
    double phi[3] = {0.,180.*CGS_U_rad,270.*CGS_U_rad};
    vec3_t<double> A {1.,0.,0.};
    
    EXPECT_LT(fabs(toolkit::par2los(A,theta[0],phi[0]) - 0.),1e-10);
    EXPECT_LT(fabs(toolkit::par2los(A,theta[1],phi[1]) + 1.),1e-10);
    EXPECT_LT(fabs(toolkit::par2los(A,theta[2],phi[2]) - 0.),1e-10);
}

TEST(toolkit, intr_pol_ang){
    double theta[3] = {0.,90.*CGS_U_rad,90.*CGS_U_rad};
    double phi[3] = {0.,180.*CGS_U_rad,270.*CGS_U_rad};
    vec3_t<double> A {1.,0.,0.};
    
    EXPECT_LT(fabs(toolkit::intr_pol_ang(A,theta[0],phi[0]) + 90.*CGS_U_rad),1e-10);
    EXPECT_LT(fabs(toolkit::intr_pol_ang(A,theta[2],phi[2]) - 180.*CGS_U_rad),1e-10);
}

TEST(toolkit, cart_coord2cyl_coord){
    vec3_t<double> A {1.,0.,0.};
    vec3_t<double> B {0.,0.,1.};
    vec3_t<double> C {0.,1.,0.};
    vec3_t<double> tmp;
    
    toolkit::cart_coord2cyl_coord(A,tmp);
    
    EXPECT_LT(fabs(tmp.x - 1.),1e-10);
    EXPECT_LT(fabs(tmp.y - 0.),1e-10);
    EXPECT_LT(fabs(tmp.z - 0.),1e-10);
    
    toolkit::cart_coord2cyl_coord(B,tmp);
    
    EXPECT_LT(fabs(tmp.x - 0.),1e-10);
    EXPECT_LT(fabs(tmp.z - 1.),1e-10);
    
    toolkit::cart_coord2cyl_coord(C,tmp);
    
    EXPECT_LT(fabs(tmp.x - 1.),1e-10);
    EXPECT_LT(fabs(tmp.y - 90.*CGS_U_rad),1e-10);
    EXPECT_LT(fabs(tmp.z - 0.),1e-10);
}

TEST(toolkit, versor){
    vec3_t<double> tmp = vec3_t<double> {0.1,0.0000001,0.0000001};
    EXPECT_LT(fabs(crossprod(tmp,toolkit::versor(tmp)).Length() - 0.),1e-10);
}

TEST(toolkit, xml_parser){
    auto doc = toolkit::loadxml("reference/example.xml");
    
    XMLElement* el = toolkit::tracexml(doc.get(),{"double"});
    EXPECT_EQ(toolkit::FetchDouble(el,"value"),3.14);
    
    el = toolkit::tracexml(doc.get(),{"integer"});
    unsigned int test_unsigned = 23;
    EXPECT_EQ(toolkit::FetchUnsigned(el,"value"),test_unsigned);
    
    el = toolkit::tracexml(doc.get(),{});
    std::string test_string = "string";
    EXPECT_EQ(toolkit::FetchString(el,"value","string"),test_string);
}
