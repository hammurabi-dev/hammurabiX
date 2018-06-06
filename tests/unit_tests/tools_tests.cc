/**
 * unit tests for namespace toolkit
 * feel free to add more rational testing blocks
 */

#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <vector>
#include <fftw3.h>
#include <namespace_toolkit.h>
#include <vec3.h>
#include <cgs_units_file.h>
#include <tinyxml2.h>

TEST(toolkit, los_versor){
    const double theta[3] = {0.,90.*CGS_U_rad,90.*CGS_U_rad};
    const double phi[3] = {0.,180.*CGS_U_rad,270.*CGS_U_rad};
    
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

TEST(toolkit, versor){
    const vec3_t<double> tmp {0.1,0.0000001,0.0000001};
    EXPECT_LT(fabs(crossprod(tmp,toolkit::versor(tmp)).Length() - 0.),1e-10);
}

TEST(toolkit, par2los){
    const double theta[3] = {0.,90.*CGS_U_rad,90.*CGS_U_rad};
    const double phi[3] = {0.,180.*CGS_U_rad,270.*CGS_U_rad};
    const vec3_t<double> A {1.,0.,0.};
    
    EXPECT_LT(fabs(toolkit::par2los(A,theta[0],phi[0]) - 0.),1e-10);
    EXPECT_LT(fabs(toolkit::par2los(A,theta[1],phi[1]) + 1.),1e-10);
    EXPECT_LT(fabs(toolkit::par2los(A,theta[2],phi[2]) - 0.),1e-10);
}

TEST(toolkit, perp2los){
    const double theta[3] = {0.,90.*CGS_U_rad,90.*CGS_U_rad};
    const double phi[3] = {0.,180.*CGS_U_rad,270.*CGS_U_rad};
    const vec3_t<double> A {1.,0.,0.};
    
    EXPECT_LT(fabs(toolkit::perp2los(A,theta[0],phi[0]) - 1.),1e-10);
    EXPECT_LT(fabs(toolkit::perp2los(A,theta[1],phi[1]) + 0.),1e-10);
    EXPECT_LT(fabs(toolkit::perp2los(A,theta[2],phi[2]) - 1.),1e-10);
}

TEST(toolkit, intr_pol_ang){
    const double theta[3] = {0.,90.*CGS_U_rad,90.*CGS_U_rad};
    const double phi[3] = {0.,180.*CGS_U_rad,270.*CGS_U_rad};
    const vec3_t<double> A {1.,0.,0.};
    
    EXPECT_LT(fabs(toolkit::intr_pol_ang(A,theta[0],phi[0]) + 90.*CGS_U_rad),1e-10);
    EXPECT_LT(fabs(toolkit::intr_pol_ang(A,theta[2],phi[2]) - 180.*CGS_U_rad),1e-10);
}

TEST(toolkit, cart_coord2cyl_coord){
    const vec3_t<double> A {1.,0.,0.};
    const vec3_t<double> B {0.,0.,1.};
    const vec3_t<double> C {0.,1.,0.};
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
    
    double tmp_r, tmp_phi, tmp_z;
    
    toolkit::cart_coord2cyl_coord(A,tmp_r,tmp_phi,tmp_z);
    
    EXPECT_LT(fabs(tmp_r - 1.),1e-10);
    EXPECT_LT(fabs(tmp_phi - 0.),1e-10);
    EXPECT_LT(fabs(tmp_z - 0.),1e-10);
    
    toolkit::cart_coord2cyl_coord(B,tmp_r,tmp_phi,tmp_z);
    
    EXPECT_LT(fabs(tmp_r - 0.),1e-10);
    EXPECT_LT(fabs(tmp_z - 1.),1e-10);
    
    toolkit::cart_coord2cyl_coord(C,tmp_r,tmp_phi,tmp_z);
    
    EXPECT_LT(fabs(tmp_r - 1.),1e-10);
    EXPECT_LT(fabs(tmp_phi - 90.*CGS_U_rad),1e-10);
    EXPECT_LT(fabs(tmp_z - 0.),1e-10);
}

TEST(toolkit, index3d){
    std::size_t test_idx {53};
    EXPECT_EQ(toolkit::Index3d(0,5,4,2,3,1),test_idx);
}

TEST(toolkit, index4d){
    std::size_t test_idx {461};
    EXPECT_EQ(toolkit::Index4d(0,4,7,3,5,1,6,2),test_idx);
}

TEST(toolkit, mean){
    const double test_array[3] = {3,4,5};
    EXPECT_DOUBLE_EQ(toolkit::Mean(test_array,3),4.);
    
    const std::vector<double> test_vector {3,4,5};
    EXPECT_DOUBLE_EQ(toolkit::Mean(test_vector),4.);
}

TEST(toolkit, variance){
    const double test_array[3] = {3,4,5};
    EXPECT_DOUBLE_EQ(toolkit::Variance(test_array,3),2./3.);
    
    const std::vector<double> test_vector {3,4,5};
    EXPECT_DOUBLE_EQ(toolkit::Variance(test_vector),2./3.);
}

TEST(toolkit, covariance){
    const double test_array1[3] = {1,2,3};
    const double test_array2[3] = {3,4,5};
    EXPECT_DOUBLE_EQ(toolkit::Covariance(test_array1,test_array2,3),2./3.);
    
    const std::vector<double> test_vector1 {1,2,3};
    const std::vector<double> test_vector2 {3,4,5};
    EXPECT_DOUBLE_EQ(toolkit::Covariance(test_vector1,test_vector2),2./3.);
}

TEST(toolkit, complex_stripping){
    const auto test_complex = fftw_alloc_complex(2);
    test_complex[0][0] = 1.; test_complex[0][1] = 2.;
    test_complex[1][0] = 3.; test_complex[1][1] = 4.;
    double test_real[2];
    double test_imag[2];
    
    toolkit::complex2real(test_complex,test_real,2);
    EXPECT_DOUBLE_EQ(test_real[0],test_complex[0][0]);
    EXPECT_DOUBLE_EQ(test_real[1],test_complex[1][0]);
    
    toolkit::complex2imag(test_complex,test_imag,2);
    EXPECT_DOUBLE_EQ(test_imag[0],test_complex[0][1]);
    EXPECT_DOUBLE_EQ(test_imag[1],test_complex[1][1]);
    
    toolkit::complex2rni(test_complex,test_real,test_imag,2);
    EXPECT_DOUBLE_EQ(test_real[0],test_complex[0][0]);
    EXPECT_DOUBLE_EQ(test_real[1],test_complex[1][0]);
    EXPECT_DOUBLE_EQ(test_imag[0],test_complex[0][1]);
    EXPECT_DOUBLE_EQ(test_imag[1],test_complex[1][1]);
    
    fftw_free(test_complex);
}

TEST(toolkit, xml_parser){
    auto doc = toolkit::loadxml("reference/example.xml");
    
    tinyxml2::XMLElement* el = toolkit::tracexml(doc.get(),{"double"});
    EXPECT_EQ(toolkit::FetchDouble(el,"value"),3.14);
    
    el = toolkit::tracexml(doc.get(),{"integer"});
    unsigned int test_unsigned = 23;
    EXPECT_EQ(toolkit::FetchUnsigned(el,"value"),test_unsigned);
    
    el = toolkit::tracexml(doc.get(),{});
    std::string test_string = "string";
    EXPECT_EQ(toolkit::FetchString(el,"value","string"),test_string);
    
    el = toolkit::tracexml(doc.get(),{});
    bool test_bool = true;
    EXPECT_EQ(toolkit::FetchBool(el,"value","bool"),test_bool);
}
