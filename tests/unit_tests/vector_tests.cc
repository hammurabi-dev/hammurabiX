// unit tests for our vector class
// feel free to add more rational testing blocks

#include <gtest/gtest.h>

#include <iostream>
#include <cmath>
#include <memory>
#include <hvec.h>

TEST (vector, basic){
    
    // argument ctor
    hvec<1,float> vec1 (0.0);
    EXPECT_EQ (vec1[0],0.);
    
    // implicit ctor
    vec1 = 0.1;
    EXPECT_EQ (vec1[0],float(0.1));
    
    // cp assign, operator==
    hvec<1,float> vec1_cpa (0.2);
    vec1 = vec1_cpa;
    EXPECT_EQ (vec1[0],float(0.2));
    
    // cp ctor
    hvec<1,float> vec1_cpc (vec1);
    EXPECT_EQ (vec1_cpc[0],vec1[0]);
    
    // mv ctor
    hvec<1,float> vec1_mvc = std::move(vec1_cpc);
    EXPECT_EQ (vec1_mvc[0],vec1[0]);
    
    // list ctor
    hvec<2,double> vec2 {0.3,0.4};
    EXPECT_EQ (vec2[0],double(0.3));
    EXPECT_EQ (vec2[1],double(0.4));
    
    // operator +
    auto vecp = vec2 + hvec<2,double> {0.4,0.8};
    EXPECT_EQ (vecp[0],vec2[0]+double(0.4));
    EXPECT_EQ (vecp[1],vec2[1]+double(0.8));
    
    // operator +=
    vec2 += vec2;
    EXPECT_EQ (vec2[0],double(0.3)+double(0.3));
    EXPECT_EQ (vec2[1],double(0.4)+double(0.4));
    
    // resign
    vec2 = {0.3,0.4};
    EXPECT_EQ (vec2[0],double(0.3));
    EXPECT_EQ (vec2[1],double(0.4));
    
    // operator -
    auto vecm = vec2 - hvec<2,float> {0.3,0.8};
    EXPECT_EQ (vecm[0],double(0.3)-static_cast<double>(float(0.3)));
    EXPECT_EQ (vecm[1],double(0.4)-static_cast<double>(float(0.8)));
    
    // operator -=
    vec2 = {0.3,0.4};
    vec2 -= vec2;
    EXPECT_EQ (vec2[0],double(0.0));
    EXPECT_EQ (vec2[1],double(0.0));
    
    // operator *
    vec2 = {0.3,0.4};
    vec2 = vec2*3;
    EXPECT_EQ (vec2[0],double(0.3)*double(3));
    EXPECT_EQ (vec2[1],double(0.4)*double(3));
    
    // operator *=
    vec2 = {0.3,0.4};
    vec2 *= 0.3;
    EXPECT_EQ (vec2[0],double(0.3)*double(0.3));
    EXPECT_EQ (vec2[1],double(0.4)*double(0.3));
    
    // operator /
    vec2 = {0.27,0.36};
    vec2 = vec2/0.3;
    EXPECT_EQ (vec2[0],double(0.27)/double(0.3));
    EXPECT_EQ (vec2[1],double(0.36)/double(0.3));
    
    // operator /=
    vec2 = {0.9,1.2};
    vec2 /= 3;
    EXPECT_EQ (vec2[0],double(0.9)/double(3));
    EXPECT_EQ (vec2[1],double(1.2)/double(3));
    
    // operator !=
    auto test = hvec<2,double>(1.,1.);
    EXPECT_TRUE (vec2!=test);
    
    hvec<3,int> vec3 {1,2,3};
    
    // function lengthsq
    EXPECT_EQ (vec3.lengthsq(),double(14.0));
    
    // function length
    EXPECT_EQ (vec3.length(),std::sqrt(double(14.0)));
    
    // function versor
    auto norm = vec3.versor();
    EXPECT_EQ (norm[0],double(1./std::sqrt(14.)));
    EXPECT_EQ (norm[1],double(2./std::sqrt(14.)));
    EXPECT_EQ (norm[2],double(3./std::sqrt(14.)));
    
    // function flip
    auto flipt = hvec<3,int>(-1,-2,-3);
    vec3.flip();
    EXPECT_TRUE (vec3==flipt);
    
    hvec<3,double> vec3a {0.1,0.2,0.3};
    hvec<3,double> vec3b {0.4,0.5,0.6};
    
    // function dotprod
    auto rslt = vec3a.dotprod (vec3b);
    EXPECT_EQ (rslt,0.32);
    
    // function crossprod
    hvec<3,double> prodt (0,0,0);
    prodt[0] = double(0.2)*double(0.6)-double(0.3)*double(0.5);
    prodt[1] = double(0.3)*double(0.4)-double(0.1)*double(0.6);
    prodt[2] = double(0.1)*double(0.5)-double(0.2)*double(0.4);
    auto vec3c = vec3a.crossprod (vec3b);
    EXPECT_TRUE (vec3c==prodt);
}
