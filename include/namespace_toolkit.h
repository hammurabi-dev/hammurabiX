/*
 *@file: namespace_toolkit.h
 *@brief: numerical tools
 */
#ifndef HAMMURABI_TOOLKIT_H
#define HAMMURABI_TOOLKIT_H

#include <vec3.h>
#include <vector>
#include <string>
#include <cmath>
#include "cgs_units_file.h"

namespace toolkit {
    /*@get_perp2LOS/per2LOS
     * calculate perpendicular/parallel components w.r.t. LOS direction
     */
    double get_perp2LOS (const vec3_t<double> &,const double &,const double &);
    double get_par2LOS (const vec3_t<double> &,const double &,const double &);
    /*@get_intr_pol_ang
     * intrinsic polarization angle
     */
    double get_intr_pol_ang(const vec3_t<double> &,const double &,const double &);
    /*@get_LOS_unit_vec
     * get the unit vec, given the LOS direction
     */
    inline vec3_t<double> get_LOS_unit_vec(const double &the_los,const double &phi_los){
        return vec3_t<double> {cos(phi_los)*sin(the_los),
            sin(phi_los)*sin(the_los),
            cos(the_los)};
    }
    /*@cart_coord2cyl_coord
     * convert coordinate between cartesian and sylindrical frame
     */
    void cart_coord2cyl_coord(const vec3_t<double> &, double &, double &, double &);
    void cart_coord2cyl_coord(const vec3_t<double> &, vec3_t<double> &);
    /*@Cyl2Cart
     * map a vector from cylindrical frame into cartesian frame
     * two frames share the same origin point
     */
    inline void Cyl2Cart(const double &phi,const vec3_t<double> &input, vec3_t<double> &output){
        output = vec3_t<double> {cos(phi)*input.x - sin(phi)*input.y,
            sin(phi)*input.x + cos(phi)*input.y,
            input.z};
    }
    inline void Cyl2Cart(const double &phi,const double input[3],double output[3]){
        output[0] = cos(phi)*input[0] - sin(phi)*input[1];
        output[1] = sin(phi)*input[0] + cos(phi)*input[1];
        output[2] = input[2];
    }
    /*@versor
     * get versor of given vector
     */
    vec3_t<double> versor(const vec3_t<double> &);
    
    /*@Index3d
     * find index of 3D grid
     */
    inline std::size_t Index3d(const unsigned int &/* n1 */,const unsigned int &n2,const unsigned int &n3,const unsigned int &i,const unsigned int &j,const unsigned int &l){
        return (i*n2*n3 + j*n3 + l);
    }
    /*@Index4d
     * find index for 4D grid
     */
    inline std::size_t Index4d(const unsigned int &/* n1 */,const unsigned int &n2,const unsigned int &n3,const unsigned int &n4,const unsigned int &e,const unsigned int &i,const unsigned int &j,const unsigned int &l){
        return (e*n2*n3*n4 + i*n3*n4 + j*n4 + l);
    }
    
    double Mean (const double *,const std::size_t &); //for array
    double Mean(const std::vector<double> &);
    
    double Variance (const double *,const std::size_t &);
    double Variance(const std::vector<double> &);
    
    double Covariance (const double *,const double *,const std::size_t &);
    double Covariance(const std::vector<double> &,const std::vector<double> &);
    
    /*@Rank
     * convert an array/vector into rank array/vector
     */
    void Rank(double *,const std::size_t &);
    void Rank(std::vector<double> &);
    
    void Cart2LOS(const double &,const double &,const vec3_t<double> &,vec3_t<double> &);
    void LOS2Cart(const double &,const double &,const vec3_t<double> &,vec3_t<double> &);
    
    /*@temp_convert
     * converting brightness temp into thermal temp
     * with T_0 = 2.725K
     * formula get from Prog Theor Exp Phys (2014) 2014 (6): 06B109.
     */
    inline double temp_convert(const double &temp_br,const double &freq) {
        const double p {CGS_U_h_planck*freq/(CGS_U_kB*2.725)};
        return temp_br*(exp(p)-1.)*(exp(p)-1.)/(p*p*exp(p));;
    }
    /*@warp
     * convert cartesian coordiante into frame with galacitc warp
     * fields built in warpped frame as in ordinary frame
     */
    vec3_t<double> warp(const vec3_t<double> &);
    
    std::size_t random_seed(void);
    
    double timestamp(void);
}

#endif
