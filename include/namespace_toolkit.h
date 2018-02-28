///
/// numerical tools
///
#ifndef HAMMURABI_TOOLKIT_H
#define HAMMURABI_TOOLKIT_H

#include <vec3.h>
#include <vector>
#include <string>
#include <cmath>
#include <cgs_units_file.h>

namespace toolkit {
    ///
    /// perpendicular component of a vector wrt LOS direction
    /// 1st argument: vector in Cartesian frame
    /// 2nd argument: polar angle (in rad)
    /// 3rd argument: azimuthal angle (in rad)
    ///
    double get_perp2LOS (const vec3_t<double> &,const double &,const double &);
    ///
    /// (signed) parallel component of a vector wrt LOS direction
    /// 1st argument: vector in Cartesian frame
    /// 2nd argument: polar angle (in rad)
    /// 3rd argument: azimuthal angle (in rad)
    ///
    double get_par2LOS (const vec3_t<double> &,const double &,const double &);
    ///
    /// intrinsic polarization angle (in rad)
    /// 1st argument: magnetic field in Cartesian frame
    /// 2nd argument: polar angle (in rad) of LOS direction
    /// 3rd argument: azimuthal angle (in rad) of LOS direction
    /// use with caution, since vector can be parallel to LOS direction
    ///
    double get_intr_pol_ang(const vec3_t<double> &,const double &,const double &);
    ///
    /// Carteisan unit vector of given LOS direction
    /// first argument: polar angle (in rad)
    /// second argument: azimuthal angle (in rad)
    ///
    inline vec3_t<double> get_LOS_unit_vec(const double &the_los,const double &phi_los){
        return vec3_t<double> {cos(phi_los)*sin(the_los),
            sin(phi_los)*sin(the_los),
            cos(the_los)};
    }
    ///
    /// convert coordinate from cartesian to cylindrical frame
    ///
    void cart_coord2cyl_coord(const vec3_t<double> &, double &, double &, double &);
    ///
    /// convert coordinate from cylindrical to cartesian frame
    ///
    void cart_coord2cyl_coord(const vec3_t<double> &, vec3_t<double> &);
    ///
    /// map a vector from Galactic centric cylindrical frame into Galactic centric cartesian frame
    ///
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
    ///
    /// get versor of given vector
    ///
    vec3_t<double> versor(const vec3_t<double> &);
    
    ///
    /// find index of 3D grid
    ///
    inline std::size_t Index3d(const std::size_t &/* n1 */,const std::size_t &n2,const std::size_t &n3,const std::size_t &i,const std::size_t &j,const std::size_t &l){
        return (i*n2*n3 + j*n3 + l);
    }
    ///
    /// find index for 4D grid
    ///
    inline std::size_t Index4d(const std::size_t &/* n1 */,const std::size_t &n2,const std::size_t &n3,const std::size_t &n4,const std::size_t &e,const std::size_t &i,const std::size_t &j,const std::size_t &l){
        return (e*n2*n3*n4 + i*n3*n4 + j*n4 + l);
    }
    
    double Mean (const double *,const std::size_t &); //for array
    double Mean(const std::vector<double> &);
    
    double Variance (const double *,const std::size_t &);
    double Variance(const std::vector<double> &);
    
    double Covariance (const double *,const double *,const std::size_t &);
    double Covariance(const std::vector<double> &,const std::vector<double> &);
    
    ///
    /// convert an array/vector into rank array/vector
    ///
    void Rank(double *,const std::size_t &);
    void Rank(std::vector<double> &);
    
    void Cart2LOS(const double &,const double &,const vec3_t<double> &,vec3_t<double> &);
    void LOS2Cart(const double &,const double &,const vec3_t<double> &,vec3_t<double> &);
    
    ///
    /// converting brightness temp into thermal temp with T_0 = 2.725K, Prog.Theor.Exp.Phys. (2014) 2014 (6): 06B109.
    ///
    inline double temp_convert(const double &temp_br,const double &freq) {
        const double p {CGS_U_h_planck*freq/(CGS_U_kB*2.725)};
        return temp_br*(exp(p)-1.)*(exp(p)-1.)/(p*p*exp(p));;
    }
    ///
    /// convert cartesian coordiante into frame with galacitc warp, while fields built in warpped frame as in ordinary frame
    ///
    vec3_t<double> warp(const vec3_t<double> &);
    
    ///
    /// use given seed number or generate random seed according to thread and clock
    ///
    std::size_t random_seed(const int &);
    
    double timestamp(void);
}

#endif
