///
/// numerical tools
///
#ifndef HAMMURABI_TOOLKIT_H
#define HAMMURABI_TOOLKIT_H

#include <vec3.h>
#include <vector>
#include <string>
#include <cmath>
#include <fftw3.h>
#include <cgs_units_file.h>
#include <tinyxml2.h>
using namespace tinyxml2;

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
    /// 1st argument: polar angle (in rad)
    /// 2nd argument: azimuthal angle (in rad)
    ///
    inline vec3_t<double> get_LOS_unit_vec(const double &the_los,const double &phi_los){
        return vec3_t<double> {cos(phi_los)*sin(the_los),
            sin(phi_los)*sin(the_los),
            cos(the_los)};
    }
    ///
    /// convert coordinate from Cartesian to cylindrical frame
    /// 1st argument: coordinate in Cartesian frame
    /// 2nd argument: r in cylindrical frame
    /// 3rd argument: phi in cylindrical frame
    /// 4th argument: z in cylindrical frame
    ///
    void cart_coord2cyl_coord(const vec3_t<double> &, double &, double &, double &);
    ///
    /// convert coordinate from Cartesian to cylindrical frame
    /// 1st argument: coordinate in Cartesian frame
    /// 2nd argument: vec3_t{r,phi,z} in cylindrical frame
    ///
    void cart_coord2cyl_coord(const vec3_t<double> &, vec3_t<double> &);
    ///
    /// get versor of given vector
    /// 1st argument: vector in an orthogonal frame
    /// return: versor of given vector
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
    double Mean (const double *,const std::size_t &);
    double Mean(const std::vector<double> &);
    double Variance (const double *,const std::size_t &);
    double Variance(const std::vector<double> &);
    double Covariance (const double *,const double *,const std::size_t &);
    double Covariance(const std::vector<double> &,const std::vector<double> &);
    ///
    /// convert an array into rank array
    /// 1st argument: array pointer
    /// 2nd argument: array size
    ///
    void Rank(double *,const std::size_t &);
    ///
    /// convert a vector into rank vector
    ///
    void Rank(std::vector<double> &);
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
    ///
    /// simple time stamp
    /// return time in ms
    ///
    double timestamp(void);
    ///
    /// substract real part of a fftw complex array
    /// 1st argument: input fftw complex array
    /// 2nd argument: output double array
    /// 3rd argument: (output) array size
    ///
    void complex2real(const fftw_complex *,double *,const std::size_t &);
    // auxiliary functions for class Grid and Param
    std::string FetchString(XMLElement *,std::string);
    int FetchInt(XMLElement *,std::string);
    unsigned int FetchUnsigned(XMLElement *,std::string);
    bool FetchBool(XMLElement *,std::string);
    double FetchDouble(XMLElement *,std::string);
}

#endif
