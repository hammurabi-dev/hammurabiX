#include <vector>
#include <vec3.h>
#include <cmath>
#include <sys/time.h>
#include <namespace_toolkit.h>
#include <cgs_units_file.h>
#include <thread> // for random seed generator
#include <sstream>
#include <tinyxml2.h>
#include <ap_err.h>
using namespace std;
using namespace tinyxml2;

namespace toolkit {
    // calculate the perpendicular to LOS component of a vector
    double get_perp2LOS (const vec3_t<double> &input,const double &the_los,const double &phi_los){
        const vec3_t<double> unit_vec {get_LOS_unit_vec(the_los, phi_los)};
        const vec3_t<double> perp_vec {crossprod(unit_vec,input)};
        return perp_vec.Length();
    }
    // calculate the parallel to LOS component of a vector
    double get_par2LOS (const vec3_t<double> &input,const double &the_los,const double &phi_los){
        const vec3_t<double> unit_vec {get_LOS_unit_vec(the_los, phi_los)};
        return dotprod(unit_vec,input);
    }
    // calculate intrinsic polarization angle
    double get_intr_pol_ang(const vec3_t<double> &input,const double &the_ec,const double &phi_ec){
        vec3_t<double> sph_unit_v_the;
        vec3_t<double> sph_unit_v_phi;
        sph_unit_v_the = vec3_t<double> {cos(the_ec)*cos(phi_ec),
            cos(the_ec)*sin(phi_ec),
            -sin(the_ec)};
        sph_unit_v_phi = vec3_t<double> {-sin(phi_ec),
            cos(phi_ec),
            0.};
        // IAU convention
        const double y_component {-dotprod(sph_unit_v_the, input)};
        const double x_component {-dotprod(sph_unit_v_phi, input)};
        double result {atan2(y_component, x_component)};
        return result;
    }
    // from Cartesian coordiante to cylindrical coordinate
    void cart_coord2cyl_coord(const vec3_t<double> &input, double &r, double &phi, double &z){
        r = sqrt(input.x*input.x+input.y*input.y);
        phi = atan2(input.y,input.x);
        //if(phi<0.) {phi+=2.*CGS_U_pi;} // may not be necessary
        z = input.z;
    }
    // overload for vec3 to vec3
    void cart_coord2cyl_coord(const vec3_t<double> &input, vec3_t<double> &cyl_vec){
        cyl_vec.x = sqrt(input.x*input.x+input.y*input.y);
        cyl_vec.y = atan2(input.y,input.x);
        //if(cyl_vec.y<0.) {cyl_vec.y+=2.*CGS_U_pi;} // may not be necessary
        cyl_vec.z = input.z;
    }
    // get versor of a vector
    vec3_t<double> versor(const vec3_t<double> &b){
        if(b.Length()==0.) {return vec3_t<double> {0.,0.,0.};}
        vec3_t<double> V = b;
        V.Normalize();
        return V;
    }
    // Mean for array
    double Mean(const double *arr,const std::size_t &size){
        double avg {0};
        for(std::size_t i=0;i!=size;++i) {
            avg += arr[i];
        }
        avg/=size;
        return avg;
    }
    // Mean for vector
    double Mean(const std::vector<double> &vect){
#ifdef DEBUG
        if(vect.empty()){
            ap_err("empty vector");
            exit(1);
        }
#endif
        double avg {0};
        for(auto &i : vect) {
            avg += i;
        }
        avg/=vect.size();
        return avg;
    }
    // Variance for array
    double Variance(const double *arr,const std::size_t &size){
        const double avg {Mean(arr,size)};
        double var {0.};
        for(std::size_t i=0;i!=size;++i){
            var += pow((arr[i]-avg),2.);
        }
        var/=size;
        return var;
    }
    // Variance for vector
    double Variance(const std::vector<double> &vect){
#ifdef DEBUG
        if(vect.empty()){
            ap_err("empty vector");
            exit(1);
        }
#endif
        const double avg {Mean(vect)};
        double var {0.};
        for(auto &i : vect){
            var += pow((i-avg),2.);
        }
        var/=vect.size();
        return var;
    }
    // cov for array
    double Covariance (const double *arr1,const double *arr2,const std::size_t &size){
        double avg1 {Mean(arr1,size)};
        double avg2 {Mean(arr2,size)};
        double covar {0.};
        for(std::size_t m=0;m!=size;++m){
            covar += (arr1[m]-avg1)*(arr2[m]-avg2);
        }
        covar /= size;
        return covar;
    }
    // cov for vector
    double Covariance (const vector<double> &vect1,const vector<double> &vect2){
#ifdef DEBUG
        if(vect1.size()!=vect2.size()){
            ap_err("unequal vector size");
            exit(1);
        }
#endif
        double avg1 {Mean(vect1)};
        double avg2 {Mean(vect2)};
        double covar {0.};
        for(decltype(vect1.size())i=0;i!=vect1.size();++i){
            covar += (vect1[i]-avg1)*(vect2[i]-avg2);
        }
        covar /= vect1.size();
        return covar;
    }
    // get ranked array
    void Rank(double *arr,const size_t &size){
        // get max and min
        double max {arr[0]}; double min {arr[0]};
        for(std::size_t i=0;i!=size;++i){
            if(arr[i]>max) max=arr[i];
            else if(arr[i]<min) min=arr[i];
        }
#ifdef DEBUG
        if(max==min){
            ap_err("invalid array");
            exit(1);
        }
#endif
        // get elements ranked
        for(std::size_t i=0;i!=size;++i){
            arr[i] = (arr[i]-min)/(max-min);
        }
    }
    // get ranked vector
    void Rank(vector<double> &vect){
#ifdef DEBUG
        if(vect.empty()){
            ap_err("empty vector");
            exit(1);
        }
#endif
        // get max and min
        double max=vect[0];double min=vect[0];
        for(auto &i : vect){
            if(i>max) max=i;
            else if(i<min) min=i;
        }
        // get elements ranked
        for(auto &i : vect){
            i = (i-min)/(max-min);
        }
    }
    // galactic warp of cartesian coordinate
    // with this module, we still do modelling in ordinary flat frame.
    // input, Cartesian gc_pos in cgs units, output, warp z in cgs units
    vec3_t<double> warp(const vec3_t<double> &gc_pos){
        vec3_t<double> wpos {gc_pos};
        const double R_w {8.4*CGS_U_kpc};
        const double sr {sqrt( gc_pos.x*gc_pos.x + gc_pos.y*gc_pos.y )};
        if (sr >= R_w) {
            const double gamma_w {0.14};//using cgs units
            const double theta {atan2(gc_pos.y,gc_pos.x)};
            wpos.z -= gamma_w*(sr-R_w)*sin(theta);
        }
        return wpos;
    }
    // offer random seed
    std::size_t random_seed(const int &s){
        if(s<0){
            struct timeval tv;
            gettimeofday(&tv,nullptr);
            // casting thread id into unsinged long
            stringstream ss;
            ss << this_thread::get_id();
            auto th_id = stoul(ss.str());
            return (th_id + tv.tv_sec + tv.tv_usec);
        }
        return s;
    }
    // record time in ms
    double timestamp(void){
        struct timeval tv;
        gettimeofday(&tv,nullptr);
        return tv.tv_sec*1.e+3 + tv.tv_usec*1e-3;
    }
    // auxiliary functions for class Grid and Param
    std::string FetchString(XMLElement* el, string obj){
#ifdef DEBUG
        try{
#endif
            return el->FirstChildElement(obj.c_str())->Attribute("value");
#ifdef DEBUG
        }catch(...){
            ap_err("fail");
            exit(1);
        }
#endif
    }
    
    int FetchInt(XMLElement* el, string obj){
#ifdef DEBUG
        try{
#endif
        return el->FirstChildElement(obj.c_str())->IntAttribute("value");
#ifdef DEBUG
        }catch(...){
            ap_err("fail");
            exit(1);
        }
#endif
    }
    
    unsigned int FetchUnsigned(XMLElement* el, string obj){
#ifdef DEBUG
        try{
#endif
        return el->FirstChildElement(obj.c_str())->UnsignedAttribute("value");
#ifdef DEBUG
        }catch(...){
            ap_err("fail");
            exit(1);
        }
#endif
    }
    
    bool FetchBool(XMLElement* el, string obj){
#ifdef DEBUG
        try{
#endif
        return el->FirstChildElement(obj.c_str())->BoolAttribute("cue");
#ifdef DEBUG
        }catch(...){
            ap_err("fail");
            exit(1);
        }
#endif
    }
    
    double FetchDouble(XMLElement* el, string obj){
#ifdef DEBUG
        try{
#endif
        return el->FirstChildElement(obj.c_str())->DoubleAttribute("value");
#ifdef DEBUG
        }catch(...){
            ap_err("fail");
            exit(1);
        }
#endif
    }
}
// END
