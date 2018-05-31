#include <vector>
#include <vec3.h>
#include <cmath>
#include <fftw3.h>
#include <ctime>
#include <chrono>
#include <cassert>
#include <memory>
#include <omp.h>
#include <thread> // for random seed generator
#include <sstream>
#include <namespace_toolkit.h>
#include <cgs_units_file.h>
#include <tinyxml2.h>

using namespace tinyxml2;

namespace toolkit {
    // calculate the perpendicular to LOS component of a vector
    double perp2los (const vec3_t<double> &input,const double &the_los,const double &phi_los){
        const vec3_t<double> unit_vec {los_versor(the_los, phi_los)};
        const vec3_t<double> perp_vec {crossprod(unit_vec,input)};
        return perp_vec.Length();
    }
    // calculate the parallel to LOS component of a vector
    double par2los (const vec3_t<double> &input,const double &the_los,const double &phi_los){
        const vec3_t<double> unit_vec {los_versor(the_los, phi_los)};
        return dotprod(unit_vec,input);
    }
    // calculate intrinsic polarization angle
    double intr_pol_ang(const vec3_t<double> &input,const double &the_ec,const double &phi_ec){
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
        return atan2(y_component, x_component);
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
        const double L = 1./sqrt(b.x*b.x+b.y*b.y+b.z*b.z);
        return vec3_t<double> {b.x*L,b.y*L,b.z*L};
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
        assert(!vect.empty());
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
        assert(!vect.empty());
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
    double Covariance (const std::vector<double> &vect1,const std::vector<double> &vect2){
        assert(vect1.size()==vect2.size());
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
    void Rank(double *arr,const std::size_t &size){
        // get max and min
        double max {arr[0]}; double min {arr[0]};
        for(std::size_t i=0;i!=size;++i){
            if(arr[i]>max) max=arr[i];
            else if(arr[i]<min) min=arr[i];
        }
        assert(max!=min);
        // get elements ranked
        for(std::size_t i=0;i!=size;++i){
            arr[i] = (arr[i]-min)/(max-min);
        }
    }
    // get ranked vector
    void Rank(std::vector<double> &vect){
        assert(!vect.empty());
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
        assert(s>=0);
        if(s==0){
            auto p = std::chrono::system_clock::now();
            // valid until 19 January, 2038 03:14:08 UTC
            time_t today_time = std::chrono::system_clock::to_time_t(p);
            // casting thread id into unsinged long
            std::stringstream ss;
            ss << std::this_thread::get_id();
            auto th_id = std::stoul(ss.str());
            // precision in (thread,second)
            return (th_id + today_time);
        }
        return s;
    }
    
    // auxiliary functions for parsing parameters
    std::unique_ptr<XMLDocument> loadxml(const std::string& filename){
        std::unique_ptr<XMLDocument> doc = std::make_unique<XMLDocument>();
        doc->LoadFile(filename.c_str());
        assert(!doc->Error());
        return move(doc);
    }
    //
    XMLElement* tracexml(XMLDocument *doc,const std::vector<std::string>& keychain){
        XMLElement* el {doc->FirstChildElement("root")};
        if(!keychain.empty()){
            for(auto key: keychain){
#ifndef NDEBUG
                std::cout<<"key: "<<key<<std::endl;
#endif
                el = el->FirstChildElement(key.c_str());
            }
        }
        return el;
    }
    //
    std::string FetchString(XMLElement* el,const std::string& att_type,const std::string& key){
#ifndef NDEBUG
        std::cout<<"key: "<<key<<" attrib: "<<att_type<<std::endl;
#endif
        return el->FirstChildElement(key.c_str())->Attribute(att_type.c_str());
    }
    //
    std::string FetchString(XMLElement* el,const std::string& att_type){
#ifndef NDEBUG
        std::cout<<"attrib: "<<att_type<<std::endl;
#endif
        return el->Attribute(att_type.c_str());
    }
    //
    int FetchInt(XMLElement* el,const std::string& att_type,const std::string& key){
#ifndef NDEBUG
        std::cout<<"key: "<<key<<" attrib: "<<att_type<<std::endl;
#endif
        return el->FirstChildElement(key.c_str())->IntAttribute(att_type.c_str());
    }
    //
    int FetchInt(XMLElement* el,const std::string& att_type){
#ifndef NDEBUG
        std::cout<<"attrib: "<<att_type<<std::endl;
#endif
        return el->IntAttribute(att_type.c_str());
    }
    //
    unsigned int FetchUnsigned(XMLElement* el,const std::string& att_type,const std::string& key){
#ifndef NDEBUG
        std::cout<<"key: "<<key<<" attrib: "<<att_type<<std::endl;
#endif
        return el->FirstChildElement(key.c_str())->UnsignedAttribute(att_type.c_str());
    }
    //
    unsigned int FetchUnsigned(XMLElement* el,const std::string& att_type){
#ifndef NDEBUG
        std::cout<<"attrib: "<<att_type<<std::endl;
#endif
        return el->UnsignedAttribute(att_type.c_str());
    }
    //
    bool FetchBool(XMLElement* el,const std::string& att_type,const std::string& key){
#ifndef NDEBUG
        std::cout<<"key: "<<key<<" attrib: "<<att_type<<std::endl;
#endif
        return el->FirstChildElement(key.c_str())->BoolAttribute(att_type.c_str());
    }
    //
    bool FetchBool(XMLElement* el,const std::string& att_type){
#ifndef NDEBUG
        std::cout<<"attrib: "<<att_type<<std::endl;
#endif
        return el->BoolAttribute(att_type.c_str());
    }
    //
    double FetchDouble(XMLElement* el,const std::string& att_type,const std::string& key){
#ifndef NDEBUG
        std::cout<<"key: "<<key<<" attrib: "<<att_type<<std::endl;
#endif
        return el->FirstChildElement(key.c_str())->DoubleAttribute(att_type.c_str());
    }
    //
    double FetchDouble(XMLElement* el,const std::string& att_type){
#ifndef NDEBUG
        std::cout<<"attrib: "<<att_type<<std::endl;
#endif
        return el->DoubleAttribute(att_type.c_str());
    }
    // get real components from fftw_complex arrays
    void complex2real(const fftw_complex *input,double *output,const std::size_t &size){
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        // DO NOT CHANGE SCHEDULE TYPE
        for(std::size_t i=0;i<size;++i){
            output[i] = input[i][0];
        }
    }
    //
    void complex2imag(const fftw_complex *input,double *output,const std::size_t &size){
#ifdef _OPENMP
#pragma omp parallel for schedule(static) // DO NOT CHANGE SCHEDULE TYPE
#endif
        for(std::size_t i=0;i<size;++i){
            output[i] = input[i][1];
        }
    }
    //
    void complex2rni(const fftw_complex *input,double *realout,double *imagout,const std::size_t &size){
#ifdef _OPENMP
#pragma omp parallel for schedule(static) // DO NOT CHANGE SCHEDULE TYPE
#endif
        for(std::size_t i=0;i<size;++i){
            realout[i] = input[i][0];
            imagout[i] = input[i][1];
        }
    }
}// end of namespace
// END
