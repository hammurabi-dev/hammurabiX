#include <cassert>
#include <memory>
#include <omp.h>
#include <thread> // for random seed generator
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <iostream>

#include <hvec.h>
#include <fftw3.h>
#include <tinyxml2.h>

#include <namespace_toolkit.h>
#include <cgs_units_file.h>

namespace toolkit {
    
    // calculate the perpendicular to LOS component of a vector
    double perp2los (const hvec<3,double>& input,
                     const double& the_los,
                     const double& phi_los){
        const hvec<3,double> unit_vec {los_versor(the_los, phi_los)};
        const hvec<3,double> perp_vec {unit_vec.crossprod(input)};
        return perp_vec.length();
    }
    
    // calculate the parallel to LOS component of a vector
    double par2los (const hvec<3,double>& input,
                    const double& the_los,
                    const double& phi_los){
        const hvec<3,double> unit_vec {los_versor(the_los, phi_los)};
        return unit_vec.dotprod (input);
    }
    
    // calculate intrinsic polarization angle
    double intr_pol_ang (const hvec<3,double>& input,
                         const double& the_ec,
                         const double& phi_ec){
        hvec<3,double> sph_unit_v_the;
        hvec<3,double> sph_unit_v_phi;
        sph_unit_v_the = hvec<3,double> {cos(the_ec)*cos(phi_ec),
            cos(the_ec)*sin(phi_ec),
            -sin(the_ec)};
        sph_unit_v_phi = hvec<3,double> {-sin(phi_ec),
            cos(phi_ec),
            0.};
        // IAU convention
        const double y_component {-sph_unit_v_the.dotprod (input)};
        const double x_component {-sph_unit_v_phi.dotprod (input)};
        return atan2(y_component, x_component);
    }
    
    // from Cartesian coordiante to cylindrical coordinate
    void cart_coord2cyl_coord (const hvec<3,double>& input,
                               double& r,
                               double& phi,
                               double& z){
        r = sqrt(input[0]*input[0]+input[1]*input[1]);
        phi = atan2(input[1],input[0]);
        //if(phi<0.) {phi+=2.*CGS_U_pi;} // may not be necessary
        z = input[2];
    }
    
    // overload for vec3 to vec3
    void cart_coord2cyl_coord (const hvec<3,double>& input,
                               hvec<3,double>& cyl_vec){
        cyl_vec[0] = sqrt(input[0]*input[0]+input[1]*input[1]);
        cyl_vec[1] = atan2(input[1],input[0]);
        //if(cyl_vec[1]<0.) {cyl_vec[1]+=2.*CGS_U_pi;} // may not be necessary
        cyl_vec[2] = input[2];
    }
    
    // mean for array
    double mean (const double* arr,
                 const std::size_t& size){
        double avg {0};
        for(std::size_t i=0;i!=size;++i) {
            avg += arr[i];
        }
        avg/=size;
        return avg;
    }
    
    // mean for vector
    double mean (const std::vector<double>& vect){
        assert(!vect.empty());
        double avg {0};
        for(auto& i : vect) {
            avg += i;
        }
        avg/=vect.size();
        return avg;
    }
    
    // variance for array
    double variance (const double* arr,
                     const std::size_t& size){
        const double avg {mean(arr,size)};
        double var {0.};
        for(std::size_t i=0;i!=size;++i){
            var += (arr[i]-avg)*(arr[i]-avg);
        }
        var/=size;
        return var;
    }
    
    // variance for vector
    double variance (const std::vector<double>& vect){
        assert(!vect.empty());
        const double avg {mean(vect)};
        double var {0.};
        for(auto& i : vect){
            var += (i-avg)*(i-avg);
        }
        var/=vect.size();
        return var;
    }
    
    // cov for array
    double covariance (const double* arr1,
                       const double* arr2,
                       const std::size_t& size){
        double avg1 {mean(arr1,size)};
        double avg2 {mean(arr2,size)};
        double covar {0.};
        for(std::size_t m=0;m!=size;++m){
            covar += (arr1[m]-avg1)*(arr2[m]-avg2);
        }
        covar /= size;
        return covar;
    }
    
    // cov for vector
    double covariance (const std::vector<double>& vect1,
                       const std::vector<double>& vect2){
        assert(vect1.size()==vect2.size());
        double avg1 {mean(vect1)};
        double avg2 {mean(vect2)};
        double covar {0.};
        for(decltype(vect1.size())i=0;i!=vect1.size();++i){
            covar += (vect1[i]-avg1)*(vect2[i]-avg2);
        }
        covar /= vect1.size();
        return covar;
    }
    
    // offer random seed
    std::size_t random_seed (const int& s){
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
    std::unique_ptr<tinyxml2::XMLDocument> loadxml(const std::string& filename){
        auto doc = std::make_unique<tinyxml2::XMLDocument> ();
        doc->LoadFile (filename.c_str());
        assert (!doc->Error());
        return std::move (doc);
    }
    
    tinyxml2::XMLElement* tracexml (tinyxml2::XMLDocument* doc,
                                    const std::vector<std::string>& keychain){
        tinyxml2::XMLElement* el {doc->FirstChildElement("root")};
        if (!keychain.empty()){
            for(auto key: keychain){
#ifdef VERBOSE
                std::cout<<"key: "<<key<<std::endl;
#endif
                el = el->FirstChildElement(key.c_str());
            }
        }
        return el;
    }
    
    std::string fetchstring (tinyxml2::XMLElement* el,
                             const std::string& att_type,
                             const std::string& key){
#ifdef VERBOSE
        std::cout<<"key: "<<key<<" attrib: "<<att_type<<std::endl;
#endif
        return el->FirstChildElement(key.c_str())->Attribute(att_type.c_str());
    }
    
    std::string fetchstring (tinyxml2::XMLElement* el,
                             const std::string& att_type){
#ifdef VERBOSE
        std::cout<<"attrib: "<<att_type<<std::endl;
#endif
        return el->Attribute(att_type.c_str());
    }
    
    int fetchint (tinyxml2::XMLElement* el,
                  const std::string& att_type,
                  const std::string& key,
                  const int& dft){
#ifdef VERBOSE
        std::cout<<"key: "<<key<<" attrib: "<<att_type<<std::endl;
#endif
        return el->FirstChildElement(key.c_str())->IntAttribute(att_type.c_str(),dft);
    }
    
    int fetchint (tinyxml2::XMLElement* el,
                  const std::string& att_type,
                  const int& dft){
#ifdef VERBOSE
        std::cout<<"attrib: "<<att_type<<std::endl;
#endif
        return el->IntAttribute(att_type.c_str(),dft);
    }
    
    unsigned int fetchunsigned (tinyxml2::XMLElement* el,
                                const std::string& att_type,
                                const std::string& key,
                                const unsigned& dft){
#ifdef VERBOSE
        std::cout<<"key: "<<key<<" attrib: "<<att_type<<std::endl;
#endif
        return el->FirstChildElement(key.c_str())->UnsignedAttribute(att_type.c_str(),dft);
    }
    
    unsigned int fetchunsigned (tinyxml2::XMLElement* el,
                                const std::string& att_type,
                                const unsigned& dft){
#ifdef VERBOSE
        std::cout<<"attrib: "<<att_type<<std::endl;
#endif
        return el->UnsignedAttribute(att_type.c_str(),dft);
    }
    
    bool fetchbool (tinyxml2::XMLElement* el,
                    const std::string& att_type,
                    const std::string& key,
                    const int& dft){
#ifdef VERBOSE
        std::cout<<"key: "<<key<<" attrib: "<<att_type<<std::endl;
#endif
        return el->FirstChildElement(key.c_str())->BoolAttribute(att_type.c_str(),dft);
    }
    
    bool fetchbool (tinyxml2::XMLElement* el,
                    const std::string& att_type,
                    const int& dft){
#ifdef VERBOSE
        std::cout<<"attrib: "<<att_type<<std::endl;
#endif
        return el->BoolAttribute(att_type.c_str(),dft);
    }
    
    double fetchdouble (tinyxml2::XMLElement* el,
                        const std::string& att_type,
                        const std::string& key,
                        const double& dft){
#ifdef VERBOSE
        std::cout<<"key: "<<key<<" attrib: "<<att_type<<std::endl;
#endif
        return el->FirstChildElement(key.c_str())->DoubleAttribute(att_type.c_str(),dft);
    }
    
    double fetchdouble (tinyxml2::XMLElement* el,
                        const std::string& att_type,
                        const double& dft){
#ifdef VERBOSE
        std::cout<<"attrib: "<<att_type<<std::endl;
#endif
        return el->DoubleAttribute(att_type.c_str(),dft);
    }
    
    // get real components from fftw_complex arrays
    void complex2real (const fftw_complex* input,
                       double* output,
                       const std::size_t& size){
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        // DO NOT CHANGE SCHEDULE TYPE
        for(std::size_t i=0;i<size;++i){
            output[i] = input[i][0];
        }
    }
    
    void complex2imag (const fftw_complex* input,
                       double* output,
                       const std::size_t& size){
#ifdef _OPENMP
#pragma omp parallel for schedule(static) // DO NOT CHANGE SCHEDULE TYPE
#endif
        for(std::size_t i=0;i<size;++i){
            output[i] = input[i][1];
        }
    }
    
    void complex2rni (const fftw_complex* input,
                      double* realout,
                      double* imagout,
                      const std::size_t& size){
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
