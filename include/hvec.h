// multi dimensional vector class based on std::vector

#ifndef HAMMURABI_VECTOR_H
#define HAMMURABI_VECTOR_H

#include <iostream>
#include <cmath>
#include <vector>
#include <cassert>

template<int dim,typename T>
class hvec{
protected:
    std::vector<T> ele;
public:
    // default ctor
    hvec () {
        switch(dim){
            case 1: this->ele = {T(0)}; break;
            case 2: this->ele = {T(0),T(0)}; break;
            case 3: this->ele = {T(0),T(0),T(0)}; break;
            default: std::cerr<<"unsupported dimension"; break;
        }
    }
    virtual ~hvec () = default;
    // 1D
    hvec<dim,T> (const T& x){
        assert (dim==1);
        this->ele.push_back (x);
    }
    // 2D
    hvec<dim,T> (const T& x,
                 const T& y){
        assert (dim==2);
        this->ele.push_back (x);
        this->ele.push_back (y);
    }
    // 3D
    hvec<dim,T> (const T& x,
                 const T& y,
                 const T& z){
        assert (dim==3);
        this->ele.push_back (x);
        this->ele.push_back (y);
        this->ele.push_back (z);
    }
    // push back
    void push_back (const T& s){
        this->ele.push_back (s);
    }
    // copy ctor
    hvec<dim,T> (const hvec<dim,T>& v){
        this->ele = v.content();
    }
    // move ctor
    hvec<dim,T> (hvec<dim,T>&& v)
    : ele(std::move(v.content())){
    }
    // copy assign
    hvec<dim,T>& operator= (const hvec<dim,T>& v) noexcept{
        this->ele = std::move(v.content());
        return *this;
    }
    // move assign
    hvec<dim,T>& operator= (hvec<dim,T>&& v) noexcept{
        this->ele = std::move(v.content());
        return *this;
    }
    // operator []
    T operator[] (const int& i) const{
        return this->ele[i];
    }
    // operator []
    T& operator[] (const int& i){
        return this->ele[i];
    }
    // get std::vector<T> ele
    const std::vector<T> content () const{
        return this->ele;
    }
    // get std::vector<T> ele
    std::vector<T>& content (){
        return this->ele;
    }
    // operator +
    template<typename R>
    hvec<dim,T> operator+ (const hvec<dim,R>& v) const{
        hvec<dim,T> tmp(*this);
        for (unsigned int i=0;i<dim;++i){
            tmp[i] += static_cast<T>(v[i]);
        }
        return tmp;
    }
    // operator +=
    // works with different vector data types
    // play with caution
    template<typename R>
    hvec<dim,T>& operator+= (const hvec<dim,R>& v){
        for (unsigned int i=0;i<dim;++i){
            this->ele[i] += static_cast<T>(v[i]);
        }
        return *this;
    }
    // operator -
    // works with different vector data types
    // play with caution
    template<typename R>
    hvec<dim,T> operator- (const hvec<dim,R>& v) const{
        hvec<dim,T> tmp(*this);
        for (unsigned int i=0;i<dim;++i){
            tmp[i] -= static_cast<T>(v[i]);
        }
        return tmp;
    }
    // operator -=
    // works with different vector data types
    // play with caution
    template<typename R>
    hvec<dim,T>& operator-= (const hvec<dim,R>& v){
        for (unsigned int i=0;i<dim;++i){
            this->ele[i] -= static_cast<T>(v[i]);
        }
        return *this;
    }
    // operator *
    // works with different vector data types
    // play with caution
    template<typename R>
    hvec<dim,T> operator* (const R& s) const{
        hvec<dim,T> tmp(*this);
        for (unsigned int i=0;i<dim;++i){
            tmp[i] *= static_cast<T>(s);
        }
        return tmp;
    }
    // operaotr *=
    // works with different vector data types
    // play with caution
    template<typename R>
    hvec<dim,T>& operator*= (const R& s){
        for (unsigned int i=0;i<dim;++i){
            this->ele[i] *= static_cast<T>(s);
        }
        return *this;
    }
    // operator /
    // works with different vector data types
    // play with caution
    template<typename R>
    hvec<dim,T> operator/ (const R& s) const{
        hvec<dim,T> tmp(*this);
        for (unsigned int i=0;i<dim;++i){
            tmp[i] /= static_cast<T>(s);
        }
        return tmp;
    }
    // operator /=
    // works with different vector data types
    // play with caution
    template<typename R>
    hvec<dim,T>& operator/= (const R& s){
        assert(s!=0);
        for (unsigned int i=0;i<dim;++i){
            this->ele[i] /= static_cast<T>(s);
        }
        return *this;
    }
    // operator ==
    bool operator== (const hvec<dim,T>& v){
        for (unsigned int i=0;i<dim;++i){
            if (this->ele[i]!=v[i]){
                return false;
            }
        }
        return true;
    }
    // operator !=
    bool operator!= (const hvec<dim,T>& v){
        for (unsigned int i=0;i<dim;++i){
            if (this->ele[i]!=v[i]){
                return true;
            }
        }
        return false;
    }
    // vector length
    double length () const{
        double tmp {0};
        for (unsigned int i=0;i<dim;++i){
            const double cache = static_cast<double>(this->ele[i]);
            tmp += cache*cache;
        }
        return std::sqrt(tmp);
    }
    // vector squared length
    double lengthsq () const{
        double tmp {0};
        for (unsigned int i=0;i<dim;++i){
            const double cache = static_cast<double>(this->ele[i]);
            tmp += cache*cache;
        }
        return tmp;
    }
    // flip sign
    void flip (){
        for(unsigned int i=0;i<dim;++i){
            this->ele[i] *= static_cast<T>(-1.0);
        }
    }
    // versor
    hvec<dim,double> versor () const{
        hvec<dim,double> tmp;
        for (unsigned int i=0;i<dim;++i)
            tmp[i] = static_cast<double>(this->ele[i]);
        if (tmp.lengthsq() == 0)
            return tmp;
        tmp /= this->length();
        return tmp;
    }
    // inner product
    // works with different vector data type
    // play with caution
    template<typename R>
    double dotprod (const hvec<dim,R>& v) const{
        double tmp {0};
        for (unsigned int i=0;i<dim;++i){
            tmp += this->ele[i]*static_cast<T>(v[i]);
        }
        return tmp;
    }
    // cross product, works in 3D only
    // works with different vector data type
    // play with caution
    template<typename R>
    hvec<dim,double> crossprod (const hvec<dim,R>& v) const{
        assert (dim==3);
        hvec<dim,double> tmp;
        for (unsigned int i=0;i<dim;++i){
            tmp[i] = (static_cast<double>(this->ele[(i+1)%3])*static_cast<double>(v[(i+2)%3])
                      - static_cast<double>(this->ele[(i+2)%3])*static_cast<double>(v[(i+1)%3]));
        }
        return tmp;
    }
    // osteam function
    friend std::ostream& operator<< (std::ostream& os,
                                     const hvec<dim,T>& v){
        for (unsigned int i=0;i<dim;++i){
            os << v[i] << "\t";
        }
        os << std::endl;
        return os;
    }
};

#endif
