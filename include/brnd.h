/**
 * (an)iso-tropic (Gaussian) turbulent/random GMF generators
 */
#ifndef HAMMURABI_BRND_H
#define HAMMURABI_BRND_H

#include <memory>
#include <fftw3.h>
#include <vec3.h>
#include <vector>
#include <gsl/gsl_integration.h>
#include <param.h>
#include <grid.h>
#include <cgs_units_file.h>
#include <breg.h>
#include <namespace_toolkit.h>

/**
 * base class with read_grid implemented
 * get_brnd is invoked when no specific derived class object is instantiated
 */
class Brnd{
public:
    Brnd(void) = default;
    virtual ~Brnd(void) = default;
    /**
     * get random field vector
     * return zero field when read_grid is not invoked
     * 1st argument: Galactic centric Cartesian frame position
     * 2nd argument: random GMF grid object
     */
    virtual vec3_t<double> get_brnd(const vec3_t<double> &,
                                    Grid_brnd *);
    /**
     * read field from grid with trilinear interpolation
     * user has to call write_grid ahead in main routine
     * 1st argument: Galactic centric Cartesian frame position
     * 2nd argument: random GMF grid object
     */
    virtual vec3_t<double> read_grid(const vec3_t<double> &,
                                     Grid_brnd *);
    /**
     * write field to grid (model dependent)
     * user can export_grid to binary file with Grid_xxx::export_grid
     * for dynamic binding only, implemented in derived class
     * 1st argument: parameter class object
     * 2nd argument: regular GMF class object
     * 3rd argument: regular GMF grid class object
     * 4th argument: random GMF grid class object
     */
    virtual void write_grid(Param *,
                            Breg *,
                            Grid_breg *,
                            Grid_brnd *);
};

/**
 * global (an)isotropic turbulent GMF
 * this class is treated as a covering class for specific methods
 */
class Brnd_global : public Brnd{
public:
    Brnd_global(void) = default;
    virtual ~Brnd_global(void) = default;
};

/**
 * local (an)isotropic turbulent GMF
 * this class is treated as a covering class for specific methods
 */
class Brnd_local : public Brnd{
public:
    Brnd_local(void) = default;
    virtual ~Brnd_local(void) = default;
};

/**
 * Ensslin-Steininger method of global (an)isotropic turbulent GMF
 */
class Brnd_es final : public Brnd_global{
public:
    Brnd_es(void) = default;
    virtual ~Brnd_es(void) = default;
    /**
     * use triple Fourier transform scheme
     * check technical report for details
     */
    void write_grid(Param *,
                    Breg *,
                    Grid_breg *,
                    Grid_brnd *) override;
private:
    /**
     * isotropic power-spectrum
     * 1st argument: isotropic wave-vector magnitude
     * 2nd argument: parameter class object
     */
    virtual double spec(const double &,
                        Param *);
    /**
     * field energy density rescaling factor
     * 1st argument: Galactic centric Cartesian frame position
     * 2nd argument: parameter class object
     */
    virtual double rescal(const vec3_t<double> &,
                          Param *);
    /**
     * anisotropy factor, and anisotropy direction
     * check technical report for details
     * 1st argument: Galactic centric Cartesian frame position
     * 2nd argument: (in/output) anisotropy direction
     * 3rd argument: parameter class object
     * 4th argument: regular GMF class object
     * 5th argument: regular GMF grid class object
     */
    double anisotropy(const vec3_t<double> &,
                      vec3_t<double> &,
                      Param *,
                      Breg *,
                      Grid_breg *);
    /**
     * Gram-Schmidt orthogonalization process
     * 1st argument: wave-vector
     * 2nd arugment: input GMF vector (in Fourier space)
     * remark: real and imagine parts of complex GMF vector in Fourier space handled separately
     */
    vec3_t<double> gramschmidt(const vec3_t<double> &,
                               const vec3_t<double> &);
};

/**
 * Jaffe method of global (an)isotropic turbulent GMF
 */
//class Brnd_Jaffe final : public Brnd_global

/**
 * local anisotropic turbulent GMF
 * in compressive MHD plasma
 * check technical report for details
 */
class Brnd_mhd final : public Brnd_local{
public:
    Brnd_mhd(void) = default;
    virtual ~Brnd_mhd(void) = default;
    /**
     * use vector field decomposition scheme
     * use regular GMF at position of the Sun
     */
    void write_grid(Param *,
                    Breg *,
                    Grid_breg *,
                    Grid_brnd *) override;
    
private:
    /**
     * dynamo number
     * 1st argument: plasma beta
     * 2nd argument: cosine of k-B pitch angle
     */
    inline double dynamo(const double &beta,
                         const double &cosa){
        return (1+0.5*beta)*(1+0.5*beta) - 2.*beta*cosa*cosa;
    }
    /**
     * fast mode anisotropic tensor structure
     * 1st argument: plasma beta
     * 2nd argument: cosine of k-B pitch angle
     */
    double hf(const double &,
              const double &);
    /**
     * slow mode anisotropic tensor structure
     * 1st argument: plasma beta
     * 2nd argument: cosine of k-B pitch angle
     */
    double hs(const double &,
              const double &);
    /**
     * fast mode anisotropy power factor
     * 1st argument: plasma March number
     * 2nd argument: cosine of k-B pitch angle
     */
    inline double fa(const double &ma,
                     const double &cosa){
        return exp( -pow(ma,-1.33333333)*cosa*cosa/pow(1-cosa*cosa,0.66666667) );
    }
    /**
     * slow mode anisotropy power factor
     * 1st argument: plasma March number
     * 2nd argument: cosine of k-B pitch angle
     */
    inline double fs(const double &ma,
                     const double &cosa){
        return exp( -pow(ma,-1.33333333)*cosa*cosa/pow(1-cosa*cosa,0.66666667) );
    }
    
    /**
     * cosine of pitch angle between wavevector and regular GMF
     * 1st argument: field vector
     * 2nd argument: wave vector
     */
    inline double cosa(const vec3_t<double> &b,
                       const vec3_t<double> &k){
        return dotprod(toolkit::versor(b),toolkit::versor(k));
    }
    /**
     * direction of Alfven mode
     * 1st argument: regualr GMF vector
     * 2nd argument: wave vector
     */
    vec3_t<double> eplus(const vec3_t<double> &,
                         const vec3_t<double> &);
    /**
     * direction of slow and fast modes
     * 1st argument: regualr GMF vector
     * 2nd argument: wave vector
     */
    vec3_t<double> eminus(const vec3_t<double> &,
                          const vec3_t<double> &);
    
    /**
     * isotropic part of power spectrum of Alfvenic mode
     * 1st argument: wave vector magnitude
     * 2nd argument: parameter class object
     */
    double speca(const double &,
                 Param *);
    /**
     * isotropic part of power spectrum of fast mode
     * 1st argument: wave vector magnitude
     * 2nd argument: parameter class object
     */
    double specf(const double &,
                 Param *);
    /**
     * isotropic part of power spectrum of slow mode
     * 1st argument: wave vector magnitude
     * 2nd argument: parameter class object
     */
    double specs(const double &,
                 Param *);
};

#endif

// END
