///
/// iso/aniso-tropic turbulent Galactic magnetic field generators
///
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

///
/// base class with \p read_grid implemented,
/// \p get_brnd is invoked when no specific derived class object is instantiated
///
class Brnd{
public:
    Brnd(void) = default;
    virtual ~Brnd(void) = default;
    ///
    /// return empty field without invoking \p read_grid
    ///
    virtual vec3_t<double> get_brnd(const vec3_t<double> &,Grid_brnd *);
    ///
    /// read field from grid with trilinear interpolation,
    /// user has to call \p write_grid ahead
    ///
    virtual vec3_t<double> read_grid(const vec3_t<double> &,Grid_brnd *);
    ///
    /// write field to grid (model dependent),
    /// user can \p export_grid to binary file with \p Grid_xxx::export_grid
    /// for dynamic binding only, implemented in derived class
    ///
    virtual void write_grid(Param *,Breg *,Grid_breg *,Grid_brnd *);
};

///
/// global (an)isotropic turbulent GMF,
///
class Brnd_global final : public Brnd{
public:
    Brnd_global(void) = default;
    virtual ~Brnd_global(void) = default;
    vec3_t<double> get_brnd(const vec3_t<double> &,Grid_brnd *) override;
    ///
    /// triple Fourier transform scheme
    ///
    void write_grid(Param *,Breg *,Grid_breg *,Grid_brnd *) override;
    
private:
    ///
    /// isotropic power-spectrum
    ///
    virtual double spec(const double &,Param *);
    ///
    /// field energy density rescaling factor
    ///
    virtual double rescal(const vec3_t<double> &,Param *);
    ///
    /// anisotropy factor at given point
    ///
    double anisotropy(const vec3_t<double> &,vec3_t<double> &,Param *,Breg *,Grid_breg *);
    ///
    /// get real part of each element from a complex array
    ///
    void complex2real(const fftw_complex *,double *,const std::size_t &);
    ///
    /// Gram-Schmidt orthogonalization process
    ///
    vec3_t<double> gramschmidt(const vec3_t<double> &,const vec3_t<double> &);
};

///
/// local anisotropic turbulent GMF, in compressive MHD plasma
///
class Brnd_local final : public Brnd{
public:
    Brnd_local(void) = default;
    virtual ~Brnd_local(void) = default;
    vec3_t<double> get_brnd(const vec3_t<double> &,Grid_brnd *) override;
    ///
    /// vector field decomposition scheme,
    /// using regular GMF at Sun position
    ///
    void write_grid(Param *,Breg *,Grid_breg *,Grid_brnd *) override;
    
private:
    ///
    /// dynamo number
    ///
    double dynamo(const double &,const double &);
    ///
    /// fast mode anisotropic tensor structure
    ///
    double hf(const double &,const double &);
    ///
    /// slow mode anisotropic tensor structure
    ///
    double hs(const double &,const double &);
    ///
    /// fast mode anisotropy power factor
    ///
    double fa(const double &,const double &);
    ///
    /// slow mode anisotropy power factor
    ///
    double fs(const double &,const double &);
    
    ///
    /// cosine alpha, where alpha is pitch angle between wavevector and regular GMF
    ///
    double cosa(const vec3_t<double> &,const vec3_t<double> &);
    ///
    /// direction of Alfven mode
    ///
    vec3_t<double> eplus(const vec3_t<double> &,const vec3_t<double> &);
    ///
    /// direction of slow and fast modes
    ///
    vec3_t<double> eminus(const vec3_t<double> &,const vec3_t<double> &);
    
    ///
    /// isotropic part of power spectrum of Alfvenic mode
    ///
    double speca(const double &,Param *);
    ///
    /// isotropic part of power spectrum of fast mode
    ///
    double specf(const double &,Param *);
    ///
    /// isotropic part of power spectrum of slow mode
    ///
    double specs(const double &,Param *);
    
    void complex2real(const fftw_complex *,double *,const std::size_t &);
};

#endif

// END
