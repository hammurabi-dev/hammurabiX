///
/// regular GMF generators
///
#ifndef HAMMURABI_BREG_H
#define HAMMURABI_BREG_H

#include <fftw3.h>
#include <vec3.h>
#include <vector>
#include "param.h"
#include "grid.h"

///
/// base class of GMF generator,
/// \p read_grid and \p write_grid are implemented here
///
class Breg{
public:
    Breg(void) = default;
    virtual ~Breg(void) = default;
    ///
    /// fetch regular magnetic field,
    /// inovke \p read_grid regardless of field type if permitted, otherwise invoke \p breg
    ///
    virtual vec3_t<double> get_breg(const vec3_t<double> &,Param *,Grid_breg *);
    ///
    /// field assembler, specified only in derived class
    ///
    virtual vec3_t<double> breg(const vec3_t<double> &,Param *);
    ///
    /// read from field grid with trilinear interpolation
    ///
    virtual vec3_t<double> read_grid(const vec3_t<double> &,Grid_breg *);
    ///
    /// write to field grid
    ///
    virtual void write_grid(Param *,Grid_breg *);
};

///
/// designed for verification
///
class Breg_verify final : public Breg{
public:
    Breg_verify(void) = default;
    virtual ~Breg_verify(void) = default;
    vec3_t<double> breg(const vec3_t<double> &,Param *) override;
};

///
/// WMAP LSA modeling
///
class Breg_wmap final : public Breg {
public:
    Breg_wmap(void) = default;
    virtual ~Breg_wmap(void) = default;
    vec3_t<double> breg(const vec3_t<double> &,Param *) override;
};

///
/// Jaffe modeling
///
class Breg_jaffe final : public Breg{
public:
    Breg_jaffe(void) = default;
    virtual ~Breg_jaffe(void) = default;
    vec3_t<double> breg(const vec3_t<double> &,Param *) override;
private:
    // field direction
    vec3_t<double> versor(const vec3_t<double> &,Param *);
    // field amplitude radial scaling
    double radial_scaling(const vec3_t<double> &,Param *);
    // spiral arm height scaling
    inline double arm_scaling(const vec3_t<double> &pos,Param *par){
        return 1./(cosh(pos.z/par->breg_jaffe.arm_z0)*cosh(pos.z/par->breg_jaffe.arm_z0));
    }
    // disk height scaling
    inline double disk_scaling(const vec3_t<double> &pos,Param *par){
        return 1./(cosh(pos.z/par->breg_jaffe.disk_z0)*cosh(pos.z/par->breg_jaffe.disk_z0));
    }
    // halo height scaling
    inline double halo_scaling(const vec3_t<double> &pos,Param *par){
        return 1./(cosh(pos.z/par->breg_jaffe.halo_z0)*cosh(pos.z/par->breg_jaffe.halo_z0));
    }
    // spiral arm compression factor, for each arm
    std::vector<double> arm_compress(const vec3_t<double> &,Param *);
    // spiral arm compression factor for dust, for each arm
    std::vector<double> arm_compress_dust(const vec3_t<double> &,Param *);
    // distance to each spiral arm
    std::vector<double> dist2arm(const vec3_t<double> &,Param *);
};

#endif

// END
