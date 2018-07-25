/**
 * regular GMF generators
 */
#ifndef HAMMURABI_BREG_H
#define HAMMURABI_BREG_H

#include <fftw3.h>
#include <vec3.h>
#include <vector>
#include <param.h>
#include <grid.h>

/**
 * base class of GMF generator,
 * read_grid and write_grid are implemented here
 */
class Breg{
public:
    Breg () = default;
    virtual ~Breg () = default;
    /**
     * fetch regular magnetic field
     * 1st argument: Galactic centric Cartesian frame position
     * 2nd argument: parameter class object
     * 3rd argument: regular GMF grid class object
     * inovke read_grid regardless of field type if permitted, otherwise invoke breg
     */
    virtual vec3_t<double> get_breg (const vec3_t<double> &,
                                     const Param *,
                                     const Grid_breg *);
    /**
     * GMF assembler, specified only in derived class
     * 1st argument: Galactic centric Cartesian frame position
     * 2nd argument: parameter class object
     */
    virtual vec3_t<double> breg (const vec3_t<double> &,
                                 const Param *);
    /**
     * read from field grid with trilinear interpolation
     * 1st argument: Galactic centric Cartesian frame position
     * 2nd argument: regular GMF grid class object
     */
    virtual vec3_t<double> read_grid (const vec3_t<double> &,
                                      const Grid_breg *);
    /**
     * write to field grid
     * 1st argument: parameter class object
     * 2nd argument: regular GMF grid class object
     */
    virtual void write_grid (const Param *,
                             const Grid_breg *);
};

/**
 * designed for verification
 */
class Breg_verify final : public Breg{
public:
    Breg_verify () = default;
    virtual ~Breg_verify () = default;
    vec3_t<double> breg (const vec3_t<double> &,
                         const Param *) override;
};

/**
 * WMAP LSA modeling
 * http://iopscience.iop.org/article/10.1086/513699/meta
 * with errata for GMF modeling
 * https://lambda.gsfc.nasa.gov/product/map/dr2/pub_papers/threeyear/polarization/errata.cfm
 */
class Breg_wmap final : public Breg {
public:
    Breg_wmap () = default;
    virtual ~Breg_wmap () = default;
    vec3_t<double> breg (const vec3_t<double> &,
                         const Param *) override;
};

/**
 * Jaffe modeling
 * https://academic.oup.com/mnras/article/401/2/1013/1150693
 * https://www.aanda.org/articles/aa/abs/2016/12/aa28033-15/aa28033-15.html
 */
class Breg_jaffe final : public Breg{
public:
    Breg_jaffe () = default;
    virtual ~Breg_jaffe  () = default;
    vec3_t<double> breg (const vec3_t<double> &,
                         const Param *) override;
private:
    /**
     * field orientation
     * 1st argument: Galactic centric Cartesian frame position
     * 2nd argument: parameter class object
     */
    vec3_t<double> orientation (const vec3_t<double> &,
                                const Param *);
    /**
     * field amplitude radial scaling
     * 1st argument: Galactic centric Cartesian frame position
     * 2nd argument: parameter class object
     */
    double radial_scaling (const vec3_t<double> &,
                           const Param *);
    /**
     * spiral arm height scaling
     * 1st argument: Galactic centric Cartesian frame position
     * 2nd argument: parameter class object
     * remark: inlined function
     */
    inline double arm_scaling (const vec3_t<double> &pos,
                               const Param *par){
        return 1./(cosh(pos.z/par->breg_jaffe.arm_z0)*cosh(pos.z/par->breg_jaffe.arm_z0));
    }
    /**
     * disk height scaling
     * 1st argument: Galactic centric Cartesian frame position
     * 2nd argument: parameter class object
     * remark: inlined function
     */
    inline double disk_scaling (const vec3_t<double> &pos,
                                const Param *par){
        return 1./(cosh(pos.z/par->breg_jaffe.disk_z0)*cosh(pos.z/par->breg_jaffe.disk_z0));
    }
    /**
     * halo height scaling
     * 1st argument: Galactic centric Cartesian frame position
     * 2nd argument: parameter class object
     * remark: inlined function
     */
    inline double halo_scaling (const vec3_t<double> &pos,
                                const Param *par){
        return 1./(cosh(pos.z/par->breg_jaffe.halo_z0)*cosh(pos.z/par->breg_jaffe.halo_z0));
    }
    /**
     * spiral arm compression factor, for each arm
     * 1st argument: Galactic centric Cartesian frame position
     * 2nd argument: parameter class object
     */
    std::vector<double> arm_compress (const vec3_t<double> &,
                                      const Param *);
    /**
     * spiral arm compression factor for dust, for each arm
     * 1st argument: Galactic centric Cartesian frame position
     * 2nd argument: parameter class object
     */
    std::vector<double> arm_compress_dust (const vec3_t<double> &,
                                           const Param *);
    /**
     * distance to each spiral arm
     * 1st argument: Galactic centric Cartesian frame position
     * 2nd argument: parameter class object
     */
    std::vector<double> dist2arm (const vec3_t<double> &,
                                  const Param *);
};

#endif

// END
