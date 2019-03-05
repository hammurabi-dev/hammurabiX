// free electron density distribution
#ifndef HAMMURABI_FE_H
#define HAMMURABI_FE_H

#include <iostream>
#include <hvec.h>
#include <param.h>
#include <grid.h>

class FEreg{
public:
    FEreg () = default;
    FEreg (const FEreg &) = delete;
    FEreg (FEreg &&) = delete;
    FEreg& operator= (const FEreg &) = delete;
    FEreg& operator= (FEreg &&) = delete;
    virtual ~FEreg () = default;
    // get free electron density,
    // read from grid if granted, otherwise,
    // calculate directly from density function
    virtual double get_density (const hvec<3,double> &,
                                const Param *,
                                const Grid_fereg *) const;
    // read from grid with trilinear interpolation
    virtual double read_grid (const hvec<3,double> &,
                              const Param *,
                              const Grid_fereg *) const;
    // write to grid
    virtual void write_grid (const Param *,
                             Grid_fereg *) const;
    // assemble free electron density at given position
    virtual double density (const hvec<3,double> &,
                            const Param *) const;
    // gaussian blur free electron density (grid elemental scale as FWHM),
    // computing time consuming function based on \p density
    virtual double density_blur (const hvec<3,double> &,
                                 const Param *) const;
};

// uniform field modeling
class FEreg_unif final : public FEreg{
public:
    FEreg_unif () = default;
    FEreg_unif (const FEreg_unif &) = delete;
    FEreg_unif (FEreg_unif &&) = delete;
    FEreg_unif& operator= (const FEreg_unif &) = delete;
    FEreg_unif& operator= (FEreg_unif &&) = delete;
    virtual ~FEreg_unif () = default;
    double density (const hvec<3,double> &,
                    const Param *) const override;
};

// YMW16 modeling (ignore Fermi Bubble due to lack of observation)
class FEreg_ymw16 final : public FEreg{
public:
    FEreg_ymw16 () = default;
    FEreg_ymw16 (const FEreg_ymw16 &) = delete;
    FEreg_ymw16 (FEreg_ymw16 &&) = delete;
    FEreg_ymw16& operator= (const FEreg_ymw16 &) = delete;
    FEreg_ymw16& operator= (FEreg_ymw16 &&) = delete;
    virtual ~FEreg_ymw16 () = default;
    double density (const hvec<3,double> &,
                    const Param *) const override;
#ifdef NDEBUG
private:
#endif
    // thick disk
    double thick (const double &,
                  const double &,
                  const Param *) const;
    // thin disk
    double thin (const double &,
                 const double &,
                 const Param *) const;
    // spiral arms
    double spiral (const double &,
                   const double &,
                   const double &,
                   const double &,
                   const Param *) const;
    // galactic center
    double galcen (const double &,
                   const double &,
                   const double &,
                   const Param *) const;
    // gum nebula
    double gum (const double &,
                const double &,
                const double &,
                const Param *) const;
    // local bubble
    double localbubble (const double &,
                        const double &,
                        const double &,
                        const double &,
                        const double &,
                        const Param *) const;
    // northern polar spurs
    double nps (const double &,
                const double &,
                const double &,
                const Param *) const;
};

#endif

// END
