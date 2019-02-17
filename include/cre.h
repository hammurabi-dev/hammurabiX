// synchrotron emissivity calculator

#ifndef HAMMURABI_CRE_H
#define HAMMURABI_CRE_H

#include <hvec.h>

#include <param.h>
#include <grid.h>

// base class, all functions are implemented in derived class
class CRE{
public:
    CRE () = default;
    CRE (const CRE &) = delete;
    CRE (CRE &&) = delete;
    CRE& operator= (const CRE &) = delete;
    CRE& operator= (CRE &&) = delete;
    virtual ~CRE () = default;
    virtual double get_emissivity_t (const hvec<3,double> &,
                                     const Param *,
                                     const Grid_cre *,
                                     const double &) const;
    virtual double get_emissivity_p (const hvec<3,double> &,
                                     const Param *,
                                     const Grid_cre *,
                                     const double &) const;
    // read CRE flux from grid at given position
    // (E_index, sylindrical_r, sylindrical_z) with {r,z} in cgs units,
    // actual value of E is calculated from {E_index,Ek_min,Ek_fact}
    // in get_emissivity automatically select bi/trilinear interpolation
    // according to 2+1/3+1 spatial-spectral CRE flux grid
    virtual double read_grid (const std::size_t &,
                              const hvec<3,double> &,
                              const Param *,
                              const Grid_cre *) const;
    virtual void write_grid (const Param *,
                             Grid_cre *) const;
    // CRE flux at given CRE energy,
    // input CRE energy at CGS units,
    // output at [GeV m^2 s sr]^-1 units
    virtual double flux (const hvec<3,double> &,
                         const Param *,
                         const double &) const;
};

// designed for testing
#ifndef NDEBUG
class CRE_test : public CRE{
public:
    CRE_test () = default;
    CRE_test (const CRE_test &) = delete;
    CRE_test (CRE_test &&) = delete;
    CRE_test& operator= (const CRE_test &) = delete;
    CRE_test& operator= (CRE_test &&) = delete;
    virtual ~CRE_test () = default;
    double get_emissivity_t (const hvec<3,double> &,
                             const Param *,
                             const Grid_cre *,
                             const double &) const override;
    double get_emissivity_p (const hvec<3,double> &,
                             const Param *,
                             const Grid_cre *,
                             const double &) const override;
    double flux (const hvec<3,double> &,
                 const Param *,
                 const double &) const override;
    // flux normalization at given position
    double flux_norm (const hvec<3,double> &,
                      const Param *) const;
    // flux index at given position
    double flux_idx (const hvec<3,double> &,
                     const Param *) const;
    // spatial CRE flux rescaling
    double rescal (const hvec<3,double> &,
                   const Param *) const;
};
#endif

// Analytical CRE modeling
class CRE_ana : public CRE{
public:
    CRE_ana () = default;
    CRE_ana (const CRE_ana &) = delete;
    CRE_ana (CRE_ana &&) = delete;
    CRE_ana& operator= (const CRE_ana &) = delete;
    CRE_ana& operator= (CRE_ana &&) = delete;
    virtual ~CRE_ana () = default;
    double get_emissivity_t (const hvec<3,double> &,
                             const Param *,
                             const Grid_cre *,
                             const double &) const override;
    double get_emissivity_p (const hvec<3,double> &,
                             const Param *,
                             const Grid_cre *,
                             const double &) const override;
    double flux (const hvec<3,double> &,
                 const Param *,
                 const double &) const override;
#ifdef NDEBUG
protected:
#endif
    // flux normalization at given position
    double flux_norm (const hvec<3,double> &,
                      const Param *) const;
    // flux index at given position
    double flux_idx (const hvec<3,double> &,
                     const Param *) const;
    // spatial CRE flux rescaling
    double rescal (const hvec<3,double> &,
                   const Param *) const;
};

// use Numerical CRE flux
class CRE_num final: public CRE{
public:
    CRE_num () = default;
    CRE_num (const CRE_num &) = delete;
    CRE_num (CRE_num &&) = delete;
    CRE_num& operator= (const CRE_num &) = delete;
    CRE_num& operator= (CRE_num &&) = delete;
    virtual ~CRE_num () = default;
    double get_emissivity_t (const hvec<3,double> &,
                             const Param *,
                             const Grid_cre *,
                             const double &) const override;
    double get_emissivity_p (const hvec<3,double> &,
                             const Param *,
                             const Grid_cre *,
                             const double &) const override;
};
#endif

// END
