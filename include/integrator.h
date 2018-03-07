///
/// line-of-sight integration calculator
///
#ifndef HAMMURABI_INT_H
#define HAMMURABI_INT_H

#include <omp.h>
#include <string>
#include <healpix_map.h>
#include <memory>

#include <breg.h>
#include <brnd.h>
#include <cre.h>
#include <fereg.h>
#include <fernd.h>
#include <param.h>
#include <cgs_units_file.h>

class Integrator {
public:
    Integrator(void) = default;
    virtual ~Integrator(void) = default;
    ///
    /// assmebling pixels/shells into Healpix map
    ///
    void write_grid(Breg *,Brnd *,FEreg *,FErnd *,CRE *,Grid_breg *,Grid_brnd *,Grid_fereg *,Grid_fernd *,Grid_cre *,Grid_int *,Param *);
#ifndef DEBUG
private:
#endif
    ///
    /// to hold temporary information of observables
    ///
    struct struct_observables{
        double Is, Qs, Us;
        double fd;
        double dm;
        double ff;
    };
    ///
    /// to hold temporary information of spherical shells
    ///
    struct struct_shell{
        std::size_t shell_num;
        double d_start;
        double d_stop;
        double delta_d;
        std::size_t step;
        std::vector<double> dist;
    };
    ///
    /// conduct LOS integration in one pixel at given shell
    ///
    void radial_integration(struct_shell *,pointing &,struct_observables &,Breg *,Brnd *,FEreg *,FErnd *,CRE *,Grid_breg *,Grid_brnd *,Grid_fereg *,Grid_fernd *,Grid_cre *,Grid_int *,Param *);
    // auxiliary functions
    ///
    /// max/min shell radius takes 3 arguments
    /// which are, in order, current shell number, total shell number, maximum shell raidus
    ///
    double get_max_shell_radius(const std::size_t &,const std::size_t &,const double &) const;
    double get_min_shell_radius(const std::size_t &,const std::size_t &,const double &) const;
    ///
    /// general upper boundary check
    /// return false if 1st argument is larger than 2nd
    ///
    inline bool check_simulation_upper_limit(const double &value, const double &limit) const {
        return (value>limit);
    }
    ///
    /// general lower boundary check
    /// return false if 1st argument is smaller than 2nd
    ///
    inline bool check_simulation_lower_limit(const double &value, const double &limit) const {
        return (value<limit);
    }
    ///
    /// assembling shell_ref
    /// assembling shell_ref structure
    /// this part may introduce precision loss
    /// wrap into a subroutine for unit-test
    ///
    void assemble_shell_ref(struct_shell *,const Grid_int *,const std::size_t &);
};
#endif
