/**
 * allocating physical/observable fields
 */
#ifndef HAMMURABI_GRID_H
#define HAMMURABI_GRID_H

#include <fftw3.h>
#include <string>
#include <array>
#include <vector>
#include <tinyxml2.h>
#include <healpix_map.h>
#include <memory>

using namespace tinyxml2;

class Grid{
public:
    Grid(void) = default;
    virtual ~Grid(void) = default;
    /**
     * build up grid and allocate memory
     */
    virtual void build_grid(XMLDocument *);
    /**
     * export grid to file
     */
    virtual void export_grid(void);
    /**
     * import file to grid
     */
    virtual void import_grid(void);
};

/**
 * regular GMF grid
 */
class Grid_breg final : public Grid{
public:
    Grid_breg(const std::string &);
    virtual ~Grid_breg(void) = default;
    void build_grid(XMLDocument *) override;
    void export_grid(void) override;
    void import_grid(void) override;
    std::unique_ptr<double[]> bx, by, bz;
    std::string filename;
    bool read_permission, write_permission;
    double x_max, x_min, y_max, y_min, z_max, z_min;
    std::size_t nx, ny, nz, full_size;
};

/**
 * turbulent GMF grid
 */
class Grid_brnd final : public Grid{
public:
    Grid_brnd(const std::string &);
    virtual ~Grid_brnd(void) {
        if(build_permission or read_permission){
            fftw_destroy_plan(plan_c0_bw);
            fftw_destroy_plan(plan_c1_bw);
            fftw_destroy_plan(plan_c0_fw);
            fftw_destroy_plan(plan_c1_fw);
            fftw_free(c0);
            fftw_free(c1);
#ifdef _OPENMP
            fftw_cleanup_threads();
#else
            fftw_cleanup();
#endif
        }
    };
    void build_grid(XMLDocument *) override;
    void export_grid(void) override;
    void import_grid(void) override;
    // spatial space
    std::unique_ptr<double[]> bx, by, bz;
    // Fourier space
    fftw_complex *c0, *c1;
    // for/backward plans
    fftw_plan plan_c0_bw, plan_c1_bw;
    fftw_plan plan_c0_fw, plan_c1_fw;
    std::string filename;
    bool read_permission, write_permission, build_permission;
    double x_max, x_min, y_max, y_min, z_max, z_min;
    std::size_t nx, ny, nz, full_size;
};

/**
 * regular free electron field grid
 */
class Grid_fereg final : public Grid{
public:
    Grid_fereg(const std::string &);
    virtual ~Grid_fereg(void) = default;
    void build_grid(XMLDocument *) override;
    void export_grid(void) override;
    void import_grid(void) override;
    std::unique_ptr<double[]> fe;
    std::string filename;
    bool read_permission, write_permission;
    double x_max, x_min, y_max, y_min, z_max, z_min;
    std::size_t nx, ny, nz, full_size;
};

/**
 * turbulent free electron field grid
 */
class Grid_fernd final : public Grid{
public:
    Grid_fernd(const std::string &);
    virtual ~Grid_fernd(void) {
        if(build_permission or read_permission){
            fftw_destroy_plan(plan_fe_bw);
            fftw_free(fe_k);
#ifdef _OPENMP
            fftw_cleanup_threads();
#else
            fftw_cleanup();
#endif
        }
    };
    void build_grid(XMLDocument *) override;
    void export_grid(void) override;
    void import_grid(void) override;
    std::unique_ptr<double[]> fe;
    fftw_complex *fe_k;
    fftw_plan plan_fe_bw;
    std::string filename;
    bool read_permission, write_permission, build_permission;
    double x_max, x_min, y_max, y_min, z_max, z_min;
    std::size_t nx, ny, nz, full_size;
};

/**
 * CRE grid
 */
class Grid_cre final : public Grid{
public:
    Grid_cre(const std::string &);
    virtual ~Grid_cre(void) = default;
    void build_grid(XMLDocument *) override;
    void export_grid(void) override;
    void import_grid(void) override;
    std::unique_ptr<double[]> cre_flux;
    std::string filename;
    bool read_permission, write_permission;
    
    std::size_t nE, nz, nx, ny;
    std::size_t cre_size;
    double x_max, x_min, y_max, y_min, z_max, z_min;
    double E_min, E_max, E_fact;
};

/**
 * observable field grid
 */
class Grid_int final : public Grid{
public:
    Grid_int(const std::string &);
    Grid_int(void) = default;
    virtual ~Grid_int(void) = default;
    void build_grid(XMLDocument *) override;
    void export_grid(void) override;
    Healpix_Map<double> dm_map;
    Healpix_Map<double> Is_map;
    Healpix_Map<double> Qs_map;
    Healpix_Map<double> Us_map;
    Healpix_Map<double> fd_map;
    // shell parameters
    std::size_t nside_sim, npix_sim, total_shell;
    std::vector<std::size_t> nside_shell;
    // shell boundary
    double gc_r_max, ec_r_max, gc_z_max, radial_res, lat_lim;
    // switches
    bool do_dm, do_sync, do_fd;
    std::string sim_sync_name;
    std::string sim_fd_name;
    std::string sim_dm_name;
};

#endif

// END
