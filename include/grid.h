///
/// allocating physical/observable fields
///
#ifndef GENERIC_GRID_H
#define GENERIC_GRID_H

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
    ///
    /// build up grid and allocate memory
    ///
    virtual void build_grid(XMLDocument *);
    ///
    /// export grid to file
    ///
    virtual void export_grid(void);
    ///
    /// import file to grid
    ///
    virtual void import_grid(void);
    
protected:
    //auxiliary functions
    std::string FetchString(XMLElement *,std::string);
    int FetchInt(XMLElement *,std::string);
    unsigned int FetchUnsigned(XMLElement *,std::string);
    bool FetchBool(XMLElement *,std::string);
    double FetchDouble(XMLElement *,std::string);
};

///
/// regular GMF grid
///
class Grid_breg final : public Grid{
public:
    Grid_breg(std::string);
    virtual ~Grid_breg(void) = default;
    void build_grid(XMLDocument *) override;
    void export_grid(void) override;
    void import_grid(void) override;
    std::unique_ptr<double[]> reg_b_x, reg_b_y, reg_b_z; ///< 3D regular GMF arrays
    
    std::string filename;
    bool read_permission, write_permission;
    
    double x_max, x_min, y_max, y_min, z_max, z_min;
    unsigned int nx, ny, nz;
    std::size_t full_size;
    
};

///
/// turbulent GMF grid
///
class Grid_brnd final : public Grid{
public:
    Grid_brnd(std::string);
    virtual ~Grid_brnd(void) {
        if(build_permission or read_permission){
            fftw_destroy_plan(fftw_px_bw);
            fftw_destroy_plan(fftw_py_bw);
            fftw_destroy_plan(fftw_pz_bw);
            fftw_destroy_plan(fftw_px_fw);
            fftw_destroy_plan(fftw_py_fw);
            fftw_destroy_plan(fftw_pz_fw);
            fftw_free(fftw_b_kx);
            fftw_free(fftw_b_ky);
            fftw_free(fftw_b_kz);
        }
    };
    void build_grid(XMLDocument *) override;
    void export_grid(void) override;
    void import_grid(void) override;
    // spatial space
    std::unique_ptr<double[]> fftw_b_x, fftw_b_y, fftw_b_z;
    // Fourier space
    fftw_complex *fftw_b_kx, *fftw_b_ky, *fftw_b_kz;
    // for/backward plans
    fftw_plan fftw_px_bw, fftw_py_bw, fftw_pz_bw;
    fftw_plan fftw_px_fw, fftw_py_fw, fftw_pz_fw;
    
    std::string filename;
    bool read_permission, write_permission, build_permission;
    
    double x_max, x_min, y_max, y_min, z_max, z_min;
    unsigned int nx, ny, nz;
    std::size_t full_size;
};

///
/// regular free electron field grid
///
class Grid_fereg final : public Grid{
public:
    Grid_fereg(std::string);
    virtual ~Grid_fereg(void) = default;
    void build_grid(XMLDocument *) override;
    void export_grid(void) override;
    void import_grid(void) override;
    
    std::unique_ptr<double[]> fe;
    
    std::string filename;
    bool read_permission, write_permission;
    
    double x_max, x_min, y_max, y_min, z_max, z_min;
    unsigned int nx, ny, nz;
    std::size_t full_size;
};

///
/// turbulent free electron field grid
///
class Grid_fernd final : public Grid{
public:
    Grid_fernd(std::string);
    virtual ~Grid_fernd(void) {
        if(build_permission or read_permission){
            fftw_destroy_plan(fftw_p);
            fftw_free(fftw_fe_k);
        }
    };
    void build_grid(XMLDocument *) override;
    void export_grid(void) override;
    void import_grid(void) override;
    
    std::unique_ptr<double[]> fftw_fe;
    
    fftw_complex *fftw_fe_k;
    fftw_plan fftw_p;
    
    std::string filename;
    bool read_permission, write_permission, build_permission;
    
    double x_max, x_min, y_max, y_min, z_max, z_min;
    unsigned int nx, ny, nz;
    std::size_t full_size;
};

///
/// CRE grid
///
class Grid_cre final : public Grid{
public:
    Grid_cre(std::string);
    virtual ~Grid_cre(void) = default;
    void build_grid(XMLDocument *) override;
    void export_grid(void) override;
    void import_grid(void) override;
    
    std::unique_ptr<double[]> cre_flux;
    
    std::string filename;
    bool read_permission, write_permission;
    
    // 2-D spatial 1-D spectral grid
    unsigned int nE, nr, nz;
    std::size_t cre_size;
    double r_max, z_max, z_min, E_min, E_max, E_fact;
    // or 3+1 dimension grid
    unsigned int nx,ny;
    double x_max, x_min, y_max, y_min;
};

///
/// observable field grid
///
class Grid_int final : public Grid{
public:
    Grid_int(std::string);
    virtual ~Grid_int(void) = default;
    void build_grid(XMLDocument *) override;
    void export_grid(void) override;
    
    Healpix_Map<double> dm_map; ///< dispersion measure
    Healpix_Map<double> Is_map; ///< synchrotron total intensity
    Healpix_Map<double> Qs_map; ///< synchrotron Sotkes Q
    Healpix_Map<double> Us_map; ///< synchrotron Stokes U
    Healpix_Map<double> fd_map; ///< Faraday depth
    
    // shell parameters
    unsigned int nside_sim, npix_sim, total_shell;
    std::vector<unsigned int> nside_shell;
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
