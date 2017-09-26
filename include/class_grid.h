/*
 *@file: class_grid.h
 *@author: Jiaxin Wang
 *@email: jiwang@sissa.it
 *@brief: allocating/storing/reading/writing physical fields/observables
 * we use TinyXML2 as parameter parser
 
 *@special note: <bool>permission is cue for allocating and im/exporting grid
 * in random field cases, permission is cue for im/exporting only
 */
#ifndef GENERIC_GRID_H
#define GENERIC_GRID_H

#include <fftw3.h>
#include <string>
#include <array>
#include <tinyxml2.h>
#include <healpix_map.h>

using namespace tinyxml2;

/* base */
class Grid{
    public:
    Grid(void) = default;
    virtual ~Grid(void) = default;
    /*@build_grid
     * build up grid and allocate memory
     */
    virtual void build_grid(XMLDocument *);
    /*@clean_grid
     * clean up grid memory
     */
    virtual void clean_grid(void);
    /*@export_grid
     * export grid to file
     */
    virtual void export_grid(void);
    /*@import_grid
     * import file to grid
     */
    virtual void import_grid(void);
    
    protected:
    //auxiliary functions
    std::string FetchString(XMLElement *,std::string);
    int FetchInt(XMLElement *,std::string);
    unsigned int FetchUnsigned(XMLElement *,std::string);
    bool FetchBool(XMLElement *,std::string);
    double FetchDouble(XMLElement *,std::string);
};

/* regular magnetic field */
class Grid_breg final : public Grid{
    public:
    Grid_breg(std::string);
    void build_grid(XMLDocument *) override;
    void clean_grid(void) override;
    void export_grid(void) override;
    void import_grid(void) override;
    
    // 3D regular field arrays
    double *reg_b_x, *reg_b_y, *reg_b_z;
    
    std::string filename;
    bool read_permission, write_permission, ec_frame;
    
    double lx, ly, lz;
    unsigned int nx, ny, nz;
    unsigned long int full_size;
    
};

/* random magnetic field */
class Grid_brnd final : public Grid{
    public:
    Grid_brnd(std::string);
    void build_grid(XMLDocument *) override;
    void clean_grid(void) override;
    void export_grid(void) override;
    void import_grid(void) override;
    // spatial space
    double *fftw_b_x, *fftw_b_y, *fftw_b_z;
    // Fourier space
    fftw_complex *fftw_b_kx, *fftw_b_ky, *fftw_b_kz;
    // for/backward plans
    fftw_plan fftw_px_bw, fftw_py_bw, fftw_pz_bw;
    fftw_plan fftw_px_fw, fftw_py_fw, fftw_pz_fw;
    
    std::string filename;
    bool read_permission, write_permission, ec_frame;
    
    double lx, ly, lz;
    unsigned int nx, ny, nz;
    unsigned long int full_size;
};

/* free electron density field */
class Grid_fe final : public Grid{
    public:
    Grid_fe(std::string);
    void build_grid(XMLDocument *) override;
    void clean_grid(void) override;
    void export_grid(void) override;
    void import_grid(void) override;
    
    double *fe;
    
    std::string filename;
    bool read_permission, write_permission, ec_frame;
    
    double lx, ly, lz;
    unsigned int nx, ny, nz;
    unsigned long int full_size;
};

/* random free electron density field */
class Grid_fernd final : public Grid{
    public:
    Grid_fernd(std::string);
    void build_grid(XMLDocument *) override;
    void clean_grid(void) override;
    void export_grid(void) override;
    void import_grid(void) override;
    
    double *fftw_fe;
    
    fftw_complex *fftw_fe_k;
    fftw_plan fftw_p;
    
    std::string filename;
    bool read_permission, write_permission, ec_frame;
    
    double lx, ly, lz;
    unsigned int nx, ny, nz;
    unsigned long int full_size;
};

/* cosmic ray electron flux field */
class Grid_cre final : public Grid{
    public:
    Grid_cre(std::string);
    void build_grid(XMLDocument *) override;
    void clean_grid(void) override;
    void import_grid(void) override;
    
    double *cre_flux;
    
    std::string filename;
    bool read_permission, write_permission, ec_frame;
    
    // 2-D spatial 1-D spectral grid
    unsigned int nE, nr, nz;
    unsigned long int cre_size;
    // lr is radius not diameter
    // while lz = zmax - zmin
    double lr,lz, Ekmin, Ekmax, Ekfact;
    
    // or 3+1 dimension grid
    unsigned int nx,ny;
    double lx,ly;
};

/* line of signt integrator */
class Grid_int final : public Grid{
    public:
    Grid_int(std::string);
    void build_grid(XMLDocument *) override;
    void export_grid(void) override;
    
    Healpix_Map<double> dm_map;
    Healpix_Map<double> Is_map;
    Healpix_Map<double> Qs_map;
    Healpix_Map<double> Us_map;
    Healpix_Map<double> fd_map;
    
    // shell parameters
    unsigned int nside_sim, nside_min, npix_sim, total_shell, bin_num;
    // shell boundary
    double gc_r_max, ec_r_max, gc_z_max, lat_lim;
    // switches
    bool do_dm, do_sync, do_fd;
    
    std::string sim_sync_name;
    std::string sim_fd_name;
    std::string sim_dm_name;
};

#endif
// END

