// allocating physical/observable fields

#ifndef HAMMURABI_GRID_H
#define HAMMURABI_GRID_H

#include <string>
#include <array>
#include <vector>
#include <memory>
#include <cassert>

#include <fftw3.h>
#include <tinyxml2.h>
#include <healpix_map.h>
#include <param.h>

class Grid{
public:
    Grid () = default;
    Grid (const Grid&) = delete;
    Grid (const Grid&&) = delete;
    Grid& operator= (const Grid &) = delete;
    Grid& operator= (Grid &&) = delete;
    virtual ~Grid () = default;
    // build up grid and allocate memory
    virtual void build_grid (const Param *);
    // export grid to file
    virtual void export_grid (const Param *);
    // import file to grid
    virtual void import_grid (const Param *);
};

// regular GMF grid
class Grid_breg final : public Grid{
public:
    Grid_breg () = default;
    Grid_breg (const Param *);
    Grid_breg (const Grid_breg&) = delete;
    Grid_breg (const Grid_breg&&) = delete;
    Grid_breg& operator= (const Grid_breg &) = delete;
    Grid_breg& operator= (Grid_breg &&) = delete;
    virtual ~Grid_breg () = default;
    void build_grid (const Param *) override;
    void export_grid (const Param *) override;
    void import_grid (const Param *) override;
    std::unique_ptr<double[]> bx, by, bz;
};

// turbulent GMF grid
class Grid_brnd final : public Grid{
public:
    Grid_brnd () = default;
    Grid_brnd (const Param *);
    Grid_brnd (const Grid_brnd&) = delete;
    Grid_brnd (const Grid_brnd&&) = delete;
    Grid_brnd& operator= (const Grid_brnd &) = delete;
    Grid_brnd& operator= (Grid_brnd &&) = delete;
    virtual ~Grid_brnd (){
        if (clean_switch){
            fftw_destroy_plan (plan_c0_bw);
            fftw_destroy_plan (plan_c1_bw);
            fftw_destroy_plan (plan_c0_fw);
            fftw_destroy_plan (plan_c1_fw);
            fftw_free (c0);
            fftw_free (c1);
#ifdef _OPENMP
            fftw_cleanup_threads();
#else
            fftw_cleanup();
#endif
        }
    };
    void build_grid (const Param *) override;
    void export_grid (const Param *) override;
    void import_grid (const Param *) override;
    // spatial space
    std::unique_ptr<double[]> bx, by, bz;
    // Fourier space
    fftw_complex *c0, *c1;
    // for/backward plans
    fftw_plan plan_c0_bw, plan_c1_bw, plan_c0_fw, plan_c1_fw;
    // for destructor
    bool clean_switch = false;
};

// regular free electron field grid
class Grid_fereg final : public Grid{
public:
    Grid_fereg () = default;
    Grid_fereg (const Param *);
    Grid_fereg (const Grid_fereg&) = delete;
    Grid_fereg (const Grid_fereg&&) = delete;
    Grid_fereg& operator= (const Grid_fereg &) = delete;
    Grid_fereg& operator= (Grid_fereg &&) = delete;
    virtual ~Grid_fereg () = default;
    void build_grid (const Param *) override;
    void export_grid (const Param *) override;
    void import_grid (const Param *) override;
    std::unique_ptr<double[]> fe;
};

// turbulent free electron field grid
class Grid_fernd final : public Grid{
public:
    Grid_fernd () = default;
    Grid_fernd (const Param *);
    Grid_fernd (const Grid_fernd&) = delete;
    Grid_fernd (const Grid_fernd&&) = delete;
    Grid_fernd& operator= (const Grid_fernd &) = delete;
    Grid_fernd& operator= (Grid_fernd &&) = delete;
    virtual ~Grid_fernd (){
        if (clean_switch){
            fftw_destroy_plan (plan_fe_bw);
            fftw_free (fe_k);
#ifdef _OPENMP
            fftw_cleanup_threads();
#else
            fftw_cleanup();
#endif
        }
    };
    void build_grid (const Param *) override;
    void export_grid (const Param *) override;
    void import_grid (const Param *) override;
    std::unique_ptr<double[]> fe;
    fftw_complex *fe_k;
    fftw_plan plan_fe_bw;
    bool clean_switch = false;
};

// CRE grid
class Grid_cre final : public Grid{
public:
    Grid_cre () = default;
    Grid_cre (const Param *);
    Grid_cre (const Grid_cre&) = delete;
    Grid_cre (const Grid_cre&&) = delete;
    Grid_cre& operator= (const Grid_cre &) = delete;
    Grid_cre& operator= (Grid_cre &&) = delete;
    virtual ~Grid_cre () = default;
    void build_grid (const Param *) override;
    void export_grid (const Param *) override;
    void import_grid (const Param *) override;
    std::unique_ptr<double[]> cre_flux;
};

// observable field grid
class Grid_int final : public Grid{
public:
    Grid_int () = default;
    Grid_int (const Param *);
    Grid_int (const Grid_int&) = delete;
    Grid_int (const Grid_int&&) = delete;
    Grid_int& operator= (const Grid_int &) = delete;
    Grid_int& operator= (Grid_int &&) = delete;
    virtual ~Grid_int () = default;
    void build_grid (const Param *) override;
    void export_grid (const Param *) override;
    std::unique_ptr<Healpix_Map<double>> dm_map, Is_map, Qs_map, Us_map, fd_map;
    // temporary map cache for each shell
    std::unique_ptr<Healpix_Map<double>> tmp_dm_map, tmp_Is_map, tmp_Qs_map, tmp_Us_map, tmp_fd_map;
};

#endif

// END
