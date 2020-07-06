// grid for allocating physical/observable fields
//
// class design:
//
// Grid
//  |-- Grid_b (magnetic field grid)
//  |-- Grid_te (thermal electron field grid)
//  |-- Grid_cre (cosmic ray electron flux field grid)
//  |-- Grid_obs (observable grid)
//
// There are three major functions in Grid class
// ``build_grid`` is in charge of preparing data structure and memory allocation
// ``import_grid`` and ``export_grid`` interface with external data sotrage

#ifndef HAMMURABI_GRID_H
#define HAMMURABI_GRID_H

#include <array>
#include <cassert>
#include <fftw3.h>
#include <hamdis.h>
#include <hamsk.h>
#include <hamtype.h>
#include <memory>
#include <param.h>
#include <string>
#include <tinyxml2.h>
#include <vector>

class Grid {
public:
  Grid() = default;
  Grid(const Grid &) = delete;
  Grid(const Grid &&) = delete;
  Grid &operator=(const Grid &) = delete;
  Grid &operator=(Grid &&) = delete;
  virtual ~Grid() = default;
  // build up grid and allocate memory
  virtual void build_grid(const Param *);
  // export grid to file
  virtual void export_grid(const Param *);
  // import file to grid
  virtual void import_grid(const Param *);
};

// regular magnetic vector field grid
class Grid_b final : public Grid {
public:
  Grid_b() = default;
  Grid_b(const Param *);
  Grid_b(const Grid_b &) = delete;
  Grid_b(const Grid_b &&) = delete;
  Grid_b &operator=(const Grid_b &) = delete;
  Grid_b &operator=(Grid_b &&) = delete;
  virtual ~Grid_b() {
    if (clean_fft) {
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
  void build_grid(const Param *) override;
  void export_grid(const Param *) override;
  void import_grid(const Param *) override;
  std::vector<ham_float> *export_grid();
  void import_grid(const std::vector<ham_float> *);
  // spatial domain magnetic field
  std::unique_ptr<std::vector<ham_float>> bx, by, bz;
  // Fourier domain magnetic field
  fftw_complex *c0, *c1;
  // for/backward FFT plans
  fftw_plan plan_c0_bw, plan_c1_bw, plan_c0_fw, plan_c1_fw;
  // for destructor
  bool clean_fft = false;
};

// regular thermal electron density field grid
class Grid_te final : public Grid {
public:
  Grid_te() = default;
  Grid_te(const Param *);
  Grid_te(const Grid_te &) = delete;
  Grid_te(const Grid_te &&) = delete;
  Grid_te &operator=(const Grid_te &) = delete;
  Grid_te &operator=(Grid_te &&) = delete;
  virtual ~Grid_te() {
    if (clean_fft) {
      fftw_destroy_plan(plan_te_bw);
      fftw_free(te_k);
#ifdef _OPENMP
      fftw_cleanup_threads();
#else
      fftw_cleanup();
#endif
    }
  };
  void build_grid(const Param *) override;
  void export_grid(const Param *) override;
  void import_grid(const Param *) override;
  // spatial domain thermal electron field
  std::unique_ptr<ham_float[]> te;
  // Fourier domain thermal electron field
  fftw_complex *te_k;
  // backward FFT plan
  fftw_plan plan_te_bw;
  // for destructor
  bool clean_fft = false;
};

// cosmic ray electron flux phase-space density grid
class Grid_cre final : public Grid {
public:
  Grid_cre() = default;
  Grid_cre(const Param *);
  Grid_cre(const Grid_cre &) = delete;
  Grid_cre(const Grid_cre &&) = delete;
  Grid_cre &operator=(const Grid_cre &) = delete;
  Grid_cre &operator=(Grid_cre &&) = delete;
  virtual ~Grid_cre() = default;
  void build_grid(const Param *) override;
  void export_grid(const Param *) override;
  void import_grid(const Param *) override;
  // phase-space domain CRE flux field
  std::unique_ptr<ham_float[]> cre_flux;
};

// observable field grid
class Grid_obs final : public Grid {
public:
  Grid_obs() = default;
  Grid_obs(const Param *);
  Grid_obs(const Grid_obs &) = delete;
  Grid_obs(const Grid_obs &&) = delete;
  Grid_obs &operator=(const Grid_obs &) = delete;
  Grid_obs &operator=(Grid_obs &&) = delete;
  virtual ~Grid_obs() = default;
  void build_grid(const Param *) override;
  void export_grid(const Param *) override;
  // Hampix map for observables
  // dm_map: dispersion measure
  // is_map: synchrotron Stokes I
  // qs_map: synchrotron Stokes Q
  // us_map: synchrotron Stokes U
  // fd_map: Faraday depth
  std::unique_ptr<Hampix<ham_float>> dm_map, is_map, qs_map, us_map, fd_map;
  // temporary HEALPix map for each shell
  std::unique_ptr<Hampix<ham_float>> tmp_dm_map, tmp_is_map, tmp_qs_map,
      tmp_us_map, tmp_fd_map;
  // mask map
  std::unique_ptr<Hampisk<ham_float>> mask_map;
};

#endif
