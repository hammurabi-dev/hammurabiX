// observable field grid

#include <array>
#include <cassert>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <fitshandle.h>
#include <fitsio.h>
#include <healpix_map_fitsio.h>

#include <cgs_units.h>
#include <grid.h>
#include <param.h>

// line of sight integrator
Grid_obs::Grid_obs(const Param *par) { build_grid(par); }

void Grid_obs::build_grid(const Param *par) {
  if (par->grid_obs.do_dm) {
    dm_map = std::make_unique<Healpix_Map<double>>();
    dm_map->SetNside(par->grid_obs.nside_dm, RING);
    tmp_dm_map = std::make_unique<Healpix_Map<double>>();
  }
  if (par->grid_obs.do_sync.back()) {
    is_map = std::make_unique<Healpix_Map<double>>();
    is_map->SetNside(par->grid_obs.nside_sync.back(), RING);
    tmp_is_map = std::make_unique<Healpix_Map<double>>();
    qs_map = std::make_unique<Healpix_Map<double>>();
    qs_map->SetNside(par->grid_obs.nside_sync.back(), RING);
    tmp_qs_map = std::make_unique<Healpix_Map<double>>();
    us_map = std::make_unique<Healpix_Map<double>>();
    us_map->SetNside(par->grid_obs.nside_sync.back(), RING);
    tmp_us_map = std::make_unique<Healpix_Map<double>>();
    fd_map = std::make_unique<Healpix_Map<double>>();
    fd_map->SetNside(par->grid_obs.nside_sync.back(), RING);
    tmp_fd_map = std::make_unique<Healpix_Map<double>>();
  }
  if (par->grid_obs.do_fd) {
    fd_map = std::make_unique<Healpix_Map<double>>();
    fd_map->SetNside(par->grid_obs.nside_fd, RING);
    tmp_fd_map = std::make_unique<Healpix_Map<double>>();
  }
}

void Grid_obs::export_grid(const Param *par) {
  if (par->grid_obs.do_dm) {
    // in units pc/cm^3, conventional units
    dm_map->Scale(cgs_ccm / cgs_pc);
    fitshandle out_dm;
    if (std::ifstream(par->grid_obs.sim_dm_name.c_str())) {
      remove(par->grid_obs.sim_dm_name.c_str());
    }
    out_dm.create(par->grid_obs.sim_dm_name);
    write_Healpix_map_to_fits(out_dm, *dm_map, PLANCK_FLOAT64);
    out_dm.close();
  }
  if (par->grid_obs.do_sync.back()) {
    fitshandle out_sc;
    if (std::ifstream(par->grid_obs.sim_sync_name.back().c_str())) {
      remove(par->grid_obs.sim_sync_name.back().c_str());
    }
    out_sc.create(par->grid_obs.sim_sync_name.back());
    write_Healpix_map_to_fits(out_sc, *is_map, *qs_map, *us_map,
                              PLANCK_FLOAT64);
    out_sc.close();
  }
  if (par->grid_obs.do_fd) {
    // FD units is rad*m^(-2) in our calculation
    fd_map->Scale(cgs_m * cgs_m);
    fitshandle out_fd;
    if (std::ifstream(par->grid_obs.sim_fd_name.c_str())) {
      remove(par->grid_obs.sim_fd_name.c_str());
    }
    out_fd.create(par->grid_obs.sim_fd_name);
    write_Healpix_map_to_fits(out_fd, *fd_map, PLANCK_FLOAT64);
    out_fd.close();
  }
}
