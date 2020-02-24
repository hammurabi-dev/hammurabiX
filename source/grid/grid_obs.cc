// observable field grid

#include <array>
#include <cassert>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <grid.h>
#include <hamdis.h>
#include <hamio.h>
#include <hamtype.h>
#include <hamunits.h>
#include <param.h>

// line of sight integrator
Grid_obs::Grid_obs(const Param *par) { build_grid(par); }

void Grid_obs::build_grid(const Param *par) {
  if (par->grid_obs.do_dm) {
    dm_map = std::make_unique<Hampix<ham_float>>(par->grid_obs.nside_dm);
    tmp_dm_map = std::make_unique<Hampix<ham_float>>();
  }
  if (par->grid_obs.do_sync.back()) {
    is_map =
        std::make_unique<Hampix<ham_float>>(par->grid_obs.nside_sync.back());
    tmp_is_map = std::make_unique<Hampix<ham_float>>();
    qs_map =
        std::make_unique<Hampix<ham_float>>(par->grid_obs.nside_sync.back());
    tmp_qs_map = std::make_unique<Hampix<ham_float>>();
    us_map =
        std::make_unique<Hampix<ham_float>>(par->grid_obs.nside_sync.back());
    tmp_us_map = std::make_unique<Hampix<ham_float>>();
    fd_map =
        std::make_unique<Hampix<ham_float>>(par->grid_obs.nside_sync.back());
    tmp_fd_map = std::make_unique<Hampix<ham_float>>();
  }
  if (par->grid_obs.do_fd) {
    fd_map = std::make_unique<Hampix<ham_float>>(par->grid_obs.nside_fd);
    tmp_fd_map = std::make_unique<Hampix<ham_float>>();
  }
  if (par->grid_obs.do_mask) {
    Hamio<ham_float> maskio(par->grid_obs.mask_name);
    auto map_host =
        std::make_unique<Hampix<ham_float>>(par->grid_obs.nside_mask);
    maskio.load(*map_host);
    mask_map = std::make_unique<Hampisk<ham_float>>(*map_host);
  }
}

void Grid_obs::export_grid(const Param *par) {
  Hamio<ham_float> expio;
  if (par->grid_obs.do_dm) {
    // in units pc/cm^3, conventional units
    dm_map->rescale(cgs::ccm / cgs::pc);
    expio.filename(par->grid_obs.sim_dm_name);
    expio.dump(*dm_map);
  }
  if (par->grid_obs.do_sync.back()) {
    std::string name_tag_i = "I_";
    std::string name_tag_q = "Q_";
    std::string name_tag_u = "U_";
    name_tag_i.append(par->grid_obs.sim_sync_name.back());
    name_tag_q.append(par->grid_obs.sim_sync_name.back());
    name_tag_u.append(par->grid_obs.sim_sync_name.back());
    expio.filename(name_tag_i);
    expio.dump(*is_map);
    expio.filename(name_tag_q);
    expio.dump(*qs_map);
    expio.filename(name_tag_u);
    expio.dump(*us_map);
  }
  if (par->grid_obs.do_fd) {
    // FD units is rad*m^(-2) in our calculation
    fd_map->rescale(cgs::m * cgs::m);
    expio.filename(par->grid_obs.sim_fd_name);
    expio.dump(*fd_map);
  }
}
