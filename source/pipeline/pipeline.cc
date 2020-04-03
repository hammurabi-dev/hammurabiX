#include <bfield.h>
#include <crefield.h>
#include <cstdlib>
#include <grid.h>
#include <hamtype.h>
#include <integrator.h>
#include <iostream>
#include <memory>
#include <param.h>
#include <pipeline.h>
#include <stdexcept>
#include <string>
#include <tefield.h>
#include <toolkit.h>
#include <vector>

// constructor
Pipeline::Pipeline(const std::string &filename) {
  par = std::make_unique<Param>(filename);
}

void Pipeline::assemble_grid() {
  grid_b = std::make_unique<Grid_b>(par.get());
  grid_te = std::make_unique<Grid_te>(par.get());
  grid_cre = std::make_unique<Grid_cre>(par.get());
  grid_obs = std::make_unique<Grid_obs>(par.get());
}

// regular magnetic field
void Pipeline::assemble_b() {
  bfield = std::make_unique<Bfield>(par.get());
  if (par->grid_b.read_permission or par->grid_b.write_permission or
      par->grid_b.build_permission) {
    grid_b->build_grid(par.get());
    // allows adding models to numerical input
    bfield->write_field(par.get(), grid_b.get());
    if (par->grid_b.write_permission)
      grid_b->export_grid(par.get());
  }
}

// regular thermel electron field
void Pipeline::assemble_te() {
  tefield = std::make_unique<TEfield>(par.get());
  if (par->grid_te.read_permission or par->grid_te.write_permission or
      par->grid_te.build_permission) {
    grid_te->build_grid(par.get());
    // allows adding models to numerical input
    tefield->write_field(par.get(), grid_te.get());
    if (par->grid_te.write_permission)
      grid_te->export_grid(par.get());
  }
}

// cre flux field
void Pipeline::assemble_cre() {
  crefield = std::make_unique<CREfield>(par.get());
  if (par->grid_cre.read_permission or par->grid_cre.write_permission or
      par->grid_cre.build_permission) {
    grid_cre->build_grid(par.get());
    // allows adding models to numerical input
    crefield->write_field(par.get(), grid_cre.get());
    if (par->grid_cre.write_permission)
      grid_cre->export_grid(par.get());
  }
}

// LoS integration for observables
void Pipeline::assemble_obs() {
  integrator = std::make_unique<Integrator>();
  if (par->grid_obs.write_permission) {
    const auto repeat = par->grid_obs.do_sync.size();
    for (ham_uint i = 0; i < repeat; ++i) {
      if (i > 0) {
        par->grid_obs.do_dm = false;
        par->grid_obs.do_fd = false;
        // need to rebuild integration grid
        grid_obs->build_grid(par.get());
      }
      // connect fields
      integrator->fields.b = bfield.get();
      integrator->fields.te = tefield.get();
      integrator->fields.cre = crefield.get();
      // connect grids
      integrator->grids.b = grid_b.get();
      integrator->grids.te = grid_te.get();
      integrator->grids.cre = grid_cre.get();
      integrator->grids.obs = grid_obs.get();
      // write and export
      integrator->write_grid(par.get());
      grid_obs->export_grid(par.get());
      // delete obsolete parameters
      if (par->grid_obs.do_sync.back()) {
        par->grid_obs.nside_sync.pop_back();
        par->grid_obs.do_sync.pop_back();
        par->grid_obs.sim_sync_freq.pop_back();
        par->grid_obs.sim_sync_name.pop_back();
      }
    }
  }
}
