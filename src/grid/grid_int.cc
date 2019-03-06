#include <fstream>
#include <memory>
#include <array>
#include <string>
#include <vector>
#include <cassert>

#include <fitshandle.h>
#include <healpix_map_fitsio.h>
#include <fitsio.h>

#include <grid.h>
#include <cgs_units_file.h>
#include <namespace_toolkit.h>
#include <param.h>

// line of sight integrator
Grid_int::Grid_int (const Param *par){
    build_grid (par);
}

void Grid_int::build_grid (const Param *par){
    if (par->grid_int.do_dm) {
        dm_map = std::make_unique<Healpix_Map<double>>();
        dm_map->SetNside (par->grid_int.nside_dm,
                          RING);
    }
    if (par->grid_int.do_sync.back()) {
        Is_map = std::make_unique<Healpix_Map<double>>();
        Is_map->SetNside (par->grid_int.nside_sync.back(),
                          RING);
        Qs_map = std::make_unique<Healpix_Map<double>>();
        Qs_map->SetNside (par->grid_int.nside_sync.back(),
                          RING);
        Us_map = std::make_unique<Healpix_Map<double>>();
        Us_map->SetNside (par->grid_int.nside_sync.back(),
                          RING);
        fd_map = std::make_unique<Healpix_Map<double>>();
        fd_map->SetNside (par->grid_int.nside_sync.back(),
                          RING);
    }
    if (par->grid_int.do_fd) {
        fd_map = std::make_unique<Healpix_Map<double>>();
        fd_map->SetNside (par->grid_int.nside_fd,
                          RING);
    }
}

void Grid_int::export_grid (const Param *par){
    if (par->grid_int.do_dm){
        // in units pc/cm^3, conventional units
        dm_map->Scale (CGS_U_ccm/CGS_U_pc);
        fitshandle out_dm;
        if (std::ifstream (par->grid_int.sim_dm_name.c_str())){
            remove (par->grid_int.sim_dm_name.c_str());
        }
        out_dm.create (par->grid_int.sim_dm_name);
        write_Healpix_map_to_fits (out_dm,
                                   *dm_map,
                                   PLANCK_FLOAT64);
        out_dm.close();
    }
    if (par->grid_int.do_sync.back()){
        fitshandle out_sc;
        if (std::ifstream (par->grid_int.sim_sync_name.back().c_str())){
            remove (par->grid_int.sim_sync_name.back().c_str());
        }
        out_sc.create (par->grid_int.sim_sync_name.back());
        write_Healpix_map_to_fits (out_sc,
                                   *Is_map,
                                   *Qs_map,
                                   *Us_map,
                                   PLANCK_FLOAT64);
        out_sc.close();
    }
    if (par->grid_int.do_fd){
        // FD units is rad*m^(-2) in our calculation
        fd_map->Scale (CGS_U_m*CGS_U_m);
        fitshandle out_fd;
        if (std::ifstream (par->grid_int.sim_fd_name.c_str())){
            remove (par->grid_int.sim_fd_name.c_str());
        }
        out_fd.create (par->grid_int.sim_fd_name);
        write_Healpix_map_to_fits (out_fd,
                                   *fd_map,
                                   PLANCK_FLOAT64);
        out_fd.close();
    }
}
