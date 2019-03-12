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
        tmp_dm_map = std::make_unique<Healpix_Map<double>>();
        tmp_dm_map->SetNside (par->grid_int.nside_sim,
                              RING);
    }
    if (par->grid_int.do_sync.back()) {
        Is_map = std::make_unique<Healpix_Map<double>>();
        Is_map->SetNside (par->grid_int.nside_sync.back(),
                          RING);
        tmp_Is_map = std::make_unique<Healpix_Map<double>>();
        tmp_Is_map->SetNside (par->grid_int.nside_sim,
                              RING);
        Qs_map = std::make_unique<Healpix_Map<double>>();
        Qs_map->SetNside (par->grid_int.nside_sync.back(),
                          RING);
        tmp_Qs_map = std::make_unique<Healpix_Map<double>>();
        tmp_Qs_map->SetNside (par->grid_int.nside_sim,
                              RING);
        Us_map = std::make_unique<Healpix_Map<double>>();
        Us_map->SetNside (par->grid_int.nside_sync.back(),
                          RING);
        tmp_Us_map = std::make_unique<Healpix_Map<double>>();
        tmp_Us_map->SetNside (par->grid_int.nside_sim,
                              RING);
        fd_map = std::make_unique<Healpix_Map<double>>();
        fd_map->SetNside (par->grid_int.nside_sync.back(),
                          RING);
        tmp_fd_map = std::make_unique<Healpix_Map<double>>();
        tmp_fd_map->SetNside (par->grid_int.nside_sim,
                              RING);
    }
    if (par->grid_int.do_fd) {
        fd_map = std::make_unique<Healpix_Map<double>>();
        fd_map->SetNside (par->grid_int.nside_fd,
                          RING);
        tmp_fd_map = std::make_unique<Healpix_Map<double>>();
        tmp_fd_map->SetNside (par->grid_int.nside_sim,
                              RING);
    }
}

void Grid_int::export_grid (const Param *par){
    if (par->grid_int.do_dm){
        if (par->grid_int.nside_dm > par->grid_int.nside_sim){
            dm_map->Import_upgrade (*tmp_dm_map);
        }
        else if (par->grid_int.nside_dm < par->grid_int.nside_sim){
            dm_map->Import_degrade (*tmp_dm_map,true);
        }
        else{
            dm_map->Import_nograde (*tmp_dm_map);
        }
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
        if (par->grid_int.nside_sync.back() > par->grid_int.nside_sim){
            Is_map->Import_upgrade (*tmp_Is_map);
            Qs_map->Import_upgrade (*tmp_Qs_map);
            Us_map->Import_upgrade (*tmp_Us_map);
        }
        else if (par->grid_int.nside_sync.back() < par->grid_int.nside_sim){
            Is_map->Import_degrade (*tmp_Is_map,true);
            Qs_map->Import_degrade (*tmp_Qs_map,true);
            Us_map->Import_degrade (*tmp_Us_map,true);
        }
        else{
            Is_map->Import_nograde (*tmp_Is_map);
            Qs_map->Import_nograde (*tmp_Qs_map);
            Us_map->Import_nograde (*tmp_Us_map);
        }
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
        if (par->grid_int.nside_fd > par->grid_int.nside_sim){
            fd_map->Import_upgrade (*tmp_fd_map);
        }
        else if (par->grid_int.nside_fd < par->grid_int.nside_sim){
            fd_map->Import_degrade (*tmp_fd_map,true);
        }
        else{
            fd_map->Import_nograde (*tmp_fd_map);
        }
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
