// map masking module for hammurabi
// it relies on the HEALPix support
//
// it imports external mask map (at the pivot resolution)
// in the HEALPix form (in RING ordering)
// we need the upgraded or downgraded resolutions
// and use copies of the same mask but at various resolutions
// to meet such application requirement, we define a structure
// with multiple binary arrays defined by HEALPix Nside
//
// following the HEALPix definition,
// we distinguish the opti/passi-mistic mask map downgrading
// passimistic downgrading means if a low-res pixel is not fully masked
// then we set it as masked
// optimastic downgrading means if a low-res pixel is not fully masked
// then we leave it as UN-masked
// we import the mask map as it is (with its resolution)
// during checking masking regions for each LoS shell
// we use optismatic downgrading so we won't miss any information at the boundaries
//
// alternatively, with non-HEALPix pixelization
// we provide a list of pointing directions
// in this case there is no need to introduce any mask
// since the masked "pixel"s should have not been included in the first place

#ifndef HAMMURABI_MASK_H
#define HAMMURABI_MASK_H

#include <memory>
#include <utility>
#include <vector>
#include <map>
#include <fstream>
#include <stdexcept>
#include <cassert>

#include <healpix_map.h>
#include <fitshandle.h>
#include <fitsio.h>
#include <healpix_map_fitsio.h>

#include <param.h>

class hampixma{
public:
  hampixma() = default;
  // disable copy
  hampixma(const hampixma &) = delete;
  hampixma &operator=(const hampixma &) = delete;
  virtual ~hampixma() = default;
  // make a copy at the new resolution
  // from the pivot(input)-resolution mask
  // 1st argument: target Nside
  virtual void duplicate(const std::size_t &);
  // import from HEALPix map array in RING ordering
  // 1st argument: Param object of hammurabi
  virtual void import(const Param *);
  // return the mask info at given resolution and index
  // 1st argument: target mask HEALPix resolution
  // 2nd argument: target map index
  inline double info(const std::size_t &nside, const std::size_t &idx) const{
    return (maps->at(nside))[idx];
  }
#ifdef NDEBUG
protected:
#endif
  // a map collection which can host mask copies at various resolutions
  std::unique_ptr<std::map<std::size_t, Healpix_Map<double>>> maps;
  // pivot mask map resolution
  std::size_t pivot_nside = 0;
};

void hampixma::import(const Param *par){
  // double check if a mask map is required
  if (not par->grid_obs.do_mask) {
    throw std::runtime_error("no masking");
  }
  if (not std::ifstream(par->grid_obs.mask_name.c_str())) {
    throw std::runtime_error("no mask file");
  }
  // read in mask map in HEALPix fits format
  auto tmp = std::make_unique<Healpix_Map<double>>();
  read_Healpix_map_from_fits(par->grid_obs.mask_name, *tmp);
  // initialize the maps
  maps = std::make_unique<std::map<std::size_t, Healpix_Map<double>>> ();
  pivot_nside = tmp->Nside();
  maps->insert({pivot_nside,*tmp});
}

void hampixma::duplicate(const std::size_t &nside){
  auto tmp = std::make_unique<Healpix_Map<double>>();
  tmp->SetNside(nside, RING);
  if ( maps->find(nside) != maps->end() ) return; // avoid existed copy
  if(pivot_nside>nside) {
    tmp->Import_degrade(maps->at(pivot_nside));
    const std::size_t npix {static_cast<std::size_t>(tmp->Npix())};
    for (std::size_t i=0;i<npix;++i) {
      (*tmp)[i] = double((*tmp)[i]>0.0);
    }
  }
  else if (pivot_nside<nside)
    tmp->Import_upgrade(maps->at(pivot_nside));
  maps->insert({nside,*tmp});
}

#endif
