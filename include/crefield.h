// cosmic ray electron flux field generator/interface
//
// class design:
//
// CREfield
//
// method design:
//
//    (model collection)
//       |   |
//       |   | write_field
//       |   |
//       |   |--> (grid)
//       |          |
//       |          | read_field
//       |----------------------> (grid position/index)
//
// CREmodel
//    |-- CREmodel_xxx
//    |-- CREmodel_yyy
//
// method design:
//
//    (model) -----------> (grid)
//       |      write_grid    |
//       |                    |
//       |                    | read_grid
//       | read_model         |
//       ------------------------> (grid position/index)

#ifndef HAMMURABI_CRE_H
#define HAMMURABI_CRE_H

#include <cassert>
#include <grid.h>
#include <hamtype.h>
#include <hamvec.h>
#include <param.h>
#include <stdexcept>
#include <string>

// base class, all functions are implemented in derived class
class CREmodel {
public:
  CREmodel() = default;
  CREmodel(const CREmodel &) = delete;
  CREmodel(CREmodel &&) = delete;
  CREmodel &operator=(const CREmodel &) = delete;
  CREmodel &operator=(CREmodel &&) = delete;
  virtual ~CREmodel() = default;
  // read CRE flux from grid at given grid index
  // 1st argument: CRE grid index
  // 2nd argument: CRE grid
  virtual ham_float read_grid(const ham_uint &, const Grid_cre *) const;
  // read CRE flux from grid at given spatial position and energy coordinate
  // index
  // for reading CRE flux at arbitrary energy, use base class read_field
  // actual value of E is calculated from {E_index,Ek_min,Ek_fact}
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: index in energy coordinate
  // 3rd argument: parameter set
  // 4th argument: CRE grid
  virtual ham_float read_grid(const Hamvec<3, ham_float> &, const ham_uint &,
                              const Param *, const Grid_cre *) const;
  // read CRE flux from grid at given spatial position and energy
  // for reading CRE flux at arbitrary energy, use base class read_field
  // actual value of E is calculated from {E_index,Ek_min,Ek_fact}
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: CRE energy
  // 3rd argument: parameter set
  // 4th argument: CRE grid
  virtual ham_float read_grid(const Hamvec<3, ham_float> &,
                              const Hamvec<1, ham_float> &, const Param *,
                              const Grid_cre *) const;
  // fill the grid with CRE flux distribution
  // 1st argument: parameter set
  // 2nd argument: CRE grid
  virtual void write_grid(const Param *, Grid_cre *) const;
  // assemble CRE flux at given position
  // 1st argument: galactic centric Cartesian spatial position
  // 2nd argument: CRE energy in GeV
  // 3rd argument: parameter set
  virtual ham_float read_model(const Hamvec<3, ham_float> &, const ham_float &,
                               const Param *) const;
  // flux normalization at given position
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter set
  virtual ham_float flux_norm(const Hamvec<3, ham_float> &,
                              const Param *) const;
  // flux index at given position
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter set
  virtual ham_float flux_idx(const Hamvec<3, ham_float> &, const Param *) const;

protected:
  // spatial CRE flux reprofiling
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter set
  virtual ham_float spatial_profile(const Hamvec<3, ham_float> &,
                                    const Param *) const;
};

// uniform CRE flux
class CREmodel_unif : public CREmodel {
public:
  CREmodel_unif() = default;
  CREmodel_unif(const CREmodel_unif &) = delete;
  CREmodel_unif(CREmodel_unif &&) = delete;
  CREmodel_unif &operator=(const CREmodel_unif &) = delete;
  CREmodel_unif &operator=(CREmodel_unif &&) = delete;
  virtual ~CREmodel_unif() = default;
  ham_float read_model(const Hamvec<3, ham_float> &, const ham_float &,
                       const Param *) const override;
  // flux normalization at given position
  ham_float flux_norm(const Hamvec<3, ham_float> &,
                      const Param *) const override;
  // flux index at given position
  ham_float flux_idx(const Hamvec<3, ham_float> &,
                     const Param *) const override;

protected:
  // spatial CRE flux reprofiling
  ham_float spatial_profile(const Hamvec<3, ham_float> &,
                            const Param *) const override;
};

// analytic CRE flux
class CREmodel_ana : public CREmodel {
public:
  CREmodel_ana() = default;
  CREmodel_ana(const CREmodel_ana &) = delete;
  CREmodel_ana(CREmodel_ana &&) = delete;
  CREmodel_ana &operator=(const CREmodel_ana &) = delete;
  CREmodel_ana &operator=(CREmodel_ana &&) = delete;
  virtual ~CREmodel_ana() = default;
  ham_float read_model(const Hamvec<3, ham_float> &, const ham_float &,
                       const Param *) const override;
  // flux normalization at given position
  ham_float flux_norm(const Hamvec<3, ham_float> &,
                      const Param *) const override;
  // flux index at given position
  ham_float flux_idx(const Hamvec<3, ham_float> &,
                     const Param *) const override;

protected:
  // spatial CRE flux reprofiling
  ham_float spatial_profile(const Hamvec<3, ham_float> &,
                            const Param *) const override;
};

// magnetic field collection class
class CREfield {
public:
  CREfield() = default;
  CREfield(const CREfield &) = delete;
  CREfield(CREfield &&) = delete;
  CREfield &operator=(const CREfield &) = delete;
  CREfield &operator=(CREfield &&) = delete;
  virtual ~CREfield() = default;
  // initialize magnetic field model collection
  // model ordering arranged during parameter parsing
  CREfield(Param *par) {
    if (par->cremodel_list.empty()) {
      model_list.push_back(std::make_unique<CREmodel>());
    } else {
      for (auto &name : par->cremodel_list) {
        if (name == "unif") {
          model_list.push_back(std::make_unique<CREmodel_unif>());
        } else if (name == "pwrlaw") {
          model_list.push_back(std::make_unique<CREmodel_ana>());
        } else {
          throw std::runtime_error("unknown CRE model");
        }
      }
    }
  }
  // read CRE flux at given grid index
  // 1st argument: CRE grid index
  // 2nd argument: CRE grid
  ham_float read_field(const ham_uint &idx, const Grid_cre *grid) const {
    return grid->cre_flux[idx];
  }
  // read CRE flux at given position (without energy interpolation)
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: index in energy coordinate
  // 3rd argument: parameter set
  // 4th argument: CRE grid
  ham_float read_field(const Hamvec<3, ham_float> &pos, const ham_uint &idx,
                       const Param *par, const Grid_cre *grid) const {
    ham_float tmp{0.};
    // get field content from grid
    // model ordering arranged during parameter parsing
    if (par->grid_cre.build_permission or par->grid_cre.read_permission) {
      tmp = model_list[0]->read_grid(pos, idx, par, grid);
    }
    // get field content directly from model
    else {
      throw std::runtime_error("un-implemented");
    }
    return tmp;
  }
  // read CRE flux at given position (with energy interpolation)
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: CRE energy
  // 3rd argument: parameter set
  // 4th argument: CRE grid
  ham_float read_field(const Hamvec<3, ham_float> &pos,
                       const Hamvec<1, ham_float> &eng, const Param *par,
                       const Grid_cre *grid) const {
    ham_float tmp{0.};
    // get field content from grid
    // model ordering arranged during parameter parsing
    if (par->grid_cre.build_permission or par->grid_cre.read_permission) {
      // read_grid defined in base class
      tmp = model_list[0]->read_grid(pos, eng, par, grid);
    }
    // get field content directly from model
    else {
      for (auto &model : model_list) {
        tmp += model->read_model(pos, eng[0], par);
      }
    }
    return tmp;
  }
  // write magnetic field to grid before reading
  void write_field(const Param *par, Grid_cre *grid) {
    // write field content to grid
    // model ordering arranged during parameter parsing
    assert(par->grid_cre.build_permission or par->grid_cre.write_permission or
           par->grid_cre.read_permission);
    for (auto &model : model_list) {
      model->write_grid(par, grid);
    }
  }
  // flux normalization at given position
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter set
  ham_float flux_norm(const Hamvec<3, ham_float> &pos, const Param *par) const {
    assert(model_list.size() == 1);
    return model_list[0]->flux_norm(pos, par);
  }
  // flux index at given position
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter set
  ham_float flux_idx(const Hamvec<3, ham_float> &pos, const Param *par) const {
    assert(model_list.size() == 1);
    return model_list[0]->flux_idx(pos, par);
  }

protected:
  // list of magnetic field model pointers
  std::vector<std::unique_ptr<CREmodel>> model_list;
};

#endif
