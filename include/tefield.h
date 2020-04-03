// thermal electron field generators/interface
//
// class design:
//
// TEfield
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
//
// TEmodel
//   |-- TEmodel_xxx
//   |-- TEmodel_yyy
//
//
// method design:
//
//    (model) -----------> (grid)
//       |      write_grid    |
//       |                    |
//       |                    | read_grid
//       | read_model         |
//       ------------------------> (grid position/index)

#ifndef HAMMURABI_TE_H
#define HAMMURABI_TE_H

#include <cassert>
#include <grid.h>
#include <hamtype.h>
#include <hamvec.h>
#include <param.h>
#include <stdexcept>
#include <string>

class TEmodel {
public:
  TEmodel() = default;
  TEmodel(const TEmodel &) = delete;
  TEmodel(TEmodel &&) = delete;
  TEmodel &operator=(const TEmodel &) = delete;
  TEmodel &operator=(TEmodel &&) = delete;
  virtual ~TEmodel() = default;
  // assemble thermal electron field
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter set
  virtual ham_float read_model(const Hamvec<3, ham_float> &,
                               const Param *) const;
  // read from grid with given grid index
  // 1st argument: thermal electron field grid index
  // 2nd argument: thermal electron field grid
  virtual ham_float read_grid(const ham_uint &, const Grid_te *) const;
  // read from grid with linear interpolation
  // 1st argument: galactic centric Cartesian frame position
  // 2nd argument: parameter set
  // 3rd argument: electron field grid
  virtual ham_float read_grid(const Hamvec<3, ham_float> &, const Param *,
                              const Grid_te *) const;
  // write thermal electron field to grid
  // 1st argument: parameter set
  // 2nd argument: thermal electron field grid
  virtual void write_grid(const Param *, Grid_te *) const;
};

// uniform field modeling
class TEmodel_unif final : public TEmodel {
public:
  TEmodel_unif() = default;
  TEmodel_unif(const TEmodel_unif &) = delete;
  TEmodel_unif(TEmodel_unif &&) = delete;
  TEmodel_unif &operator=(const TEmodel_unif &) = delete;
  TEmodel_unif &operator=(TEmodel_unif &&) = delete;
  virtual ~TEmodel_unif() = default;
  ham_float read_model(const Hamvec<3, ham_float> &,
                       const Param *) const override;
};

// YMW16 modeling (ignore Fermi Bubble due to lack of observation)
class TEmodel_ymw16 final : public TEmodel {
public:
  TEmodel_ymw16() = default;
  TEmodel_ymw16(const TEmodel_ymw16 &) = delete;
  TEmodel_ymw16(TEmodel_ymw16 &&) = delete;
  TEmodel_ymw16 &operator=(const TEmodel_ymw16 &) = delete;
  TEmodel_ymw16 &operator=(TEmodel_ymw16 &&) = delete;
  virtual ~TEmodel_ymw16() = default;
  ham_float read_model(const Hamvec<3, ham_float> &,
                       const Param *) const override;
#ifdef NDEBUG
private:
#endif
  // thick disk
  ham_float thick(const ham_float &, const ham_float &, const Param *) const;
  // thin disk
  ham_float thin(const ham_float &, const ham_float &, const Param *) const;
  // spiral arms
  ham_float spiral(const ham_float &, const ham_float &, const ham_float &,
                   const ham_float &, const Param *) const;
  // galactic center
  ham_float galcen(const ham_float &, const ham_float &, const ham_float &,
                   const Param *) const;
  // gum nebula
  ham_float gum(const ham_float &, const ham_float &, const ham_float &,
                const Param *) const;
  // local bubble
  ham_float localbubble(const ham_float &, const ham_float &, const ham_float &,
                        const ham_float &, const ham_float &,
                        const Param *) const;
  // northern polar spurs
  ham_float nps(const ham_float &, const ham_float &, const ham_float &,
                const Param *) const;
};

// default method of global random FE
class TEmodel_dft final : public TEmodel {
public:
  TEmodel_dft() = default;
  TEmodel_dft(const TEmodel_dft &) = delete;
  TEmodel_dft(TEmodel_dft &&) = delete;
  TEmodel_dft &operator=(const TEmodel_dft &) = delete;
  TEmodel_dft &operator=(TEmodel_dft &&) = delete;
  virtual ~TEmodel_dft() = default;
  // trivial Fourier transform, with rescaling applied in spatial space
  void write_grid(const Param *, Grid_te *) const override;
#ifndef NDEBUG
protected:
#endif
  // isotropic turubulent power spectrum
  virtual ham_float spectrum(const ham_float &, const Param *) const;
  // density variance rescaling factor
  virtual ham_float spatial_profile(const Hamvec<3, ham_float> &,
                                    const Param *) const;
};

// thermal electron field collection class
class TEfield {
public:
  TEfield() = default;
  TEfield(const TEfield &) = delete;
  TEfield(TEfield &&) = delete;
  TEfield &operator=(const TEfield &) = delete;
  TEfield &operator=(TEfield &&) = delete;
  virtual ~TEfield() = default;
  // initialize thermal electron field model collection
  // model ordering arranged during parameter parsing
  TEfield(Param *par) {
    if (par->temodel_list.empty()) {
      model_list.push_back(std::make_unique<TEmodel>());
    } else {
      for (auto &name : par->temodel_list) {
        if (name == "unif") {
          model_list.push_back(std::make_unique<TEmodel_unif>());
        } else if (name == "ymw16") {
          model_list.push_back(std::make_unique<TEmodel_ymw16>());
        } else if (name == "dft") {
          model_list.push_back(std::make_unique<TEmodel_dft>());
        } else
          throw std::runtime_error("unknown TE model");
      }
    }
  }
  // read thermal electron field at given grid index
  // 1st argument: thermal electron field grid index
  // 2nd argument: thermal electron field grid
  ham_float read_field(const ham_uint &idx, const Grid_te *grid) const {
    return grid->te[idx];
  }
  // read thermal electron field at given position (with interpolation)
  // 1st argument: Cartesian frame Galactic centric position
  // 2nd argument: parameter set
  // 34d argument: thermal electron field grid
  ham_float read_field(const Hamvec<3, ham_float> &pos, const Param *par,
                       const Grid_te *grid) const {
    ham_float tmp{0};
    // get field content from grid
    // model ordering arranged during parameter parsing
    if (par->grid_te.build_permission or par->grid_te.read_permission) {
      tmp = model_list[0]->read_grid(pos, par, grid);
    }
    // get field content directly from model
    // model ordering arranged during parameter parsing
    else {
      for (auto &model : model_list) {
        tmp += model->read_model(pos, par);
      }
    }
    return tmp;
  }
  // write thermal electron field to grid before reading
  // 1st argument: parameter set
  // 2nd argument: thermal electron field grid
  void write_field(const Param *par, Grid_te *grid) {
    // write field content to grid
    // model ordering arranged during parameter parsing
    assert(par->grid_te.build_permission or par->grid_te.write_permission or
           par->grid_te.read_permission);
    for (auto &model : model_list) {
      model->write_grid(par, grid);
    }
  }

protected:
  // list of thermal electron field model pointers
  std::vector<std::unique_ptr<TEmodel>> model_list;
};

#endif
