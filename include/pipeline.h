#ifndef HAMMURABI_PIPELINE_H
#define HAMMURABI_PIPELINE_H

#include <bfield.h>
#include <crefield.h>
#include <grid.h>
#include <integrator.h>
#include <param.h>
#include <string>
#include <tefield.h>
#include <toolkit.h>

class Pipeline {
public:
  Pipeline() = default;
  Pipeline(const std::string &);
  Pipeline(const Pipeline &) = delete;
  Pipeline(Pipeline &&) = delete;
  Pipeline &operator=(const Pipeline &) = delete;
  Pipeline &operator=(Pipeline &&) = delete;
  virtual ~Pipeline() = default;
  virtual void assemble_grid();
  virtual void assemble_b();
  virtual void assemble_te();
  virtual void assemble_cre();
  virtual void assemble_obs();

protected:
  std::unique_ptr<Param> par;
  std::unique_ptr<Grid_te> grid_te;
  std::unique_ptr<Grid_b> grid_b;
  std::unique_ptr<Grid_cre> grid_cre;
  std::unique_ptr<Grid_obs> grid_obs;
  std::unique_ptr<Bfield> bfield;
  std::unique_ptr<TEfield> tefield;
  std::unique_ptr<CREfield> crefield;
  std::unique_ptr<Integrator> integrator;
};

#endif
