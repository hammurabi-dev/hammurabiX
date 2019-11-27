#ifndef HAMMURABI_PIPELINE
#define HAMMURABI_PIPELINE

#include <string>

#include <bfield.h>
#include <crefield.h>
#include <grid.h>
#include <integrator.h>
#include <param.h>
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
  virtual void assemble_tereg();
  virtual void assemble_breg();
  virtual void assemble_ternd();
  virtual void assemble_brnd();
  virtual void assemble_cre();
  virtual void assemble_obs();

protected:
  std::unique_ptr<Param> par;
  std::unique_ptr<Grid_tereg> grid_tereg;
  std::unique_ptr<Grid_breg> grid_breg;
  std::unique_ptr<Grid_brnd> grid_brnd;
  std::unique_ptr<Grid_ternd> grid_ternd;
  std::unique_ptr<Grid_cre> grid_cre;
  std::unique_ptr<Grid_obs> grid_obs;
  std::unique_ptr<TEreg> tereg;
  std::unique_ptr<Breg> breg;
  std::unique_ptr<TErnd> ternd;
  std::unique_ptr<Brnd> brnd;
  std::unique_ptr<CREfield> cre;
  std::unique_ptr<Integrator> intobj;
};

#endif
