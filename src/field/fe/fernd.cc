#include <cmath>
#include <cassert>
#include <stdexcept>

#include <vec3.h>

#include <param.h>
#include <grid.h>
#include <fernd.h>
#include <fereg.h>
#include <namespace_toolkit.h>

double FErnd::get_fernd (const hvec<3,double> &pos,
                         const Param *par,
                         const Grid_fernd *grid) const{
    if (par->grid_fernd.read_permission or par->grid_fernd.build_permission){
        return read_grid (pos,
                          par,
                          grid);
    }
    // if no specific turbulent field model is called
    // base class will return zero vector
    else {
        return 0.;
    }
}

double FErnd::read_grid (const hvec<3,double> &pos,
                         const Param *par,
                         const Grid_fernd *grid) const{
    double tmp {(par->grid_fernd.nx-1)*(pos[0]-par->grid_fernd.x_min)/(par->grid_fernd.x_max-par->grid_fernd.x_min)};
    if (tmp<0 or tmp>par->grid_fernd.nx-1) { return 0.;}
    decltype(par->grid_fernd.nx) xl {(std::size_t)std::floor(tmp)};
    const double xd {tmp - xl};
    
    tmp = (par->grid_fernd.ny-1)*(pos[1]-par->grid_fernd.y_min)/(par->grid_fernd.y_max-par->grid_fernd.y_min);
    if (tmp<0 or tmp>par->grid_fernd.ny-1) { return 0.;}
    decltype(par->grid_fernd.nx) yl {(std::size_t)std::floor(tmp)};
    const double yd {tmp - yl};
    
    tmp = (par->grid_fernd.nz-1)*(pos[2]-par->grid_fernd.z_min)/(par->grid_fernd.z_max-par->grid_fernd.z_min);
    if (tmp<0 or tmp>par->grid_fernd.nz-1) { return 0.;}
    decltype(par->grid_fernd.nx) zl {(std::size_t)std::floor(tmp)};
    const double zd {tmp - zl};
    assert(xd>=0 and yd>=0 and zd>=0 and xd<=1 and yd<=1 and zd<=1);
    double density;
    //trilinear interpolation
    if (xl+1<par->grid_fernd.nx and yl+1<par->grid_fernd.ny and zl+1<par->grid_fernd.nz) {
        //interpolate along z direction, there are four interpolated vectors
        std::size_t idx1 {toolkit::index3d(par->grid_fernd.nx,par->grid_fernd.ny,par->grid_fernd.nz,xl,yl,zl)};
        std::size_t idx2 {toolkit::index3d(par->grid_fernd.nx,par->grid_fernd.ny,par->grid_fernd.nz,xl,yl,zl+1)};
        double i1 {grid->fe[idx1]*(1.-zd) + grid->fe[idx2]*zd};
        
        idx1 = toolkit::index3d(par->grid_fernd.nx,par->grid_fernd.ny,par->grid_fernd.nz,xl,yl+1,zl);
        idx2 = toolkit::index3d(par->grid_fernd.nx,par->grid_fernd.ny,par->grid_fernd.nz,xl,yl+1,zl+1);
        double i2 {grid->fe[idx1]*(1.-zd) + grid->fe[idx2]*zd};
        
        idx1 = toolkit::index3d(par->grid_fernd.nx,par->grid_fernd.ny,par->grid_fernd.nz,xl+1,yl,zl);
        idx2 = toolkit::index3d(par->grid_fernd.nx,par->grid_fernd.ny,par->grid_fernd.nz,xl+1,yl,zl+1);
        double j1 {grid->fe[idx1]*(1.-zd) + grid->fe[idx2]*zd};
        
        idx1 = toolkit::index3d(par->grid_fernd.nx,par->grid_fernd.ny,par->grid_fernd.nz,xl+1,yl+1,zl);
        idx2 = toolkit::index3d(par->grid_fernd.nx,par->grid_fernd.ny,par->grid_fernd.nz,xl+1,yl+1,zl+1);
        double j2 {grid->fe[idx1]*(1.-zd) + grid->fe[idx2]*zd};
        // interpolate along y direction, two interpolated vectors
        double w1 {i1*(1.-yd) + i2*yd};
        double w2 {j1*(1.-yd) + j2*yd};
        // interpolate along x direction
        density = w1*(1.-xd) + w2*xd;
    }
    // on the boundary
    else {
        std::size_t idx{toolkit::index3d(par->grid_fernd.nx,par->grid_fernd.ny,par->grid_fernd.nz,xl,yl,zl)};
        density = grid->fe[idx];
    }
    return density;
}

void FErnd::write_grid (const Param *,
                        Grid_fernd *) const{
    throw std::runtime_error("wrong inheritance");
}

// END
