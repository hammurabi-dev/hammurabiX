#include <iostream>
#include <vec3.h>
#include <vector>
#include <array>
#include <cmath>
#include <gsl/gsl_sf_synchrotron.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_integration.h>
#include <cre.h>
#include <param.h>
#include <grid.h>
#include <cgs_units_file.h>
#include <namespace_toolkit.h>
#include <cassert>
using namespace std;

double CRE::flux(const vec3_t<double> &,Param *,const double &){
    assert(false);
    return 0;
}

double CRE::get_emissivity_t(const vec3_t<double> &,Param *,Grid_cre *,const double &){
    assert(false);
    return 0;
}

double CRE::get_emissivity_p(const vec3_t<double> &,Param *,Grid_cre *,const double &){
    assert(false);
    return 0;
}

// use bilinear/trilinear interpolationi according to the dimension of CRE flux grid
double CRE::read_grid(const std::size_t &Eidx, const vec3_t<double> &pos,Grid_cre *grid){
    //trilinear interpolation
    double tmp {(grid->nx-1)*(pos.x-grid->x_min)/(grid->x_max-grid->x_min)};
    if (tmp<0 or tmp>grid->nx-1) { return 0.;}
    decltype(grid->nx) xl {(std::size_t)floor(tmp)};
    const double xd {tmp - xl};
    tmp = (grid->ny-1)*(pos.y-grid->y_min)/(grid->y_max-grid->y_min);
    if (tmp<0 or tmp>grid->ny-1) { return 0.;}
    decltype(grid->nx) yl {(std::size_t)floor(tmp)};
    const double yd {tmp - yl};
    tmp = (grid->nz-1)*(pos.z-grid->z_min)/(grid->z_max-grid->z_min);
    if (tmp<0 or tmp>grid->nz-1) { return 0.;}
    decltype(grid->nx) zl {(std::size_t)floor(tmp)};
    const double zd {tmp - zl};
    assert(xd>=0 and yd>=0 and zd>=0 and xd<=1 and yd<=1 and zd<=1);
    double cre;
    if (xl+1<grid->nx and yl+1<grid->ny and zl+1<grid->nz) {
        std::size_t idx1 {toolkit::Index4d(grid->nE,grid->nx,grid->ny,grid->nz,Eidx,xl,yl,zl)};
        std::size_t idx2 {toolkit::Index4d(grid->nE,grid->nx,grid->ny,grid->nz,Eidx,xl,yl,zl+1)};
        const double i1 {grid->cre_flux[idx1]*(1.-zd) + grid->cre_flux[idx2]*zd};
        idx1 = toolkit::Index4d(grid->nE,grid->nx,grid->ny,grid->nz,Eidx,xl,yl+1,zl);
        idx2 = toolkit::Index4d(grid->nE,grid->nx,grid->ny,grid->nz,Eidx,xl,yl+1,zl+1);
        const double i2 {grid->cre_flux[idx1]*(1-zd) + grid->cre_flux[idx2]*zd};
        idx1 = toolkit::Index4d(grid->nE,grid->nx,grid->ny,grid->nz,Eidx,xl+1,yl,zl);
        idx2 = toolkit::Index4d(grid->nE,grid->nx,grid->ny,grid->nz,Eidx,xl+1,yl,zl+1);
        const double j1 {grid->cre_flux[idx1]*(1-zd) + grid->cre_flux[idx2]*zd};
        idx1 = toolkit::Index4d(grid->nE,grid->nx,grid->ny,grid->nz,Eidx,xl+1,yl+1,zl);
        idx2 = toolkit::Index4d(grid->nE,grid->nx,grid->ny,grid->nz,Eidx,xl+1,yl+1,zl+1);
        const double j2 {grid->cre_flux[idx1]*(1-zd) + grid->cre_flux[idx2]*zd};
        const double w1 {i1*(1-yd)+i2*yd};
        const double w2 {j1*(1-yd)+j2*yd};
        cre = (w1*(1-xd)+w2*xd);
        
    }
    else {
        std::size_t idx1 {toolkit::Index4d(grid->nE,grid->nx,grid->ny,grid->nz,Eidx,xl,yl,zl)};
        cre = grid->cre_flux[idx1];
    }
    assert(cre>=0);
    return cre;
}

// writing out CRE DIFFERENTIAL density flux, [GeV m^2 s sr]^-1
void CRE::write_grid(Param *par, Grid_cre *grid){
    assert(grid->write_permission);
    vec3_t<double> gc_pos {0.,0.,0.};
    
    double lx {grid->x_max-grid->x_min};
    double ly {grid->y_max-grid->y_min};
    double lz {grid->z_max-grid->z_min};
    for(decltype(grid->nE) i=0;i!=grid->nE;++i){
        double E {grid->E_min*exp(i*grid->E_fact)};
        for(decltype(grid->nx) j=0;j!=grid->nx;++j){
            gc_pos.x = lx*j/(grid->nx-1) + grid->x_min;
            for(decltype(grid->ny) k=0;k!=grid->ny;++k){
                gc_pos.y = ly*k/(grid->ny-1) + grid->y_min;
                for(decltype(grid->nz) m=0;m!=grid->nz;++m){
                    gc_pos.z = lz*m/(grid->nz-1) + grid->z_min;
                    std::size_t idx {toolkit::Index4d(grid->nE,grid->nx,grid->ny,grid->nz,i,j,k,m)};
                    grid->cre_flux[idx] = flux(gc_pos,par,E);
                }
            }
        }
    }
}

// END
