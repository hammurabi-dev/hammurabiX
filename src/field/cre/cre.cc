#include <cmath>
#include <cassert>
#include <stdexcept>

#include <hvec.h>

#include <cre.h>
#include <param.h>
#include <grid.h>
#include <namespace_toolkit.h>

double CRE::flux (const hvec<3,double> &,
                  const Param *,
                  const double &) const{
    throw std::runtime_error("wrong inheritance");
}

double CRE::get_emissivity_t (const hvec<3,double> &,
                              const Param *,
                              const Grid_cre *,
                              const double &) const{
    throw std::runtime_error("wrong inheritance");
}

double CRE::get_emissivity_p (const hvec<3,double> &,
                              const Param *,
                              const Grid_cre *,
                              const double &) const{
    throw std::runtime_error("wrong inheritance");
}

// use bilinear/trilinear interpolationi according to the dimension of CRE flux grid
double CRE::read_grid (const std::size_t &Eidx,
                       const hvec<3,double> &pos,
                       const Param *par,
                       const Grid_cre *grid) const{
    //trilinear interpolation
    double tmp {(par->grid_cre.nx-1)*(pos[0]-par->grid_cre.x_min)/(par->grid_cre.x_max-par->grid_cre.x_min)};
    if (tmp<0 or tmp>par->grid_cre.nx-1) { return 0.;}
    decltype(par->grid_cre.nx) xl {(std::size_t)std::floor(tmp)};
    const double xd {tmp - xl};
    tmp = (par->grid_cre.ny-1)*(pos[1]-par->grid_cre.y_min)/(par->grid_cre.y_max-par->grid_cre.y_min);
    if (tmp<0 or tmp>par->grid_cre.ny-1) { return 0.;}
    decltype(par->grid_cre.nx) yl {(std::size_t)std::floor(tmp)};
    const double yd {tmp - yl};
    tmp = (par->grid_cre.nz-1)*(pos[2]-par->grid_cre.z_min)/(par->grid_cre.z_max-par->grid_cre.z_min);
    if (tmp<0 or tmp>par->grid_cre.nz-1) { return 0.;}
    decltype(par->grid_cre.nx) zl {(std::size_t)std::floor(tmp)};
    const double zd {tmp - zl};
    assert(xd>=0 and yd>=0 and zd>=0 and xd<=1 and yd<=1 and zd<=1);
    double cre;
    if (xl+1<par->grid_cre.nx and yl+1<par->grid_cre.ny and zl+1<par->grid_cre.nz) {
        std::size_t idx1 {toolkit::index4d(par->grid_cre.nE,par->grid_cre.nx,par->grid_cre.ny,par->grid_cre.nz,Eidx,xl,yl,zl)};
        std::size_t idx2 {toolkit::index4d(par->grid_cre.nE,par->grid_cre.nx,par->grid_cre.ny,par->grid_cre.nz,Eidx,xl,yl,zl+1)};
        const double i1 {grid->cre_flux[idx1]*(1.-zd) + grid->cre_flux[idx2]*zd};
        idx1 = toolkit::index4d(par->grid_cre.nE,par->grid_cre.nx,par->grid_cre.ny,par->grid_cre.nz,Eidx,xl,yl+1,zl);
        idx2 = toolkit::index4d(par->grid_cre.nE,par->grid_cre.nx,par->grid_cre.ny,par->grid_cre.nz,Eidx,xl,yl+1,zl+1);
        const double i2 {grid->cre_flux[idx1]*(1-zd) + grid->cre_flux[idx2]*zd};
        idx1 = toolkit::index4d(par->grid_cre.nE,par->grid_cre.nx,par->grid_cre.ny,par->grid_cre.nz,Eidx,xl+1,yl,zl);
        idx2 = toolkit::index4d(par->grid_cre.nE,par->grid_cre.nx,par->grid_cre.ny,par->grid_cre.nz,Eidx,xl+1,yl,zl+1);
        const double j1 {grid->cre_flux[idx1]*(1-zd) + grid->cre_flux[idx2]*zd};
        idx1 = toolkit::index4d(par->grid_cre.nE,par->grid_cre.nx,par->grid_cre.ny,par->grid_cre.nz,Eidx,xl+1,yl+1,zl);
        idx2 = toolkit::index4d(par->grid_cre.nE,par->grid_cre.nx,par->grid_cre.ny,par->grid_cre.nz,Eidx,xl+1,yl+1,zl+1);
        const double j2 {grid->cre_flux[idx1]*(1-zd) + grid->cre_flux[idx2]*zd};
        const double w1 {i1*(1-yd)+i2*yd};
        const double w2 {j1*(1-yd)+j2*yd};
        cre = (w1*(1-xd)+w2*xd);
        
    }
    else {
        std::size_t idx1 {toolkit::index4d(par->grid_cre.nE,par->grid_cre.nx,par->grid_cre.ny,par->grid_cre.nz,Eidx,xl,yl,zl)};
        cre = grid->cre_flux[idx1];
    }
    assert(cre>=0);
    return cre;
}

// writing out CRE DIFFERENTIAL density flux, [GeV m^2 s sr]^-1
void CRE::write_grid (const Param *par,
                      Grid_cre *grid) const{
    assert(par->grid_cre.write_permission);
    hvec<3,double> gc_pos;
    
    double lx {par->grid_cre.x_max-par->grid_cre.x_min};
    double ly {par->grid_cre.y_max-par->grid_cre.y_min};
    double lz {par->grid_cre.z_max-par->grid_cre.z_min};
    for(decltype(par->grid_cre.nE) i=0;i!=par->grid_cre.nE;++i){
        double E {par->grid_cre.E_min*std::exp(i*par->grid_cre.E_fact)};
        for(decltype(par->grid_cre.nx) j=0;j!=par->grid_cre.nx;++j){
            gc_pos[0] = lx*j/(par->grid_cre.nx-1) + par->grid_cre.x_min;
            for(decltype(par->grid_cre.ny) k=0;k!=par->grid_cre.ny;++k){
                gc_pos[1] = ly*k/(par->grid_cre.ny-1) + par->grid_cre.y_min;
                for(decltype(par->grid_cre.nz) m=0;m!=par->grid_cre.nz;++m){
                    gc_pos[2] = lz*m/(par->grid_cre.nz-1) + par->grid_cre.z_min;
                    std::size_t idx {toolkit::index4d(par->grid_cre.nE,par->grid_cre.nx,par->grid_cre.ny,par->grid_cre.nz,i,j,k,m)};
                    grid->cre_flux[idx] = flux (gc_pos,par,E);
                }
            }
        }
    }
}

// END
