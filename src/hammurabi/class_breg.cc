#include <iostream>
#include <vec3.h>
#include <vector>
#include <array>
#include <cmath>

#include "class_pond.h"
#include "class_grid.h"
#include "class_breg.h"
#include "cgs_units_file.h"
#include "namespace_toolkit.h"

using namespace std;

vec3 Breg::get_breg(const vec3 &pos, Pond *par,Grid_breg *grid){
    if(grid->permission){
        return read_grid(pos,grid,par);
    }
    else {
        return breg(pos,par);
    }
}

vec3 Breg::breg(const vec3 &pos, Pond *par){
    cerr<<"ERR:"<<__FILE__
    <<" : in function "<<__func__<<endl
    <<" at line "<<__LINE__<<endl
    <<"DYNAMIC BINDING FAILURE"<<endl;
    exit(1);
    return vec3(0.,0.,0.);
}

vec3 Breg::read_grid(const vec3 &pos, Grid_breg *grid, Pond *par){
    
    // get position in sim box
    decltype(grid->nx) xl, yl, zl;
    
    double tmp = (grid->nx-1)*(pos.x/grid->lx + 0.5);
    // if use solar-centric frame
    if(grid->ec_frame) tmp -= (grid->nx-1)*(par->SunPosition.x/grid->lx);
    if (tmp<0 or tmp>grid->nx-1) { return vec3(0.,0.,0.);}
    else xl = floor(tmp);
    const double xd = tmp - xl;
    
    tmp = (grid->ny-1)*(pos.y/grid->ly + 0.5);
    if (tmp<0 or tmp>grid->ny-1) { return vec3(0.,0.,0.);}
    else yl = floor(tmp);
    const double yd = tmp - yl;
    
    tmp = (grid->nz-1)*(pos.z/grid->lz + 0.5);
    if (tmp<0 or tmp>grid->nz-1) { return vec3(0.,0.,0.);}
    else zl = floor(tmp);
    const double zd = tmp - zl;
#ifndef NDEBUG
    if(xd<0 or yd<0 or zd<0 or xd>1 or yd>1 or zd>1){
        cerr<<"ERR:"<<__FILE__
        <<" : in function "<<__func__<<endl
        <<" at line "<<__LINE__<<endl
        <<"WRONG VALUE: "<<endl;
        exit(1);
    }
#endif
    vec3 b_vec3;
    // trilinear interpolation
    if (xl+1<grid->nx and yl+1<grid->ny and zl+1<grid->nz) {
        // for storing vectors
        vec3 i1,i2,j1,j2,w1,w2;
        
        // interpolate along z direction, there are four interpolated vectors
        auto idx1 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl,yl,zl);
        auto idx2 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl,yl,zl+1);
        i1 = vec3(grid->reg_b_x[idx1]*(1.-zd) + grid->reg_b_x[idx2]*zd,
                  grid->reg_b_y[idx1]*(1.-zd) + grid->reg_b_y[idx2]*zd,
                  grid->reg_b_z[idx1]*(1.-zd) + grid->reg_b_z[idx2]*zd );
        
        idx1 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl,yl+1,zl);
        idx2 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl,yl+1,zl+1);
        i2 = vec3(grid->reg_b_x[idx1]*(1.-zd) + grid->reg_b_x[idx2]*zd,
                  grid->reg_b_y[idx1]*(1.-zd) + grid->reg_b_y[idx2]*zd,
                  grid->reg_b_z[idx1]*(1.-zd) + grid->reg_b_z[idx2]*zd);
        
        idx1 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl+1,yl,zl);
        idx2 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl+1,yl,zl+1);
        j1 = vec3(grid->reg_b_x[idx1]*(1.-zd) + grid->reg_b_x[idx2]*zd,
                  grid->reg_b_y[idx1]*(1.-zd) + grid->reg_b_y[idx2]*zd,
                  grid->reg_b_z[idx1]*(1.-zd) + grid->reg_b_z[idx2]*zd);
        
        idx1 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl+1,yl+1,zl);
        idx2 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl+1,yl+1,zl+1);
        j2 = vec3(grid->reg_b_x[idx1]*(1.-zd) + grid->reg_b_x[idx2]*zd,
                  grid->reg_b_y[idx1]*(1.-zd) + grid->reg_b_y[idx2]*zd,
                  grid->reg_b_z[idx1]*(1.-zd) + grid->reg_b_z[idx2]*zd);
        // interpolate along y direction, two interpolated vectors
        w1 = i1*(1.-yd) + i2*yd;
        w2 = j1*(1.-yd) + j2*yd;
        // interpolate along x direction
        b_vec3 = w1*(1.-xd) + w2*xd;
    }
    // no interpolation
    else {
        auto idx = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl,yl,zl);
        b_vec3 = vec3(grid->reg_b_x[idx],grid->reg_b_y[idx],grid->reg_b_z[idx]);
    }
#ifndef NDEBUG
    if (b_vec3.Length()>50.*CGS_U_muGauss) {
        cerr<<"WAR:"<<__FILE__
        <<" : in function "<<__func__<<endl
        <<" at line "<<__LINE__<<endl
        <<" too strong field at xl="<<xl<<", yl="<<yl<<", zl="<<zl<<endl
        <<" field amplitude: "<< b_vec3.Length()/CGS_U_muGauss <<" microGauss"<<endl;
        exit(1);
    }
#endif
    return b_vec3;
}

void Breg::write_grid(Pond *par, Grid_breg *grid){
    if(!grid->permission){
        cerr<<"ERR:"<<__FILE__
        <<" : in function "<<__func__<<endl
        <<" at line "<<__LINE__<<endl
        <<"NO PERMISSION"<<endl;
        exit(1);
    }
    
    vec3 gc_pos, tmp_vec;
    for(decltype(grid->nx) i=0;i!=grid->nx;++i){
        for(decltype(grid->ny) j=0;j!=grid->ny;++j){
            for(decltype(grid->nz) k=0;k!=grid->nz;++k){
                
                gc_pos.x = grid->lx*(double(i)/(grid->nx-1) - 0.5);
                gc_pos.y = grid->ly*(double(j)/(grid->ny-1) - 0.5);
                gc_pos.z = grid->lz*(double(k)/(grid->nz-1) - 0.5);
                
                if(grid->ec_frame){
                    gc_pos.x += par->SunPosition.x;
                }
                
                tmp_vec = breg(gc_pos,par);
                auto idx = toolkit::Index3d(grid->nx,grid->ny,grid->nz,i,j,k);
                
                grid->reg_b_x[idx] = tmp_vec.x;
                grid->reg_b_y[idx] = tmp_vec.y;
                grid->reg_b_z[idx] = tmp_vec.z;
            }
        }
    }
}

/* wmap-3yr */
/*
Bwmap::Bwmap(const vector<double> &in, Pond *par){
    if(in.size()!=par->bwmap.size()){
        cerr<<"WAR:"<<__FILE__
        <<" : in function "<<__func__<<endl
        <<" at line "<<__LINE__<<endl
        <<"Bwmap: NOT ASSIGNING FULL PAR-SET"<<endl;
    }
    for(decltype(in.size()) i=0;i!=in.size();++i){
        par->bwmap[i] = in[i];
    }
}
*/
vec3 Bwmap::breg(const vec3 &pos,Pond *par){
    vec3 b_vec3(0.,0.,0.);
    const double r = sqrt(pos.x*pos.x + pos.y*pos.y);
    if (r>(20.*CGS_U_kpc) or r<(3.*CGS_U_kpc )) {
        return b_vec3;
    }
    // units
    const double b0 = par->bwmap[0]*CGS_U_muGauss;
    const double psi0 = par->bwmap[1]*CGS_U_pi/180.;
    const double psi1 = par->bwmap[2]*CGS_U_pi/180.;
    const double chi0 = par->bwmap[3]*CGS_U_pi/180.;
    
    const double phi = atan2(pos.y,pos.x);
    const double psi = psi0 + psi1*log(r/(8.*CGS_U_kpc));
    const double chi = chi0*tanh(pos.z/(1.*CGS_U_kpc));
    const vec3 b_cyl(
                     b0*sin(psi)*cos(chi),
                     b0*cos(psi)*cos(chi),
                     b0*sin(chi)
                     );
    toolkit::Cyl2Cart(phi,b_cyl,b_vec3);
    return b_vec3;
}

/* local */
/*
Blocal::Blocal(const vector<double> &in, Pond *par){
    if(in.size()!=par->blocal.size()){
        cerr<<"WAR:"<<__FILE__
        <<" : in function "<<__func__<<endl
        <<" at line "<<__LINE__<<endl
        <<"Blocal: NOT ASSIGNING FULL PAR-SET"<<endl;
    }
    for(decltype(in.size()) i=0;i!=in.size();++i){
        par->blocal[i] = in[i];
    }
}
*/
vec3 Blocal::breg(const vec3 &pos, Pond *par){
    // units
    const double bd = par->blocal[0]*CGS_U_muGauss;
    const double l0 = par->blocal[1]*CGS_U_pi/180.;
    const double z0 = par->blocal[2]*CGS_U_kpc;
    
    const double bn = par->blocal[3]*CGS_U_muGauss;
    const double bs = par->blocal[4]*CGS_U_muGauss;
    
    // sylindrical coordinates
    double r_cyl, phi, z;
    toolkit::cart_coord2cyl_coord(pos,r_cyl,phi,z);
    
    // disk-halo scaling
    const double L = exp(-fabs(z)/z0);
    // total field
    vec3 b_tot(0.,0.,0.);
    
    // disk
    // we ignore global r scaling of disk field
    b_tot.x = bd*L*cos(l0);
    b_tot.y = bd*L*sin(l0);
    
    // poloidal
    // constant field, 2 param
    
    if(z>0){
        b_tot.z += bn;
    }
    else{
        b_tot.z += bs;
    }
    
    //toroidal
    /*
    if(z>0){
        b_tot.x += -bn*sin(phi)*(1.-L);
        b_tot.y += bn*cos(phi)*(1.-L);
    }
    else{
        b_tot.x += -bs*sin(phi)*(1.-L);
        b_tot.y += bs*cos(phi)*(1.-L);
    }
    */
    return b_tot;
}

// END
