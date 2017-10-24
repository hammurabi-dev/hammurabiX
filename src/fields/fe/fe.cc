#include <iostream>
#include <cmath>
#include <vec3.h>
#include <omp.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <string>
#include <fstream>
#include "fe.h"
#include "pond.h"
#include "grid.h"
#include "namespace_toolkit.h"
#include "cgs_units_file.h"

using namespace std;

double FE::get_density(const vec3 &pos, Pond *par, Grid_fe *grid){
    if(grid->read_permission){
        return read_grid(pos,grid);
    }
    else {
        return density(pos,par);
    }
}

// not recommended to use without enough computing source
// recommend to use this once (replace density in write_grid) if no
// free parameters in FE
double FE::density_blur(const vec3 &pos, Pond *par, Grid_fe *grid){
    double ne_blur {0.};
    // element value in each thread
    double element {0.};
    // sampling point number
    unsigned int step {1000};
    // gaussian blur scale
    double blur_scale_x {(grid->x_max-grid->x_min)/(grid->nx*CGS_U_kpc)};
    double blur_scale_y {(grid->y_max-grid->y_min)/(grid->ny*CGS_U_kpc)};
    double blur_scale_z {(grid->z_max-grid->z_min)/(grid->nz*CGS_U_kpc)};
    // sample position
    vec3 pos_s;
    gsl_rng *r {gsl_rng_alloc(gsl_rng_taus)};
    gsl_rng_set(r, toolkit::random_seed());
#pragma omp parallel sections private(element) reduction(+:ne_blur)
    {
#pragma omp section
        {
            for(decltype(step)i=0;i!=step;++i){
                pos_s = pos + vec3 {gsl_ran_gaussian(r,(blur_scale_x/2.355))*CGS_U_kpc,
                    gsl_ran_gaussian(r,(blur_scale_y/2.355))*CGS_U_kpc,
                    gsl_ran_gaussian(r,(blur_scale_z/2.355))*CGS_U_kpc};
                element = density(pos_s,par);
                ne_blur += element;
            }
        }// omp section
    }//parallel
    gsl_rng_free(r);
    return ne_blur/step;
}

double FE::density(const vec3 &, Pond *){
    cerr<<"ERR:"<<__FILE__
    <<" : in function "<<__func__<<endl
    <<" at line "<<__LINE__<<endl
    <<"DYNAMIC BINDING FAILURE"<<endl;
    exit(1);
    return 0.;
}

double FE::read_grid(const vec3 &pos, Grid_fe *grid){
    double tmp {(grid->nx-1)*(pos.x-grid->x_min)/(grid->x_max-grid->x_min)};
    if (tmp<1 or tmp>grid->nx-1) { return 0.;}
    decltype(grid->nx) xl {(unsigned int)floor(tmp)};
    const double xd = tmp - xl;
    
    tmp = (grid->ny-1)*(pos.y-grid->y_min)/(grid->y_max-grid->y_min);
    if (tmp<1 or tmp>grid->ny-1) { return 0.;}
    decltype(grid->nx) yl {(unsigned int)floor(tmp)};
    const double yd = tmp - yl;
    
    tmp = (grid->nz-1)*(pos.z-grid->z_min)/(grid->z_max-grid->z_min);
    if (tmp<1 or tmp>grid->nz-1) { return 0.;}
    decltype(grid->nx) zl {(unsigned int)floor(tmp)};
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
    double fe;
    if (xl+1<grid->nx and yl+1<grid->ny and zl+1<grid->nz) {
        unsigned long int idx1 {toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl,yl,zl)};
        unsigned long int idx2 {toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl,yl,zl+1)};
        double i1 {grid->fe[idx1]*(1.-zd) + grid->fe[idx2]*zd};
        idx1 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl,yl+1,zl);
        idx2 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl,yl+1,zl+1);
        double i2 {grid->fe[idx1]*(1-zd) + grid->fe[idx2]*zd};
        idx1 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl+1,yl,zl);
        idx2 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl+1,yl,zl+1);
        double j1 {grid->fe[idx1]*(1-zd) + grid->fe[idx2]*zd};
        idx1 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl+1,yl+1,zl);
        idx2 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl+1,yl+1,zl+1);
        double j2 {grid->fe[idx1]*(1-zd) + grid->fe[idx2]*zd};
        double w1 {i1*(1-yd)+i2*yd};
        double w2 {j1*(1-yd)+j2*yd};
        fe = (w1*(1-xd)+w2*xd);
    }
    else {
        unsigned long int idx1 {toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl,yl,zl)};
        fe = grid->fe[idx1];
    }
#ifndef NDEBUG
    if(fe<0){
        cerr<<"WAR:"<<__FILE__
        <<" : in function "<<__func__<<endl
        <<" at line "<<__LINE__<<endl
        <<"WRONG VALUE"<<endl;
        exit(1);
    }
#endif
    return fe;
}

void FE::write_grid(Pond *par, Grid_fe *grid){
    if(!grid->write_permission){
        cerr<<"ERR:"<<__FILE__
        <<" : in function "<<__func__<<endl
        <<" at line "<<__LINE__<<endl
        <<"NO PERMISSION"<<endl;
        exit(1);
    }
    cout<<"...FE: WRITING OUTPUT..."<<endl;
    vec3 gc_pos;
    double lx {grid->x_max-grid->x_min};
    double ly {grid->y_max-grid->y_min};
    double lz {grid->z_max-grid->z_min};
    for(decltype(grid->nx) i=0;i!=grid->nx;++i){
        gc_pos.x = lx*i/(grid->nx-1) + grid->x_min;
        for(decltype(grid->ny) j=0;j!=grid->ny;++j){
            gc_pos.y = ly*j/(grid->ny-1) + grid->y_min;
            for(decltype(grid->nz) k=0;k!=grid->nz;++k){
                unsigned long int idx {toolkit::Index3d(grid->nx,grid->ny,grid->nz,i,j,k)};
                gc_pos.z = lz*k/(grid->nz-1) + grid->z_min;
                // two solutions
                //grid->fe[idx] = density_blur(gc_pos, par, grid);
                grid->fe[idx] = density(gc_pos,par);
            }
        }
    }
}

// END
