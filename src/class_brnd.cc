#include <iostream>
#include <vec3.h>
#include <vector>
#include <array>
#include <cmath>
#include <fftw3.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_integration.h>
#include "class_pond.h"
#include "class_grid.h"
#include "class_brnd.h"
#include "class_breg.h"
#include "cgs_units_file.h"
#include "namespace_toolkit.h"

using namespace std;

vec3 Brnd::get_brnd(const vec3 &pos, Grid_brnd *grid){
    if(grid->read_permission){
        return read_grid(pos,grid);
    }
    // if no specific random field model is called
    // base class will return zero vector
    else{
        return vec3 {0.,0.,0.};
    }
}

// since we are using rms normalization
// p0 is hidden and not affecting anything
double Brnd::b_spec(const double &k, Pond *par){
    //units fixing, wave vector in 1/kpc units
    const double p0 {par->brnd_iso[0]*pow(CGS_U_muGauss,2)};
    const double k0 {par->brnd_iso[1]};
    const double a0 {par->brnd_iso[2]};
    const double unit = 1./(4*CGS_U_pi*k*k);
    // avoid nan
    if(k<=0.){
        return 0.;
    }
    // Giacalone-Jokipii
    //const double P = 2*p0/(1.+pow(k/k0,a0));
    // power law
    double P{0.};
    if(k>k0){
        P = p0/pow(k/k0,a0);
    }
    return P*unit;
}

// galactic scaling of random field energy density
// set to 1 at observer's place
double Brnd::rescal_fact(const vec3 &pos, Pond *par){
    const double r_cyl {(sqrt(pos.x*pos.x+pos.y*pos.y) - fabs(par->SunPosition.x))/CGS_U_kpc};
    const double z {(fabs(pos.z) - fabs(par->SunPosition.z))/CGS_U_kpc};
    const double r0 {par->brnd_scal[0]};
    const double z0 {par->brnd_scal[1]};
    if(r_cyl==0. or z==0){return 1.;}
    else{
    	return exp(-r_cyl/r0)*exp(-z/z0);
    }
}

vec3 Brnd::read_grid(const vec3 &pos, Grid_brnd *grid){
    double tmp {(grid->nx-1)*(pos.x-grid->x_min)/(grid->x_max-grid->x_min)};
    if (tmp<0 or tmp>grid->nx-1) { return vec3 {0.,0.,0.};}
    decltype(grid->nx) xl {(unsigned int)floor(tmp)};
    const double xd {tmp - xl};
    
    tmp = (grid->ny-1)*(pos.y-grid->y_min)/(grid->y_max-grid->y_min);
    if (tmp<0 or tmp>grid->ny-1) { return vec3 {0.,0.,0.};}
    decltype(grid->nx) yl {(unsigned int)floor(tmp)};
    const double yd {tmp - yl};
    
    tmp = (grid->nz-1)*(pos.z-grid->z_min)/(grid->z_max-grid->z_min);
    if (tmp<0 or tmp>grid->nz-1) { return vec3 {0.,0.,0.};}
    decltype(grid->nx) zl {(unsigned int)floor(tmp)};
    const double zd {tmp - zl};
#ifndef NDEBUG
    if(xd<0 or yd<0 or zd<0 or xd>1 or yd>1 or zd>1){
        cerr<<"ERR:"<<__FILE__
        <<" : in function "<<__func__<<endl
        <<" at line "<<__LINE__<<endl
        <<"WRONG VALUE"<<endl;
        exit(1);
    }
#endif
    vec3 b_vec3;
    //trilinear interpolation
    if (xl+1<grid->nx and yl+1<grid->ny and zl+1<grid->nz) {
        unsigned long int idx1 {toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl,yl,zl)};
        unsigned long int idx2 {toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl,yl,zl+1)};
        vec3 i1 {grid->fftw_b_x[idx1]*(1.-zd) + grid->fftw_b_x[idx2]*zd,
            grid->fftw_b_y[idx1]*(1.-zd) + grid->fftw_b_y[idx2]*zd,
            grid->fftw_b_z[idx1]*(1.-zd) + grid->fftw_b_z[idx2]*zd};
        
        idx1 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl,yl+1,zl);
        idx2 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl,yl+1,zl+1);
        vec3 i2 {grid->fftw_b_x[idx1]*(1.-zd) + grid->fftw_b_x[idx2]*zd,
            grid->fftw_b_y[idx1]*(1.-zd) + grid->fftw_b_y[idx2]*zd,
            grid->fftw_b_z[idx1]*(1.-zd) + grid->fftw_b_z[idx2]*zd};
        
        idx1 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl+1,yl,zl);
        idx2 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl+1,yl,zl+1);
        vec3 j1 {grid->fftw_b_x[idx1]*(1.-zd) + grid->fftw_b_x[idx2]*zd,
            grid->fftw_b_y[idx1]*(1.-zd) + grid->fftw_b_y[idx2]*zd,
            grid->fftw_b_z[idx1]*(1.-zd) + grid->fftw_b_z[idx2]*zd};
        
        idx1 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl+1,yl+1,zl);
        idx2 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl+1,yl+1,zl+1);
        vec3 j2 {grid->fftw_b_x[idx1]*(1.-zd) + grid->fftw_b_x[idx2]*zd,
            grid->fftw_b_y[idx1]*(1.-zd) + grid->fftw_b_y[idx2]*zd,
            grid->fftw_b_z[idx1]*(1.-zd) + grid->fftw_b_z[idx2]*zd};
        // interpolate along y direction, two interpolated vectors
        vec3 w1 {i1*(1.-yd) + i2*yd};
        vec3 w2 {j1*(1.-yd) + j2*yd};
        // interpolate along x direction
        b_vec3 = w1*(1.-xd) + w2*xd;
    }
    // on the boundary
    else {
        unsigned long int idx {toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl,yl,zl)};
        b_vec3 = vec3 {grid->fftw_b_x[idx],grid->fftw_b_y[idx],grid->fftw_b_z[idx]};
    }
#ifndef NDEBUG
    if (b_vec3.Length()>50.*CGS_U_muGauss) {
        cerr<<"WAR:"<<__FILE__
        <<" : in function "<<__func__<<endl
        <<" at line "<<__LINE__<<endl
        <<" too strong field at xl="<<xl<<", yl="<<yl<<", zl="<<zl<<endl
        <<" b_strength = "<< b_vec3.Length()/CGS_U_muGauss <<" microGauss"<<endl;
        exit(1);
    }
#endif
    return b_vec3;
}

// base class does not write out to grid
void Brnd::write_grid_iso(Pond *,Grid_brnd *){
    cerr<<"WAR:"<<__FILE__
    <<" : in function "<<__func__<<endl
    <<" at line "<<__LINE__<<endl
    <<"DYNAMIC BINDING FAILURE"<<endl;
    exit(1);
}

void Brnd::write_grid_ani(Pond *, Breg *, Grid_breg *, Grid_brnd *){
    cerr<<"WAR:"<<__FILE__
    <<" : in function "<<__func__<<endl
    <<" at line "<<__LINE__<<endl
    <<"DYNAMIC BINDING FAILURE"<<endl;
    exit(1);
}

// isotropic turbulent field
vec3 Brnd_iso::get_brnd(const vec3 &pos, Grid_brnd *grid){
    // interpolate written grid to given position
    // check if you have called ::write_grid
    return read_grid(pos,grid);
}

void Brnd_iso::write_grid_iso(Pond *par, Grid_brnd *grid){
    // PHASE I
    // GENERATE GAUSSIAN RANDOM FROM SPECTRUM
    // initialize random seed
    gsl_rng *r {gsl_rng_alloc(gsl_rng_taus)};
    gsl_rng_set(r, toolkit::random_seed());
    // start Fourier space filling
    double lx {grid->x_max-grid->x_min};
    double ly {grid->y_max-grid->y_min};
    double lz {grid->z_max-grid->z_min};
    for (decltype(grid->nx) i=0;i!=grid->nx;++i) {
        for (decltype(grid->ny) j=0;j!=grid->ny;++j) {
            for (decltype(grid->nz) l=0;l!=grid->nz;++l) {
                unsigned long int idx {toolkit::Index3d(grid->nx,grid->ny,grid->nz,i,j,l)};
                // FFT expects up to n/2 positive while n/2 to n negative
                // physical k in 1/kpc dimension
                double kx {i/(lx/CGS_U_kpc)};
                double ky {j/(ly/CGS_U_kpc)};
                double kz {l/(lz/CGS_U_kpc)};
                if(i>=grid->nx/2) kx -= grid->nx/(lx/CGS_U_kpc);
                if(j>=grid->ny/2) ky -= grid->ny/(ly/CGS_U_kpc);
                if(l>=grid->nz/2) kz -= grid->nz/(lz/CGS_U_kpc);
                const double k {sqrt(kx*kx + ky*ky + kz*kz)};
                // physical dk^3
                const double dk3 {pow(CGS_U_kpc,3.)/(lx*ly*lz)};
                // simpson's rule
                double element {2.*b_spec(k,par)/3.};
                const double halfdk {0.5*sqrt(pow(CGS_U_kpc/lx,2)+pow(CGS_U_kpc/ly,2)+pow(CGS_U_kpc/lz,2))};
                element += b_spec(k+halfdk,par)/6.;
                element += b_spec(k-halfdk,par)/6.;
                // amplitude, dividing by two because equal allocation to Re and Im parts
                const double sigma {sqrt(element*dk3/2.0)};
                grid->fftw_b_kx[idx][0] = gsl_ran_gaussian(r,sigma);
                grid->fftw_b_ky[idx][0] = gsl_ran_gaussian(r,sigma);
                grid->fftw_b_kz[idx][0] = gsl_ran_gaussian(r,sigma);
                grid->fftw_b_kx[idx][1] = gsl_ran_gaussian(r,sigma);
                grid->fftw_b_ky[idx][1] = gsl_ran_gaussian(r,sigma);
                grid->fftw_b_kz[idx][1] = gsl_ran_gaussian(r,sigma);
                if(i==0 and j==0 and l==0) {
                    grid->fftw_b_kx[idx][0] = 0.;
                    grid->fftw_b_ky[idx][0] = 0.;
                    grid->fftw_b_kz[idx][0] = 0.;
                    grid->fftw_b_kx[idx][1] = 0.;
                    grid->fftw_b_ky[idx][1] = 0.;
                    grid->fftw_b_kz[idx][1] = 0.;
                }
            }// l
        }// j
    }// i
    // free random memory
    gsl_rng_free(r);
    // no Hermiticity fixing, complex 2 complex
    // execute DFT backward plan
    fftw_execute(grid->fftw_px_bw);
    fftw_execute(grid->fftw_py_bw);
    fftw_execute(grid->fftw_pz_bw);
    // PHASE II
    // RESCALING FIELD PROFILE IN REAL SPACE
    double b_var {toolkit::Variance(grid->fftw_b_kx[0], grid->full_size)};
    for (decltype(grid->nx) i=0;i!=grid->nx;++i) {
        for (decltype(grid->ny) j=0;j!=grid->ny;++j) {
            for (decltype(grid->nz) l=0;l!=grid->nz;++l) {
                // get physical position
                vec3 pos {i*lx/(grid->nx-1) + grid->x_min,
                    j*ly/(grid->ny-1) + grid->y_min,
                    l*lz/(grid->nz-1) + grid->z_min};
                // get rescaling factor
                double ratio {sqrt(rescal_fact(pos,par))*par->brnd_iso[0]*CGS_U_muGauss/sqrt(3.*b_var)};
                // add anisotropic field to random one
                unsigned long int idx {toolkit::Index3d(grid->nx,grid->ny,grid->nz,i,j,l)};
                grid->fftw_b_kx[idx][0] *= ratio;
                grid->fftw_b_ky[idx][0] *= ratio;
                grid->fftw_b_kz[idx][0] *= ratio;
                grid->fftw_b_kx[idx][1] *= ratio;
                grid->fftw_b_ky[idx][1] *= ratio;
                grid->fftw_b_kz[idx][1] *= ratio;
            } //l
        } //j
    } //i
    // execute DFT forward plan
    fftw_execute(grid->fftw_px_fw);
    fftw_execute(grid->fftw_py_fw);
    fftw_execute(grid->fftw_pz_fw);
    // PHASE III
    // RE-ORTHOGONALIZING IN FOURIER SPACE
    // Gram-Schmidt process
    for (decltype(grid->nx) i=0;i!=grid->nx;++i) {
        for (decltype(grid->ny) j=0;j!=grid->ny;++j) {
            for (decltype(grid->nz) l=0;l!=grid->nz;++l) {
                unsigned long int idx {toolkit::Index3d(grid->nx,grid->ny,grid->nz,i,j,l)};
                // FFT expects up to n/2 positive while n/2 to n negative
                // physical k in cgs dimension
                vec3 tmp_k {i/(lx/CGS_U_kpc),
                    j/(ly/CGS_U_kpc),
                    l/(lz/CGS_U_kpc)};
                if(i>=grid->nx/2) tmp_k.x -= grid->nx/(lx/CGS_U_kpc);
                if(j>=grid->ny/2) tmp_k.y -= grid->ny/(ly/CGS_U_kpc);
                if(l>=grid->nz/2) tmp_k.z -= grid->nz/(lz/CGS_U_kpc);

                const vec3 tmp_b_re {grid->fftw_b_kx[idx][0],grid->fftw_b_ky[idx][0],grid->fftw_b_kz[idx][0]};
                const vec3 tmp_b_im {grid->fftw_b_kx[idx][1],grid->fftw_b_ky[idx][1],grid->fftw_b_kz[idx][1]};
                
                const vec3 free_b_re {gramschmidt(tmp_k,tmp_b_re)};
                const vec3 free_b_im {gramschmidt(tmp_k,tmp_b_im)};
                
                grid->fftw_b_kx[idx][0] = free_b_re.x;
                grid->fftw_b_ky[idx][0] = free_b_re.y;
                grid->fftw_b_kz[idx][0] = free_b_re.z;
                grid->fftw_b_kx[idx][1] = free_b_im.x;
                grid->fftw_b_ky[idx][1] = free_b_im.y;
                grid->fftw_b_kz[idx][1] = free_b_im.z;
                
                if(i==0 and j==0 and l==0) {
                    grid->fftw_b_kx[idx][0] = 0.;
                    grid->fftw_b_ky[idx][0] = 0.;
                    grid->fftw_b_kz[idx][0] = 0.;
                    grid->fftw_b_kx[idx][1] = 0.;
                    grid->fftw_b_ky[idx][1] = 0.;
                    grid->fftw_b_kz[idx][1] = 0.;
                }
            }// l
        }// j
    }// i
    // execute DFT backward plan
    fftw_execute(grid->fftw_px_bw);
    fftw_execute(grid->fftw_py_bw);
    fftw_execute(grid->fftw_pz_bw);
    // get real elements, use auxiliary function
    complex2real(grid->fftw_b_kx, grid->fftw_b_x, grid->full_size);
    complex2real(grid->fftw_b_ky, grid->fftw_b_y, grid->full_size);
    complex2real(grid->fftw_b_kz, grid->fftw_b_z, grid->full_size);
    // according to FFTW3 manual
    // transform forward followed by backword scale up array by nx*ny*nz
    for(unsigned long int i=0;i!=grid->full_size;++i){
        grid->fftw_b_x[i] /= grid->full_size;
        grid->fftw_b_y[i] /= grid->full_size;
        grid->fftw_b_z[i] /= grid->full_size;
    }

#ifndef NDEBUG
    // check RMS
    b_var = toolkit::Variance(grid->fftw_b_x, grid->full_size);
    b_var +=toolkit::Variance(grid->fftw_b_y, grid->full_size);
    b_var +=toolkit::Variance(grid->fftw_b_z, grid->full_size);
    cout<< "BRND: Numerical RMS: "<<sqrt(b_var)/CGS_U_muGauss<<" microG"<<endl;
    // check average divergence
    double div {0};
    for(decltype(grid->nx) i=2;i!=grid->nx-2;++i) {
        for(decltype(grid->ny) j=2;j!=grid->ny-2;++j) {
            for(decltype(grid->nz) k=2;k!=grid->nz-2;++k) {
                unsigned long int idxf {toolkit::Index3d(grid->nx,grid->ny,grid->nz,i+1,j,k)};
                unsigned long int idxb {toolkit::Index3d(grid->nx,grid->ny,grid->nz,i-1,j,k)};
                unsigned long int idyf {toolkit::Index3d(grid->nx,grid->ny,grid->nz,i,j+1,k)};
                unsigned long int idyb {toolkit::Index3d(grid->nx,grid->ny,grid->nz,i,j-1,k)};
                unsigned long int idzf {toolkit::Index3d(grid->nx,grid->ny,grid->nz,i,j,k+1)};
                unsigned long int idzb {toolkit::Index3d(grid->nx,grid->ny,grid->nz,i,j,k-1)};
                div += fabs( (grid->fftw_b_x[idxf]-grid->fftw_b_x[idxb]) + (grid->fftw_b_y[idyf]-grid->fftw_b_y[idyb]) + (grid->fftw_b_z[idzf]-grid->fftw_b_z[idzb]) );
            }
        }
    }
    cout<<"BRND: Averaged divergence: "<<div/(CGS_U_muGauss*grid->nx*grid->ny*grid->nz)<<" microG/grid"<<endl;
#endif
}

// get real components from fftw_complex arrays
void Brnd_iso::complex2real(const fftw_complex *input,double *output,const unsigned long int &size) {
    for(unsigned long int i=0;i!=size;++i){
        output[i] = input[i][0];
    }
}

// Gram-Schimdt, rewritten using Healpix vec3 library
// tiny error caused by machine is inevitable
vec3 Brnd_iso::gramschmidt(const vec3 &k,const vec3 &b){
    if(k.Length()==0 or b.Length()==0){
        return vec3 {0,0,0};
    }
    vec3 b_free {b.x - k.x*dotprod(k,b)/k.SquaredLength(),
        b.y - k.y*dotprod(k,b)/k.SquaredLength(),
        b.z - k.z*dotprod(k,b)/k.SquaredLength()};
    
    b_free = toolkit::versor(b_free)*b.Length();
    return b_free;
}


// global anisotropic turbulent field
double Brnd_anig::anisotropy(const vec3 &pos,vec3 &H,Pond *par,Breg *breg,Grid_breg *gbreg){
    // the simplest case, const.
    return par->brnd_anig[0];
    // H, direction of anisotropy
    H = toolkit::versor(breg->get_breg(pos,par,gbreg));
}

void Brnd_anig::write_grid_ani(Pond *par, Breg *breg, Grid_breg *gbreg, Grid_brnd *grid){
    // PHASE I
    // GENERATE GAUSSIAN RANDOM FROM SPECTRUM
    // initialize random seed
    gsl_rng *r{gsl_rng_alloc(gsl_rng_taus)};
    gsl_rng_set(r, toolkit::random_seed());
    // start Fourier space filling
    double lx {grid->x_max-grid->x_min};
    double ly {grid->y_max-grid->y_min};
    double lz {grid->z_max-grid->z_min};
    for (decltype(grid->nx) i=0;i!=grid->nx;++i) {
        for (decltype(grid->ny) j=0;j!=grid->ny;++j) {
            for (decltype(grid->nz) l=0;l!=grid->nz;++l) {
                unsigned long int idx {toolkit::Index3d(grid->nx,grid->ny,grid->nz,i,j,l)};
                // FFT expects up to n/2 positive while n/2 to n negative
                // physical k in 1/kpc dimension
                double kx {i/(lx/CGS_U_kpc)};
                double ky {j/(ly/CGS_U_kpc)};
                double kz {l/(lz/CGS_U_kpc)};
                if(i>=grid->nx/2) kx -= grid->nx/(lx/CGS_U_kpc);
                if(j>=grid->ny/2) ky -= grid->ny/(ly/CGS_U_kpc);
                if(l>=grid->nz/2) kz -= grid->nz/(lz/CGS_U_kpc);
                const double k {sqrt(kx*kx + ky*ky + kz*kz)};
                // physical dk^3
                const double dk3 {pow(CGS_U_kpc,3.)/(lx*ly*lz)};
                // simpson's rule
                double element {2.*b_spec(k,par)/3.};
                const double halfdk {0.5*sqrt(pow(CGS_U_kpc/lx,2)+pow(CGS_U_kpc/ly,2)+pow(CGS_U_kpc/lz,2))};
                element += b_spec(k+halfdk,par)/6.;
                element += b_spec(k-halfdk,par)/6.;
                // amplitude, dividing by two because equal allocation to Re and Im parts
                const double sigma {sqrt(element*dk3/2.0)};
                grid->fftw_b_kx[idx][0] = gsl_ran_gaussian(r,sigma);
                grid->fftw_b_ky[idx][0] = gsl_ran_gaussian(r,sigma);
                grid->fftw_b_kz[idx][0] = gsl_ran_gaussian(r,sigma);
                grid->fftw_b_kx[idx][1] = gsl_ran_gaussian(r,sigma);
                grid->fftw_b_ky[idx][1] = gsl_ran_gaussian(r,sigma);
                grid->fftw_b_kz[idx][1] = gsl_ran_gaussian(r,sigma);
                if(i==0 and j==0 and l==0) {
                    grid->fftw_b_kx[idx][0] = 0.;
                    grid->fftw_b_ky[idx][0] = 0.;
                    grid->fftw_b_kz[idx][0] = 0.;
                    grid->fftw_b_kx[idx][1] = 0.;
                    grid->fftw_b_ky[idx][1] = 0.;
                    grid->fftw_b_kz[idx][1] = 0.;
                }
            }// l
        }// j
    }// i
    // free random memory
    gsl_rng_free(r);
    // no Hermiticity fixing, complex 2 complex
    // execute DFT backward plan
    fftw_execute(grid->fftw_px_bw);
    fftw_execute(grid->fftw_py_bw);
    fftw_execute(grid->fftw_pz_bw);
    // PHASE II
    // RESCALING FIELD PROFILE IN REAL SPACE
    // isotropicity has been checked
    double b_var {toolkit::Variance(grid->fftw_b_kx[0], grid->full_size)};
    for (decltype(grid->nx) i=0;i!=grid->nx;++i) {
        for (decltype(grid->ny) j=0;j!=grid->ny;++j) {
            for (decltype(grid->nz) l=0;l!=grid->nz;++l) {
                // get physical position
                vec3 pos {i*lx/(grid->nx-1) + grid->x_min,
                    j*ly/(grid->ny-1) + grid->y_min,
                    l*lz/(grid->nz-1) + grid->z_min};
                // assemble b_Re and b_Im
                unsigned long int idx {toolkit::Index3d(grid->nx,grid->ny,grid->nz,i,j,l)};
                vec3 b_re {grid->fftw_b_kx[idx][0],grid->fftw_b_ky[idx][0],grid->fftw_b_kz[idx][0]};
                vec3 b_im {grid->fftw_b_kx[idx][1],grid->fftw_b_ky[idx][1],grid->fftw_b_kz[idx][1]};
                // get rescaling factor
                double ratio {sqrt(rescal_fact(pos,par))*par->brnd_iso[0]*CGS_U_muGauss/sqrt(3.*b_var)};
                // impose anisotropy
                vec3 H_versor;
                double rho {anisotropy(pos,H_versor,par,breg,gbreg)};
#ifndef NDEBUG
                if(rho<0. or rho>1.){
                    cerr<<"ERR:"<<__FILE__
                    <<" : in function "<<__func__<<endl
                    <<" at line "<<__LINE__<<endl
                    <<"WRONG VALUE"<<endl;
                    exit(1);
                }
#endif
                if(H_versor.Length()!=0){
                    vec3 b_re_par {H_versor*dotprod(H_versor,b_re)};
                    vec3 b_re_perp {b_re - b_re_par};
                    b_re = toolkit::versor(b_re_par*rho + b_re_perp*(1-rho))*(b_re.Length()*ratio);
                    vec3 b_im_par {H_versor*dotprod(H_versor,b_im)};
                    vec3 b_im_perp {b_im - b_im_par};
                    b_im = toolkit::versor(b_im_par*rho + b_im_perp*(1-rho))*(b_im.Length()*ratio);
                }
                // add anisotropic field to random one
                grid->fftw_b_kx[idx][0] = b_re.x;
                grid->fftw_b_ky[idx][0] = b_re.y;
                grid->fftw_b_kz[idx][0] = b_re.z;
                grid->fftw_b_kx[idx][1] = b_im.x;
                grid->fftw_b_ky[idx][1] = b_im.y;
                grid->fftw_b_kz[idx][1] = b_im.z;
            } //l
        } //j
    } //i
    // execute DFT forward plan
    fftw_execute(grid->fftw_px_fw);
    fftw_execute(grid->fftw_py_fw);
    fftw_execute(grid->fftw_pz_fw);
    // PHASE III
    // RE-ORTHOGONALIZING IN FOURIER SPACE
    // Gram-Schmidt process
    for (decltype(grid->nx) i=0;i!=grid->nx;++i) {
        for (decltype(grid->ny) j=0;j!=grid->ny;++j) {
            for (decltype(grid->nz) l=0;l!=grid->nz;++l) {
                unsigned long int idx {toolkit::Index3d(grid->nx,grid->ny,grid->nz,i,j,l)};
                // FFT expects up to n/2 positive while n/2 to n negative
                // physical k in 1/kpc
                vec3 tmp_k {i/(lx/CGS_U_kpc),
                    j/(ly/CGS_U_kpc),
                    l/(lz/CGS_U_kpc)};
                if(i>=grid->nx/2) tmp_k.x -= grid->nx/(lx/CGS_U_kpc);
                if(j>=grid->ny/2) tmp_k.y -= grid->ny/(ly/CGS_U_kpc);
                if(l>=grid->nz/2) tmp_k.z -= grid->nz/(lz/CGS_U_kpc);
                
                const vec3 tmp_b_re {grid->fftw_b_kx[idx][0],grid->fftw_b_ky[idx][0],grid->fftw_b_kz[idx][0]};
                const vec3 tmp_b_im {grid->fftw_b_kx[idx][1],grid->fftw_b_ky[idx][1],grid->fftw_b_kz[idx][1]};
                
                const vec3 free_b_re {gramschmidt(tmp_k,tmp_b_re)};
                const vec3 free_b_im {gramschmidt(tmp_k,tmp_b_im)};
                
                grid->fftw_b_kx[idx][0] = free_b_re.x;
                grid->fftw_b_ky[idx][0] = free_b_re.y;
                grid->fftw_b_kz[idx][0] = free_b_re.z;
                grid->fftw_b_kx[idx][1] = free_b_im.x;
                grid->fftw_b_ky[idx][1] = free_b_im.y;
                grid->fftw_b_kz[idx][1] = free_b_im.z;
                
                if(i==0 and j==0 and l==0) {
                    grid->fftw_b_kx[idx][0] = 0.;
                    grid->fftw_b_ky[idx][0] = 0.;
                    grid->fftw_b_kz[idx][0] = 0.;
                    grid->fftw_b_kx[idx][1] = 0.;
                    grid->fftw_b_ky[idx][1] = 0.;
                    grid->fftw_b_kz[idx][1] = 0.;
                }
            }// l
        }// j
    }// i
    // execute DFT backward plan
    fftw_execute(grid->fftw_px_bw);
    fftw_execute(grid->fftw_py_bw);
    fftw_execute(grid->fftw_pz_bw);
    // get real elements, use auxiliary function
    complex2real(grid->fftw_b_kx, grid->fftw_b_x, grid->full_size);
    complex2real(grid->fftw_b_ky, grid->fftw_b_y, grid->full_size);
    complex2real(grid->fftw_b_kz, grid->fftw_b_z, grid->full_size);
    // according to FFTW3 manual
    // transform forward followed by backword scale up array by nx*ny*nz
    for(unsigned long int i=0;i!=grid->full_size;++i){
        grid->fftw_b_x[i] /= grid->full_size;
        grid->fftw_b_y[i] /= grid->full_size;
        grid->fftw_b_z[i] /= grid->full_size;
    }
    
#ifndef NDEBUG
    // check RMS
    b_var = toolkit::Variance(grid->fftw_b_x, grid->full_size);
    b_var +=toolkit::Variance(grid->fftw_b_y, grid->full_size);
    b_var +=toolkit::Variance(grid->fftw_b_z, grid->full_size);
    cout<< "BRND: Numerical RMS: "<<sqrt(b_var)/CGS_U_muGauss<<" microG"<<endl;
    // check average divergence
    double div=0;
    for(decltype(grid->nx) i=2;i!=grid->nx-2;++i) {
        for(decltype(grid->ny) j=2;j!=grid->ny-2;++j) {
            for(decltype(grid->nz) k=2;k!=grid->nz-2;++k) {
                unsigned long int idxf {toolkit::Index3d(grid->nx,grid->ny,grid->nz,i+1,j,k)};
                unsigned long int idxb {toolkit::Index3d(grid->nx,grid->ny,grid->nz,i-1,j,k)};
                unsigned long int idyf {toolkit::Index3d(grid->nx,grid->ny,grid->nz,i,j+1,k)};
                unsigned long int idyb {toolkit::Index3d(grid->nx,grid->ny,grid->nz,i,j-1,k)};
                unsigned long int idzf {toolkit::Index3d(grid->nx,grid->ny,grid->nz,i,j,k+1)};
                unsigned long int idzb {toolkit::Index3d(grid->nx,grid->ny,grid->nz,i,j,k-1)};
                div += fabs( (grid->fftw_b_x[idxf]-grid->fftw_b_x[idxb]) + (grid->fftw_b_y[idyf]-grid->fftw_b_y[idyb]) + (grid->fftw_b_z[idzf]-grid->fftw_b_z[idzb]) );
            }
        }
    }
    cout<<"BRND: Averaged divergence: "<<div/(CGS_U_muGauss*grid->nx*grid->ny*grid->nz)<<" microG/grid"<<endl;
#endif
}

//local anisotropic turbulent field
void Brnd_anil::anisotropy(double *tensor, Pond *par, Breg *breg, Grid_breg *gbreg){
    
}

void Brnd_anil::write_grid_ani(Pond *par, Breg *breg, Grid_breg *gbreg, Grid_brnd *grid){
    // initialize random seed
    gsl_rng *r {gsl_rng_alloc(gsl_rng_taus)};
    gsl_rng_set(r, toolkit::random_seed());
    
    for (decltype(grid->nx) i=0;i!=grid->nx;++i) {
        for (decltype(grid->ny) j=0;j!=grid->ny;++j) {
            for (decltype(grid->nz) l=0;l!=grid->nz;++l) {
                unsigned long int idx {toolkit::Index3d(grid->nx,grid->ny,grid->nz,i,j,l)};
                // FFT expects up to n/2 positive while n/2 to n negative
                // physical k in 1/kpc dimension
                double kx {i/(lx/CGS_U_kpc)};
                double ky {j/(ly/CGS_U_kpc)};
                double kz {l/(lz/CGS_U_kpc)};
                if(i>=grid->nx/2) kx -= grid->nx/(lx/CGS_U_kpc);
                if(j>=grid->ny/2) ky -= grid->ny/(ly/CGS_U_kpc);
                if(l>=grid->nz/2) kz -= grid->nz/(lz/CGS_U_kpc);
                const double k {sqrt(kx*kx + ky*ky + kz*kz)};
                // physical dk^3
                const double dk3 {pow(CGS_U_kpc,3.)/(lx*ly*lz)};
                
                // cheating simpson evaluation
                const double halfdk {0.5*sqrt(pow(CGS_U_kpc/lx,2)+pow(CGS_U_kpc/ly,2)+pow(CGS_U_kpc/lz,2))};
                double element {2.*b_spec(k,par)/3.};
                element += b_spec(k+halfdk,par)/6.;
                element += b_spec(k-halfdk,par)/6.;
                
                /* UNFINISHED */
                }
            }// l
        }// j
    }// i
    // free random memory
    gsl_rng_free(r);
    
    // execute DFT plan
    fftw_execute(grid->fftw_px);
    fftw_execute(grid->fftw_py);
    fftw_execute(grid->fftw_pz);
    
    // get real elements, use auxiliary function
    complex2real(grid->fftw_b_kx, grid->fftw_b_x, grid->full_size);
    complex2real(grid->fftw_b_ky, grid->fftw_b_y, grid->full_size);
    complex2real(grid->fftw_b_kz, grid->fftw_b_z, grid->full_size);
    
#ifndef NDEBUG
    double b_var;
    b_var = toolkit::Variance(grid->fftw_b_x, grid->full_size);
    b_var +=toolkit::Variance(grid->fftw_b_y, grid->full_size);
    b_var +=toolkit::Variance(grid->fftw_b_z, grid->full_size);
    cout<< "BRND: Numerical RMS: "<<sqrt(b_var)/CGS_U_muGauss<<" microG"<<endl;
    // check average divergence
    double div=0;
    for(decltype(grid->nx) i=2;i!=grid->nx-2;++i) {
        for(decltype(grid->ny) j=2;j!=grid->ny-2;++j) {
            for(decltype(grid->nz) k=2;k!=grid->nz-2;++k) {
                auto idxf = toolkit::Index3d(grid->nx,grid->ny,grid->nz,i+1,j,k);
                auto idxb = toolkit::Index3d(grid->nx,grid->ny,grid->nz,i-1,j,k);
                auto idyf = toolkit::Index3d(grid->nx,grid->ny,grid->nz,i,j+1,k);
                auto idyb = toolkit::Index3d(grid->nx,grid->ny,grid->nz,i,j-1,k);
                auto idzf = toolkit::Index3d(grid->nx,grid->ny,grid->nz,i,j,k+1);
                auto idzb = toolkit::Index3d(grid->nx,grid->ny,grid->nz,i,j,k-1);
                div += fabs( (grid->fftw_b_x[idxf]-grid->fftw_b_x[idxb]) + (grid->fftw_b_y[idyf]-grid->fftw_b_y[idyb]) + (grid->fftw_b_z[idzf]-grid->fftw_b_z[idzb]) );
            }
        }
    }
    cout<<"BRND: Averaged divergence: "<<div/(CGS_U_muGauss*grid->nx*grid->ny*grid->nz)<<" microG/grid"<<endl;
#endif
}

// END
