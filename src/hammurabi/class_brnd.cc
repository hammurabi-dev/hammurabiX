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

/* base class */

vec3 Brnd::get_brnd(const vec3 &pos, Pond *par, Grid_brnd *grid){
    // if no specific random field model is called
    // base class will return zero vector
    return vec3(0.,0.,0.);
}

double Brnd::b_spec(const double &k, Pond *par){
    //units fixing
    const double p0 = par->bgrnd[0]*pow(CGS_U_muGauss,2)*CGS_U_kpc;
    const double k0 = par->bgrnd[1]/CGS_U_kpc;
    const double a0 = par->bgrnd[2];
    const double unit = 1./(4*CGS_U_pi*k*k);
    // avoid nan
    if(k<=0.){
        return 0.;
    }
    // modified Giacalone-Jokipii model
    //const double P = 2*p0*pow(k/k0,a0)/(1.+pow(k/k0,a0-a1));
    
    // power law
    
    double P = 0.;
    if(k>k0){
        P = p0*pow(k/k0,a0);
    }
    return P*unit;
}

// theoretical integration from k_min in simulation to k=1/1pc of P(k)pi^2
double Brnd::get_energy(Pond *par, Grid_brnd *grid){
    // minimum physical k in DFT
    const double k_start = sqrt(pow(1./grid->lx,2)+pow(1./grid->ly,2)+pow(1./grid->lz,2));
    // end at 1pc
    const double k_end = 1000./CGS_U_kpc;
    // steps
    const unsigned int N = 3000;
    const double k_lapse = (k_end-k_start)/N;
    double result = 0.;
    // using simpson's rule
    for(unsigned int i=0;i!=N;++i){
        double k0 = k_start + i*k_lapse;
        double k1 = k0 + 0.5*k_lapse;
        double k2 = k0 + k_lapse;
        result += (k0*k0*b_spec(k0,par)+4.*k1*k1*b_spec(k1,par)+k2*k2*b_spec(k2,par))*k_lapse/12.;
    }
    return result;
}

// theoretical integration from k_max in simulation to k=1/1pc of P(k)pi^2
double Brnd::get_missing(Pond *par, Grid_brnd *grid){
    // maximum physical k in DFT
    double k_start = sqrt(3.)*(grid->nx/grid->lx)/2.;
    // end at 1pc
    const double k_end = 1000./CGS_U_kpc;
    // steps
    const unsigned int N = 1000;
    const double k_lapse = (k_end-k_start)/N;
    double result = 0.;
    // using simpson's rule
    for(unsigned int i=0;i!=N;++i){
        double k0 = k_start + i*k_lapse;
        double k1 = k0 + 0.5*k_lapse;
        double k2 = k0 + k_lapse;
        result += (k0*k0*b_spec(k0,par)+4.*k1*k1*b_spec(k1,par)+k2*k2*b_spec(k2,par))*k_lapse/12.;
    }
    return result;
}

double Brnd::rescal_fact(const vec3 &pos, Pond *par){
    const double r_cyl = sqrt(pos.x*pos.x+pos.y*pos.y)/CGS_U_kpc;
    const double z = fabs(pos.z/CGS_U_kpc);
    const double r0 = par->bscal[0];
    const double z0 = par->bscal[1];
    if(r_cyl==0. or z==0){return 1.;}
    else{
    	return exp(-r_cyl/r0)*exp(-z/z0);
    }
}

vec3 Brnd::read_grid(const vec3 &pos, Grid_brnd *grid, Pond *par){
    
    decltype(grid->nx) xl, yl, zl;
    
    double tmp = (grid->nx-1)*(pos.x/grid->lx + 0.5);
    // if use solar-centric grid
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
        <<"WRONG VALUE"<<endl;
        exit(1);
    }
#endif
    vec3 b_vec3;
    //trilinear interpolation
    if (xl+1<grid->nx and yl+1<grid->ny and zl+1<grid->nz) {
        //for storing vectors
        vec3 i1,i2,j1,j2,w1,w2;
        
        //interpolate along z direction, there are four interpolated vectors
        auto idx1 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl,yl,zl);
        auto idx2 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl,yl,zl+1);
        i1 = vec3(grid->fftw_b_x[idx1]*(1.-zd) + grid->fftw_b_x[idx2]*zd,
                  grid->fftw_b_y[idx1]*(1.-zd) + grid->fftw_b_y[idx2]*zd,
                  grid->fftw_b_z[idx1]*(1.-zd) + grid->fftw_b_z[idx2]*zd );
        
        idx1 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl,yl+1,zl);
        idx2 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl,yl+1,zl+1);
        i2 = vec3(grid->fftw_b_x[idx1]*(1.-zd) + grid->fftw_b_x[idx2]*zd,
                  grid->fftw_b_y[idx1]*(1.-zd) + grid->fftw_b_y[idx2]*zd,
                  grid->fftw_b_z[idx1]*(1.-zd) + grid->fftw_b_z[idx2]*zd);
        
        idx1 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl+1,yl,zl);
        idx2 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl+1,yl,zl+1);
        j1 = vec3(grid->fftw_b_x[idx1]*(1.-zd) + grid->fftw_b_x[idx2]*zd,
                  grid->fftw_b_y[idx1]*(1.-zd) + grid->fftw_b_y[idx2]*zd,
                  grid->fftw_b_z[idx1]*(1.-zd) + grid->fftw_b_z[idx2]*zd);
        
        idx1 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl+1,yl+1,zl);
        idx2 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl+1,yl+1,zl+1);
        j2 = vec3(grid->fftw_b_x[idx1]*(1.-zd) + grid->fftw_b_x[idx2]*zd,
                  grid->fftw_b_y[idx1]*(1.-zd) + grid->fftw_b_y[idx2]*zd,
                  grid->fftw_b_z[idx1]*(1.-zd) + grid->fftw_b_z[idx2]*zd);
        // interpolate along y direction, two interpolated vectors
        w1 = i1*(1.-yd) + i2*yd;
        w2 = j1*(1.-yd) + j2*yd;
        // interpolate along x direction
        b_vec3 = w1*(1.-xd) + w2*xd;
    }
    // on the boundary
    else {
        auto idx = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl,yl,zl);
        b_vec3 = vec3(grid->fftw_b_x[idx],grid->fftw_b_y[idx],grid->fftw_b_z[idx]);
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
void Brnd::write_grid(Pond *par,Grid_brnd *grid){
    cerr<<"WAR:"<<__FILE__
    <<" : in function "<<__func__<<endl
    <<" at line "<<__LINE__<<endl
    <<"DYNAMIC BINDING FAILURE"<<endl;
    exit(1);
}
void Brnd::write_grid_plus(Pond *par, Breg *breg, Grid_breg *gbreg, Grid_brnd *gbrnd){
    cerr<<"WAR:"<<__FILE__
    <<" : in function "<<__func__<<endl
    <<" at line "<<__LINE__<<endl
    <<"DYNAMIC BINDING FAILURE"<<endl;
    exit(1);
}

/* gaussian random field */

Bgrnd::Bgrnd(Pond *par, Grid_brnd *grid){
    get_energy_rslt = get_energy(par,grid);
    get_missing_rslt = get_missing(par,grid);
}

vec3 Bgrnd::get_brnd(const vec3 &pos, Pond *par, Grid_brnd *grid){
    // interpolate written grid to given position
    // check if you have called ::write_grid
    return read_grid(pos,grid,par);
}

void Bgrnd::write_grid(Pond *par, Grid_brnd *grid){
    // initialize random seed
    gsl_rng *r = gsl_rng_alloc(gsl_rng_taus);
    // from tools
    gsl_rng_set(r, toolkit::random_seed());
    
    for (decltype(grid->nx) i=0;i!=grid->nx;++i) {
        for (decltype(grid->ny) j=0;j!=grid->ny;++j) {
            for (decltype(grid->nz) l=0;l!=grid->nz;++l) {
                
                // point out where we are
                auto idx = toolkit::Index3d(grid->nx,grid->ny,grid->nz,i,j,l);
                
                // FFT expects up to n/2 positive while n/2 to n negative
                // physical k in cgs dimension
                double kx = double(i)/grid->lx;
                double ky = double(j)/grid->ly;
                double kz = double(l)/grid->lz;
                if(i>=grid->nx/2) kx -= double(grid->nx)/grid->lx;
                if(j>=grid->ny/2) ky -= double(grid->ny)/grid->ly;
                if(l>=grid->nz/2) kz -= double(grid->nz)/grid->lz;
                
                const double k = sqrt(kx*kx + ky*ky + kz*kz);
                
                // kx=ky=0 elements, since sqrt(kx^2+ky^2) in above
                if(i==0 and j==0) {
                    grid->fftw_b_kx[idx][0]=0.; grid->fftw_b_kx[idx][1]=0.;
                    grid->fftw_b_ky[idx][0]=0.; grid->fftw_b_ky[idx][1]=0.;
                    grid->fftw_b_kz[idx][0]=0.; grid->fftw_b_kz[idx][1]=0.;
                    continue;
                }
                
                // now we express k space base {e1,e2,\hat{k}} by {x,y,z} base
                // the following form is one of infinite combinations
                vec3 ep,em; // ep and em are actually e+ and e- in draft
                ep.x = kx*kz/(k*sqrt(kx*kx+ky*ky));
                ep.y = ky*kz/(k*sqrt(kx*kx+ky*ky));
                ep.z = -sqrt(kx*kx+ky*ky)/k;
                
                em.x = -ky/sqrt(kx*kx+ky*ky);
                em.y = kx/sqrt(kx*kx+ky*ky);
                em.z = 0.;
                // the idea of ep and em is that, following divergenceless requirement
                // we decompose b(k) into two vector modes
                // the pre-normalized bases of two vector modes are ep and em
                // bp and bm are two complex coeffcients in two directions
                // which satisfy <bp^2 + bm^2> = P(k)(2\pi)^3
                // in terms of Re and Im components in bp and bm, they can be just random
                // no constraints on these two degrees of freedom
                
                // ratio of plus mode variance
                // rho should go just random within [0,1]
                const double rho = gsl_rng_uniform(r);
                // angles
                // angles also go random within [0, 2pi]
                const double angp = 2*CGS_U_pi*gsl_rng_uniform(r);
                const double angm = 2*CGS_U_pi*gsl_rng_uniform(r);
                // physical dk^3
                const double dk3 = 1./(grid->lx*grid->ly*grid->lz);
                
                // cheating simpson evaluation
                double element = 2.*b_spec(k,par)/3.;
                const double halfdk = 0.5*sqrt(pow(1./grid->lx,2)+pow(1./grid->ly,2)+pow(1./grid->lz,2));
                element += b_spec(k+halfdk,par)/6.;
                element += b_spec(k-halfdk,par)/6.;
                // sigma+ and sigma-
                const double sigmap = sqrt(element*dk3*rho);
                const double sigmam = sqrt(element*dk3*(1-rho));
                
                // amplitude sampling
                const double bp = gsl_ran_gaussian(r,sigmap);
                const double bm = gsl_ran_gaussian(r,sigmam);
                
                grid->fftw_b_kx[idx][0] = bp*cos(angp)*ep.x + bm*cos(angm)*em.x;
                grid->fftw_b_ky[idx][0] = bp*cos(angp)*ep.y + bm*cos(angm)*em.y;
                grid->fftw_b_kz[idx][0] = bp*cos(angp)*ep.z;
                
                // fix zero imaginary elements
                if((i==0 or i==grid->nx/2) and (j==0 or j==grid->ny/2) and (l==0 or l==grid->nz/2)) {
                    grid->fftw_b_kx[idx][1] = 0.;
                    grid->fftw_b_ky[idx][1] = 0.;
                    grid->fftw_b_kz[idx][1] = 0.;
                }else {
                    grid->fftw_b_kx[idx][1] = bp*sin(angp)*ep.x + bm*sin(angm)*em.x;
                    grid->fftw_b_ky[idx][1] = bp*sin(angp)*ep.y + bm*sin(angm)*em.y;
                    grid->fftw_b_kz[idx][1] = bp*sin(angp)*ep.z;
                }
            }// l
        }// j
    }// i
    // free random memory
    gsl_rng_free(r);
    
    // fix Hermiticity, use auxiliary function
    hermiticity(grid->fftw_b_kx, grid->nx, grid->ny, grid->nz);
    hermiticity(grid->fftw_b_ky, grid->nx, grid->ny, grid->nz);
    hermiticity(grid->fftw_b_kz, grid->nx, grid->ny, grid->nz);
    
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
    cout<<"BRND: Analytical RMS: "<<sqrt(get_energy_rslt*8.*CGS_U_pi)/CGS_U_muGauss<<" microG"<<endl;
    cout<<"BRND: Missing RMS: "<<sqrt(get_missing_rslt*8.*CGS_U_pi)/CGS_U_muGauss<<" microG"<<endl;
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

// hermiticity fixer
void Bgrnd::hermiticity(fftw_complex *a_c,const unsigned int &N1,const unsigned int &N2,const unsigned int &N3){
    /* hermiticity: complete volume for k < N3/2 */
    for(unsigned int i=1;i!=N1;i++){
        for(unsigned int j=1;j!=N2;j++){
            for(unsigned int k=1;k!=N3/2;k++){
                auto idx = toolkit::Index3d(N1,N2,N3,i,j,k);
                auto idxc = toolkit::Index3d(N1,N2,N3,N1-i,N2-j,N3-k);
                a_c[idxc][0] = a_c[idx][0];
                a_c[idxc][1] = -a_c[idx][1];
            }
        }
    }
    /* hermiticity: plane for k = N3/2, j<N2/2 */
    for(unsigned int i=1;i!=N1;i++){
        for(unsigned int j=1;j!=N2/2;j++){
            unsigned int k=N3/2;
            auto idx = toolkit::Index3d(N1,N2,N3,i,j,k);
            auto idxc = toolkit::Index3d(N1,N2,N3,N1-i,N2-j,k);
            a_c[idxc][0] = a_c[idx][0];
            a_c[idxc][1] = -a_c[idx][1];
        }
    }
    /* hermiticity: plane for k = 0, j<N2/2 */
    for(unsigned int i=1;i!=N1;i++){
        for(unsigned int j=1;j!=N2/2;j++){
            unsigned int k=0;
            auto idx = toolkit::Index3d(N1,N2,N3,i,j,k);
            auto idxc = toolkit::Index3d(N1,N2,N3,N1-i,N2-j,k);
            a_c[idxc][0] = a_c[idx][0];
            a_c[idxc][1] = -a_c[idx][1];
        }
    }
    /* hermiticity: plane for j = 0, k<N3/2 */
    for(unsigned int i=1;i!=N1;i++){
        for(unsigned int k=1;k!=N3/2;k++){
            unsigned int j=0;
            auto idx = toolkit::Index3d(N1,N2,N3,i,j,k);
            auto idxc = toolkit::Index3d(N1,N2,N3,N1-i,j,N3-k);
            a_c[idxc][0] = a_c[idx][0];
            a_c[idxc][1] = -a_c[idx][1];
        }
    }
    /* hermiticity: plane for i = 0, k<N3/2 */
    for(unsigned int j=1;j!=N2;j++){
        for(unsigned int k=1;k!=N3/2;k++){
            unsigned int i=0;
            auto idx = toolkit::Index3d(N1,N2,N3,i,j,k);
            auto idxc = toolkit::Index3d(N1,N2,N3,i,N2-j,N3-k);
            a_c[idxc][0] = a_c[idx][0];
            a_c[idxc][1] = -a_c[idx][1];
        }
    }
    /* hermiticity: line for k = N3/2, j=N2/2 */
    for(unsigned int i=1;i!=N1/2;i++){
        unsigned int j=N2/2;
        unsigned int k=N3/2;
        auto idx = toolkit::Index3d(N1,N2,N3,i,j,k);
        auto idxc = toolkit::Index3d(N1,N2,N3,N1-i,j,k);
        a_c[idxc][0] = a_c[idx][0];
        a_c[idxc][1] = -a_c[idx][1];
    }
    /* hermiticity: line for k = 0, j=N2/2 */
    for(unsigned int i=1;i!=N1/2;i++){
        unsigned int j=N2/2;
        unsigned int k=0;
        auto idx = toolkit::Index3d(N1,N2,N3,i,j,k);
        auto idxc = toolkit::Index3d(N1,N2,N3,N1-i,j,k);
        a_c[idxc][0] = a_c[idx][0];
        a_c[idxc][1] = -a_c[idx][1];
    }
    /* hermiticity: line for k = N3/2, j=0 */
    for(unsigned int i=1;i!=N1/2;i++){
        unsigned int k=N3/2;
        unsigned int j=0;
        auto idx = toolkit::Index3d(N1,N2,N3,i,j,k);
        auto idxc = toolkit::Index3d(N1,N2,N3,N1-i,j,k);
        a_c[idxc][0] = a_c[idx][0];
        a_c[idxc][1] = -a_c[idx][1];
    }
    /* hermiticity: line for k = N3/2, i=0 */
    for(unsigned int j=1;j!=N2/2;j++){
        unsigned int k=N3/2;
        unsigned int i=0;
        auto idx = toolkit::Index3d(N1,N2,N3,i,j,k);
        auto idxc = toolkit::Index3d(N1,N2,N3,i,N2-j,k);
        a_c[idxc][0] = a_c[idx][0];
        a_c[idxc][1] = -a_c[idx][1];
    }
    /* hermiticity: line for k = 0, j=0 */
    for(unsigned int i=1;i!=N1/2;i++){
        unsigned int k=0;
        unsigned int j=0;
        auto idx = toolkit::Index3d(N1,N2,N3,i,j,k);
        auto idxc = toolkit::Index3d(N1,N2,N3,N1-i,j,k);
        a_c[idxc][0] = a_c[idx][0];
        a_c[idxc][1] = -a_c[idx][1];
    }
    /* hermiticity: line for k = 0, i=0 */
    for(unsigned int j=1;j!=N2/2;j++){
        unsigned int k=0;
        unsigned int i=0;
        auto idx = toolkit::Index3d(N1,N2,N3,i,j,k);
        auto idxc = toolkit::Index3d(N1,N2,N3,i,N2-j,k);
        a_c[idxc][0] = a_c[idx][0];
        a_c[idxc][1] = -a_c[idx][1];
    }
    /* hermiticity: line for i = 0, j=0 */
    for(unsigned int k=1;k!=N3/2;k++){
        unsigned int i=0;
        unsigned int j=0;
        auto idx = toolkit::Index3d(N1,N2,N3,i,j,k);
        auto idxc = toolkit::Index3d(N1,N2,N3,i,j,N3-k);
        a_c[idxc][0] = a_c[idx][0];
        a_c[idxc][1] = -a_c[idx][1];
    }
}

// get real components from fftw_complex arrays
void Bgrnd::complex2real(const fftw_complex *input,double *output,const unsigned long int &size) {
    for(unsigned long int i=0;i!=size;++i){
        output[i] = input[i][0];
    }
}

/* gaussian random + anisotropic random */
Bfrnd::Bfrnd(Pond *par, Grid_brnd *grid){
    get_energy_rslt = get_energy(par,grid);
    get_missing_rslt = get_missing(par,grid);
}

void Bfrnd::write_grid_plus(Pond *par, Breg *breg, Grid_breg *gbreg, Grid_brnd *grid){
    // first, evaluate Gaussian random field to grid
    // inherit from Bgrnd
    write_grid(par,grid);
    // second, add anisotropic value
    // initialize random seed
    gsl_rng *r = gsl_rng_alloc(gsl_rng_taus);
    // from tools
    gsl_rng_set(r, toolkit::random_seed());
    // last parameter is beta (the scaling factor wrt regular field strength)
    const double beta = par->bgrnd[par->bgrnd.size()-1];
    // position, regualr field
    vec3 pos, regular;
    for (decltype(grid->nx) i=0;i!=grid->nx;++i) {
        for (decltype(grid->ny) j=0;j!=grid->ny;++j) {
            for (decltype(grid->nz) l=0;l!=grid->nz;++l) {
                // get regular field vector at given grid point
                pos.x = grid->lx*(double(i)/(grid->nx-1) - 0.5);
                pos.y = grid->ly*(double(j)/(grid->ny-1) - 0.5);
                pos.z = grid->lz*(double(l)/(grid->nz-1) - 0.5);
                // if use solar-centric grid
                if(grid->ec_frame){
                    pos.x += par->SunPosition.x;
                }
                regular = breg->get_breg(pos,par,gbreg);
                // add anisotropic field to random one
                auto idx = toolkit::Index3d(grid->nx,grid->ny,grid->nz,i,j,l);
                regular *= sqrt(beta)*gsl_rng_uniform(r);
                grid->fftw_b_x[idx] += regular.x;
                grid->fftw_b_y[idx] += regular.y;
                grid->fftw_b_z[idx] += regular.z;
            }
        }
    }
    // free random memory
    gsl_rng_free(r);
}
// END
