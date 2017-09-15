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
#include "class_fernd.h"
#include "class_fe.h"
#include "cgs_units_file.h"
#include "namespace_toolkit.h"

using namespace std;

/* base class */

double FErnd::get_fernd(const vec3 &pos, Pond *par, Grid_fernd *grid){
    if(grid->read_permission){
        return read_grid(pos,grid,par);
    }
    else{
        // if no specific random field model is called
        // base class will return zero vector
        return 0.;
    }
}

double FErnd::fe_spec(const double &k, Pond *par){
    //units fixing
    const double p0 = par->fernd_iso[0]*CGS_U_kpc;
    const double k0 = par->fernd_iso[1]/CGS_U_kpc;
    const double a0 = par->fernd_iso[2];
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
double FErnd::get_variance(Pond *par, Grid_fernd *grid){
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
        result += (k0*k0*fe_spec(k0,par)+4.*k1*k1*fe_spec(k1,par)+k2*k2*fe_spec(k2,par))*k_lapse/6.;
    }
    return result*4.*CGS_U_pi;
}

// theoretical integration from k_max in simulation to k=1/1pc of P(k)pi^2
double FErnd::get_missing(Pond *par, Grid_fernd *grid){
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
        result += (k0*k0*fe_spec(k0,par)+4.*k1*k1*fe_spec(k1,par)+k2*k2*fe_spec(k2,par))*k_lapse/6.;
    }
    return result*4.*CGS_U_pi;
}

double FErnd::rescal_fact(const vec3 &pos, Pond *par){
    const double r_cyl = sqrt(pos.x*pos.x+pos.y*pos.y)/CGS_U_kpc;
    const double z = fabs(pos.z/CGS_U_kpc);
    const double r0 = par->fernd_scal[0];
    const double z0 = par->fernd_scal[1];
    if(r_cyl==0. or z==0) {return 1.;}
    else{
	return exp(-r_cyl/r0)*exp(-z/z0);
    }
}

double FErnd::read_grid(const vec3 &pos, Grid_fernd *grid, Pond *par){
    
    decltype(grid->nx) xl, yl, zl;
    
    double tmp = (grid->nx-1)*(pos.x/grid->lx + 0.5);
    // if use solar-centric grid
    if(grid->ec_frame) tmp -= (grid->nx-1)*(par->SunPosition.x/grid->lx);
    if (tmp<0 or tmp>grid->nx-1) { return 0.;}
    else xl = floor(tmp);
    const double xd = tmp - xl;
    
    tmp = (grid->ny-1)*(pos.y/grid->ly + 0.5);
    if (tmp<0 or tmp>grid->ny-1) { return 0.;}
    else yl = floor(tmp);
    const double yd = tmp - yl;
    
    tmp = (grid->nz-1)*(pos.z/grid->lz + 0.5);
    if (tmp<0 or tmp>grid->nz-1) { return 0.;}
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
    double density;
    //trilinear interpolation
    if (xl+1<grid->nx and yl+1<grid->ny and zl+1<grid->nz) {
        //for storing vectors
        double i1,i2,j1,j2,w1,w2;
        
        //interpolate along z direction, there are four interpolated vectors
        auto idx1 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl,yl,zl);
        auto idx2 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl,yl,zl+1);
        i1 = grid->fftw_fe[idx1]*(1.-zd) + grid->fftw_fe[idx2]*zd;
        
        idx1 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl,yl+1,zl);
        idx2 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl,yl+1,zl+1);
        i2 = grid->fftw_fe[idx1]*(1.-zd) + grid->fftw_fe[idx2]*zd;
        
        idx1 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl+1,yl,zl);
        idx2 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl+1,yl,zl+1);
        j1 = grid->fftw_fe[idx1]*(1.-zd) + grid->fftw_fe[idx2]*zd;
        
        idx1 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl+1,yl+1,zl);
        idx2 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl+1,yl+1,zl+1);
        j2 = grid->fftw_fe[idx1]*(1.-zd) + grid->fftw_fe[idx2]*zd;
        // interpolate along y direction, two interpolated vectors
        w1 = i1*(1.-yd) + i2*yd;
        w2 = j1*(1.-yd) + j2*yd;
        // interpolate along x direction
        density = w1*(1.-xd) + w2*xd;
    }
    // on the boundary
    else {
        auto idx = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl,yl,zl);
        density = grid->fftw_fe[idx];
    }
    /*
#ifndef NDEBUG
    if (density<0) {
        cerr<<"WAR:"<<__FILE__
        <<" : in function "<<__func__<<endl
        <<" at line "<<__LINE__<<endl
        <<"NEGATIVE FREE ELECTRON DENSITY"<<endl;
        exit(1);
    }
#endif
     */
    return density;
}

// base class does not write out to grid
void FErnd::write_grid(Pond *par,Grid_fernd *grid){
    cerr<<"WAR:"<<__FILE__
    <<" : in function "<<__func__<<endl
    <<" at line "<<__LINE__<<endl
    <<"DYNAMIC BINDING FAILURE"<<endl;
    exit(1);
}

/* gaussian random field */
FEgrnd::FEgrnd(Pond *par, Grid_fernd *grid){
    get_variance_rslt = get_variance(par,grid);
    get_missing_rslt = get_missing(par,grid);
}

double FEgrnd::get_fernd(const vec3 &pos, Pond *par, Grid_fernd *grid){
    // interpolate written grid to given position
    // check if you have called ::write_grid
    return read_grid(pos,grid,par);
}

void FEgrnd::write_grid(Pond *par, Grid_fernd *grid){
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
                
                if(i==0 and j==0 and l==0) {
                    grid->fftw_fe_k[idx][0]=0.;
                    grid->fftw_fe_k[idx][1]=0.;
                    continue;
                }
                // angle
                // angle also go random within [0, 2pi]
                const double ang = 2*CGS_U_pi*gsl_rng_uniform(r);
                // physical dk^3
                const double dk3 = 1./(grid->lx*grid->ly*grid->lz);
                
                // cheating simpson evaluation
                double element = 2.*fe_spec(k,par)/3.;
                const double halfdk = 0.5*sqrt(pow(1./grid->lx,2)+pow(1./grid->ly,2)+pow(1./grid->lz,2));
                element += fe_spec(k+halfdk,par)/6.;
                element += fe_spec(k-halfdk,par)/6.;
                // sigma+ and sigma-
                const double sigma = sqrt(element*dk3);
                
                // amplitude sampling
                const double fe = gsl_ran_gaussian(r,sigma);
                
                grid->fftw_fe_k[idx][0] = fe*cos(ang);
                // fix zero imaginary elements
                if((i==0 or i==grid->nx/2) and (j==0 or j==grid->ny/2) and (l==0 or l==grid->nz/2)) {
                    grid->fftw_fe_k[idx][1] = 0.;
                }else {
                    grid->fftw_fe_k[idx][1] = fe*sin(ang);
                }
            }// l
        }// j
    }// i
    // free random memory
    gsl_rng_free(r);
    
    // fix Hermiticity, use auxiliary function
    hermiticity(grid->fftw_fe_k, grid->nx, grid->ny, grid->nz);
    
    // execute DFT plan
    fftw_execute(grid->fftw_p);
    
    // get real elements, use auxiliary function
    complex2real(grid->fftw_fe_k, grid->fftw_fe, grid->full_size);
    
#ifndef NDEBUG
    double fe_var;
    fe_var = toolkit::Variance(grid->fftw_fe, grid->full_size);
    cout<< "FERND: Numerical RMS: "<<sqrt(fe_var)<<" pccm"<<endl;
    cout<<"FERND: Analytical RMS: "<<sqrt(get_variance_rslt)<<" pccm"<<endl;
    cout<<"FERND: Missing RMS: "<<sqrt(get_missing_rslt)<<" pccm"<<endl;
    //double fe_mean = toolkit::Mean(grid->fftw_fe,grid->full_size);
    //cout<<"FERND: Numrerical MEAN: "<<fe_mean<<"pccm"<<endl;
#endif
}

// hermiticity fixer
void FEgrnd::hermiticity(fftw_complex *a_c,const unsigned int &N1,const unsigned int &N2,const unsigned int &N3){
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
void FEgrnd::complex2real(const fftw_complex *input,double *output,const unsigned long int &size) {
    for(unsigned long int i=0;i!=size;++i){
        output[i] = input[i][0];
    }
}

// END
