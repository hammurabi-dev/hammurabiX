#include <iostream>
#include <vec3.h>
#include <vector>
#include <array>
#include <cmath>
#include <gsl/gsl_sf_synchrotron.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_integration.h>

#include "class_cre.h"
#include "class_pond.h"
#include "class_grid.h"
#include "cgs_units_file.h"
#include "namespace_toolkit.h"

using namespace std;

double CRE::get_emissivity(const vec3 &pos,Pond *par,Grid_cre *grid,const double &Bper,const bool &cue){
    cerr<<"ERR:"<<__FILE__
    <<" : in function "<<__func__<<endl
    <<" at line "<<__LINE__<<endl
    <<"DYNAMIC BINDING FAILURE"<<endl;
    exit(1);
}

double CRE::read_grid(const unsigned int &n, const vec3 &pos,Grid_cre *grid){
    cerr<<"ERR:"<<__FILE__
    <<" : in function "<<__func__<<endl
    <<" at line "<<__LINE__<<endl
    <<"DYNAMIC BINDING FAILURE"<<endl;
    exit(1);
}

void CRE::write_grid(Pond *par,Grid_cre *grid){
    cerr<<"ERR:"<<__FILE__
    <<" : in function "<<__func__<<endl
    <<" at line "<<__LINE__<<endl
    <<"DYNAMIC BINDING FAILURE"<<endl;
    exit(1);
}

/* analytical CRE flux */
// give values to spectral index and norm factor
// we drag this part out to meet users need of changing analytical model of CRE
void CRE_ana::flux_param(const vec3 &pos,Pond *par,const double &Bper, double &index,double &norm){
    // units
    const double alpha = par->creana[0];
    const double beta = par->creana[1];
    const double theta = par->creana[2];
    const double hr = par->creana[3]*CGS_U_kpc;
    const double hz = par->creana[4]*CGS_U_kpc;
    const double je = par->creana[5];
    // R0 as it should be
    const double R0 = sqrt(par->SunPosition.x*par->SunPosition.x+par->SunPosition.y*par->SunPosition.y);
    // sylindrical position
    const double r = sqrt(pos.x*pos.x+pos.y*pos.y);
    const double z = fabs(pos.z);
    // gamma,beta at 10GeV
    const double cre_gamma_10 = 10.*CGS_U_GeV/CGS_U_MEC2;
    const double cre_beta_10 = sqrt(1.-1./cre_gamma_10);
    // from flux to density unit convertion
    const double unitfactor = (4.*CGS_U_pi*CGS_U_MEC/cre_beta_10)/(CGS_U_GeV*10000.*CGS_U_cm*CGS_U_cm*CGS_U_sec);
    // coefficients which do not attend integration
    const double forefactor = pow(CGS_U_qe,2.5)*sqrt(fabs(Bper)*2.*CGS_U_MEC*2.*CGS_U_pi*par->sim_freq)/(4.*CGS_U_pi*CGS_U_MEC2);
    
    /* MODEL DEPENDENT PARAMETERS */
    // CRE flux normalizaton factor at earth, model dependent
    const double C_earth = je*unitfactor*pow(cre_gamma_10,alpha-beta*R0);
    const double normfactor = C_earth*exp(R0/hr);
    // for scaling of CRE density at spatial position, same as wmap3yr model
    // this is changeable by users
    const double scalfactor = exp(-r/hr)*(1./pow(cosh(z/hz),2.));
    // this is changeable by users
    index = -alpha+beta*r+theta*z;
    
    norm = forefactor*normfactor*scalfactor;
}

double CRE_ana::get_emissivity(const vec3 &pos,Pond *par,Grid_cre *grid,const double &Bper,const bool &cue){
    if(grid->read_permission){
        cerr<<"ERR:"<<__FILE__
        <<" : in function "<<__func__<<endl
        <<" at line "<<__LINE__<<endl
        <<"WRONG MODULE"<<endl;
        exit(1);
    }
    double J = 0.;
        // allocating values to index, norm according to user defined model
        // user may consider building derived class from CRE_ana
        double index, norm;
        flux_param(pos,par,Bper,index,norm);
        
        // synchrotron integration
        //double abserr; //switch on with gsl integration
        struct int_pars inner;
        inner.A = sqrt(2.*CGS_U_MEC)/sqrt(3.*CGS_U_qe*abs(Bper));
        inner.omega = 2.*CGS_U_pi*par->sim_freq;
        inner.index = index;
        
        /* TWO STATIC FUNCTIONS ARE DEFINED IN HEADER FILE
         double gF(double x, void * pars) {
         return pow(A*sqrt(omega/x), index ) * gsl_sf_synchrotron_1(x) * pow(x, -1.5);
         }
         double gG(double x, void * pars) {
         return pow(A*sqrt(omega/x), index ) * gsl_sf_synchrotron_2(x) * pow(x, -1.5);
         }
         */
	if(cue){
        // test with gsl integration
        /*
        gsl_function f;
        f.function = &gF;
        f.params = &inner;
        gsl_integration_workspace *wf = gsl_integration_workspace_alloc(100);
        gsl_integration_qag(&f,0,10,0,1.0e-2,100,1,wf,&J,&abserr);
        gsl_integration_workspace_free(wf);
        */
        // analytical solution
        double mu = -(3.+inner.index)/2.0;
        J = pow(inner.A*sqrt(inner.omega),inner.index)*pow(2,mu+1)*gsl_sf_gamma(0.5*mu+7./3.)*gsl_sf_gamma(0.5*mu+2./3.)/(mu+2.);
        
        return norm*J/(4.*CGS_U_pi);
    }
	else{
        // test with gsl integration
        /*
        gsl_function g;
        g.function = &gG;
        g.params = &inner;
        gsl_integration_workspace *wg = gsl_integration_workspace_alloc(100);
        gsl_integration_qag(&g,0,30,0,1.0e-2,100,1,wg,&J,&abserr);
        gsl_integration_workspace_free(wg);
        */
        // analytical solution
        double mu = -(3.+inner.index)/2.0;
        J = pow(inner.A*sqrt(inner.omega),inner.index)*pow(2,mu)*gsl_sf_gamma(0.5*mu+4./3.)*gsl_sf_gamma(0.5*mu+2./3.);
        
        return norm*J/(4.*CGS_U_pi);
	}
        /* the last 4pi comes from solid-angle integration/deviation,
         check eq(6.16) in Ribiki-Lightman's where Power is defined,
         we need isotropic power which means we need a 1/4pi factor!
         */
}

/* numerical CRE flux */
// use bilinear/trilinear interpolationi according to the dimension of CRE flux grid
double CRE_num::read_grid(const unsigned int &Eidx, const vec3 &pos,Grid_cre *grid){
    // if grid in spatial 2D
    if(grid->nx==0){
        // sylindrical galactic centric position
        const double r = sqrt(pos.x*pos.x+pos.y*pos.y);
        const double z = pos.z;
        // bilinear interpolation
        decltype(grid->nr) rl, zl;
        // notice that lr is radius
        double tmp = (grid->nr-1)*r/grid->lr;
        if(tmp<0 or tmp>grid->nr-1) {return 0.;}
        else rl = floor(tmp);
        const double rd = tmp - rl;
        
        tmp = (grid->nz-1)*(z/grid->lz + 0.5);
        if(tmp<0 or tmp>grid->nz-1) {return 0.;}
        else zl = floor(tmp);
        const double zd = tmp - zl;
        
        if(rl+1<grid->nr and zl+1<grid->nz){
            auto idx1 = toolkit::Index3d(grid->nE,grid->nr,grid->nz,Eidx,rl,zl);
            auto idx2 = toolkit::Index3d(grid->nE,grid->nr,grid->nz,Eidx,rl+1,zl);
            double i1 = grid->cre_flux[idx1]*(1-rd) + grid->cre_flux[idx2]*rd;
            
            idx1 = toolkit::Index3d(grid->nE,grid->nr,grid->nz,Eidx,rl,zl+1);
            idx2 = toolkit::Index3d(grid->nE,grid->nr,grid->nz,Eidx,rl+1,zl+1);
            double i2 = grid->cre_flux[idx1]*(1-rd) + grid->cre_flux[idx2]*rd;
#ifndef NDEBUG
            if(i1<0 or i2<0){
                cerr<<"ERR:"<<__FILE__
                <<" : in function "<<__func__<<endl
                <<" at line "<<__LINE__<<endl
                <<"NEGATIVE CRE FLUX"<<endl;
                exit(1);
            }
#endif
            return (i1*(1-zd)+i2*zd);
        }
        else{
            auto idx1 = toolkit::Index3d(grid->nE,grid->nr,grid->nz,Eidx,rl,zl);
            return grid->cre_flux[idx1];
        }
    }
    // if grid in spatial 3D
    else if(grid->nr==0){
        //trilinear interpolation
        decltype(grid->nx) xl, yl, zl;
        
        double tmp = (grid->nx-1)*(pos.x/grid->lx + 0.5);
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
            <<"WRONG VALUE: "<<endl;
            exit(1);
        }
#endif
        if (xl+1<grid->nx and yl+1<grid->ny and zl+1<grid->nz) {
            double i1,i2,j1,j2,w1,w2;
            
            auto idx1 = toolkit::Index4d(grid->nE,grid->nx,grid->ny,grid->nz,Eidx,xl,yl,zl);
            auto idx2 = toolkit::Index4d(grid->nE,grid->nx,grid->ny,grid->nz,Eidx,xl,yl,zl+1);
            i1=grid->cre_flux[idx1]*(1.-zd) + grid->cre_flux[idx2]*zd;
            
            idx1 = toolkit::Index4d(grid->nE,grid->nx,grid->ny,grid->nz,Eidx,xl,yl+1,zl);
            idx2 = toolkit::Index4d(grid->nE,grid->nx,grid->ny,grid->nz,Eidx,xl,yl+1,zl+1);
            i2=grid->cre_flux[idx1]*(1-zd) + grid->cre_flux[idx2]*zd;
            
            idx1 = toolkit::Index4d(grid->nE,grid->nx,grid->ny,grid->nz,Eidx,xl+1,yl,zl);
            idx2 = toolkit::Index4d(grid->nE,grid->nx,grid->ny,grid->nz,Eidx,xl+1,yl,zl+1);
            j1=grid->cre_flux[idx1]*(1-zd) + grid->cre_flux[idx2]*zd;
            
            idx1 = toolkit::Index4d(grid->nE,grid->nx,grid->ny,grid->nz,Eidx,xl+1,yl+1,zl);
            idx2 = toolkit::Index4d(grid->nE,grid->nx,grid->ny,grid->nz,Eidx,xl+1,yl+1,zl+1);
            j2=grid->cre_flux[idx1]*(1-zd) + grid->cre_flux[idx2]*zd;
            
            w1=i1*(1-yd)+i2*yd;
            w2=j1*(1-yd)+j2*yd;
            
            return (w1*(1-xd)+w2*xd);
            
        }// if not edge
        else {
            auto idx1 = toolkit::Index4d(grid->nE,grid->nx,grid->ny,grid->nz,Eidx,xl,yl,zl);
            return grid->cre_flux[idx1];
        }
    }
    else{
        cerr<<"ERR:"<<__FILE__
        <<" : in function "<<__func__<<endl
        <<" at line "<<__LINE__<<endl
        <<"WRONG CRE GRID DIMENSION"<<endl;
        exit(1);
    }
}

double CRE_num::get_emissivity(const vec3 &pos,Pond *par,Grid_cre *grid,const double &Bper,const bool &cue){
    double J = 0.;
    if(!grid->read_permission){
        cerr<<"ERR:"<<__FILE__
        <<" : in function "<<__func__<<endl
        <<" at line "<<__LINE__<<endl
        <<"NO INPUT"<<endl;
        exit(1);
    }
    
    // allocate energy grid
    double KE[grid->nE] = {0.};
    // we want F(x[E]) and G(x[E]) in spectral integration
    // also we need beta, for getting correct density
    double x[grid->nE] = {0.};
    double beta[grid->nE] = {0.};
    // consts used in loop, using cgs units
    const double shift_fact = (2.*CGS_U_MEC*CGS_U_MEC2*CGS_U_MEC2*2.*CGS_U_pi*par->sim_freq)/(3.*CGS_U_qe*Bper*CGS_U_GeV*CGS_U_GeV);
    // KE in GeV
    const double gamma_fact = CGS_U_MEC2/CGS_U_GeV;
    for(decltype(grid->nE) i=0;i!=grid->nE;++i){
        KE[i] = exp(log(grid->Ekmin) + (double)i*log(grid->Ekfact));
        x[i] = shift_fact/(KE[i]*KE[i]);
        beta[i] = fabs(sqrt(1-gamma_fact/KE[i]));
    }
    
    // do energy spectrum integration at given position
    // unit_factor for DIFFERENTIAL density flux, [GeV m^2 s sr]^-1
    // n(E,pos) = \phi(E,pos)*(4\pi/\beta*c), the relatin between flux \phi and density n
    // ref: "Cosmic rays n' particle physics", A3
    // we integrate over E instead of x
    // we drop 1/GeV since dE in GeV
    const double unitfactor = 4.*CGS_U_pi/(CGS_U_C_light*100.*CGS_U_cm*100.*CGS_U_cm*CGS_U_sec);
    
    for(decltype(grid->nE) i=0;i!=grid->nE-1;++i){
        double xv = (x[i+1]+x[i])/2.;
        // we control value of xv to avoid underflow in gsl functions
        if(xv>10){continue;}
        double dE = fabs(KE[i+1]-KE[i]);
        // we use beta here
        double de = unitfactor*(read_grid(i+1,pos,grid)/beta[i+1]+read_grid(i,pos,grid)/beta[i])/2.;
#ifndef NDEBUG
        if(de<0){
            cerr<<"ERR:"<<__FILE__
            <<" : in function "<<__func__<<endl
            <<" at line "<<__LINE__<<endl
            <<"NEGATIVE CRE DENSITY"<<endl;
            exit(1);
        }
#endif
	if(cue){
        J += gsl_sf_synchrotron_1(xv)*de*dE;
	}
	else{
        J += gsl_sf_synchrotron_2(xv)*de*dE;
	}
    }
    // fore factor
    const double fore_factor = sqrt(3.)*pow(CGS_U_qe,3.)*abs(Bper)/(2.*CGS_U_pi*CGS_U_MEC2);
    return fore_factor*J/(4.*CGS_U_pi);
    // 4\pi here explained in type1 function
}

// END