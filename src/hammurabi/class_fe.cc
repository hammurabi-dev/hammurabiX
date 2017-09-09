#include <iostream>
#include <cmath>
#include <vec3.h>
#include <omp.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <string>
#include <fstream>

#include "class_fe.h"
#include "class_pond.h"
#include "class_grid.h"
#include "namespace_toolkit.h"
#include "cgs_units_file.h"

using namespace std;

double FE::get_density(const vec3 &pos, Pond *par, Grid_fe *grid){
    if(grid->permission){
        return read_grid(pos,grid,par);
    }
    else {
        return density(pos,par);
    }
}

// not recommended to use without enough computing source
// recommend to use this once (replace density in write_grid) if no
// free parameters in FE
double FE::density_blur(const vec3 &pos, Pond *par, Grid_fe *grid){
    double ne_blur = 0.;
    // element value in each thread
    double element = 0.;
    // sampling point number
    unsigned int step = 1000;
    // gaussian blur scale
    double blur_scale_x = grid->lx/(grid->nx*CGS_U_kpc);
    double blur_scale_y = grid->ly/(grid->ny*CGS_U_kpc);
    double blur_scale_z = grid->lz/(grid->nz*CGS_U_kpc);
    // sample position
    vec3 pos_s;
    gsl_rng *r = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(r, toolkit::random_seed());
#pragma omp parallel sections private(element) reduction(+:ne_blur)
    {
#pragma omp section
        {
            for(decltype(step)i=0;i!=step;++i){
                pos_s.x = pos.x + gsl_ran_gaussian(r,(blur_scale_x/2.355))*CGS_U_kpc;
                pos_s.y = pos.y + gsl_ran_gaussian(r,(blur_scale_y/2.355))*CGS_U_kpc;
                pos_s.z = pos.z + gsl_ran_gaussian(r,(blur_scale_z/2.355))*CGS_U_kpc;
                element = density(pos_s,par);
                ne_blur += element;
            }
        }// omp section
    }//parallel
    gsl_rng_free(r);
    return ne_blur/step;
}

double FE::density(const vec3 &pos, Pond *par){
    cerr<<"ERR:"<<__FILE__
    <<" : in function "<<__func__<<endl
    <<" at line "<<__LINE__<<endl
    <<"DYNAMIC BINDING FAILURE"<<endl;
    exit(1);
    return 0.;
}

double FE::read_grid(const vec3 &pos, Grid_fe *grid, Pond *par){
    
    decltype(grid->nx) xl, yl, zl;
    
    double tmp = (grid->nx-1)*(pos.x/grid->lx + 0.5);
    // if use solar-centric grid
    if(grid->ec_frame) tmp -= (grid->nx-1)*(par->SunPosition.x/grid->lx);
    if (tmp<1 or tmp>grid->nx-1) { return 0.;}
    else xl = floor(tmp);
    const double xd = tmp - xl;
    
    tmp = (grid->ny-1)*(pos.y/grid->ly + 0.5);
    if (tmp<1 or tmp>grid->ny-1) { return 0.;}
    else yl = floor(tmp);
    const double yd = tmp - yl;
    
    tmp = (grid->nz-1)*(pos.z/grid->lz + 0.5);
    if (tmp<1 or tmp>grid->nz-1) { return 0.;}
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
        
        auto idx1 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl,yl,zl);
        auto idx2 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl,yl,zl+1);
        i1=grid->fe[idx1]*(1.-zd) + grid->fe[idx2]*zd;
        
        idx1 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl,yl+1,zl);
        idx2 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl,yl+1,zl+1);
        i2=grid->fe[idx1]*(1-zd) + grid->fe[idx2]*zd;
        
        idx1 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl+1,yl,zl);
        idx2 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl+1,yl,zl+1);
        j1=grid->fe[idx1]*(1-zd) + grid->fe[idx2]*zd;
        
        idx1 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl+1,yl+1,zl);
        idx2 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl+1,yl+1,zl+1);
        j2=grid->fe[idx1]*(1-zd) + grid->fe[idx2]*zd;
        
        w1=i1*(1-yd)+i2*yd;
        w2=j1*(1-yd)+j2*yd;
        
        return (w1*(1-xd)+w2*xd);
        
    }
    else {
        auto idx1 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl,yl,zl);
        return grid->fe[idx1];
    }
    
}

void FE::write_grid(Pond *par, Grid_fe *grid){
    if(!grid->permission){
        cerr<<"ERR:"<<__FILE__
        <<" : in function "<<__func__<<endl
        <<" at line "<<__LINE__<<endl
        <<"NO PERMISSION"<<endl;
        exit(1);
    }
    cout<<"...FE: WRITING OUTPUT..."<<endl;
    vec3 gc_pos;
    for(decltype(grid->nx) i=0;i!=grid->nx;++i){
        for(decltype(grid->ny) j=0;j!=grid->ny;++j){
            for(decltype(grid->nz) k=0;k!=grid->nz;++k){
                auto idx = toolkit::Index3d(grid->nx,grid->ny,grid->nz,i,j,k);
                
                gc_pos.x = grid->lx*(double(i)/(grid->nx-1) - 0.5);
                gc_pos.y = grid->ly*(double(j)/(grid->ny-1) - 0.5);
                gc_pos.z = grid->lz*(double(k)/(grid->nz-1) - 0.5);
                
                if(grid->ec_frame){
                    gc_pos.x += par->SunPosition.x;
                }
                // two solutions
                //grid->fe[idx] = density_blur(gc_pos, par, grid);
                grid->fe[idx] = density(gc_pos,par);
            }
        }
    }
}


/* YMW16 */

double YMW16::density(const vec3 &pos, Pond *par){
    // YMW16 using a different Cartesian frame from our default one
    vec3 gc_pos;
    gc_pos.x = pos.y;
    gc_pos.y = -pos.x;
    gc_pos.z = pos.z;
    // sylindrical r
    double r_cyl = sqrt(gc_pos.x*gc_pos.x+gc_pos.y*gc_pos.y);
    // warp
    if(r_cyl>=par->R_warp){
        double theta_warp = atan2(gc_pos.y,gc_pos.x);
        gc_pos.z -= par->t0_Gamma_w*(r_cyl-par->R_warp)*cos(theta_warp);
    }
    
    if (gc_pos.Length()>25*CGS_U_kpc){
        return 0.;
    }
    else{
        double ne = 0.;
        double NE[8] = {0.};
        double WLB = 0., WGN = 0., WLI = 0.;
        // longitude, in deg
        const double ec_l = atan2(gc_pos.x,par->R0-gc_pos.y)*180./CGS_U_pi;
        
        //call structure functions
        //since in YMW16, Fermi Bubble is not actually contributing, we ignore FB for thick disk
        NE[1] = thick(gc_pos.z,r_cyl,par);
        NE[2] = thin(gc_pos.z,r_cyl,par);
        NE[3] = spiral(gc_pos.x,gc_pos.y,gc_pos.z,r_cyl,par);
        NE[4] = galcen(gc_pos.x,gc_pos.y,gc_pos.z,par);
        NE[5] = gum(gc_pos.x,gc_pos.y,gc_pos.z,par);
        // localbubble boundary
        const double RLB = 110.*CGS_U_pc;
        NE[6] = localbubble(gc_pos.x,gc_pos.y,gc_pos.z,ec_l,RLB,par);
        NE[7] = nps(gc_pos.x,gc_pos.y,gc_pos.z,par);
        
        //adding up rules
        NE[0] = NE[1]+max(NE[2],NE[3]);
        // distance to local bubble
        const double rlb = sqrt( pow(((gc_pos.y-8.34*CGS_U_kpc)*0.94-0.34*gc_pos.z),2.) + pow(gc_pos.x,2.) );
        if(rlb<RLB){ //inside local bubble
            NE[0]=par->t6_J_LB*NE[1]+max(NE[2],NE[3]);
            if(NE[6]>NE[0]) {WLB=1;}
        }
        else{ //outside local bubble
            if(NE[6]>NE[0] and NE[6]>NE[5]) {WLB=1;}
        }
        if(NE[7]>NE[0]) {WLI=1;}
        if(NE[5]>NE[0]) {WGN = 1;}
        // final density
        ne = (1-WLB)*( (1-WGN)*( (1-WLI)*(NE[0]+NE[4]) + WLI*NE[7] )  + WGN*NE[5] ) + (WLB)*(NE[6]);
        
        return ne;
    }
    
}
// thick disk
double YMW16::thick(const double &zz, const double &rr, Pond *par){
    double gd = 1.;
    if(rr>par->t1_Bd){
        gd = pow(1/cosh((rr-par->t1_Bd)/par->t1_Ad),2.);
    }
    return par->t1_n1*gd*pow(1/cosh(zz/par->t1_H1),2.);
}
// thin disk
double YMW16::thin(const double &zz,const double &rr, Pond *par){
    double gd = 1.;
    if(rr>par->t1_Bd){
        gd = pow(1/cosh((rr-par->t1_Bd)/par->t1_Ad),2.);
    }
    // z scaling, K_2*H in ref
    double H = par->t2_K2*(32*CGS_U_pc+1.6e-3*rr+(4.e-7/CGS_U_pc)*pow(rr,2.));
    
    return par->t2_n2*gd*pow(1/cosh((rr-par->t2_B2)/par->t2_A2),2)*pow(1/cosh(zz/H),2.);
}
// spiral arms
double YMW16::spiral(const double &xx,const double &yy,const double &zz,const double &rr,Pond *par){
    double detrr=1e10*CGS_U_pc;
    double armr, armr1, armr2, smin;
    
    double ga=0;
    //for calculating Carina-Sagittarius correction
    double RAD = 180./CGS_U_pi;
    double gd = 1.;
    if(rr>par->t1_Bd){
        gd = pow(1/cosh((rr-par->t1_Bd)/par->t1_Ad),2.);
    }
    // z scaling, K_a*H in ref
    const double H = par->t3_Ka*(32*CGS_U_pc+1.6e-3*rr+(4.e-7/CGS_U_pc)*pow(rr,2.));
    // R_ai
    const double rmin[5] = {3.35*CGS_U_kpc,3.707*CGS_U_kpc,3.56*CGS_U_kpc,3.670*CGS_U_kpc,8.21*CGS_U_kpc};
    // phi_ai
    const double thmin[5] = {44.4*CGS_U_pi/180,120*CGS_U_pi/180,218.6*CGS_U_pi/180,330.3*CGS_U_pi/180,55.1*CGS_U_pi/180};
    // tan(psi_ai)
    const double tpitch[5] = {0.20200,0.17300,0.18300,0.18600,0.048300};
    // cos(psi_ai)
    const double cspitch[5] = {0.98000,0.98500,0.98360,0.98300,0.99880};
    double theta=atan2(yy,xx);
    if(theta<0) theta += 2*CGS_U_pi;
    
    // 普通角度的计算 (general angle calculation)
    double ne3s=0.;
    // looping through arms
    for(int i=0;i!=4;++i){
        ga=0;
        //Norma-Outer
        if(i==0){
            if(theta>=0 and theta<thmin[i]){
                armr=rmin[i]*exp( (theta+2*CGS_U_pi-thmin[i])*tpitch[i] );//R_a
                detrr=fabs(rr-armr);//radial distance
            }
            if(theta>=thmin[i] and theta<2*CGS_U_pi){
                armr1=rmin[i]*exp( (theta-thmin[i])*tpitch[i] );
                armr2=rmin[i]*exp( (theta+2*CGS_U_pi-thmin[i])*tpitch[i] );
                detrr=min(fabs(rr-armr1), fabs(rr-armr2));
            }
        }
        //Perseus
        if(i==1){
            if(theta>=0 and theta<thmin[i]){
                armr=rmin[i]*exp( (theta+2*CGS_U_pi-thmin[i])*tpitch[i] );
                detrr=fabs(rr-armr);
            }
            if(theta>=thmin[i] and theta<2*CGS_U_pi){
                armr1=rmin[i]*exp( (theta-thmin[i])*tpitch[i] );
                armr2=rmin[i]*exp( (theta+2*CGS_U_pi-thmin[i])*tpitch[i] );
                detrr=min(fabs(rr-armr1), fabs(rr-armr2));
            }
        }
        //Carina-Sagittarius
        if(i==2){
            if(theta>=0 and theta<thmin[i]){
                armr1=rmin[i]*exp( (theta+2*CGS_U_pi-thmin[i])*tpitch[i] );
                armr2=rmin[i]*exp( (theta+4*CGS_U_pi-thmin[i])*tpitch[i] );
                detrr=min(fabs(rr-armr1), fabs(rr-armr2));
            }
            if(theta>=thmin[i] and theta<2*CGS_U_pi){
                armr1=rmin[i]*exp( (theta-thmin[i])*tpitch[i] );
                armr2=rmin[i]*exp( (theta+2*CGS_U_pi-thmin[i])*tpitch[i] );
                detrr=min(fabs(rr-armr1), fabs(rr-armr2));
            }
        }
        //Crux_Scutum
        if(i==3){
            if(theta>=0 and theta<thmin[i]){
                armr1=rmin[i]*exp( (theta+2*CGS_U_pi-thmin[i])*tpitch[i] );
                armr2=rmin[i]*exp( (theta+4*CGS_U_pi-thmin[i])*tpitch[i] );
                detrr=min(fabs(rr-armr1), fabs(rr-armr2));
            }
            if(theta>=thmin[i] and theta<2*CGS_U_pi){
                armr1=rmin[i]*exp( (theta-thmin[i])*tpitch[i] );
                armr2=rmin[i]*exp( (theta+2*CGS_U_pi-thmin[i])*tpitch[i] );
                detrr=min(fabs(rr-armr1), fabs(rr-armr2));
            }
        }
        //Local
        if(i==4){
            if(theta>=0 and theta<thmin[i]){
                detrr=1e10*CGS_U_pc;//a very big value
            }
            if(theta>=thmin[i] and theta<2){
                armr=rmin[i]*exp((theta-thmin[i])*tpitch[i]);
                detrr=fabs(rr-armr);
            }
            if(theta>=2 and theta<2*CGS_U_pi){
                detrr=1e10*CGS_U_pc;//a very big value
            }
        }
        smin=detrr*cspitch[i];//s_ai
        // correction for Carina-Sagittarius
        if(i==2){
            ga=(1-(par->t3_nsg)*(exp(-pow((theta*RAD-par->t3_thetasg)/par->t3_wsg,2.))))*(1+par->t3_ncn*exp(-pow((theta*RAD-par->t3_thetacn)/par->t3_wcn,2.)))*pow(1/cosh(smin/par->t3_warm[i]),2.);
            
            if(rr>6*CGS_U_kpc and theta*RAD>par->t3_thetacn) {
                ga=(1-(par->t3_nsg)*(exp(-pow((theta*RAD-par->t3_thetasg)/par->t3_wsg,2.))))*(1+par->t3_ncn)*pow(1/cosh(smin/par->t3_warm[i]),2.);
            }
        }
        // for other three arms
        else {
            ga=pow(1/cosh(smin/par->t3_warm[i]),2.);
        }
        
        ne3s += par->t3_narm[i]*ga*pow(1/cosh((rr-par->t3_B2s)/par->t3_Aa),2.)*gd*pow(1/cosh(zz/H),2.);
    }// end of looping through arms
    return ne3s;
}
// galactic center
double YMW16::galcen(const double &xx,const double &yy,const double &zz,Pond *par){
    //pos of center
    const double Xgc = 50.*CGS_U_pc;
    const double Ygc = 0.;
    const double Zgc = -7.*CGS_U_pc;
    const double Rgc = sqrt((xx-Xgc)*(xx-Xgc)+(yy-Ygc)*(yy-Ygc));
    const double Ar = exp(-pow(Rgc/par->t4_Agc,2.));
    const double Az = pow(1/cosh((zz-Zgc)/par->t4_Hgc),2.);
    
    return par->t4_ngc*Ar*Az;
}
// gum nebula
double YMW16::gum(const double &xx,const double &yy,const double &zz,Pond *par){
    double Dmin=1e5*CGS_U_pc;
    //center of Gum Nebula
    const double lc = 264.*CGS_U_pi/180.;
    const double bc = -4.*CGS_U_pi/180.;
    const double dc = 450.*CGS_U_pc;
    
    const double xc = dc*cos(bc)*sin(lc);
    const double yc = par->R0-dc*cos(bc)*cos(lc);
    const double zc = dc*sin(bc);
    
    const double theta = fabs(atan((zz-zc)/sqrt((xx-xc)*(xx-xc)+(yy-yc)*(yy-yc))));
    double zp = ((par->t5_Agn)*(par->t5_Agn)*(par->t5_Kgn))/sqrt(((par->t5_Agn)*(par->t5_Agn))+((par->t5_Agn)*(par->t5_Agn)*(par->t5_Kgn)*(par->t5_Kgn))/(tan(theta)*tan(theta)));
    const double xyp = zp/tan(theta);
    const double alpha = -atan((-(par->t5_Agn)*(par->t5_Kgn)*xyp)/((par->t5_Agn)*sqrt((par->t5_Agn)*(par->t5_Agn)-xyp*xyp)));
    const double RR = sqrt((xx-xc)*(xx-xc)+(yy-yc)*(yy-yc)+(zz-zc)*(zz-zc));
    const double rp = sqrt((zp)*(zp)+(xyp)*(xyp));
    Dmin = fabs((RR-rp)*sin(theta+alpha));
    
    return par->t5_ngn*exp(-pow(Dmin/par->t5_Wgn,2));
}
// local bubble
double YMW16::localbubble(const double &xx,const double &yy,const double &zz,const double &ll,const double &Rlb,Pond *par){
    // r_LB in ref
    const double rLB = sqrt(pow(((yy-8.34*CGS_U_kpc)*0.94-0.34*zz),2)+pow(xx,2));
    // l-l_LB1 in ref
    const double dl1 = min(fabs(ll+360.-par->t6_thetalb1),fabs(par->t6_thetalb1-(ll)));
    const double nelb1 = par->t6_nlb1*pow(1/cosh(dl1/par->t6_detlb1),2)*pow(1/cosh((rLB-Rlb)/par->t6_wlb1),2)*pow(1/cosh(zz/par->t6_hlb1),2);
    // l-l_LB2 in ref
    const double dl2 = min(fabs(ll+360-par->t6_thetalb2),fabs(par->t6_thetalb2-(ll)));
    const double nelb2 = par->t6_nlb2*pow(1/cosh(dl2/par->t6_detlb2),2)*pow(1/cosh((rLB-Rlb)/par->t6_wlb2),2)*pow(1/cosh(zz/par->t6_hlb2),2);
    
    return nelb1+nelb2;
}
// north spur
double YMW16::nps(const double &xx,const double &yy,const double &zz,Pond *par){
    const double theta_LI = (par->t7_thetaLI)*CGS_U_pi/180.;
    const double x_c = -10.156*CGS_U_pc;
    const double y_c = 8106.207*CGS_U_pc;
    const double z_c = 10.467*CGS_U_pc;
    // r_LI in ref
    const double rLI = sqrt((xx-x_c)*(xx-x_c)+(yy-y_c)*(yy-y_c)+(zz-z_c)*(zz-z_c));
    const double theta = acos(((xx-x_c)*(cos(theta_LI))+(zz-z_c)*(sin(theta_LI)))/rLI)*180./CGS_U_pi;
    
    return (par->t7_nLI)*exp(-pow((rLI-par->t7_RLI)/par->t7_WLI,2.))*exp(-pow(theta/par->t7_detthetaLI,2.));
}

/* NE2001 */

// END
