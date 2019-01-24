#include <cassert>
#include <omp.h>
#include <cmath>

#include <vec3.h>
#include <fftw3.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_integration.h>

#include <param.h>
#include <grid.h>
#include <brnd.h>
#include <breg.h>
#include <cgs_units_file.h>
#include <namespace_toolkit.h>

// global anisotropic turbulent field
vec3_t<double> Brnd_es::anisotropy_direction (const vec3_t<double> &pos,
                                              const Param *par,
                                              const Breg *breg,
                                              const Grid_breg *gbreg) const{
    return toolkit::versor (breg->get_breg (pos,par,gbreg));
}

// global anisotropic turbulent field
double Brnd_es::anisotropy_ratio (const vec3_t<double> &,
                                  const Param *par,
                                  const Breg *,
                                  const Grid_breg *) const{
    // the simplest case, const.
    return par->brnd_es.rho;
}

// since we are using rms normalization
// p0 is hidden and not affecting anything
double Brnd_es::spec(const double &k,
                     const Param *par) const{
    //units fixing, wave vector in 1/kpc units
    const double p0 {par->brnd_es.rms*CGS_U_muGauss};
    const double kr {k/par->brnd_es.k0};
    const double a0 {par->brnd_es.a0};
    const double unit = 1./(4*CGS_U_pi*k*k);
    // power law
    double P {0.};
    if (kr>1){
        P = p0/std::pow(kr,a0);
    }
    return P*unit;
}

// galactic scaling of random field energy density
// set to 1 at observer's place
double Brnd_es::rescal (const vec3_t<double> &pos,
                        const Param *par) const{
    const double r_cyl {std::sqrt(pos.x*pos.x+pos.y*pos.y) - std::fabs(par->SunPosition.x)};
    const double z {std::fabs(pos.z) - std::fabs(par->SunPosition.z)};
    const double r0 {par->brnd_es.r0};
    const double z0 {par->brnd_es.z0};
    return std::exp(-r_cyl/r0)*std::exp(-z/z0);
}

// Gram-Schimdt, rewritten using Healpix vec3 library
// tiny error caused by machine is inevitable
vec3_t<double> Brnd_es::gramschmidt (const vec3_t<double> &k,
                                     const vec3_t<double> &b) const{
    if (k.SquaredLength()==0 or b.SquaredLength()==0){
        return vec3_t<double> {0,0,0};
    }
    const double inv_k_mod = 1./k.SquaredLength();
    vec3_t<double> b_free {b.x - k.x*dotprod(k,b)*inv_k_mod,
        b.y - k.y*dotprod(k,b)*inv_k_mod,
        b.z - k.z*dotprod(k,b)*inv_k_mod};
    
    b_free = toolkit::versor(b_free)*b.Length();
    return b_free;
}

void Brnd_es::write_grid (const Param *par,
                          const Breg *breg,
                          const Grid_breg *gbreg,
                          Grid_brnd *grid) const{
    // PHASE I
    // GENERATE GAUSSIAN RANDOM FROM SPECTRUM
    // initialize random seed
#ifdef _OPENMP
    gsl_rng **threadvec = new gsl_rng *[omp_get_max_threads()];
    for (int b=0;b<omp_get_max_threads();++b){
        threadvec[b] = gsl_rng_alloc(gsl_rng_taus);
        gsl_rng_set (threadvec[b],b+toolkit::random_seed(par->brnd_seed));
    }
#else
    gsl_rng *r {gsl_rng_alloc(gsl_rng_taus)};
    gsl_rng_set (r, toolkit::random_seed(par->brnd_seed));
#endif
    // start Fourier space filling, physical k in 1/kpc dimension
    const double lx {par->grid_brnd.x_max-par->grid_brnd.x_min};
    const double ly {par->grid_brnd.y_max-par->grid_brnd.y_min};
    const double lz {par->grid_brnd.z_max-par->grid_brnd.z_min};
    // physical dk^3
    const double dk3 {CGS_U_kpc*CGS_U_kpc*CGS_U_kpc/(lx*ly*lz)};
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (decltype(par->grid_brnd.nx) i=0;i<par->grid_brnd.nx;++i){
#ifdef _OPENMP
        auto seed_id = threadvec[omp_get_thread_num()];
#else
        auto seed_id = r;
#endif
        double kx {CGS_U_kpc*i/lx};
        if (i>=(par->grid_brnd.nx+1)/2) kx -= CGS_U_kpc*par->grid_brnd.nx/lx;
        // it's better to calculate indeces manually
        // just for reference, how indeces are calculated
        // const size_t idx {toolkit::index3d(par->grid_brnd.nx,par->grid_brnd.ny,par->grid_brnd.nz,i,j,l)};
        const std::size_t idx_lv1 {i*par->grid_brnd.ny*par->grid_brnd.nz};
        for (decltype(par->grid_brnd.ny) j=0;j<par->grid_brnd.ny;++j){
            double ky {CGS_U_kpc*j/ly};
            if (j>=(par->grid_brnd.ny+1)/2) ky -= CGS_U_kpc*par->grid_brnd.ny/ly;
            const std::size_t idx_lv2 {idx_lv1+j*par->grid_brnd.nz};
            for (decltype(par->grid_brnd.nz) l=0;l<par->grid_brnd.nz;++l){
                // 0th term is fixed to zero in allocation
                if (i==0 and j==0 and l==0) continue;
                double kz {CGS_U_kpc*l/lz};
                if (l>=(par->grid_brnd.nz+1)/2) kz -= CGS_U_kpc*par->grid_brnd.nz/lz;
                const double ks {std::sqrt(kx*kx + ky*ky + kz*kz)};
                const std::size_t idx {idx_lv2+l};
                // turbulent power is shared in following pattern
                // P ~ (bx^2 + by^2 + bz^2)
                // c0^2 ~ c1^2 ~ (bx^2 + by^2) ~ P*2/3
                // as renormalization comes in PHASE II,
                // 2/3, P0 in spec, dk3 are numerically redundant
                const double sigma {std::sqrt(0.66666667*spec(ks,par)*dk3)};
                grid->c0[idx][0] = sigma*gsl_ran_ugaussian(seed_id);
                grid->c0[idx][1] = sigma*gsl_ran_ugaussian(seed_id);
                grid->c1[idx][0] = sigma*gsl_ran_ugaussian(seed_id);
                grid->c1[idx][1] = sigma*gsl_ran_ugaussian(seed_id);
            }// l
        }// j
    }// i
    // ks=0 should be automatically addressed in P(k)
    // execute DFT backward plan
    fftw_execute_dft (grid->plan_c0_bw,grid->c0,grid->c0);
    fftw_execute_dft (grid->plan_c1_bw,grid->c1,grid->c1);
    // free random memory
#ifdef _OPENMP
    for (int b=0;b<omp_get_max_threads();++b)
        gsl_rng_free(threadvec[b]);
    delete [] threadvec;
#else
    // free random memory
    gsl_rng_free(r);
#endif
    // PHASE II
    // RESCALING FIELD PROFILE IN REAL SPACE
    // 1./std::sqrt(3*bi_var)
    // after 1st Fourier transformation, c0_R = bx, c0_I = by, c1_I = bz
    const double b_var_invsq {1./std::sqrt(3.*toolkit::variance(grid->c0[0],par->grid_brnd.full_size))};
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (decltype(par->grid_brnd.nx) i=0;i<par->grid_brnd.nx;++i){
        vec3_t<double> pos {i*lx/(par->grid_brnd.nx-1) + par->grid_brnd.x_min,0,0};
        const std::size_t idx_lv1 {i*par->grid_brnd.ny*par->grid_brnd.nz};
        for (decltype(par->grid_brnd.ny) j=0;j<par->grid_brnd.ny;++j){
            pos.y = j*ly/(par->grid_brnd.ny-1) + par->grid_brnd.y_min;
            const std::size_t idx_lv2 {idx_lv1+j*par->grid_brnd.nz};
            for (decltype(par->grid_brnd.nz) l=0;l<par->grid_brnd.nz;++l){
                // get physical position
                pos.z = l*lz/(par->grid_brnd.nz-1) + par->grid_brnd.z_min;
                // get rescaling factor
                double ratio {std::sqrt(rescal(pos,par))*par->brnd_es.rms*b_var_invsq};
                const std::size_t idx {idx_lv2+l};
                // assemble b_Re
                // after 1st Fourier transformation, c0_R = bx, c0_I = by, c1_I = bz
                vec3_t<double> b_re {grid->c0[idx][0]*ratio,
                    grid->c0[idx][1]*ratio,
                    grid->c1[idx][1]*ratio};
                // impose anisotropy
                vec3_t<double> H_versor = anisotropy_direction (pos,par,breg,gbreg);
                double rho {anisotropy_ratio (pos,par,breg,gbreg)};
                assert (rho>=0. and rho<=1.);
                if (H_versor.SquaredLength()<1e-10)// zero regular field, no prefered anisotropy
                    continue;
                vec3_t<double> b_re_par {H_versor*dotprod(H_versor,b_re)};
                vec3_t<double> b_re_perp {b_re - b_re_par};
                b_re = toolkit::versor(b_re_par*rho + b_re_perp*(1-rho))*b_re.Length();
                // push b_re back to c0 and c1
                grid->c0[idx][0] = b_re.x;
                grid->c0[idx][1] = b_re.y;
                grid->c1[idx][0] = b_re.y;
                grid->c1[idx][1] = b_re.z;
            } //l
        } //j
    } //i
    // execute DFT forward plan
    fftw_execute_dft (grid->plan_c0_fw,grid->c0,grid->c0);
    fftw_execute_dft (grid->plan_c1_fw,grid->c1,grid->c1);
    // PHASE III
    // RE-ORTHOGONALIZING IN FOURIER SPACE
    // Gram-Schmidt process
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (decltype(par->grid_brnd.nx) i=0;i<par->grid_brnd.nx;++i){
        decltype(par->grid_brnd.nx) i_sym {par->grid_brnd.nx-i};// apply Hermitian symmetry
        if (i==0) i_sym = i;
        vec3_t<double> tmp_k {CGS_U_kpc*i/lx,0,0};
        // it's better to calculate indeces manually
        // just for reference, how indeces are calculated
        // const std::size_t idx {toolkit::index3d(par->grid_brnd.nx,par->grid_brnd.ny,par->grid_brnd.nz,i,j,l)};
        // const std::size_t idx_sym {toolkit::index3d(par->grid_brnd.nx,par->grid_brnd.ny,par->grid_brnd.nz,i_sym,j_sym,l_sym)};
        const std::size_t idx_lv1 {i*par->grid_brnd.ny*par->grid_brnd.nz};
        const std::size_t idx_sym_lv1 {i_sym*par->grid_brnd.ny*par->grid_brnd.nz};
        if (i>=(par->grid_brnd.nx+1)/2) tmp_k.x -= CGS_U_kpc*par->grid_brnd.nx/lx;
        for (decltype(par->grid_brnd.ny) j=0;j<par->grid_brnd.ny;++j){
            decltype(par->grid_brnd.ny) j_sym {par->grid_brnd.ny-j};// apply Hermitian symmetry
            if (j==0) j_sym = j;
            tmp_k.y = CGS_U_kpc*j/ly;
            if (j>=(par->grid_brnd.ny+1)/2) tmp_k.y -= CGS_U_kpc*par->grid_brnd.ny/ly;
            const std::size_t idx_lv2 {idx_lv1+j*par->grid_brnd.nz};
            const std::size_t idx_sym_lv2 {idx_sym_lv1+j_sym*par->grid_brnd.nz};
            for (decltype(par->grid_brnd.nz) l=0;l<par->grid_brnd.nz;++l){
                decltype(par->grid_brnd.nz) l_sym {par->grid_brnd.nz-l};// apply Hermitian symmetry
                if (l==0) l_sym = l;
                tmp_k.z = CGS_U_kpc*l/lz;
                if (l>=(par->grid_brnd.nz+1)/2) tmp_k.z -= CGS_U_kpc*par->grid_brnd.nz/lz;
                const std::size_t idx {idx_lv2+l}; //k
                const std::size_t idx_sym {idx_sym_lv2+l_sym}; //-k
                // reconstruct bx,by,bz from c0,c1,c*0,c*1
                // c0(k) = bx(k) + i by(k)
                // c*0(-k) = bx(k) - i by(k)
                // c1(k) = by(k) + i bz(k)
                // c1*1(-k) = by(k) - i bz(k)
                const vec3_t<double> tmp_b_re {0.5*(grid->c0[idx][0]+grid->c0[idx_sym][0]),
                    0.5*(grid->c1[idx][0]+grid->c1[idx_sym][0]),
                    0.5*(grid->c1[idx_sym][1]+grid->c1[idx][1])};
                
                const vec3_t<double> tmp_b_im {0.5*(grid->c0[idx][1]-grid->c0[idx_sym][1]),
                    0.5*(grid->c1[idx][1]-grid->c1[idx_sym][1]),
                    0.5*(grid->c1[idx_sym][0]-grid->c1[idx][0])};
                
                const vec3_t<double> free_b_re {gramschmidt(tmp_k,tmp_b_re)};
                const vec3_t<double> free_b_im {gramschmidt(tmp_k,tmp_b_im)};
                // reassemble c0,c1 from bx,by,bz
                // c0(k) = bx(k) + i by(k)
                // c1(k) = by(k) + i bz(k)
                grid->c0[idx][0] = free_b_re.x-free_b_im.y;
                grid->c0[idx][1] = free_b_im.x+free_b_re.y;
                grid->c1[idx][0] = free_b_re.y-free_b_im.z;
                grid->c1[idx][1] = free_b_im.y+free_b_re.z;
            }// l
        }// j
    }// i
    // execute DFT backward plan
    fftw_execute_dft (grid->plan_c0_bw,grid->c0,grid->c0);
    fftw_execute_dft (grid->plan_c1_bw,grid->c1,grid->c1);
    // according to FFTW convention
    // transform forward followed by backword scale up array by nx*ny*nz
    double inv_grid_size = 1.0/par->grid_brnd.full_size;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t i=0;i<par->grid_brnd.full_size;++i){
        grid->bx[i] = grid->c0[i][0]*inv_grid_size;
        grid->by[i] = grid->c0[i][1]*inv_grid_size;
        grid->bz[i] = grid->c1[i][1]*inv_grid_size;
    }
}
