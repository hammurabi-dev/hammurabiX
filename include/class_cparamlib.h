/*
 cparamlib is developed by:
 Niklas Karlsson
 Johann Cohen-Tanugi
 published in:
 https://github.com/niklask/cparamlib
 The Astrophysical Journal, 647:692–708, 2006 August 10
 
 version updated to Jun.2017,
 Include impelmentation of Kamae model in Dragon/Kamae.cc
 
 we merge functions into a single file under c++11 std
 */
#ifndef CPARAMLIB_H
#define CPARAMLIB_H

#include <iostream>
#include <cmath>

class Cparamlib{
    public:
    /* Definition of particle ids */
    enum PARTICLE_IDS{
        ID_GAMMA,    /**< Id for gamma ray */
        ID_ELECTRON, /**< Id for electron */
        ID_POSITRON, /**< Id for positron */
        ID_NUE,      /**< Id for electron neutrino */
        ID_NUMU,     /**< Id for muon neutrino */
        ID_ANTINUE,  /**< Id for electron antineutrino */
        ID_ANTINUMU  /**< Id for muon antineutrino */
    };
    
    /**
     * c struct containing the four double arrays used for calculating inclusive
     * cross sections, dsigma/dlog(E). Each double array corresponds to a parameter
     * array (a, b, c and d) in Kamae et al. (2006).
     */
    struct PARAMSET{
        /**
         * Double array (sizeof = 9) which stores parameters
         * a<sub>0</sub>,...,a<sub>8</sub> describing the inclusive cross section
         * for secondary particles produced by the nondiffractive interaction.
         */
        double a[9];
        /**
         * Double array (sizeof = 8) which stores parameters
         * b<sub>0</sub>,...,b<sub>7</sub> describing the inclusive cross section
         * for secondary particles produced by diffraction dissociation.
         */
        double b[8];
        /**
         * Double array (sizeof = 5) which stores parameters
         * c<sub>0</sub>,...,c<sub>4</sub> describing the inclusive cross section
         * for secondary particles produced by the Delta(1232) resonance.
         */
        double c[5];
        /**
         * Double array (sizeof = 5) which stores parameters
         * d<sub>0</sub>,...,d<sub>4</sub> describing the inclusive cross section
         * for secondary particles produced by the res(1600) resonance.
         */
        double d[5];
    };
    
    /**
     * c struct containing doubles and double arrays used for calculating pT
     * distributions, i.e. the differential cross section
     * d<sup>2</sup>sigma/dlog(E)dp<sub>t</sub>,
     * as defined in Karlsson and Kamae (2008).
     *
     * Doubles correspond to parameters describing the differential cross
     * section as functions of the secondary particle energy (E) and and
     * transverse momentum (pT). Double arrays describe above parameters as
     * functions of proton kinetic energy (Tp).
     */
    struct PARAMSET_PT{
        /**
         * Double array (sizeof = 6) which stores parameters a<sub>0</sub>,
         * a<sub>1</sub> and a<sub>10</sub>,...,a<sub>13</sub> (in that order)
         * describing the p<sub>t</sub> distribution for gamma rays produced by
         * nonresonance interaction.
         */
        double a[6];
        /**
         * Array of three double arrays (sizeof = 3x4) which stores parameters
         * b<sub>0</sub>, b<sub>1</sub>, b<sub>2</sub>,
         * b<sub>10</sub>,...,b<sub>13</sub> and b<sub>20</sub>,...,b<sub>23</sub>
         * describing the p<sub>t</sub> distribution for gamma rays produced by
         * the Delta(1232) resonance.
         *
         * Layout:
         * - b[0][..] = b<sub>0</sub>, b<sub>1</sub>, b<sub>2</sub>
         * - b[1][..] = b<sub>10</sub>,...,b<sub>13</sub>
         * - b[2][..] = b<sub>20</sub>,...,b<sub>23</sub>
         *
         * The element b[0][3] is not used and always zero.
         */
        double b[3][4];
        /**
         * Array of three double arrays (sizeof = 3x4) which stores parameters
         * c<sub>0</sub>, c<sub>1</sub>, c<sub>2</sub>,
         * c<sub>10</sub>,...,c<sub>13</sub> and c<sub>20</sub>,...,c<sub>23</sub>
         * describing the p<sub>t</sub> distribution for gamma rays produced by
         * the Delta(1232) resonance.
         *
         * Layout:
         * - c[0][..] = c<sub>0</sub>, c<sub>1</sub>, c<sub>2</sub>
         * - c[1][..] = c<sub>10</sub>,...,c<sub>13</sub>
         * - c[2][..] = c<sub>20</sub>,...,c<sub>23</sub>
         *
         * The element c[0][3] is not used and always zero.
         */
        double c[3][4];
        
        PARAMSET params;
    };
    
    Cparamlib(void) = default;
    ~Cparamlib(void) = default;
    /* Implementation of Kamea model for for secondary electron/positron production in pp interactions. */
    double GetSigma(const double &Eg,const double &Ekinbullet, PARTICLE_IDS par);
    
    private:
    /* Electron antineutrino */
    void antinue_param(const double &Tp, PARAMSET* params);
    void antinue_param_nd(const double &Tp, PARAMSET* params);
    void antinue_param_diff(const double &Tp, PARAMSET* params);
    void antinue_param_delta(const double &Tp, PARAMSET* params);
    void antinue_param_res(const double &Tp, PARAMSET* params);
    /* Muon antineutrinos */
    void antinumu_param(const double &Tp, PARAMSET* params);
    void antinumu_param_nd(const double &Tp, PARAMSET* params);
    void antinumu_param_diff(const double &Tp, PARAMSET* params);
    void antinumu_param_delta(const double &Tp, PARAMSET* params);
    void antinumu_param_res(const double &Tp, PARAMSET* params);
    /* Electrons */
    void elec_param(const double &Tp, PARAMSET* params);
    void elec_param_nd(const double &Tp, PARAMSET* params);
    void elec_param_diff(const double &Tp, PARAMSET* params);
    void elec_param_delta(const double &Tp, PARAMSET* params);
    void elec_param_res(const double &Tp, PARAMSET* params);
    /* Gamma rays */
    void gamma_param(const double &Tp, PARAMSET* params);
    void gamma_pt_param(const double &E,const double &Tp, PARAMSET_PT* pt_params,const bool &flag);
    void gamma_param_nd(const double &Tp, PARAMSET* params);
    void gamma_param_diff(const double &Tp, PARAMSET* params);
    void gamma_param_delta(const double &Tp, PARAMSET* params);
    void gamma_param_res(const double &Tp, PARAMSET* params);
    void gamma_pt_param_nr(const double &E,const double &Tp, PARAMSET_PT* pt_params,const bool &flag);
    void gamma_pt_param_delta(const double &E,const double &Tp, PARAMSET_PT* pt_params,const bool &flag);
    void gamma_pt_param_res(const double &E,const double &Tp, PARAMSET_PT* pt_params,const bool &flag);
    /* Muon neutrinos */
    void numu_param(const double &Tp, PARAMSET* params);
    void numu_param_nd(const double &Tp, PARAMSET* params);
    void numu_param_diff(const double &Tp, PARAMSET* params);
    void numu_param_delta(const double &Tp, PARAMSET* params);
    void numu_param_res(const double &Tp, PARAMSET* params);
    /* Electron neutrinos */
    void nue_param(const double &Tp, PARAMSET* params);
    void nue_param_nd(const double &Tp, PARAMSET* params);
    void nue_param_diff(const double &Tp, PARAMSET* params);
    void nue_param_delta(const double &Tp, PARAMSET* params);
    void nue_param_res(const double &Tp, PARAMSET* params);
    /* Positrons */
    void posi_param(const double &Tp, PARAMSET* params);
    void posi_param_nd(const double &Tp, PARAMSET* params);
    void posi_param_diff(const double &Tp, PARAMSET* params);
    void posi_param_delta(const double &Tp, PARAMSET* params);
    void posi_param_res(const double &Tp, PARAMSET* params);
    /*
     * Calculation of the inclusive cross section, dsigma/dlogE, for any given secondary
     * particle and any given energy of that secondary particle
     *
     * Proton kinetic energy Tp is passed as well as a pointer to a struct holding the
     * parameters (at that Tp)
     */
    double sigma_incl_delta(const int &particle,const double &E,const double &Tp, PARAMSET* params);
    double sigma_incl_diff(const int &particle,const double &E,const double &Tp, PARAMSET* params);
    double sigma_incl_nd(const int &particle,const double &E,const double &Tp, PARAMSET* params);
    double sigma_incl_res(const int &particle,const double &E,const double &Tp, PARAMSET* params);
    double sigma_incl_tot(const int &particle,const double &E,const double &Tp, PARAMSET* params);
    /*
     * Calculation of the p-p cross section for a given proton momentum
     */
    double sigma_pp_delta(const double &Pp);
    double sigma_pp_diff(const double &Pp);
    double sigma_pp_nd(const double &Pp);
    double sigma_pp_res(const double &Pp);
    /*
     * Calculation of the differential cross section dsigma/dlogEdpT, for any gamma-ray energy
     *
     * Proton kinetic energy Tp is passed as well as a pointer to a struct keeping all the
     * parameters (at given Tp and log(E))
     */
    double sigma_pt_delta(const int &particle,const double &pT,const double &E,const double &Tp, PARAMSET_PT* pt_params);
    double sigma_pt_nr(const int &particle,const double &pT,const double &E,const double &Tp, PARAMSET_PT* pt_params);
    double sigma_pt_res(const int &particle,const double &pT,const double &E,const double &Tp, PARAMSET_PT* pt_params);
    double sigma_pt_tot(const int &particle,const double &pT,const double &E,const double &Tp, PARAMSET_PT* pt_params);
};

#endif
