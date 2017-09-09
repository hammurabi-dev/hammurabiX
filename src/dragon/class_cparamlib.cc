#include <iostream>
#include <cmath>

#include "class_cparamlib.h"

using namespace std;

double Cparamlib::GetSigma(const double &Eg, const double &Ekinbullet, PARTICLE_IDS par){
    PARAMSET *parameters;
    switch (par) {
        case ID_GAMMA :
            gamma_param(Ekinbullet, parameters);
            break;
        case ID_ELECTRON :
            elec_param(Ekinbullet, parameters);
            break;
        case ID_POSITRON :
            posi_param(Ekinbullet, parameters);
            break;
        default:
            cerr<<"ERR:"<<__FILE__
            <<" : in function "<<__func__<<endl
            <<" at line "<<__LINE__<<endl
            << "WRONG PARTICLE ID" << endl;
            exit(WRONGID);
    }
    
    return 1.e-3*sigma_incl_tot(par,Eg, Ekinbullet, parameters)/Eg;
}

/*Electron antineutrinos*/

/**
 * Calculates parameters a<sub>0</sub>,...,a<sub>8</sub> describing the
 * non-diffraction interaction electron antineutrino inclusive cross section as
 * a function of the proton kinetic energy T<sub>p</sub>.
 *
 * @param Tp     Proton kinetic energy in GeV.
 * @param params Pointer to a ::PARAMSET struct where the calculated parameters
 *               will be stored.
 */
void Cparamlib::antine_param_nd(const double &Tp, PARAMSET* params){
    double y, z;
    
    /* check whether params is a null pointer or not */
    if (params == NULL)
        return;
    
    y = log10(Tp*0.001);
    
    /* 06/06/06: removed unneccessary use of pow() to increase performance
     also added use of z = y + constant */
    if ((Tp > 0.487) && (Tp < 512000.1)) {
        z = y + 3.31;
        params->a[0] = 0.0013113 + z*(0.36538 + z*(1.5178 + z*(-0.20668 + 0.024255*z)));
        params->a[1] = -4.7833e-6 + 4.5837e-5*exp(-0.42980*(y + 3.4)) + 6.1559e-6/(y + 4.1731) + 1.1928e-6*y;
        params->a[2] = -245.22 + y*(73.223 + y*(-19.652 + y*(0.83138 + 0.71564*y)));
        params->a[3] = 0.45232 + y*(0.52934 + y*(0.010078 - 0.0017092*y));
        z = y + 3.32;
        params->a[4] = -0.0025734 + z*(0.38424 + z*(1.5517 + z*(0.17336 + z*(-0.17160 + 0.021059*z))));
        params->a[5] = 4.7673e-5 + 5.4936e-5*log10(0.0067905*(y + 4.3)) + 0.00020740/(y + 4.9772);
        params->a[6] = -270.30 - 114.47*log10(0.34352*(y + 3.4)) + y*(80.085 - 7.9240*y);
        params->a[7] = 3271.9 - 2.9161e5/(y + 87.847) - 6.2330*y*y;
        params->a[8] = -0.17787 + y*(0.36771 + y*(-0.025397 + y*(0.0019238 + 0.0032725*y)));
    } else {
        for (int i = 0; i != 9; ++i)
            params->a[i] = 0;
    }
}

/**
 * Calculates parameters b<sub>0</sub>,...,b<sub>7</sub> describing the
 * diffraction dissociation electron antineutrino inclusive cross section as a
 * function of the proton kinetic energy T<sub>p</sub>.
 *
 * @param Tp     Proton kinetic energy in GeV.
 * @param params Pointer to a ::PARAMSET struct where the calculated parameters
 *               will be stored.
 */
void Cparamlib::antinue_param_diff(const double &Tp, PARAMSET* params){
    double y, z1, z2, pow;
    
    /* check whether params is a null pointer or not */
    if (params == NULL)
        return;
    
    y = log10(Tp*0.001);
    
    /* 06/06/06: removed unneccessary use of pow() to increase performance
     also added use of z = y + constant and pow = <expression> */
    if ((Tp > 1.94) && (Tp < 512000.1)) {
        if (Tp > 11.0) {
            z1 = y + 0.55505;
            z2 = y + 9.2685;
            params->b[0] = 41.307*tanh(-0.37411*(y + 2.2)) - 4.1223*z1*z1 + 0.0042652*z2*z2*z2*z2;
            pow = (y + 1.9196)/(1.0 + 11.530*(y + 1.9196));
            params->b[1] = -132.50 + 142.12*exp(-8.0289*pow*pow);
            z1 = y + 250.77;
            params->b[2] = -17.223 + 0.011285*tanh(69.746*(y + 1.9)) - 0.048233*y + 0.00025881*z1*z1;
            pow = (y + 1.9292)/(1.0 + 16.682*(y + 1.9292));
            params->b[3] = 8.1991 - 9.6437*exp(-45.261*pow*pow);
        } else {
            params->b[0] = 0;
            params->b[1] = 0;
            params->b[2] = 0;
            params->b[3] = 0;
        }
        z1 = y + 2.95;
        pow = (y + 2.2) + 0.43867*(y + 2.2)*(y + 2.2);
        params->b[4] = 0.55919 + z1*(0.36647 + 0.056194*z1) + 0.49957*exp(-5.5317*pow*pow);
        params->b[5] = 1.2544 - 0.52362*tanh(2.7638*(y + 1.9)) - 0.055837*(y - 17.638);
        params->b[6] = 1.4788 + y*(1.0278 + y*(-0.092852 + y*(-0.0062734 + 0.011920*y)));
        z1 = y - 2.7889;
        params->b[7] = 5.1651 + 5.7398*tanh(-0.37356*(y + 2.1)) - 0.22234*z1*z1;
    } else {
        for (int i = 0; i != 8; ++i)
            params->b[i] = 0;
    }
}

/**
 * Calculates parameters c<sub>0</sub>,...,c<sub>4</sub> describing the
 * Delta(1232) electron antineutrino inclusive cross section as a function of
 * the proton kinetic energy T<sub>p</sub>. Since the negative charged
 * Delta(1232) resonance is not formed, no electron antineutrinos are produced
 * and thus this function only sets the ::PARAMSET::c array to zero.
 *
 * @param Tp     Proton kinetic energy in GeV.
 * @param params Pointer to a ::PARAMSET struct where the calculated parameters
 *               will be stored.
 */
void Cparamlib::antinue_param_delta(const double &Tp, PARAMSET* params){
    
    /* check whether params is a null pointer or not */
    if (params == NULL)
        return;
    
    for (int i = 0; i != 5; ++i)
        params->c[i] = 0;
}

/**
 * Calculates parameters d<sub>0</sub>,...,d<sub>4</sub> describing the
 * red(1600) electron antineutrino inclusive cross section as a function of
 * the proton kinetic energy T<sub>p</sub>.
 *
 * @param Tp     Proton kinetic energy in GeV.
 * @param params Pointer to a ::PARAMSET struct where the calculated parameters
 *               will be stored.
 */
void Cparamlib::antinue_param_res(const double &Tp, PARAMSET* params){
    double y, pow;
    
    /* check whether params is a null pointer or not */
    if (params == NULL)
        return;
    
    y = log10(Tp*0.001);
    
    /* 06/06/06: removed unneccessary use of pow() to increase performance
     also added use of pow = <expression> */
    if ((Tp < 0.69) || (Tp > 2.76)) {
        for (int i = 0; i != 5; ++i)
            params->d[i] = 0;
    } else {
        pow = (y + 2.9537)/(1.0 + 1.4320*(y + 2.9537));
        params->d[0] = 0.36459*exp(-58.210*pow*pow) - (0.11283 + 0.046244*y);
        params->d[1] = -9.5066 - y*(5.4655 + 0.31769*y);
        params->d[2] = -7.1831 - 7.1551*tanh(30.354*(y + 2.1)) + 0.33757*y;
        params->d[3] = 2.7938 + y*(1.6992 + 0.20161*y);
        params->d[4] = 0.61878 + y*(0.62371 + y*(0.18913 + 0.019118*y));
    }
}


/**
 * Calculates all parameters a<sub>0</sub>,...,a<sub>8</sub>,
 * b<sub>0</sub>,...,b<sub>7</sub>, c<sub>0</sub>,...,c<sub>4</sub> and
 * d<sub>0</sub>,...,d<sub>4</sub> describing electron antineutrino inclusive
 * cross sections as a function of the proton kinetic energy.
 *
 * @param Tp     Proton kinetic energy in GeV.
 * @param params Pointer to a ::PARAMSET struct where the calculated parameters
 *               will be stored.
 */
void Cparamlib::antinue_param(const double &Tp, PARAMSET* params){
    /* check whether params is a null pointer or not */
    if (params == NULL)
        return;
    
    antinue_param_nd(Tp, params);
    antinue_param_diff(Tp, params);
    antinue_param_delta(Tp, params);
    antinue_param_res(Tp, params);
}

/* Muon antineutrinos */

/**
 * Calculates parameters a<sub>0</sub>,...,a<sub>8</sub> describing the
 * non-diffraction interaction muon antineutrino inclusive cross section as
 * a function of the proton kinetic energy T<sub>p</sub>.
 *
 * @param Tp     Proton kinetic energy in GeV.
 * @param params Pointer to a ::PARAMSET struct where the calculated parameters
 *               will be stored.
 */
void Cparamlib::antinumu_param_nd(const double &Tp, PARAMSET* params){
    double y, z;
    
    /* check whether params is a null pointer or not */
    if (params == NULL)
        return;
    
    y = log10(Tp*0.001);
    
    /* 06/06/06: removed unneccessary use of pow() to increase performance
     also added use of z = y + constant */
    if ((Tp > 0.487) && (Tp < 512000.1)) {
        z = y + 3.3;
        params->a[0] = z*(-1.5243 + z*(10.107 + z*(-4.3126 + z*(0.80081 - 0.048724*z))));
        params->a[1] = -2.6297e-5 + 9.3858e-5*exp(-0.32384*(y + 3.4)) + 7.7821e-6/(y + 4.0560) + 7.6149e-6*y - 8.4091e-7*y*y;
        params->a[2] = -243.62 + 59.374*y - 5.7356*y*y + 1.9815*y*y*y - 1.0478*y*y*y*y;
        params->a[3] = 0.50807 + 0.60221*y + 0.0034120*y*y - 0.011139*y*y*y;
        z = y + 3.32;
        params->a[4] = z*(2.6483 + z*(4.4585 + z*(-1.2744 + z*(0.11659 + 0.0030477*z))));
        z = y + 4.7707;
        params->a[5] = 9.1101e-7 + 1.3880e-6*log10(0.016998*(y + 4.3)) + 0.00012583/(z*z);
        params->a[6] = -272.11 + 53.477*log10(0.35531*(y + 3.4)) + y*(56.041 - 6.0876*y);
        params->a[7] = 6431.8 + 893.92*log10(5.7013e-9*(y + 3.9)) + 2103.6/(y + 5.6740) - 6.1125*y*y;
        params->a[8] = -0.11120 + y*(0.38144 + y*(-0.040128 + y*(0.0047484 + 0.0054707*y)));
    } else {
        for (int i = 0; i != 9; i++)
            params->a[i] = 0;
    }
}

/**
 * Calculates parameters b<sub>0</sub>,...,b<sub>7</sub> describing the
 * diffraction dissociation muon antineutrino inclusive cross section as a
 * function of the proton kinetic energy T<sub>p</sub>.
 *
 * @param Tp     Proton kinetic energy in GeV.
 * @param params Pointer to a ::PARAMSET struct where the calculated parameters
 *               will be stored.
 */
void Cparamlib::antinumu_param_diff(const double &Tp, PARAMSET* params){
    double y, z1, z2, pow;
    
    /* check whether params is a null pointer or not */
    if (params == NULL)
        return;
    
    y = log10(Tp*0.001);
    
    /* 06/06/06: removed unneccessary use of pow() to increase performance
     also added use of z = y + constant and pow = <expression> */
    if ((Tp > 1.94) && (Tp < 512000.1)) {
        if (Tp > 11.0) {
            z1 = y + 0.52273;
            z2 = y + 9.5266;
            params->b[0] = 70.430*tanh(-0.35816*(y + 2.2)) - 6.6796*z1*z1 + 0.0065659*z2*z2*z2*z2;
            pow = (y + 2.2190)/(1.0 + 81.105*(y + 2.2190));
            params->b[1] = -8.1145 + 7686.0*exp(-44046*pow*pow);
            z1 = y - 1.8683;
            params->b[2] = -1.3095 + 0.071270*tanh(0.0075463*(y + 1.9)) + 0.067759*(y + 5.3433) - 0.0044205*z1*z1;
            pow = (y + 2.8363)/(1.0 + 7.0976*(y + 2.8363));
            params->b[3] = 0.082149 - 2190.1*exp(-533.75*pow*pow);
        } else {
            params->b[0] = 0;
            params->b[1] = 0;
            params->b[2] = 0;
            params->b[3] = 0;
        }
        z1 = y + 2.95;
        pow = y + 2.28 - 0.18922*(y + 2.2)*(y + 2.2);
        params->b[4] = 2.7540 + 0.33859*z1*z1 - 0.0034274*z1*z1*z1*z1 + 1.1679*exp(-10.408*pow*pow);
        params->b[5] = 2.1817 - 0.59584*tanh(2.7054*(y + 1.9)) - 0.010909*(y - 14.880);
        params->b[6] = 1.4591 + y*(1.0275 + y*(-0.074949 + y*(-0.0060396 + 0.0097568*y)));
        z1 = y - 2.7653;
        params->b[7] = 3.7609 + 4.2843*tanh(-0.37148*(y + 2.1)) - 0.16479*z1*z1;
    } else {
        for (int i = 0; i != 8; i++)
            params->b[i] = 0;
    }
}

/**
 * Calculates parameters d<sub>0</sub>,...,d<sub>4</sub> describing the
 * red(1600) muon antineutrino inclusive cross section as a function of
 * the proton kinetic energy T<sub>p</sub>.
 *
 * @param Tp     Proton kinetic energy in GeV.
 * @param params Pointer to a ::PARAMSET struct where the calculated parameters
 *               will be stored.
 */
void Cparamlib::antinumu_param_delta(const double &Tp, PARAMSET* params){
    double y, p;
    
    /* check whether params is a null pointer or not */
    if (params == NULL)
        return;
    
    y = log10(Tp*0.001);
    
    /* 06/06/06: removed unneccessary use of pow() to increase performance
     also added use of pow = <expression> */
    if ((Tp < 0.488) || (Tp > 1.95)) {
        for (int i = 0; i != 5; i++)
            params->c[i] = 0;
    } else {
        p = (y + 3.1250)/(1.0 - 0.47567*(y + 3.1250));
        params->c[0] = 2.8262*exp(-62.894*p*p) + 5.6845 + 13.409/y - 0.097296*y*y;
        params->c[1] = 16.721 + 11.750*y + 2.4637*y*y;
        params->c[2] = -6.0557 - 6.3378*tanh(21.984*(y + 2.1)) + 0.43173*y;
        params->c[3] = 0.37009 + 0.27706*y;
        params->c[4] = 0.047507 + 0.061570*y + 0.0070117*y*y;
    }
}

/**
 * Calculates parameters d<sub>0</sub>,...,d<sub>4</sub> describing the
 * red(1600) muon antineutrino inclusive cross section as a function of
 * the proton kinetic energy T<sub>p</sub>.
 *
 * @param Tp     Proton kinetic energy in GeV.
 * @param params Pointer to a ::PARAMSET struct where the calculated parameters
 *               will be stored.
 */
void Cparamlib::antinumu_param_res(const double &Tp, PARAMSET* params){
    double y, pow;
    
    /* check whether params is a null pointer or not */
    if (params == NULL)
        return;
    
    y = log10(Tp*0.001);
    
    /* 06/06/06: removed unneccessary use of pow() to increase performance
     also added use of pow = <expression> */
    if ((Tp < 0.69) || (Tp > 2.76)) {
        for (int i = 0; i != 5; i++)
            params->d[i] = 0;
    } else {
        pow = (y + 2.9492)/(1.0 + 1.2994*(y + 2.9492));
        params->d[0] = 2.2400*exp(-57.159*pow*pow) - (0.66521 + 0.27554*y);
        params->d[1] = -7.0650 - 4.2773*y - 0.17648*y*y;
        params->d[2] = -7.0410 - 7.1977*tanh(31.095*(y + 2.1)) + 0.40238*y;
        params->d[3] = -1.2354 - 0.87581*y - 0.20829*y*y;
        params->d[4] = -0.11395 + 0.34418*y + 0.27103*y*y + 0.050248*y*y*y;
    }
}

/**
 * Calculates all parameters a<sub>0</sub>,...,a<sub>8</sub>,
 * b<sub>0</sub>,...,b<sub>7</sub>, c<sub>0</sub>,...,c<sub>4</sub> and
 * d<sub>0</sub>,...,d<sub>4</sub> describing muon antineutrino inclusive
 * cross sections as a function of the proton kinetic energy.
 *
 * @param Tp     Proton kinetic energy in GeV.
 * @param params Pointer to a ::PARAMSET struct where the calculated parameters
 *               will be stored.
 */
void Cparamlib::antinumu_param(const double &Tp, PARAMSET* params){
    /* check whether params is a null pointer or not */
    if (params == NULL)
        return;
    
    antinumu_param_nd(Tp, params);
    antinumu_param_diff(Tp, params);
    antinumu_param_delta(Tp, params);
    antinumu_param_res(Tp, params);
}

/* Electrons */
/**
 * Calculates parameters a<sub>0</sub>,...,a<sub>8</sub> describing the
 * nondiffraction interaction electron inclusive cross section as a function
 * of the proton kinetic energy T<sub>p</sub>.
 *
 * @param Tp     Proton kinetic energy in GeV.
 * @param params Pointer to a ::PARAMSET struct where the calculated parameters
 *               will be stored.
 */
void Cparamlib::elec_param_nd(const double &Tp, PARAMSET* params){
    double y, z;
    
    /* check whether params is a null pointer or not */
    if (params == NULL)
        return;
    
    y = log10(Tp*0.001);
    
    /* 06/06/06: removed unneccessary use of pow() to increase performance
     also added use of z = y + constant */
    if ((Tp > 0.487) && (Tp < 512000.1)) {
        z = y + 3.3;
        params->a[0] = z*(-0.018639 + z*(2.4315 + z*(-0.57719 + 0.063435*z)));
        params->a[1] = 7.1827e-6 + y*(-3.5067e-6 + y*(1.3264e-6 + y*(-3.3481e-7 + y*(2.3551e-8 + 3.4297e-9*y))));
        z = y + 7.9031;
        params->a[2] = 563.91 - 362.18*log10(2.7187*(y + 3.4)) - 2.8924e4/(z*z);
        params->a[3] = 0.52684 + y*(0.57717 + y*(0.0045336 - 0.0089066*y));
        z = y + 3.32;
        params->a[4] = z*(0.36108 + z*(1.6963 + z*(-0.074456 + z*(-0.071455 + 0.010473*z))));
        params->a[5] = 9.7387e-5 + 7.8573e-5*log10(0.0036055*(y + 4.3)) + 0.00024660/(y + 4.9390) - 3.8097e-7*y*y;
        params->a[6] = -273.00 - 106.22*log10(0.34100*(y + 3.4)) + 89.037*y - 12.546*y*y;
        z = y + 8.5518;
        params->a[7] = 432.53 - 883.99*log10(0.19737*(y + 3.9)) - 4.1938e4/(z*z);
        params->a[8] = -0.12756 + y*(0.43478 + y*(-0.0027797 - 0.0083074*y));
    } else {
        for (int i = 0; i != 9; i++)
            params->a[i] = 0;
    }
}

/**
 * Calculates parameters b<sub>0</sub>,...,b<sub>7</sub> describing the
 * diffraction dissociation electron inclusive cross section as a function
 * of the proton kinetic energy T<sub>p</sub>.
 *
 * @param Tp     Proton kinetic energy in GeV.
 * @param params Pointer to a ::PARAMSET struct where the calculated parameters
 *               will be stored.
 */
void Cparamlib::elec_param_diff(const double &Tp, PARAMSET* params){
    double y, z1, z2, pow;
    
    /* check whether params is a null pointer or not */
    if (params == NULL)
        return;
    
    y = log10(Tp*0.001);
    
    /* 06/06/06: removed unneccessary use of pow() to increase performance
     also added use of z = y + constant and pow = <expression> */
    if ((Tp > 1.94) && (Tp < 512000.1)){
        if (Tp > 5.51) {
            z1 = y + 1.6878;
            z2 = y + 9.6400;
            params->b[0] = 0.20463*tanh(-6.2370*(y + 2.2)) - 0.16362*z1*z1 + 3.5183e-4*z2*z2*z2*z2;
            pow = (y + 2.0154)/(1.0 + 0.62779*(y + 2.0154));
            params->b[1] = 1.6537 + 3.8530*exp(-3.2027*pow*pow);
            z1 = y + 256.63;
            params->b[2] = -10.722 + 0.082672*tanh(1.8879*(y + 2.1)) + 0.00014895*z1*z1;
            pow = (y + 1.9877)/(1.0 + 0.40300*(y + 1.988));
            params->b[3] = -0.023752 - 0.51734*exp(-3.3087*pow*pow);
        } else {
            params->b[0] = 0;
            params->b[1] = 0;
            params->b[2] = 0;
            params->b[3] = 0;
        }
        z1 = y + 2.9;
        params->b[4] = 0.94921 + 0.12280*z1*z1 - 7.1585e-4*z1*z1*z1*z1 + 0.52130*log10(z1);
        params->b[5] = -4.2295 - 1.0025*tanh(9.0733*(y + 1.9)) - 0.11452*(y - 62.382);
        params->b[6] = 1.4862 + y*(0.99544 + y*(-0.042763 + y*(-0.0040065 + 0.0057987*y)));
        z1 = y - 2.8542;
        params->b[7] = 6.2629 + 6.9517*tanh(-0.36480*(y + 2.1)) - 0.26033*z1*z1;
    } else {
        for (int i = 0; i != 8; i++)
            params->b[i] = 0;
    }
}

/**
 * Calculates parameters c<sub>0</sub>,...,c<sub>4</sub> describing the
 * Delta(1232) electron inclusive cross section as a function of the proton
 * kinetic energy T<sub>p</sub>. Since the negative charged Delta(1232)
 * resonance is not formed, no electrons are produced and thus this function
 * only sets the ::PARAMSET::c array to zero.
 *
 * @param Tp     Proton kinetic energy in GeV.
 * @param params Pointer to a ::PARAMSET struct where the calculated parameters
 *               will be stored.
 */
void Cparamlib::elec_param_delta(const double &Tp, PARAMSET* params){
    
    /* check whether params is a null pointer or not */
    if (params == NULL)
        return;
    
    for (int i = 0; i != 5; i++)
        params->c[i] = 0;
}

/**
 * Calculates parameters d<sub>0</sub>,...,d<sub>4</sub> describing the
 * res(1600) electron inclusive cross section as a function of the proton
 * kinetic energy T<sub>p</sub>.
 *
 * @param Tp     Proton kinetic energy in GeV.
 * @param params Pointer to a ::PARAMSET struct where the calculated parameters
 *               will be stored.
 */
void Cparamlib::elec_param_res(const double &Tp, PARAMSET* params){
    double y, pow;
    
    /* check whether params is a null pointer or not */
    if (params == NULL)
        return;
    
    y = log10(Tp*0.001);
    
    /* 06/06/06: removed unneccessary use of pow() to increase performance
     also added use of pow = <expression> */
    if ((Tp < 0.69) || (Tp > 2.76)) {
        for (int i = 0; i != 5; i++)
            params->d[i] = 0;
    } else {
        pow = (y + 2.9537)/(1.0 + 1.5221*(y + 2.9537));
        params->d[0] = 0.37790*exp(-56.826*pow*pow) - 0.059458 + 0.0096583*y*y;
        params->d[1] = -5.5135 - 3.3988*y;
        params->d[2] = -7.1209 - 7.1850*tanh(30.801*(y + 2.1)) + 0.35108*y;
        params->d[3] = -6.7841 - 4.8385*y - 0.91523*y*y;
        params->d[4] = -134.03 - 139.63*y - 48.316*y*y - 5.5526*y*y*y;
    }
}

/**
 * Calculates all parameters a<sub>0</sub>,...,a<sub>8</sub>,
 * b<sub>0</sub>,...,b<sub>7</sub>, c<sub>0</sub>,...,c<sub>4</sub> and
 * d<sub>0</sub>,...,d<sub>4</sub> describing electron inclusive cross
 * sections as a function of the proton kinetic energy.
 *
 * @param Tp     Proton kinetic energy in GeV.
 * @param params Pointer to a ::PARAMSET struct where the calculated parameters
 *               will be stored.
 */
void Cparamlib::elec_param(const double &Tp, PARAMSET* params){
    /* check whether params is a null pointer or not */
    if (params == NULL)
        return;
    
    elec_param_nd(Tp, params);
    elec_param_diff(Tp, params);
    elec_param_delta(Tp, params);
    elec_param_res(Tp, params);
}

/* Gamma rays */

/**
 * Calculates parameters a<sub>0</sub>,...,a<sub>8</sub> describing the
 * nondiffraction interaction gamma-ray inclusive cross section as a function
 * of the proton kinetic energy T<sub>p</sub>.
 *
 * @param Tp     Proton kinetic energy in GeV.
 * @param params Pointer to a ::PARAMSET struct where the calculated parameters
 *               will be stored.
 */
void Cparamlib::gamma_param_nd(const double &Tp, PARAMSET* params){
    double y, z;
    
    /* check whether params is a null pointer or not */
    if (params == NULL)
        return;
    
    y = log10(Tp*0.001);
    
    /* 06/06/06: removed unneccessary use of pow() to increase performance
     also added use of z = y + constant */
    if ((Tp > 0.487) && (Tp < 512000.1)) {
        z = y + 3.3;
        params->a[0] = z*(-0.51187 + z*(7.6179 + z*(-2.1332 + 0.22184*z)));
        params->a[1] = -1.2592e-5 + 1.4439e-5*exp(-0.29360*(y + 3.4)) + 5.9363e-5/(y + 4.1485) + y*(2.2640e-6 - 3.3723e-7*y);
        params->a[2] = -174.83 + 152.78*log10(1.5682*(y + 3.4)) - 808.74/(y + 4.6157);
        params->a[3] = 0.81177 + y*(0.56385 + y*(0.0040031 + y*(-0.0057658 + 0.00012057*y)));
        z = y + 3.32;
        params->a[4] = z*(0.68631 + z*(10.145 + z*(-4.6176 + z*(0.86824 - 0.053741*z))));
        z = y + 4.7171;
        params->a[5] = 9.0466e-7 + 1.4539e-6*log10(0.015204*(y + 3.4)) + 0.00013253/(z*z) + y*(-4.1228e-7 + 2.2036e-7*y);
        params->a[6] = -339.45 + 618.73*log10(0.31595*(y + 3.9)) + 250.20/((y + 4.4395)*(y + 4.4395));
        params->a[7] = -35.105 + y*(36.167 + y*(-9.3575 + 0.33717*y));
        params->a[8] = 0.17554 + y*(0.37300 + y*(-0.014938 + y*(0.0032314 + 0.0025579*y)));
    } else {
        for (int i = 0; i != 9; i++)
            params->a[i] = 0;
    }
}

/**
 * Calculates parameters b<sub>0</sub>,...,b<sub>7</sub> describing the
 * diffraction dissociation gamma-ray inclusive cross section as a function
 * of the proton kinetic energy T<sub>p</sub>.
 *
 * @param Tp     Proton kinetic energy in GeV.
 * @param params Pointer to a ::PARAMSET struct where the calculated parameters
 *               will be stored.
 */
void Cparamlib::gamma_param_diff(const double &Tp, PARAMSET* params){
    double y, z1, z2, pow;
    
    /* check whether params is a null pointer or not */
    if (params == NULL)
        return;
    
    y = log10(Tp*0.001);
    
    /* 06/06/06: removed unneccessary use of pow() to increase performance
     also added use of z = y + constant and pow = <expression> */
    if ((Tp > 1.94) && (Tp < 512000.1)) {
        if (Tp > 5.51) {
            z1 = y + 0.59913;
            z2 = y + 9.4773;
            params->b[0] = 60.142*tanh(-0.37555*(y + 2.2)) - 5.9564*z1*z1 + 0.0060162*z2*z2*z2*z2;
            z1 = y + 369.13;
            params->b[1] = 35.322 + 3.8026*tanh(-2.4979*(y + 1.9)) - 0.00021870*z1*z1;
            z1 = y + 252.43;
            params->b[2] = -15.732 - 0.082064*tanh(-1.9621*(y + 2.1)) + 0.00023355*z1*z1;
            pow = (y + 1.0444)/(1.0 + 0.27437*(y + 1.0444));
            params->b[3] = -0.086827 + 0.37646*exp(-0.53053*pow*pow);
        } else {
            params->b[0] = 0;
            params->b[1] = 0;
            params->b[2] = 0;
            params->b[3] = 0;
        }
        z1 = y + 2.95;
        pow = (y + 2.45) - 0.19717*(y + 2.45)*(y + 2.45);
        params->b[4] = 2.5982 + 0.39131*z1*z1 - 0.0049693*z1*z1*z1*z1 + 0.94131*exp(-24.347*pow*pow);
        z1 = (y - 0.83562)/(1.0 + 0.33933*(y - 0.83562));
        params->b[5] = 0.11198 + y*(-0.64582 + 0.16114*y) + 2.2853*exp(-0.0032432*z1*z1);
        params->b[6] = 1.7843 + y*(0.91914 + y*(0.050118 + y*(0.038096 + y*(-0.027334 + y*(-0.0035556 + 0.0025742*y)))));
        z1 = y + 1.8441;
        params->b[7] = -0.19870 + y*(-0.071003 + 0.019328*y) - 0.28321*exp(-6.0516*z1*z1);
    } else {
        for (int i = 0; i != 8; i++)
            params->b[i] = 0;
    }
}

/**
 * Calculates parameters c<sub>0</sub>,...,c<sub>4</sub> describing the
 * Delta(1232) gamma-ray inclusive cross section as a function of the proton
 * kinetic energy T<sub>p</sub>.
 *
 * @param Tp     Proton kinetic energy in GeV.
 * @param params Pointer to a ::PARAMSET struct where the calculated parameters
 *               will be stored.
 */
void Cparamlib::gamma_param_delta(const double &Tp, PARAMSET* params){
    double y, pow;
    
    /* check whether params is a null pointer or not */
    if (params == NULL)
        return;
    
    y = log10(Tp*0.001);
    
    /* 06/06/06: removed unneccessary use of pow() to increase performance
     also added use of pow = <expression> */
    if ((Tp < 0.488) || (Tp > 1.95))
        for (int i = 0; i != 5; i++)
            params->c[i] = 0;
    else {
        pow = ((y + 3.1301)/(1.0 + 0.14921*(y + 3.1301)));
        params->c[0] = 2.4316*exp(-69.484*pow*pow) - (6.3003 + 9.5349/y - 0.38121*y*y);
        params->c[1] = 56.872 + y*(40.627 + 7.7528*y);
        params->c[2] = -5.4918 - 6.7872*tanh(4.7128*(y + 2.1)) + 0.68048*y;
        params->c[3] = -0.36414 + 0.039777*y;
        params->c[4] = -0.72807 + y*(-0.48828 - 0.092876*y);
    }
}

/**
 * Calculates parameters d<sub>0</sub>,...,d<sub>4</sub> describing the
 * res(1600) gamma-ray inclusive cross section as a function of the proton
 * kinetic energy T<sub>p</sub>.
 *
 * @param Tp     Proton kinetic energy in GeV.
 * @param params Pointer to a ::PARAMSET struct where the calculated parameters
 *               will be stored.
 */
void Cparamlib::gamma_param_res(const double &Tp, PARAMSET* params){
    double y, pow;
    
    /* check whether params is a null pointer or not */
    if (params == NULL)
        return;
    
    y = log10(Tp*0.001);
    
    /* 06/06/06: removed unneccessary use of pow() to increase performance
     also added use of pow = <expression> */
    if ((Tp < 0.69) || (Tp > 2.76)) {
        for (int i = 0; i != 5; i++)
            params->d[i] = 0;
    } else {
        pow = ((y + 2.9507)/(1.0 + 1.2912*(y + 2.9507)));
        params->d[0] = 3.2433*exp(-57.133*pow*pow) - (1.0640 + 0.43925*y);
        params->d[1] = 16.901 + y*(5.9539 + y*(-2.1257 - 0.92057*y));
        params->d[2] = -6.6638 - 7.5010*tanh(30.322*(y + 2.1)) + 0.54662*y;
        params->d[3] = -1.50648 + y*(-0.87211 - 0.17097*y);
        params->d[4] = 0.42795 + y*(0.55136 + y*(0.20707 + 0.027552*y));
    }
}

/**
 * Calculates all parameters a<sub>0</sub>,...,a<sub>8</sub>,
 * b<sub>0</sub>,...,b<sub>7</sub>, c<sub>0</sub>,...,c<sub>4</sub> and
 * d<sub>0</sub>,...,d<sub>4</sub> describing gamma-ray inclusive cross
 * sections as a function of the proton kinetic energy.
 *
 * @param Tp     Proton kinetic energy in GeV.
 * @param params Pointer to a ::PARAMSET struct where the calculated parameters
 *               will be stored.
 */
void Cparamlib::gamma_param(const double &Tp, PARAMSET* params){
    /* check whether params is a null pointer or not */
    if (params == NULL)
        return;
    
    gamma_param_nd(Tp, params);
    gamma_param_diff(Tp, params);
    gamma_param_delta(Tp, params);
    gamma_param_res(Tp, params);
}

/**
 * Calculates parameters a<sub>0</sub>, a<sub>1</sub> and
 * a<sub>10</sub>,...,a<sub>13</sub> describing the nonresonance gamma-ray
 * p<sub>t</sub> distribution as a function of the proton kinetic energy
 * T<sub>p</sub> and the gamma-ray energy E.
 *
 * @param E         Gamma-ray energy in GeV.
 * @param Tp        Proton kinetic energy in GeV.
 * @param pt_params Pointer to a ::PARAMSET_PT struct where the calculated
 *                  parameters will be stored.
 * @param flag      Flag if Tp has changed since last call to this function. If
 *                  set to 0 then only a<sub>0</sub> and a<sub>1</sub> are
 *                  recalculated using the new gamma-ray energy E.
 */
void Cparamlib::gamma_pt_param_nr(const double &E,const double &Tp, PARAMSET_PT* pt_params,const bool &flag){
    double a14;
    double x, xp, xa, pow;
    double y;
    double sigma_incl;
    
    /* check whether pt_params is a null pointer or not */
    if (pt_params == NULL)
        return;
    
    x = log10(E);
    y = log10(Tp*0.001);
    
    if ((Tp > 0.579) && (Tp < 512000.1)) {
        /* first calculate a10,...,a14 using functions from Table 1
         but only do this if flag = 1 */
        if (flag) {
            pt_params->a[2] = 0.043775 + 0.010271*exp(-0.55808*y);
            pt_params->a[3] = 0.8;
            pt_params->a[4] = 0.34223 + 0.027134*y - 0.0089229*y*y + 4.9996e-4*y*y*y;
            pt_params->a[5] = -0.20480 + 0.013372*y + 0.13087*exp(0.0044021*(y - 11.467)*(y - 11.467));
            
            gamma_param_nd(Tp, &(pt_params->params));
            gamma_param_diff(Tp, &(pt_params->params));
        }
        
        /* then calculate a1 using a10,...,a14; when x > -0.75 coefficent a14 = a1(x = -0.75) */
        if (x > -0.75) {
            xa = -0.75 + pt_params->a[4];
            a14 = pt_params->a[2]*exp(-pt_params->a[3]*xa*xa);
            
            pt_params->a[1] = pt_params->a[5]*(x + 0.75) + a14;
        } else {
            xa = x + pt_params->a[4];
            pt_params->a[1] = pt_params->a[2]*exp(-pt_params->a[3]*xa*xa);
        }
        
        sigma_incl = sigma_incl_nd(ID_GAMMA, E, Tp, &(pt_params->params));
        sigma_incl += sigma_incl_diff(ID_GAMMA, E, Tp, &(pt_params->params));
        
        /* finally calculate a0 from a1 and sigma_incl*/
        if (pt_params->a[1] != 0) {
            pt_params->a[0] = 1/(pt_params->a[1]*pt_params->a[1])*sigma_incl;
        } else {
            pt_params->a[0] = 0;
        }
    } else {
        for (int i = 0; i != 6; i++) {
            pt_params->a[i] = 0;
        }
    }
}

/**
 * Calculates parameters b<sub>0</sub>, b<sub>1</sub>, b<sub>2</sub>,
 * b<sub>10</sub>,...,b<sub>13</sub> and b<sub>20</sub>,...,b<sub>23</sub>
 * describing the Delta(1232) gamma-ray p<sub>t</sub> distribution as a
 * function of the proton kinetic energy T<sub>p</sub> and the gamma-ray
 * energy E.
 *
 * @param E         Gamma-ray energy in GeV.
 * @param Tp        Proton kinetic energy in GeV.
 * @param pt_params Pointer to a ::PARAMSET_PT struct where the calculated
 *                  parameters will be stored.
 * @param flag      Flag if Tp has changed since last call to this function. If
 *                  set to 0 then only b<sub>0</sub>,...,b<sub>2</sub> are
 *                  recalculated using the new gamma-ray energy E.
 */
void Cparamlib::gamma_pt_param_delta(const double &E,const double &Tp, PARAMSET_PT* pt_params,const bool &flag){
    double x, xb;
    double y;
    double A;
    double pow;
    double sigma_incl;
    
    /* check whether pt_params is a null pointer or not */
    if (pt_params == NULL)
        return;
    
    x = log10(E);
    y = log10(Tp*0.001);
    
    if ((Tp < 0.488) || (Tp > 1.95)) {
        for (int j = 0; j != 3; j++) {
            for (int i = 0; i != 4; i++) {
                pt_params->b[j][i] = 0;
                pt_params->b[j][i] = 0;
            }
        }
    } else {
        /* first calculate b10,...,b14 and b20,...,b24 using functions from Table 1
         but only do this if flag = 1 */
        if (flag) {
            pt_params->b[1][0] = 18.712 + 18.030*y + 5.8239*y*y + 0.62728*y*y*y;
            pt_params->b[1][1] = 612.61 + 404.80*y + 67.406*y*y;
            pt_params->b[1][2] = 98.639 + 96.741*y + 31.597*y*y + 3.4567*y*y*y;
            pt_params->b[1][3] = -208.38 - 183.65*y - 53.283*y*y - 5.0470*y*y*y;
            
            pt_params->b[2][0] = 0.21977 + 0.064073*y;
            pt_params->b[2][1] = 3.3187e3 + 3.4634e3*y + 1.1982e3*y*y + 136.71*y*y*y;
            pt_params->b[2][2] = 91.410 + 91.613*y + 30.621*y*y + 3.4296*y*y*y;
            pt_params->b[2][3] = -521.40 - 529.06*y - 178.49*y*y - 19.975*y*y*y;
            
            gamma_param_delta(Tp, &(pt_params->params));
        }
        
        /* then calculate b1 and b2 using b10,...,b24
         when x > 0.5 bi = 0; i=0,1,2 */
        A = 0.81*(y + 3.32) - 0.5;
        if (x < A) {
            xb = x - pt_params->b[1][2];
            pow = xb/(1.0 + pt_params->b[1][3]*xb);
            pt_params->b[0][1] = pt_params->b[1][0]*exp(-pt_params->b[1][1]*pow*pow);
            
            xb = x - pt_params->b[2][2];
            pow = xb/(1.0 + pt_params->b[2][3]*xb);
            pt_params->b[0][2] = pt_params->b[2][0]*exp(-pt_params->b[2][1]*pow*pow);
            
            sigma_incl = sigma_incl_delta(ID_GAMMA, E, Tp, &(pt_params->params));
            
            /* calculate b0 from b1 and b2 */
            if ((pt_params->b[0][1] != 0) && (pt_params->b[0][2] != 0)) {
                pt_params->b[0][0] = 2*sigma_incl/(pt_params->b[0][1]*sqrt(M_PI*pt_params->b[0][2])*(erf(pt_params->b[0][1]/sqrt(pt_params->b[0][2])) + 1) +
                                                   pt_params->b[0][2]*exp(-pt_params->b[0][1]*pt_params->b[0][1]/pt_params->b[0][2]));
            } else {
                pt_params->b[0][0] = 0;
            }
        } else {
            for (int i = 0; i != 4; i++) {
                pt_params->b[0][i] = 0;
            }
        }
    }
}

/**
 * Calculates parameters c<sub>0</sub>, c<sub>1</sub>, c<sub>2</sub>,
 * c<sub>10</sub>,...,c<sub>13</sub> and c<sub>10</sub>,...,c<sub>13</sub>
 * describing the res(1600) gamma-ray p<sub>t</sub> distribution as a
 * function of the proton kinetic energy T<sub>p</sub> and the gamma-ray
 * energy.
 *
 * @param E         Gamma-ray energy in GeV.
 * @param Tp        Proton kinetic energy in GeV.
 * @param pt_params Pointer to a ::PARAMSET_PT struct where the calculated
 *                  parameters will be stored.
 * @param flag      Flag if Tp has changed since last call to this function. If
 *                  set to 0 then only c<sub>0</sub>,...,c<sub>2</sub> are
 *                  recalculated using the new gamma-ray energy E.
 */
void Cparamlib::gamma_pt_param_res(const double &E,const double &Tp, PARAMSET_PT* pt_params,const bool &flag){
    double x, xc;
    double y;
    double A;
    double pow;
    double sigma_incl;
    
    /* check whether pt_params is a null pointer or not */
    if (pt_params == NULL)
        return;
    
    x = log10(E);
    y = log10(Tp*0.001);
    
    if ((Tp < 0.69) || (Tp > 2.76)) {
        for (int j = 0; j != 3; j++) {
            for (int i = 0; i != 4; i++) {
                pt_params->c[j][i] = 0;
                pt_params->c[j][i] = 0;
            }
        }
    } else {
        /* first calculate c10,...,c14 and c20,...,c24 using functions from Table 1
         but only do this if flag = 1 */
        if (flag) {
            pt_params->c[1][0] = -1.5013 - 1.1281*y - 0.19813*y*y;
            pt_params->c[1][1] = -33.179 - 22.496*y - 3.3108*y*y;
            pt_params->c[1][2] = 116.44 + 122.11*y + 42.594*y*y + 4.9609*y*y*y;
            pt_params->c[1][3] = -545.77 - 574.80*y - 201.25*y*y - 23.400*y*y*y;
            
            pt_params->c[2][0] = 0.68849 + 0.36438*y + 0.047958*y*y;
            pt_params->c[2][1] = -1.6871e4 - 1.7412e4*y - 5.9648e3*y*y - 679.27*y*y*y;
            pt_params->c[2][2] = -88.565 - 94.034*y - 33.014*y*y - 3.8205*y*y*y;
            pt_params->c[2][3] = 1.5141e3 + 1.5757e3*y + 544.20*y*y + 62.446*y*y*y;
            
            gamma_param_res(Tp, &(pt_params->params));
        }
        
        /* then calculate c1 and c2 using c10,...,c24
         when x > 0.5 ci = 0; i=0,1,2 */
        A = 0.82*(x + 3.17) - 0.25;
        if (x < 0.5) {
            xc = x - pt_params->c[1][2];
            pow = xc/(1.0 + pt_params->c[1][3]*xc);
            pt_params->c[0][1] = pt_params->c[1][0]*exp(-pt_params->c[1][1]*pow*pow);
            
            xc = x - pt_params->c[2][2];
            pow = xc/(1.0 + pt_params->c[2][3]*xc);
            pt_params->c[0][2] = pt_params->c[2][0]*exp(-pt_params->c[2][1]*pow*pow);
            
            sigma_incl = sigma_incl_res(ID_GAMMA, E, Tp, &(pt_params->params));
            
            /* calculate c0 from c1 and c2 */
            if ((pt_params->c[0][1] != 0) && (pt_params->c[0][2] != 0)) {
                pt_params->c[0][0] = 2*sigma_incl/(pt_params->c[0][1]*sqrt(M_PI*pt_params->c[0][2])*(erf(pt_params->c[0][1]/sqrt(pt_params->c[0][2])) + 1) +
                                                   pt_params->c[0][2]*exp(-pt_params->c[0][1]*pt_params->c[0][1]/pt_params->c[0][2]));
            } else {
                pt_params->c[0][0] = 0;
            }
        } else {
            for (int i = 0; i != 4; i++) {
                pt_params->c[0][i] = 0;
            }
        }
    }
}

/**
 * Calculates all parameters describing the gamma-ray p<sub>t</sub> distribution
 * as a function of the proton kinetic energy T<sub>p</sub> and the gamma-ray
 * energy.
 *
 * @param E         Gamma-ray energy in GeV.
 * @param Tp        Proton kinetic energy in GeV.
 * @param pt_params Pointer to a ::PARAMSET_PT struct where the calculated
 *                  parameters will be stored.
 * @param flag      Flag if Tp has changed since last call to this function.
 */
void Cparamlib::gamma_pt_param(const double &E,const double &Tp, PARAMSET_PT* pt_params,const bool &flag){
    /* check whether params is a null pointer or not */
    if (pt_params == NULL)
        return;
    
    gamma_pt_param_nr(E, Tp, pt_params, flag);
    gamma_pt_param_delta(E, Tp, pt_params, flag);
    gamma_pt_param_res(E, Tp, pt_params, flag);
}

/* Electron neutrinos */
/**
 * Calculates parameters a<sub>0</sub>,...,a<sub>8</sub> describing the
 * non-diffraction interaction electron neutrino inclusive cross section as a
 * function of the proton kinetic energy T<sub>p</sub>.
 *
 * @param Tp     Proton kinetic energy in GeV.
 * @param params Pointer to a ::PARAMSET struct where the calculated parameters
 *               will be stored.
 */
void Cparamlib::nue_param_nd(const double &Tp, PARAMSET* params){
    double y, z;
    
    /* check whether params is a null pointer or not */
    if (params == NULL)
        return;
    
    y = log10(Tp*0.001);
    
    /* 06/06/06: removed unneccessary use of pow() to increase performance
     also added use of z = y + constant */
    if ((Tp > 0.487) && (Tp < 512000.1)) {
        z = y + 3.31;
        params->a[0] = 0.0074087 + z*(2.9161 + z*(0.99061 + z*(-0.28694 + 0.038799*z)));
        params->a[1] = -3.2480e-5 + 7.1944e-5*exp(-0.21814*(y + 3.4)) + 2.0467e-5/(y + 4.1640) + y*(5.6954e-6 - 3.4105e-7*y);
        params->a[2] = -230.50 + y*(58.802 + y*(-9.9393 + y*(1.2473 - 0.26322*y)));
        params->a[3] = 0.45064 + y*(0.56930 + y*(0.012428 - 0.0070889*y));
        z = y + 3.32;
        params->a[4] = -0.011883 + z*(1.7992 + z*(3.5264 + z*(-1.7478 + z*(0.32077 - 0.017667*z))));
        z = y + 4.8229;
        params->a[5] = -1.6238e-7 + 1.8116e-6*exp(-0.30111*(y + 3.4)) + 9.6112e-5/(z*z);
        params->a[6] = -261.30 - 43.351*log10(0.35298*(y + 3.4)) + y*(70.925 - 8.7147*y);
        params->a[7] = 184.45 - 1473.6/(y + 6.8788) - 4.0536*y*y;
        params->a[8] = -0.24019 + y*(0.38504 + y*(0.0096869 - 0.0015046*y));
    } else {
        for (int i = 0; i != 9; i++)
            params->a[i] = 0;
    }
}

/**
 * Calculates parameters b<sub>0</sub>,...,b<sub>7</sub> describing the
 * diffraction dissociation electron neutrino inclusive cross section as a
 * function of the proton kinetic energy T<sub>p</sub>.
 *
 * @param Tp     Proton kinetic energy in GeV.
 * @param params Pointer to a ::PARAMSET struct where the calculated parameters
 *               will be stored.
 */
void Cparamlib::nue_param_diff(const double &Tp, PARAMSET* params){
    double y, z1, z2, pow;
    
    /* check whether params is a null pointer or not */
    if (params == NULL)
        return;
    
    y = log10(Tp*0.001);
    
    /* 06/06/06: removed unneccessary use of pow() to increase performance
     also added use of z = y + constant and pow = <expression> */
    if ((Tp > 1.94) && (Tp < 512000.1)) {
        if (Tp > 11.0) {
            z1 = y + 0.76010;
            z2 = y + 8.5075;
            params->b[0] = 53.809*tanh(-0.41421*(y + 2.2)) - 6.7538*z1*z1 + 0.0088080*z2*z2*z2*z2;
            pow = (y + 1.8901)/(1.0 + 4.4440*(y + 1.8901));
            params->b[1] = -50.211 + 55.131*exp(-1.3651*pow*pow);
            z1 = y + 250.68;
            params->b[2] = -17.231 + 0.041100*tanh(7.9638*(y + 1.9)) - 0.055449*y + 0.00025866*z1*z1;
            pow = (y + 1.8998)/(1.0 + 5.5969*(y + 1.8998));
            params->b[3] = 12.335 - 12.893*exp(-1.4412*pow*pow);
        } else {
            params->b[0] = 0;
            params->b[1] = 0;
            params->b[2] = 0;
            params->b[3] = 0;
        }
        z1 = y + 2.95;
        pow = y + 2.2 + 4.6121*(y + 2.2)*(y + 2.2);
        params->b[4] = 1.3558 + 0.46601*z1 + 0.052978*z1*z1 + 0.79575*exp(-5.4007*pow*pow);
        params->b[5] = 1.8756 - 0.42169*tanh(1.6100*(y + 1.9)) - 0.051026*(y - 3.9573);
        params->b[6] = 1.5016 + 1.0118*y - 0.072787*y*y - 0.0038858*y*y*y + 0.0093650*y*y*y*y;
        z1 = y - 2.8604;
        params->b[7] = 4.9735 + 5.5674*tanh(-0.36249*(y + 2.1)) - 0.20660*z1*z1;
    } else {
        for (int i = 0; i != 8; i++)
            params->b[i] = 0;
    }
}

/**
 * Calculates parameters c<sub>0</sub>,...,c<sub>4</sub> describing the
 * Delta(1232) electron neutrino inclusive cross section as a function of the
 * proton kinetic energy T<sub>p</sub>.
 *
 * @param Tp     Proton kinetic energy in GeV.
 * @param params Pointer to a ::PARAMSET struct where the calculated parameters
 *               will be stored.
 */
void Cparamlib::nue_param_delta(const double &Tp, PARAMSET* params){
    double y, pow;
    
    /* check whether params is a null pointer or not */
    if (params == NULL)
        return;
    
    y = log10(Tp*0.001);
    
    /* 06/06/06: removed unneccessary use of pow() to increase performance
     also added use of pow = <expression> */
    if ((Tp < 0.488) || (Tp > 1.95)) {
        for (int i = 0; i != 5; i++)
            params->c[i] = 0;
    } else {
        pow = ((y + 3.1282)/(1.0 + 0.48420*(y + 3.1282)));
        params->c[0] = 2.8290*exp(-71.339*pow*pow) - (9.6339 + 15.733/y - 0.52413*y*y);
        params->c[1] = -24.571 - 15.831*y - 2.1200*y*y;
        params->c[2] = -5.9593 - 6.4695*tanh(4.7225*(y + 2.1)) + 0.50003*y;
        params->c[3] = 0.26022 + 0.24545*y;
        params->c[4] = 0.076498 + 0.061678*y + 0.0040028*y*y;
    }
}

/**
 * Calculates parameters d<sub>0</sub>,...,d<sub>4</sub> describing the
 * res(1600) electron neutrino inclusive cross section as a function of the
 * proton kinetic energy T<sub>p</sub>.
 *
 * @param Tp     Proton kinetic energy in GeV.
 * @param params Pointer to a ::PARAMSET struct where the calculated parameters
 *               will be stored.
 */
void Cparamlib::nue_param_res(const double &Tp, PARAMSET* params){
    double y, pow;
    
    /* check whether params is a null pointer or not */
    if (params == NULL)
        return;
    
    y = log10(Tp*0.001);
    
    /* 06/06/06: removed unneccessary use of pow() to increase performance
     also added use of pow = <expression> */
    if ((Tp < 0.69) || (Tp > 2.76)) {
        for (int i = 0; i != 5; i++)
            params->d[i] = 0;
    } else {
        pow = (y + 2.9509)/(1.0 + 1.4101*(y + 2.9509));
        params->d[0] = 1.7951*exp(-57.260*pow*pow) - (0.58604 + 0.23868*y);
        params->d[1] = -2.6395 - 1.5105*y + 0.22174*y*y;
        params->d[2] = -7.0512 - 7.1970*tanh(31.074*(y + 2.1)) + 0.39007*y;
        params->d[3] = -1.4271 - 1.0399*y - 0.24179*y*y;
        params->d[4] = 0.74875 + 0.63616*y + 0.17396*y*y + 0.017636*y*y*y;
    }
}

/**
 * Calculates all parameters a<sub>0</sub>,...,a<sub>8</sub>,
 * b<sub>0</sub>,...,b<sub>7</sub>, c<sub>0</sub>,...,c<sub>4</sub> and
 * d<sub>0</sub>,...,d<sub>4</sub> describing electron neutrino inclusive cross
 * sections as a function of the proton kinetic energy.
 *
 * @param Tp     Proton kinetic energy in GeV.
 * @param params Pointer to a ::PARAMSET struct where the calculated parameters
 *               will be stored.
 */
void Cparamlib::nue_param(const double &Tp, PARAMSET* params){
    /* check whether params is a null pointer or not */
    if (params == NULL)
        return;
    
    nue_param_nd(Tp, params);
    nue_param_diff(Tp, params);
    nue_param_delta(Tp, params);
    nue_param_res(Tp, params);
}

/* Muon neutrinos */
/**
 * Calculates parameters a<sub>0</sub>,...,a<sub>8</sub> describing the
 * non-diffraction interaction muon neutrino inclusive cross section as a
 * function of the proton kinetic energy T<sub>p</sub>.
 *
 * @param Tp     Proton kinetic energy in GeV.
 * @param params Pointer to a ::PARAMSET struct where the calculated parameters
 *               will be stored.
 */
void Cparamlib::numu_param_nd(const double &Tp, PARAMSET* params){
    double y, z;
    
    /* check whether params is a null pointer or not */
    if (params == NULL)
        return;
    
    y = log10(Tp*0.001);
    
    /* 06/06/06: removed unneccessary use of pow() to increase performance
     also added use of z = y + constant */
    if ((Tp > 0.487) && (Tp < 512000.1)) {
        z = y + 3.3;
        params->a[0] = z*(-0.63611 + z*(9.9015 + z*(-4.5897 + z*(0.91778 - 0.060724*z))));
        params->a[1] = 6.8700e-6 + y*(-2.8245e-6 + y*(7.6032e-7 + y*(-3.2953e-7 + 7.4292e-8*y)));
        params->a[2] = -240.46 + y*(58.405 + y*(-9.8556 + y*(3.1401 - 0.88932*y)));
        params->a[3] = 0.49935 + y*(0.60919 + y*(0.0024963 - 0.0099910*y));
        z = y + 3.32;
        params->a[4] = z*(2.5094 + z*(4.1350 + z*(-0.89534 + z*(- 2.7577e-3 + 0.014511*z))));
        z = y + 4.7136;
        params->a[5] = 8.2046e-7 + 1.4085e-6*log10(0.016793*(y + 4.3)) + 0.00013340/(z*z);
        params->a[6] = -267.55 - 0.21018*log10(0.35217*(y + 3.4)) + y*(69.586 - 9.9930*y);
        params->a[7] = 2741.8 + 222.01*log10(9.7401e-13*(y + 3.9)) - 4772.5/(y + 19.773) - 6.1001*y*y;
        params->a[8] = -0.11857 + y*(0.39072 + y*(-0.037813 + y*(0.0022265 + 0.0046931*y)));
    } else {
        for (int i = 0; i != 9; i++)
            params->a[i] = 0;
    }
}

/**
 * Calculates parameters b<sub>0</sub>,...,b<sub>7</sub> describing the
 * diffraction dissociation muon neutrino inclusive cross section as a
 * function of the proton kinetic energy T<sub>p</sub>.
 *
 * @param Tp     Proton kinetic energy in GeV.
 * @param params Pointer to a ::PARAMSET struct where the calculated parameters
 *               will be stored.
 */
void Cparamlib::numu_param_diff(const double &Tp, PARAMSET* params){
    double y, z1, z2, pow;
    
    /* check whether params is a null pointer or not */
    if (params == NULL)
        return;
    
    y = log10(Tp*0.001);
    
    /* 06/06/06: removed unneccessary use of pow() to increase performance
     also added use of z = y + constant and pow = <expression> */
    if ((Tp > 1.94) && (Tp < 512000.1)) {
        if (Tp > 11.0) {
            z1 = y + 0.44754;
            z2 = y + 9.9165;
            params->b[0] = 64.682*tanh(-0.34313*(y + 2.2)) - 5.5955*z1*z1 + 0.0050117*z2*z2*z2*z2;
            pow = (y + 2.3066)/(1.0 + 41.612*(y + 2.3066));
            params->b[1] = -7.6016 + 3042.7*exp(-1.0134e4*pow*pow);
            z1 = y - 1.8861;
            params->b[2] = -1.4978 - 0.58163*tanh(-0.36488*(y + 1.9)) + 0.031825*(y + 2.8097) + 0.022796*z1*z1;
            pow = (y + 3.8835)/(1.0 + 0.53343*(y + 3.8835));
            params->b[3] = -0.0061483 - 65.799*exp(-4.8239*pow*pow);
        } else {
            params->b[0] = 0;
            params->b[1] = 0;
            params->b[2] = 0;
            params->b[3] = 0;
        }
        z1 = y + 2.95;
        pow = y + 2.28 - 0.19149*(y + 2.28)*(y + 2.28);
        params->b[4] = 2.8009 + z1*z1*(0.35351 - 0.0039779*z1*z1) + 1.3012*exp(-10.592*pow*pow);
        params->b[5] = 1.8016 - 0.69847*tanh(2.8627*(y + 1.9)) - 0.015722*(y - 45.410);
        params->b[6] = 1.4617 + y*(1.0167 + y*(-0.078617 + y*(-0.0038336 + 0.010141*y)));
        z1 = y - 2.4209;
        params->b[7] = 3.5599 + 4.0041*tanh(-0.41889*(y + 2.1)) - 0.18182*z1*z1;
    } else {
        for (int i = 0; i != 8; i++)
            params->b[i] = 0;
    }
}

/**
 * Calculates parameters c<sub>0</sub>,...,c<sub>4</sub> describing the
 * Delta(1232) muon neutrino inclusive cross section as a function of the proton
 * kinetic energy T<sub>p</sub>.
 *
 * @param Tp     Proton kinetic energy in GeV.
 * @param params Pointer to a ::PARAMSET struct where the calculated parameters
 *               will be stored.
 */
void Cparamlib::numu_param_delta(const double &Tp, PARAMSET* params){
    double y, pow;
    
    /* check whether params is a null pointer or not */
    if (params == NULL)
        return;
    
    y = log10(Tp*0.001);
    
    /* 06/06/06: removed unneccessary use of pow() to increase performance
     also added use of pow = <expression> */
    if ((Tp < 0.488) || (Tp > 1.95)) {
        for (int i = 0; i != 5; i++)
            params->c[i] = 0;
    } else {
        pow = (y + 3.1278)/(1.0 - 0.19497*(y + 3.1278));
        params->c[0] = 3.6052*exp(-60.914*pow*pow) - (0.92514 - 2.1315/y - 0.23548*y*y);
        params->c[1] = 95.310 + y*(70.497 + 13.636*y);
        params->c[2] = -6.2158 - 6.2939*tanh(21.592*(y + 2.1)) + 0.37440*y;
        params->c[3] = 2.7485 + 1.1692*y;
        params->c[4] = -2.7568 + y*(-1.8461 - 0.31376*y);
    }
}

/**
 * Calculates parameters d<sub>0</sub>,...,d<sub>4</sub> describing the
 * res(1600) muon neutrino inclusive cross section as a function of the proton
 * kinetic energy T<sub>p</sub>.
 *
 * @param Tp     Proton kinetic energy in GeV.
 * @param params Pointer to a ::PARAMSET struct where the calculated parameters
 *               will be stored.
 */
void Cparamlib::numu_param_res(double Tp, PARAMSET* params){
    double y, pow;
    
    /* check whether params is a null pointer or not */
    if (params == NULL)
        return;
    
    y = log10(Tp*0.001);
    
    /* 06/06/06: removed unneccessary use of pow() to increase performance
     also added use of pow = <expression> */
    if ((Tp < 0.69) || (Tp > 2.76)) {
        for (int i = 0; i != 5; i++)
            params->d[i] = 0;
    } else {
        pow = (y + 2.9509)/(1.0 + 1.3154*(y + 2.9509));
        params->d[0] = 2.5489*exp(-58.488*pow*pow) - (0.83039 + 0.34412*y);
        params->d[1] = 88.173 + 65.148*y + 12.585*y*y;
        params->d[2] = -7.0962 - 7.1690*tanh(30.890*(y + 2.1)) + 0.38032*y;
        params->d[3] = -4.1440 + y*(-3.2717 - 0.70537*y);
        params->d[4] = 2.2624 + y*(1.1806 + y*(0.0043450 - 0.043020*y));
    }
}

/**
 * Calculates all parameters a<sub>0</sub>,...,a<sub>8</sub>,
 * b<sub>0</sub>,...,b<sub>7</sub>, c<sub>0</sub>,...,c<sub>4</sub> and
 * d<sub>0</sub>,...,d<sub>4</sub> describing muon neutrino inclusive cross
 * sections as a function of the proton kinetic energy.
 *
 * @param Tp     Proton kinetic energy in GeV.
 * @param params Pointer to a ::PARAMSET struct where the calculated parameters
 *               will be stored.
 */
void Cparamlib::numu_param(const double &Tp, PARAMSET* params){
    /* check whether params is a null pointer or not */
    if (params == NULL)
        return;
    
    numu_param_nd(Tp, params);
    numu_param_diff(Tp, params);
    numu_param_delta(Tp, params);
    numu_param_res(Tp, params);
}

/* Positrons */
/**
 * Calculates parameters a<sub>0</sub>,...,a<sub>8</sub> describing the
 * non-diffraction interaction positron inclusive cross section as a function
 * of the proton kinetic energy T<sub>p</sub>.
 *
 * @param Tp     Proton kinetic energy in GeV.
 * @param params Pointer to a ::PARAMSET struct where the calculated parameters
 *               will be stored.
 */
void Cparamlib::posi_param_nd(const double &Tp, PARAMSET* params){
    double y, z;
    
    /* check whether params is a null pointer or not */
    if (params == NULL)
        return;
    
    y = log10(Tp*0.001);
    
    /* 06/06/06: removed unneccessary use of pow() to increase performance
     also added use of z = y + constant */
    if ((Tp > 0.487) && (Tp < 512000.1)) {
        z = y + 3.3;
        params->a[0] = z*(-0.79606 + z*(7.7496 + z*(-3.9326 + z*(0.80202 - 0.054994*z))));
        params->a[1] = 6.7943e-6 + y*(-3.5345e-6 + y*(6.0927e-7 + y*(2.0219e-7 + y*(5.1005e-8 - 4.2622e-8*y))));
        params->a[2] = 44.827 + 81.378*log10(0.027733*(y + 3.5)) - 1.3886e4/((y + 8.4417)*(y + 8.4417));
        params->a[3] = 0.52010 + y*(0.59336 + y*(0.012032 - 0.0064242*y));
        z = y + 3.32;
        params->a[4] = z*(2.1361 + z*(1.8514 + z*(-0.47872 + z*(0.0032043 + 0.0082955*z))));
        params->a[5] = 1.0845e-6 + 1.4336e-6*log10(0.0077255*(y + 4.3)) + 0.00013018/((y + 4.8188)*(y + 4.8188)) + 9.3601e-8*y;
        params->a[6] = -267.74 + 14.175*log10(0.35391*(y + 3.4)) + y*(64.669 - 7.7036*y);
        params->a[7] = 138.26 - 529.84*log10(0.12467*(y + 3.9)) - 1.9869e4/((y + 7.6884)*(y + 7.6884)) + 1.0675*y*y;
        params->a[8] = -0.14707 + y*(0.40135 + y*(0.0039899 - 0.0016602*y));
    } else {
        for (int i = 0; i != 9; i++)
            params->a[i] = 0;
    }
}

/**
 * Calculates parameters b<sub>0</sub>,...,b<sub>7</sub> describing the
 * diffraction dissociation positron inclusive cross section as a function
 * of the proton kinetic energy T<sub>p</sub>.
 *
 * @param Tp     Proton kinetic energy in GeV.
 * @param params Pointer to a ::PARAMSET struct where the calculated parameters
 *               will be stored.
 */
void Cparamlib::posi_param_diff(const double &Tp, PARAMSET* params){
    double y, z1, z2, pow;
    
    /* check whether params is a null pointer or not */
    if (params == NULL)
        return;
    
    y = log10(Tp*0.001);
    
    /* 06/06/06: removed unneccessary use of pow() to increase performance
     also added use of z = y + constant and pow = <expression> */
    if ((Tp > 1.94) && (Tp < 512000.1)) {
        if (Tp > 11.0) {
            z1 = y + 0.67500;
            z2 = y + 9.0824;
            params->b[0] = 29.192*tanh(-0.37879*(y + 2.2)) - 3.2196*z1*z1 + 3.6687e-3*z2*z2*z2*z2;
            pow = (y + 1.8781)/(1.0 + 3.8389*(y + 1.8781));
            params->b[1] = -142.97 + 147.86*exp(-0.37194*pow*pow);
            z1 = y + 234.65;
            params->b[2] = -14.487 - 4.2223*tanh(-13.546*(y + 2.2)) + 0.00016988*z1*z1;
            pow = (y + 1.8194)/(1.0 + 0.99946*(y + 1.8194));
            params->b[3] = -0.0036974 - 0.41976*exp(-6.1527*pow*pow);
        } else {
            params->b[0] = 0;
            params->b[1] = 0;
            params->b[2] = 0;
            params->b[3] = 0;
        }
        z1 = y + 2.95;
        pow = y + 2.29 - 0.18967*(y + 2.29);
        params->b[4] = 1.8108 + z1*z1*(0.18545 - 2.0049e-3*z1*z1) + 0.85084*exp(-14.987*pow*pow);
        params->b[5] = 2.0404 - 0.51548*tanh(2.2758*(y + 1.9)) - 0.035009*(y - 6.6555);
        params->b[6] = 1.5258 + y*(1.0132 + y*(-0.064388 + y*(-0.0040209 + 0.0082772*y)));
        z1 = y - 2.7718;
        params->b[7] = 3.0551 + 3.5240*tanh(-0.36739*(y + 2.1)) - 0.13382*z1*z1;
    } else {
        for (int i = 0; i < 8; i++)
            params->b[i] = 0;
    }
}

/**
 * Calculates parameters c<sub>0</sub>,...,c<sub>4</sub> describing the
 * Delta(1232) positron inclusive cross section as a function of the proton
 * kinetic energy T<sub>p</sub>.
 *
 * @param Tp     Proton kinetic energy in GeV.
 * @param params Pointer to a ::PARAMSET struct where the calculated parameters
 *               will be stored.
 */
void Cparamlib::posi_param_delta(const double &Tp, PARAMSET* params){
    double y, pow;
    
    /* check whether params is a null pointer or not */
    if (params == NULL)
        return;
    
    y = log10(Tp*0.001);
    
    /* 06/06/06: removed unneccessary use of pow() to increase performance
     also added use of pow = <expression> */
    if ((Tp < 0.488) || (Tp > 1.95))
        for (int i = 0; i != 5; i++)
            params->c[i] = 0;
    else {
        pow = ((y + 3.1272)/(1.0 + 0.22831*(y + 3.1272)));
        params->c[0] = 2.9841*exp(-67.857*pow*pow) - (6.5855 + 9.6984/y - 0.41256*y*y);
        params->c[1] = 6.8276 + 5.2236*y + 1.4630*y*y;
        params->c[2] = -6.0291 - 6.4581*tanh(5.0830*(y + 2.1)) + 0.46352*y;
        params->c[3] = 0.59300 + 0.36093*y;
        params->c[4] = 0.77368 + 0.44776*y + 0.056409*y*y;
    }
}

/**
 * Calculates parameters d<sub>0</sub>,...,d<sub>4</sub> describing the
 * res(1600) positron inclusive cross section as a function of the proton
 * kinetic energy T<sub>p</sub>.
 *
 * @param Tp     Proton kinetic energy in GeV.
 * @param params Pointer to a ::PARAMSET struct where the calculated parameters
 *               will be stored.
 */
void Cparamlib::posi_param_res(const double &Tp, PARAMSET* params){
    double y, pow;
    
    /* check whether params is a null pointer or not */
    if (params == NULL)
        return;
    
    y = log10(Tp*0.001);
    
    /* 06/06/06: removed unneccessary use of pow() to increase performance
     also added use of pow = <expression> */
    if ((Tp < 0.69) || (Tp >= 2.76)) {
        for (int i = 0; i != 5; i++)
            params->d[i] = 0;
    } else {
        pow = ((y + 2.9485)/(1.0 + 1.2892*(y + 2.9485)));
        params->d[0] = 1.9186*exp(-56.544*pow*pow) - (0.23720 - 0.041315*y*y);
        params->d[1] = -4.9866 - 3.1435*y;
        params->d[2] = -7.0550 - 7.2165*tanh(31.033*(y + 2.1)) + 0.38541*y;
        params->d[3] = -2.8915 - 2.1495*y - 0.45006*y*y;
        params->d[4] = -1.2970 - 0.13947*y + 0.41197*y*y + 0.10641*y*y*y;
    }
}

/**
 * Calculates all parameters a<sub>0</sub>,...,a<sub>8</sub>,
 * b<sub>0</sub>,...,b<sub>7</sub>, c<sub>0</sub>,...,c<sub>4</sub> and
 * d<sub>0</sub>,...,d<sub>4</sub> describing positron inclusive cross
 * sections as a function of the proton kinetic energy.
 *
 * @param Tp     Proton kinetic energy in GeV.
 * @param params Pointer to a ::PARAMSET struct where the calculated parameters
 *               will be stored.
 */
void Cparamlib::posi_param(const double &Tp, PARAMSET* params){
    /* check whether params is a null pointer or not */
    if (params == NULL)
        return;
    
    posi_param_nd(Tp, params);
    posi_param_diff(Tp, params);
    posi_param_delta(Tp, params);
    posi_param_res(Tp, params);
}

/* Differential Cross Section calculaor */
/*
 *
 * @brief Main part of cparamlib; functions to calculate of differential
 * cross sections.
 *
 * File provides functions for the calculation of inclusive cross sections for
 * stable secondary particles (gamma rays, electrons, positrons, electron
 * neutrinos, electron antineutrinos, muon neutrinos and muon antineturinos)
 * and transverse momentum distributions for gamma rays, as described in
 * <a href="http://www.arxiv.org/abs/astro-ph/0605581">Kamae et al. (2006)</a>
 * and
 * <a href="http://www.arxiv.org/abs/0709.0233">Karlsson and Kamae (2008)</a>.
 * It also provides functions for the calculation of the inelastic
 * proton-proton cross section given in
 * <a href="http://www.arxiv.org/abs/astro-ph/0605581">Kamae et al. (2006)</a>.
 *
 *
 * Coefficients for the kinematic cutoff functions as given in Table 2 of
 * Kamae et al. (2006)
 *
 */
const double L_MAX[7] = {0.96, 0.96, 0.94, 0.98, 0.98, 0.94, 0.98};
const double W_NDL[7] = {15.0, 20.0, 15.0, 15.0, 15.0, 20.0, 15.0};
const double W_NDH[7] = {44.0, 45.0, 47.0, 42.0, 40.0, 45.0, 40.0};

/**
 * Calculates the inclusive cross section for the non-diffractive interaction
 * (equation 6 of Kamae et al. 2006). The inclusive cross section is defined as
 * dsigma/dlog(E) and returned in the units of mb.
 *
 * @param particle Secondary particle id number (see cparamlib.h).
 * @param E        Secondary particle energy in GeV.
 * @param Tp       Proton kinetic energy in GeV.
 * @param params   Pointer to a ::PARAMSET struct. The struct should be
 *                 initialized before being passed to this function.
 * @return         Inclusive cross section dsigma/dlogE in mb.
 *
 */
double Cparamlib::sigma_incl_nd(const int &particle,const double &E,const double &Tp, PARAMSET* params){
    double Wl, Wh, Lmin, Lmax;
    double x, y;
    double xa3, xa8;
    double pow1, pow2;
    double sigma;
    double cutoff, r_factor;
    
    /* check whether params is a null pointer or not */
    if (params == NULL)
        return 0;
    
    /* calculate log(E) and log(Tp) */
    x = log10(E);
    y = log10(Tp*0.001);
    
    /* init some variables, given in table 2 */
    Lmin = -2.6;
    Lmax = L_MAX[particle]*(y + 3.0);
    Wl = W_NDL[particle];
    Wh = W_NDH[particle];
    
    /* calculate the flux due to non-diffractive process for given gamma-ray energy */
    /* 06/12/06: replaced code involving pow() */
    xa3 = x - params->a[3];
    pow1 = xa3*(1 + params->a[2]*xa3);
    xa8 = x - params->a[8];
    pow2 = xa8*(1 + xa8*(params->a[6] + params->a[7]*xa8));
    sigma = params->a[0]*exp(-params->a[1]*pow1*pow1) + params->a[4]*exp(-params->a[5]*pow2*pow2);
    
    /* cutoff is the kinematic limit function as in the paper */
    cutoff = (1/(1 + exp(Wl*(Lmin - x))))*(1/(1 + exp(Wh*(x - Lmax))));
    sigma = sigma*cutoff;
    
    if (sigma < 0)
        sigma = 0;
    
    /* renormalization
     this is different for each particle, thus we must use if statements
     */
    r_factor = 1;
    switch (particle) {
            /* gamma */
        case ID_GAMMA:
            if (Tp <= 1.95) {
                pow1 = (y + 3.25)/(1 + 8.08*(y + 3.25));
                r_factor = 3.05*exp(-107.0*pow1*pow1);
            } else
                r_factor = 1.01;
            break;
            /* electron */
        case ID_ELECTRON:
            if (Tp <= 15.6) {
                pow1 = (y + 3.26)/(1 + 9.21*(y + 3.26));
                r_factor = 3.63*exp(-106*pow1*pow1) + y*(-0.182 - 0.175*y);
            } else
                r_factor = 1.01;
            break;
            /* positron */
        case ID_POSITRON:
            if (Tp <= 5.52) {
                pow1 = (y + 3.25)/(1 + 10.4*(y + 3.25));
                r_factor = 2.22*exp(-98.9*pow1*pow1);
            }
            break;
            /* electron neutrino */
        case ID_NUE:
            if (Tp <= 7.81) {
                pow1 = (y + 3.26)/(1 + 6.56*(y + 3.26));
                r_factor = 0.329*exp(-247*pow1*pow1) + y*(-0.959 - 0.229*y);
            }
            break;
            /* muon neutrino */
        case ID_NUMU:
            if (Tp <= 15.6) {
                pow1 = (y + 3.25)/(1 + 8.38*(y + 3.25));
                r_factor = 2.23*exp(-93.4*pow1*pow1) + y*(-0.376 - 0.121*y);
            }
            break;
            /* electron antineutrino */
        case ID_ANTINUE:
            if (Tp <= 15.6) {
                pow1 = (y + 3.27)/(1 + 6.59*(y + 3.27));
                r_factor = 2.67*exp(-45.7*pow1*pow1) + y*(-0.301 - 0.208*y);
            }
            break;
            /* muon antineutrino */
        case ID_ANTINUMU:
            if (Tp <= 15.6) {
                pow1 = (y + 3.25)/(1 + 8.34*(y + 3.25));
                r_factor = 2.56*exp(-107*pow1*pow1) + y*(-0.385 - 0.125*y);
            }
            break;
    }
    sigma = sigma*r_factor;
    
    return sigma;
}

/**
 * Calculates the inclusive cross section for the diffraction dissociation
 * process (equation 9 of Kamae et al. 2006). The inclusive cross section is
 * defined as dsigma/dlog(E) and is returned in units of mb.
 *
 * @param particle Secondary particle id number (see cparamlib.h).
 * @param E        Secondary particle energy in GeV.
 * @param Tp       Proton kinetic energy in GeV.
 * @param params   Pointer to a ::PARAMSET struct. The struct should be
 *                 initialized before being passed to this function.
 * @return         Inclusive cross section dsigma/dlogE in mb.
 *
 */
double Cparamlib::sigma_incl_diff(const int &particle,const double &E,const double &Tp, PARAMSET* params){
    double Wdiff, Lmax;
    double x, y;
    double pow1, pow2;
    double sigma;
    double cutoff;
    
    /* check whether params is a null pointer or not */
    if (params == NULL)
        return 0;
    
    /* calculate log(E) and log(Tp) */
    x = log10(E);
    y = log10(Tp*0.001);
    
    /* init some variables */
    Lmax = y + 3.0;
    Wdiff = 75.0;
    
    /* calculate the sigma due to diffractive process for given gamma-ray energy */
    /* 06/12/06: replaced code involving pow() */
    pow1 = (x - params->b[2])/(1 + params->b[3]*(x - params->b[2]));
    pow2 = (x - params->b[6])/(1 + params->b[7]*(x - params->b[6]));
    sigma = params->b[0]*exp(-params->b[1]*pow1*pow1) + params->b[4]*exp(-params->b[5]*pow2*pow2);
    
    /* cutoff is the kinematic limit function as in the paper */
    cutoff = 1/(1 + exp(Wdiff*(x - Lmax)));
    sigma = sigma*cutoff;
    
    if (sigma < 0)
        sigma = 0;
    
    return sigma;
}

/**
 * Calculates the inclusive cross section for the Delta(1232) resonance
 * (equation 12 of Kamae et al. 2006). The inclusive cross section is defined as
 * dsigma/dlog(E) and is returned in units of mb.
 *
 * @param particle Secondary particle id number (see cparamlib.h)
 * @param E        Secondary particle energy in GeV
 * @param Tp       Proton kinetic energy in GeV
 * @param params   Pointer to a ::PARAMSET struct. The struct should be
 *                 initialized before being passed to this function.
 * @return         Inclusive cross section dsigma/dlogE in mb
 *
 */
double Cparamlib::sigma_incl_delta(const int &particle,const double &E,const double &Tp, PARAMSET* params){
    double Wdiff, Lmax;
    double x, y;
    double xc2;
    double pow;
    double sigma;
    double cutoff;
    
    /* check whether params is a null pointer or not */
    if (params == NULL)
        return 0;
    
    /* calculate log(E) and log(Tp) */
    x = log10(E);
    y = log10(Tp*0.001);
    
    /* init some variables */
    Lmax = y + 3.0;
    Wdiff = 75.0;
    
    /* calculate the sigma due to resonance process for given gamma-ray energy */
    /* 06/12/06: replaced code involving pow() */
    xc2 = x - params->c[2];
    pow = xc2/(1 + xc2*(params->c[3] + params->c[4]*xc2));
    sigma = params->c[0]*exp(-params->c[1]*pow*pow);
    
    /* cutoff is the kinematic limit function as in the paper */
    cutoff = 1/(1 + exp(Wdiff*(x - Lmax)));
    sigma = sigma*cutoff;
    
    if (sigma < 0)
        sigma = 0;
    
    return sigma;
}

/**
 * Calculates the inclusive cross section for the res(1600) resonance
 * (equation 12 of Kamae et al. 2006). The inclusive cross section is defined as
 * dsigma/dlog(E) and is returned in units of mb.
 *
 * @param particle Secondary particle id number (see cparamlib.h).
 * @param E        Secondary particle energy in GeV.
 * @param Tp       Proton kinetic energy in GeV.
 * @param params   Pointer to a ::PARAMSET struct. The struct should be
 *                 initialized before being passed to this function.
 * @return         Inclusive cross section dsigma/dlogE in mb.
 *
 */
double Cparamlib::sigma_incl_res(const int &particle,const double &E,const double &Tp, PARAMSET* params){
    double Wdiff, Lmax;
    double x, y;
    double xd2;
    double pow;
    double sigma;
    double cutoff;
    
    /* check whether params is a null pointer or not */
    if (params == NULL)
        return 0;
    
    /* calculate log(E) and log(Tp) */
    x = log10(E);
    y = log10(Tp*0.001);
    
    /* init some variables */
    Lmax = y + 3.0;
    Wdiff = 75.0;
    
    /* calculate the sigma due to resonance process for given gamma-ray energy */
    /* 06/12/06: replaced code involving pow() */
    xd2 = x - params->d[2];
    pow = xd2/(1 + xd2*(params->d[3] + params->d[4]*xd2));
    sigma = params->d[0]*exp(-params->d[1]*pow*pow);
    
    /* cutoff is the kinematic limit function as in the paper */
    cutoff = 1/(1 + exp(Wdiff*(x - Lmax)));
    sigma = sigma*cutoff;
    
    if (sigma < 0)
        sigma = 0;
    
    return sigma;
}

/**
 * Calculates the total inclusive cross section, i.e. the sum of all components
 * (non-diffractive interaction, diffraction dissociation, Delta(1232) and
 * res(1600) resonances). The inclusive cross section is defined as
 * dsigma/dlog(E) and is returned in units of mb.
 *
 * @param particle Secondary particle id number (see cparamlib.h).
 * @param E        Secondary particle energy in GeV.
 * @param Tp       Proton kinetic energy in GeV.
 * @param params   Pointer to a ::PARAMSET struct. The struct should be
 *                 initialized before being passed to this function.
 * @return         Inclusive cross section dsigma/dlogE in mb.
 *
 */
double Cparamlib::sigma_incl_tot(const int &particle,const double &E,const double &Tp, PARAMSET* params){
    double s_tot;
    double s_nd, s_diff, s_delta, s_res;
    
    /* check whether params is a null pointer or not */
    if (params == NULL)
        return 0;
    
    s_tot = s_nd = s_diff = s_delta = s_res = 0;
    
    /* calculate sigmas */
    s_nd = sigma_incl_nd(particle, E, Tp, params);
    s_diff = sigma_incl_diff(particle, E, Tp, params);
    s_delta = sigma_incl_delta(particle, E, Tp, params);
    s_res = sigma_incl_res(particle, E, Tp, params);
    
    s_tot = s_nd + s_diff + s_delta + s_res;
    
    return s_tot;
}

/**
 * Calculates the pT distribution for non-resonance (non-diffractive plus
 * diffraction dissociation). The pT distribution is defined as the differential
 * cross section d^2sigma/dlog(E)dpT and is returned in units of mb/(GeV/c).
 *
 * @param particle  Secondary particle id number (see cparamlib.h).
 * @param pT        Secondary particle transverse momentum in GeV/c.
 * @param E         Secondary particle energy in GeV.
 * @param Tp        Proton kinetic energy in GeV.
 * @param pt_params Pointer to a ::PARAMSET_PT struct. The struct should be
 *                  initialized before being passed to this function.
 * @return          Inclusive cross section d^2sigma/dlog(E)dpT in mb/(GeV/c).
 *
 */
double Cparamlib::sigma_pt_nr(const int &particle,const double &pT,const double &E,const double &Tp, PARAMSET_PT* pt_params){
    double x;
    double x1, Lp, W;
    double sigma;
    double cutoff;
    
    /* check whether params is a null pointer or not */
    if (pt_params == NULL)
        return 0;
    
    sigma = 0;
    W = 75.0;
    Lp = 0;
    
    x = log10(E);
    
    if (particle == ID_GAMMA) {
        sigma = pt_params->a[0]*pT*exp(-pT/pt_params->a[1]);
        
        /* cutoff is the kinematic limit function as in the paper */
        if (x > 0.5) {
            Lp = 2.5;
        } else if (x > -1.0) {
            x1 = x + 1.0;
            Lp = -0.793 + exp(0.271*x1 + 0.363*x1*x1);
        } else {
            Lp = 0.09755 + 0.670*exp(1.81*x);
        }
        cutoff = 1/(1 + exp(W*(pT - Lp)));
        sigma = sigma*cutoff;
        
        if (sigma < 0)
            sigma = 0;
    }
    
    return sigma;
}

/**
 * Calculates the pT distribution for the Delta(1232) resonance. The pT
 * distribution is defined as the differential cross section d^2sigma/dlog(E)dpT
 * and is returned in units of mb/(GeV/c).
 *
 * @param particle  Secondary particle id number (see cparamlib.h).
 * @param pT        Secondary particle transverse momentum in GeV/c.
 * @param E         Secondary particle energy in GeV.
 * @param Tp        Proton kinetic energy in GeV.
 * @param pt_params Pointer to a ::PARAMSET_PT struct. The struct should be
 *                  initialized before being passed to this function.
 * @return          Inclusive cross section d^2sigma/dlog(E)dpT in mb/(GeV/c).
 *
 */
double Cparamlib::sigma_pt_delta(const int &particle,const double &pT,const double &E,const double &Tp, PARAMSET_PT* pt_params){
    double x;
    double pow;
    double x1, Lp, W;
    double sigma;
    double cutoff;
    
    /* check whether params is a null pointer or not */
    if (pt_params == NULL)
        return 0;
    
    sigma = 0;
    W = 75.0;
    Lp = 0;
    
    x = log10(E);
    
    if (particle == ID_GAMMA) {
        pow = pT - pt_params->b[0][1];
        sigma = pt_params->b[0][0]*pT*exp(-pow*pow/pt_params->b[0][2]);
        
        /* cutoff is the kinematic limit function as in the paper */
        if (x > 0.5) {
            Lp = 2.5;
        } else if (x > -1.0) {
            x1 = x + 1.0;
            Lp = -0.793 + exp(0.271*x1 + 0.363*x1*x1);
        } else {
            Lp = 0.09755 + 0.670*exp(1.81*x);
        }
        cutoff = 1/(1 + exp(W*(pT - Lp)));
        
        sigma = sigma*cutoff;
        
        if (sigma < 0)
            sigma = 0;
    }
    
    return sigma;
}

/**
 * Calculates the pT distribution for the res(1600) resonance. The pT
 * distribution is defined as the differential cross section d^2sigma/dlog(E)dpT
 * and is returned in units of mb/(GeV/c).
 *
 * @param particle  Secondary particle id number (see cparamlib.h).
 * @param pT        Secondary particle transverse momentum in GeV/c.
 * @param E         Secondary particle energy in GeV.
 * @param Tp        Proton kinetic energy in GeV.
 * @param pt_params Pointer to a ::PARAMSET_PT struct. The struct should be
 *                  initialized before being passed to this function.
 * @return          Inclusive cross section d^2sigma/dlog(E)dpT in mb/(GeV/c).
 *
 */
double Cparamlib::sigma_pt_res(const int &particle,const double &pT,const double &E,const double &Tp, PARAMSET_PT* pt_params){
    double x;
    double pow;
    double x1, Lp, W;
    double sigma;
    double cutoff;
    
    /* check whether params is a null pointer or not */
    if (pt_params == NULL)
        return 0;
    
    sigma = 0;
    W = 75.0;
    Lp = 0;
    
    x = log10(E);
    
    if (particle == ID_GAMMA) {
        pow = pT - pt_params->c[0][1];
        sigma = pt_params->c[0][0]*pT*exp(-pow*pow/pt_params->c[0][2]);
        
        /* cutoff is the kinematic limit function as in the paper */
        if (x > 0.5) {
            Lp = 2.5;
        } else if (x > -1.0) {
            x1 = x + 1.0;
            Lp = -0.793 + exp(0.271*x1 + 0.363*x1*x1);
        } else {
            Lp = 0.09755 + 0.670*exp(1.81*x);
        }
        cutoff = 1/(1 + exp(W*(pT - Lp)));
        sigma = sigma*cutoff;
        
        if (sigma < 0)
            sigma = 0;
    }
    
    return sigma;
}

/**
 * Calculates the total pT distribution, i.e. sum of all components,
 * non-resonance, Delta(1232) and res(1600) resonances. The pT distribution is
 * defined as the differential cross section d^2sigma/dlog(E)dpT and is
 * returned in units of mb/(GeV/c).
 *
 * @param particle  Secondary particle id number (see cparamlib.h).
 * @param pT        Secondary particle transverse momentum in GeV/c.
 * @param E         Secondary particle energy in GeV.
 * @param Tp        Proton kinetic energy in GeV.
 * @param pt_params Pointer to a ::PARAMSET_PT struct. The struct should be
 *                  initialized before being passed to this function.
 * @return          Inclusive cross section d^2sigma/dlog(E)dpT in mb/(GeV/c).
 *
 */
double Cparamlib::sigma_pt_tot(const int &particle,const double &pT,const double &E,const double &Tp, PARAMSET_PT* pt_params){
    double s_tot;
    double s_nr, s_delta, s_res;
    
    /* check whether params is a null pointer or not */
    if (pt_params == NULL)
        return 0;
    
    s_tot = s_nr = s_delta = s_res = 0;
    
    /* calculate sigmas */
    s_nr = sigma_pt_nr(particle, pT, E, Tp, pt_params);
    s_delta = sigma_pt_delta(particle, pT, E, Tp, pt_params);
    s_res = sigma_pt_res(particle, pT, E, Tp, pt_params);
    
    s_tot = s_nr + s_delta + s_res;
    
    return s_tot;
}

/**
 * Calculates the inelastic proton-proton cross section for the non-diffractive
 * interaction component, as given by equation 3 in Kamae et al. (2006). The
 * cross section sigma(pp) is returned in units of mb.
 *
 * @param Pp       Proton momentum in GeV/c.
 * @return         Inelastic cross section sigma(pp) in mb.
 *
 */
double Cparamlib::sigma_pp_nd(const double &Pp){
    const double a[9] = {0.57, 0.1176, 0.3829, 23.10, 6.454, -5.764, -23.63, 94.75, 0.02667};
    const double b[2] = {11.34, 23.72};
    const double c[3] = {28.5, -6.133, 1.464};
    double x, sigma;
    
    x = log10(Pp);
    sigma = 0;
    
    if ((Pp >= 1) && (Pp < 1.3)) {
        sigma = a[0]*pow(x/a[1], 1.2)*(a[3] + x*x*(a[4] + a[5]*x) + a[6]*exp(-a[7]*(x + a[8])*(x + a[8])));
    } else if ((Pp >= 1.3) && (Pp < 2.4)) {
        sigma = (b[0]*fabs(a[2] - x) + b[1]*fabs(a[1] - x))/(a[2] - a[1]);
    } else if ((Pp >= 2.4) && (Pp <= 10.0)) {
        sigma = a[3] + x*x*(a[4] + a[5]*x) + a[6]*exp(-a[7]*(x + a[8])*(x + a[8]));
    } else if (Pp > 10.0) {
        sigma = c[0] + c[1]*x + c[2]*x*x;
    }
    
    return sigma;
}

/**
 * Calculates the inelastic proton-proton cross section for the diffraction
 * dissociation component, as given by equation 3 in Kamae et al. (2006). The
 * cross section sigma(pp) is returned in units of mb.
 *
 * @param Pp       Proton momentum in GeV/c.
 * @return         Inelastic cross section sigma(pp) in mb.
 *
 */
double Cparamlib::sigma_pp_diff(const double &Pp){
    const double d[7] = {0.3522, 0.1530, 1.498, 2.0, 30.0, 3.155, 1.042};
    const double e[2] = {5.922, 1.632};
    
    double x, sigma;
    
    x = log10(Pp);
    sigma = 0;
    
    if ((Pp >= 2.25) && (Pp < 3.2)) {
        sigma = sqrt((x-d[0])/d[1])*(d[2] + d[3]*log10(d[4]*(x - 0.25)) + x*x*(d[5] - d[6]*x));
    } else if ((Pp >= 3.2) && (Pp <= 100.0)) {
        sigma = d[2] + d[3]*log10(d[4]*(x - 0.25)) + x*x*(d[5] - d[6]*x);
    } else if (Pp > 100.0) {
        sigma = e[0] + e[1]*x;
    }
    
    return sigma;
}

/**
 * Calculates the inelastic proton-proton cross section for the Delta(1232)
 * component, as given by equation 3 in Kamae et al. (2006). The cross section
 * sigma(pp) is returned in units of mb.
 *
 * @param Pp       Proton momentum in GeV/c.
 * @return         Inelastic cross section sigma(pp) in mb.
 *
 */
double Cparamlib::sigma_pp_delta(const double &Pp){
    const double f[5] = {0.0834, 9.5, -5.5, 1.68, 3134.0};
    
    double Ep, sigma;
    
    /* proton mass 0.938 GeV/c^2 gives mass squared  0.879844 (GeV/c^2)^2 */
    Ep = sqrt(Pp*Pp + 0.879844);
    sigma = 0;
    
    if ((Ep >= 1.4) && (Ep < 1.6)) {
        sigma = f[0]*pow(Ep, 10);
    } else if ((Ep >= 1.6) && (Ep < 1.8)) {
        sigma = f[1]*exp(-f[2]*(Ep - f[3])*(Ep - f[3]));
    } else if ((Ep >= 1.8) && (Ep <= 10.0)) {
        sigma = f[4]*pow(Ep, -10);
    }
    
    return sigma;
}

/**
 * Calculates the inelastic proton-proton cross section for the res(1600)
 * component, as given by equation 4 in Kamae et al. (2006). The cross section
 * sigma(pp) is returned in units of mb.
 *
 * @param Pp       Proton momentum in GeV/c.
 * @return         Inelastic cross section sigma(pp) in mb.
 *
 */
double Cparamlib::sigma_pp_res(const double &Pp){
    const double g[5] = {0.0004257, 4.5, -7.0, 2.1, 503.5};
    
    double Ep, sigma;
    
    /* proton mass 0.938 GeV/c^2 gives mass squared  0.879844 (GeV/c^2)^2 */
    Ep = sqrt(Pp*Pp + 0.879844);
    sigma = 0;
    
    if ((Ep >= 1.6) && (Ep < 1.9)) {
        sigma = g[0]*pow(Ep, 14);
    } else if ((Ep >= 1.9) && (Ep < 2.3)) {
        sigma = g[1]*exp(-g[2]*(Ep - g[3])*(Ep - g[3]));
    } else if ((Ep >= 2.3) && (Ep <= 20.0)) {
        sigma = g[4]*pow(Ep, -6);
    }
    
    return sigma;
}
