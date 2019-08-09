// constants in CGS units

#define cgs_erg 1.
#define cgs_cm 1.
#define cgs_sec 1.
#define cgs_Kelvin 1.

#define cgs_gram (cgs_erg * cgs_sec * cgs_sec / (cgs_cm * cgs_cm))
#define cgs_joule (1.e7 * cgs_erg)
#define cgs_watt (cgs_joule / cgs_sec)

#define cgs_pi 3.14159265358979
#define cgs_rad (cgs_pi / 180.)
#define cgs_arcmin (2. * cgs_pi / (360. * 60.))
#define cgs_arcsec (2. * cgs_pi / (360. * 60. * 60.))
#define cgs_sterad (4. * cgs_pi)
#define cgs_esu sqrt(cgs_erg *cgs_cm)
#define cgs_Coulomb (2997924579.99957 * cgs_esu)
#define cgs_year (365.25 * 24. * 3600. * cgs_sec)
#define cgs_m (1.e2 * cgs_cm)
#define cgs_km (1.e5 * cgs_cm)
#define cgs_pc (3.0856775806e+18 * cgs_cm)
#define cgs_kpc (3.0856775806e+21 * cgs_cm)
#define cgs_Mpc (1.e3 * cgs_kpc)
#define cgs_Gpc (1.e6 * cgs_kpc)
#define cgs_eV (1.60217733e-12 * cgs_erg)
#define cgs_keV (1.e3 * cgs_eV)
#define cgs_MeV (1.e6 * cgs_eV)
#define cgs_GeV (1.e9 * cgs_eV)
#define cgs_Hz (1. / cgs_sec)
#define cgs_kHz (1.e3 * cgs_Hz)
#define cgs_MHz (1.e6 * cgs_Hz)
#define cgs_GHz (1.e9 * cgs_Hz)
#define cgs_barn (1.e-24 * cgs_cm * cgs_cm)
#define cgs_mbarn (1.e-27 * cgs_cm * cgs_cm)
#define cgs_Jansky (1.e-23 * cgs_erg / cgs_sec / (cgs_cm * cgs_cm) / cgs_Hz)
#define cgs_Gauss sqrt(cgs_erg / (cgs_cm * cgs_cm * cgs_cm))
#define cgs_muGauss (1.e-6 * cgs_Gauss)
#define cgs_ccm (cgs_cm * cgs_cm * cgs_cm)

#define cgs_c_light (2.99792458e+10 * cgs_cm / cgs_sec)
#define cgs_h_planck (6.626075540e-27 * cgs_erg * cgs_sec) ///< Planck constant
#define cgs_hq (1.05457266e-27 * cgs_erg * cgs_sec) ///< Plancks constant/(2pi)
#define cgs_qe (4.80320425e-10 * cgs_esu)
#define cgs_mec2 (0.51099907e-3 * cgs_GeV) ///< Electron Mass times c^2
#define cgs_mpc2 (938.272310e-3 * cgs_GeV) ///< Proton Mass times c^2
#define cgs_me (cgs_mec2 / (cgs_c_light * cgs_c_light))
#define cgs_mec (cgs_mec2 / cgs_c_light)
#define cgs_mp (cgs_mpc2 / (cgs_c_light * cgs_c_light))
#define cgs_sigmaT (0.66524616 * cgs_barn)
#define cgs_kB (1.380622e-16 * cgs_erg / cgs_Kelvin)
#define cgs_re (2.81794092e-13 * cgs_cm) ///< classical electron radius

#define cgs_GV (cgs_GeV / cgs_qe) ///< rigidity for cosmic-rays
