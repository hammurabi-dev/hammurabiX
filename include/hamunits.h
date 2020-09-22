// constants in CGS units

#ifndef HAMMURABI_UNITS_H
#define HAMMURABI_UNITS_H

#include <cmath>

namespace cgs {

const double erg = 1.;
const double cm = 1.;
const double sec = 1.;
const double Kelvin = 1.;

const double gram = (erg * sec * sec / (cm * cm));
const double joule = (1.e7 * erg);
const double watt = (joule / sec);

const double pi = 3.141592653589793238462643383279502884197;
const double halfpi = 1.570796326794896619231321691639751442099;
const double onethird = 0.333333333333333333333333333333333333333;
const double twothirds = 0.666666666666666666666666666666666666667;
const double twopi = 6.283185307179586476925286766559005768394;
const double sqrtwo = 1.4142135623730950488016887242096980785697;
const double sqrthree = 1.7320508075688772935274463415058723669428;
const double rad = (pi / 180.);
const double arcmin = (2. * pi / (360. * 60.));
const double arcsec = (2. * pi / (360. * 60. * 60.));
const double sterad = (4. * pi);
const double esu = std::sqrt(erg * cm);
const double Coulomb = (2997924579.99957 * esu);
const double year = (365.25 * 24. * 3600. * sec);
const double m = (1.e2 * cm);
const double km = (1.e5 * cm);
const double pc = (3.0856775806e+18 * cm);
const double kpc = (3.0856775806e+21 * cm);
const double Mpc = (1.e3 * kpc);
const double Gpc = (1.e6 * kpc);
const double eV = (1.60217733e-12 * erg);
const double keV = (1.e3 * eV);
const double MeV = (1.e6 * eV);
const double GeV = (1.e9 * eV);
const double Hz = (1. / sec);
const double kHz = (1.e3 * Hz);
const double MHz = (1.e6 * Hz);
const double GHz = (1.e9 * Hz);
const double barn = (1.e-24 * cm * cm);
const double mbarn = (1.e-27 * cm * cm);
const double Jansky = (1.e-23 * erg / sec / (cm * cm) / Hz);
const double Gauss = std::sqrt(erg / (cm * cm * cm));
const double muGauss(1.e-6 * Gauss);
const double ccm(cm *cm *cm);

const double c_light = (2.99792458e+10 * cm / sec);
const double h_planck = (6.626075540e-27 * erg * sec); ///< Planck constant
const double hq = (1.05457266e-27 * erg * sec); ///< Plancks constant/(2pi)
const double qe = (4.80320425e-10 * esu);
const double mec2 = (0.51099907e-3 * GeV); ///< Electron Mass times c^2
const double mpc2 = (938.272310e-3 * GeV); ///< Proton Mass times c^2
const double me = (mec2 / (c_light * c_light));
const double mec = (mec2 / c_light);
const double mp = (mpc2 / (c_light * c_light));
const double sigmaT = (0.66524616 * barn);
const double kB = (1.380622e-16 * erg / Kelvin);
const double re = (2.81794092e-13 * cm); ///< classical electron radius

const double GV = (GeV / qe); ///< rigidity for cosmic-rays
} // namespace cgs

#endif
