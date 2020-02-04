// reforged HEALPix pointing class
// spherical frame described by (theta,phi) in radian
// associated default Cartesian frame described by (x,y,z)
// and the default relation is
// theta = 0, z direction
// phi = 0, x direction
// phi = 90, y direction

#ifndef HAMMURABI_POINTING_H
#define HAMMURABI_POINTING_H

#include <cgsunits.h>
#include <cmath>
#include <hamvec.h>

class hampixp {
public:
#ifdef NDEBUG
protected:
#endif
  double Theta{0};
  double Phi{0};
  // normalize the angles so that 0<=theta<=pi and 0<=phi<2*pi
  void norm();

public:
  hampixp() = default;
  // copy
  hampixp(const hampixp &p) : Theta(p.Theta), Phi(p.Phi) {}
  // copy
  hampixp &operator=(const hampixp &p) noexcept {
    Theta = p.theta();
    Phi = p.phi();
    return *this;
  }
  // explicit
  hampixp(const double &t, const double &p) : Theta(t), Phi(p) { norm(); }
  // from Caretesian vector (not necessarily a versor)
  hampixp(const hamvec<3, double> &vec) {
    const double x{vec[0]};
    const double y{vec[1]};
    const double z{vec[2]};
    Theta = std::fmod(std::atan2(sqrt(x * x + y * y), z), cgs::pi);
    Theta += (Theta < 0) * cgs::pi;
    Phi = std::fmod(std::atan2(y, x), 2. * cgs::pi);
    Phi += (Phi < 0) * 2. * cgs::pi;
  }
  ~hampixp() = default;

  // read theta value
  void theta(const double &t) {
    Theta = std::fmod(t, cgs::pi);
    Theta += (Theta < 0) * cgs::pi;
  }
  // read phi value
  void phi(const double &p) {
    Phi = std::fmod(p, 2. * cgs::pi);
    Phi += (Phi < 0) * 2. * cgs::pi;
  }
  // return theta value
  double theta() const { return Theta; }
  // return phi value
  double phi() const { return Phi; }
  // return to Cartesian versor
  hamvec<3, double> versor() const;
};

void hampixp::norm() {
  Theta = std::fmod(Theta, cgs::pi);
  Theta += (Theta < 0) * cgs::pi;
  Phi = std::fmod(Phi, 2. * cgs::pi);
  Phi += (Phi < 0) * 2. * cgs::pi;
}

hamvec<3, double> hampixp::versor() const {
  const double sth{std::sin(Theta)};
  return hamvec<3, double>{sth * std::cos(Phi), sth * std::sin(Phi),
                           std::cos(Theta)};
}

#endif
