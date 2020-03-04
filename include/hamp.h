// reforged HEALPix pointing class
// spherical frame described by (theta,phi) in radian
// associated default Cartesian frame described by (x,y,z)
// and the default relation is
// theta = 0, z direction
// phi = 0, x direction
// phi = 90, y direction

#ifndef HAMMURABI_POINTING_H
#define HAMMURABI_POINTING_H

#include <cmath>

#include <hamtype.h>
#include <hamunits.h>
#include <hamvec.h>

class Hamp {
protected:
  ham_float Theta{0};
  ham_float Phi{0};
  // normalize the angles so that 0<=theta<=pi and 0<=phi<2*pi
  void norm() {
    this->Theta = std::fmod(this->Theta, cgs::pi);
    this->Theta += (this->Theta < 0) * cgs::pi;
    this->Phi = std::fmod(this->Phi, cgs::twopi);
    this->Phi += (this->Phi < 0) * cgs::twopi;
  }

public:
  Hamp() = default;
  // copy constr
  Hamp(const Hamp &p) noexcept {
    this->Theta = p.theta();
    this->Phi = p.phi();
  }
  // copy assign
  Hamp &operator=(const Hamp &p) noexcept {
    this->Theta = p.theta();
    this->Phi = p.phi();
    return *this;
  }
  // move constr
  Hamp(Hamp &&p) noexcept {
    this->Theta = std::move(p.Theta);
    this->Phi = std::move(p.Phi);
  }
  // move assign
  Hamp &operator=(Hamp &&p) noexcept {
    this->Theta = std::move(p.Theta);
    this->Phi = std::move(p.Phi);
    return *this;
  }
  // explicit
  Hamp(const ham_float &t, const ham_float &p) {
    this->Theta = t;
    this->Phi = p;
    norm();
  }
  // from Caretesian vector (not necessarily a versor)
  Hamp(const Hamvec<3, ham_float> &vec) {
    const ham_float x{vec[0]};
    const ham_float y{vec[1]};
    const ham_float z{vec[2]};
    this->Theta = std::fmod(std::atan2(sqrt(x * x + y * y), z), cgs::pi);
    this->Theta += (this->Theta < 0) * cgs::pi;
    this->Phi = std::fmod(std::atan2(y, x), cgs::twopi);
    this->Phi += (this->Phi < 0) * cgs::twopi;
  }
  ~Hamp() = default;

  // read in theta value
  void theta(const ham_float &t) {
    this->Theta = std::fmod(t, cgs::pi);
    this->Theta += (this->Theta < 0) * cgs::pi;
  }
  // read in phi value
  void phi(const ham_float &p) {
    this->Phi = std::fmod(p, cgs::twopi);
    this->Phi += (this->Phi < 0) * cgs::twopi;
  }
  // return theta value
  ham_float theta() const { return this->Theta; }
  // return phi value
  ham_float phi() const { return this->Phi; }
  // return to Cartesian versor
  Hamvec<3, ham_float> versor() const {
    const ham_float sth{std::sin(this->Theta)};
    return Hamvec<3, ham_float>{sth * std::cos(this->Phi),
                                sth * std::sin(this->Phi),
                                std::cos(this->Theta)};
  }
};

#endif
