///
/// regular GMF generators
///
#ifndef HAMMURABI_BREG_H
#define HAMMURABI_BREG_H

#include <fftw3.h>
#include <vec3.h>
#include <vector>
#include "pond.h"
#include "grid.h"

///
/// base class of GMF generator,
/// \p read_grid and \p write_grid are implemented here
///
class Breg{
public:
    Breg(void) = default;
    virtual ~Breg(void) = default;
    ///
    /// fetch regular magnetic field,
    /// inovke \p read_grid regardless of field type if permitted, otherwise invoke \p breg
    ///
    virtual vec3_t<double> get_breg(const vec3_t<double> &,Pond *,Grid_breg *);
    ///
    /// field assembler, specified only in derived class
    ///
    virtual vec3_t<double> breg(const vec3_t<double> &,Pond *);
    ///
    /// read from field grid with trilinear interpolation
    ///
    virtual vec3_t<double> read_grid(const vec3_t<double> &,Grid_breg *);
    ///
    /// write to field grid
    ///
    virtual void write_grid(Pond *,Grid_breg *);
};

///
/// designed for verification
///
class Breg_verify final : public Breg{
public:
    Breg_verify(void) = default;
    virtual ~Breg_verify(void) = default;
    vec3_t<double> breg(const vec3_t<double> &,Pond *) override;
};

///
/// WMAP-3yr field modeling
///
class Breg_wmap final : public Breg {
public:
    Breg_wmap(void) = default;
    virtual ~Breg_wmap(void) = default;
    vec3_t<double> breg(const vec3_t<double> &,Pond *) override;
};

#endif

// END
