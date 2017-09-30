/*
 *@file: class_breg.h
 *@author: Jiaxin Wang
 *@email: jiwang@sissa.it
 *@brief: regular magnetic field generator
 */
#ifndef HAMMURABI_BREG_H
#define HAMMURABI_BREG_H

#include <fftw3.h>
#include <vec3.h>
#include <vector>

#include "class_pond.h"
#include "class_grid.h"

/* base class */
class Breg{
    public:
    Breg(void) = default;
    virtual ~Breg(void) = default;
    /* 
     *@get_breg
     * get regular magnetic field at gc position
     * read_grid regardless of field type if permitted
     * or calculate directly
     */
    virtual vec3 get_breg(const vec3 &,Pond *,Grid_breg *);
    /*@breg
     * regular field builder specified in derived class
     */
    virtual vec3 breg(const vec3 &,Pond *);
    /*@read_grid
     * read from grid with trilinear interpolation
     */
    virtual vec3 read_grid(const vec3 &, Grid_breg *);
    /*@write_grid
     * write to grid
     */
    virtual void write_grid(Pond *, Grid_breg *);
};

/* WMAP-3yr */
class Bwmap final : public Breg {
    public:
    Bwmap(void) = default;
    virtual ~Bwmap(void) = default;
    /*@Bwmap(not void)
     * reassign breg parameters in pond
     * with given vector of values
     */
    //Bwmap(const std::vector<double>&, Pond *);
    vec3 breg(const vec3 &, Pond *) override;
};

/* Local */
class Blocal final : public Breg {
    public:
    Blocal(void) = default;
    virtual ~Blocal(void) = default;
    //Blocal(const std::vector<double>&, Pond *);
    vec3 breg(const vec3 &,Pond *) override;
};

#endif
// END
