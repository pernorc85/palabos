#ifndef MOVING_LUMEN_2D_H
#define MOVING_LUMEN_2D_H

#include "core/globalDefs.h"
#include "core/array.h"
//#include "atomicBlock/blockLattice2D.h"
#include "particles/particle2D.h"
#include "newBulkParticles/movingObject2D.h"

#include <vector>

namespace plb {

template<typename T, template<typename U> class Descriptor> class BulkParticleField2D;

template<typename T, template<typename U> class Descriptor>
class movingLumen2D: public movingObject2D<T,Descriptor> {
public:
    movingLumen2D();
    movingLumen2D(plint tag_, T x0, T x1, T y0, T y1, Array<T,2> velocity_);
    virtual ~movingLumen2D() { }

    virtual void reset(T x0_, T x1_, T y0_, T y1_);
    virtual void advance();
    virtual movingObject2D<T,Descriptor>* clone() const;
    virtual void rescale(int dxScale, int dtScale);
public:
    void computeBoundaryLinks(Box2D domain);
    bool inside(T x, T y){}
    void objectToParticle(BulkParticle2D<T,Descriptor>* particle);

protected:
    T x0, x1, y0, y1;
    T h;//deform_at_center
    //b = y0+h
    //a = -4*h/(x0-x1)^2
    //y = a(x-(x0+x1)/2)^2 + b

};

}  // namespace plb

#endif  // MOVING_OBJECT_2D_H
