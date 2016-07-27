#ifndef MOVING_RECT_2D_H
#define MOVING_RECT_2D_H

#include "core/globalDefs.h"
#include "core/array.h"
//#include "atomicBlock/blockLattice2D.h"
#include "particles/particle2D.h"
#include "newBulkParticles/movingObject2D.h"

#include <vector>

namespace plb {

template<typename T, template<typename U> class Descriptor> class BulkParticleField2D;

template<typename T, template<typename U> class Descriptor>
class movingRect2D: public movingObject2D<T,Descriptor> {
public:
    movingRect2D();
    movingRect2D(plint tag_, T x0_, T x1_, T y0_, T y1_, Array<T,2> const& velocity_);
    ~movingRect2D(){}

    virtual void reset(T x0_, T x1_, T y0_, T y1_);
    virtual void advance();
    virtual movingObject2D<T,Descriptor>* clone() const;
    virtual void rescale(int dxScale, int dtScale);
public:
    virtual void computeBoundaryLinks(Box2D domain);

protected:
    T x0, x1, y0, y1;
};

}  // namespace plb

#endif  // MOVING_RECT_2D_H
