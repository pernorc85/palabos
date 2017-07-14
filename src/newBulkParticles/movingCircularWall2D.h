#ifndef MOVING_CIRCULAR_WALL_2D_H
#define MOVING_CIRCULAR_WALL_2D_H

#include "core/globalDefs.h"
#include "core/array.h"
//#include "atomicBlock/blockLattice2D.h"
#include "particles/particle2D.h"
#include "newBulkParticles/movingObject2D.h"

#include <vector>

namespace plb {

template<typename T, template<typename U> class Descriptor>
class movingCircularWall2D: public movingObject2D<T,Descriptor> {
public:
    movingCircularWall2D();
    movingCircularWall2D(plint tag_, T x0_, T y0_, T innerRadius_, T outerRadius_, T angularVelocity_);
    virtual ~movingCircularWall2D() { }

    virtual void reset(T x0_, T y0_, T innerRadius_, T outerRadius_);
    virtual T get_x0(){return x0;}
    virtual T get_y0(){return y0;}
    virtual T get_innerRadius(){return innerRadius;}
    virtual T get_outerRadius(){return outerRadius;}
    virtual void advance();
    virtual movingObject2D<T,Descriptor>* clone() const;
    virtual void rescale(int dxScale, int dtScale);
public:
    void computeBoundaries(Box2D domain);
    void computeBoundaryLinks(Box2D domain);
    bool inside(T x, T y){
        T dx = x-this->x0;
        T dy = y-this->y0;
        if(dx*dx + dy*dy >= innerRadius*innerRadius && dx*dx + dy*dy <= outerRadius*outerRadius){
            return true;
        }
        else{ return false; }
    }
    void objectToParticle(BulkParticle2D<T,Descriptor>* particle){}
    void objectToFluid(BlockLattice2D<T,Descriptor>& fluid);
    void objectToFluidSlip(BlockLattice2D<T,Descriptor>& fluid, T zeta);
    void objectToFluidSlip2(BlockLattice2D<T,Descriptor>& fluid, T zeta);

protected:
    T x0, y0, innerRadius, outerRadius, angularVelocity;
    std::vector<Array<int,2> > fluidBoundary;
    std::vector<Array<int,2> > solidBoundary;
};

}  // namespace plb

#endif  // MOVING_VERTICAL_WALL_2D_H

