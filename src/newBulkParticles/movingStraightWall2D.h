#ifndef MOVING_STRAIGHT_WALL_2D_H
#define MOVING_STRAIGHT_WALL_2D_H

#include "core/globalDefs.h"
#include "core/array.h"
//#include "atomicBlock/blockLattice2D.h"
#include "particles/particle2D.h"
#include "newBulkParticles/movingObject2D.h"

#include <vector>

namespace plb {

template<typename T, template<typename U> class Descriptor>
class movingVerticalWall2D: public movingObject2D<T,Descriptor> {
public:
    movingVerticalWall2D();
    movingVerticalWall2D(plint tag_, T x0_, T y0_, T y1_, Array<T,2> velocity_);
    virtual ~movingVerticalWall2D() { }

    virtual void reset(T x0_, T y0_, T y1_);
    virtual T get_x0(){return x0;}
    virtual T get_y0(){return y0;}
    virtual T get_y1(){return y1;}
    virtual void advance();
    virtual movingObject2D<T,Descriptor>* clone() const;
    virtual void rescale(int dxScale, int dtScale);
public:
    void computeBoundaryLinks(Box2D domain);
    bool inside(T x, T y){}
    void objectToParticle(BulkParticle2D<T,Descriptor>* particle){};

protected:
    T x0, y0, y1;
};

template<typename T, template<typename U> class Descriptor>
class movingInclinedWall2D: public movingObject2D<T,Descriptor> {
public:
    movingInclinedWall2D();
    movingInclinedWall2D(plint tag_, T x0_, T y0_, T y1_, T k_, Array<T,2> velocity_);
    virtual ~movingInclinedWall2D() { }

    virtual void reset(T x0_, T y0_, T y1_, T k);
    virtual T get_x0(){return x0;}
    virtual T get_y0(){return y0;}
    virtual T get_y1(){return y1;}
    virtual T get_k(){return k;}
    virtual void advance();
    virtual movingObject2D<T,Descriptor>* clone() const;
    virtual void rescale(int dxScale, int dtScale);
public:
    void computeBoundaryLinks(Box2D domain);
    bool inside(T x, T y){}
    void objectToParticle(BulkParticle2D<T,Descriptor>* particle){};

protected:
    T x0, y0, y1, k;
};

}  // namespace plb

#endif  // MOVING_VERTICAL_WALL_2D_H
