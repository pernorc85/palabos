#ifndef MOVING_OBJECT_2D_H
#define MOVING_OBJECT_2D_H

#include "core/globalDefs.h"
#include "core/array.h"
//#include "atomicBlock/blockLattice2D.h"
#include "particles/particle2D.h"

#include <vector>

namespace plb {

template<typename T, template<typename U> class Descriptor> class BulkParticle2D;

template<typename T, template<typename U> class Descriptor>
class movingObject2D {
public:
    movingObject2D();
    movingObject2D(plint tag_, Array<T,2> const& velocity_);
    virtual ~movingObject2D(){};

    void objectToFluid(BlockLattice2D<T,Descriptor>& fluid);
    void objectToFluidSlip(BlockLattice2D<T,Descriptor>& fluid,T zeta);
    virtual void objectToParticle(BulkParticle2D<T,Descriptor>* particle)=0;
    virtual void advance()=0;
    Array<T,2>& getVelocity() { return velocity; }
    virtual void setVelocity(Array<T,2> velocity_){velocity = velocity_;}
    int getId() const {return id; }
    plint getTag() const {return tag; }
    virtual movingObject2D<T,Descriptor>* clone() const=0;
//    virtual bool getScalar(plint whichScalar, T& scalar) const;
//    virtual bool getVector(plint whichVector, Array<T,2>& vector) const;
//    virtual bool getTensor(plint whichVector, Array<T,SymmetricTensorImpl<T,2>::n>& tensor) const;
//    virtual bool setScalars(std::vector<T> const& scalars);
//    virtual bool setVectors(std::vector<Array<T,2> > const& vectors);
//    virtual bool setTensors(std::vector<Array<T,SymmetricTensorImpl<T,2>::n> > const& tensors);
    virtual void rescale(int dxScale, int dtScale)=0;
public:
    virtual void computeBoundaryLinks(Box2D domain)=0;
    virtual bool inside(T x, T y)=0;


protected:
    Array<T,2> velocity;
    int tag;
    int id;
    std::vector<Array<int,4> > LeftRightLink;
    std::vector<Array<int,4> > UpDownLink;
    std::vector<Array<int,4> > DiagonalLink;
    std::vector<Array<int,4> > antiDiagonalLink;
};

}  // namespace plb

#endif  // MOVING_OBJECT_2D_H
