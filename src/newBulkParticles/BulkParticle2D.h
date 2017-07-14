#ifndef BULK_PARTICLE_2D_H
#define BULK_PARTICLE_2D_H

#include "core/globalDefs.h"
#include "core/array.h"
//#include "atomicBlock/blockLattice2D.h"
#include "particles/particle2D.h"

#include <vector>

namespace plb {
template<typename T, template<typename U> class Descriptor> class solidBlockLattice2D;
template<typename T, template<typename U> class Descriptor> class BulkParticleField2D;
template<typename T, template<typename U> class Descriptor> class BulkParticleChunk2D;

template<typename T, template<typename U> class Descriptor>
class BulkParticle2D : public Particle2D<T,Descriptor>{
public:
    BulkParticle2D();
    /// The tag does not need to be unique and exists only for the
    ///   user's convenience.
    BulkParticle2D(plint id_, plint tag_, Array<T,2> const& position_, T AngularPosition_,
                                          Array<T,2> const& velocity_, T AngularVelocity_, T radius_, T density_);

    virtual ~BulkParticle2D() { }
    virtual void velocityToParticle(TensorField2D<T,2>& velocityField, T scaling=1.);
    virtual void rhoBarJtoParticle(NTensorField2D<T>& rhoBarJfield, bool velIsJ, T scaling=1.);
    virtual void fluidToParticle(solidBlockLattice2D<T,Descriptor>& fluid, T scaling=1.);
    virtual void particleToFluid(BlockLattice2D<T,Descriptor>& fluid);
    virtual void particleToFluidSlip(BlockLattice2D<T,Descriptor>& fluid, T zeta);
    virtual void particleToFluidSlip2(BlockLattice2D<T,Descriptor>& fluid, T zeta);
    virtual void particleToParticleMagnetic(BulkParticleField2D<T,Descriptor>& bulkParticleField_);
    virtual void particleToParticle(BulkParticleField2D<T,Descriptor>& bulkParticleField_);
    virtual void particleToChunk(BulkParticleChunk2D<T,Descriptor>& bulkParticleChunk_);
    virtual void wallToParticle(BlockLattice2D<T,Descriptor>& fluid);
    virtual void wallToParticle(solidBlockLattice2D<T,Descriptor>& fluid);
    virtual void advance();
    Array<T,2> const& getPosition() const { return this->position; }
    Array<T,2>& getPosition() { return this->position; }
    Array<T,2>& getVelocity() { return velocity; }
    void setVelocity(Array<T,2> velocity_){ velocity = velocity_; }
    T getAngularVelocity() { return AngularVelocity; }
    T getRadius() {return radius; }
    T getDensity() {return density; }
    virtual void reset(Array<T,2> const& position_, T AngularPosition_);
    virtual void serialize(HierarchicSerializer& serializer) const;
    virtual void unserialize(HierarchicUnserializer& unserializer);
    /// Return a unique ID for this class.
    virtual int getId() const {return id; }
    /// Return this particle's tag (not necessarily unique).
    virtual plint getTag() const {return this->tag; }
    virtual BulkParticle2D<T,Descriptor>* clone() const;
//    virtual bool getScalar(plint whichScalar, T& scalar) const;
//    virtual bool getVector(plint whichVector, Array<T,2>& vector) const;
//    virtual bool getTensor(plint whichVector, Array<T,SymmetricTensorImpl<T,2>::n>& tensor) const;
//    virtual bool setScalars(std::vector<T> const& scalars);
//    virtual bool setVectors(std::vector<Array<T,2> > const& vectors);
//    virtual bool setTensors(std::vector<Array<T,SymmetricTensorImpl<T,2>::n> > const& tensors);
    virtual void rescale(int dxScale, int dtScale);
public:
       //Compute fluid boundary, solid boundary and internal solid points
    void computeBoundaries(Box2D domain);
    void computeBoundaryLinks(Box2D domain);
    bool inside(T x, T y){
        if((x-this->position[0])*(x-this->position[0]) + (y-this->position[1])*(y-this->position[1]) <= radius*radius )
            return 1;
        else return 0;
    }


protected:
    Array<T,2> velocity;
    T AngularVelocity;
    //Array<T,2> position; // inherited from base class
    T AngularPosition;
    T radius;
    T density;
    int id;
    std::vector<Array<int,2> > fluidBoundary;
    std::vector<Array<int,2> > solidBoundary;
    std::vector<Array<int,2> > internalSolid;
public:
    std::vector<Array<int,4> > LeftRightLink;
    std::vector<Array<int,4> > UpDownLink;
    std::vector<Array<int,4> > DiagonalLink;
    std::vector<Array<int,4> > antiDiagonalLink;


    template<typename T_, template<typename U_> class Descriptor_> friend class solidBlockLattice2D;
};

/// Serialize the whole particle chain into a byte-stream, and append it to
///   the vector data.
template<typename T, template<typename U> class Descriptor>
void serialize(Particle2D<T,Descriptor> const& particle, std::vector<char>& data);

/// Unserialize all data into newly generated particle objects.
template<typename T, template<typename U> class Descriptor>
void generateAndUnserializeParticles (
        std::vector<char> const& data,
        std::vector<Particle2D<T,Descriptor>*>& particle );


}  // namespace plb

#endif  // BULK_PARTICLE_2D_H
