#ifndef BULK_PARTICLE_3D_H
#define BULK_PARTICLE_3D_H

#include "core/globalDefs.h"
#include "core/array.h"
//#include "atomicBlock/blockLattice2D.h"
#include "particles/particle2D.h"

#include <vector>

namespace plb {

template<typename T, template<typename U> class Descriptor> class solidBlockLattice3D;
template<typename T, template<typename U> class Descriptor> class BulkParticleField3D;

template<typename T, template<typename U> class Descriptor>
class BulkParticle3D : public Particle3D<T,Descriptor>{
public:
    BulkParticle3D();
    /// The tag does not need to be unique and exists only for the
    ///   user's convenience.
    BulkParticle3D(plint tag_, Array<T,3> const& position_, T AngularPosition_,
                               Array<T,3> const& velocity_, T AngularVelocity_, T radius_, T density_);

    virtual ~BulkParticle3D() { }
    virtual void velocityToParticle(TensorField3D<T,3>& velocityField, T scaling=1.);
    virtual void rhoBarJtoParticle(NTensorField3D<T>& rhoBarJfield, bool velIsJ, T scaling=1.);
    virtual void fluidToParticle(solidBlockLattice3D<T,Descriptor>& fluid, T scaling=1.);
    virtual void particleToFluid(BlockLattice3D<T,Descriptor>& fluid);
    virtual void particleToParticle(BulkParticleField3D<T,Descriptor>& bulkParticleField_);
    virtual void wallToParticle(solidBlockLattice3D<T,Descriptor>& fluid);
    virtual void advance();
    Array<T,3> const& getPosition() const { return this->position; }
    Array<T,3>& getPosition() { return this->position; }
    Array<T,3>& getVelocity() { return velocity; }
    T getAngularVelocity() { return AngularVelocity; }
    T getRadius() {return radius; }
    T getDensity() {return density; }
    virtual void reset(Array<T,3> const& position_, T AngularPosition_);
    virtual void serialize(HierarchicSerializer& serializer) const;
    virtual void unserialize(HierarchicUnserializer& unserializer);
    /// Return a unique ID for this class.
    virtual int getId() const {return id; }
    /// Return this particle's tag (not necessarily unique).
    virtual plint getTag() const {return this->tag; }
    virtual BulkParticle3D<T,Descriptor>* clone() const;
//    virtual bool getScalar(plint whichScalar, T& scalar) const;
//    virtual bool getVector(plint whichVector, Array<T,2>& vector) const;
//    virtual bool getTensor(plint whichVector, Array<T,SymmetricTensorImpl<T,2>::n>& tensor) const;
//    virtual bool setScalars(std::vector<T> const& scalars);
//    virtual bool setVectors(std::vector<Array<T,2> > const& vectors);
//    virtual bool setTensors(std::vector<Array<T,SymmetricTensorImpl<T,2>::n> > const& tensors);
    virtual void rescale(int dxScale, int dtScale);
public:
       //Compute fluid boundary, solid boundary and internal solid points
    void computeBoundaryLinks(Box3D domain);
    bool inside(T x, T y, T z){
        if((x-this->position[0])*(x-this->position[0]) + (y-this->position[1])*(y-this->position[1])
                                                       + (z-this->position[2])*(z-this->position[2]) ) <= radius*radius*radius )
            return 1;
        else return 0;
    }


protected:
    Array<T,3> velocity;
    T AngularVelocity;
    //Array<T,2> position; // inherited from base class
    T AngularPosition;
    T radius;
    T density;
    int id;
    std::vector<Array<int,2> > fluidBoundary;
    std::vector<Array<int,2> > solidBoundary;
    std::vector<Array<int,2> > internalSolid;
    std::vector<Array<int,6> > LeftRightLink;
    std::vector<Array<int,6> > UpDownLink;
    std::vector<Array<int,6> > DiagonalLink;
    std::vector<Array<int,6> > antiDiagonalLink;
    
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
