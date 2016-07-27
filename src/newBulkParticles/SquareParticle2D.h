#ifndef SQAURE_PARTICLE_2D_H
#define SQAURE_PARTICLE_2D_H

#include "core/globalDefs.h"
#include "core/array.h"
//#include "atomicBlock/blockLattice2D.h"
#include "particles/particle2D.h"

#include <vector>

namespace plb {

template<typename T, template<typename U> class Descriptor> class solidBlockLattice2D;
template<typename T, template<typename U> class Descriptor> class SquareParticleField2D;

template<typename T, template<typename U> class Descriptor>
class SquareParticle2D : public Particle2D<T,Descriptor>{
public:
    SquareParticle2D();
    /// The tag does not need to be unique and exists only for the
    ///   user's convenience.
    SquareParticle2D(plint tag_, Array<T,2> const& position_, T AngularPosition_,
                               Array<T,2> const& velocity_, T AngularVelocity_, T length_, T density_);

    virtual ~SquareParticle2D() { }
    virtual void velocityToParticle(TensorField2D<T,2>& velocityField, T scaling=1.);
    virtual void rhoBarJtoParticle(NTensorField2D<T>& rhoBarJfield, bool velIsJ, T scaling=1.);
    virtual void fluidToParticle(solidBlockLattice2D<T,Descriptor>& fluid, T scaling=1.);
    virtual void advance();
    Array<T,2> const& getPosition() const { return this->position; }
    Array<T,2>& getPosition() { return this->position; }
    Array<T,2>& getVelocity() { return velocity; }
    T getAngularVelocity() { return AngularVelocity; }
    T getLength() {return length; }
    T getDensity() {return density; }
    virtual void reset(Array<T,2> const& position_, T AngularPosition_);
    virtual void serialize(HierarchicSerializer& serializer) const;
    virtual void unserialize(HierarchicUnserializer& unserializer);
    /// Return a unique ID for this class.
    virtual int getId() const {return id; }
    /// Return this particle's tag (not necessarily unique).
    virtual plint getTag() const {return this->tag; }
    virtual SquareParticle2D<T,Descriptor>* clone() const;
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
    //void computeSolidFraction(Array<int,2> cellPosition, double& solidFrac);
    bool line1(T x, T y);
    bool line2(T x, T y);
    bool line3(T x, T y);
    bool line4(T x, T y);
    T line1_getX(T y);  T line1_getY(T x);
    T line2_getX(T y);  T line2_getY(T x);
    T line3_getX(T y);  T line3_getY(T x);
    T line4_getX(T y);  T line4_getY(T x);
    bool inside(T x, T y){
        T theta = AngularPosition;
        T x0 = this->position[0];
        T y0 = this->position[1];
        T l = length;
        //y - (y0 + l/2*cos(theta) ) =  tan(theta)  *(x - (x0 - l/2*sin(theta) ) )
        //y - (y0 - l/2*cos(theta) ) =  tan(theta)  *(x - (x0 + l/2*sin(theta) ) )
        //y - (y0 + l/2*sin(theta) ) = -1/tan(theta)*(x - (x0 + l/2*cos(theta) ) )
        //y - (y0 - 1/2*sin(theta) ) = -1/tan(theta)*(x - (x0 - l/2*cos(theta) ) )
        if( (y - (y0 + l/2*cos(theta) ) -  tan(theta)  *(x - (x0 - l/2*sin(theta) ) ) )*
            (y - (y0 - l/2*cos(theta) ) -  tan(theta)  *(x - (x0 + l/2*sin(theta) ) ) ) <= 0
           &&
            (y - (y0 + l/2*sin(theta) ) +  1/tan(theta)*(x - (x0 + l/2*cos(theta) ) ) )*
            (y - (y0 - l/2*sin(theta) ) +  1/tan(theta)*(x - (x0 - l/2*cos(theta) ) ) ) <= 0
          )
            return 1;
        else return 0;
    }


protected:
    Array<T,2> velocity;
    T AngularVelocity;
    //Array<T,2> position; // inherited from base class
    T AngularPosition;
    T length;
    T density;
    int id;
    std::vector<Array<int,2> > fluidBoundary;
    std::vector<Array<int,2> > solidBoundary;
    std::vector<Array<int,2> > internalSolid;

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
