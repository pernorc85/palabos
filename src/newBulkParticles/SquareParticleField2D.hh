#ifndef SQUARE_PARTICLE_FIELD_2D_HH
#define SQUARE_PARTICLE_FIELD_2D_HH

#include "core/globalDefs.h"
#include "particles/particleField2D.h"
#include "newBulkParticles/SquareParticleField2D.h"
#include "newBulkParticles/solidBlockLattice2D.h"
#include <utility>

namespace plb {

/* *************** class BulkParticleDataTransfer2D ************************ */

template<typename T, template<typename U> class Descriptor>
SquareParticleDataTransfer2D<T,Descriptor>::SquareParticleDataTransfer2D (
        SquareParticleField2D<T,Descriptor>& squareParticleField_)
    : squareParticleField(squareParticleField_)
{ }

template<typename T, template<typename U> class Descriptor>
plint SquareParticleDataTransfer2D<T,Descriptor>::staticCellSize() const {
    return 0;  // Particle containers have only dynamic data.
}

template<typename T, template<typename U> class Descriptor>
void SquareParticleDataTransfer2D<T,Descriptor>::send (
        Box2D domain, std::vector<char>& buffer, modif::ModifT kind ) const
{
    buffer.clear();
    // Particles, by definition, are dynamic data, and they need to
    //   be reconstructed in any case. Therefore, the send procedure
    //   is run whenever kind is one of the dynamic types.
    if ( (kind==modif::dynamicVariables) ||
         (kind==modif::allVariables) ||
         (kind==modif::dataStructure) )
    {
        std::vector<Particle2D<T,Descriptor>*> foundParticles;
        squareParticleField.findParticles(domain, foundParticles);
        for (pluint iParticle=0; iParticle<foundParticles.size(); ++iParticle) {
            // The serialize function automatically reallocates memory for buffer.
            serialize(*foundParticles[iParticle], buffer);
        }
    }
}

template<typename T, template<typename U> class Descriptor>
void SquareParticleDataTransfer2D<T,Descriptor>::receive (
        Box2D domain, std::vector<char> const& buffer, modif::ModifT kind )
{
    PLB_PRECONDITION(contained(domain, particleField.getBoundingBox()));
    // Clear the existing data before introducing the new data.
    squareParticleField.removeParticles(domain);
    // Particles, by definition, are dynamic data, and they need to
    //   be reconstructed in any case. Therefore, the receive procedure
    //   is run whenever kind is one of the dynamic types.
    if ( (kind==modif::dynamicVariables) ||
         (kind==modif::allVariables) ||
         (kind==modif::dataStructure) )
    {
        pluint posInBuffer = 0;
        while (posInBuffer < buffer.size()) {
            // 1. Generate dynamics object, and unserialize dynamic data.
            HierarchicUnserializer unserializer(buffer, posInBuffer);
            Particle2D<T,Descriptor>* newParticle =
                meta::particleRegistration2D<T,Descriptor>().generate(unserializer);
            posInBuffer = unserializer.getCurrentPos();
            squareParticleField.addParticle(domain, newParticle);
        }
    }
}

template<typename T, template<typename U> class Descriptor>
void SquareParticleDataTransfer2D<T,Descriptor>::receive (
        Box2D domain, std::vector<char> const& buffer, modif::ModifT kind, Dot2D absoluteOffset )
{
    if (absoluteOffset.x == 0 && absoluteOffset.y == 0) {
        receive(domain, buffer, kind);
        return;
    }
    PLB_PRECONDITION(contained(domain, particleField.getBoundingBox()));
    Array<T,2> realAbsoluteOffset((T)absoluteOffset.x, (T)absoluteOffset.y);
    // Clear the existing data before introducing the new data.
    squareParticleField.removeParticles(domain);
    // Particles, by definition, are dynamic data, and they need to
    //   be reconstructed in any case. Therefore, the receive procedure
    //   is run whenever kind is one of the dynamic types.
    if ( (kind==modif::dynamicVariables) ||
         (kind==modif::allVariables) ||
         (kind==modif::dataStructure) )
    {
        pluint posInBuffer = 0;
//         while (posInBuffer < buffer.size()) {
//             // 1. Generate dynamics object, and unserialize dynamic data.
//             HierarchicUnserializer unserializer(buffer, posInBuffer);
//             Particle2D<T,Descriptor>* newParticle =
//                 meta::particleRegistration2D<T,Descriptor>().generate(unserializer);
//             posInBuffer = unserializer.getCurrentPos();
//             newParticle -> getPosition() += realAbsoluteOffset;
//             particleField.addParticle(domain, newParticle);
//         }
        while (posInBuffer < buffer.size()) {
            // 1. Generate dynamics object, and unserialize dynamic data.
            HierarchicUnserializer unserializer(buffer, posInBuffer);
            Particle2D<T,Descriptor>* newParticle =
                meta::particleRegistration2D<T,Descriptor>().generate(unserializer);
            posInBuffer = unserializer.getCurrentPos();
            newParticle -> getPosition() += realAbsoluteOffset;
            Array<T,2> pos(newParticle->getPosition());
            Dot2D location(squareParticleField.getLocation());

            squareParticleField.addParticle(domain, newParticle);
        }
    }
}

template<typename T, template<typename U> class Descriptor>
void SquareParticleDataTransfer2D<T,Descriptor>::attribute (
        Box2D toDomain, plint deltaX, plint deltaY,
        AtomicBlock2D const& from, modif::ModifT kind )
{
    Box2D fromDomain(toDomain.shift(deltaX,deltaY));
    std::vector<char> buffer;
    SquareParticleField2D<T,Descriptor> const& fromParticleField =
        dynamic_cast<SquareParticleField2D<T,Descriptor>const &>(from);
    fromParticleField.getDataTransfer().send(fromDomain, buffer, kind);
    receive(toDomain, buffer, kind);
}

template<typename T, template<typename U> class Descriptor>
void SquareParticleDataTransfer2D<T,Descriptor>::attribute (
        Box2D toDomain, plint deltaX, plint deltaY,
        AtomicBlock2D const& from, modif::ModifT kind, Dot2D absoluteOffset )
{
    Box2D fromDomain(toDomain.shift(deltaX,deltaY));
    std::vector<char> buffer;
    SquareParticleField2D<T,Descriptor> const& fromParticleField =
        dynamic_cast<SquareParticleField2D<T,Descriptor>const &>(from);
    fromParticleField.getDataTransfer().send(fromDomain, buffer, kind);
    receive(toDomain, buffer, kind, absoluteOffset);
}

/* *************** class BulkParticleField2D ************************************ */

template<typename T, template<typename U> class Descriptor>
SquareParticleField2D<T,Descriptor>::SquareParticleField2D(plint nx, plint ny)
    : ParticleField2D<T,Descriptor>(nx,ny),
      dataTransfer(*this)
{ }

template<typename T, template<typename U> class Descriptor>
SquareParticleField2D<T,Descriptor>::~SquareParticleField2D()
{
    particles.clear();
}



template<typename T, template<typename U> class Descriptor>
SquareParticleField2D<T,Descriptor>::SquareParticleField2D(SquareParticleField2D const& rhs)
    : dataTransfer(*this)
{
    particles.resize(rhs.particles.size());
    for (pluint iParticle=0; iParticle<particles.size(); ++iParticle) {
        particles[iParticle] = rhs.particles[iParticle].clone();
    }
}


template<typename T, template<typename U> class Descriptor>
SquareParticleField2D<T,Descriptor>&
    SquareParticleField2D<T,Descriptor>::operator=(SquareParticleField2D<T,Descriptor> const& rhs)
{
    SquareParticleField2D<T,Descriptor>(rhs).swap(*this);
    return *this;
}


template<typename T, template<typename U> class Descriptor>
SquareParticleField2D<T,Descriptor>*
    SquareParticleField2D<T,Descriptor>::clone() const
{
    return new SquareParticleField2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void SquareParticleField2D<T,Descriptor>::swap(SquareParticleField2D<T,Descriptor>& rhs) {
    particles.swap(rhs.particles);
}

template<typename T, template<typename U> class Descriptor>
void SquareParticleField2D<T,Descriptor>::addParticle(Box2D domain,  SquareParticle2D<T,Descriptor>* particle) {
    plint iX, iY;
    computeGridPosition(particle->getPosition(), iX, iY);
    bool overlap = 0;

    if( contained(iX,iY, domain) )
    {
        particles.push_back(particle);
        std::cout<<"particle added."<<std::endl;
    }
    else {
        delete particle;
        std::cout<<"particle deleted."<<std::endl;
    }
}

template<typename T, template<typename U> class Descriptor>
void SquareParticleField2D<T,Descriptor>::removeParticles(Box2D domain) {
    plint iX, iY;
    typename std::vector<SquareParticle2D<T,Descriptor>*>::iterator it = particles.begin();
    for( ; it != particles.end(); ){
        computeGridPosition((*it)->getPosition(), iX, iY);
        if( contained(iX,iY, domain) ){
            delete *it;
            it = particles.erase(it);
        }
        else { ++it; }
    }
}



template<typename T, template<typename U> class Descriptor>
void SquareParticleField2D<T,Descriptor>::removeParticles(Box2D domain, plint tag) {
    plint iX, iY;
    typename std::vector<SquareParticle2D<T,Descriptor>*>::iterator it = particles.begin();
    for( ; it != particles.end(); ) {
        computeGridPosition((*it)->getPosition(), iX, iY);
        if ( contained(iX,iY, domain) && (*it)->getTag() == tag) {
            delete *it;
            it = particles.erase(it);
        }
        else { ++it; }
    }
}



template<typename T, template<typename U> class Descriptor>
void SquareParticleField2D<T,Descriptor>::findParticles (
        Box2D domain, std::vector<Particle2D<T,Descriptor>*>& found )
{
    found.clear();
    //PLB_ASSERT( contained(domain, particleGrid.getBoundingBox()) );
    plint iX, iY;
    typename std::vector<SquareParticle2D<T,Descriptor>*>::iterator it = particles.begin();
    for( ; it != particles.end(); ) {
        computeGridPosition((*it)->getPosition(), iX, iY);
        if ( contained(iX,iY, domain))  found.push_back(*it);
    }
}

template<typename T, template<typename U> class Descriptor>
void SquareParticleField2D<T,Descriptor>::findParticles (
        Box2D domain, std::vector<Particle2D<T,Descriptor> const*>& found ) const
{
    found.clear();
    //PLB_ASSERT( contained(domain, particleGrid.getBoundingBox()) );
    plint iX, iY;
    for(plint i=0 ; i < particles.size(); ++i) {
        computeGridPosition(particles[i]->getPosition(), iX, iY);
        //if ( contained(iX,iY, domain))  found.push_back(*(particles[i]));
    }
}



template<typename T, template<typename U> class Descriptor>
void SquareParticleField2D<T,Descriptor>::velocityToParticleCoupling (
        Box2D domain, TensorField2D<T,2>& velocityField, T scaling )
{
    Box2D finalDomain;
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (pluint iParticle=0; iParticle<particles.size(); ++iParticle) {
                    particles[iParticle]->velocityToParticle(velocityField, scaling);
                }
            }
        }
}



template<typename T, template<typename U> class Descriptor>
void SquareParticleField2D<T,Descriptor>::rhoBarJtoParticleCoupling (
        Box2D domain, NTensorField2D<T>& rhoBarJfield, bool velIsJ, T scaling )
{
    Box2D finalDomain;
        for (plint iX=finalDomain.x0; iX<=finalDomain.x1; ++iX) {
            for (plint iY=finalDomain.y0; iY<=finalDomain.y1; ++iY) {
                for (pluint iParticle=0; iParticle<particles.size(); ++iParticle) {
                    particles[iParticle]->rhoBarJtoParticle(rhoBarJfield, velIsJ, scaling);
                }
            }
        }
}



template<typename T, template<typename U> class Descriptor>
void SquareParticleField2D<T,Descriptor>::fluidToParticleCoupling (
        Box2D domain, solidBlockLattice2D<T,Descriptor>& lattice )
{
    plint iX, iY;
    typename std::vector<SquareParticle2D<T,Descriptor>* >::iterator it = particles.begin();
    for( ; it != particles.end(); it++) {
        computeGridPosition((*it)->getPosition(), iX, iY);
        if ( contained(iX,iY, domain)) (*it)->fluidToParticle(lattice, 0.3);
    }
}



template<typename T, template<typename U> class Descriptor>
void SquareParticleField2D<T,Descriptor>::advanceParticles(Box2D domain, T cutOffValue) {

    //PLB_ASSERT( contained(domain, particleGrid.getBoundingBox()) );
    for(int it = 0; it < particles.size(); it++){
        particles[it]->advance();
        std::cout<<"particle advance."<<std::endl;
    }
}


template<typename T, template<typename U> class Descriptor>
SquareParticleDataTransfer2D<T,Descriptor>& SquareParticleField2D<T,Descriptor>::getDataTransfer() {
    return dataTransfer;
}

template<typename T, template<typename U> class Descriptor>
SquareParticleDataTransfer2D<T,Descriptor> const& SquareParticleField2D<T,Descriptor>::getDataTransfer() const {
    return dataTransfer;
}

template<typename T, template<typename U> class Descriptor>
std::string SquareParticleField2D<T,Descriptor>::getBlockName() {
    return std::string("BulkParticleField2D");
}

template<typename T, template<typename U> class Descriptor>
std::string SquareParticleField2D<T,Descriptor>::basicType() {
    return std::string(NativeType<T>::getName());
}

template<typename T, template<typename U> class Descriptor>
std::string SquareParticleField2D<T,Descriptor>::descriptorType() {
    return std::string(Descriptor<T>::name);
}

}  // namespace plb

#endif  // SQUARE_PARTICLE_FIELD_2D_HH
