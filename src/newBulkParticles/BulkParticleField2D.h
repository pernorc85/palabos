#ifndef BULK_PARTICLE_FIELD_2D_H
#define BULK_PARTICLE_FIELD_2D_H

#include "core/globalDefs.h"
#include "atomicBlock/atomicBlock2D.h"
#include "atomicBlock/dataField2D.h"
#include "newBulkParticles/BulkParticle2D.h"
#include "particles/particleField2D.h"
#include <vector>

namespace plb {
template<typename T, template<typename U> class Descriptor> class BulkParticleField2D;

template<typename T, template<typename U> class Descriptor> class solidBlockLattice2D;

template<typename T, template<typename U> class Descriptor>
class BulkParticleDataTransfer2D : public BlockDataTransfer2D {
public:
    BulkParticleDataTransfer2D(BulkParticleField2D<T,Descriptor>& particleField_);
    virtual plint staticCellSize() const;
    virtual void send(Box2D domain, std::vector<char>& buffer, modif::ModifT kind) const;
    virtual void receive(Box2D domain, std::vector<char> const& buffer, modif::ModifT kind);
    virtual void receive(Box2D domain, std::vector<char> const& buffer, modif::ModifT kind, Dot2D absoluteOffset);
    virtual void receive( Box2D domain, std::vector<char> const& buffer,
                          modif::ModifT kind, std::map<int,std::string> const& foreignIds )
    {
        receive(domain, buffer, kind);
    }
    virtual void attribute(Box2D toDomain, plint deltaX, plint deltaY,
                           AtomicBlock2D const& from, modif::ModifT kind);
    virtual void attribute(Box2D toDomain, plint deltaX, plint deltaY,
                           AtomicBlock2D const& from, modif::ModifT kind, Dot2D absoluteOffset);
private:
    BulkParticleField2D<T,Descriptor>& bulkParticleField;
};

template<typename T, template<typename U> class Descriptor>
class BulkParticleField2D : public ParticleField2D<T,Descriptor> {
public:
    BulkParticleField2D(plint nx, plint ny);
    virtual ~BulkParticleField2D();
    BulkParticleField2D(BulkParticleField2D<T,Descriptor> const& rhs);
    BulkParticleField2D<T,Descriptor>& operator=(BulkParticleField2D<T,Descriptor> const& rhs);
    BulkParticleField2D<T,Descriptor>* clone() const;
    void swap(BulkParticleField2D<T,Descriptor>& rhs);
    /// Add a particle if it is part of the domain, else delete it.
    /** This method with domain-argument is provided here exclusively,
     *    because it may not be easy for an outside instance to decide
     *    autonomously whether a particle is inside a domain or not
     *    (because domains are enlarged by half a lattice site).
     **/
    virtual void addParticle(Box2D domain, Particle2D<T,Descriptor>* particle) {};
    virtual void addParticle(Box2D domain, BulkParticle2D<T,Descriptor>* particle);
    /// Remove all particles found in the indicated domain.
    virtual void removeParticles(Box2D domain);
    /// Remove all particles of a certain tag found in the indicated domain.
    virtual void removeParticles(Box2D domain, plint tag);
    /// Return all particles found in the indicated domain.
    virtual void findParticles(Box2D domain,
                               std::vector<Particle2D<T,Descriptor>*>& found);
    /// Return all particles found in the indicated domain; const version
    virtual void findParticles(Box2D domain,
                               std::vector<Particle2D<T,Descriptor> const*>& found) const;
    /// Execute velocity-particle interaction for all particles contained in the domain.
    virtual void velocityToParticleCoupling(Box2D domain, TensorField2D<T,2>& velocity, T scaling=0.);
    /// Execute velocity-particle interaction for all particles contained in the domain.
    virtual void rhoBarJtoParticleCoupling(Box2D domain, NTensorField2D<T>& rhoBarJ, bool velIsJ, T scaling=0.);
    /// Execute fluid-particle interaction for all particles contained in the domain.
    virtual void fluidToParticleCoupling(Box2D domain, solidBlockLattice2D<T,Descriptor>& lattice);
    virtual void particleToFluidCoupling(BlockLattice2D<T,Descriptor>& lattice);
    /// Advance all particles contained in the domain. When the speed of a particle drops
    ///   below sqrt(cutOffValue), the particle is eliminated. Negative cutOffValue means
    ///   no cutoff.
    virtual void advanceParticles(Box2D domain, T cutOffValue=-1.);
    virtual identifiers::BlockId getBlockId() const { return identifiers::ParticleId; }
public:
    virtual BulkParticleDataTransfer2D<T,Descriptor>& getDataTransfer();
    virtual BulkParticleDataTransfer2D<T,Descriptor> const& getDataTransfer() const;
    static std::string getBlockName();
    static std::string basicType();
    static std::string descriptorType();
public:
    static plint nearestCell(T pos);
    std::vector<BulkParticle2D<T,Descriptor>*>  particles;
    BulkParticleDataTransfer2D<T,Descriptor> dataTransfer;
};

}  // namespace plb

#endif  // BULK_PARTICLE_FIELD_2D_H

