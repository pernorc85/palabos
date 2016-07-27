#ifndef BULK_PARTICLE_CHUNK_2D_H
#define BULK_PARTICLE_CHUNK_2D_H

#include "core/globalDefs.h"
#include "atomicBlock/atomicBlock2D.h"
#include "atomicBlock/dataField2D.h"
#include "newBulkParticles/BulkParticle2D.h"
#include "particles/particleField2D.h"
#include "newBulkParticles/BulkParticleField2D.h"
#include <vector>

namespace plb {
template<typename T, template<typename U> class Descriptor> class BulkParticleChunk2D;

template<typename T, template<typename U> class Descriptor>
class BulkParticleChunk2D : public ParticleField2D<T,Descriptor> {
public:
    BulkParticleChunk2D(plint nx, plint ny);
    virtual ~BulkParticleChunk2D();
    BulkParticleChunk2D(BulkParticleChunk2D<T,Descriptor> const& rhs);
    BulkParticleChunk2D<T,Descriptor>& operator=(BulkParticleChunk2D<T,Descriptor> const& rhs);
    BulkParticleChunk2D<T,Descriptor>* clone() const;
    void swap(BulkParticleChunk2D<T,Descriptor>& rhs);
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
    virtual void removeParticles(Box2D domain, plint tag){}
    /// Return all particles found in the indicated domain.
    virtual void findParticles(Box2D domain,
                               std::vector<Particle2D<T,Descriptor>*>& found){}
    /// Return all particles found in the indicated domain; const version
    virtual void findParticles(Box2D domain,
                               std::vector<Particle2D<T,Descriptor> const*>& found) const{}
    ///
    void releaseParticle(plint i, BulkParticleField2D<T,Descriptor>* particleField);
    void releaseParticles(std::vector<int>& idToRelease, BulkParticleField2D<T,Descriptor>* particleField);
    void rigidize();
    /// Execute velocity-particle interaction for all particles contained in the domain.
    virtual void velocityToParticleCoupling(Box2D domain, TensorField2D<T,2>& velocity, T scaling=0.);
    /// Execute velocity-particle interaction for all particles contained in the domain.
    virtual void rhoBarJtoParticleCoupling(Box2D domain, NTensorField2D<T>& rhoBarJ, bool velIsJ, T scaling=0.);
    /// Execute fluid-particle interaction for all particles contained in the domain.
    virtual void fluidToParticleCoupling(Box2D domain, solidBlockLattice2D<T,Descriptor>& lattice);
    virtual void particleToFluidCoupling(BlockLattice2D<T,Descriptor>& lattice);
    virtual void chunkToFluid(BlockLattice2D<T,Descriptor>& lattice);
    virtual void fixedChunkToFluid(BlockLattice2D<T,Descriptor>& fluid,
                                   BulkParticleField2D<T,Descriptor>* particleField);
    /// Advance all particles contained in the domain. When the speed of a particle drops
    ///   below sqrt(cutOffValue), the particle is eliminated. Negative cutOffValue means
    ///   no cutoff.
    virtual void advanceParticles(Box2D domain, T cutOffValue=-1.){}
    virtual void advanceChunk(Box2D domain, T cutOffValue=-1.);
    virtual identifiers::BlockId getBlockId() const { return identifiers::ParticleId; }
    void computeMassCenter();
public:
    virtual BulkParticleDataTransfer2D<T,Descriptor>& getDataTransfer();
    virtual BulkParticleDataTransfer2D<T,Descriptor> const& getDataTransfer() const;
    static std::string getBlockName();
    static std::string basicType();
    static std::string descriptorType();
public:
    static plint nearestCell(T pos);
    std::vector<BulkParticle2D<T,Descriptor>*>  particles;
    std::vector<std::vector<int> > rigidBonds;
    Array<T,2>  massCenter;
};

}  // namespace plb

#endif  // BULK_PARTICLE_FIELD_2D_H

