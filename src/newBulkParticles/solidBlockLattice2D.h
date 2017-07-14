#ifndef SOLID_BLOCK_LATTICE_2D_H
#define SOLID_BLOCK_LATTICE_2D_H

#include "core/globalDefs.h"
#include "core/plbDebug.h"
#include "core/cell.h"
#include "atomicBlock/dataField2D.h"
#include "core/blockLatticeBase2D.h"
#include "atomicBlock/atomicBlock2D.h"
#include "core/blockIdentifiers.h"
#include "atomicBlock/blockLattice2D.h"
#include "newBulkParticles/BulkParticleField2D.h"
#include "newBulkParticles/SquareParticleField2D.h"
#include "newBulkParticles/solidCell.h"
#include <vector>
#include <map>


namespace plb {

template<typename T, template<typename U> class Descriptor> struct Dynamics;
template<typename T, template<typename U> class Descriptor> class BlockLattice2D;


/** A block lattice contains a regular array of Cell objects and
 * some useful methods to execute the LB dynamics on the lattice.
 *
 * This class is not intended to be derived from.
 */
template<typename T, template<typename U> class Descriptor> class BlockLatticeDataTransfer2D;

template<typename T, template<typename U> class Descriptor>
class solidBlockLattice2D : public BlockLattice2D<T,Descriptor> {
public:
    /// Construction of an nx_ by ny_ lattice
    solidBlockLattice2D(plint nx_, plint ny_, Dynamics<T,Descriptor>* backgroundDynamics);
    /// Constructino according to a bulkParticleField
    solidBlockLattice2D(BulkParticleField2D<T,Descriptor>& myBulkParticleField_, Dynamics<T,Descriptor>* backgroundDynamics);
    /// Destruction of the lattice
    ~solidBlockLattice2D();
    /// Copy construction
    solidBlockLattice2D(solidBlockLattice2D<T,Descriptor> const& rhs);
    /// Copy assignment
    solidBlockLattice2D& operator=(solidBlockLattice2D<T,Descriptor> const& rhs);
    /// Swap the content of two BlockLattices
    void swap(solidBlockLattice2D& rhs);
public:
    /// Read/write access to lattice cells

    virtual solidCell<T,Descriptor>& getSolidCell(plint iX, plint iY) {
        PLB_PRECONDITION(iX<this->getNx());
        PLB_PRECONDITION(iY<this->getNy());
        return grid[iX][iY];
    }

    /// Read only access to lattice cells
    virtual solidCell<T,Descriptor> const& getSolidCell(plint iX, plint iY) const {
        PLB_PRECONDITION(iX<this->getNx());
        PLB_PRECONDITION(iY<this->getNy());
        return grid[iX][iY];
    }

    /// Specify wheter statistics measurements are done on given rect. domain
    //virtual void specifyStatisticsStatus(Box2D domain, bool status);
    virtual void collide(Box2D domain);
    //virtual void collide();
//    virtual void stream(Box2D domain);
//    virtual void stream();
    /// Apply first collision, then streaming step to a rectangular domain
    virtual void collideAndStream(Box2D domain);
    /// Apply first collision, then streaming step to the whole domain
    virtual void collideAndStream();
    /// Increment time counter
    //virtual void incrementTime();
    /// Get access to data transfer between blocks
    //virtual BlockLatticeDataTransfer2D<T,Descriptor>& getDataTransfer();
    /// Get access to data transfer between blocks (const version)
    //virtual BlockLatticeDataTransfer2D<T,Descriptor> const& getDataTransfer() const;
public:
    /// Attribute dynamics to a cell.
    void solidAttributeDynamics(plint iX, plint iY, Dynamics<T,Descriptor>* dynamics);
    void attributeSolidFracAndDensity(BulkParticleField2D<T,Descriptor>& bulkParticleField_);
    void attributeSolidVelocity(BulkParticleField2D<T,Descriptor>& bulkParticleField_);
    
    void attributeSolidFracAndDensity(SquareParticleField2D<T,Descriptor>& squareParticleField_);
    void attributeSolidVelocity(SquareParticleField2D<T,Descriptor>& squareParticleField_);
//    Dynamics<T,Descriptor>& getBackgroundDynamics();
//    Dynamics<T,Descriptor> const& getBackgroundDynamics() const;
    /// Apply streaming step to bulk (non-boundary) cells
    void bulkStream(Box2D domain);
    /// Apply streaming step to boundary cells
    //void boundaryStream(Box2D bound, Box2D domain);
private:
    /// Helper method for memory allocation
    void allocateAndInitialize();
    void releaseMemory();
    //void implementPeriodicity();
private:
    void periodicDomain(Box2D domain);
private:
    solidCell<T,Descriptor>     *rawData;
    solidCell<T,Descriptor>    **grid;
    BlockLatticeDataTransfer2D<T,Descriptor> dataTransfer;
public:
    static CachePolicy2D& cachePolicy();
template<typename T_, template<typename U_> class Descriptor_>
    friend class ExternalRhoJcollideAndStream2D;
};

template<typename T, template<typename U> class Descriptor>
double getStoredAverageDensity(BlockLattice2D<T,Descriptor> const& blockLattice);

template<typename T, template<typename U> class Descriptor>
double getStoredAverageEnergy(BlockLattice2D<T,Descriptor> const& blockLattice);

template<typename T, template<typename U> class Descriptor>
double getStoredAverageVelocity(BlockLattice2D<T,Descriptor> const& blockLattice);

}  // namespace plb

#endif  // SOLID_BLOCK_LATTICE_2D_H
