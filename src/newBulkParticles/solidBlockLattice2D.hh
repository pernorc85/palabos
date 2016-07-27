#ifndef SOLID_BLOCK_LATTICE_2D_HH
#define SOLID_BLOCK_LATTICE_2D_HH

#include "newBulkParticles/solidBlockLattice2D.h"
#include "core/dynamics.h"
#include "core/cell.h"
#include "latticeBoltzmann/latticeTemplates.h"
#include "latticeBoltzmann/indexTemplates.h"
#include "core/util.h"
#include "core/latticeStatistics.h"
#include "core/dynamicsIdentifiers.h"
#include "core/plbProfiler.h"
#include <algorithm>
#include <typeinfo>
#include <cmath>

namespace plb {

////////////////////// Class BlockLattice2D /////////////////////////

/** \param nx_ lattice width (first index)
 *  \param ny_ lattice height (second index)
 */
template<typename T, template<typename U> class Descriptor>
solidBlockLattice2D<T,Descriptor>::solidBlockLattice2D (
        plint nx_, plint ny_,
        Dynamics<T,Descriptor>* backgroundDynamics_ )
    : BlockLattice2D<T,Descriptor>(nx_, ny_, backgroundDynamics_),
      dataTransfer(*this)
{
    plint nx = this->getNx();
    plint ny = this->getNy();
    // Allocate memory and attribute dynamics.
    allocateAndInitialize();
    this->backgroundDynamics->report();
    for (plint iX=0; iX<nx; ++iX) {
        for (plint iY=0; iY<ny; ++iY) {
            solidAttributeDynamics(iX, iY, this->backgroundDynamics);
            grid[iX][iY].defineSolidFrac(0);
            grid[iX][iY].defineBoundaryType(-1);
        }
    }
}




/** During destruction, the memory for the lattice and the contained
 * cells is released. However, the dynamics objects pointed to by
 * the cells must be deleted manually by the user.
 */

template<typename T, template<typename U> class Descriptor>
solidBlockLattice2D<T,Descriptor>::~solidBlockLattice2D()
{
    releaseMemory();
}


/** The whole data of the lattice is duplicated. This includes
 * both particle distribution function and external fields.
 * \warning The dynamics objects and internalProcessors are not copied
 * \param rhs the lattice to be duplicated
 */

template<typename T, template<typename U> class Descriptor>
solidBlockLattice2D<T,Descriptor>::solidBlockLattice2D(solidBlockLattice2D<T,Descriptor> const& rhs)
    : BlockLattice2D<T,Descriptor>(rhs)
{
    plint nx = this->getNx();
    plint ny = this->getNy();
    for (plint iX=0; iX<nx; ++iX) {
        for (plint iY=0; iY<ny; ++iY) {
            grid[iX][iY].defineSolidFrac(rhs.grid[iX][iY].solidFrac);
            grid[iX][iY].defineSolidDensity(rhs.grid[iX][iY].solidDensity);
            grid[iX][iY].defineSolidVelocity(rhs.grid[iX][iY].solidVelocity);
        }
    }
}



/** The current lattice is deallocated, then the lattice from the rhs
 * is duplicated. This includes both particle distribution function
 * and external fields.
 * \warning The dynamics objects and internalProcessors are not copied
 * \param rhs the lattice to be duplicated
 */
template<typename T, template<typename U> class Descriptor>
solidBlockLattice2D<T,Descriptor>& solidBlockLattice2D<T,Descriptor>::operator= (
        solidBlockLattice2D<T,Descriptor> const& rhs )
{
    solidBlockLattice2D<T,Descriptor> tmp(rhs);
    swap(tmp);
    return *this;
}

/** The swap is efficient, in the sense that only pointers to the
 * lattice are copied, and not the lattice itself.
 */
/*
template<typename T, template<typename U> class Descriptor>
void BlockLattice2D<T,Descriptor>::swap(BlockLattice2D& rhs) {
    BlockLatticeBase2D<T,Descriptor>::swap(rhs);
    AtomicBlock2D::swap(rhs);
    std::swap(backgroundDynamics, rhs.backgroundDynamics);
    std::swap(rawData, rhs.rawData);
    std::swap(grid, rhs.grid);
}
*/

/*
template<typename T, template<typename U> class Descriptor>
void BlockLattice2D<T,Descriptor>::specifyStatisticsStatus (
        Box2D domain, bool status )
{
    // Make sure domain is contained within current lattice
    PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            grid[iX][iY].specifyStatisticsStatus(status);
        }
    }
}
*/


template<typename T, template<typename U> class Descriptor>
void solidBlockLattice2D<T,Descriptor>::collide(Box2D domain) {
    PLB_PRECONDITION( (plint)Descriptor<T>::q==(plint)Descriptor<T>::numPop );
    // Make sure domain is contained within current lattice
    std::cout<<"solidBlockLattice2D collide(Box2D)"<<std::endl;
    PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            grid[iX][iY].collide(this->getInternalStatistics());
            grid[iX][iY].revert();
        }
    }
}



/** \sa collide(int,int,int,int) */
/*
template<typename T, template<typename U> class Descriptor>
void BlockLattice2D<T,Descriptor>::collide() {
    collide(this->getBoundingBox());
}
*/

/** The distribution functions never leave the rectangular domain. On the
 * domain boundaries, the (outgoing) distribution functions that should
 * be streamed outside are simply left untouched.
 * The finalization of an iteration step is not automatically executed,
 * as it is in the method stream(). If you want it to be executed, you
 * must explicitly call the methods finalizeIteration() and
 * executeInternalProcessors().
 * \sa stream()
 */
/*
template<typename T, template<typename U> class Descriptor>
void BlockLattice2D<T,Descriptor>::stream(Box2D domain) {
    // Make sure domain is contained within current lattice
    PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );

    static const plint vicinity = Descriptor<T>::vicinity;

    bulkStream( Box2D(domain.x0+vicinity,domain.x1-vicinity,
                      domain.y0+vicinity,domain.y1-vicinity) );

    boundaryStream(domain, Box2D(domain.x0,domain.x0+vicinity-1,
                                 domain.y0,domain.y1));
    boundaryStream(domain, Box2D(domain.x1-vicinity+1,domain.x1,
                                 domain.y0,domain.y1));
    boundaryStream(domain, Box2D(domain.x0+vicinity,domain.x1-vicinity,
                                 domain.y0,domain.y0+vicinity-1));
    boundaryStream(domain, Box2D(domain.x0+vicinity,domain.x1-vicinity,
                                 domain.y1-vicinity+1,domain.y1));
}
*/

/** At the end of this method, the methods finalizeIteration() and
 * executeInternalProcessors() are automatically invoked.
 * \sa stream(int,int,int,int)
 */
/*
template<typename T, template<typename U> class Descriptor>
void BlockLattice2D<T,Descriptor>::stream()
{
    stream(this->getBoundingBox());

    implementPeriodicity();

    this->executeInternalProcessors();
    this->evaluateStatistics();
    this->incrementTime();
}
*/

/** This operation is more efficient than a successive application of
 * collide(int,int,int,int) and stream(int,int,int,int), because memory
 * is traversed only once instead of twice.
 * The finalization of an iteration step is not automatically invoked by this
 * method, as it is in the method stream(). If you want it to be executed, you
 * must explicitly call the methods finalizeIteration() and
 * executeInternalProcessors().
 * \sa collideAndStream()
 */
template<typename T, template<typename U> class Descriptor>
void solidBlockLattice2D<T,Descriptor>::collideAndStream(Box2D domain)
{
    PLB_PRECONDITION( (plint)Descriptor<T>::q==(plint)Descriptor<T>::numPop );
    // Make sure domain is contained within current lattice
    PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );
    std::cout<<"Enter solidBlockLattice2D collideAndStream."<<std::endl;
    //global::profiler().start("collStream");
    //global::profiler().increment("collStreamCells", domain.nCells());

    static const plint vicinity = Descriptor<T>::vicinity;

    // First, do the collision on cells within a boundary envelope of width
    // equal to the range of the lattice vectors (e.g. 1 for D2Q9)
    this->collide(Box2D(domain.x0,domain.x0+vicinity-1, domain.y0,domain.y1));
    this->collide(Box2D(domain.x1-vicinity+1,domain.x1, domain.y0,domain.y1));
    this->collide(Box2D(domain.x0+vicinity,domain.x1-vicinity, domain.y0,domain.y0+vicinity-1));
    this->collide(Box2D(domain.x0+vicinity,domain.x1-vicinity, domain.y1-vicinity+1,domain.y1));

    std::cout<<"boundary collide finished"<<std::endl;

    // Then, do the efficient collideAndStream algorithm in the bulk,
    // excluding the envelope (this is efficient because there is no
    // if-then-else statement within the loop, given that the boundary
    // region is excluded)

    this->collide(Box2D (domain.x0+vicinity, domain.x1-vicinity,
                   domain.y0+vicinity, domain.y1-vicinity));
    std::cout<<"bulk collide finished"<<std::endl;
    system("PAUSE");

    bulkStream(Box2D (domain.x0+vicinity, domain.x1-vicinity,
                      domain.y0+vicinity, domain.y1-vicinity));
    std::cout<<"bulkStream finished"<<std::endl;
//    bulkCollideAndStream(Box2D(domain.x0+vicinity,domain.x1-vicinity,
//                               domain.y0+vicinity,domain.y1-vicinity));

    // Finally, do streaming in the boundary envelope to conclude the
    // collision-stream cycle
    this->boundaryStream(domain, Box2D(domain.x0,domain.x0+vicinity-1,
                                 domain.y0,domain.y1));
    this->boundaryStream(domain, Box2D(domain.x1-vicinity+1,domain.x1,
                                 domain.y0,domain.y1));
    this->boundaryStream(domain, Box2D(domain.x0+vicinity,domain.x1-vicinity,
                                 domain.y0,domain.y0+vicinity-1));
    this->boundaryStream(domain, Box2D(domain.x0+vicinity,domain.x1-vicinity,
                                 domain.y1-vicinity+1,domain.y1));
    std::cout<<"boundaryStream finished"<<std::endl;
    system("PAUSE");
    //global::profiler().stop("collStream");
}

/** At the end of this method, the methods finalizeIteration() and
 * executeInternalProcessors() are automatically invoked.
 * \sa collideAndStream(int,int,int,int) */

template<typename T, template<typename U> class Descriptor>
void solidBlockLattice2D<T,Descriptor>::collideAndStream() {
    collideAndStream(this->getBoundingBox());

    this->implementPeriodicity();

    this->executeInternalProcessors();
    this->evaluateStatistics();
    this->incrementTime();
}


/*
template<typename T, template<typename U> class Descriptor>
void BlockLattice2D<T,Descriptor>::incrementTime() {
    this->getTimeCounter().incrementTime();
}
*/

template<typename T, template<typename U> class Descriptor>
void solidBlockLattice2D<T,Descriptor>::allocateAndInitialize() {
    this->getInternalStatistics().subscribeAverage(); // Subscribe average rho-bar
    this->getInternalStatistics().subscribeAverage(); // Subscribe average uSqr
    this->getInternalStatistics().subscribeMax();     // Subscribe max uSqr
    plint nx = this->getNx();
    plint ny = this->getNy();
    rawData = new solidCell<T,Descriptor> [nx*ny];
    grid    = new solidCell<T,Descriptor>* [nx];
    for (plint iX=0; iX<nx; ++iX) {
        grid[iX] = rawData + iX*ny;
    }
}


template<typename T, template<typename U> class Descriptor>
void solidBlockLattice2D<T,Descriptor>::attributeSolidFracAndDensity(BulkParticleField2D<T,Descriptor>& bulkParticleField_) {
    double solidFrac;

    std::cout<<"AttributeSolidFrac"<<std::endl;
//    for(plint i=0; i<this->getNx(); i++){
//        for(plint j=0; j<this->getNy(); j++){
//            this->getSolidCell(i,j).defineSolidFrac(0);
//            this->getSolidCell(i,j).defineSolidDensity(1);
//            this->getSolidCell(i,j).defineBoundaryType(-1);
//        }
//    }

    for(pluint iParticle = 0; iParticle < bulkParticleField_.particles.size(); iParticle++){
        //std::cout<<"fluidBoundary"<<std::endl;
        for(std::vector<Array<int,2> >::iterator i= bulkParticleField_.particles[iParticle]->fluidBoundary.begin();
                                                 i!=bulkParticleField_.particles[iParticle]->fluidBoundary.end(); ++i){
            Array<int,2> cellPosition = (*i);
            solidFrac = 0;
            
            if(contained(cellPosition[0], cellPosition[1], this->getBoundingBox())){
                computeSolidFraction(cellPosition, bulkParticleField_.particles[iParticle]->position,
                                 bulkParticleField_.particles[iParticle]->radius, solidFrac);
                //std::cout<<"solidFrac="<<solidFrac<<std::endl;
                this->getSolidCell((*i)[0], (*i)[1]).defineSolidFrac(solidFrac);
                this->getSolidCell((*i)[0], (*i)[1]).defineSolidDensity(bulkParticleField_.particles[iParticle]->getDensity());
                this->getSolidCell((*i)[0], (*i)[1]).defineBoundaryType(2);
            }
        }
        //std::cout<<"solidBoundary"<<std::endl;
        for(std::vector<Array<int,2> >::iterator i= bulkParticleField_.particles[iParticle]->solidBoundary.begin();
                                                 i!=bulkParticleField_.particles[iParticle]->solidBoundary.end(); ++i){
            Array<int,2> cellPosition = (*i);
            solidFrac = 0;
            if(contained(cellPosition[0], cellPosition[1], this->getBoundingBox())){
                computeSolidFraction(cellPosition, bulkParticleField_.particles[iParticle]->position,
                                 bulkParticleField_.particles[iParticle]->radius, solidFrac);
                //std::cout<<"solidFrac="<<solidFrac<<std::endl;
                this->getSolidCell((*i)[0], (*i)[1]).defineSolidFrac(solidFrac);
                this->getSolidCell((*i)[0], (*i)[1]).defineSolidDensity(bulkParticleField_.particles[iParticle]->getDensity());
                this->getSolidCell((*i)[0], (*i)[1]).defineBoundaryType(1);
            }
        }
        for(std::vector<Array<int,2> >::iterator i= bulkParticleField_.particles[iParticle]->internalSolid.begin();
                                                 i!=bulkParticleField_.particles[iParticle]->internalSolid.end(); ++i){
            Array<int,2> cellPosition = (*i);
            if(contained(cellPosition[0], cellPosition[1], this->getBoundingBox())){
                 this->getSolidCell((*i)[0], (*i)[1]).defineSolidFrac(1);
                 this->getSolidCell((*i)[0], (*i)[1]).defineSolidDensity(bulkParticleField_.particles[iParticle]->getDensity());
                 this->getSolidCell((*i)[0], (*i)[1]).defineBoundaryType(0);
            }
        }
    }
}


template<typename T, template<typename U> class Descriptor>
void solidBlockLattice2D<T,Descriptor>::attributeSolidVelocity(BulkParticleField2D<T,Descriptor>& bulkParticleField_) {

    std::cout<<"AttributeSolidVelocity"<<std::endl;
//    for(pluint iX = 0; iX < this->getNx(); iX++){
//        for(pluint iY = 0; iY < this->getNy(); iY++){
//            iniCellAtEquilibrium(this->getSolidCell(iX, iY), 1., Array<T,2> (10,0));
//        }
//    }

    for(pluint iParticle = 0; iParticle < bulkParticleField_.particles.size(); iParticle++){
        for(std::vector<Array<int,2> >::iterator i= bulkParticleField_.particles[iParticle]->fluidBoundary.begin();
                                                 i!=bulkParticleField_.particles[iParticle]->fluidBoundary.end(); ++i){
            Array<int,2> cellPosition = (*i);
            T AngularVelocity = bulkParticleField_.particles[iParticle]->getAngularVelocity();
            Array<T,2> radius;
            radius[0] = (T)cellPosition[0] - bulkParticleField_.particles[iParticle]->getPosition()[0];
            radius[1] = (T)cellPosition[1] - bulkParticleField_.particles[iParticle]->getPosition()[1];
            Array<T,2> additionalVelocity;
            additionalVelocity[0] = radius[1] * AngularVelocity;
            additionalVelocity[1] =-radius[0] * AngularVelocity;
            if(contained(cellPosition[0], cellPosition[1], this->getBoundingBox())){
                this->getSolidCell((*i)[0], (*i)[1]).defineSolidVelocity(additionalVelocity
                                        + bulkParticleField_.particles[iParticle]->getVelocity() );
            }
        }
        for(std::vector<Array<int,2> >::iterator i= bulkParticleField_.particles[iParticle]->solidBoundary.begin();
                                                 i!=bulkParticleField_.particles[iParticle]->solidBoundary.end(); ++i){
            Array<int,2> cellPosition = (*i);
            T AngularVelocity = bulkParticleField_.particles[iParticle]->getAngularVelocity();
            Array<T,2> radius;
            radius[0] = (T)cellPosition[0] - bulkParticleField_.particles[iParticle]->getPosition()[0];
            radius[1] = (T)cellPosition[1] - bulkParticleField_.particles[iParticle]->getPosition()[1];
            Array<T,2> additionalVelocity;
            additionalVelocity[0] = radius[1] * AngularVelocity;
            additionalVelocity[1] =-radius[0] * AngularVelocity;
            if(contained(cellPosition[0], cellPosition[1], this->getBoundingBox())){
                this->getSolidCell((*i)[0], (*i)[1]).defineSolidVelocity(additionalVelocity
                                        + bulkParticleField_.particles[iParticle]->getVelocity() );
            }
        }
        for(std::vector<Array<int,2> >::iterator i= bulkParticleField_.particles[iParticle]->internalSolid.begin();
                                                 i!=bulkParticleField_.particles[iParticle]->internalSolid.end(); ++i){
            Array<int,2> cellPosition = (*i);
            T AngularVelocity = bulkParticleField_.particles[iParticle]->getAngularVelocity();
            Array<T,2> radius;
            radius[0] = (T)cellPosition[0] - bulkParticleField_.particles[iParticle]->getPosition()[0];
            radius[1] = (T)cellPosition[1] - bulkParticleField_.particles[iParticle]->getPosition()[1];
            Array<T,2> additionalVelocity;
            additionalVelocity[0] = radius[1] * AngularVelocity;
            additionalVelocity[1] =-radius[0] * AngularVelocity;
            if(contained(cellPosition[0], cellPosition[1], this->getBoundingBox())){
                this->getSolidCell((*i)[0], (*i)[1]).defineSolidVelocity(additionalVelocity
                                        + bulkParticleField_.particles[iParticle]->getVelocity() );
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
void solidBlockLattice2D<T,Descriptor>::attributeSolidFracAndDensity(SquareParticleField2D<T,Descriptor>& squareParticleField_) {
    double solidFrac;

    std::cout<<"AttributeSolidFrac from SquareParticleField2D"<<std::endl;
    for(plint i=0; i<this->getNx(); i++){
        for(plint j=0; j<this->getNy(); j++){
            this->getSolidCell(i,j).defineSolidFrac(0);
            this->getSolidCell(i,j).defineSolidDensity(1);
            this->getSolidCell(i,j).defineBoundaryType(-1);
        }
    }

    for(pluint iParticle = 0; iParticle < squareParticleField_.particles.size(); iParticle++){
        for(std::vector<Array<int,2> >::iterator i= squareParticleField_.particles[iParticle]->fluidBoundary.begin();
                                                 i!=squareParticleField_.particles[iParticle]->fluidBoundary.end(); ++i){
            Array<int,2> cellPosition = (*i);
            solidFrac = 0;
            if(contained(cellPosition[0], cellPosition[1], this->getBoundingBox())){
                this->getSolidCell((*i)[0], (*i)[1]).defineSolidFrac(solidFrac);
                this->getSolidCell((*i)[0], (*i)[1]).defineSolidDensity(squareParticleField_.particles[iParticle]->getDensity());
                this->getSolidCell((*i)[0], (*i)[1]).defineBoundaryType(2);
            }
        }
        for(std::vector<Array<int,2> >::iterator i= squareParticleField_.particles[iParticle]->solidBoundary.begin();
                                                 i!=squareParticleField_.particles[iParticle]->solidBoundary.end(); ++i){
            Array<int,2> cellPosition = (*i);
            solidFrac = 1;
            if(contained(cellPosition[0], cellPosition[1], this->getBoundingBox())){
                this->getSolidCell((*i)[0], (*i)[1]).defineSolidFrac(solidFrac);
                this->getSolidCell((*i)[0], (*i)[1]).defineSolidDensity(squareParticleField_.particles[iParticle]->getDensity());
                this->getSolidCell((*i)[0], (*i)[1]).defineBoundaryType(1);
            }
        }
        for(std::vector<Array<int,2> >::iterator i= squareParticleField_.particles[iParticle]->internalSolid.begin();
                                                 i!=squareParticleField_.particles[iParticle]->internalSolid.end(); ++i){
            Array<int,2> cellPosition = (*i);
            if(contained(cellPosition[0], cellPosition[1], this->getBoundingBox())){
                 this->getSolidCell((*i)[0], (*i)[1]).defineSolidFrac(1);
                 this->getSolidCell((*i)[0], (*i)[1]).defineSolidDensity(squareParticleField_.particles[iParticle]->getDensity());
                 this->getSolidCell((*i)[0], (*i)[1]).defineBoundaryType(0);
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
void solidBlockLattice2D<T,Descriptor>::attributeSolidVelocity(SquareParticleField2D<T,Descriptor>& squareParticleField_) {

    std::cout<<"AttributeSolidVelocity from SquareParticleField2D"<<std::endl;
//    for(pluint iX = 0; iX < this->getNx(); iX++){
//        for(pluint iY = 0; iY < this->getNy(); iY++){
//            iniCellAtEquilibrium(this->getSolidCell(iX, iY), 1., Array<T,2> (10,0));
//        }
//    }

    for(pluint iParticle = 0; iParticle < squareParticleField_.particles.size(); iParticle++){
        for(std::vector<Array<int,2> >::iterator i= squareParticleField_.particles[iParticle]->fluidBoundary.begin();
                                                 i!=squareParticleField_.particles[iParticle]->fluidBoundary.end(); ++i){
            Array<int,2> cellPosition = (*i);
            T AngularVelocity = squareParticleField_.particles[iParticle]->getAngularVelocity();
            Array<T,2> radius;
            radius[0] = (T)cellPosition[0] - squareParticleField_.particles[iParticle]->getPosition()[0];
            radius[1] = (T)cellPosition[1] - squareParticleField_.particles[iParticle]->getPosition()[1];
            Array<T,2> additionalVelocity;
            additionalVelocity[0] = radius[1] * AngularVelocity;
            additionalVelocity[1] =-radius[0] * AngularVelocity;
            if(contained(cellPosition[0], cellPosition[1], this->getBoundingBox())){
                this->getSolidCell((*i)[0], (*i)[1]).defineSolidVelocity(additionalVelocity
                                        + squareParticleField_.particles[iParticle]->getVelocity() );
            }
            //Don't initialize this cell since it is fluid boundary, it should still possess the fluid velocity.
            //iniCellAtEquilibrium(this->getSolidCell((*i)[0], (*i)[1]), 1.,
            //                (additionalVelocity + bulkParticleField_.particles[iParticle]->getVelocity() )
            //                 * (this->getSolidCell((*i)[0], (*i)[1]).getSolidFrac()) );
        }
        for(std::vector<Array<int,2> >::iterator i= squareParticleField_.particles[iParticle]->solidBoundary.begin();
                                                 i!=squareParticleField_.particles[iParticle]->solidBoundary.end(); ++i){
            Array<int,2> cellPosition = (*i);
            T AngularVelocity = squareParticleField_.particles[iParticle]->getAngularVelocity();
            Array<T,2> radius;
            radius[0] = (T)cellPosition[0] - squareParticleField_.particles[iParticle]->getPosition()[0];
            radius[1] = (T)cellPosition[1] - squareParticleField_.particles[iParticle]->getPosition()[1];
            Array<T,2> additionalVelocity;
            additionalVelocity[0] = radius[1] * AngularVelocity;
            additionalVelocity[1] =-radius[0] * AngularVelocity;
            if(contained(cellPosition[0], cellPosition[1], this->getBoundingBox())){
                this->getSolidCell((*i)[0], (*i)[1]).defineSolidVelocity(additionalVelocity
                                        + squareParticleField_.particles[iParticle]->getVelocity() );
            }
            //Don't initialize this cell since it is fluid boundary, it should still possess the fluid velocity.
            //iniCellAtEquilibrium(this->getSolidCell((*i)[0], (*i)[1]), 1.,
            //                (additionalVelocity + bulkParticleField_.particles[iParticle]->getVelocity() )
            //                 * (this->getSolidCell((*i)[0], (*i)[1]).getSolidFrac()) );
        }
        for(std::vector<Array<int,2> >::iterator i= squareParticleField_.particles[iParticle]->internalSolid.begin();
                                                 i!=squareParticleField_.particles[iParticle]->internalSolid.end(); ++i){
            Array<int,2> cellPosition = (*i);
            T AngularVelocity = squareParticleField_.particles[iParticle]->getAngularVelocity();
            Array<T,2> radius;
            radius[0] = (T)cellPosition[0] - squareParticleField_.particles[iParticle]->getPosition()[0];
            radius[1] = (T)cellPosition[1] - squareParticleField_.particles[iParticle]->getPosition()[1];
            Array<T,2> additionalVelocity;
            additionalVelocity[0] = radius[1] * AngularVelocity;
            additionalVelocity[1] =-radius[0] * AngularVelocity;
            if(contained(cellPosition[0], cellPosition[1], this->getBoundingBox())){
                this->getSolidCell((*i)[0], (*i)[1]).defineSolidVelocity(additionalVelocity
                                        + squareParticleField_.particles[iParticle]->getVelocity() );
            }
            //iniCellAtEquilibrium(this->getSolidCell((*i)[0], (*i)[1]), 1.,
            //                additionalVelocity + bulkParticleField_.particles[iParticle]->getVelocity() );
        }
    }
}

template<typename T, template<typename U> class Descriptor>
void solidBlockLattice2D<T,Descriptor>::releaseMemory() {
    plint nx = this->getNx();
    plint ny = this->getNy();
    for (plint iX=0; iX<nx; ++iX) {
        for (plint iY=0; iY<ny; ++iY) {
            Dynamics<T,Descriptor>* dynamics = &grid[iX][iY].getDynamics();
            if (dynamics != this->backgroundDynamics) {
                delete dynamics;
            }
        }
    }
    delete this->backgroundDynamics;
    delete [] rawData;
    delete [] grid;
}


template<typename T, template<typename U> class Descriptor>
void solidBlockLattice2D<T,Descriptor>::solidAttributeDynamics(plint iX, plint iY, Dynamics<T,Descriptor>* dynamics) {
    Dynamics<T,Descriptor>* previousDynamics = &grid[iX][iY].getDynamics();
    if (previousDynamics != this->backgroundDynamics) {
        delete previousDynamics;
    }
    grid[iX][iY].solidAttributeDynamics(dynamics);
}

/*
template<typename T, template<typename U> class Descriptor>
Dynamics<T,Descriptor>& BlockLattice2D<T,Descriptor>::getBackgroundDynamics() {
    return *backgroundDynamics;
}
*/

/*
template<typename T, template<typename U> class Descriptor>
Dynamics<T,Descriptor> const& BlockLattice2D<T,Descriptor>::getBackgroundDynamics() const {
    return *backgroundDynamics;
}
*/

/** This method is slower than bulkStream(int,int,int,int), because one needs
 * to verify which distribution functions are to be kept from leaving
 * the domain.
 * \sa stream(int,int,int,int)
 * \sa stream()
 */
/*
template<typename T, template<typename U> class Descriptor>
void BlockLattice2D<T,Descriptor>::boundaryStream(Box2D bound, Box2D domain) {
    // Make sure bound is contained within current lattice
    PLB_PRECONDITION( contained(bound, this->getBoundingBox()) );
    // Make sure domain is contained within bound
    PLB_PRECONDITION( contained(domain, bound) );

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iPop=1; iPop<=Descriptor<T>::q/2; ++iPop) {
                plint nextX = iX + Descriptor<T>::c[iPop][0];
                plint nextY = iY + Descriptor<T>::c[iPop][1];
                if (nextX>=bound.x0 && nextX<=bound.x1 && nextY>=bound.y0 && nextY<=bound.y1) {
                    std::swap(grid[iX][iY][iPop+Descriptor<T>::q/2],
                              grid[nextX][nextY][iPop]);
                }
            }
        }
    }
}
*/

/** This method is faster than boundaryStream(int,int,int,int), but it
 * is erroneous when applied to boundary cells.
 * \sa stream(int,int,int,int)
 * \sa stream()
 */
template<typename T, template<typename U> class Descriptor>
void solidBlockLattice2D<T,Descriptor>::bulkStream(Box2D domain) {
    // Make sure domain is contained within current lattice
    PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {

            for (plint iPop=1; iPop<=Descriptor<T>::q/2; ++iPop) {
                plint nextX = iX + Descriptor<T>::c[iPop][0];
                plint nextY = iY + Descriptor<T>::c[iPop][1];
                    //if(this->getSolidCell(nextX,nextY).getSolidFrac() < 1 && this->getSolidCell(iX,iY).getSolidFrac() < 1)
                    std::swap(grid[iX][iY][iPop+Descriptor<T>::q/2],
                          grid[nextX][nextY][iPop]);
            }


        }
    }
}


/*
template<typename T, template<typename U> class Descriptor>
void BlockLattice2D<T,Descriptor>::periodicDomain(Box2D domain) {
    plint nx = this->getNx();
    plint ny = this->getNy();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iPop=1; iPop<Descriptor<T>::q; ++iPop) {
                plint prevX = iX - Descriptor<T>::c[iPop][0];
                plint prevY = iY - Descriptor<T>::c[iPop][1];
                if ( (prevX>=0 && prevX<nx) &&
                     (prevY>=0 && prevY<ny) )
                {
                    plint nextX = (iX+nx)%nx;
                    plint nextY = (iY+ny)%ny;
                    std::swap (
                        grid[prevX][prevY][indexTemplates::opposite<Descriptor<T> >(iPop)],
                        grid[nextX][nextY][iPop] );
                }
            }
        }
    }
}
*/

/*
template<typename T, template<typename U> class Descriptor>
void BlockLattice2D<T,Descriptor>::implementPeriodicity() {
    static const plint vicinity = Descriptor<T>::vicinity;
    plint maxX = this->getNx()-1;
    plint maxY = this->getNy()-1;
    // Periodicity of edges orthogonal to x-axis.
    periodicDomain(Box2D(-vicinity,-1,0,maxY));
    // Periodicity of edges orthogonal to y-axis.
    periodicDomain(Box2D(0,maxX,-vicinity,-1));
    // Periodicity between (-1,-1) and (+1,+1) corner.
    periodicDomain(Box2D(-vicinity,-1,-vicinity,-1));
    // Periodicity between (-1,+1) and (+1,-1) corner.
    periodicDomain(Box2D(-vicinity,-1,maxY+1,maxY+vicinity));
}
*/


}  // namespace plb

#endif  // BLOCK_LATTICE_2D_HH
