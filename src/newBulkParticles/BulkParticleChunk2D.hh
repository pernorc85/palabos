#ifndef BULK_PARTICLE_CHUNK_2D_HH
#define BULK_PARTICLE_CHUNK_2D_HH

#include <cmath>
#include "core/globalDefs.h"
#include "particles/particleField2D.h"
#include "newBulkParticles/BulkParticleChunk2D.h"
#include "newBulkParticles/solidBlockLattice2D.h"
#include <utility>

namespace plb {
/* *************** class BulkParticleChunk2D ************************************ */

template<typename T, template<typename U> class Descriptor>
BulkParticleChunk2D<T,Descriptor>::BulkParticleChunk2D(plint nx, plint ny)
    : ParticleField2D<T,Descriptor>(nx,ny)
{ }

template<typename T, template<typename U> class Descriptor>
BulkParticleChunk2D<T,Descriptor>::~BulkParticleChunk2D()
{
    particles.clear();
}



template<typename T, template<typename U> class Descriptor>
BulkParticleChunk2D<T,Descriptor>::BulkParticleChunk2D(BulkParticleChunk2D const& rhs)
{
    particles.resize(rhs.particles.size());
    for (pluint iParticle=0; iParticle<particles.size(); ++iParticle) {
        particles[iParticle] = rhs.particles[iParticle].clone();
    }
}


template<typename T, template<typename U> class Descriptor>
BulkParticleChunk2D<T,Descriptor>&
    BulkParticleChunk2D<T,Descriptor>::operator=(BulkParticleChunk2D<T,Descriptor> const& rhs)
{
    BulkParticleChunk2D<T,Descriptor>(rhs).swap(*this);
    return *this;
}


template<typename T, template<typename U> class Descriptor>
BulkParticleChunk2D<T,Descriptor>*
    BulkParticleChunk2D<T,Descriptor>::clone() const
{
    return new BulkParticleChunk2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void BulkParticleChunk2D<T,Descriptor>::swap(BulkParticleChunk2D<T,Descriptor>& rhs) {
    particles.swap(rhs.particles);
}

template<typename T, template<typename U> class Descriptor>
void BulkParticleChunk2D<T,Descriptor>::addParticle(Box2D domain,  BulkParticle2D<T,Descriptor>* particle) {
    plint iX, iY;
    computeGridPosition(particle->getPosition(), iX, iY);
    bool overlap = 0;

    for(plint i = 0; i < particles.size(); i++) {
        BulkParticle2D<T,Descriptor>& particle2 = *(particles[i]);
        if( pow(particle->getPosition()[0] - particle2.getPosition()[0],2.0) +
            pow(particle->getPosition()[1] - particle2.getPosition()[1],2.0) <
            pow(particle->getRadius() + particle2.getRadius(), 2.0) )
            overlap = 1;
    }

    if( contained(iX,iY, domain) && overlap == 0)//modify
    {
        particles.push_back(particle);
        std::cout<<"particle added to chunk."<<std::endl;
    }
    else {
        delete particle;
        std::cout<<"particle deleted."<<std::endl;
    }
}


template<typename T, template<typename U> class Descriptor>
void BulkParticleChunk2D<T,Descriptor>::removeParticles(Box2D domain) {
    plint iX, iY;
    typename std::vector<BulkParticle2D<T,Descriptor>*>::iterator it = particles.begin();
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
void BulkParticleChunk2D<T,Descriptor>::releaseParticle(plint i, BulkParticleField2D<T,Descriptor>* particleField){
    typename std::vector<BulkParticle2D<T,Descriptor>*>::iterator it = particles.begin()+i-1;
    typename std::vector<std::vector<int> >::iterator itB = rigidBonds.begin()+i-1;
    particleField->addParticle(Box2D(0,1200-1,0,600-1), (*it));
    particles.erase(it);
    (*itB).clear();
    rigidBonds.erase(itB);
}

template<typename T, template<typename U> class Descriptor>
void BulkParticleChunk2D<T,Descriptor>::releaseParticles(std::vector<int>& idToRelease,
                                                         BulkParticleField2D<T,Descriptor>* particleField){
    if(idToRelease.size() == 0)return;
    typename std::vector<BulkParticle2D<T,Descriptor>*>::iterator it = particles.begin();
    typename std::vector<std::vector<int> >::iterator itB = rigidBonds.begin();
    typename std::vector<int>::iterator itRelease = idToRelease.begin();
    printf("in releaseParticles");
    printf("first to release id = %d\n", (*itRelease));
    for(; it != particles.end(); ){
        printf("processing particle with id%d\n", (*it)->getId());
        if((*it)->getId() == (*itRelease)){
            printf("particle with id%d is removed from chunk\n",(*it)->getId());
            system("PAUSE");
            particleField->addParticle(Box2D(0,1200-1,0,600-1), (*it));
            it = particles.erase(it);
            (*itB).clear();
            itB = rigidBonds.erase(itB);
            if(itRelease != idToRelease.end() )itRelease++;
            printf("in if\n");
        }
        else{
            it++;
            itB++;
            printf("in else\n");
        }
    }
    for(it = particles.begin(); it != particles.end(); it++){
        printf("%d\t",(*it)->getId() );
    }
    printf("\n");
}

template<typename T, template<typename U> class Descriptor>
void BulkParticleChunk2D<T,Descriptor>::rigidize(){
//after rigidize, when the particle experience a big torque, it will determine if
//to break off from the chunk, by comparing the external torque with internal torque
//given by adjacent particles.

//Therefore, rigidization involves identification of adjacent particles.

    T epsilon = 0.1;
    typename std::vector<BulkParticle2D<T,Descriptor>*>::iterator it = particles.begin(), it2;
    rigidBonds.clear();
    std::vector<int> tmpRigidBonds;
    for( ; it != particles.end(); it++){
        tmpRigidBonds.clear();
        for(it2 = particles.begin(); it2 != particles.end(); it2++){
            if( sqrt( pow( (*it)->getPosition()[0] - (*it2)->getPosition()[0], 2.0) +
                      pow( (*it)->getPosition()[1] - (*it2)->getPosition()[1], 2.0) )
                <= (*it)->getRadius() + (*it2)->getRadius() - epsilon ){
                tmpRigidBonds.push_back((*it2)->getId());
            }
        }
        rigidBonds.push_back(tmpRigidBonds);
    }

}


template<typename T, template<typename U> class Descriptor>
void BulkParticleChunk2D<T,Descriptor>::velocityToParticleCoupling (
        Box2D domain, TensorField2D<T,2>& velocityField, T scaling )
{
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (pluint iParticle=0; iParticle<particles.size(); ++iParticle) {
                    particles[iParticle]->velocityToParticle(velocityField, scaling);
                }
            }
        }
}



template<typename T, template<typename U> class Descriptor>
void BulkParticleChunk2D<T,Descriptor>::rhoBarJtoParticleCoupling (
        Box2D domain, NTensorField2D<T>& rhoBarJfield, bool velIsJ, T scaling )
{
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (pluint iParticle=0; iParticle<particles.size(); ++iParticle) {
                    particles[iParticle]->rhoBarJtoParticle(rhoBarJfield, velIsJ, scaling);
                }
            }
        }
}



template<typename T, template<typename U> class Descriptor>
void BulkParticleChunk2D<T,Descriptor>::fluidToParticleCoupling (
        Box2D domain, solidBlockLattice2D<T,Descriptor>& lattice )
{
    plint iX, iY;
    typename std::vector<BulkParticle2D<T,Descriptor>* >::iterator it = particles.begin();
    for( ; it != particles.end(); it++) {
        computeGridPosition((*it)->getPosition(), iX, iY);
        if ( contained(iX,iY, domain)) (*it)->fluidToParticle(lattice, 0.3);
    }
}

template<typename T, template<typename U> class Descriptor>
void BulkParticleChunk2D<T,Descriptor>::particleToFluidCoupling(BlockLattice2D<T,Descriptor>& lattice)
{
     typename std::vector<BulkParticle2D<T,Descriptor>* >::iterator it = particles.begin();
     for( ; it != particles.end(); it++) {
         /*(*it)->particleToFluid(lattice); */
         (*it)->particleToFluidSlip2(lattice, 0.6);
     }
}

template<typename T, template<typename U> class Descriptor>
void BulkParticleChunk2D<T,Descriptor>::chunkToFluid(BlockLattice2D<T,Descriptor>& fluid)
{
//    std::vector<Array<int,4> > LeftRightLink;
//    std::vector<Array<int,4> > UpDownLink;
//    std::vector<Array<int,4> > DiagonalLink;
//    std::vector<Array<int,4> > antiDiagonalLink;
    double tmp1, tmp2;
    double Mx=0, My=0, Zxx=0, Zxy=0, Zyy=0;
    double dx, dy;
    typename std::vector<BulkParticle2D<T,Descriptor>* >::iterator it = particles.begin();
    typename std::vector<std::vector<int> >::iterator itRigidBonds = rigidBonds.begin();
    for( ; it != particles.end(); it++) {
        itRigidBonds++;
        for(std::vector<Array<int,4> >::iterator i=(*it)->LeftRightLink.begin();   i!=(*it)->LeftRightLink.end(); ++i){
            tmp1=fluid.get( (*i)[0],(*i)[1] ).f[6];
            tmp2=fluid.get( (*i)[2],(*i)[3] ).f[2];
            fluid.get( (*i)[0],(*i)[1] ).f[6] = tmp2 + 2./3.*(*it)->getVelocity()[0];
            fluid.get( (*i)[2],(*i)[3] ).f[2] = tmp1 - 2./3.*(*it)->getVelocity()[0];
            Mx += 2*(tmp1-tmp2);
            Zxx += 2*2./3.;
        }
        for(std::vector<Array<int,4> >::iterator i=(*it)->UpDownLink.begin();      i!=(*it)->UpDownLink.end();    ++i){
            tmp1=fluid.get( (*i)[0],(*i)[1] ).f[8];
            tmp2=fluid.get( (*i)[2],(*i)[3] ).f[4];
            fluid.get( (*i)[0],(*i)[1] ).f[8] = tmp2 + 2./3.*((*it)->getVelocity()[1]);
            fluid.get( (*i)[2],(*i)[3] ).f[4] = tmp1 - 2./3.*((*it)->getVelocity()[1]);
            My += 2*(tmp1-tmp2);
            Zyy += 2*2./3.;
        }
        for(std::vector<Array<int,4> >::iterator i=(*it)->DiagonalLink.begin();    i!=(*it)->DiagonalLink.end();  ++i){
            tmp1=fluid.get( (*i)[0],(*i)[1] ).f[7];
            tmp2=fluid.get( (*i)[2],(*i)[3] ).f[3];
            fluid.get( (*i)[0],(*i)[1] ).f[7] = tmp2 + 1./6.*((*it)->getVelocity()[0] + (*it)->getVelocity()[1]);
            fluid.get( (*i)[2],(*i)[3] ).f[3] = tmp1 - 1./6.*((*it)->getVelocity()[0] + (*it)->getVelocity()[1]);
            Mx += 2*(tmp1-tmp2);
            My += 2*(tmp1-tmp2);
            Zxx += 2*1./6.;
            Zxy += 2*1./6.;
            Zyy += 2*1./6.;
        }
        for(std::vector<Array<int,4> >::iterator i=(*it)->antiDiagonalLink.begin();i!=(*it)->antiDiagonalLink.end(); ++i){
            tmp1=fluid.get( (*i)[0],(*i)[1] ).f[1];
            tmp2=fluid.get( (*i)[2],(*i)[3] ).f[5];
            fluid.get( (*i)[0],(*i)[1] ).f[1] = tmp2 + 1./6.*(-(*it)->getVelocity()[0] + (*it)->getVelocity()[1]);
            fluid.get( (*i)[2],(*i)[3] ).f[5] = tmp1 - 1./6.*(-(*it)->getVelocity()[0] + (*it)->getVelocity()[1]);
            Mx += -2*(tmp1-tmp2);
            My +=  2*(tmp1-tmp2);
            Zxx += 2*1./6.;
            Zxy += -2*1./6.;
            Zyy += 2*1./6.;
        }
    }
    double zxx, zxy, zyy;
    double alp = 1;
    double M = 0;
    for(it = particles.begin() ; it != particles.end(); it++) {
        M += 3.141592654*(*it)->getRadius()*(*it)->getRadius()*(*it)->getDensity();
    }
//This should be consistent with
    zxx = (M+alp*Zyy)/((M+alp*Zxx)*(M+alp*Zyy)-alp*alp*Zxy*Zxy);
    zxy =    alp*Zxy /((M+alp*Zxx)*(M+alp*Zyy)-alp*alp*Zxy*Zxy);
    zyy = (M+alp*Zxx)/((M+alp*Zxx)*(M+alp*Zyy)-alp*alp*Zxy*Zxy);


    for(it = particles.begin() ; it != particles.end(); it++) {
        double tmpX = (*it)->getVelocity()[0], tmpY = (*it)->getVelocity()[1];
        (*it)->getVelocity()[0] += zxx*(Mx-Zxx*tmpX-Zxy*tmpY)+zxy*(My-Zxy*tmpX-Zyy*tmpY);
        (*it)->getVelocity()[1] += zxy*(Mx-Zxx*tmpX-Zxy*tmpY)+zyy*(My-Zxy*tmpX-Zyy*tmpY);
    }
}

template<typename T, template<typename U> class Descriptor>
void BulkParticleChunk2D<T,Descriptor>::fixedChunkToFluid(BlockLattice2D<T,Descriptor>& fluid,
                                                          BulkParticleField2D<T,Descriptor>* particleField)
{
//    std::vector<Array<int,4> > LeftRightLink;
//    std::vector<Array<int,4> > UpDownLink;
//    std::vector<Array<int,4> > DiagonalLink;
//    std::vector<Array<int,4> > antiDiagonalLink;
    double tmp1, tmp2;
    double Mx=0, My=0, Zxx=0, Zxy=0, Zyy=0;
    double Torque, thresholdTorque = 0.05;
    double dx, dy;
    typename std::vector<BulkParticle2D<T,Descriptor>* >::iterator it = particles.begin();
    typename std::vector<std::vector<int> >::iterator itRigidBonds = rigidBonds.begin();
    std::vector<int> idToRelease;
    printf("in fixedChunkToFluid\n");
    for( ; it != particles.end(); it++) {
        printf("processing particle with id%d\n", (*it)->getId());
        Torque = 0;
        for(std::vector<Array<int,4> >::iterator i=(*it)->LeftRightLink.begin();   i!=(*it)->LeftRightLink.end(); ++i){
            tmp1=fluid.get( (*i)[0],(*i)[1] ).f[6];
            tmp2=fluid.get( (*i)[2],(*i)[3] ).f[2];
            fluid.get( (*i)[0],(*i)[1] ).f[6] = tmp2 + 2./3.*(*it)->getVelocity()[0];
            fluid.get( (*i)[2],(*i)[3] ).f[2] = tmp1 - 2./3.*(*it)->getVelocity()[0];
            Torque += 2*(tmp1-tmp2)*((*i)[1] - (*it)->getPosition()[1]);
        }
        for(std::vector<Array<int,4> >::iterator i=(*it)->UpDownLink.begin();      i!=(*it)->UpDownLink.end();    ++i){
            tmp1=fluid.get( (*i)[0],(*i)[1] ).f[8];
            tmp2=fluid.get( (*i)[2],(*i)[3] ).f[4];
            fluid.get( (*i)[0],(*i)[1] ).f[8] = tmp2 + 2./3.*((*it)->getVelocity()[1]);
            fluid.get( (*i)[2],(*i)[3] ).f[4] = tmp1 - 2./3.*((*it)->getVelocity()[1]);
            Torque -= 2*(tmp1-tmp2)*((*i)[0] - (*it)->getPosition()[0]);
        }
        for(std::vector<Array<int,4> >::iterator i=(*it)->DiagonalLink.begin();    i!=(*it)->DiagonalLink.end();  ++i){
            tmp1=fluid.get( (*i)[0],(*i)[1] ).f[7];
            tmp2=fluid.get( (*i)[2],(*i)[3] ).f[3];
            fluid.get( (*i)[0],(*i)[1] ).f[7] = tmp2 + 1./6.*((*it)->getVelocity()[0] + (*it)->getVelocity()[1]);
            fluid.get( (*i)[2],(*i)[3] ).f[3] = tmp1 - 1./6.*((*it)->getVelocity()[0] + (*it)->getVelocity()[1]);
            Torque += 2*(tmp1-tmp2)*cos(PI/4)*((*i)[1]-(*i)[0] - (*it)->getPosition()[1] + (*it)->getPosition()[0]);
        }
        for(std::vector<Array<int,4> >::iterator i=(*it)->antiDiagonalLink.begin();i!=(*it)->antiDiagonalLink.end(); ++i){
            tmp1=fluid.get( (*i)[0],(*i)[1] ).f[1];
            tmp2=fluid.get( (*i)[2],(*i)[3] ).f[5];
            fluid.get( (*i)[0],(*i)[1] ).f[1] = tmp2 + 1./6.*(-(*it)->getVelocity()[0] + (*it)->getVelocity()[1]);
            fluid.get( (*i)[2],(*i)[3] ).f[5] = tmp1 - 1./6.*(-(*it)->getVelocity()[0] + (*it)->getVelocity()[1]);
            Torque -= 2*(tmp1-tmp2)*cos(PI/4)*((*i)[0]+(*i)[1] - (*it)->getPosition()[0] - (*it)->getPosition()[1]);
        }
/*
        for(std::vector<int>::iterator i=(*itRigidBonds).begin(); i!=(*itRigidBonds).end(); i++){
            if(Torque > 0)Torque -= 0.5;
            else if(Torque < 0)Torque += 0.5;
        }*/
        printf("Torque = %f\n", Torque);
        if(Torque > thresholdTorque || Torque < -thresholdTorque){
            idToRelease.push_back((*it)->getId() );
            printf("ParticleToRelease added.\n");
        }
        itRigidBonds++;
    }

    printf("number of particles to release %d\n", idToRelease.size());
    for(int i=0;i<idToRelease.size();i++)printf("idToRelease=%d\n",idToRelease[i]);
    if(particles.size() > 1)releaseParticles(idToRelease, particleField);
    if(idToRelease.size() > 0)system("PAUSE");
    printf("exit fixedChunkToFluid\n");
//==============================================================================================
    //decide which particle should be released or which bond should be broken.
    //there are two ways to do this.
    //1) all the bonds associated with a particle are broken concurrently and the particle is released right away.
    //   This is easier to implement. Compare the external torque with the combined internal one.
    //2) the bonds break one by one, and the particle is not released until the last bond is broken.
    //   This will involve some chemical process.
    //We choose to implement 1) due to its simplicity.


}



template<typename T, template<typename U> class Descriptor>
void BulkParticleChunk2D<T,Descriptor>::advanceChunk(Box2D domain, T cutOffValue) {

    //PLB_ASSERT( contained(domain, particleGrid.getBoundingBox()) );
    for(int it = 0; it < particles.size(); it++){
        particles[it]->advance();
    }
    std::cout<<"Chunk advance."<<std::endl;
}


template<typename T, template<typename U> class Descriptor>
BulkParticleDataTransfer2D<T,Descriptor>& BulkParticleChunk2D<T,Descriptor>::getDataTransfer() {
}

template<typename T, template<typename U> class Descriptor>
BulkParticleDataTransfer2D<T,Descriptor> const& BulkParticleChunk2D<T,Descriptor>::getDataTransfer() const {
}

template<typename T, template<typename U> class Descriptor>
std::string BulkParticleChunk2D<T,Descriptor>::getBlockName() {
    return std::string("BulkParticleChunk2D");
}

template<typename T, template<typename U> class Descriptor>
std::string BulkParticleChunk2D<T,Descriptor>::basicType() {
    return std::string(NativeType<T>::getName());
}

template<typename T, template<typename U> class Descriptor>
std::string BulkParticleChunk2D<T,Descriptor>::descriptorType() {
    return std::string(Descriptor<T>::name);
}

template<typename T, template<typename U> class Descriptor>
void BulkParticleChunk2D<T,Descriptor>::computeMassCenter(){
    T weighted_positionX = 0, weighted_positionY = 0;
    T weight;
    T accu_weight=0;
    typename std::vector<BulkParticle2D<T,Descriptor>* >::iterator it = particles.begin();
    for(; it != particles.end(); it++){
        weight = 3.1415926535798932*(*it)->getDensity()*(*it)->getRadius()*(*it)->getRadius();
        accu_weight += weight;
        weighted_positionX += weight*(*it)->getPosition()[0];
        weighted_positionY += weight*(*it)->getPosition()[1];
    }
    weighted_positionX /= accu_weight;
    weighted_positionY /= accu_weight;
    massCenter[0] = weighted_positionX;
    massCenter[1] = weighted_positionY;
    return;
}

}  // namespace plb

#endif  // BULK_PARTICLE_FIELD_2D_HH
