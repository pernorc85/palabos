#ifndef BULK_PARTICLE_3D_HH
#define BULK_PARTICLE_3D_HH

#include "core/globalDefs.h"
#include "finiteDifference/interpolations3D.h"
#include "particles/particleIdentifiers3D.h"
#include "newBulkParticles/BulkParticle3D.h"
#include "newBulkParticles/solidBlockLattice3D.h"
#include "newBulkParticles/BulkParticleField3D.h"
#include <cmath>
#include <vector>
#define TRUE 1
#define FALSE 0

namespace plb {

/* *************** class BulkParticle2D ***************************************** */

template<typename T, template<typename U> class Descriptor>
BulkParticle3D<T,Descriptor>::BulkParticle3D()
    : Particle3D<T,Descriptor>(),
      AngularPosition(0),
      radius(1)
{ }

template<typename T, template<typename U> class Descriptor>
BulkParticle3D<T,Descriptor>::BulkParticle3D(plint tag_, Array<T,3> const& position_, T AngularPosition_,
                                                         Array<T,3> const& velocity_, T AngularVelocity_, T radius_, T density_)
    : Particle3D<T,Descriptor>(tag_, position_),
      AngularPosition(AngularPosition_),
      velocity(velocity_),
      AngularVelocity(AngularVelocity_),
      radius(radius_),
      density(density_)
{
}

template<typename T, template<typename U> class Descriptor>
void BulkParticle3D<T,Descriptor>::reset(Array<T,3> const& position_, T AngularPosition_)
{
    this->position = position_;
    AngularPosition = AngularPosition_;
}


template<typename T, template<typename U> class Descriptor>
void BulkParticle3D<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
{
    serializer.addValue(this->tag);
    serializer.addValues<T,3>(this->position);
    serializer.addValue(AngularPosition);
    serializer.addValue(radius);
}

template<typename T, template<typename U> class Descriptor>
void BulkParticle3D<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer)
{
    unserializer.readValue(this->tag);
    unserializer.readValues<T,3>(this->position);
    unserializer.readValue(AngularPosition);
    unserializer.readValue(radius);
}


template<typename T, template<typename U> class Descriptor>
BulkParticle3D<T,Descriptor>* BulkParticle3D<T,Descriptor>::clone() const {
    return new BulkParticle3D<T,Descriptor>(*this);
}
/*
template<typename T, template<typename U> class Descriptor>
bool Particle2D<T,Descriptor>::getVector(plint whichVector, Array<T,2>& vector) const {
    return false;
}

template<typename T, template<typename U> class Descriptor>
bool Particle2D<T,Descriptor>::getScalar(plint whichScalar, T& scalar) const {
    return false;
}

template<typename T, template<typename U> class Descriptor>
bool Particle2D<T,Descriptor>::getTensor(plint whichVector, Array<T,SymmetricTensorImpl<T,2>::n>& tensor) const
{
    return false;
}

template<typename T, template<typename U> class Descriptor>
bool Particle2D<T,Descriptor>::setScalars(std::vector<T> const& scalars)
{
    return false;
}

template<typename T, template<typename U> class Descriptor>
bool Particle2D<T,Descriptor>::setVectors(std::vector<Array<T,2> > const& vectors)
{
    return false;
}

template<typename T, template<typename U> class Descriptor>
bool Particle2D<T,Descriptor>::setTensors (
        std::vector<Array<T,SymmetricTensorImpl<T,2>::n> > const& tensors )
{
    return false;
}
*/

template<typename T, template<typename U> class Descriptor>
void BulkParticle3D<T,Descriptor>::rescale(int dxScale, int dtScale) {
    int dimDx = 1;
    int dimDt = 0;
    T scaleFactor = scaleFromReference(dxScale, dimDx, dtScale, dimDt);
    this->position *= scaleFactor;
    radius *= scaleFactor;
}



template<typename T, template<typename U> class Descriptor>
void BulkParticle3D<T,Descriptor>::computeBoundaryLinks(Box3D domain)
{
        LeftRightLink.clear();
        UpDownLink.clear();
        FrontBackLink.clear();
        DiagonalLink.clear();
        antiDiagonalLink.clear();
        DiagonalLink.clear();
        antiDiagonalLink.clear();
        DiagonalLink.clear();
        antiDiagonalLink.clear();
        double dx, dy, dz, R3;
        R3 = radius*radius*radius;
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                dx = iX - this->position[0];
                dy = iY - this->position[1];
                if( (dx*dx+dy*dy-R2)*((dx+1)*(dx+1)+dy*dy-R2)<=0 )
                    LeftRightLink.push_back(Array<int,4>(iX, iY, iX+1, iY));
                if( (dx*dx+dy*dy-R2)*(dx*dx+(dy+1)*(dy+1)-R2)<=0 )
                    UpDownLink.push_back(Array<int,4>(iX, iY, iX, iY+1));
                if( (dx*dx+dy*dy-R2)*((dx+1)*(dx+1)+(dy+1)*(dy+1)-R2)<=0 )
                    DiagonalLink.push_back(Array<int,4>(iX, iY, iX+1, iY+1));
                if( (dx*dx+(dy+1)*(dy+1)-R2)*((dx+1)*(dx+1)+dy*dy-R2)<=0 )
                    antiDiagonalLink.push_back(Array<int,4>(iX+1, iY, iX, iY+1));
            }
        }
}


template<typename T, template<typename U> class Descriptor>
void BulkParticle3D<T,Descriptor>::velocityToParticle(TensorField2D<T,3>& velocityField, T scaling)
{
    velocity = predictorCorrectorTensorField<T,3>(velocityField, this->getPosition(), scaling);
}

template<typename T, template<typename U> class Descriptor>
void BulkParticle3D<T,Descriptor>::rhoBarJtoParticle (
        NTensorField3D<T>& rhoBarJfield, bool velIsJ, T scaling )
{

}



template<typename T, template<typename U> class Descriptor>
void BulkParticle3D<T,Descriptor>::fluidToParticle(solidBlockLattice3D<T,Descriptor>& fluid, T scaling)
{
//    std::vector<Array<T,2> > fluidBoundary;
//    std::vector<Array<T,2> > solidBoundary;
//    std::vector<Array<T,2> > internalSolid;
    double solidFrac;
    Array<T,2> Force(0,0);
    T Torque=0;
    Array<T,2> cellForce;

    std::cout<<"fluidToParticle"<<std::endl;

    std::cout<<"fluidBoundary:"<<std::endl;
    for(std::vector<Array<int,2> >::iterator i=fluidBoundary.begin(); i!=fluidBoundary.end(); ++i){
        Array<int,2> cellPosition = (*i);
        computeSolidFraction(cellPosition, this->position, radius, solidFrac);
        T B = solidFrac;
        //std::cout<<"solidFrac="<<solidFrac<<std::endl;
        Array<T,Descriptor<T>::q> additional_collision;
        cellForce.resetToZero();
        if (contained(cellPosition[0],cellPosition[1], fluid.getBoundingBox())) {
            for(plint iPop=0; iPop<Descriptor<T>::q; iPop++){
                plint opposite_iPop;
                if(iPop == 0)opposite_iPop = 0;
                else if(iPop <= 4)opposite_iPop = iPop+Descriptor<T>::q/2;
                else opposite_iPop = iPop-Descriptor<T>::q/2;
                T rho = fluid.get(cellPosition[0],cellPosition[1]).computeDensity();
                Array<T,2> j = rho*fluid.getSolidCell(cellPosition[0],cellPosition[1]).getSolidVelocity();
                T jSqr = j[0]*j[0] + j[1]*j[1];
/*
omega_i = f_opposite_iPop(x,t) - feq_opposite_iPop(rho,u)
        + feq_iPop(rho,u_solid) - f_iPop(x,t);
*/
                additional_collision[iPop]  = fluid.get(cellPosition[0],cellPosition[1])[opposite_iPop];
                additional_collision[iPop] -= fluid.get(cellPosition[0],cellPosition[1]).computeEquilibrium(opposite_iPop, rho, j, jSqr, 0);
                additional_collision[iPop] -= fluid.get(cellPosition[0],cellPosition[1])[iPop];
                additional_collision[iPop] += fluid.get(cellPosition[0],cellPosition[1]).computeEquilibrium(iPop, rho, j, jSqr, 0);
                cellForce[0] += additional_collision[iPop]*Descriptor<T>::c[iPop][0];
                cellForce[1] += additional_collision[iPop]*Descriptor<T>::c[iPop][1];
            }
        }
        else {
            cellForce.resetToZero();
        }
        //std::cout<<"cellForce"<<cellForce[0]<<","<<cellForce[1]<<std::endl;
        Force[0] += cellForce[0]*scaling*B;
        Force[1] += cellForce[1]*scaling*B;
        Torque +=-(this->position[1]-cellPosition[1])*cellForce[0]*scaling*solidFrac
                 +(this->position[0]-cellPosition[0])*cellForce[1]*scaling*solidFrac;
    }

    system("PAUSE");
    std::cout<<"SolidBoundary:"<<std::endl;
    for(std::vector<Array<int,2> >::iterator i=solidBoundary.begin(); i!=solidBoundary.end(); ++i){
        Array<int,2> cellPosition = (*i);
        computeSolidFraction(cellPosition, this->position, radius, solidFrac);
        T B = solidFrac;
        //std::cout<<"solidFrac="<<solidFrac<<std::endl;
        Array<T,Descriptor<T>::q> additional_collision;
        cellForce.resetToZero();
        if (contained(cellPosition[0],cellPosition[1], fluid.getBoundingBox())) {
            for(plint iPop=0; iPop<Descriptor<T>::q; iPop++){
/*
omega_i = f_opposite_iPop(x,t) - feq_opposite_iPop(rho,u)
        + feq_iPop(rho,u_solid) - f_iPop(x,t);
*/
                plint opposite_iPop;
                if(iPop == 0)opposite_iPop = 0;
                else if(iPop <= 4)opposite_iPop = iPop+Descriptor<T>::q/2;
                else opposite_iPop = iPop-Descriptor<T>::q/2;
                T rho = fluid.get(cellPosition[0],cellPosition[1]).computeDensity();
                Array<T,2> j = rho*fluid.getSolidCell(cellPosition[0],cellPosition[1]).getSolidVelocity();
                T jSqr = j[0]*j[0] + j[1]*j[1];
//additional collision operator = (f_oppi - f_eq_oppi) - (f_i - f_eq_i);
                additional_collision[iPop]  = fluid.get(cellPosition[0],cellPosition[1])[opposite_iPop];
                additional_collision[iPop] -= fluid.get(cellPosition[0],cellPosition[1]).computeEquilibrium(opposite_iPop, rho, j, jSqr, 0);
                additional_collision[iPop] -= fluid.get(cellPosition[0],cellPosition[1])[iPop];
                additional_collision[iPop] += fluid.get(cellPosition[0],cellPosition[1]).computeEquilibrium(iPop, rho, j, jSqr, 0);
                cellForce[0] += additional_collision[iPop]*Descriptor<T>::c[iPop][0];
                cellForce[1] += additional_collision[iPop]*Descriptor<T>::c[iPop][1];
            }
        }
        else {
            cellForce.resetToZero();
        }
        //std::cout<<"cellForce"<<cellForce[0]<<","<<cellForce[1]<<std::endl;
        Force[0] += cellForce[0]*scaling*B;
        Force[1] += cellForce[1]*scaling*B;
        Torque +=-(this->position[1]-cellPosition[1])*cellForce[0]*scaling*solidFrac
                 +(this->position[0]-cellPosition[0])*cellForce[1]*scaling*solidFrac;
    }

    std::cout<<"internalSolid:"<<std::endl;
    for(std::vector<Array<int,2> >::iterator i=internalSolid.begin(); i!=internalSolid.end(); ++i){
        Array<int,2> cellPosition = (*i);
        //std::cout<<"solidFrac="<<solidFrac<<std::endl;
        Array<T,Descriptor<T>::q> additional_collision;
        cellForce.resetToZero();
        if (contained(cellPosition[0],cellPosition[1], fluid.getBoundingBox())) {
            for(plint iPop=0; iPop<Descriptor<T>::q; iPop++){
/*
omega_i = f_opposite_iPop(x,t) - feq_opposite_iPop(rho,u)
        + feq_iPop(rho,u_solid) - f_iPop(x,t);
*/
                plint opposite_iPop;
                if(iPop == 0)opposite_iPop = 0;
                else if(iPop <= 4)opposite_iPop = iPop+Descriptor<T>::q/2;
                else opposite_iPop = iPop-Descriptor<T>::q/2;
                T rho = fluid.get(cellPosition[0],cellPosition[1]).computeDensity();
                Array<T,2> j = rho*fluid.getSolidCell(cellPosition[0],cellPosition[1]).getSolidVelocity();
                T jSqr = j[0]*j[0] + j[1]*j[1];
//additional collision operator = (f_oppi - f_eq_oppi) - (f_i - f_eq_i);
                additional_collision[iPop]  = fluid.get(cellPosition[0],cellPosition[1])[opposite_iPop];
                additional_collision[iPop] -= fluid.get(cellPosition[0],cellPosition[1]).computeEquilibrium(opposite_iPop, rho, j, jSqr, 0);
                additional_collision[iPop] -= fluid.get(cellPosition[0],cellPosition[1])[iPop];
                additional_collision[iPop] += fluid.get(cellPosition[0],cellPosition[1]).computeEquilibrium(iPop, rho, j, jSqr, 0);
                cellForce[0] += additional_collision[iPop]*Descriptor<T>::c[iPop][0];
                cellForce[1] += additional_collision[iPop]*Descriptor<T>::c[iPop][1];
            }
        }
        else {
            cellForce.resetToZero();
        }
        //std::cout<<"cellForce"<<cellForce[0]<<","<<cellForce[1]<<std::endl;
        Force[0] += cellForce[0]*scaling;
        Force[1] += cellForce[1]*scaling;
        Torque +=-(this->position[1]-cellPosition[1])*cellForce[0]*scaling
                 +(this->position[0]-cellPosition[0])*cellForce[1]*scaling;
    }

    std::cout<<"Force"<<Force[0]<<","<<Force[1]<<std::endl;
    double mass = 3.141592654*getRadius()*getRadius()*getDensity();
    double moment_of_inertia = 2*3.141592654*getRadius()*getRadius()*getRadius()*getRadius()/4*getDensity();
    velocity[0] += -Force[0]/mass;
    velocity[1] += -Force[1]/mass;
    AngularVelocity += -Torque/moment_of_inertia;
    std::cout<<"After fluidToParticle acceleration, velocity="<<velocity[0]<<","<<velocity[1]<<std::endl;
    std::cout<<"AngularVelocity="<<AngularVelocity<<std::endl;

}

template<typename T, template<typename U> class Descriptor>
void BulkParticle3D<T,Descriptor>::particleToFluid(BlockLattice3D<T,Descriptor>& fluid)
{

//    std::vector<Array<int,4> > LeftRightLink;
//    std::vector<Array<int,4> > UpDownLink;
//    std::vector<Array<int,4> > DiagonalLink;
//    std::vector<Array<int,4> > antiDiagonalLink;
    double tmp1, tmp2;
    double Mx, My, Zxx, Zxy, Zyy;
    double dx, dy;
    Array<T,2> r;
    Array<T,2> additionalVelocity;
    for(std::vector<Array<int,6> >::iterator i=LeftRightLink.begin();   i!=LeftRightLink.end(); ++i){
        tmp1=fluid.get( (*i)[0],(*i)[1], (*i)[2] ).f[6];
        tmp2=fluid.get( (*i)[3],(*i)[4], (*i)[5] ).f[2];
        dx = (*i)[0] - this->position[0];
        dy = (*i)[1] - this->position[1];
        dz = (*i)[2] - this->position[2];
        if(dx*dx+dy*dy+dz*dz<=radius*radius*radius){
            r[0] = (T)((*i)[0]) - this->position[0];
            r[1] = (T)((*i)[1]) - this->position[1];
            r[2] = (T)((*i)[2]) - this->position[2];
        }
        else {
            r[0] = (T)((*i)[3]) - this->position[0];
            r[1] = (T)((*i)[4]) - this->position[1];
            r[2] = (T)((*i)[5]) - this->position[2];
        }
//        additionalVelocity[0] = r[1] * AngularVelocity;
//        additionalVelocity[1] =-r[0] * AngularVelocity;
    
        fluid.get( (*i)[0],(*i)[1] ).f[6] = tmp2 - 2./3.*(this->getVelocity()[0]+additionalVelocity[0]);
        fluid.get( (*i)[2],(*i)[3] ).f[2] = tmp1 + 2./3.*(this->getVelocity()[0]+additionalVelocity[0]);
        Mx += 2*(tmp2-tmp1);
        Zxx += 2*2./3.;
    }
    for(std::vector<Array<int,4> >::iterator i=UpDownLink.begin();      i!=UpDownLink.end();    ++i){
        tmp1=fluid.get( (*i)[0],(*i)[1] ).f[8];
        tmp2=fluid.get( (*i)[2],(*i)[3] ).f[4];
        dx = (*i)[0] - this->position[0];
        dy = (*i)[1] - this->position[1];
        if(dx*dx+dy*dy<=radius*radius){
            r[0] = (T)((*i)[0]) - this->position[0];
            r[1] = (T)((*i)[1]) - this->position[1];
        }
        else {
            r[0] = (T)((*i)[2]) - this->position[0];
            r[1] = (T)((*i)[3]) - this->position[1];
        }
        additionalVelocity[0] = r[1] * AngularVelocity;
        additionalVelocity[1] =-r[0] * AngularVelocity;
        fluid.get( (*i)[0],(*i)[1] ).f[8] = tmp2 - 2./3.*(this->getVelocity()[1]+additionalVelocity[1]);
        fluid.get( (*i)[2],(*i)[3] ).f[4] = tmp1 + 2./3.*(this->getVelocity()[1]+additionalVelocity[1]);
        My += 2*(tmp2-tmp1);
        Zyy += 2*2./3.;
    }
    for(std::vector<Array<int,4> >::iterator i=DiagonalLink.begin();    i!=DiagonalLink.end();  ++i){
        tmp1=fluid.get( (*i)[0],(*i)[1] ).f[7];
        tmp2=fluid.get( (*i)[2],(*i)[3] ).f[3];
        dx = (*i)[0] - this->position[0];
        dy = (*i)[1] - this->position[1];
        if(dx*dx+dy*dy<=radius*radius){
            r[0] = (T)((*i)[0]) - this->position[0];
            r[1] = (T)((*i)[1]) - this->position[1];
        }
        else {
            r[0] = (T)((*i)[2]) - this->position[0];
            r[1] = (T)((*i)[3]) - this->position[1];
        }
        additionalVelocity[0] = r[1] * AngularVelocity;
        additionalVelocity[1] =-r[0] * AngularVelocity;
        fluid.get( (*i)[0],(*i)[1] ).f[7] = tmp2 - 1./6.*(this->getVelocity()[0]+additionalVelocity[0] + this->getVelocity()[1]+additionalVelocity[1]);
        fluid.get( (*i)[2],(*i)[3] ).f[3] = tmp1 + 1./6.*(this->getVelocity()[0]+additionalVelocity[0] + this->getVelocity()[1]+additionalVelocity[1]);
        Mx += 2*(tmp2-tmp1);
        My += 2*(tmp2-tmp1);
        Zxx += 2*1./6.;
        Zxy += 2*1./6.;
        Zyy += 2*1./6.;
    }
    for(std::vector<Array<int,4> >::iterator i=antiDiagonalLink.begin();i!=antiDiagonalLink.end(); ++i){
        tmp1=fluid.get( (*i)[0],(*i)[1] ).f[1];
        tmp2=fluid.get( (*i)[2],(*i)[3] ).f[5];
        dx = (*i)[0] - this->position[0];
        dy = (*i)[1] - this->position[1];
        if(dx*dx+dy*dy<=radius*radius){
            r[0] = (T)((*i)[0]) - this->position[0];
            r[1] = (T)((*i)[1]) - this->position[1];
        }
        else {
            r[0] = (T)((*i)[2]) - this->position[0];
            r[1] = (T)((*i)[3]) - this->position[1];
        }
        additionalVelocity[0] = r[1] * AngularVelocity;
        additionalVelocity[1] =-r[0] * AngularVelocity;
        fluid.get( (*i)[0],(*i)[1] ).f[1] = tmp2 - 1./6.*(-this->getVelocity()[0]-additionalVelocity[0] + this->getVelocity()[1]+additionalVelocity[1]);
        fluid.get( (*i)[2],(*i)[3] ).f[5] = tmp1 + 1./6.*(-this->getVelocity()[0]-additionalVelocity[0] + this->getVelocity()[1]+additionalVelocity[1]);
        Mx += -2*(tmp2-tmp1);
        My +=  2*(tmp2-tmp1);
        Zxx += 2*1./6.;
        Zxy += -2*1./6.;
        Zyy += 2*1./6.;
    }
    double zxx, zxy, zyy;
    double alp = 1;
    double M = 3.141592654*getRadius()*getRadius()*getDensity();
    zxx = (M+alp*Zyy)/((M+alp*Zxx)*(M+alp*Zyy)+alp*alp*Zxy*Zxy);
    zxy =    alp*Zxy /((M+alp*Zxx)*(M+alp*Zyy)+alp*alp*Zxy*Zxy);
    zyy = (M+alp*Zxx)/((M+alp*Zxx)*(M+alp*Zyy)+alp*alp*Zxy*Zxy);

//    if(this->tag!=3){
        double tmpX = velocity[0], tmpY = velocity[1];
        velocity[0] += zxx*(Mx-Zxx*tmpX-Zxy*tmpY)+zxy*(My-Zxy*tmpX-Zyy*tmpY);
        velocity[1] += zxy*(Mx-Zxx*tmpX-Zxy*tmpY)+zyy*(My-Zxy*tmpX-Zyy*tmpY);
//    }

}

template<typename T, template<typename U> class Descriptor>
void BulkParticle2D<T,Descriptor>::particleToParticle(BulkParticleField2D<T,Descriptor>& bulkParticleField_)
{
    std::cout<<"particleToParticle:"<<std::endl;
    Array<T,2> Force;
    Force[0] = 0, Force[1] = 0;
    T epsilon = 10.0;

    for(pluint iParticle = 0; iParticle < bulkParticleField_.particles.size(); iParticle++){
        if(bulkParticleField_.particles[iParticle]->getTag() != this->getTag()){
            T distance = sqrt( pow(bulkParticleField_.particles[iParticle]->getPosition()[0] - this->getPosition()[0], 2.0) +
                               pow(bulkParticleField_.particles[iParticle]->getPosition()[1] - this->getPosition()[1], 2.0) );
            T theOtherRadius = bulkParticleField_.particles[iParticle]->getRadius();
            if(     distance > theOtherRadius + this->getRadius() + epsilon )
                Force += 0;
            else if(distance < theOtherRadius + this->getRadius() + epsilon &&
                    distance > theOtherRadius + this->getRadius() ){
                std::cout<<"approaching..."<<std::endl;
                Force[0] += 100*pow(distance - theOtherRadius - this->getRadius() - epsilon, 4.0) *
                           (this->getPosition()[0] - bulkParticleField_.particles[iParticle]->getPosition()[0])/distance;
                Force[1] += 100*pow(distance - theOtherRadius - this->getRadius() - epsilon, 4.0) *
                           (this->getPosition()[1] - bulkParticleField_.particles[iParticle]->getPosition()[1])/distance;
            }
            else if(distance < theOtherRadius + this->getRadius() && distance > 0 ){
                 std::cout<<"colliding..."<<std::endl;
                 Force[0] +=100*( pow(distance - theOtherRadius - this->getRadius() - epsilon, 4.0) +
                             (theOtherRadius + this->getRadius() - distance)/epsilon ) *
                           (this->getPosition()[0] - bulkParticleField_.particles[iParticle]->getPosition()[0])/distance;
                 Force[1] +=100*( pow(distance - theOtherRadius - this->getRadius() - epsilon, 4.0) +
                             (theOtherRadius + this->getRadius() - distance)/epsilon ) *
                           (this->getPosition()[1] - bulkParticleField_.particles[iParticle]->getPosition()[1])/distance;
            }
        }
    }
    std::cout<<"Force:"<<Force[0]<<","<<Force[1]<<std::endl;
    //system("PAUSE");
    double mass = 3.141592654*getRadius()*getRadius()*getDensity();
    velocity[0] += Force[0]/mass;
    velocity[1] += Force[1]/mass;
    std::cout<<"After particleToParticle acceleration, velocity="<<velocity[0]<<","<<velocity[1]<<std::endl;
}



template<typename T, template<typename U> class Descriptor>
void BulkParticle2D<T,Descriptor>::wallToParticle(solidBlockLattice2D<T,Descriptor>& fluid)
{
    std::cout<<"wallToParticle:"<<std::endl;
    Array<T,2> Force;
    Force[0] = 0;
    Force[1] = 0;
    T epsilon = 2.0;

    Array<T,2> theOtherPosition;//fictitious particle
    for(plint i=1; i<=4; i++){
        if(i==1){
            theOtherPosition[0] = -this->getPosition()[0];
            theOtherPosition[1] =  this->getPosition()[1];
        }
        else if(i==2){
            theOtherPosition[0] = 2*fluid.getNx() - this->getPosition()[0];
            theOtherPosition[1] =                   this->getPosition()[1];
        }
        else if(i==3){
            theOtherPosition[0] =  this->getPosition()[0];
            theOtherPosition[1] = -this->getPosition()[1];
        }
        else if(i==4){
            theOtherPosition[0] =                   this->getPosition()[0];
            theOtherPosition[1] = 2*fluid.getNy() - this->getPosition()[1];
        }

        T distance = sqrt( pow(theOtherPosition[0] - this->getPosition()[0],2.0) +
                           pow(theOtherPosition[1] - this->getPosition()[1],2.0) );

        if( distance > 2*this->getRadius() + epsilon )
            Force += 0;
        else if(distance < 2*this->getRadius() + epsilon && distance > 2*this->getRadius() ){
            std::cout<<"approaching wall..."<<std::endl;
            system("PAUSE");
            Force[0] += 200*pow(distance - 2*this->getRadius() - epsilon, 2.0)/pow(epsilon, 2.0) *
                               (this->getPosition()[0] - theOtherPosition[0])/distance;
            Force[1] += 200*pow(distance - 2*this->getRadius() - epsilon, 2.0)/pow(epsilon, 2.0) *
                               (this->getPosition()[1] - theOtherPosition[1])/distance;
        }
        else if(distance < 2*this->getRadius() && distance > 0){
            std::cout<<"colliding with wall..."<<std::endl;
            system("PAUSE");
            Force[0] += 200*( pow(distance - 2*this->getRadius() - epsilon, 2.0)/pow(epsilon, 2.0) +
                               (2*this->getRadius() - distance)/epsilon ) *
                               (this->getPosition()[0] - theOtherPosition[0])/distance;
            Force[1] += 200*( pow(distance - 2*this->getRadius() - epsilon, 2.0)/pow(epsilon, 2.0) +
                               (2*this->getRadius() - distance)/epsilon ) *
                               (this->getPosition()[1] - theOtherPosition[1])/distance;
        }
    }
    std::cout<<"Force:"<<Force[0]<<","<<Force[1]<<std::endl;
    double mass = 3.141592654*getRadius()*getRadius()*getDensity();
    velocity[0] += Force[0]/mass;
    velocity[1] += Force[1]/mass;
    std::cout<<"after wallToParticle:"<<velocity[0]<<","<<velocity[1]<<std::endl;
}




template<typename T, template<typename U> class Descriptor>
void BulkParticle2D<T,Descriptor>::advance() {
    PLB_ASSERT( norm(velocity)<1. );
    std::cout<<this->position[0]<<","<<this->position[1]<<std::endl;
    std::cout<<this->velocity[0]<<","<<this->velocity[1]<<std::endl;
    this->position[0] += velocity[0];
    this->position[1] += velocity[1];
    AngularPosition += AngularVelocity;
}
