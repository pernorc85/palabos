#ifndef BULK_PARTICLE_2D_HH
#define BULK_PARTICLE_2D_HH

#include "core/globalDefs.h"
#include "finiteDifference/interpolations2D.h"
#include "particles/particleIdentifiers2D.h"
#include "newBulkParticles/BulkParticle2D.h"
#include "newBulkParticles/solidBlockLattice2D.h"
#include "newBulkParticles/BulkParticleField2D.h"
#include <cmath>
#include <cstdio>
#include <vector>
#define TRUE 1
#define FALSE 0

namespace plb {

/* *************** class BulkParticle2D ***************************************** */

template<typename T, template<typename U> class Descriptor>
BulkParticle2D<T,Descriptor>::BulkParticle2D()
    : Particle2D<T,Descriptor>(),
      AngularPosition(0),
      radius(1)
{ }

template<typename T, template<typename U> class Descriptor>
BulkParticle2D<T,Descriptor>::BulkParticle2D(plint id_, plint tag_, Array<T,2> const& position_, T AngularPosition_,
                                                                    Array<T,2> const& velocity_, T AngularVelocity_, T radius_, T density_)
    : Particle2D<T,Descriptor>(tag_, position_),
      id(id_),
      AngularPosition(AngularPosition_),
      velocity(velocity_),
      AngularVelocity(AngularVelocity_),
      radius(radius_),
      density(density_)
{
}

template<typename T, template<typename U> class Descriptor>
void BulkParticle2D<T,Descriptor>::reset(Array<T,2> const& position_, T AngularPosition_)
{
    this->position = position_;
    AngularPosition = AngularPosition_;
}


template<typename T, template<typename U> class Descriptor>
void BulkParticle2D<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
{
    serializer.addValue(this->tag);
    serializer.addValues<T,2>(this->position);
    serializer.addValue(AngularPosition);
    serializer.addValue(radius);
}

template<typename T, template<typename U> class Descriptor>
void BulkParticle2D<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer)
{
    unserializer.readValue(this->tag);
    unserializer.readValues<T,2>(this->position);
    unserializer.readValue(AngularPosition);
    unserializer.readValue(radius);
}


template<typename T, template<typename U> class Descriptor>
BulkParticle2D<T,Descriptor>* BulkParticle2D<T,Descriptor>::clone() const {
    return new BulkParticle2D<T,Descriptor>(*this);
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
void BulkParticle2D<T,Descriptor>::rescale(int dxScale, int dtScale) {
    int dimDx = 1;
    int dimDt = 0;
    T scaleFactor = scaleFromReference(dxScale, dimDx, dtScale, dimDt);
    this->position *= scaleFactor;
    radius *= scaleFactor;
}

template<typename T, template<typename U> class Descriptor>
void BulkParticle2D<T,Descriptor>::computeBoundaries(Box2D domain) {
        solidBoundary.clear();
        internalSolid.clear();
        fluidBoundary.clear();
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if( inside(iX, iY)==TRUE &&
                                 (inside(iX-1,iY)==FALSE ||
                                  inside(iX+1,iY)==FALSE ||
                                  inside(iX,iY-1)==FALSE ||
                                  inside(iX,iY+1)==FALSE ||
                                  inside(iX-1,iY-1)==FALSE ||
                                  inside(iX-1,iY+1)==FALSE ||
                                  inside(iX+1,iY-1)==FALSE ||
                                  inside(iX+1,iY+1)==FALSE) ){
                    solidBoundary.push_back(Array<int,2>(iX, iY));
                }
                else if( inside(iX, iY)==TRUE) {
                    internalSolid.push_back(Array<int,2>(iX, iY));
                }
            }
        }
//        pcout<<"solidBoundary:\n";
//        for(plint it = 0; it < solidBoundary.size(); it++)
//            pcout<<solidBoundary[it][0]<<" "<<solidBoundary[it][1]<<"\n";
//        system("PAUSE");
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if( inside(iX, iY)==FALSE &&
                                (find(solidBoundary.begin(),solidBoundary.end(),Array<int,2>(iX-1,iY))!=solidBoundary.end() ||
                                 find(solidBoundary.begin(),solidBoundary.end(),Array<int,2>(iX+1,iY))!=solidBoundary.end() ||
                                 find(solidBoundary.begin(),solidBoundary.end(),Array<int,2>(iX,iY-1))!=solidBoundary.end() ||
                                 find(solidBoundary.begin(),solidBoundary.end(),Array<int,2>(iX,iY+1))!=solidBoundary.end() ||
                                 find(solidBoundary.begin(),solidBoundary.end(),Array<int,2>(iX-1,iY-1))!=solidBoundary.end() ||
                                 find(solidBoundary.begin(),solidBoundary.end(),Array<int,2>(iX+1,iY-1))!=solidBoundary.end() ||
                                 find(solidBoundary.begin(),solidBoundary.end(),Array<int,2>(iX-1,iY+1))!=solidBoundary.end() ||
                                 find(solidBoundary.begin(),solidBoundary.end(),Array<int,2>(iX+1,iY+1))!=solidBoundary.end()) ){
                    fluidBoundary.push_back(Array<T,2>(iX, iY));
                }
            }
        }
//        pcout<<"fluidBoundary:\n";
//        for(plint i = 0; i < fluidBoundary.size(); i++)
//            pcout<<fluidBoundary[i][0]<<" "<<fluidBoundary[i][1]<<"\n";

}

template<typename T, template<typename U> class Descriptor>
void BulkParticle2D<T,Descriptor>::computeBoundaryLinks(Box2D domain)
{
        LeftRightLink.clear();
        UpDownLink.clear();
        DiagonalLink.clear();
        antiDiagonalLink.clear();
        double dx, dy, R2;
        R2 = radius*radius;
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
void BulkParticle2D<T,Descriptor>::velocityToParticle(TensorField2D<T,2>& velocityField, T scaling)
{
    velocity = predictorCorrectorTensorField<T,2>(velocityField, this->getPosition(), scaling);
}

template<typename T, template<typename U> class Descriptor>
void BulkParticle2D<T,Descriptor>::rhoBarJtoParticle (
        NTensorField2D<T>& rhoBarJfield, bool velIsJ, T scaling )
{

}



template<typename T, template<typename U> class Descriptor>
void BulkParticle2D<T,Descriptor>::fluidToParticle(solidBlockLattice2D<T,Descriptor>& fluid, T scaling)
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
void BulkParticle2D<T,Descriptor>::particleToFluid(BlockLattice2D<T,Descriptor>& fluid)
{
//    std::vector<Array<int,4> > LeftRightLink;
//    std::vector<Array<int,4> > UpDownLink;
//    std::vector<Array<int,4> > DiagonalLink;
//    std::vector<Array<int,4> > antiDiagonalLink;
    double tmp1, tmp2;
    double Mx=0, My=0, Zxx=0, Zxy=0, Zyy=0;
    double dx, dy;
    Array<T,2> r;
    Array<T,2> additionalVelocity;
    for(std::vector<Array<int,4> >::iterator i=LeftRightLink.begin();   i!=LeftRightLink.end(); ++i){
        tmp1=fluid.get( (*i)[0],(*i)[1] ).f[6];
        tmp2=fluid.get( (*i)[2],(*i)[3] ).f[2];
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
        fluid.get( (*i)[0],(*i)[1] ).f[6] = tmp2 + 2./3.*(this->getVelocity()[0]+additionalVelocity[0]);
        fluid.get( (*i)[2],(*i)[3] ).f[2] = tmp1 - 2./3.*(this->getVelocity()[0]+additionalVelocity[0]);
        //tmp1=fluid.get( (*i)[0],(*i)[1] ).f[6];
        //tmp2=fluid.get( (*i)[2],(*i)[3] ).f[2];
        //==>
        //fluid.get( (*i)[0],(*i)[1] ).f[6] = fluid.get( (*i)[2],(*i)[3] ).f[2] + 2./3.*this->getVelocity()[0];
        //fluid.get( (*i)[2],(*i)[3] ).f[2] = fluid.get( (*i)[0],(*i)[1] ).f[6] - 2./3.*this->getVelocity()[0];
        //==>
        //fluid.get( (*i)[2],(*i)[3] ).f[6] = fluid.get( (*i)[2],(*i)[3] ).f[2] + 2./3.*this->getVelocity()[0];
        //fluid.get( (*i)[0],(*i)[1] ).f[6] = fluid.get( (*i)[2],(*i)[3] ).f[6];//and will propogate back next step
        //fluid.get( (*i)[0],(*i)[1] ).f[2] = fluid.get( (*i)[0],(*i)[1] ).f[6] - 2./3.*this->getVelocity()[0];
        //fluid.get( (*i)[2],(*i)[3] ).f[2] = fluid.get( (*i)[0],(*i)[1] ).f[2];//and will propagate back next step

        Mx += 2*(tmp1-tmp2);
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
        fluid.get( (*i)[0],(*i)[1] ).f[8] = tmp2 + 2./3.*(this->getVelocity()[1]+additionalVelocity[1]);
        fluid.get( (*i)[2],(*i)[3] ).f[4] = tmp1 - 2./3.*(this->getVelocity()[1]+additionalVelocity[1]);
        My += 2*(tmp1-tmp2);
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
        fluid.get( (*i)[0],(*i)[1] ).f[7] = tmp2 + 1./6.*(this->getVelocity()[0]+additionalVelocity[0] + this->getVelocity()[1]+additionalVelocity[1]);
        fluid.get( (*i)[2],(*i)[3] ).f[3] = tmp1 - 1./6.*(this->getVelocity()[0]+additionalVelocity[0] + this->getVelocity()[1]+additionalVelocity[1]);
        Mx += 2*(tmp1-tmp2);
        My += 2*(tmp1-tmp2);
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
        fluid.get( (*i)[0],(*i)[1] ).f[1] = tmp2 + 1./6.*(-this->getVelocity()[0]-additionalVelocity[0] + this->getVelocity()[1]+additionalVelocity[1]);
        fluid.get( (*i)[2],(*i)[3] ).f[5] = tmp1 - 1./6.*(-this->getVelocity()[0]-additionalVelocity[0] + this->getVelocity()[1]+additionalVelocity[1]);
        Mx += -2*(tmp1-tmp2);
        My +=  2*(tmp1-tmp2);
        Zxx += 2*1./6.;
        Zxy += -2*1./6.;
        Zyy += 2*1./6.;
    }
    double zxx, zxy, zyy;
    double alp = 1;
    double M = 3.141592654*getRadius()*getRadius()*getDensity();
//calculation of determinant is specious
//checking result: changing the sign almost has no effect
    zxx = (M+alp*Zyy)/((M+alp*Zxx)*(M+alp*Zyy)-alp*alp*Zxy*Zxy);
    zxy =    alp*Zxy /((M+alp*Zxx)*(M+alp*Zyy)-alp*alp*Zxy*Zxy);
    zyy = (M+alp*Zxx)/((M+alp*Zxx)*(M+alp*Zyy)-alp*alp*Zxy*Zxy);

//    if(this->tag!=2){
//        double tmpX = velocity[0], tmpY = velocity[1];
//        velocity[0] += zxx*(Mx-Zxx*tmpX-Zxy*tmpY)+zxy*(My-Zxy*tmpX-Zyy*tmpY);
//        velocity[1] += zxy*(Mx-Zxx*tmpX-Zxy*tmpY)+zyy*(My-Zxy*tmpX-Zyy*tmpY);
//    }
}


template<typename T, template<typename U> class Descriptor>
void BulkParticle2D<T,Descriptor>::particleToFluidSlip(BlockLattice2D<T,Descriptor>& fluid, T zeta)
{
//    std::vector<Array<int,4> > LeftRightLink;
//    std::vector<Array<int,4> > UpDownLink;
//    std::vector<Array<int,4> > DiagonalLink;
//    std::vector<Array<int,4> > antiDiagonalLink;
    double tmp1, tmp2;
    double Mx=0, My=0, Zxx=0, Zxy=0, Zyy=0;
    double dx, dy;
    printf("LeftRightLink adjustment.\n");
    for(std::vector<Array<int,4> >::iterator i=LeftRightLink.begin();   i!=LeftRightLink.end(); ++i){
        //right boundary
        if(pow((*i)[0] - this->position[0],2.0) + pow((*i)[1] - this->position[1],2.0) >= pow(this->radius,2.0) ){
            if(find(    DiagonalLink.begin(),     DiagonalLink.end(), Array<int,4>((*i)[0],(*i)[1],(*i)[0]+1, (*i)[1]+1)) !=     DiagonalLink.end() &&
               find(antiDiagonalLink.begin(), antiDiagonalLink.end(), Array<int,4>((*i)[0]+1,(*i)[1]-1,(*i)[0], (*i)[1])) != antiDiagonalLink.end() ){
                //right boundary is vertical, exchange 1<-->3
                printf("swap 1--3\n");
                //1  8|  7   full-slip : 1<-->7  3<-->5 if added in stream step, since during collide(revert), 1<-->5, 3<-->7. we only need to do 1<-->3
                //2   |  6   non-slip :  1<-->5  3<-->7 if added in stream step, since during collide(revert),...we do not need to do anything
                //3  4|  5
                tmp1=fluid.get( (*i)[0],(*i)[1] ).f[6];
                tmp2=fluid.get( (*i)[2],(*i)[3] ).f[2];
                fluid.get( (*i)[0], (*i)[1] ).f[1] = (1-zeta)*(fluid.get( (*i)[0], (*i)[1] ).f[5] + 1/6*this->getVelocity()[1]);
                fluid.get( (*i)[0], (*i)[1] ).f[1] += zeta*fluid.get( (*i)[0], (*i)[1] ).f[7];
                fluid.get( (*i)[0], (*i)[1] ).f[1] -= 1/6*this->getVelocity()[0];
                fluid.get( (*i)[0], (*i)[1] ).f[3] = (1-zeta)*(fluid.get( (*i)[0], (*i)[1] ).f[7] - 1/6*this->getVelocity()[1]);
                fluid.get( (*i)[0], (*i)[1] ).f[3] += zeta*fluid.get( (*i)[0], (*i)[1] ).f[5];
                fluid.get( (*i)[0], (*i)[1] ).f[3] -= 1/6*this->getVelocity()[0];

                fluid.get( (*i)[0], (*i)[1] ).f[2] = fluid.get( (*i)[0],(*i)[1] ).f[6] - 2./3.*this->getVelocity()[0];
                //swap(fluid.get( (*i)[0],(*i)[1] ).f[1], fluid.get( (*i)[0],(*i)[1] ).f[3]);
                fluid.get( (*i)[2], (*i)[3] ).f[5] = (1-zeta)*(fluid.get( (*i)[2], (*i)[3] ).f[1] - 1/6*this->getVelocity()[1]);
                fluid.get( (*i)[2], (*i)[3] ).f[5] += zeta*fluid.get( (*i)[2], (*i)[3] ).f[3];
                fluid.get( (*i)[2], (*i)[3] ).f[5] += 1/6*this->getVelocity()[0];
                fluid.get( (*i)[2], (*i)[3] ).f[7] = (1-zeta)*(fluid.get( (*i)[2], (*i)[3] ).f[3] + 1/6*this->getVelocity()[1]);
                fluid.get( (*i)[2], (*i)[3] ).f[7] += zeta*fluid.get( (*i)[2], (*i)[3] ).f[1];
                fluid.get( (*i)[2], (*i)[3] ).f[7] += 1/6*this->getVelocity()[0];

                fluid.get( (*i)[2], (*i)[3] ).f[6] = fluid.get( (*i)[2],(*i)[3] ).f[2] + 2./3.*this->getVelocity()[0];
                Mx += 2*(tmp1-tmp2);
                Zxx += 2*2./3.;
            }
            else if(find(    DiagonalLink.begin(),     DiagonalLink.end(), Array<int,4>((*i)[0],(*i)[1],(*i)[0]+1, (*i)[1]+1)) == DiagonalLink.end() &&
                    find(antiDiagonalLink.begin(), antiDiagonalLink.end(), Array<int,4>((*i)[0]+1,(*i)[1]-1,(*i)[0], (*i)[1])) != antiDiagonalLink.end() ){
                //right boundary is like "/", exchange 1<-->2
                printf("swap 1--2\n");
                swap(fluid.get( (*i)[0],(*i)[1] ).f[1], fluid.get( (*i)[0],(*i)[1] ).f[2]);
            }
            else if(find(    DiagonalLink.begin(),     DiagonalLink.end(), Array<int,4>((*i)[0],(*i)[1],(*i)[0]+1, (*i)[1]+1)) != DiagonalLink.end() &&
                    find(antiDiagonalLink.begin(), antiDiagonalLink.end(), Array<int,4>((*i)[0]+1,(*i)[1]-1,(*i)[0], (*i)[1])) == antiDiagonalLink.end() ){
                //right boundary is like "\", exchange 2<-->3
                printf("swap 2--3\n");
                swap(fluid.get( (*i)[0],(*i)[1] ).f[2], fluid.get( (*i)[0],(*i)[1] ).f[3]);
            }
            else{}
        }
        //left boundary
        if(pow((*i)[2] - this->position[0],2.0) + pow((*i)[3] - this->position[1],2.0) >= pow(this->radius,2.0) ){
            if(find(    DiagonalLink.begin(),     DiagonalLink.end(), Array<int,4>((*i)[2]-1,(*i)[3]-1,(*i)[2], (*i)[3])) !=     DiagonalLink.end() &&
               find(antiDiagonalLink.begin(), antiDiagonalLink.end(), Array<int,4>((*i)[2],(*i)[3],(*i)[2]-1, (*i)[3]+1)) != antiDiagonalLink.end() ){
                //left boundary is vertical, exchange 5<-->7
                printf("swap 5--7\n");
                tmp1=fluid.get( (*i)[0],(*i)[1] ).f[6];
                tmp2=fluid.get( (*i)[2],(*i)[3] ).f[2];
                fluid.get( (*i)[0], (*i)[1] ).f[1] = (1-zeta)*(fluid.get( (*i)[0], (*i)[1] ).f[5] + 1/6*this->getVelocity()[1]);
                fluid.get( (*i)[0], (*i)[1] ).f[1] += zeta*fluid.get( (*i)[0], (*i)[1] ).f[7];
                fluid.get( (*i)[0], (*i)[1] ).f[1] -= 1/6*this->getVelocity()[0];
                fluid.get( (*i)[0], (*i)[1] ).f[3] = (1-zeta)*(fluid.get( (*i)[0], (*i)[1] ).f[7] - 1/6*this->getVelocity()[1]);
                fluid.get( (*i)[0], (*i)[1] ).f[3] += zeta*fluid.get( (*i)[0], (*i)[1] ).f[5];
                fluid.get( (*i)[0], (*i)[1] ).f[3] -= 1/6*this->getVelocity()[0];

                fluid.get( (*i)[0], (*i)[1] ).f[2] = fluid.get( (*i)[0],(*i)[1] ).f[6] - 2./3.*this->getVelocity()[0];
                //swap(fluid.get( (*i)[0],(*i)[1] ).f[1], fluid.get( (*i)[0],(*i)[1] ).f[3]);
                fluid.get( (*i)[2], (*i)[3] ).f[5] = (1-zeta)*(fluid.get( (*i)[2], (*i)[3] ).f[1] - 1/6*this->getVelocity()[1]);
                fluid.get( (*i)[2], (*i)[3] ).f[5] += zeta*fluid.get( (*i)[2], (*i)[3] ).f[3];
                fluid.get( (*i)[2], (*i)[3] ).f[5] += 1/6*this->getVelocity()[0];
                fluid.get( (*i)[2], (*i)[3] ).f[7] = (1-zeta)*(fluid.get( (*i)[2], (*i)[3] ).f[3] + 1/6*this->getVelocity()[1]);
                fluid.get( (*i)[2], (*i)[3] ).f[7] += zeta*fluid.get( (*i)[2], (*i)[3] ).f[1];
                fluid.get( (*i)[2], (*i)[3] ).f[7] += 1/6*this->getVelocity()[0];

                fluid.get( (*i)[2], (*i)[3] ).f[6] = fluid.get( (*i)[2],(*i)[3] ).f[2] + 2./3.*this->getVelocity()[0];
                Mx += 2*(tmp1-tmp2);
                Zxx += 2*2./3.;
            }
            else if(find(    DiagonalLink.begin(),     DiagonalLink.end(), Array<int,4>((*i)[2]-1,(*i)[3]-1,(*i)[2], (*i)[3])) == DiagonalLink.end() &&
                find(antiDiagonalLink.begin(), antiDiagonalLink.end(), Array<int,4>((*i)[2],(*i)[3],(*i)[2]-1, (*i)[3]+1)) != antiDiagonalLink.end() ){
                //left boundary is like "/", exchange 5<-->6
                printf("swap 5--6\n");
                swap(fluid.get( (*i)[2],(*i)[3] ).f[5], fluid.get( (*i)[2],(*i)[3] ).f[6]);
            }
            else if(find(    DiagonalLink.begin(),     DiagonalLink.end(), Array<int,4>((*i)[2]-1,(*i)[3]-1,(*i)[2], (*i)[3])) != DiagonalLink.end() &&
                find(antiDiagonalLink.begin(), antiDiagonalLink.end(), Array<int,4>((*i)[2],(*i)[3],(*i)[2]-1, (*i)[3]+1)) == antiDiagonalLink.end() ){
                //left boundary is like "\", exchange 6<-->7
                printf("swap 6--7\n");
                swap(fluid.get( (*i)[2],(*i)[3] ).f[6], fluid.get( (*i)[2],(*i)[3] ).f[7]);
            }
            else{}
        }
    }
//==================================================================
    printf("UpDownLink adjustment.\n");
    for(std::vector<Array<int,4> >::iterator i=UpDownLink.begin();   i!=UpDownLink.end(); ++i){
        //bottom boundary
        if(pow((*i)[0] - this->position[0],2.0) + pow((*i)[1] - this->position[1],2.0) >= pow(this->radius,2.0) ){
            if(find(    DiagonalLink.begin(),     DiagonalLink.end(), Array<int,4>((*i)[2]-1,(*i)[3]-1,(*i)[2], (*i)[3])) !=     DiagonalLink.end() &&
               find(antiDiagonalLink.begin(), antiDiagonalLink.end(), Array<int,4>((*i)[2]+1,(*i)[3]-1,(*i)[2], (*i)[3])) != antiDiagonalLink.end() ){
                //bottom boundary is horizontal, exchange 1<-->7
                printf("swap 1--7\n");
                tmp1=fluid.get( (*i)[0],(*i)[1] ).f[8];
                tmp2=fluid.get( (*i)[2],(*i)[3] ).f[4];
                fluid.get( (*i)[2], (*i)[3] ).f[1] *= (1-zeta)*(fluid.get( (*i)[2], (*i)[3] ).f[5] - 1/6*this->getVelocity()[0]);
                fluid.get( (*i)[2], (*i)[3] ).f[1] +=     zeta*fluid.get( (*i)[2], (*i)[3] ).f[3];
                fluid.get( (*i)[2], (*i)[3] ).f[1] += 1/6*this->getVelocity()[1];
                fluid.get( (*i)[2], (*i)[3] ).f[7] *= (1-zeta)*(fluid.get( (*i)[2], (*i)[3] ).f[3] + 1/6*this->getVelocity()[0]);
                fluid.get( (*i)[2], (*i)[3] ).f[7] +=     zeta*fluid.get( (*i)[2], (*i)[3] ).f[5];
                fluid.get( (*i)[2], (*i)[3] ).f[7] += 1/6*this->getVelocity()[1];

                fluid.get( (*i)[2], (*i)[3] ).f[8] = fluid.get( (*i)[2], (*i)[3] ).f[4] + 2./3.*this->getVelocity()[1];
                //swap(fluid.get( (*i)[0],(*i)[1] ).f[1], fluid.get( (*i)[0],(*i)[1] ).f[7]);
                fluid.get( (*i)[0], (*i)[1] ).f[3] *= (1-zeta)*(fluid.get( (*i)[0], (*i)[1] ).f[7] - 1/6*this->getVelocity()[0]);
                fluid.get( (*i)[0], (*i)[1] ).f[3] +=     zeta*fluid.get( (*i)[0], (*i)[1] ).f[1];
                fluid.get( (*i)[0], (*i)[1] ).f[3] -= 1/6*this->getVelocity()[1];
                fluid.get( (*i)[0], (*i)[1] ).f[5] *= (1-zeta)*(fluid.get( (*i)[0], (*i)[1] ).f[1] + 1/6*this->getVelocity()[0]);
                fluid.get( (*i)[0], (*i)[1] ).f[5] +=     zeta*fluid.get( (*i)[0], (*i)[1] ).f[7];
                fluid.get( (*i)[0], (*i)[1] ).f[5] -= 1/6*this->getVelocity()[1];

                fluid.get( (*i)[0], (*i)[1] ).f[4] = fluid.get( (*i)[0], (*i)[1] ).f[8] - 2./3.*this->getVelocity()[1];
                My += 2*(tmp1-tmp2);
                Zyy += 2*2./3.;
            }
            else if(find(    DiagonalLink.begin(),     DiagonalLink.end(), Array<int,4>((*i)[0]-1,(*i)[1]-1,(*i)[0], (*i)[1])) == DiagonalLink.end() &&
                    find(antiDiagonalLink.begin(), antiDiagonalLink.end(), Array<int,4>((*i)[0]+1,(*i)[1]-1,(*i)[0], (*i)[1])) != antiDiagonalLink.end() ){
                //bottom boundary is like "", exchange 1<-->8
                printf("swap 1--8\n");
                swap(fluid.get( (*i)[0],(*i)[1] ).f[1], fluid.get( (*i)[0],(*i)[1] ).f[8]);
            }
            else if(find(    DiagonalLink.begin(),     DiagonalLink.end(), Array<int,4>((*i)[0]-1,(*i)[1]-1,(*i)[0], (*i)[1])) != DiagonalLink.end() &&
                    find(antiDiagonalLink.begin(), antiDiagonalLink.end(), Array<int,4>((*i)[0]+1,(*i)[1]-1,(*i)[0], (*i)[1])) == antiDiagonalLink.end() ){
                //right boundary is like "\", exchange 2<-->3
                printf("swap 7--8\n");
                swap(fluid.get( (*i)[0],(*i)[1] ).f[7], fluid.get( (*i)[0],(*i)[1] ).f[8]);
            }
            else{}
        }
        //boundary on top
        if(pow((*i)[2] - this->position[0],2.0) + pow((*i)[3] - this->position[1],2.0) >= pow(this->radius,2.0) ){
            if(find(    DiagonalLink.begin(),     DiagonalLink.end(), Array<int,4>((*i)[0],(*i)[1],(*i)[0]+1, (*i)[1]+1)) !=     DiagonalLink.end() &&
               find(antiDiagonalLink.begin(), antiDiagonalLink.end(), Array<int,4>((*i)[0],(*i)[1],(*i)[0]-1, (*i)[1]+1)) != antiDiagonalLink.end() ){
                printf("swap 3--5\n");
                tmp1=fluid.get( (*i)[0],(*i)[1] ).f[8];
                tmp2=fluid.get( (*i)[2],(*i)[3] ).f[4];
                fluid.get( (*i)[2], (*i)[3] ).f[1] *= (1-zeta)*(fluid.get( (*i)[2], (*i)[3] ).f[5] - 1/6*this->getVelocity()[0]);
                fluid.get( (*i)[2], (*i)[3] ).f[1] +=     zeta*fluid.get( (*i)[2], (*i)[3] ).f[3];
                fluid.get( (*i)[2], (*i)[3] ).f[1] += 1/6*this->getVelocity()[1];
                fluid.get( (*i)[2], (*i)[3] ).f[7] *= (1-zeta)*(fluid.get( (*i)[2], (*i)[3] ).f[3] + 1/6*this->getVelocity()[0]);
                fluid.get( (*i)[2], (*i)[3] ).f[7] +=     zeta*fluid.get( (*i)[2], (*i)[3] ).f[5];
                fluid.get( (*i)[2], (*i)[3] ).f[7] += 1/6*this->getVelocity()[1];

                fluid.get( (*i)[2], (*i)[3] ).f[8] = fluid.get( (*i)[2], (*i)[3] ).f[4] + 2./3.*this->getVelocity()[1];
                //swap(fluid.get( (*i)[0],(*i)[1] ).f[1], fluid.get( (*i)[0],(*i)[1] ).f[7]);
                fluid.get( (*i)[0], (*i)[1] ).f[3] *= (1-zeta)*(fluid.get( (*i)[0], (*i)[1] ).f[7] - 1/6*this->getVelocity()[0]);
                fluid.get( (*i)[0], (*i)[1] ).f[3] +=     zeta*fluid.get( (*i)[0], (*i)[1] ).f[1];
                fluid.get( (*i)[0], (*i)[1] ).f[3] -= 1/6*this->getVelocity()[1];
                fluid.get( (*i)[0], (*i)[1] ).f[5] *= (1-zeta)*(fluid.get( (*i)[0], (*i)[1] ).f[1] + 1/6*this->getVelocity()[0]);
                fluid.get( (*i)[0], (*i)[1] ).f[5] +=     zeta*fluid.get( (*i)[0], (*i)[1] ).f[7];
                fluid.get( (*i)[0], (*i)[1] ).f[5] -= 1/6*this->getVelocity()[1];

                fluid.get( (*i)[0], (*i)[1] ).f[4] = fluid.get( (*i)[0], (*i)[1] ).f[8] - 2./3.*this->getVelocity()[1];
                //swap(fluid.get( (*i)[2],(*i)[3] ).f[3], fluid.get( (*i)[2],(*i)[3] ).f[5]);
                My += 2*(tmp1-tmp2);
                Zyy += 2*2./3.;
            }
            else if(find(    DiagonalLink.begin(),     DiagonalLink.end(), Array<int,4>((*i)[2],(*i)[3],(*i)[2]+1, (*i)[3]+1)) == DiagonalLink.end() &&
                    find(antiDiagonalLink.begin(), antiDiagonalLink.end(), Array<int,4>((*i)[2],(*i)[3],(*i)[2]-1, (*i)[3]+1)) != antiDiagonalLink.end() ){
                printf("swap 4--5\n");
                swap(fluid.get( (*i)[2],(*i)[3] ).f[4], fluid.get( (*i)[2],(*i)[3] ).f[5]);
            }
            else if(find(    DiagonalLink.begin(),     DiagonalLink.end(), Array<int,4>((*i)[2],(*i)[3],(*i)[2]+1, (*i)[3]+1)) != DiagonalLink.end() &&
                    find(antiDiagonalLink.begin(), antiDiagonalLink.end(), Array<int,4>((*i)[2],(*i)[3],(*i)[2]-1, (*i)[3]+1)) == antiDiagonalLink.end() ){
                //left boundary is like "\", exchange 6<-->7
                printf("swap 3--4\n");
                swap(fluid.get( (*i)[2],(*i)[3] ).f[3], fluid.get( (*i)[2],(*i)[3] ).f[4]);
            }
            else{}
        }
    }
//=================================================================================
    printf("DiagonalLink adjustment.\n");
    for(std::vector<Array<int,4> >::iterator i=DiagonalLink.begin();   i!=DiagonalLink.end(); ++i){
        //up-right boundary
        tmp1=fluid.get( (*i)[0],(*i)[1] ).f[7];
        tmp2=fluid.get( (*i)[2],(*i)[3] ).f[3];
        if(pow((*i)[0] - this->position[0],2.0) + pow((*i)[1] - this->position[1],2.0) >= pow(this->radius,2.0) ){
            if(find(LeftRightLink.begin(), LeftRightLink.end(), Array<int,4>((*i)[0],(*i)[1],(*i)[0]+1, (*i)[1])) != LeftRightLink.end() &&
               find(UpDownLink.begin(),    UpDownLink.end(),    Array<int,4>((*i)[0],(*i)[1]+1,(*i)[0], (*i)[1])) != UpDownLink.end() ){
                printf("swap 2--4\n");
                fluid.get( (*i)[0], (*i)[1] ).f[2] = (1-zeta)*(fluid.get( (*i)[0], (*i)[1] ).f[6] - 2./3.*this->getVelocity()[0]);
                fluid.get( (*i)[0], (*i)[1] ).f[2] += zeta*(fluid.get( (*i)[0], (*i)[1] ).f[8] - 1/3*(this->getVelocity()[0]+this->getVelocity()[1]));

                fluid.get( (*i)[0], (*i)[1] ).f[4] = (1-zeta)*(fluid.get( (*i)[0], (*i)[1] ).f[8] = 2./3.*this->getVelocity()[1]);
                fluid.get( (*i)[0], (*i)[1] ).f[4] += zeta*(fluid.get( (*i)[0], (*i)[1] ).f[6] - 1/3*(this->getVelocity()[0]+this->getVelocity()[1]));

                fluid.get( (*i)[0], (*i)[1] ).f[3] = fluid.get( (*i)[0], (*i)[1] ).f[7] - 1/6*(this->getVelocity()[0]+this->getVelocity()[1]);

                fluid.get( (*i)[2], (*i)[3] ).f[7] = fluid.get( (*i)[2], (*i)[3] ).f[3] + 1/6*(this->getVelocity()[0]+this->getVelocity()[1]);
            }
        }
        //bottom-left boundary
        if(pow((*i)[2] - this->position[0],2.0) + pow((*i)[3] - this->position[1],2.0) >= pow(this->radius,2.0) ){
            if(find(LeftRightLink.begin(), LeftRightLink.end(), Array<int,4>((*i)[2]-1,(*i)[3],(*i)[2], (*i)[3])) != LeftRightLink.end() &&
               find(UpDownLink.begin(),    UpDownLink.end(),    Array<int,4>((*i)[2],(*i)[3],(*i)[2], (*i)[3]-1)) != UpDownLink.end() ){
                printf("swap 6--8\n");
                fluid.get( (*i)[2], (*i)[3] ).f[6] = (1-zeta)*(fluid.get( (*i)[2], (*i)[3] ).f[2] + 2./3.*this->getVelocity()[0]);
                fluid.get( (*i)[2], (*i)[3] ).f[6] += zeta*(fluid.get( (*i)[2], (*i)[3] ).f[4] + 1/3*(this->getVelocity()[0]+this->getVelocity()[1]));

                fluid.get( (*i)[2], (*i)[3] ).f[8] = (1-zeta)*(fluid.get( (*i)[2], (*i)[3] ).f[4] + 2./3.*this->getVelocity()[1]);
                fluid.get( (*i)[2], (*i)[3] ).f[8] += zeta*(fluid.get( (*i)[2], (*i)[3] ).f[2] + 1/3*(this->getVelocity()[0]+this->getVelocity()[1]));

                fluid.get( (*i)[2], (*i)[3] ).f[7] = fluid.get( (*i)[2], (*i)[3] ).f[3] + 1/6*(this->getVelocity()[0]+this->getVelocity()[1]);

                fluid.get( (*i)[0], (*i)[1] ).f[3] = fluid.get( (*i)[0], (*i)[1] ).f[7] - 1/6*(this->getVelocity()[0]+this->getVelocity()[1]);
            }
        }
        Mx += 2*(tmp1-tmp2);
        My += 2*(tmp1-tmp2);
        Zxx += 2*1./6.;
        Zxy += 2*1./6.;
        Zyy += 2*1./6.;
    }
//=================================================================================
    printf("antiDiagonalLink adjustment.\n");
    for(std::vector<Array<int,4> >::iterator i=antiDiagonalLink.begin();   i!=antiDiagonalLink.end(); ++i){
        //up-left boundary
        tmp1=fluid.get( (*i)[0],(*i)[1] ).f[1];
        tmp2=fluid.get( (*i)[2],(*i)[3] ).f[5];
        if(pow((*i)[0] - this->position[0],2.0) + pow((*i)[1] - this->position[1],2.0) >= pow(this->radius,2.0) ){
            if(find(LeftRightLink.begin(), LeftRightLink.end(), Array<int,4>((*i)[0]-1,(*i)[1],(*i)[0], (*i)[1])) != LeftRightLink.end() &&
               find(UpDownLink.begin(),    UpDownLink.end(),    Array<int,4>((*i)[0],(*i)[1]+1,(*i)[0], (*i)[1])) != UpDownLink.end() ){
                printf("swap 4--6\n");
                fluid.get( (*i)[0], (*i)[1] ).f[4] = (1-zeta)*(fluid.get( (*i)[0], (*i)[1] ).f[8] - 2./3.*this->getVelocity()[1]);
                fluid.get( (*i)[0], (*i)[1] ).f[4] += zeta*(fluid.get( (*i)[0], (*i)[1] ).f[2] + 1/3*(this->getVelocity()[0]-this->getVelocity()[1]));

                fluid.get( (*i)[0], (*i)[1] ).f[6] = (1-zeta)*(fluid.get( (*i)[0], (*i)[1] ).f[2] + 2./3.*this->getVelocity()[0]);
                fluid.get( (*i)[0], (*i)[1] ).f[6] += zeta*(fluid.get( (*i)[0], (*i)[1] ).f[8] + 1/3*(this->getVelocity()[0]-this->getVelocity()[1]));

                fluid.get( (*i)[0], (*i)[1] ).f[5] = fluid.get( (*i)[0], (*i)[1] ).f[1] + 1/6*(this->getVelocity()[0]-this->getVelocity()[1]);

                fluid.get( (*i)[2], (*i)[3] ).f[1] = fluid.get( (*i)[2], (*i)[3] ).f[5] - 1/6*(this->getVelocity()[0]-this->getVelocity()[1]);
            }
        }
        //bottom-right boundary
        if(pow((*i)[2] - this->position[0],2.0) + pow((*i)[3] - this->position[1],2.0) >= pow(this->radius,2.0) ){
            if(find(LeftRightLink.begin(), LeftRightLink.end(), Array<int,4>((*i)[2],(*i)[3],(*i)[2]+1, (*i)[3])) != LeftRightLink.end() &&
               find(UpDownLink.begin(),    UpDownLink.end(),    Array<int,4>((*i)[2],(*i)[3],(*i)[2], (*i)[3]-1)) != UpDownLink.end() ){
                printf("swap 2--8\n");
                fluid.get( (*i)[2], (*i)[3] ).f[2] = (1-zeta)*(fluid.get( (*i)[2], (*i)[3] ).f[6] - 2./3.*this->getVelocity()[1]);
                fluid.get( (*i)[2], (*i)[3] ).f[2] += zeta*(fluid.get( (*i)[2], (*i)[3] ).f[4] - 1/3*(this->getVelocity()[0]-this->getVelocity()[1]));

                fluid.get( (*i)[2], (*i)[3] ).f[8] = (1-zeta)*(fluid.get( (*i)[2], (*i)[3] ).f[4] + 2./3.*this->getVelocity()[1]);
                fluid.get( (*i)[2], (*i)[3] ).f[8] += zeta*(fluid.get( (*i)[2], (*i)[3] ).f[6] - 1/3*(this->getVelocity()[0]-this->getVelocity()[1]));

                fluid.get( (*i)[2], (*i)[3] ).f[1] = fluid.get( (*i)[2], (*i)[3] ).f[5] - 1/6*(this->getVelocity()[0] - this->getVelocity()[1]);

                fluid.get( (*i)[0], (*i)[1] ).f[5] = fluid.get( (*i)[0], (*i)[1] ).f[1] + 1/6*(this->getVelocity()[0] - this->getVelocity()[1]);
            }
        }
        Mx += -2*(tmp1-tmp2);
        My +=  2*(tmp1-tmp2);
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

    if(this->tag!=2){
        double tmpX = velocity[0], tmpY = velocity[1];
        velocity[0] += zxx*(Mx-Zxx*tmpX-Zxy*tmpY)+zxy*(My-Zxy*tmpX-Zyy*tmpY);
        velocity[1] += zxy*(Mx-Zxx*tmpX-Zxy*tmpY)+zyy*(My-Zxy*tmpX-Zyy*tmpY);
    }
    return;
}

template<typename T, template<typename U> class Descriptor>
void BulkParticle2D<T,Descriptor>::particleToFluidSlip2(BlockLattice2D<T,Descriptor>& fluid, T zeta){
    zeta = 0.8;
    double tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
    double Mx=0, My=0, Zxx=0, Zxy=0, Zyy=0;
    double dx, dy;
    Array<T,2> r;
    Array<T,2> additionalVelocity, tangent_velocity, normal_velocity;
    T alpha;
    Array<T,2> slip_tangent_velocity;
    for(std::vector<Array<int,2> >::iterator i=this->solidBoundary.begin();   i!=this->solidBoundary.end(); ++i){
        //calculate tangent velocity of solid boundary.
        tmp1=fluid.get( (*i)[0],(*i)[1] ).f[1];
        tmp2=fluid.get( (*i)[0],(*i)[1] ).f[2];
        tmp3=fluid.get( (*i)[0],(*i)[1] ).f[3];
        tmp4=fluid.get( (*i)[0],(*i)[1] ).f[4];
        tmp5=fluid.get( (*i)[0],(*i)[1] ).f[5];
        tmp6=fluid.get( (*i)[0],(*i)[1] ).f[6];
        tmp7=fluid.get( (*i)[0],(*i)[1] ).f[7];
        tmp8=fluid.get( (*i)[0],(*i)[1] ).f[8];
        r[0] = (T)((*i)[0]) - this->position[0];
        r[1] = (T)((*i)[1]) - this->position[1];

        //calculate tangent velocity of boundary cells after slip
        if(r[0] > 0.0 )alpha = atan(r[1]/r[0]);
        else if(r[0] < 0.0 && r[1] >= 0.0)alpha = atan(r[1]/r[0]) + PI;
        else if(r[0] < 0.0 && r[1] < 0.0)alpha = atan(r[1]/r[0]) - PI;
        else if(r[0] == 0.0 && r[1] > 0.0)alpha = PI/2;
        else if(r[0] == 0.0 && r[1] < 0.0)alpha = -PI/2;
        else{printf("error"); system("PAUSE");}
        tangent_velocity[0]=(-velocity[0]*sin(alpha) + velocity[1]*cos(alpha) )*(-sin(alpha));
        tangent_velocity[1]=(-velocity[0]*sin(alpha) + velocity[1]*cos(alpha) )*cos(alpha);
        normal_velocity[0] =( velocity[0]*cos(alpha) + velocity[1]*sin(alpha) )*cos(alpha);
        normal_velocity[1] =( velocity[0]*cos(alpha) + velocity[1]*sin(alpha) )*sin(alpha);
        slip_tangent_velocity[0] = tangent_velocity[0] - 3*(1-zeta)/zeta*(
                                            (fluid.get( (*i)[0],(*i)[1] ).f[3]
                                            +fluid.get( (*i)[0],(*i)[1] ).f[7]
                                            -fluid.get( (*i)[0],(*i)[1] ).f[1]
                                            -fluid.get( (*i)[0],(*i)[1] ).f[5])*cos(2.0*alpha)
                                                            +
                                            (fluid.get( (*i)[0],(*i)[1] ).f[2]
                                            +fluid.get( (*i)[0],(*i)[1] ).f[6]
                                            -fluid.get( (*i)[0],(*i)[1] ).f[4]
                                            -fluid.get( (*i)[0],(*i)[1] ).f[8])*sin(2.0*alpha)*0.5
                                                            )*(-sin(alpha));
        printf("tangent_velocity[0] = %f\n", tangent_velocity[0]);
        printf("slip_tangent_velocity[0] = %f\n", slip_tangent_velocity[0]);
        slip_tangent_velocity[1] = tangent_velocity[1] - 3*(1-zeta)/zeta*(
                                            (fluid.get( (*i)[0],(*i)[1] ).f[3]
                                            +fluid.get( (*i)[0],(*i)[1] ).f[7]
                                            -fluid.get( (*i)[0],(*i)[1] ).f[1]
                                            -fluid.get( (*i)[0],(*i)[1] ).f[5])*cos(2.0*alpha)
                                                            +
                                            (fluid.get( (*i)[0],(*i)[1] ).f[2]
                                            +fluid.get( (*i)[0],(*i)[1] ).f[6]
                                            -fluid.get( (*i)[0],(*i)[1] ).f[4]
                                            -fluid.get( (*i)[0],(*i)[1] ).f[8])*sin(2.0*alpha)*0.5
                                                            )*(cos(alpha));
        printf("tangent_velocity[1] = %f\n", tangent_velocity[1]);
        printf("slip_tangent_velocity[1] = %f\n", slip_tangent_velocity[1]);
        printf("slip 0/1 = %f\t", slip_tangent_velocity[0]/slip_tangent_velocity[1]);
        printf("-tan(alpha) = %f\n", -tan(alpha));
        iniCellAtEquilibrium(fluid.get( (*i)[0],(*i)[1] ), 1.0, slip_tangent_velocity+normal_velocity);
        //calculate the force from fluid back to Object

        Mx -= 2*(fluid.get( (*i)[0],(*i)[1] ).f[7] - tmp7
                +fluid.get( (*i)[0],(*i)[1] ).f[6] - tmp6
                +fluid.get( (*i)[0],(*i)[1] ).f[5] - tmp5
                -fluid.get( (*i)[0],(*i)[1] ).f[1] + tmp1
                -fluid.get( (*i)[0],(*i)[1] ).f[2] + tmp2
                -fluid.get( (*i)[0],(*i)[1] ).f[3] + tmp3);
        My -= 2*(fluid.get( (*i)[0],(*i)[1] ).f[1] - tmp1
                +fluid.get( (*i)[0],(*i)[1] ).f[7] - tmp7
                +fluid.get( (*i)[0],(*i)[1] ).f[8] - tmp8
                -fluid.get( (*i)[0],(*i)[1] ).f[3] + tmp3
                -fluid.get( (*i)[0],(*i)[1] ).f[4] + tmp4
                -fluid.get( (*i)[0],(*i)[1] ).f[5] + tmp5);
/*
        if(tmp2 != tmp6){
            Zxx += 2*2./3.*abs(fluid.get( (*i)[0],(*i)[1] ).f[6] - fluid.get( (*i)[0],(*i)[1] ).f[2])/abs(tmp6-tmp2);
            printf("difference ratio %f\n", (fluid.get( (*i)[0],(*i)[1] ).f[6] - fluid.get( (*i)[0],(*i)[1] ).f[2])/(tmp6-tmp2));
        }
        if(tmp4 != tmp8){
            Zyy += 2*2./3.*abs(fluid.get( (*i)[0],(*i)[1] ).f[8] - fluid.get( (*i)[0],(*i)[1] ).f[4])/abs(tmp8-tmp4);
            printf("difference ratio %f\n", (fluid.get( (*i)[0],(*i)[1] ).f[8] - fluid.get( (*i)[0],(*i)[1] ).f[4])/(tmp8-tmp4));
        }
        if(tmp1 != tmp5){
            Zxx += 2*1./6.*abs(fluid.get( (*i)[0],(*i)[1] ).f[5] - fluid.get( (*i)[0],(*i)[1] ).f[1])/abs(tmp5-tmp1);
            Zyy += 2*1./6.*abs(fluid.get( (*i)[0],(*i)[1] ).f[1] - fluid.get( (*i)[0],(*i)[1] ).f[5])/abs(tmp1-tmp5);
            Zxy +=-2*1./6.*abs(fluid.get( (*i)[0],(*i)[1] ).f[1] - fluid.get( (*i)[0],(*i)[1] ).f[5])/abs(tmp1-tmp5);
            printf("difference ratio %f\n", (fluid.get( (*i)[0],(*i)[1] ).f[5] - fluid.get( (*i)[0],(*i)[1] ).f[1])/(tmp5-tmp1));
        }
        if(tmp3 != tmp7){
            Zxx += 2*1./6.*abs(fluid.get( (*i)[0],(*i)[1] ).f[7] - fluid.get( (*i)[0],(*i)[1] ).f[3])/abs(tmp7-tmp3);
            Zyy += 2*1./6.*abs(fluid.get( (*i)[0],(*i)[1] ).f[7] - fluid.get( (*i)[0],(*i)[1] ).f[3])/abs(tmp7-tmp3);
            Zxy += 2*1./6.*abs(fluid.get( (*i)[0],(*i)[1] ).f[7] - fluid.get( (*i)[0],(*i)[1] ).f[3])/abs(tmp7-tmp3);
            printf("difference ratio %f\n", (fluid.get( (*i)[0],(*i)[1] ).f[7] - fluid.get( (*i)[0],(*i)[1] ).f[3])/(tmp7-tmp3));
        }
*/
    }

    for(std::vector<Array<int,2> >::iterator i=this->fluidBoundary.begin();   i!=this->fluidBoundary.end(); ++i){
        //calculate tangent velocity of solid boundary.
        tmp1=fluid.get( (*i)[0],(*i)[1] ).f[1];
        tmp2=fluid.get( (*i)[0],(*i)[1] ).f[2];
        tmp3=fluid.get( (*i)[0],(*i)[1] ).f[3];
        tmp4=fluid.get( (*i)[0],(*i)[1] ).f[4];
        tmp5=fluid.get( (*i)[0],(*i)[1] ).f[5];
        tmp6=fluid.get( (*i)[0],(*i)[1] ).f[6];
        tmp7=fluid.get( (*i)[0],(*i)[1] ).f[7];
        tmp8=fluid.get( (*i)[0],(*i)[1] ).f[8];
        r[0] = (T)((*i)[0]) - this->position[0];
        r[1] = (T)((*i)[1]) - this->position[1];
        //calculate tangent velocity of boundary cells after slip
        if(r[0] > 0.0 )alpha = atan(r[1]/r[0]);
        else if(r[0] < 0.0 && r[1] >= 0.0)alpha = atan(r[1]/r[0]) + PI;
        else if(r[0] < 0.0 && r[1] < 0.0)alpha = atan(r[1]/r[0]) - PI;
        else if(r[0] == 0.0 && r[1] > 0.0)alpha = PI/2;
        else if(r[0] == 0.0 && r[1] < 0.0)alpha = -PI/2;
        else{printf("error"); system("PAUSE");}
        tangent_velocity[0]=(-velocity[0]*sin(alpha) + velocity[1]*cos(alpha) )*sin(-alpha);
        tangent_velocity[1]=(-velocity[0]*sin(alpha) + velocity[1]*cos(alpha) )*cos(alpha);
        normal_velocity[0] =( velocity[0]*cos(alpha) + velocity[1]*sin(alpha) )*cos(alpha);
        normal_velocity[1] =( velocity[0]*cos(alpha) + velocity[1]*sin(alpha) )*sin(alpha);
        slip_tangent_velocity[0] = tangent_velocity[0] - 3*(1-zeta)/zeta*(
                                            (fluid.get( (*i)[0],(*i)[1] ).f[3]
                                            +fluid.get( (*i)[0],(*i)[1] ).f[7]
                                            -fluid.get( (*i)[0],(*i)[1] ).f[1]
                                            -fluid.get( (*i)[0],(*i)[1] ).f[5])*cos(2.0*alpha)
                                                            +
                                            (fluid.get( (*i)[0],(*i)[1] ).f[2]
                                            +fluid.get( (*i)[0],(*i)[1] ).f[6]
                                            -fluid.get( (*i)[0],(*i)[1] ).f[4]
                                            -fluid.get( (*i)[0],(*i)[1] ).f[8])*sin(2.0*alpha)*0.5
                                                            )*(-sin(alpha));
        slip_tangent_velocity[1] = tangent_velocity[1] - 3*(1-zeta)/zeta*(
                                            (fluid.get( (*i)[0],(*i)[1] ).f[3]
                                            +fluid.get( (*i)[0],(*i)[1] ).f[7]
                                            -fluid.get( (*i)[0],(*i)[1] ).f[1]
                                            -fluid.get( (*i)[0],(*i)[1] ).f[5])*cos(2.0*alpha)
                                                            +
                                            (fluid.get( (*i)[0],(*i)[1] ).f[2]
                                            +fluid.get( (*i)[0],(*i)[1] ).f[6]
                                            -fluid.get( (*i)[0],(*i)[1] ).f[4]
                                            -fluid.get( (*i)[0],(*i)[1] ).f[8])*sin(2.0*alpha)*0.5
                                                            )*(cos(alpha));
        iniCellAtEquilibrium(fluid.get( (*i)[0],(*i)[1] ), 1.0, slip_tangent_velocity+normal_velocity);
        //calculate the force from fluid back to Object
        Mx -= 2*(fluid.get( (*i)[0],(*i)[1] ).f[7] - tmp7
                +fluid.get( (*i)[0],(*i)[1] ).f[6] - tmp6
                +fluid.get( (*i)[0],(*i)[1] ).f[5] - tmp5
                -fluid.get( (*i)[0],(*i)[1] ).f[1] + tmp1
                -fluid.get( (*i)[0],(*i)[1] ).f[2] + tmp2
                -fluid.get( (*i)[0],(*i)[1] ).f[3] + tmp3);
        My -= 2*(fluid.get( (*i)[0],(*i)[1] ).f[1] - tmp1
                +fluid.get( (*i)[0],(*i)[1] ).f[7] - tmp7
                +fluid.get( (*i)[0],(*i)[1] ).f[8] - tmp8
                -fluid.get( (*i)[0],(*i)[1] ).f[3] + tmp3
                -fluid.get( (*i)[0],(*i)[1] ).f[4] + tmp4
                -fluid.get( (*i)[0],(*i)[1] ).f[5] + tmp5);
/*
        if(tmp2 != tmp6){
            Zxx += 2*2./3.*abs(fluid.get( (*i)[0],(*i)[1] ).f[6] - fluid.get( (*i)[0],(*i)[1] ).f[2])/abs(tmp6-tmp2);
        }
        if(tmp4 != tmp8){
            Zyy += 2*2./3.*abs(fluid.get( (*i)[0],(*i)[1] ).f[8] - fluid.get( (*i)[0],(*i)[1] ).f[4])/abs(tmp8-tmp4);
        }
        if(tmp1 != tmp5){
            Zxx += 2*1./6.*abs(fluid.get( (*i)[0],(*i)[1] ).f[5] - fluid.get( (*i)[0],(*i)[1] ).f[1])/abs(tmp5-tmp1);
            Zyy += 2*1./6.*abs(fluid.get( (*i)[0],(*i)[1] ).f[1] - fluid.get( (*i)[0],(*i)[1] ).f[5])/abs(tmp1-tmp5);
            Zxy +=-2*1./6.*abs(fluid.get( (*i)[0],(*i)[1] ).f[1] - fluid.get( (*i)[0],(*i)[1] ).f[5])/abs(tmp1-tmp5);
        }
        if(tmp3 != tmp7){
            Zxx += 2*1./6.*abs(fluid.get( (*i)[0],(*i)[1] ).f[7] - fluid.get( (*i)[0],(*i)[1] ).f[3])/abs(tmp7-tmp3);
            Zyy += 2*1./6.*abs(fluid.get( (*i)[0],(*i)[1] ).f[7] - fluid.get( (*i)[0],(*i)[1] ).f[3])/abs(tmp7-tmp3);
            Zxy += 2*1./6.*abs(fluid.get( (*i)[0],(*i)[1] ).f[7] - fluid.get( (*i)[0],(*i)[1] ).f[3])/abs(tmp7-tmp3);
        }
*/
    }
    double zxx, zxy, zyy;
    double alp = 0;
    double M = 3.141592654*getRadius()*getRadius()*getDensity();
    zxx = (M+alp*Zyy)/((M+alp*Zxx)*(M+alp*Zyy)+alp*alp*Zxy*Zxy);
    zxy =    alp*Zxy /((M+alp*Zxx)*(M+alp*Zyy)+alp*alp*Zxy*Zxy);
    zyy = (M+alp*Zxx)/((M+alp*Zxx)*(M+alp*Zyy)+alp*alp*Zxy*Zxy);

    double tmpX = velocity[0], tmpY = velocity[1];
    velocity[0] += zxx*(Mx-Zxx*tmpX-Zxy*tmpY)+zxy*(My-Zxy*tmpX-Zyy*tmpY);
    velocity[1] += zxy*(Mx-Zxx*tmpX-Zxy*tmpY)+zyy*(My-Zxy*tmpX-Zyy*tmpY);
    return;
}


template<typename T, template<typename U> class Descriptor>
void BulkParticle2D<T,Descriptor>::particleToParticleMagnetic(BulkParticleField2D<T,Descriptor>& bulkParticleField_)
{
    std::cout<<"particleToParticleMagnetic:"<<std::endl;
    Array<T,2> Force;
    Force[0] = 0, Force[1] = 0;
    T epsilon = 10.0;

    //repulsive force
    for(pluint iParticle = 0; iParticle < bulkParticleField_.particles.size(); iParticle++){
        if(bulkParticleField_.particles[iParticle]->getId() != this->getId()){
            BulkParticle2D<T,Descriptor> *theOtherParticle = bulkParticleField_.particles[iParticle];
            T distance = sqrt( pow(theOtherParticle->getPosition()[0] - this->getPosition()[0], 2.0) +
                               pow(theOtherParticle->getPosition()[1] - this->getPosition()[1], 2.0) );
            T theOtherRadius = theOtherParticle->getRadius();
            if(     distance > theOtherRadius + this->getRadius() + epsilon )
                Force += 0;
            else if(distance < theOtherRadius + this->getRadius() + epsilon &&
                    distance > theOtherRadius + this->getRadius() ){
                std::cout<<"approaching..."<<std::endl;
                Force[0] += 100*pow(distance - theOtherRadius - this->getRadius() - epsilon, 4.0) *
                           (this->getPosition()[0] - theOtherParticle->getPosition()[0])/distance;
                Force[1] += 100*pow(distance - theOtherRadius - this->getRadius() - epsilon, 4.0) *
                           (this->getPosition()[1] - theOtherParticle->getPosition()[1])/distance;
            }
            else if(distance < theOtherRadius + this->getRadius() && distance > 0 ){
                 std::cout<<"colliding..."<<std::endl;
                 Force[0] +=100*( pow(distance - theOtherRadius - this->getRadius() - epsilon, 4.0) +
                             (theOtherRadius + this->getRadius() - distance)/epsilon ) *
                           (this->getPosition()[0] - theOtherParticle->getPosition()[0])/distance;
                 Force[1] +=100*( pow(distance - theOtherRadius - this->getRadius() - epsilon, 4.0) +
                             (theOtherRadius + this->getRadius() - distance)/epsilon ) *
                           (this->getPosition()[1] - theOtherParticle->getPosition()[1])/distance;
            }
        }
    }

    //magnetic attractive force
    for(pluint iParticle = 0; iParticle < bulkParticleField_.particles.size(); iParticle++){
        if(bulkParticleField_.particles[iParticle]->getId() != this->getId() &&
           bulkParticleField_.particles[iParticle]->getTag() == this->getTag() ){
            T distance = sqrt( pow(bulkParticleField_.particles[iParticle]->getPosition()[0] - this->getPosition()[0], 2.0) +
                               pow(bulkParticleField_.particles[iParticle]->getPosition()[1] - this->getPosition()[1], 2.0) );
            T theOtherRadius = bulkParticleField_.particles[iParticle]->getRadius();
            if(distance > theOtherRadius + this->getRadius() + epsilon ){
                std::cout<<"approaching..."<<std::endl;
                Force[0] += -10000/pow(distance, 2.0) *
                           (this->getPosition()[0] - bulkParticleField_.particles[iParticle]->getPosition()[0])/distance;
                Force[1] += -10000/pow(distance, 2.0) *
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
void BulkParticle2D<T,Descriptor>::particleToParticle(BulkParticleField2D<T,Descriptor>& bulkParticleField_)
{
    std::cout<<"particleToParticle:"<<std::endl;
    Array<T,2> Force;
    Force[0] = 0, Force[1] = 0;
    T epsilon = 10.0;

    for(pluint iParticle = 0; iParticle < bulkParticleField_.particles.size(); iParticle++){
        if(bulkParticleField_.particles[iParticle]->getId() != this->getId()){
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
void BulkParticle2D<T,Descriptor>::particleToChunk(BulkParticleChunk2D<T,Descriptor>& bulkParticleChunk_){
    std::cout<<"particleToChunk:"<<std::endl;
    Array<T,2> Force;
    Force[0] = 0, Force[1] = 0;
    T epsilon = 10.0;

    for(pluint iParticle = 0; iParticle < bulkParticleChunk_.particles.size(); iParticle++){
        if(bulkParticleChunk_.particles[iParticle]->getId() != this->getId()){
            T distance = sqrt( pow(bulkParticleChunk_.particles[iParticle]->getPosition()[0] - this->getPosition()[0], 2.0) +
                               pow(bulkParticleChunk_.particles[iParticle]->getPosition()[1] - this->getPosition()[1], 2.0) );
            T theOtherRadius = bulkParticleChunk_.particles[iParticle]->getRadius();
            if(     distance > theOtherRadius + this->getRadius() + epsilon )
                Force += 0;
            else if(distance < theOtherRadius + this->getRadius() + epsilon &&
                    distance > theOtherRadius + this->getRadius() ){
                std::cout<<"approaching..."<<std::endl;
                Force[0] += 100*pow(distance - theOtherRadius - this->getRadius() - epsilon, 4.0) *
                           (this->getPosition()[0] - bulkParticleChunk_.particles[iParticle]->getPosition()[0])/distance;
                Force[1] += 100*pow(distance - theOtherRadius - this->getRadius() - epsilon, 4.0) *
                           (this->getPosition()[1] - bulkParticleChunk_.particles[iParticle]->getPosition()[1])/distance;
            }
            else if(distance < theOtherRadius + this->getRadius() && distance > 0 ){
                 std::cout<<"colliding..."<<std::endl;
                 Force[0] +=100*( pow(distance - theOtherRadius - this->getRadius() - epsilon, 4.0) +
                             (theOtherRadius + this->getRadius() - distance)/epsilon ) *
                           (this->getPosition()[0] - bulkParticleChunk_.particles[iParticle]->getPosition()[0])/distance;
                 Force[1] +=100*( pow(distance - theOtherRadius - this->getRadius() - epsilon, 4.0) +
                             (theOtherRadius + this->getRadius() - distance)/epsilon ) *
                           (this->getPosition()[1] - bulkParticleChunk_.particles[iParticle]->getPosition()[1])/distance;
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
    T epsilon = 10.0;

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
            Force[0] += 100*pow(distance - 2*this->getRadius() - epsilon, 4.0) *
                               (this->getPosition()[0] - theOtherPosition[0])/distance;
            Force[1] += 100*pow(distance - 2*this->getRadius() - epsilon, 4.0) *
                               (this->getPosition()[1] - theOtherPosition[1])/distance;
        }
        else if(distance < 2*this->getRadius() && distance > 0){
            std::cout<<"colliding with wall..."<<std::endl;
            system("PAUSE");
            Force[0] += 100*( pow(distance - 2*this->getRadius() - epsilon, 4.0) +
                               (2*this->getRadius() - distance)/epsilon ) *
                               (this->getPosition()[0] - theOtherPosition[0])/distance;
            Force[1] += 100*( pow(distance - 2*this->getRadius() - epsilon, 4.0) +
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
void BulkParticle2D<T,Descriptor>::wallToParticle(BlockLattice2D<T,Descriptor>& fluid)
{
    std::cout<<"wallToParticle:"<<std::endl;
    Array<T,2> Force;
    Force[0] = 0;
    Force[1] = 0;
    T epsilon = 10.0;

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
            Force[0] += 200*pow(distance - 2*this->getRadius() - epsilon, 4.0) *
                               (this->getPosition()[0] - theOtherPosition[0])/distance;
            Force[1] += 200*pow(distance - 2*this->getRadius() - epsilon, 4.0) *
                               (this->getPosition()[1] - theOtherPosition[1])/distance;
        }
        else if(distance < 2*this->getRadius() && distance > 0){
            std::cout<<"colliding with wall..."<<std::endl;
            system("PAUSE");
            Force[0] += 200*( pow(distance - 2*this->getRadius() - epsilon, 4.0) +
                               (2*this->getRadius() - distance)/epsilon ) *
                               (this->getPosition()[0] - theOtherPosition[0])/distance;
            Force[1] += 200*( pow(distance - 2*this->getRadius() - epsilon, 4.0) +
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


//=============================free function============================================
template<typename T>
void computeSolidFraction(Array<int,2> cellPosition,Array<T,2> particlePosition, T radius, double& solidFrac){
    //std::cout<<"cell"    <<cellPosition[0]    <<","<<cellPosition[1]    <<std::endl;
    //std::cout<<"particle"<<particlePosition[0]<<","<<particlePosition[1]<<std::endl;
    Array<T,2> UpLeft;
    Array<T,2> UpRight;
    Array<T,2> DownLeft;
    Array<T,2> DownRight;
    UpLeft[0] = (T)cellPosition[0] - 0.5;
    UpLeft[1] = (T)cellPosition[1] + 0.5;
    UpRight[0] = (T)cellPosition[0] + 0.5;
    UpRight[1] = (T)cellPosition[1] + 0.5;
    DownLeft[0] = (T)cellPosition[0] - 0.5;
    DownLeft[1] = (T)cellPosition[1] - 0.5;
    DownRight[0] = (T)cellPosition[0] + 0.5;
    DownRight[1] = (T)cellPosition[1] - 0.5;
//    std::cout<<UpLeft[0]<<","<<UpLeft[1]<<","<<UpRight[0]<<","<<UpRight[1]<<
//          DownLeft[0]<<","<<DownLeft[1]<<","<<DownRight[0]<<","<<DownRight[1]<<std::endl;
    T upSideLength, downSideLength;
    T leftSideLength, rightSideLength;
//There are twelve cases
//1)trapezium on the left; 2)trapezium on the right;
//3)trapezium on the up;   4)trapezium at the bottom;
//5)triangle on the upleft; 6)anti-triangle on the upleft;
//7)triangle on the upright; 8)anti-triangle on the upright;
//9)triangle on the downleft; 10)anti-triangle on the downleft;
//11)triangle on the downright; 12)anti-triangle on the downright.
//
    if(pow(UpLeft[0]-particlePosition[0],2.0) + pow(UpLeft[1]-particlePosition[1],2.0) >= radius*radius &&
       pow(DownLeft[0]-particlePosition[0],2.0) + pow(DownLeft[1]-particlePosition[1],2.0) >= radius*radius &&
       pow(UpRight[0]-particlePosition[0],2.0) + pow(UpRight[1]-particlePosition[1],2.0) <= radius*radius &&
       pow(DownRight[0]-particlePosition[0],2.0) + pow(DownRight[1]-particlePosition[1],2.0) <= radius*radius){
        //case 1: trapezium on the left
        //std::cout<<"case 1"<<std::endl;
        upSideLength = sqrt( pow(radius,2.0)-pow(particlePosition[1]-(cellPosition[1]+0.5),2.0) ) -
                       abs(particlePosition[0] - (cellPosition[0]+0.5) );
        downSideLength = sqrt( pow(radius,2.0)-pow(particlePosition[1]-(cellPosition[1]-0.5),2.0) ) -
                       abs(particlePosition[0] - (cellPosition[0]+0.5) );
        solidFrac = (upSideLength+downSideLength)/2;
    }
    else if(pow(UpLeft[0]-particlePosition[0],2.0) + pow(UpLeft[1]-particlePosition[1],2.0) <= radius*radius &&
       pow(DownLeft[0]-particlePosition[0],2.0) + pow(DownLeft[1]-particlePosition[1],2.0) <= radius*radius &&
       pow(UpRight[0]-particlePosition[0],2.0) + pow(UpRight[1]-particlePosition[1],2.0) >= radius*radius &&
       pow(DownRight[0]-particlePosition[0],2.0) + pow(DownRight[1]-particlePosition[1],2.0) >= radius*radius){
        //case 2: trapezium on the right
        //std::cout<<"case 2"<<std::endl;
        upSideLength = sqrt ( pow(radius,2.0)-pow(particlePosition[1]-(cellPosition[1]+0.5),2.0) ) -
                       abs(particlePosition[0] - (cellPosition[0]-0.5) );
        downSideLength = sqrt ( pow(radius,2.0)-pow(particlePosition[1]-(cellPosition[1]-0.5),2.0) ) -
                       abs(particlePosition[0] - (cellPosition[0]-0.5) );
        solidFrac = (upSideLength+downSideLength)/2;
    }
    else if(pow(UpLeft[0]-particlePosition[0],2.0) + pow(UpLeft[1]-particlePosition[1],2.0) >= radius*radius &&
       pow(UpRight[0]-particlePosition[0],2.0) + pow(UpRight[1]-particlePosition[1],2.0) >= radius*radius &&
       pow(DownLeft[0]-particlePosition[0],2.0) + pow(DownLeft[1]-particlePosition[1],2.0) <= radius*radius &&
       pow(DownRight[0]-particlePosition[0],2.0) + pow(DownRight[1]-particlePosition[1],2.0) <= radius*radius){
        //case 3: trapezium on the top
        //std::cout<<"case 3"<<std::endl;
        leftSideLength = sqrt ( pow(radius,2.0)-pow(particlePosition[0]-(cellPosition[0]-0.5),2.0) ) -
                         abs(particlePosition[1] - (cellPosition[1]-0.5) );
        rightSideLength = sqrt ( pow(radius,2.0)-pow(particlePosition[0]-(cellPosition[0]+0.5),2.0) ) -
                         abs(particlePosition[1] - (cellPosition[1]-0.5) );
        solidFrac = (leftSideLength+rightSideLength)/2;
    }
    else if(pow(UpLeft[0]-particlePosition[0],2.0) + pow(UpLeft[1]-particlePosition[1],2.0) <= radius*radius &&
       pow(UpRight[0]-particlePosition[0],2.0) + pow(UpRight[1]-particlePosition[1],2.0) <= radius*radius &&
       pow(DownLeft[0]-particlePosition[0],2.0) + pow(DownLeft[1]-particlePosition[1],2.0) >= radius*radius &&
       pow(DownRight[0]-particlePosition[0],2.0) + pow(DownRight[1]-particlePosition[1],2.0) >= radius*radius){
        //case 4: trapezium at the bottom
        //std::cout<<"case 4"<<std::endl;
        leftSideLength = sqrt ( pow(radius,2.0)-pow(particlePosition[0]-(cellPosition[0]-0.5),2.0) ) -
                         abs(particlePosition[1] - (cellPosition[1]+0.5) );
        rightSideLength = sqrt ( pow(radius,2.0)-pow(particlePosition[0]-(cellPosition[0]+0.5),2.0) ) -
                          abs(particlePosition[1] - (cellPosition[1]+0.5) );
        solidFrac = (leftSideLength+rightSideLength)/2;
    }
    else if(pow(UpLeft[0]-particlePosition[0],2.0) + pow(UpLeft[1]-particlePosition[1],2.0) > radius*radius &&
       pow(UpRight[0]-particlePosition[0],2.0) + pow(UpRight[1]-particlePosition[1],2.0) >= radius*radius &&
       pow(DownLeft[0]-particlePosition[0],2.0) + pow(DownLeft[1]-particlePosition[1],2.0) >= radius*radius &&
       pow(DownRight[0]-particlePosition[0],2.0) + pow(DownRight[1]-particlePosition[1],2.0) < radius*radius){
        //case 5: triangle at the upleft
        //std::cout<<"case 5"<<std::endl;
        rightSideLength = sqrt ( pow(radius,2.0)-pow(particlePosition[0]-(cellPosition[0]+0.5),2.0) ) -
                         abs(particlePosition[1] - (cellPosition[1]-0.5) );
        downSideLength = sqrt( pow(radius,2.0)-pow(particlePosition[1]-(cellPosition[1]-0.5),2.0) ) -
                       abs(particlePosition[0] - (cellPosition[0]+0.5) );
        solidFrac = rightSideLength*downSideLength/2;
    }
    else if(pow(UpLeft[0]-particlePosition[0],2.0) + pow(UpLeft[1]-particlePosition[1],2.0) > radius*radius &&
       pow(UpRight[0]-particlePosition[0],2.0) + pow(UpRight[1]-particlePosition[1],2.0) < radius*radius &&
       pow(DownLeft[0]-particlePosition[0],2.0) + pow(DownLeft[1]-particlePosition[1],2.0) < radius*radius &&
       pow(DownRight[0]-particlePosition[0],2.0) + pow(DownRight[1]-particlePosition[1],2.0) < radius*radius){
        //case 6: anti-triangle at the upleft
        //std::cout<<"case 6:"<<std::endl;
        leftSideLength = sqrt ( pow(radius,2.0)-pow(particlePosition[0]-(cellPosition[0]-0.5),2.0) ) -
                         abs(particlePosition[1] - (cellPosition[1]-0.5) );
        upSideLength = sqrt ( pow(radius,2.0)-pow(particlePosition[1]-(cellPosition[1]+0.5),2.0) ) -
                       abs(particlePosition[0] - (cellPosition[0]+0.5) );
        solidFrac = 1-(1-leftSideLength)*(1-upSideLength)/2;
    }
    else if(pow(UpLeft[0]-particlePosition[0],2.0) + pow(UpLeft[1]-particlePosition[1],2.0) >= radius*radius &&
       pow(UpRight[0]-particlePosition[0],2.0) + pow(UpRight[1]-particlePosition[1],2.0) > radius*radius &&
       pow(DownLeft[0]-particlePosition[0],2.0) + pow(DownLeft[1]-particlePosition[1],2.0) < radius*radius &&
       pow(DownRight[0]-particlePosition[0],2.0) + pow(DownRight[1]-particlePosition[1],2.0) >= radius*radius){
        //case 7: triangle at the upright
        //std::cout<<"case 7:"<<std::endl;
        leftSideLength = sqrt ( pow(radius,2.0)-pow(particlePosition[0]-(cellPosition[0]-0.5),2.0) ) -
                         abs(particlePosition[1] - (cellPosition[1]-0.5) );
        downSideLength = sqrt( pow(radius,2.0)-pow(particlePosition[1]-(cellPosition[1]-0.5),2.0) ) -
                         abs(particlePosition[0] - (cellPosition[0]-0.5) );
        solidFrac = leftSideLength*downSideLength/2;
    }
    else if(pow(UpLeft[0]-particlePosition[0],2.0) + pow(UpLeft[1]-particlePosition[1],2.0) < radius*radius &&
       pow(UpRight[0]-particlePosition[0],2.0) + pow(UpRight[1]-particlePosition[1],2.0) > radius*radius &&
       pow(DownLeft[0]-particlePosition[0],2.0) + pow(DownLeft[1]-particlePosition[1],2.0) < radius*radius &&
       pow(DownRight[0]-particlePosition[0],2.0) + pow(DownRight[1]-particlePosition[1],2.0) < radius*radius){
        //case 8: anti-triangle at the upright
        //std::cout<<"case 8:"<<std::endl;
        rightSideLength = sqrt ( pow(radius,2.0)-pow(particlePosition[0]-(cellPosition[0]+0.5),2.0) ) -
                         abs(particlePosition[1] - (cellPosition[1]-0.5) );
        upSideLength = sqrt ( pow(radius,2.0)-pow(particlePosition[1]-(cellPosition[1]+0.5),2.0) ) -
                       abs(particlePosition[0] - (cellPosition[0]-0.5) );
        solidFrac = 1-(1-rightSideLength)*(1-upSideLength)/2;
    }
    else if(pow(UpLeft[0]-particlePosition[0],2.0) + pow(UpLeft[1]-particlePosition[1],2.0) >= radius*radius &&
       pow(UpRight[0]-particlePosition[0],2.0) + pow(UpRight[1]-particlePosition[1],2.0) < radius*radius &&
       pow(DownLeft[0]-particlePosition[0],2.0) + pow(DownLeft[1]-particlePosition[1],2.0) > radius*radius &&
       pow(DownRight[0]-particlePosition[0],2.0) + pow(DownRight[1]-particlePosition[1],2.0) >= radius*radius){
        //case 9: triangle at the downleft
        //std::cout<<"case 9:"<<std::endl;
        rightSideLength = sqrt ( pow(radius,2.0)-pow(particlePosition[0]-(cellPosition[0]+0.5),2.0) ) -
                         abs(particlePosition[1] - (cellPosition[1]+0.5) );
        upSideLength = sqrt ( pow(radius,2.0)-pow(particlePosition[1]-(cellPosition[1]+0.5),2.0) ) -
                       abs(particlePosition[0] - (cellPosition[0]+0.5) );
        solidFrac = rightSideLength*upSideLength/2;
    }
    else if(pow(UpLeft[0]-particlePosition[0],2.0) + pow(UpLeft[1]-particlePosition[1],2.0) < radius*radius &&
       pow(UpRight[0]-particlePosition[0],2.0) + pow(UpRight[1]-particlePosition[1],2.0) < radius*radius &&
       pow(DownLeft[0]-particlePosition[0],2.0) + pow(DownLeft[1]-particlePosition[1],2.0) > radius*radius &&
       pow(DownRight[0]-particlePosition[0],2.0) + pow(DownRight[1]-particlePosition[1],2.0) < radius*radius){
        //case 10: anti-triangle at the downleft
        //std::cout<<"case 10:"<<std::endl;
        leftSideLength = sqrt ( pow(radius,2.0)-pow(particlePosition[0]-(cellPosition[0]-0.5),2.0) ) -
                         abs(particlePosition[1] - (cellPosition[1]+0.5) );
        downSideLength = sqrt( pow(radius,2.0)-pow(particlePosition[1]-(cellPosition[1]-0.5),2.0) ) -
                       abs(particlePosition[0] - (cellPosition[0]+0.5) );
        solidFrac = 1-(1-leftSideLength)*(1-downSideLength)/2;
    }
    else if(pow(UpLeft[0]-particlePosition[0],2.0) + pow(UpLeft[1]-particlePosition[1],2.0) < radius*radius &&
       pow(UpRight[0]-particlePosition[0],2.0) + pow(UpRight[1]-particlePosition[1],2.0) >= radius*radius &&
       pow(DownLeft[0]-particlePosition[0],2.0) + pow(DownLeft[1]-particlePosition[1],2.0) >= radius*radius &&
       pow(DownRight[0]-particlePosition[0],2.0) + pow(DownRight[1]-particlePosition[1],2.0) > radius*radius){
        //case 11: triangle at the downright
        //std::cout<<"case 11:"<<std::endl;
        leftSideLength = sqrt ( pow(radius,2.0)-pow(particlePosition[0]-(cellPosition[0]-0.5),2.0) ) -
                         abs(particlePosition[1] - (cellPosition[1]+0.5) );
        upSideLength = sqrt ( pow(radius,2.0)-pow(particlePosition[1]-(cellPosition[1]+0.5),2.0) ) -
                       abs(particlePosition[0] - (cellPosition[0]-0.5) );
        solidFrac = leftSideLength*upSideLength/2;
    }
    else if(pow(UpLeft[0]-particlePosition[0],2.0) + pow(UpLeft[1]-particlePosition[1],2.0) < radius*radius &&
       pow(UpRight[0]-particlePosition[0],2.0) + pow(UpRight[1]-particlePosition[1],2.0) < radius*radius &&
       pow(DownLeft[0]-particlePosition[0],2.0) + pow(DownLeft[1]-particlePosition[1],2.0) < radius*radius &&
       pow(DownRight[0]-particlePosition[0],2.0) + pow(DownRight[1]-particlePosition[1],2.0) > radius*radius){
        //case 12: anti-triangle at the downright
        //std::cout<<"case 12:"<<std::endl;
        rightSideLength = sqrt ( pow(radius,2.0)-pow(particlePosition[0]-(cellPosition[0]+0.5),2.0) ) -
                         abs(particlePosition[1] - (cellPosition[1]+0.5) );
        downSideLength = sqrt( pow(radius,2.0)-pow(particlePosition[1]-(cellPosition[1]-0.5),2.0) ) -
                       abs(particlePosition[0] - (cellPosition[0]-0.5) );
        solidFrac = 1-(1-rightSideLength)*(1-downSideLength)/2;
    }

    else if(pow(UpLeft[0]-particlePosition[0],2.0) + pow(UpLeft[1]-particlePosition[1],2.0) < radius*radius &&
       pow(UpRight[0]-particlePosition[0],2.0) + pow(UpRight[1]-particlePosition[1],2.0) < radius*radius &&
       pow(DownLeft[0]-particlePosition[0],2.0) + pow(DownLeft[1]-particlePosition[1],2.0) < radius*radius &&
       pow(DownRight[0]-particlePosition[0],2.0) + pow(DownRight[1]-particlePosition[1],2.0) < radius*radius){
        //case 13: internal
        //std::cout<<"case 13:"<<std::endl;
        solidFrac = 1;
    }

    else if(pow(UpLeft[0]-particlePosition[0],2.0) + pow(UpLeft[1]-particlePosition[1],2.0) > radius*radius &&
       pow(UpRight[0]-particlePosition[0],2.0) + pow(UpRight[1]-particlePosition[1],2.0) > radius*radius &&
       pow(DownLeft[0]-particlePosition[0],2.0) + pow(DownLeft[1]-particlePosition[1],2.0) > radius*radius &&
       pow(DownRight[0]-particlePosition[0],2.0) + pow(DownRight[1]-particlePosition[1],2.0) > radius*radius){
        //case 14: external
        //std::cout<<"case 14:"<<std::endl;
        solidFrac = 0;
    }

    else{
        std::cout<<"not included in the above 14 cases."<<std::endl;
    }
    if(solidFrac>1 || solidFrac<0){
        std::cout<<"\n\n\n\n\n\n\n\n\n\nsolidFrac="<<solidFrac<<std::endl;
        system("PAUSE");
    }
}

}
#endif

