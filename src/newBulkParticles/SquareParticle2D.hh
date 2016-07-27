#ifndef SQUARE_PARTICLE_2D_HH
#define SQUARE_PARTICLE_2D_HH

#include "core/globalDefs.h"
#include "finiteDifference/interpolations2D.h"
#include "particles/particleIdentifiers2D.h"
#include "newBulkParticles/SquareParticle2D.h"
#include "newBulkParticles/solidBlockLattice2D.h"
#include "newBulkParticles/SquareParticleField2D.h"
#include <cmath>
#include <vector>
#define TRUE 1
#define FALSE 0

namespace plb {

/* *************** class BulkParticle2D ***************************************** */

template<typename T, template<typename U> class Descriptor>
SquareParticle2D<T,Descriptor>::SquareParticle2D()
    : Particle2D<T,Descriptor>(),
      AngularPosition(0),
      length(1)
{ }

template<typename T, template<typename U> class Descriptor>
SquareParticle2D<T,Descriptor>::SquareParticle2D(plint tag_, Array<T,2> const& position_, T AngularPosition_,
                                                         Array<T,2> const& velocity_, T AngularVelocity_, T length_, T density_)
    : Particle2D<T,Descriptor>(tag_, position_),
      AngularPosition(AngularPosition_),
      velocity(velocity_),
      AngularVelocity(AngularVelocity_),
      length(length_),
      density(density_)
{
}

template<typename T, template<typename U> class Descriptor>
void SquareParticle2D<T,Descriptor>::reset(Array<T,2> const& position_, T AngularPosition_)
{
    this->position = position_;
    AngularPosition = AngularPosition_;
}


template<typename T, template<typename U> class Descriptor>
void SquareParticle2D<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
{
    serializer.addValue(this->tag);
    serializer.addValues<T,2>(this->position);
    serializer.addValue(AngularPosition);
    serializer.addValue(length);
}

template<typename T, template<typename U> class Descriptor>
void SquareParticle2D<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer)
{
    unserializer.readValue(this->tag);
    unserializer.readValues<T,2>(this->position);
    unserializer.readValue(AngularPosition);
    unserializer.readValue(length);
}


template<typename T, template<typename U> class Descriptor>
SquareParticle2D<T,Descriptor>* SquareParticle2D<T,Descriptor>::clone() const {
    return new SquareParticle2D<T,Descriptor>(*this);
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
void SquareParticle2D<T,Descriptor>::rescale(int dxScale, int dtScale) {
    int dimDx = 1;
    int dimDt = 0;
    T scaleFactor = scaleFromReference(dxScale, dimDx, dtScale, dimDt);
    this->position *= scaleFactor;
    length *= scaleFactor;
}

template<typename T, template<typename U> class Descriptor>
void SquareParticle2D<T,Descriptor>::computeBoundaries(Box2D domain) {
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
void SquareParticle2D<T,Descriptor>::velocityToParticle(TensorField2D<T,2>& velocityField, T scaling)
{
    velocity = predictorCorrectorTensorField<T,2>(velocityField, this->getPosition(), scaling);
}

template<typename T, template<typename U> class Descriptor>
void SquareParticle2D<T,Descriptor>::rhoBarJtoParticle (
        NTensorField2D<T>& rhoBarJfield, bool velIsJ, T scaling )
{

}





template<typename T, template<typename U> class Descriptor>
void SquareParticle2D<T,Descriptor>::fluidToParticle(solidBlockLattice2D<T,Descriptor>& fluid, T scaling)
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
        //computeSolidFraction(cellPosition, solidFrac);
        T B = 0;//solidFrac;
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
        Torque +=-(this->position[1]-cellPosition[1])*cellForce[0]*scaling*B
                 +(this->position[0]-cellPosition[0])*cellForce[1]*scaling*B;
    }

    system("PAUSE");
    std::cout<<"SolidBoundary:"<<std::endl;
    for(std::vector<Array<int,2> >::iterator i=solidBoundary.begin(); i!=solidBoundary.end(); ++i){
        Array<int,2> cellPosition = (*i);
        //computeSolidFraction(cellPosition, solidFrac);
        T B = 1;//solidFrac;
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
        Torque +=-(this->position[1]-cellPosition[1])*cellForce[0]*scaling*B
                 +(this->position[0]-cellPosition[0])*cellForce[1]*scaling*B;
    }

    std::cout<<"internalSolid:"<<std::endl;
    for(std::vector<Array<int,2> >::iterator i=internalSolid.begin(); i!=internalSolid.end(); ++i){
        Array<int,2> cellPosition = (*i);
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
    double mass = getLength()*getLength()*getDensity();
    double moment_of_inertia = 2*3.141592654*getLength()*getLength()*getLength()*getLength()/4*getDensity();
    velocity[0] += -Force[0]/mass;
    velocity[1] += -Force[1]/mass;
    AngularVelocity += -Torque/moment_of_inertia;
    std::cout<<"After fluidToParticle acceleration, velocity="<<velocity[0]<<","<<velocity[1]<<std::endl;
    std::cout<<"AngularVelocity="<<AngularVelocity<<std::endl;

}





template<typename T, template<typename U> class Descriptor>
void SquareParticle2D<T,Descriptor>::advance() {
    PLB_ASSERT( norm(velocity)<1. );
    std::cout<<this->position[0]<<","<<this->position[1]<<std::endl;
    std::cout<<this->velocity[0]<<","<<this->velocity[1]<<std::endl;
    this->position[0] += velocity[0];
    this->position[1] += velocity[1];
    AngularPosition += AngularVelocity;
}


/*
template<typename T, template<typename U> class Descriptor>
void SquareParticle2D<T,Descriptor>::computeSolidFraction(Array<int,2> cellPosition, double& solidFrac){
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

        T theta = AngularPosition;
        T x0 = this->position[0];
        T y0 = this->position[1];
        T l = length;
        //line1 : y - (y0 + l/2*cos(theta) ) =  tan(theta)  *(x - (x0 - l/2*sin(theta) ) )
        //line2 : y - (y0 - l/2*cos(theta) ) =  tan(theta)  *(x - (x0 + l/2*sin(theta) ) )
        //line3 : y - (y0 + l/2*sin(theta) ) = -1/tan(theta)*(x - (x0 + l/2*cos(theta) ) )
        //line4 : y - (y0 - 1/2*sin(theta) ) = -1/tan(theta)*(x - (x0 - l/2*cos(theta) ) )
        std::cout<<"cellPosition"<<cellPosition[0]<<","<<cellPosition[1]<<std::endl;

    if(inside(UpLeft[0],   UpLeft[1])   == FALSE &&
       inside(DownLeft[0], DownLeft[1]) == FALSE &&
       inside(UpRight[0],  UpRight[1])  == TRUE &&
       inside(DownRight[0],DownRight[1])== TRUE){
        //case 1: trapezium on the left
        std::cout<<"case 1:"<<std::endl;
        T xUp, xDown;
        if( line1(UpLeft[0],   UpLeft[1])   * line1(UpRight[0],   UpRight[1]) <=0  &&
            line1(DownLeft[0], DownLeft[1]) * line1(DownRight[0], DownRight[1]) <=0)
            {xUp = line1_getX(UpLeft[1]); xDown = line1_getX(DownLeft[1]);}
        else if( line2(UpLeft[0],   UpLeft[1])   * line2(UpRight[0],   UpRight[1]) <=0  &&
                 line2(DownLeft[0], DownLeft[1]) * line2(DownRight[0], DownRight[1]) <=0)
            {xUp = line2_getX(UpLeft[1]); xDown = line2_getX(DownLeft[1]);}
        else if( line3(UpLeft[0],   UpLeft[1])   * line3(UpRight[0],   UpRight[1]) <=0  &&
                 line3(DownLeft[0], DownLeft[1]) * line3(DownRight[0], DownRight[1]) <=0)
            {xUp = line3_getX(UpLeft[1]); xDown = line3_getX(DownLeft[1]);}
        else if( line4(UpLeft[0],   UpLeft[1])   * line4(UpRight[0],   UpRight[1]) <=0  &&
                 line4(DownLeft[0], DownLeft[1]) * line4(DownRight[0], DownRight[1]) <=0)
            {xUp = line4_getX(UpLeft[1]); xDown = line4_getX(DownLeft[1]);}
        std::cout<<"xUp="<<xUp<<",xDown="<<xDown<<std::endl;
        upSideLength   = UpRight[0]   - xUp;
        downSideLength = DownRight[0] - xDown;
        solidFrac = (upSideLength+downSideLength)/2;
    }
    else if(inside(UpLeft[0],   UpLeft[1])   == TRUE &&
            inside(DownLeft[0], DownLeft[1]) == TRUE &&
            inside(UpRight[0],  UpRight[1])  == FALSE &&
            inside(DownRight[0],DownRight[1]) == FALSE){
        //case 2: trapezium on the right
        std::cout<<"case 2:"<<std::endl;
        T xUp, xDown;
        if( line1(UpLeft[0],   UpLeft[1])   * line1(UpRight[0],   UpRight[1]) <=0  &&
            line1(DownLeft[0], DownLeft[1]) * line1(DownRight[0], DownRight[1]) <=0)
            {xUp = line1_getX(UpLeft[1]); xDown = line1_getX(DownLeft[1]);}
        else if( line2(UpLeft[0],   UpLeft[1])   * line2(UpRight[0],   UpRight[1]) <=0  &&
                 line2(DownLeft[0], DownLeft[1]) * line2(DownRight[0], DownRight[1]) <=0)
            {xUp = line2_getX(UpLeft[1]); xDown = line2_getX(DownLeft[1]);}
        else if( line3(UpLeft[0],   UpLeft[1])   * line3(UpRight[0],   UpRight[1]) <=0  &&
                 line3(DownLeft[0], DownLeft[1]) * line3(DownRight[0], DownRight[1]) <=0)
            {xUp = line3_getX(UpLeft[1]); xDown = line3_getX(DownLeft[1]);}
        else if( line4(UpLeft[0],   UpLeft[1])   * line4(UpRight[0],   UpRight[1]) <=0  &&
                 line4(DownLeft[0], DownLeft[1]) * line4(DownRight[0], DownRight[1]) <=0)
            {xUp = line4_getX(UpLeft[1]); xDown = line4_getX(DownLeft[1]);}
        std::cout<<"xUp="<<xUp<<",xDown="<<xDown<<std::endl;
        upSideLength   = xUp   - UpLeft[0];
        downSideLength = xDown - DownLeft[0];
        solidFrac = (upSideLength+downSideLength)/2;
    }
    else if(inside(UpLeft[0],   UpLeft[1])   == FALSE &&
            inside(UpRight[0],  UpRight[1])  == FALSE &&
            inside(DownLeft[0], DownLeft[1]) == TRUE &&
            inside(DownRight[0],DownRight[1])== TRUE){
        //case 3: trapezium on the top
        std::cout<<"case 3:"<<std::endl;
        T yLeft, yRight;
        if( line1(UpLeft[0],   UpLeft[1])   * line1(DownLeft[0],  DownLeft[1]) <=0  &&
            line1(UpRight[0],  UpRight[1])  * line1(DownRight[0], DownRight[1]) <=0)
            {yLeft = line1_getY(UpLeft[0]); yRight = line1_getY(UpRight[0]);}
        else if( line2(UpLeft[0],   UpLeft[1])   * line2(DownLeft[0],  DownLeft[1]) <=0  &&
                 line2(UpRight[0],  UpRight[1])  * line2(DownRight[0], DownRight[1]) <=0)
            {yLeft = line2_getY(UpLeft[0]); yRight = line2_getY(UpRight[0]);}
        else if( line3(UpLeft[0],   UpLeft[1])   * line3(DownLeft[0],  DownLeft[1]) <=0  &&
                 line3(UpRight[0],  UpRight[1])  * line3(DownRight[0], DownRight[1]) <=0)
            {yLeft = line3_getY(UpLeft[0]); yRight = line3_getY(UpRight[0]);}
        else if( line4(UpLeft[0],   UpLeft[1])   * line4(DownLeft[0],  DownLeft[1]) <=0  &&
                 line4(UpRight[0],  UpRight[1])  * line4(DownRight[0], DownRight[1]) <=0)
            {yLeft = line4_getY(UpLeft[0]); yRight = line4_getY(UpRight[0]);}
        std::cout<<"yLeft="<<yLeft<<",yRight="<<yRight<<std::endl;
        leftSideLength  = yLeft - DownLeft[1];
        rightSideLength = yRight- DownRight[1];
        solidFrac = (leftSideLength+rightSideLength)/2;
    }
    else if(inside(UpLeft[0],   UpLeft[1])   == TRUE  &&
            inside(UpRight[0],  UpRight[1])  == TRUE  &&
            inside(DownLeft[0], DownLeft[1]) == FALSE &&
            inside(DownRight[0],DownRight[1])== FALSE){
        //case 4: trapezium at the bottom
        std::cout<<"case 4:"<<std::endl;
        T yLeft, yRight;
        if( line1(UpLeft[0],   UpLeft[1])   * line1(DownLeft[0],  DownLeft[1]) <=0  &&
            line1(UpRight[0],  UpRight[1])  * line1(DownRight[0], DownRight[1]) <=0)
            {yLeft = line1_getY(UpLeft[0]); yRight = line1_getY(UpRight[0]);}
        else if( line2(UpLeft[0],   UpLeft[1])   * line2(DownLeft[0],  DownLeft[1]) <=0  &&
                 line2(UpRight[0],  UpRight[1])  * line2(DownRight[0], DownRight[1]) <=0)
            {yLeft = line2_getY(UpLeft[0]); yRight = line2_getY(UpRight[0]);}
        else if( line3(UpLeft[0],   UpLeft[1])   * line3(DownLeft[0],  DownLeft[1]) <=0  &&
                 line3(UpRight[0],  UpRight[1])  * line3(DownRight[0], DownRight[1]) <=0)
            {yLeft = line3_getY(UpLeft[0]); yRight = line3_getY(UpRight[0]);}
        else if( line4(UpLeft[0],   UpLeft[1])   * line4(DownLeft[0],  DownLeft[1]) <=0  &&
                 line4(UpRight[0],  UpRight[1])  * line4(DownRight[0], DownRight[1]) <=0)
            {yLeft = line4_getY(UpLeft[0]); yRight = line4_getY(UpRight[0]);}
        std::cout<<"yLeft="<<yLeft<<",yRight="<<yRight<<std::endl;
        leftSideLength  = UpLeft[1] - yLeft;
        rightSideLength = UpRight[1]- yRight;
        solidFrac = (leftSideLength+rightSideLength)/2;
    }
    else if(inside(UpLeft[0],   UpLeft[1])   == FALSE  &&
            inside(UpRight[0],  UpRight[1])  == FALSE  &&
            inside(DownLeft[0], DownLeft[1]) == FALSE &&
            inside(DownRight[0],DownRight[1])== TRUE ){
        //case 5: triangle at the upleft
        std::cout<<"case 5:"<<std::endl;
        T x, y;
        if( line1(DownRight[0], DownRight[1])   * line1(DownLeft[0],  DownLeft[1]) <=0  &&
            line1(DownRight[0], DownRight[1])   * line1(UpRight[0],   UpRight[1])  <=0  &&
            line1(DownRight[0], DownRight[1])   * line1(UpLeft[0],    UpLeft[1])   <=0 )
            {x = line1_getX(DownRight[1]); y = line1_getY(DownRight[0]);}
        else if( line2(DownRight[0], DownRight[1])   * line2(DownLeft[0],  DownLeft[1]) <=0  &&
            line2(DownRight[0], DownRight[1])   * line2(UpRight[0],   UpRight[1])  <=0  &&
            line2(DownRight[0], DownRight[1])   * line2(UpLeft[0],    UpLeft[1])   <=0 )
            {x = line2_getX(DownRight[1]); y = line2_getY(DownRight[0]);}
        else if( line3(DownRight[0], DownRight[1])   * line3(DownLeft[0],  DownLeft[1]) <=0  &&
            line3(DownRight[0], DownRight[1])   * line3(UpRight[0],   UpRight[1])  <=0  &&
            line3(DownRight[0], DownRight[1])   * line3(UpLeft[0],    UpLeft[1])   <=0 )
            {x = line3_getX(DownRight[1]); y = line3_getY(DownRight[0]);}
        else if( line4(DownRight[0], DownRight[1])   * line4(DownLeft[0],  DownLeft[1]) <=0  &&
            line4(DownRight[0], DownRight[1])   * line4(UpRight[0],   UpRight[1])  <=0  &&
            line4(DownRight[0], DownRight[1])   * line4(UpLeft[0],    UpLeft[1])   <=0 )
            {x = line4_getX(DownRight[1]); y = line4_getY(DownRight[0]);}
        std::cout<<"x="<<x<<",y="<<y<<std::endl;
        rightSideLength = y - DownRight[1];
        downSideLength  = DownRight[0] - x;
        solidFrac = rightSideLength*downSideLength/2;
    }
    else if(inside(UpLeft[0],   UpLeft[1])   == FALSE  &&
            inside(UpRight[0],  UpRight[1])  == TRUE  &&
            inside(DownLeft[0], DownLeft[1]) == TRUE &&
            inside(DownRight[0],DownRight[1])== TRUE){
        //case 6: anti-triangle at the upleft
        std::cout<<"case 6:"<<std::endl;
        T x,y;
        if( line1(UpLeft[0], UpLeft[1])   * line1(DownLeft[0],  DownLeft[1]) <=0  &&
            line1(UpLeft[0], UpLeft[1])   * line1(UpRight[0],   UpRight[1])  <=0  &&
            line1(UpLeft[0], UpLeft[1])   * line1(DownRight[0], DownRight[1])   <=0 )
            {x = line1_getX(UpLeft[1]); y = line1_getY(UpLeft[0]);}
        else if( line2(UpLeft[0], UpLeft[1])   * line2(DownLeft[0],  DownLeft[1]) <=0  &&
            line2(UpLeft[0], UpLeft[1])   * line2(UpRight[0],   UpRight[1])  <=0  &&
            line2(UpLeft[0], UpLeft[1])   * line2(DownRight[0], DownRight[1])   <=0 )
            {x = line2_getX(UpLeft[1]); y = line2_getY(UpLeft[0]);}
        else if( line3(UpLeft[0], UpLeft[1])   * line3(DownLeft[0],  DownLeft[1]) <=0  &&
            line3(UpLeft[0], UpLeft[1])   * line3(UpRight[0],   UpRight[1])  <=0  &&
            line3(UpLeft[0], UpLeft[1])   * line3(DownRight[0], DownRight[1])   <=0 )
            {x = line3_getX(UpLeft[1]); y = line3_getY(UpLeft[0]);}
        else if( line4(UpLeft[0], UpLeft[1])   * line4(DownLeft[0],  DownLeft[1]) <=0  &&
            line4(UpLeft[0], UpLeft[1])   * line4(UpRight[0],   UpRight[1])  <=0  &&
            line4(UpLeft[0], UpLeft[1])   * line4(DownRight[0], DownRight[1])   <=0 )
            {x = line4_getX(UpLeft[1]); y = line4_getY(UpLeft[0]);}
        std::cout<<"x="<<x<<",y="<<y<<std::endl;
        leftSideLength = UpLeft[1] - y;
        upSideLength   = x - UpLeft[0];
        solidFrac = 1-(1-leftSideLength)*(1-upSideLength)/2;
    }
    else if(inside(UpLeft[0],   UpLeft[1])   == FALSE  &&
            inside(UpRight[0],  UpRight[1])  == FALSE  &&
            inside(DownLeft[0], DownLeft[1]) == TRUE  &&
            inside(DownRight[0],DownRight[1])== FALSE){
        //case 7: triangle at the upright
        std::cout<<"case 7:"<<std::endl;
        T x,y;
        if( line1(DownLeft[0], DownLeft[1])   * line1(UpLeft[0],    UpLeft[1]) <=0  &&
            line1(DownLeft[0], DownLeft[1])   * line1(UpRight[0],   UpRight[1])  <=0  &&
            line1(DownLeft[0], DownLeft[1])   * line1(DownRight[0], DownRight[1])   <=0 )
            {x = line1_getX(DownLeft[1]); y = line1_getY(DownLeft[0]);}
        else if( line2(DownLeft[0], DownLeft[1])   * line2(UpLeft[0],    UpLeft[1]) <=0  &&
            line2(DownLeft[0], DownLeft[1])   * line2(UpRight[0],   UpRight[1])  <=0  &&
            line2(DownLeft[0], DownLeft[1])   * line2(DownRight[0], DownRight[1])   <=0 )
            {x = line2_getX(DownLeft[1]); y = line2_getY(DownLeft[0]);}
        else if( line3(DownLeft[0], DownLeft[1])   * line3(UpLeft[0],    UpLeft[1]) <=0  &&
            line3(DownLeft[0], DownLeft[1])   * line3(UpRight[0],   UpRight[1])  <=0  &&
            line3(DownLeft[0], DownLeft[1])   * line3(DownRight[0], DownRight[1])   <=0 )
            {x = line3_getX(DownLeft[1]); y = line3_getY(DownLeft[0]);}
        else if( line4(DownLeft[0], DownLeft[1])   * line4(UpLeft[0],    UpLeft[1]) <=0  &&
            line4(DownLeft[0], DownLeft[1])   * line4(UpRight[0],   UpRight[1])  <=0  &&
            line4(DownLeft[0], DownLeft[1])   * line4(DownRight[0], DownRight[1])   <=0 )
            {x = line4_getX(DownLeft[1]); y = line4_getY(DownLeft[0]);}
        std::cout<<"x="<<x<<",y="<<y<<std::endl;
        leftSideLength = y - DownLeft[1];
        downSideLength = x - DownLeft[0];
        solidFrac = leftSideLength*downSideLength/2;
    }
    else if(inside(UpLeft[0],   UpLeft[1])   == TRUE  &&
            inside(UpRight[0],  UpRight[1])  == FALSE  &&
            inside(DownLeft[0], DownLeft[1]) == TRUE &&
            inside(DownRight[0],DownRight[1])== TRUE){
        //case 8: anti-triangle at the upright
        std::cout<<"case 8:"<<std::endl;
        T x,y;
        if( line1(UpRight[0], UpRight[1])   * line1(DownLeft[0],  DownLeft[1]) <=0  &&
            line1(UpRight[0], UpRight[1])   * line1(DownRight[0], DownRight[1])  <=0  &&
            line1(UpRight[0], UpRight[1])   * line1(UpLeft[0],    UpLeft[1])   <=0 )
            {x = line1_getX(UpRight[1]); y = line1_getY(UpRight[0]);}
        else if( line2(UpRight[0], UpRight[1])   * line2(DownLeft[0],  DownLeft[1]) <=0  &&
            line2(UpRight[0], UpRight[1])   * line2(DownRight[0], DownRight[1])  <=0  &&
            line2(UpRight[0], UpRight[1])   * line2(UpLeft[0],    UpLeft[1])   <=0 )
            {x = line2_getX(UpRight[1]); y = line2_getY(UpRight[0]);}
        else if( line3(UpRight[0], UpRight[1])   * line3(DownLeft[0],  DownLeft[1]) <=0  &&
            line3(UpRight[0], UpRight[1])   * line3(DownRight[0], DownRight[1])  <=0  &&
            line3(UpRight[0], UpRight[1])   * line3(UpLeft[0],    UpLeft[1])   <=0 )
            {x = line3_getX(UpRight[1]); y = line3_getY(UpRight[0]);}
        else if( line4(UpRight[0], UpRight[1])   * line4(DownLeft[0],  DownLeft[1]) <=0  &&
            line4(UpRight[0], UpRight[1])   * line4(DownRight[0], DownRight[1])  <=0  &&
            line4(UpRight[0], UpRight[1])   * line4(UpLeft[0],    UpLeft[1])   <=0 )
            {x = line4_getX(UpRight[1]); y = line4_getY(UpRight[0]);}
        std::cout<<"x="<<x<<",y="<<y<<std::endl;
        rightSideLength = UpRight[1] - y;
        upSideLength    = UpRight[0] - x;
        solidFrac = 1-(1-rightSideLength)*(1-upSideLength)/2;
    }
    else if(inside(UpLeft[0],   UpLeft[1])   == FALSE  &&
            inside(UpRight[0],  UpRight[1])  == TRUE  &&
            inside(DownLeft[0], DownLeft[1]) == FALSE &&
            inside(DownRight[0],DownRight[1])== FALSE){
        //case 9: triangle at the downleft
        std::cout<<"case 9:"<<std::endl;
        T x,y;
        if( line1(UpRight[0], UpRight[1])   * line1(DownLeft[0],  DownLeft[1]) <=0  &&
            line1(UpRight[0], UpRight[1])   * line1(DownRight[0], DownRight[1])  <=0  &&
            line1(UpRight[0], UpRight[1])   * line1(UpLeft[0],    UpLeft[1])   <=0 )
            {x = line1_getX(UpRight[1]); y = line1_getY(UpRight[0]);}
        else if( line2(UpRight[0], UpRight[1])   * line2(DownLeft[0],  DownLeft[1]) <=0  &&
            line2(UpRight[0], UpRight[1])   * line2(DownRight[0], DownRight[1])  <=0  &&
            line2(UpRight[0], UpRight[1])   * line2(UpLeft[0],    UpLeft[1])   <=0 )
            {x = line2_getX(UpRight[1]); y = line2_getY(UpRight[0]);}
        else if( line3(UpRight[0], UpRight[1])   * line3(DownLeft[0],  DownLeft[1]) <=0  &&
            line3(UpRight[0], UpRight[1])   * line3(DownRight[0], DownRight[1])  <=0  &&
            line3(UpRight[0], UpRight[1])   * line3(UpLeft[0],    UpLeft[1])   <=0 )
            {x = line3_getX(UpRight[1]); y = line3_getY(UpRight[0]);}
        else if( line4(UpRight[0], UpRight[1])   * line4(DownLeft[0],  DownLeft[1]) <=0  &&
            line4(UpRight[0], UpRight[1])   * line4(DownRight[0], DownRight[1])  <=0  &&
            line4(UpRight[0], UpRight[1])   * line4(UpLeft[0],    UpLeft[1])   <=0 )
            {x = line4_getX(UpRight[1]); y = line4_getY(UpRight[0]);}
        std::cout<<"x="<<x<<",y="<<y<<std::endl;
        rightSideLength = UpRight[1] - y;
        upSideLength    = UpRight[0] - x;
        solidFrac = rightSideLength*upSideLength/2;
    }
    else if(inside(UpLeft[0],   UpLeft[1])   == TRUE  &&
            inside(UpRight[0],  UpRight[1])  == TRUE  &&
            inside(DownLeft[0], DownLeft[1]) == FALSE &&
            inside(DownRight[0],DownRight[1])== TRUE){
        //case 10: anti-triangle at the downleft
        std::cout<<"case 10:"<<std::endl;
        T x,y;
        if( line1(DownLeft[0], DownLeft[1])   * line1(UpLeft[0],    UpLeft[1]) <=0  &&
            line1(DownLeft[0], DownLeft[1])   * line1(UpRight[0],   UpRight[1])  <=0  &&
            line1(DownLeft[0], DownLeft[1])   * line1(DownRight[0], DownRight[1])   <=0 )
            {x = line1_getX(DownLeft[1]); y = line1_getY(DownLeft[0]);}
        else if( line2(DownLeft[0], DownLeft[1])   * line2(UpLeft[0],    UpLeft[1]) <=0  &&
            line2(DownLeft[0], DownLeft[1])   * line2(UpRight[0],   UpRight[1])  <=0  &&
            line2(DownLeft[0], DownLeft[1])   * line2(DownRight[0], DownRight[1])   <=0 )
            {x = line2_getX(DownLeft[1]); y = line2_getY(DownLeft[0]);}
        else if( line3(DownLeft[0], DownLeft[1])   * line3(UpLeft[0],    UpLeft[1]) <=0  &&
            line3(DownLeft[0], DownLeft[1])   * line3(UpRight[0],   UpRight[1])  <=0  &&
            line3(DownLeft[0], DownLeft[1])   * line3(DownRight[0], DownRight[1])   <=0 )
            {x = line3_getX(DownLeft[1]); y = line3_getY(DownLeft[0]);}
        else if( line4(DownLeft[0], DownLeft[1])   * line4(UpLeft[0],    UpLeft[1]) <=0  &&
            line4(DownLeft[0], DownLeft[1])   * line4(UpRight[0],   UpRight[1])  <=0  &&
            line4(DownLeft[0], DownLeft[1])   * line4(DownRight[0], DownRight[1])   <=0 )
            {x = line4_getX(DownLeft[1]); y = line4_getY(DownLeft[0]);}
        std::cout<<"x="<<x<<",y="<<y<<std::endl;
        leftSideLength = y - DownLeft[1];
        downSideLength = x - DownLeft[0];
        solidFrac = 1-(1-leftSideLength)*(1-downSideLength)/2;
    }
    else if(inside(UpLeft[0],   UpLeft[1])   == TRUE  &&
            inside(UpRight[0],  UpRight[1])  == FALSE  &&
            inside(DownLeft[0], DownLeft[1]) == FALSE &&
            inside(DownRight[0],DownRight[1])== FALSE){
        //case 11: triangle at the downright
        std::cout<<"case 11:"<<std::endl;
        T x,y;
        if( line1(UpLeft[0], UpLeft[1])   * line1(DownLeft[0],  DownLeft[1]) <=0  &&
            line1(UpLeft[0], UpLeft[1])   * line1(UpRight[0],   UpRight[1])  <=0  &&
            line1(UpLeft[0], UpLeft[1])   * line1(DownRight[0], DownRight[1])   <=0 )
            {x = line1_getX(UpLeft[1]); y = line1_getY(UpLeft[0]);}
        else if( line2(UpLeft[0], UpLeft[1])   * line2(DownLeft[0],  DownLeft[1]) <=0  &&
            line2(UpLeft[0], UpLeft[1])   * line2(UpRight[0],   UpRight[1])  <=0  &&
            line2(UpLeft[0], UpLeft[1])   * line2(DownRight[0], DownRight[1])   <=0 )
            {x = line2_getX(UpLeft[1]); y = line2_getY(UpLeft[0]);}
        else if( line3(UpLeft[0], UpLeft[1])   * line3(DownLeft[0],  DownLeft[1]) <=0  &&
            line3(UpLeft[0], UpLeft[1])   * line3(UpRight[0],   UpRight[1])  <=0  &&
            line3(UpLeft[0], UpLeft[1])   * line3(DownRight[0], DownRight[1])   <=0 )
            {x = line3_getX(UpLeft[1]); y = line3_getY(UpLeft[0]);}
        else if( line4(UpLeft[0], UpLeft[1])   * line4(DownLeft[0],  DownLeft[1]) <=0  &&
            line4(UpLeft[0], UpLeft[1])   * line4(UpRight[0],   UpRight[1])  <=0  &&
            line4(UpLeft[0], UpLeft[1])   * line4(DownRight[0], DownRight[1])   <=0 )
            {x = line4_getX(UpLeft[1]); y = line4_getY(UpLeft[0]);}
        std::cout<<"x="<<x<<",y="<<y<<std::endl;
        leftSideLength = UpLeft[1] - y;
        upSideLength =   x - UpLeft[1];
        solidFrac = leftSideLength*upSideLength/2;
    }
    else if(inside(UpLeft[0],   UpLeft[1])   == TRUE  &&
            inside(UpRight[0],  UpRight[1])  == TRUE  &&
            inside(DownLeft[0], DownLeft[1]) == TRUE  &&
            inside(DownRight[0],DownRight[1])== FALSE){
        //case 12: anti-triangle at the downright
        std::cout<<"case 12:"<<std::endl;
        T x,y;
        if( line1(DownRight[0], DownRight[1])   * line1(DownLeft[0],  DownLeft[1]) <=0  &&
            line1(DownRight[0], DownRight[1])   * line1(UpRight[0],   UpRight[1])  <=0  &&
            line1(DownRight[0], DownRight[1])   * line1(UpLeft[0],    UpLeft[1])   <=0 )
            {x = line1_getX(DownRight[1]); y = line1_getY(DownRight[0]);}
        else if( line2(DownRight[0], DownRight[1])   * line2(DownLeft[0],  DownLeft[1]) <=0  &&
            line2(DownRight[0], DownRight[1])   * line2(UpRight[0],   UpRight[1])  <=0  &&
            line2(DownRight[0], DownRight[1])   * line2(UpLeft[0],    UpLeft[1])   <=0 )
            {x = line2_getX(DownRight[1]); y = line2_getY(DownRight[0]);}
        else if( line3(DownRight[0], DownRight[1])   * line3(DownLeft[0],  DownLeft[1]) <=0  &&
            line3(DownRight[0], DownRight[1])   * line3(UpRight[0],   UpRight[1])  <=0  &&
            line3(DownRight[0], DownRight[1])   * line3(UpLeft[0],    UpLeft[1])   <=0 )
            {x = line3_getX(DownRight[1]); y = line3_getY(DownRight[0]);}
        else if( line4(DownRight[0], DownRight[1])   * line4(DownLeft[0],  DownLeft[1]) <=0  &&
            line4(DownRight[0], DownRight[1])   * line4(UpRight[0],   UpRight[1])  <=0  &&
            line4(DownRight[0], DownRight[1])   * line4(UpLeft[0],    UpLeft[1])   <=0 )
            {x = line4_getX(DownRight[1]); y = line4_getY(DownRight[0]);}
        std::cout<<"x="<<x<<",y="<<y<<std::endl;
        rightSideLength = y - DownRight[1];
        downSideLength =  x - DownRight[0];
        solidFrac = 1-(1-rightSideLength)*(1-downSideLength)/2;
    }

    else if(inside(UpLeft[0],   UpLeft[1])   == TRUE  &&
            inside(UpRight[0],  UpRight[1])  == TRUE  &&
            inside(DownLeft[0], DownLeft[1]) == TRUE &&
            inside(DownRight[0],DownRight[1])== TRUE) {
        //case 13: internal
        std::cout<<"case 13:"<<std::endl;
        solidFrac = 1;
    }

    else if(inside(UpLeft[0],   UpLeft[1])   == FALSE  &&
            inside(UpRight[0],  UpRight[1])  == FALSE  &&
            inside(DownLeft[0], DownLeft[1]) == FALSE &&
            inside(DownRight[0],DownRight[1])== FALSE){
        //case 14: external
        std::cout<<"case 14:"<<std::endl;
        solidFrac = 0;
    }

    else{
        std::cout<<"not included in the above 14 cases."<<std::endl;
    }
    std::cout<<"solidFrac="<<solidFrac<<std::endl;
    system("PAUSE");
}
*/

template<typename T, template<typename U> class Descriptor>
bool SquareParticle2D<T,Descriptor>::line1(T x, T y){
    T theta = AngularPosition;
    T x0 = this->position[0];
    T y0 = this->position[1];
    T l = length;
    return (y - (y0 + l/2*cos(theta) ) >  tan(theta)  *(x - (x0 - l/2*sin(theta) ) ) );
}

template<typename T, template<typename U> class Descriptor>
T SquareParticle2D<T,Descriptor>::line1_getY(T x){
    T theta = AngularPosition;
    T x0 = this->position[0];
    T y0 = this->position[1];
    T l = length;
    return (y0 + l/2*cos(theta) ) +  tan(theta)  *(x - (x0 - l/2*sin(theta) ) );
}

template<typename T, template<typename U> class Descriptor>
T SquareParticle2D<T,Descriptor>::line1_getX(T y){
    T theta = AngularPosition;
    T x0 = this->position[0];
    T y0 = this->position[1];
    T l = length;
    return (x0 - l/2*sin(theta) )  + 1/tan(theta)*(y - (y0 + l/2*cos(theta) ) );
}

template<typename T, template<typename U> class Descriptor>
bool SquareParticle2D<T,Descriptor>::line2(T x, T y){
    T theta = AngularPosition;
    T x0 = this->position[0];
    T y0 = this->position[1];
    T l = length;
    return (y - (y0 - l/2*cos(theta) ) >  tan(theta)  *(x - (x0 + l/2*sin(theta) ) ) );
}

template<typename T, template<typename U> class Descriptor>
T SquareParticle2D<T,Descriptor>::line2_getY(T x){
    T theta = AngularPosition;
    T x0 = this->position[0];
    T y0 = this->position[1];
    T l = length;
    return (y0 - l/2*cos(theta) ) +  tan(theta)  *(x - (x0 + l/2*sin(theta) ) );
}

template<typename T, template<typename U> class Descriptor>
T SquareParticle2D<T,Descriptor>::line2_getX(T y){
    T theta = AngularPosition;
    T x0 = this->position[0];
    T y0 = this->position[1];
    T l = length;
    return (x0 + l/2*sin(theta) )  + 1/tan(theta)*(y - (y0 - l/2*cos(theta) ) );
}

template<typename T, template<typename U> class Descriptor>
bool SquareParticle2D<T,Descriptor>::line3(T x, T y){
    T theta = AngularPosition;
    T x0 = this->position[0];
    T y0 = this->position[1];
    T l = length;
    return (y - (y0 + l/2*sin(theta) ) > -1/tan(theta)*(x - (x0 + l/2*cos(theta) ) ) );
}

template<typename T, template<typename U> class Descriptor>
T SquareParticle2D<T,Descriptor>::line3_getY(T x){
    T theta = AngularPosition;
    T x0 = this->position[0];
    T y0 = this->position[1];
    T l = length;
    return (y0 + l/2*sin(theta) ) -  1/tan(theta)*(x - (x0 + l/2*cos(theta) ) );
}

template<typename T, template<typename U> class Descriptor>
T SquareParticle2D<T,Descriptor>::line3_getX(T y){
    T theta = AngularPosition;
    T x0 = this->position[0];
    T y0 = this->position[1];
    T l = length;
    return (x0 + l/2*cos(theta) )  -   tan(theta)*(y - (y0 + l/2*sin(theta) ) );
}

template<typename T, template<typename U> class Descriptor>
bool SquareParticle2D<T,Descriptor>::line4(T x, T y){
    T theta = AngularPosition;
    T x0 = this->position[0];
    T y0 = this->position[1];
    T l = length;
    return (y - (y0 - 1/2*sin(theta) ) > -1/tan(theta)*(x - (x0 - l/2*cos(theta) ) ) );
}

template<typename T, template<typename U> class Descriptor>
T SquareParticle2D<T,Descriptor>::line4_getY(T x){
    T theta = AngularPosition;
    T x0 = this->position[0];
    T y0 = this->position[1];
    T l = length;
    return (y0 - 1/2*sin(theta) ) -  1/tan(theta)*(x - (x0 - l/2*cos(theta) ) ) ;
}

template<typename T, template<typename U> class Descriptor>
T SquareParticle2D<T,Descriptor>::line4_getX(T y){
    T theta = AngularPosition;
    T x0 = this->position[0];
    T y0 = this->position[1];
    T l = length;
    return (x0 - l/2*cos(theta) )  -   tan(theta)*(y - (y0 - 1/2*sin(theta) ) );
}

}
#endif

