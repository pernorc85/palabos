#ifndef MOVING_LUMEN_2D_HH
#define MOVING_LUMEN_2D_HH

#include "core/globalDefs.h"
#include "finiteDifference/interpolations2D.h"
#include "particles/particleIdentifiers2D.h"
#include "newBulkParticles/movingLumen2D.h"
#include "newBulkParticles/BulkParticleField2D.h"
#include <cmath>
#include <vector>
#define TRUE 1
#define FALSE 0

namespace plb {

/* *************** class BulkParticle2D ***************************************** */

template<typename T, template<typename U> class Descriptor>
movingLumen2D<T,Descriptor>::movingLumen2D()
{ }

template<typename T, template<typename U> class Descriptor>
movingLumen2D<T,Descriptor>::movingLumen2D(plint tag_, T x0_, T x1_, T y0_, T y1_, Array<T,2> velocity_):
      movingObject2D<T,Descriptor>(),
      x0(x0_), x1(x1_), y0(y0_), y1(y1_)
{
    this->tag = tag_;
    this->velocity = velocity_;
    h = 0;
}

template<typename T, template<typename U> class Descriptor>
void movingLumen2D<T,Descriptor>::reset(T x0_, T x1_, T y0_, T y1_)
{
    x0 = x0_;
    x1 = x1_;
    y0 = y0_;
    y1 = y1_;
}

template<typename T, template<typename U> class Descriptor>
void movingLumen2D<T,Descriptor>::advance() {
    PLB_ASSERT( norm(velocity)<1. );
    h += this->velocity[1];
}

template<typename T, template<typename U> class Descriptor>
movingObject2D<T,Descriptor>* movingLumen2D<T,Descriptor>::clone() const {
    return new movingLumen2D<T,Descriptor>(*this);
}


template<typename T, template<typename U> class Descriptor>
void movingLumen2D<T,Descriptor>::rescale(int dxScale, int dtScale) {
    int dimDx = 1;
    int dimDt = 0;
    T scaleFactor = scaleFromReference(dxScale, dimDx, dtScale, dimDt);
    x0 *= scaleFactor;
    x1 *= scaleFactor;
    y0 *= scaleFactor;
    y1 *= scaleFactor;
    h  *= scaleFactor;
}

template<typename T, template<typename U> class Descriptor>
void movingLumen2D<T,Descriptor>::computeBoundaryLinks(Box2D domain)
{
        this->LeftRightLink.clear();
        this->UpDownLink.clear();
        this->DiagonalLink.clear();
        this->antiDiagonalLink.clear();
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if( iY >= y0 && iY <= y1 && (iX-x0)*(iX+1-x0) <=0 )
                    this->LeftRightLink.push_back(Array<int,4>(iX, iY, iX+1, iY));
                if( iY >= y0 && iY <= y1 && (iX-x1)*(iX+1-x1) <=0 )
                    this->LeftRightLink.push_back(Array<int,4>(iX, iY, iX+1, iY));
//lower boundary
//b = y0+h
//a = -4*h/(x0-x1)^2
//y = a*(x-(x0+x1)/2)^2 + b
                double b = y0+h;
                double a = -4*h/(x0-x1)/(x0-x1);
                if( iY >= y0-abs(h)-1 && iY <= y1+abs(h) && iX>=x0 && iX<=x1 ){
                    if((iY<a*pow(iX  -(x0+x1)/2,2.0)+b && iY  >a*pow(iX+1-(x0+x1)/2,2.0)+b) ||
                       (iY>a*pow(iX  -(x0+x1)/2,2.0)+b && iY  <a*pow(iX+1-(x0+x1)/2,2.0)+b) )
                        this->LeftRightLink.push_back(Array<int,4>(iX, iY, iX+1, iY));
                    if(iY<a*pow(iX  -(x0+x1)/2,2.0)+b && iY+1>a*pow(iX  -(x0+x1)/2,2.0)+b)
                        this->UpDownLink.push_back(Array<int,4>(iX, iY, iX, iY+1));
                    if(iY<a*pow(iX  -(x0+x1)/2,2.0)+b && iY+1>a*pow(iX+1-(x0+x1)/2,2.0)+b)
                        this->DiagonalLink.push_back(Array<int,4>(iX, iY, iX+1, iY+1));
                    if(iY<a*pow(iX+1-(x0+x1)/2,2.0)+b && iY+1>a*pow(iX  -(x0+x1)/2,2.0)+b)
                        this->antiDiagonalLink.push_back(Array<int,4>(iX+1, iY, iX, iY+1));
                }
//upper boundary
//b = y1+h;
//a = -4*h/(x0-x1)^2
//y = a*(x-(x0+x1)/2)^2 + b
                b = y1+h;
                if( iY >= y0-abs(h)-1 && iY <= y1+abs(h) && iX>=x0 && iX<=x1 ){
                    if((iY<a*pow(iX  -(x0+x1)/2,2.0)+b && iY  >a*pow(iX+1-(x0+x1)/2,2.0)+b) ||
                       (iY>a*pow(iX  -(x0+x1)/2,2.0)+b && iY  <a*pow(iX+1-(x0+x1)/2,2.0)+b) )
                        this->LeftRightLink.push_back(Array<int,4>(iX, iY, iX+1, iY));
                    if(iY<a*pow(iX  -(x0+x1)/2,2.0)+b && iY+1>a*pow(iX  -(x0+x1)/2,2.0)+b)
                        this->UpDownLink.push_back(Array<int,4>(iX, iY, iX, iY+1));
                    if(iY<a*pow(iX  -(x0+x1)/2,2.0)+b && iY+1>a*pow(iX+1-(x0+x1)/2,2.0)+b)
                        this->DiagonalLink.push_back(Array<int,4>(iX, iY, iX+1, iY+1));
                    if(iY<a*pow(iX+1-(x0+x1)/2,2.0)+b && iY+1>a*pow(iX  -(x0+x1)/2,2.0)+b)
                        this->antiDiagonalLink.push_back(Array<int,4>(iX+1, iY, iX, iY+1));
                }

            }
        }
}

template<typename T, template<typename U> class Descriptor>
void movingLumen2D<T,Descriptor>::objectToParticle(BulkParticle2D<T,Descriptor>* particle){
    Array<T,2> Force;

    T epsilon = 10.0;
    std::vector<double> x_min, y_min;
    Force[0] = 0, Force[1] = 0;
    //calculate the shortest distance between the particle
    //and the boundary of movingObject
//lower boundary
//b = y0+h
//a = -4*h/(x0-x1)^2
//y = f(x) = a*(x-(x0+x1)/2)^2 + b
//f'(x) = 2a*(x-(x0+x1)/2)
//solve for:
//dD_sq/dx = 2*(x-particleX) + 2*(f(x)-particleY)*f'(x)
//get (x_min, y_min)

    double b = y0+h;
    double a = -4*h/(x0-x1)/(x0-x1);
    x_min.push_back(particle->getPosition()[0]);//trial value
    y_min.push_back(a*pow(particle->getPosition()[0]-(x0+x1)/2,2.0) + b);
//upper boundary
    b = y1+h;
    x_min.push_back(particle->getPosition()[0]);
    y_min.push_back(a*pow(particle->getPosition()[0]-(x0+x1)/2,2.0) + b);

    for(int i=0;i<x_min.size();i++){
        if(particle->getPosition()[0] >= x0+50 && particle->getPosition()[1] <= x1-50){
            T distance = sqrt(pow(particle->getPosition()[0] - x_min[i],2.0) + pow(particle->getPosition()[1] - y_min[i],2.0) );
            if( distance >= particle->getRadius() + epsilon)
                Force += 0;
            else if(distance < particle->getRadius() + epsilon &&
                    distance > particle->getRadius() ){
                Force[0] += 100*pow(distance - particle->getRadius()-epsilon,4.0) *
                        (particle->getPosition()[0] - x_min[i])/distance;
                Force[1] += 100*pow(distance - particle->getRadius()-epsilon,4.0) *
                        (particle->getPosition()[1] - y_min[i])/distance;
            }
            else if(distance < particle->getRadius() && distance > 0){
                Force[0] += 100*( pow(distance - particle->getRadius() - epsilon, 4.0) + (particle->getRadius()-distance)/epsilon ) *
                        (particle->getPosition()[0] - x_min[i])/distance;
                Force[1] += 100*( pow(distance - particle->getRadius() - epsilon, 4.0) + (particle->getRadius()-distance)/epsilon ) *
                        (particle->getPosition()[1] - y_min[i])/distance;
            }
        }
    }


//accelerate

    double mass = 3.141592654*particle->getRadius()*particle->getRadius()*particle->getDensity();
    particle->setVelocity( Array<T,2>(particle->getVelocity()[0] + Force[0]/mass, particle->getVelocity()[1] + Force[1]/mass) );
}

}
#endif
