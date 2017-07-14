#ifndef MOVING_STRAIGHT_WALL_2D_HH
#define MOVING_STRAIGHT_WALL_2D_HH

#include "core/globalDefs.h"
#include "finiteDifference/interpolations2D.h"
#include "particles/particleIdentifiers2D.h"
#include "newBulkParticles/movingStraightWall2D.h"
#include "newBulkParticles/BulkParticleField2D.h"
#include <cmath>
#include <vector>
#define TRUE 1
#define FALSE 0

namespace plb {

/* *************** class movingVerticalWall2D ***************************************** */

template<typename T, template<typename U> class Descriptor>
movingVerticalWall2D<T,Descriptor>::movingVerticalWall2D()
{ }

template<typename T, template<typename U> class Descriptor>
movingVerticalWall2D<T,Descriptor>::movingVerticalWall2D(plint tag_, T x0_, T y0_, T y1_, Array<T,2> velocity_):
      movingObject2D<T,Descriptor>(),
      x0(x0_), y0(y0_), y1(y1_)
{
    this->tag = tag_;
    this->velocity = velocity_;
}

template<typename T, template<typename U> class Descriptor>
void movingVerticalWall2D<T,Descriptor>::reset(T x0_, T y0_, T y1_)
{
    x0 = x0_;
    y0 = y0_;
    y1 = y1_;
}

template<typename T, template<typename U> class Descriptor>
void movingVerticalWall2D<T,Descriptor>::advance() {
    PLB_ASSERT( norm(velocity)<1. );
}

template<typename T, template<typename U> class Descriptor>
movingObject2D<T,Descriptor>* movingVerticalWall2D<T,Descriptor>::clone() const {
    return new movingVerticalWall2D<T,Descriptor>(*this);
}


template<typename T, template<typename U> class Descriptor>
void movingVerticalWall2D<T,Descriptor>::rescale(int dxScale, int dtScale) {
    int dimDx = 1;
    int dimDt = 0;
    T scaleFactor = scaleFromReference(dxScale, dimDx, dtScale, dimDt);
    x0 *= scaleFactor;
    y0 *= scaleFactor;
    y1 *= scaleFactor;
}

template<typename T, template<typename U> class Descriptor>
void movingVerticalWall2D<T,Descriptor>::computeBoundaryLinks(Box2D domain)
{
        this->LeftRightLink.clear();
        this->UpDownLink.clear();
        this->DiagonalLink.clear();
        this->antiDiagonalLink.clear();
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if( iY >= y0 && iY <= y1 && (iX-x0)*(iX+1-x0) <=0 ){
                    this->LeftRightLink.push_back(Array<int,4>(iX, iY, iX+1, iY));
                    this->DiagonalLink.push_back(Array<int,4>(iX, iY, iX+1, iY+1));
                    this->antiDiagonalLink.push_back(Array<int,4>(iX+1, iY, iX, iY+1));
                }
            }
        }
}

/* *************** class movingInclinedWall2D ***************************************** */

template<typename T, template<typename U> class Descriptor>
movingInclinedWall2D<T,Descriptor>::movingInclinedWall2D()
{ }

template<typename T, template<typename U> class Descriptor>
movingInclinedWall2D<T,Descriptor>::movingInclinedWall2D(plint tag_, T x0_, T y0_, T y1_, T k_, Array<T,2> velocity_):
      movingObject2D<T,Descriptor>(),
      x0(x0_), y0(y0_), y1(y1_), k(k_)
{
    this->tag = tag_;
    this->velocity = velocity_;
}

template<typename T, template<typename U> class Descriptor>
void movingInclinedWall2D<T,Descriptor>::reset(T x0_, T y0_, T y1_, T k_)
{
    x0 = x0_;
    y0 = y0_;
    y1 = y1_;
    k  = k_;
}

template<typename T, template<typename U> class Descriptor>
void movingInclinedWall2D<T,Descriptor>::advance() {
    PLB_ASSERT( norm(velocity)<1. );
}

template<typename T, template<typename U> class Descriptor>
movingObject2D<T,Descriptor>* movingInclinedWall2D<T,Descriptor>::clone() const {
    return new movingInclinedWall2D<T,Descriptor>(*this);
}


template<typename T, template<typename U> class Descriptor>
void movingInclinedWall2D<T,Descriptor>::rescale(int dxScale, int dtScale) {
    int dimDx = 1;
    int dimDt = 0;
    T scaleFactor = scaleFromReference(dxScale, dimDx, dtScale, dimDt);
    x0 *= scaleFactor;
    y0 *= scaleFactor;
    y1 *= scaleFactor;
}

template<typename T, template<typename U> class Descriptor>
void movingInclinedWall2D<T,Descriptor>::computeBoundaryLinks(Box2D domain)
{
        this->LeftRightLink.clear();
        this->UpDownLink.clear();
        this->DiagonalLink.clear();
        this->antiDiagonalLink.clear();
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if( iY >= y0 && iY <= y1 && ((iY-y0)-k*(iX-x0))*((iY-y0)-k*(iX+1-x0)) <=0 )
                    this->LeftRightLink.push_back(Array<int,4>(iX, iY, iX+1, iY));
                if( iY >= y0 && iY <= y1 && ((iY-y0)-k*(iX-x0))*((iY+1-y0)-k*(iX-x0)) <=0 )
                    this->UpDownLink.push_back(Array<int,4>(iX, iY, iX, iY+1));
                if( iY >= y0 && iY <= y1 && ((iY-y0)-k*(iX-x0))*((iY+1-y0)-k*(iX+1-x0)) <=0 )
                    this->DiagonalLink.push_back(Array<int,4>(iX, iY, iX+1, iY-1));
                if( iY >= y0 && iY <= y1 && ((iY-y0)-k*(iX-x0))*((iY-1-y0)-k*(iX+1-x0)) <=0 )
                    this->antiDiagonalLink.push_back(Array<int,4>(iX, iY, iX+1, iY-1));
            }
        }
}


}
#endif
