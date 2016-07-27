#ifndef MOVING_OVAL_2D_HH
#define MOVING_OVAL_2D_HH

#include "core/globalDefs.h"
#include "finiteDifference/interpolations2D.h"
#include "particles/particleIdentifiers2D.h"
#include "newBulkParticles/movingOval2D.h"
#include "newBulkParticles/BulkParticleField2D.h"
#include <cmath>
#include <vector>
#define TRUE 1
#define FALSE 0

namespace plb {

/* *************** class BulkParticle2D ***************************************** */

template<typename T, template<typename U> class Descriptor>
movingOval2D<T,Descriptor>::movingOval2D()
{ }

template<typename T, template<typename U> class Descriptor>
movingOval2D<T,Descriptor>::movingOval2D(plint tag_, T x0_, T x1_, T y0_, T y1_, Array<T,2> const& velocity_):
      movingObject2D<T,Descriptor>(),
      x0(x0_), x1(x1_), y0(y0_), y1(y1_)
{
    this->tag = tag_;
    this->velocity = velocity_;
}

template<typename T, template<typename U> class Descriptor>
void movingOval2D<T,Descriptor>::reset(T x0_, T x1_, T y0_, T y1_)
{
    x0 = x0_;
    x1 = x1_;
    y0 = y0_;
    y1 = y1_;
}

template<typename T, template<typename U> class Descriptor>
void movingOval2D<T,Descriptor>::advance() {
    PLB_ASSERT( norm(velocity)<1. );
    x0 += this->velocity[0];
    x1 += this->velocity[0];
    y0 += this->velocity[1];
    y1 += this->velocity[1];
}

template<typename T, template<typename U> class Descriptor>
movingObject2D<T,Descriptor>* movingOval2D<T,Descriptor>::clone() const {
    return new movingOval2D<T,Descriptor>(*this);
}


template<typename T, template<typename U> class Descriptor>
void movingOval2D<T,Descriptor>::rescale(int dxScale, int dtScale) {
    int dimDx = 1;
    int dimDt = 0;
    T scaleFactor = scaleFromReference(dxScale, dimDx, dtScale, dimDt);
    x0 *= scaleFactor;
    x1 *= scaleFactor;
    y0 *= scaleFactor;
    y1 *= scaleFactor;
}

template<typename T, template<typename U> class Descriptor>
void movingOval2D<T,Descriptor>::computeBoundaryLinks(Box2D domain)
{
        this->LeftRightLink.clear();
        this->UpDownLink.clear();
        this->DiagonalLink.clear();
        this->antiDiagonalLink.clear();
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {                    
                if( iX >= x0 && iX <= x1 && (iY-y0)*(iY+1-y0) <=0 )
                    this->UpDownLink.push_back(Array<int,4>(iX, iY, iX, iY+1));
                if( iX >= x0 && iX <= x1 && (iY-y1)*(iY+1-y1) <=0 )
                    this->UpDownLink.push_back(Array<int,4>(iX, iY, iX, iY+1));
                    
                if( iY >= y0-1 && iY <= y1 && iX<=x0 && iX>=x0-(y1-y0)/2-1 ){
                    if((pow(iX+1-x0,2.0) + pow(iY-  (y0+y1)/2,2.0) - pow((y1-y0)/2,2.0)) *
                       (pow(iX  -x0,2.0) + pow(iY-  (y0+y1)/2,2.0) - pow((y1-y0)/2,2.0)) <=0 )
                        this->LeftRightLink.push_back(Array<int,4>(iX, iY, iX+1, iY));
                    if((pow(iX  -x0,2.0) + pow(iY-  (y0+y1)/2,2.0) - pow((y1-y0)/2,2.0)) *
                       (pow(iX  -x0,2.0) + pow(iY+1-(y0+y1)/2,2.0) - pow((y1-y0)/2,2.0)) <=0 && iX<x0)
                        this->UpDownLink.push_back(   Array<int,4>(iX, iY, iX, iY+1));
                    if((pow(iX+1-x0,2.0) + pow(iY+1-(y0+y1)/2,2.0) - pow((y1-y0)/2,2.0)) *
                       (pow(iX  -x0,2.0) + pow(iY  -(y0+y1)/2,2.0) - pow((y1-y0)/2,2.0)) <=0 )
                        this->DiagonalLink.push_back(Array<int,4>(iX, iY, iX+1, iY+1));
                    if((pow(iX+1-x0,2.0) + pow(iY  -(y0+y1)/2,2.0) - pow((y1-y0)/2,2.0)) *
                       (pow(iX  -x0,2.0) + pow(iY+1-(y0+y1)/2,2.0) - pow((y1-y0)/2,2.0)) <=0 )
                        this->antiDiagonalLink.push_back(Array<int,4>(iX+1, iY, iX, iY+1));
                }
                if( iY >= y0-1 && iY <= y1 && iX>=x1-1 && iX<=x1+(y1-y0)/2 ){
                    if((pow(iX+1-x1,2.0) + pow(iY  -(y0+y1)/2,2.0) - pow((y1-y0)/2,2.0)) *
                       (pow(iX  -x1,2.0) + pow(iY  -(y0+y1)/2,2.0) - pow((y1-y0)/2,2.0)) <=0 )
                        this->LeftRightLink.push_back(Array<int,4>(iX, iY, iX+1, iY));
                    if((pow(iX  -x1,2.0) + pow(iY  -(y0+y1)/2,2.0) - pow((y1-y0)/2,2.0)) *
                       (pow(iX  -x1,2.0) + pow(iY+1-(y0+y1)/2,2.0) - pow((y1-y0)/2,2.0)) <=0 && iX>x1)
                        this->UpDownLink.push_back(   Array<int,4>(iX, iY, iX, iY+1));
                    if((pow(iX+1-x1,2.0) + pow(iY+1-(y0+y1)/2,2.0) - pow((y1-y0)/2,2.0)) *
                       (pow(iX  -x1,2.0) + pow(iY  -(y0+y1)/2,2.0) - pow((y1-y0)/2,2.0)) <=0 )
                        this->DiagonalLink.push_back(Array<int,4>(iX, iY, iX+1, iY+1));
                    if((pow(iX+1-x1,2.0) + pow(iY  -(y0+y1)/2,2.0) - pow((y1-y0)/2,2.0)) *
                       (pow(iX  -x1,2.0) + pow(iY+1-(y0+y1)/2,2.0) - pow((y1-y0)/2,2.0)) <=0 )
                        this->antiDiagonalLink.push_back(Array<int,4>(iX+1, iY, iX, iY+1));
                }
            }
        }
}

}
#endif
