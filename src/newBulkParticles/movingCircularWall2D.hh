#ifndef MOVING_CIRCULAR_WALL_2D_HH
#define MOVING_CIRCULAR_WALL_2D_HH

#include "core/globalDefs.h"
#include "finiteDifference/interpolations2D.h"
#include "particles/particleIdentifiers2D.h"
#include "newBulkParticles/movingCircularWall2D.h"
#include "newBulkParticles/BulkParticleField2D.h"
#include <cmath>
#include <vector>
#define TRUE 1
#define FALSE 0

namespace plb {

/* *************** class movingVerticalWall2D ***************************************** */

const double PI = 3.1415926535897932;
template<typename T, template<typename U> class Descriptor>
movingCircularWall2D<T,Descriptor>::movingCircularWall2D()
{ }

template<typename T, template<typename U> class Descriptor>
movingCircularWall2D<T,Descriptor>::movingCircularWall2D(plint tag_, T x0_, T y0_,
                                                         T innerRadius_, T outerRadius_,
                                                         T angularVelocity_):
      movingObject2D<T,Descriptor>(),
      x0(x0_), y0(y0_), innerRadius(innerRadius_), outerRadius(outerRadius_),
      angularVelocity(angularVelocity_)
{
    this->tag = tag_;
    this->velocity = Array<T,2>(0.0,0.0);
}

template<typename T, template<typename U> class Descriptor>
void movingCircularWall2D<T,Descriptor>::reset(T x0_, T y0_, T innerRadius_, T outerRadius_)
{
    x0 = x0_;
    y0 = y0_;
    innerRadius = innerRadius_;
    outerRadius = outerRadius_;
}

template<typename T, template<typename U> class Descriptor>
void movingCircularWall2D<T,Descriptor>::advance() {
    PLB_ASSERT( norm(velocity)<1. );
}

template<typename T, template<typename U> class Descriptor>
movingObject2D<T,Descriptor>* movingCircularWall2D<T,Descriptor>::clone() const {
    return new movingCircularWall2D<T,Descriptor>(*this);
}


template<typename T, template<typename U> class Descriptor>
void movingCircularWall2D<T,Descriptor>::rescale(int dxScale, int dtScale) {
    int dimDx = 1;
    int dimDt = 0;
    T scaleFactor = scaleFromReference(dxScale, dimDx, dtScale, dimDt);
    x0 *= scaleFactor;
    y0 *= scaleFactor;
    innerRadius *= scaleFactor;
    outerRadius *= scaleFactor;
}

template<typename T, template<typename U> class Descriptor>
void movingCircularWall2D<T,Descriptor>::computeBoundaries(Box2D domain) {
        solidBoundary.clear();
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
            }
        }

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

}


template<typename T, template<typename U> class Descriptor>
void movingCircularWall2D<T,Descriptor>::computeBoundaryLinks(Box2D domain)
{
        this->LeftRightLink.clear();
        this->UpDownLink.clear();
        this->DiagonalLink.clear();
        this->antiDiagonalLink.clear();
        double dx, dy, innerR2, outerR2;
        innerR2 = innerRadius*innerRadius;
        outerR2 = outerRadius*outerRadius;
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                dx = iX - this->x0;
                dy = iY - this->y0;
                if( (dx*dx+dy*dy-innerR2)*((dx+1)*(dx+1)+dy*dy-innerR2)<=0 ||
                    (dx*dx+dy*dy-outerR2)*((dx+1)*(dx+1)+dy*dy-outerR2)<=0 )
                    this->LeftRightLink.push_back(Array<int,4>(iX, iY, iX+1, iY));
                if( (dx*dx+dy*dy-innerR2)*(dx*dx+(dy+1)*(dy+1)-innerR2)<=0 ||
                    (dx*dx+dy*dy-outerR2)*(dx*dx+(dy+1)*(dy+1)-outerR2)<=0 )
                    this->UpDownLink.push_back(Array<int,4>(iX, iY, iX, iY+1));
                if( (dx*dx+dy*dy-innerR2)*((dx+1)*(dx+1)+(dy+1)*(dy+1)-innerR2)<=0 ||
                    (dx*dx+dy*dy-outerR2)*((dx+1)*(dx+1)+(dy+1)*(dy+1)-outerR2)<=0 )
                    this->DiagonalLink.push_back(Array<int,4>(iX, iY, iX+1, iY+1));
                if( (dx*dx+(dy+1)*(dy+1)-innerR2)*((dx+1)*(dx+1)+dy*dy-innerR2)<=0 ||
                    (dx*dx+(dy+1)*(dy+1)-outerR2)*((dx+1)*(dx+1)+dy*dy-outerR2)<=0 )
                    this->antiDiagonalLink.push_back(Array<int,4>(iX+1, iY, iX, iY+1));
            }
        }
}

template<typename T, template<typename U> class Descriptor>
void movingCircularWall2D<T,Descriptor>::objectToFluid(BlockLattice2D<T,Descriptor>& fluid)
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
    for(std::vector<Array<int,4> >::iterator i=this->LeftRightLink.begin();   i!=this->LeftRightLink.end(); ++i){
        tmp1=fluid.get( (*i)[0],(*i)[1] ).f[6];
        tmp2=fluid.get( (*i)[2],(*i)[3] ).f[2];
        dx = (*i)[0] - this->x0;
        dy = (*i)[1] - this->y0;
        if(dx*dx+dy*dy>=innerRadius*innerRadius || dx*dx+dy*dy<=outerRadius*outerRadius){
            r[0] = (T)((*i)[0]) - this->x0;
            r[1] = (T)((*i)[1]) - this->y0;
        }
        else {
            r[0] = (T)((*i)[2]) - this->x0;
            r[1] = (T)((*i)[3]) - this->y0;
        }
        additionalVelocity[0] =-r[1] * angularVelocity;
        additionalVelocity[1] = r[0] * angularVelocity;
        fluid.get( (*i)[0],(*i)[1] ).f[6] = tmp2 + 2./3.*(this->getVelocity()[0]+additionalVelocity[0]);
        fluid.get( (*i)[2],(*i)[3] ).f[2] = tmp1 - 2./3.*(this->getVelocity()[0]+additionalVelocity[0]);
        Mx += 2*(tmp1-tmp2);
        Zxx += 2*2./3.;
    }
    for(std::vector<Array<int,4> >::iterator i=this->UpDownLink.begin();      i!=this->UpDownLink.end();    ++i){
        tmp1=fluid.get( (*i)[0],(*i)[1] ).f[8];
        tmp2=fluid.get( (*i)[2],(*i)[3] ).f[4];
        dx = (*i)[0] - this->x0;
        dy = (*i)[1] - this->y0;
        if(dx*dx+dy*dy>=innerRadius*innerRadius || dx*dx+dy*dy<=outerRadius*outerRadius){
            r[0] = (T)((*i)[0]) - this->x0;
            r[1] = (T)((*i)[1]) - this->y0;
        }
        else {
            r[0] = (T)((*i)[2]) - this->x0;
            r[1] = (T)((*i)[3]) - this->y0;
        }
        additionalVelocity[0] =-r[1] * angularVelocity;
        additionalVelocity[1] = r[0] * angularVelocity;
        fluid.get( (*i)[0],(*i)[1] ).f[8] = tmp2 + 2./3.*(this->getVelocity()[1]+additionalVelocity[1]);
        fluid.get( (*i)[2],(*i)[3] ).f[4] = tmp1 - 2./3.*(this->getVelocity()[1]+additionalVelocity[1]);
        My += 2*(tmp1-tmp2);
        Zyy += 2*2./3.;
    }
    for(std::vector<Array<int,4> >::iterator i=this->DiagonalLink.begin();    i!=this->DiagonalLink.end();  ++i){
        tmp1=fluid.get( (*i)[0],(*i)[1] ).f[7];
        tmp2=fluid.get( (*i)[2],(*i)[3] ).f[3];
        dx = (*i)[0] - this->x0;
        dy = (*i)[1] - this->y0;
        if(dx*dx+dy*dy>=innerRadius*innerRadius || dx*dx+dy*dy<=outerRadius*outerRadius){
            r[0] = (T)((*i)[0]) - this->x0;
            r[1] = (T)((*i)[1]) - this->y0;
        }
        else {
            r[0] = (T)((*i)[2]) - this->x0;
            r[1] = (T)((*i)[3]) - this->y0;
        }
        additionalVelocity[0] =-r[1] * angularVelocity;
        additionalVelocity[1] = r[0] * angularVelocity;
        fluid.get( (*i)[0],(*i)[1] ).f[7] = tmp2 + 1./6.*(this->getVelocity()[0]+additionalVelocity[0] + this->getVelocity()[1]+additionalVelocity[1]);
        fluid.get( (*i)[2],(*i)[3] ).f[3] = tmp1 - 1./6.*(this->getVelocity()[0]+additionalVelocity[0] + this->getVelocity()[1]+additionalVelocity[1]);
        Mx += 2*(tmp1-tmp2);
        My += 2*(tmp1-tmp2);
        Zxx += 2*1./6.;
        Zxy += 2*1./6.;
        Zyy += 2*1./6.;
    }
    for(std::vector<Array<int,4> >::iterator i=this->antiDiagonalLink.begin();i!=this->antiDiagonalLink.end(); ++i){
        tmp1=fluid.get( (*i)[0],(*i)[1] ).f[1];
        tmp2=fluid.get( (*i)[2],(*i)[3] ).f[5];
        dx = (*i)[0] - this->x0;
        dy = (*i)[1] - this->y0;
        if(dx*dx+dy*dy>=innerRadius*innerRadius || dx*dx+dy*dy<=outerRadius*outerRadius){
            r[0] = (T)((*i)[0]) - this->x0;
            r[1] = (T)((*i)[1]) - this->y0;
        }
        else {
            r[0] = (T)((*i)[2]) - this->x0;
            r[1] = (T)((*i)[3]) - this->y0;
        }
        additionalVelocity[0] =-r[1] * angularVelocity;
        additionalVelocity[1] = r[0] * angularVelocity;
        fluid.get( (*i)[0],(*i)[1] ).f[1] = tmp2 + 1./6.*(-this->getVelocity()[0]-additionalVelocity[0] + this->getVelocity()[1]+additionalVelocity[1]);
        fluid.get( (*i)[2],(*i)[3] ).f[5] = tmp1 - 1./6.*(-this->getVelocity()[0]-additionalVelocity[0] + this->getVelocity()[1]+additionalVelocity[1]);
        Mx += -2*(tmp1-tmp2);
        My +=  2*(tmp1-tmp2);
        Zxx += 2*1./6.;
        Zxy += -2*1./6.;
        Zyy += 2*1./6.;
    }
}

template<typename T, template<typename U> class Descriptor>
void movingCircularWall2D<T,Descriptor>::objectToFluidSlip(BlockLattice2D<T,Descriptor>& fluid, T zeta)
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
    printf("LeftRightLink adjustment.\n");

    for(std::vector<Array<int,4> >::iterator i=this->LeftRightLink.begin();   i!=this->LeftRightLink.end(); ++i){
        //right boundary
        dx = (*i)[0] - this->x0;
        dy = (*i)[1] - this->y0;
        if(dx*dx + dy*dy >= outerRadius*outerRadius || dx*dx+dy*dy <= innerRadius*innerRadius ){
            //if(find(    this->DiagonalLink.begin(),     this->DiagonalLink.end(), Array<int,4>((*i)[0],(*i)[1],(*i)[0]+1, (*i)[1]+1)) !=     this->DiagonalLink.end() &&
            //   find(this->antiDiagonalLink.begin(), this->antiDiagonalLink.end(), Array<int,4>((*i)[0]+1,(*i)[1]-1,(*i)[0], (*i)[1])) != this->antiDiagonalLink.end() ){
                dx = (*i)[0] - this->x0;
                dy = (*i)[1] - this->y0;
                if(dx*dx+dy*dy>=innerRadius*innerRadius || dx*dx+dy*dy<=outerRadius*outerRadius){
                    r[0] = (T)((*i)[0]) - this->x0;
                    r[1] = (T)((*i)[1]) - this->y0;
                }
                else {
                    r[0] = (T)((*i)[2]) - this->x0;
                    r[1] = (T)((*i)[3]) - this->y0;
                }
                additionalVelocity[0] =-r[1] * angularVelocity;
                additionalVelocity[1] = r[0] * angularVelocity;
                T vx = this->getVelocity()[0] + additionalVelocity[0];
                T vy = this->getVelocity()[1] + additionalVelocity[1];
                //right boundary is vertical, exchange 1<-->3
                printf("swap 1--3\n");
                //1  8|  7   full-slip : 1<-->7  3<-->5 if added in stream step, since during collide(revert), 1<-->5, 3<-->7. we only need to do 1<-->3
                //2   |  6   non-slip :  1<-->5  3<-->7 if added in stream step, since during collide(revert),...we do not need to do anything
                //3  4|  5
                tmp1=fluid.get( (*i)[0],(*i)[1] ).f[6];
                tmp2=fluid.get( (*i)[2],(*i)[3] ).f[2];
                fluid.get( (*i)[0], (*i)[1] ).f[1] = (1-zeta)*(fluid.get( (*i)[0], (*i)[1] ).f[5] + 1/6*vy);
                fluid.get( (*i)[0], (*i)[1] ).f[1] += zeta*fluid.get( (*i)[0], (*i)[1] ).f[7];
                fluid.get( (*i)[0], (*i)[1] ).f[1] -= 1/6*vx;
                fluid.get( (*i)[0], (*i)[1] ).f[3] = (1-zeta)*(fluid.get( (*i)[0], (*i)[1] ).f[7] - 1/6*vy);
                fluid.get( (*i)[0], (*i)[1] ).f[3] += zeta*fluid.get( (*i)[0], (*i)[1] ).f[5];
                fluid.get( (*i)[0], (*i)[1] ).f[3] -= 1/6*vx;

                fluid.get( (*i)[0], (*i)[1] ).f[2] = fluid.get( (*i)[0],(*i)[1] ).f[6] - 2./3.*vx;
                //swap(fluid.get( (*i)[0],(*i)[1] ).f[1], fluid.get( (*i)[0],(*i)[1] ).f[3]);
                fluid.get( (*i)[2], (*i)[3] ).f[5] = (1-zeta)*(fluid.get( (*i)[2], (*i)[3] ).f[1] - 1/6*vy);
                fluid.get( (*i)[2], (*i)[3] ).f[5] += zeta*fluid.get( (*i)[2], (*i)[3] ).f[3];
                fluid.get( (*i)[2], (*i)[3] ).f[5] += 1/6*vx;
                fluid.get( (*i)[2], (*i)[3] ).f[7] = (1-zeta)*(fluid.get( (*i)[2], (*i)[3] ).f[3] + 1/6*vy);
                fluid.get( (*i)[2], (*i)[3] ).f[7] += zeta*fluid.get( (*i)[2], (*i)[3] ).f[1];
                fluid.get( (*i)[2], (*i)[3] ).f[7] += 1/6*vx;

                fluid.get( (*i)[2], (*i)[3] ).f[6] = fluid.get( (*i)[2],(*i)[3] ).f[2] + 2./3.*vx;
                Mx += 2*(tmp1-tmp2);
                Zxx += 2*2./3.;
            //}
            if(find(    this->DiagonalLink.begin(),     this->DiagonalLink.end(), Array<int,4>((*i)[0],(*i)[1],(*i)[0]+1, (*i)[1]+1)) == this->DiagonalLink.end() &&
                    find(this->antiDiagonalLink.begin(), this->antiDiagonalLink.end(), Array<int,4>((*i)[0]+1,(*i)[1]-1,(*i)[0], (*i)[1])) != this->antiDiagonalLink.end() ){
                //right boundary is like "/", exchange 1<-->2
                printf("swap 1--2\n");
                //swap(fluid.get( (*i)[0],(*i)[1] ).f[1], fluid.get( (*i)[0],(*i)[1] ).f[2]);
            }
            else if(find(    this->DiagonalLink.begin(),     this->DiagonalLink.end(), Array<int,4>((*i)[0],(*i)[1],(*i)[0]+1, (*i)[1]+1)) != this->DiagonalLink.end() &&
                    find(this->antiDiagonalLink.begin(), this->antiDiagonalLink.end(), Array<int,4>((*i)[0]+1,(*i)[1]-1,(*i)[0], (*i)[1])) == this->antiDiagonalLink.end() ){
                //right boundary is like "\", exchange 2<-->3
                printf("swap 2--3\n");
                //swap(fluid.get( (*i)[0],(*i)[1] ).f[2], fluid.get( (*i)[0],(*i)[1] ).f[3]);
            }
            else{}
        }
        //left boundary
        if(dx*dx + dy*dy <= outerRadius*outerRadius && dx*dx + dy*dy >= innerRadius*innerRadius ){
            //if(find(    this->DiagonalLink.begin(),     this->DiagonalLink.end(), Array<int,4>((*i)[2]-1,(*i)[3]-1,(*i)[2], (*i)[3])) !=     this->DiagonalLink.end() &&
            //   find(this->antiDiagonalLink.begin(), this->antiDiagonalLink.end(), Array<int,4>((*i)[2],(*i)[3],(*i)[2]-1, (*i)[3]+1)) != this->antiDiagonalLink.end() ){
                dx = (*i)[0] - this->x0;
                dy = (*i)[1] - this->y0;
                if(dx*dx+dy*dy>=innerRadius*innerRadius || dx*dx+dy*dy<=outerRadius*outerRadius){
                    r[0] = (T)((*i)[0]) - this->x0;
                    r[1] = (T)((*i)[1]) - this->y0;
                }
                else {
                    r[0] = (T)((*i)[2]) - this->x0;
                    r[1] = (T)((*i)[3]) - this->y0;
                }
                additionalVelocity[0] =-r[1] * angularVelocity;
                additionalVelocity[1] = r[0] * angularVelocity;
                T vx = this->getVelocity()[0] + additionalVelocity[0];
                T vy = this->getVelocity()[1] + additionalVelocity[1];
                //left boundary is vertical, exchange 5<-->7
                printf("swap 5--7\n");
                tmp1=fluid.get( (*i)[0],(*i)[1] ).f[6];
                tmp2=fluid.get( (*i)[2],(*i)[3] ).f[2];
                fluid.get( (*i)[0], (*i)[1] ).f[1] = (1-zeta)*(fluid.get( (*i)[0], (*i)[1] ).f[5] + 1/6*vy);
                fluid.get( (*i)[0], (*i)[1] ).f[1] += zeta*fluid.get( (*i)[0], (*i)[1] ).f[7];
                fluid.get( (*i)[0], (*i)[1] ).f[1] -= 1/6*vx;
                fluid.get( (*i)[0], (*i)[1] ).f[3] = (1-zeta)*(fluid.get( (*i)[0], (*i)[1] ).f[7] - 1/6*vy);
                fluid.get( (*i)[0], (*i)[1] ).f[3] += zeta*fluid.get( (*i)[0], (*i)[1] ).f[5];
                fluid.get( (*i)[0], (*i)[1] ).f[3] -= 1/6*vx;

                fluid.get( (*i)[0], (*i)[1] ).f[2] = fluid.get( (*i)[0],(*i)[1] ).f[6] - 2./3.*vx;
                //swap(fluid.get( (*i)[0],(*i)[1] ).f[1], fluid.get( (*i)[0],(*i)[1] ).f[3]);
                fluid.get( (*i)[2], (*i)[3] ).f[5] = (1-zeta)*(fluid.get( (*i)[2], (*i)[3] ).f[1] - 1/6*vy);
                fluid.get( (*i)[2], (*i)[3] ).f[5] += zeta*fluid.get( (*i)[2], (*i)[3] ).f[3];
                fluid.get( (*i)[2], (*i)[3] ).f[5] += 1/6*vx;
                fluid.get( (*i)[2], (*i)[3] ).f[7] = (1-zeta)*(fluid.get( (*i)[2], (*i)[3] ).f[3] + 1/6*vy);
                fluid.get( (*i)[2], (*i)[3] ).f[7] += zeta*fluid.get( (*i)[2], (*i)[3] ).f[1];
                fluid.get( (*i)[2], (*i)[3] ).f[7] += 1/6*vx;

                fluid.get( (*i)[2], (*i)[3] ).f[6] = fluid.get( (*i)[2],(*i)[3] ).f[2] + 2./3.*vx;
                Mx += 2*(tmp1-tmp2);
                Zxx += 2*2./3.;
            //}
            if(find(    this->DiagonalLink.begin(),     this->DiagonalLink.end(), Array<int,4>((*i)[2]-1,(*i)[3]-1,(*i)[2], (*i)[3])) == this->DiagonalLink.end() &&
                find(this->antiDiagonalLink.begin(), this->antiDiagonalLink.end(), Array<int,4>((*i)[2],(*i)[3],(*i)[2]-1, (*i)[3]+1)) != this->antiDiagonalLink.end() ){
                //left boundary is like "/", exchange 5<-->6
                printf("swap 5--6\n");
                //swap(fluid.get( (*i)[2],(*i)[3] ).f[5], fluid.get( (*i)[2],(*i)[3] ).f[6]);
            }
            else if(find(    this->DiagonalLink.begin(),     this->DiagonalLink.end(), Array<int,4>((*i)[2]-1,(*i)[3]-1,(*i)[2], (*i)[3])) != this->DiagonalLink.end() &&
                find(this->antiDiagonalLink.begin(), this->antiDiagonalLink.end(), Array<int,4>((*i)[2],(*i)[3],(*i)[2]-1, (*i)[3]+1)) == this->antiDiagonalLink.end() ){
                //left boundary is like "\", exchange 6<-->7
                printf("swap 6--7\n");
                //swap(fluid.get( (*i)[2],(*i)[3] ).f[6], fluid.get( (*i)[2],(*i)[3] ).f[7]);
            }
            else{}
        }
    }

//==================================================================

    printf("UpDownLink adjustment.\n");
    for(std::vector<Array<int,4> >::iterator i=this->UpDownLink.begin();   i!=this->UpDownLink.end(); ++i){
        //bottom boundary
        dx = (*i)[0] - this->x0;
        dy = (*i)[1] - this->y0;
        if(dx*dx + dy*dy <= outerRadius*outerRadius && dx*dx + dy*dy >= innerRadius*innerRadius ){
            //if(find(    this->DiagonalLink.begin(),     this->DiagonalLink.end(), Array<int,4>((*i)[2]-1,(*i)[3]-1,(*i)[2], (*i)[3])) !=     this->DiagonalLink.end() &&
            //   find(this->antiDiagonalLink.begin(), this->antiDiagonalLink.end(), Array<int,4>((*i)[2]+1,(*i)[3]-1,(*i)[2], (*i)[3])) != this->antiDiagonalLink.end() ){
                dx = (*i)[0] - this->x0;
                dy = (*i)[1] - this->y0;
                if(dx*dx+dy*dy>=innerRadius*innerRadius || dx*dx+dy*dy<=outerRadius*outerRadius){
                    r[0] = (T)((*i)[0]) - this->x0;
                    r[1] = (T)((*i)[1]) - this->y0;
                }
                else {
                    r[0] = (T)((*i)[2]) - this->x0;
                    r[1] = (T)((*i)[3]) - this->y0;
                }
                additionalVelocity[0] =-r[1] * angularVelocity;
                additionalVelocity[1] = r[0] * angularVelocity;
                T vx = this->getVelocity()[0] + additionalVelocity[0];
                T vy = this->getVelocity()[1] + additionalVelocity[1];
                //bottom boundary is horizontal, exchange 1<-->7
                printf("swap 1--7\n");
                tmp1=fluid.get( (*i)[0],(*i)[1] ).f[8];
                tmp2=fluid.get( (*i)[2],(*i)[3] ).f[4];
                fluid.get( (*i)[2], (*i)[3] ).f[1] *= (1-zeta)*(fluid.get( (*i)[2], (*i)[3] ).f[5] - 1/6*vx);
                fluid.get( (*i)[2], (*i)[3] ).f[1] +=     zeta*fluid.get( (*i)[2], (*i)[3] ).f[3];
                fluid.get( (*i)[2], (*i)[3] ).f[1] += 1/6*vy;
                fluid.get( (*i)[2], (*i)[3] ).f[7] *= (1-zeta)*(fluid.get( (*i)[2], (*i)[3] ).f[3] + 1/6*vx);
                fluid.get( (*i)[2], (*i)[3] ).f[7] +=     zeta*fluid.get( (*i)[2], (*i)[3] ).f[5];
                fluid.get( (*i)[2], (*i)[3] ).f[7] += 1/6*vy;

                fluid.get( (*i)[2], (*i)[3] ).f[8] = fluid.get( (*i)[2], (*i)[3] ).f[4] + 2./3.*vy;
                //swap(fluid.get( (*i)[0],(*i)[1] ).f[1], fluid.get( (*i)[0],(*i)[1] ).f[7]);
                fluid.get( (*i)[0], (*i)[1] ).f[3] *= (1-zeta)*(fluid.get( (*i)[0], (*i)[1] ).f[7] - 1/6*vx);
                fluid.get( (*i)[0], (*i)[1] ).f[3] +=     zeta*fluid.get( (*i)[0], (*i)[1] ).f[1];
                fluid.get( (*i)[0], (*i)[1] ).f[3] -= 1/6*vy;
                fluid.get( (*i)[0], (*i)[1] ).f[5] *= (1-zeta)*(fluid.get( (*i)[0], (*i)[1] ).f[1] + 1/6*vx);
                fluid.get( (*i)[0], (*i)[1] ).f[5] +=     zeta*fluid.get( (*i)[0], (*i)[1] ).f[7];
                fluid.get( (*i)[0], (*i)[1] ).f[5] -= 1/6*vy;

                fluid.get( (*i)[0], (*i)[1] ).f[4] = fluid.get( (*i)[0], (*i)[1] ).f[8] - 2./3.*vy;
                My += 2*(tmp1-tmp2);
                Zyy += 2*2./3.;
            //}
            if(find(    this->DiagonalLink.begin(),     this->DiagonalLink.end(), Array<int,4>((*i)[0]-1,(*i)[1]-1,(*i)[0], (*i)[1])) == this->DiagonalLink.end() &&
                    find(this->antiDiagonalLink.begin(), this->antiDiagonalLink.end(), Array<int,4>((*i)[0]+1,(*i)[1]-1,(*i)[0], (*i)[1])) != this->antiDiagonalLink.end() ){
                //bottom boundary is like "", exchange 1<-->8
                printf("swap 1--8\n");
                //swap(fluid.get( (*i)[0],(*i)[1] ).f[1], fluid.get( (*i)[0],(*i)[1] ).f[8]);
            }
            else if(find(    this->DiagonalLink.begin(),     this->DiagonalLink.end(), Array<int,4>((*i)[0]-1,(*i)[1]-1,(*i)[0], (*i)[1])) != this->DiagonalLink.end() &&
                    find(this->antiDiagonalLink.begin(), this->antiDiagonalLink.end(), Array<int,4>((*i)[0]+1,(*i)[1]-1,(*i)[0], (*i)[1])) == this->antiDiagonalLink.end() ){
                //right boundary is like "\", exchange 2<-->3
                printf("swap 7--8\n");
                //swap(fluid.get( (*i)[0],(*i)[1] ).f[7], fluid.get( (*i)[0],(*i)[1] ).f[8]);
            }
            else{}
        }
        //boundary on top
        if(dx*dx + dy*dy >= outerRadius*outerRadius || dx*dx + dy*dy <= innerRadius*innerRadius ){
            //if(find(    this->DiagonalLink.begin(),     this->DiagonalLink.end(), Array<int,4>((*i)[0],(*i)[1],(*i)[0]+1, (*i)[1]+1)) !=     this->DiagonalLink.end() &&
            //   find(this->antiDiagonalLink.begin(), this->antiDiagonalLink.end(), Array<int,4>((*i)[0],(*i)[1],(*i)[0]-1, (*i)[1]+1)) != this->antiDiagonalLink.end() ){
                dx = (*i)[0] - this->x0;
                dy = (*i)[1] - this->y0;
                if(dx*dx+dy*dy>=innerRadius*innerRadius || dx*dx+dy*dy<=outerRadius*outerRadius){
                    r[0] = (T)((*i)[0]) - this->x0;
                    r[1] = (T)((*i)[1]) - this->y0;
                }
                else {
                    r[0] = (T)((*i)[2]) - this->x0;
                    r[1] = (T)((*i)[3]) - this->y0;
                }
                additionalVelocity[0] =-r[1] * angularVelocity;
                additionalVelocity[1] = r[0] * angularVelocity;
                T vx = this->getVelocity()[0] + additionalVelocity[0];
                T vy = this->getVelocity()[1] + additionalVelocity[1];
                printf("swap 3--5\n");
                tmp1=fluid.get( (*i)[0],(*i)[1] ).f[8];
                tmp2=fluid.get( (*i)[2],(*i)[3] ).f[4];
                fluid.get( (*i)[2], (*i)[3] ).f[1] *= (1-zeta)*(fluid.get( (*i)[2], (*i)[3] ).f[5] - 1/6*vx);
                fluid.get( (*i)[2], (*i)[3] ).f[1] +=     zeta*fluid.get( (*i)[2], (*i)[3] ).f[3];
                fluid.get( (*i)[2], (*i)[3] ).f[1] += 1/6*vy;
                fluid.get( (*i)[2], (*i)[3] ).f[7] *= (1-zeta)*(fluid.get( (*i)[2], (*i)[3] ).f[3] + 1/6*vx);
                fluid.get( (*i)[2], (*i)[3] ).f[7] +=     zeta*fluid.get( (*i)[2], (*i)[3] ).f[5];
                fluid.get( (*i)[2], (*i)[3] ).f[7] += 1/6*vy;

                fluid.get( (*i)[2], (*i)[3] ).f[8] = fluid.get( (*i)[2], (*i)[3] ).f[4] + 2./3.*vy;
                //swap(fluid.get( (*i)[0],(*i)[1] ).f[1], fluid.get( (*i)[0],(*i)[1] ).f[7]);
                fluid.get( (*i)[0], (*i)[1] ).f[3] *= (1-zeta)*(fluid.get( (*i)[0], (*i)[1] ).f[7] - 1/6*vx);
                fluid.get( (*i)[0], (*i)[1] ).f[3] +=     zeta*fluid.get( (*i)[0], (*i)[1] ).f[1];
                fluid.get( (*i)[0], (*i)[1] ).f[3] -= 1/6*vy;
                fluid.get( (*i)[0], (*i)[1] ).f[5] *= (1-zeta)*(fluid.get( (*i)[0], (*i)[1] ).f[1] + 1/6*vx);
                fluid.get( (*i)[0], (*i)[1] ).f[5] +=     zeta*fluid.get( (*i)[0], (*i)[1] ).f[7];
                fluid.get( (*i)[0], (*i)[1] ).f[5] -= 1/6*vy;

                fluid.get( (*i)[0], (*i)[1] ).f[4] = fluid.get( (*i)[0], (*i)[1] ).f[8] - 2./3.*vy;
                //swap(fluid.get( (*i)[2],(*i)[3] ).f[3], fluid.get( (*i)[2],(*i)[3] ).f[5]);
                My += 2*(tmp1-tmp2);
                Zyy += 2*2./3.;
            //}
            if(find(    this->DiagonalLink.begin(),     this->DiagonalLink.end(), Array<int,4>((*i)[2],(*i)[3],(*i)[2]+1, (*i)[3]+1)) == this->DiagonalLink.end() &&
                    find(this->antiDiagonalLink.begin(), this->antiDiagonalLink.end(), Array<int,4>((*i)[2],(*i)[3],(*i)[2]-1, (*i)[3]+1)) != this->antiDiagonalLink.end() ){
                printf("swap 4--5\n");
                //swap(fluid.get( (*i)[2],(*i)[3] ).f[4], fluid.get( (*i)[2],(*i)[3] ).f[5]);
            }
            else if(find(    this->DiagonalLink.begin(),     this->DiagonalLink.end(), Array<int,4>((*i)[2],(*i)[3],(*i)[2]+1, (*i)[3]+1)) != this->DiagonalLink.end() &&
                    find(this->antiDiagonalLink.begin(), this->antiDiagonalLink.end(), Array<int,4>((*i)[2],(*i)[3],(*i)[2]-1, (*i)[3]+1)) == this->antiDiagonalLink.end() ){
                //left boundary is like "\", exchange 6<-->7
                printf("swap 3--4\n");
                //swap(fluid.get( (*i)[2],(*i)[3] ).f[3], fluid.get( (*i)[2],(*i)[3] ).f[4]);
            }
            else{}
        }
    }

//=================================================================================
    printf("DiagonalLink adjustment.\n");
/*
    for(std::vector<Array<int,4> >::iterator i=this->DiagonalLink.begin();   i!=this->DiagonalLink.end(); ++i){
        //up-right boundary
        tmp1=fluid.get( (*i)[0],(*i)[1] ).f[7];
        tmp2=fluid.get( (*i)[2],(*i)[3] ).f[3];
        dx = (*i)[0] - this->x0;
        dy = (*i)[1] - this->y0;
        if(dx*dx + dy*dy >= outerRadius*outerRadius || dx*dx + dy*dy <= innerRadius*innerRadius ){
            //if(find(this->LeftRightLink.begin(), this->LeftRightLink.end(), Array<int,4>((*i)[0],(*i)[1],(*i)[0]+1, (*i)[1])) != this->LeftRightLink.end() &&
            //   find(this->UpDownLink.begin(),    this->UpDownLink.end(),    Array<int,4>((*i)[0],(*i)[1]+1,(*i)[0], (*i)[1])) != this->UpDownLink.end() ){
                dx = (*i)[0] - this->x0;
                dy = (*i)[1] - this->y0;
                if(dx*dx+dy*dy>=innerRadius*innerRadius || dx*dx+dy*dy<=outerRadius*outerRadius){
                    r[0] = (T)((*i)[0]) - this->x0;
                    r[1] = (T)((*i)[1]) - this->y0;
                }
                else {
                    r[0] = (T)((*i)[2]) - this->x0;
                    r[1] = (T)((*i)[3]) - this->y0;
                }
                additionalVelocity[0] =-r[1] * angularVelocity;
                additionalVelocity[1] = r[0] * angularVelocity;
                T vx = this->getVelocity()[0] + additionalVelocity[0];
                T vy = this->getVelocity()[1] + additionalVelocity[1];
                printf("swap 2--4\n");
                fluid.get( (*i)[0], (*i)[1] ).f[2] = (1-zeta)*(fluid.get( (*i)[0], (*i)[1] ).f[6] - 2./3.*vx/4);
                fluid.get( (*i)[0], (*i)[1] ).f[2] += zeta*(fluid.get( (*i)[0], (*i)[1] ).f[8] - 1/3*(vx+vy)/4);

                fluid.get( (*i)[0], (*i)[1] ).f[4] = (1-zeta)*(fluid.get( (*i)[0], (*i)[1] ).f[8] - 2./3.*vy/4);
                fluid.get( (*i)[0], (*i)[1] ).f[4] += zeta*(fluid.get( (*i)[0], (*i)[1] ).f[6] - 1/3*(vx+vy)/4);

                fluid.get( (*i)[0], (*i)[1] ).f[3] = fluid.get( (*i)[0], (*i)[1] ).f[7] - 1/6*(vx+vy);

                fluid.get( (*i)[2], (*i)[3] ).f[7] = fluid.get( (*i)[2], (*i)[3] ).f[3] + 1/6*(vx+vy);
            //}
        }
        //bottom-left boundary
        //if the following "if" statement is commented off, solid node of top-right boundary is also considered
        if(dx*dx + dy*dy <= outerRadius*outerRadius && dx*dx + dy*dy >= innerRadius*innerRadius ){
            //if(find(this->LeftRightLink.begin(), this->LeftRightLink.end(), Array<int,4>((*i)[2]-1,(*i)[3],(*i)[2], (*i)[3])) != this->LeftRightLink.end() &&
            //   find(this->UpDownLink.begin(),    this->UpDownLink.end(),    Array<int,4>((*i)[2],(*i)[3],(*i)[2], (*i)[3]-1)) != this->UpDownLink.end() ){
                dx = (*i)[0] - this->x0;
                dy = (*i)[1] - this->y0;
                if(dx*dx+dy*dy>=innerRadius*innerRadius || dx*dx+dy*dy<=outerRadius*outerRadius){
                    r[0] = (T)((*i)[0]) - this->x0;
                    r[1] = (T)((*i)[1]) - this->y0;
                }
                else {
                    r[0] = (T)((*i)[2]) - this->x0;
                    r[1] = (T)((*i)[3]) - this->y0;
                }
                additionalVelocity[0] =-r[1] * angularVelocity;
                additionalVelocity[1] = r[0] * angularVelocity;
                T vx = this->getVelocity()[0] + additionalVelocity[0];
                T vy = this->getVelocity()[1] + additionalVelocity[1];
                printf("swap 6--8\n");
                fluid.get( (*i)[2], (*i)[3] ).f[6] = (1-zeta)*(fluid.get( (*i)[2], (*i)[3] ).f[2] + 2./3.*vx/4);
                fluid.get( (*i)[2], (*i)[3] ).f[6] += zeta*(fluid.get( (*i)[2], (*i)[3] ).f[4] + 1/3*(vx+vy)/4);

                fluid.get( (*i)[2], (*i)[3] ).f[8] = (1-zeta)*(fluid.get( (*i)[2], (*i)[3] ).f[4] + 2./3.*vy/4);
                fluid.get( (*i)[2], (*i)[3] ).f[8] += zeta*(fluid.get( (*i)[2], (*i)[3] ).f[2] + 1/3*(vx+vy)/4);

                fluid.get( (*i)[2], (*i)[3] ).f[7] = fluid.get( (*i)[2], (*i)[3] ).f[3] + 1/6*(vx+vy);

                fluid.get( (*i)[0], (*i)[1] ).f[3] = fluid.get( (*i)[0], (*i)[1] ).f[7] - 1/6*(vx+vy);
            //}
        }
        Mx += 2*(tmp1-tmp2);
        My += 2*(tmp1-tmp2);
        Zxx += 2*1./6.;
        Zxy += 2*1./6.;
        Zyy += 2*1./6.;
    }
*/
//=================================================================================
    printf("antiDiagonalLink adjustment.\n");
/*
    for(std::vector<Array<int,4> >::iterator i=this->antiDiagonalLink.begin();   i!=this->antiDiagonalLink.end(); ++i){
        //up-left boundary
        tmp1=fluid.get( (*i)[0],(*i)[1] ).f[1];
        tmp2=fluid.get( (*i)[2],(*i)[3] ).f[5];
        dx = (*i)[0] - this->x0;
        dy = (*i)[1] - this->y0;
        if(dx*dx + dy*dy >= outerRadius*outerRadius || dx*dx + dy*dy <= innerRadius*innerRadius ){
            //if(find(this->LeftRightLink.begin(), this->LeftRightLink.end(), Array<int,4>((*i)[0]-1,(*i)[1],(*i)[0], (*i)[1])) != this->LeftRightLink.end() &&
            //   find(this->UpDownLink.begin(),    this->UpDownLink.end(),    Array<int,4>((*i)[0],(*i)[1]+1,(*i)[0], (*i)[1])) != this->UpDownLink.end() ){
                dx = (*i)[0] - this->x0;
                dy = (*i)[1] - this->y0;
                if(dx*dx+dy*dy>=innerRadius*innerRadius || dx*dx+dy*dy<=outerRadius*outerRadius){
                    r[0] = (T)((*i)[0]) - this->x0;
                    r[1] = (T)((*i)[1]) - this->y0;
                }
                else {
                    r[0] = (T)((*i)[2]) - this->x0;
                    r[1] = (T)((*i)[3]) - this->y0;
                }
                additionalVelocity[0] =-r[1] * angularVelocity;
                additionalVelocity[1] = r[0] * angularVelocity;
                T vx = this->getVelocity()[0] + additionalVelocity[0];
                T vy = this->getVelocity()[1] + additionalVelocity[1];
                printf("swap 4--6\n");
                fluid.get( (*i)[0], (*i)[1] ).f[4] = (1-zeta)*(fluid.get( (*i)[0], (*i)[1] ).f[8] - 2./3.*vy/4);
                fluid.get( (*i)[0], (*i)[1] ).f[4] += zeta*(fluid.get( (*i)[0], (*i)[1] ).f[2] + 1/3*(vx-vy)/4);

                fluid.get( (*i)[0], (*i)[1] ).f[6] = (1-zeta)*(fluid.get( (*i)[0], (*i)[1] ).f[2] + 2./3.*vx/4);
                fluid.get( (*i)[0], (*i)[1] ).f[6] += zeta*(fluid.get( (*i)[0], (*i)[1] ).f[8] + 1/3*(vx-vy)/4);

                fluid.get( (*i)[0], (*i)[1] ).f[5] = fluid.get( (*i)[0], (*i)[1] ).f[1] + 1/6*(vx-vy);

                fluid.get( (*i)[2], (*i)[3] ).f[1] = fluid.get( (*i)[2], (*i)[3] ).f[5] - 1/6*(vx-vy);
            //}
        }
        //bottom-right boundary
        if(dx*dx + dy*dy <= outerRadius*outerRadius && dx*dx + dy*dy >= innerRadius*innerRadius ){
            //if(find(this->LeftRightLink.begin(), this->LeftRightLink.end(), Array<int,4>((*i)[2],(*i)[3],(*i)[2]+1, (*i)[3])) != this->LeftRightLink.end() &&
            //   find(this->UpDownLink.begin(),    this->UpDownLink.end(),    Array<int,4>((*i)[2],(*i)[3],(*i)[2], (*i)[3]-1)) != this->UpDownLink.end() ){
                dx = (*i)[0] - this->x0;
                dy = (*i)[1] - this->y0;
                if(dx*dx+dy*dy>=innerRadius*innerRadius || dx*dx+dy*dy<=outerRadius*outerRadius){
                    r[0] = (T)((*i)[0]) - this->x0;
                    r[1] = (T)((*i)[1]) - this->y0;
                }
                else {
                    r[0] = (T)((*i)[2]) - this->x0;
                    r[1] = (T)((*i)[3]) - this->y0;
                }
                additionalVelocity[0] =-r[1] * angularVelocity;
                additionalVelocity[1] = r[0] * angularVelocity;
                T vx = this->getVelocity()[0] + additionalVelocity[0];
                T vy = this->getVelocity()[1] + additionalVelocity[1];
                printf("swap 2--8\n");
                fluid.get( (*i)[2], (*i)[3] ).f[2] = (1-zeta)*(fluid.get( (*i)[2], (*i)[3] ).f[6] - 2./3.*vy/4);
                fluid.get( (*i)[2], (*i)[3] ).f[2] += zeta*(fluid.get( (*i)[2], (*i)[3] ).f[4] - 1/3*(vx-vy)/4);

                fluid.get( (*i)[2], (*i)[3] ).f[8] = (1-zeta)*(fluid.get( (*i)[2], (*i)[3] ).f[4] + 2./3.*vy/4);
                fluid.get( (*i)[2], (*i)[3] ).f[8] += zeta*(fluid.get( (*i)[2], (*i)[3] ).f[6] - 1/3*(vx-vy)/4);

                fluid.get( (*i)[2], (*i)[3] ).f[1] = fluid.get( (*i)[2], (*i)[3] ).f[5] - 1/6*(vx - vy);

                fluid.get( (*i)[0], (*i)[1] ).f[5] = fluid.get( (*i)[0], (*i)[1] ).f[1] + 1/6*(vx - vy);
            //}
        }
        Mx += -2*(tmp1-tmp2);
        My +=  2*(tmp1-tmp2);
        Zxx += 2*1./6.;
        Zxy += -2*1./6.;
        Zyy += 2*1./6.;
    }
*/
}

template<typename T, template<typename U> class Descriptor>
void movingCircularWall2D<T,Descriptor>::objectToFluidSlip2(BlockLattice2D<T,Descriptor>& fluid, T zeta){
    double tmp1, tmp2;
    double dx, dy;
    Array<T,2> r;
    Array<T,2> additionalVelocity, tangent_velocity;
    T alpha;
    Array<T,2> slip_velocity;
    for(std::vector<Array<int,2> >::iterator i=this->solidBoundary.begin();   i!=this->solidBoundary.end(); ++i){
        tmp1=fluid.get( (*i)[0],(*i)[1] ).f[6];
        tmp2=fluid.get( (*i)[0],(*i)[1] ).f[2];
        //calculate tangent velocity of solid boundary.
        r[0] = (T)((*i)[0]) - this->x0;
        r[1] = (T)((*i)[1]) - this->y0;
        additionalVelocity[0] =-r[1] * angularVelocity;
        additionalVelocity[1] = r[0] * angularVelocity;
        tangent_velocity[0]=additionalVelocity[0];
        tangent_velocity[1]=additionalVelocity[1];
        //calculate tangent velocity of boundary cells after slip
        if(r[0] > 0.0 )alpha = atan(r[1]/r[0]);
        else if(r[0] < 0.0 && r[1] >= 0.0)alpha = atan(r[1]/r[0]) + PI;
        else if(r[0] < 0.0 && r[1] < 0.0)alpha = atan(r[1]/r[0]) - PI;
        else if(r[0] == 0.0 && r[1] > 0.0)alpha = PI/2;
        else if(r[0] == 0.0 && r[1] < 0.0)alpha = -PI/2;
        else{printf("error"); system("PAUSE");}
        slip_velocity[0] = tangent_velocity[0] - 3*(1-zeta)/zeta*(
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
        slip_velocity[1] = tangent_velocity[1] - 3*(1-zeta)/zeta*(
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
        iniCellAtEquilibrium(fluid.get( (*i)[0],(*i)[1] ), 1.0, slip_velocity);
    }

    for(std::vector<Array<int,2> >::iterator i=this->fluidBoundary.begin();   i!=this->fluidBoundary.end(); ++i){
        tmp1=fluid.get( (*i)[0],(*i)[1] ).f[6];
        tmp2=fluid.get( (*i)[0],(*i)[1] ).f[2];
        //calculate tangent velocity of solid boundary.
        r[0] = (T)((*i)[0]) - this->x0;
        r[1] = (T)((*i)[1]) - this->y0;
        additionalVelocity[0] =-r[1] * angularVelocity;
        additionalVelocity[1] = r[0] * angularVelocity;
        tangent_velocity[0] = additionalVelocity[0];
        tangent_velocity[1] = additionalVelocity[1];
        //calculate tangent velocity of boundary cells after slip
        if(r[0] > 0.0 )alpha = atan(r[1]/r[0]);
        else if(r[0] < 0.0 && r[1] >= 0.0)alpha = atan(r[1]/r[0]) + PI;
        else if(r[0] < 0.0 && r[1] < 0.0)alpha = atan(r[1]/r[0]) - PI;
        else if(r[0] == 0.0 && r[1] > 0.0)alpha = PI/2;
        else if(r[0] == 0.0 && r[1] < 0.0)alpha = -PI/2;
        else{printf("error"); system("PAUSE");}
        slip_velocity[0] = tangent_velocity[0] - 3*(1-zeta)/zeta*(
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
        slip_velocity[1] = tangent_velocity[1] - 3*(1-zeta)/zeta*(
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
        iniCellAtEquilibrium(fluid.get( (*i)[0],(*i)[1] ), 1.0, slip_velocity);
    }
}

}

#endif
