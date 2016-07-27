#ifndef MOVING_OBJECT_2D_HH
#define MOVING_OBJECT_2D_HH

#include "core/globalDefs.h"
#include "finiteDifference/interpolations2D.h"
#include "particles/particleIdentifiers2D.h"
#include "newBulkParticles/movingObject2D.h"
#include "newBulkParticles/solidBlockLattice2D.h"
#include "newBulkParticles/BulkParticleField2D.h"
#include <cmath>
#include <vector>
#define TRUE 1
#define FALSE 0

namespace plb {

/* *************** class BulkParticle2D ***************************************** */

template<typename T, template<typename U> class Descriptor>
movingObject2D<T,Descriptor>::movingObject2D()
{ }

template<typename T, template<typename U> class Descriptor>
movingObject2D<T,Descriptor>::movingObject2D(plint tag_, Array<T,2> const& velocity_)
    : tag(tag_), velocity(velocity_)
{ }


template<typename T, template<typename U> class Descriptor>
movingObject2D<T,Descriptor>* movingObject2D<T,Descriptor>::clone() const {
    return new movingObject2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void movingObject2D<T,Descriptor>::objectToFluid(BlockLattice2D<T,Descriptor>& fluid)
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
        fluid.get( (*i)[0],(*i)[1] ).f[6] = tmp2 + 2./3.*(this->getVelocity()[0]);
        fluid.get( (*i)[2],(*i)[3] ).f[2] = tmp1 - 2./3.*(this->getVelocity()[0]);
        Mx += 2*(tmp2-tmp1);
        Zxx += 2*2./3.;
    }
    for(std::vector<Array<int,4> >::iterator i=UpDownLink.begin();      i!=UpDownLink.end();    ++i){
        tmp1=fluid.get( (*i)[0],(*i)[1] ).f[8];
        tmp2=fluid.get( (*i)[2],(*i)[3] ).f[4];
        fluid.get( (*i)[0],(*i)[1] ).f[8] = tmp2 + 2./3.*(this->getVelocity()[1]);
        fluid.get( (*i)[2],(*i)[3] ).f[4] = tmp1 - 2./3.*(this->getVelocity()[1]);
        My += 2*(tmp2-tmp1);
        Zyy += 2*2./3.;
    }
    for(std::vector<Array<int,4> >::iterator i=DiagonalLink.begin();    i!=DiagonalLink.end();  ++i){
        tmp1=fluid.get( (*i)[0],(*i)[1] ).f[7];
        tmp2=fluid.get( (*i)[2],(*i)[3] ).f[3];
        fluid.get( (*i)[0],(*i)[1] ).f[7] = tmp2 + 1./6.*(this->getVelocity()[0] + this->getVelocity()[1]);
        fluid.get( (*i)[2],(*i)[3] ).f[3] = tmp1 - 1./6.*(this->getVelocity()[0] + this->getVelocity()[1]);
        Mx += 2*(tmp2-tmp1);
        My += 2*(tmp2-tmp1);
        Zxx += 2*1./6.;
        Zxy += 2*1./6.;
        Zyy += 2*1./6.;
    }
    for(std::vector<Array<int,4> >::iterator i=antiDiagonalLink.begin();i!=antiDiagonalLink.end(); ++i){
        tmp1=fluid.get( (*i)[0],(*i)[1] ).f[1];
        tmp2=fluid.get( (*i)[2],(*i)[3] ).f[5];
        fluid.get( (*i)[0],(*i)[1] ).f[1] = tmp2 + 1./6.*(-this->getVelocity()[0] + this->getVelocity()[1]);
        fluid.get( (*i)[2],(*i)[3] ).f[5] = tmp1 - 1./6.*(-this->getVelocity()[0] + this->getVelocity()[1]);
        Mx += -2*(tmp2-tmp1);
        My +=  2*(tmp2-tmp1);
        Zxx += 2*1./6.;
        Zxy += -2*1./6.;
        Zyy += 2*1./6.;
    }

//    double zxx, zxy, zyy;
//    double alp = 1;
//    double M = (x1-x0)*(y1-y0);
//    zxx = (M+alp*Zyy)/((M+alp*Zxx)*(M+alp*Zyy)+alp*alp*Zxy*Zxy);
//    zxy =    alp*Zxy /((M+alp*Zxx)*(M+alp*Zyy)+alp*alp*Zxy*Zxy);
//    zyy = (M+alp*Zxx)/((M+alp*Zxx)*(M+alp*Zyy)+alp*alp*Zxy*Zxy);

//        double tmpX = velocity[0], tmpY = velocity[1];
//        velocity[0] += zxx*(Mx-Zxx*tmpX-Zxy*tmpY)+zxy*(My-Zxy*tmpX-Zyy*tmpY);
//        velocity[1] += zxy*(Mx-Zxx*tmpX-Zxy*tmpY)+zyy*(My-Zxy*tmpX-Zyy*tmpY);
}
/*
template<typename T, template<typename U> class Descriptor>
void movingObject2D<T,Descriptor>::objectToFluidSlip(BlockLattice2D<T,Descriptor>& fluid, T zeta)
{
//    std::vector<Array<int,4> > LeftRightLink;
//    std::vector<Array<int,4> > UpDownLink;
//    std::vector<Array<int,4> > DiagonalLink;
//    std::vector<Array<int,4> > antiDiagonalLink;

    printf("LeftRightLink adjustment.\n");
    for(std::vector<Array<int,4> >::iterator i=LeftRightLink.begin();   i!=LeftRightLink.end(); ++i){
        //right boundary
        if((*i)[0] < this->x0,2.0 ){
            if(find(    DiagonalLink.begin(),     DiagonalLink.end(), Array<int,4>((*i)[0],(*i)[1],(*i)[0]+1, (*i)[1]+1)) !=     DiagonalLink.end() &&
               find(antiDiagonalLink.begin(), antiDiagonalLink.end(), Array<int,4>((*i)[0]+1,(*i)[1]-1,(*i)[0], (*i)[1])) != antiDiagonalLink.end() ){
                //right boundary is vertical, exchange 1<-->3
                printf("swap 1--3\n");
                fluid.get( (*i)[0], (*i)[1] ).f[1] *= (1-zeta);
                fluid.get( (*i)[0], (*i)[1] ).f[1] += zeta*(fluid.get( (*i)[0], (*i)[1] ).f[7] - 1/6*this->getVelocity()[0]);
                fluid.get( (*i)[0], (*i)[1] ).f[3] *= (1-zeta);
                fluid.get( (*i)[0], (*i)[1] ).f[3] += zeta*(fluid.get( (*i)[0], (*i)[1] ).f[5] - 1/6*this->getVelocity()[0]);
                //swap(fluid.get( (*i)[0],(*i)[1] ).f[1], fluid.get( (*i)[0],(*i)[1] ).f[3]);
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
                fluid.get( (*i)[0], (*i)[1] ).f[5] *= (1-zeta);
                fluid.get( (*i)[0], (*i)[1] ).f[5] += zeta*(fluid.get( (*i)[0], (*i)[1] ).f[3] + 1/6*this->getVelocity()[0]);
                fluid.get( (*i)[0], (*i)[1] ).f[7] *= (1-zeta);
                fluid.get( (*i)[0], (*i)[1] ).f[7] += zeta*(fluid.get( (*i)[0], (*i)[1] ).f[1] + 1/6*this->getVelocity()[0]);
                //swap(fluid.get( (*i)[2],(*i)[3] ).f[5], fluid.get( (*i)[2],(*i)[3] ).f[7]);
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
            if(find(    DiagonalLink.begin(),     DiagonalLink.end(), Array<int,4>((*i)[0]-1,(*i)[1]-1,(*i)[0], (*i)[1])) !=     DiagonalLink.end() &&
               find(antiDiagonalLink.begin(), antiDiagonalLink.end(), Array<int,4>((*i)[0]+1,(*i)[1]-1,(*i)[0], (*i)[1])) != antiDiagonalLink.end() ){
                //bottom boundary is horizontal, exchange 1<-->7
                printf("swap 1--7\n");
                fluid.get( (*i)[0], (*i)[1] ).f[1] *= (1-zeta);
                fluid.get( (*i)[0], (*i)[1] ).f[1] += zeta*(fluid.get( (*i)[0], (*i)[1] ).f[3] + 1/6*this->getVelocity()[1]);
                fluid.get( (*i)[0], (*i)[1] ).f[7] *= (1-zeta);
                fluid.get( (*i)[0], (*i)[1] ).f[7] += zeta*(fluid.get( (*i)[0], (*i)[1] ).f[5] + 1/6*this->getVelocity()[1]);
                //swap(fluid.get( (*i)[0],(*i)[1] ).f[1], fluid.get( (*i)[0],(*i)[1] ).f[7]);
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
            if(find(    DiagonalLink.begin(),     DiagonalLink.end(), Array<int,4>((*i)[2],(*i)[3],(*i)[2]+1, (*i)[3]+1)) !=     DiagonalLink.end() &&
               find(antiDiagonalLink.begin(), antiDiagonalLink.end(), Array<int,4>((*i)[2],(*i)[3],(*i)[2]-1, (*i)[3]+1)) != antiDiagonalLink.end() ){
                printf("swap 3--5\n");
                fluid.get( (*i)[0], (*i)[1] ).f[3] *= (1-zeta);
                fluid.get( (*i)[0], (*i)[1] ).f[3] += zeta*(fluid.get( (*i)[0], (*i)[1] ).f[1] - 1/6*this->getVelocity()[1]);
                fluid.get( (*i)[0], (*i)[1] ).f[5] *= (1-zeta);
                fluid.get( (*i)[0], (*i)[1] ).f[5] += zeta*(fluid.get( (*i)[0], (*i)[1] ).f[7] - 1/6*this->getVelocity()[1]);
                //swap(fluid.get( (*i)[2],(*i)[3] ).f[3], fluid.get( (*i)[2],(*i)[3] ).f[5]);
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
        if(pow((*i)[0] - this->position[0],2.0) + pow((*i)[1] - this->position[1],2.0) >= pow(this->radius,2.0) ){
            if(find(LeftRightLink.begin(), LeftRightLink.end(), Array<int,4>((*i)[0],(*i)[1],(*i)[0]+1, (*i)[1])) != LeftRightLink.end() &&
               find(UpDownLink.begin(),    UpDownLink.end(),    Array<int,4>((*i)[0],(*i)[1]+1,(*i)[0], (*i)[1])) != UpDownLink.end() ){
                printf("swap 2--4\n");
                fluid.get( (*i)[0], (*i)[1] ).f[2] *= (1-zeta);
                fluid.get( (*i)[0], (*i)[1] ).f[2] += zeta*(fluid.get( (*i)[0], (*i)[1] ).f[8] - 1/3*(this->getVelocity()[0]+this->getVelocity()[1]));
                fluid.get( (*i)[0], (*i)[1] ).f[4] *= (1-zeta);
                fluid.get( (*i)[0], (*i)[1] ).f[4] += zeta*(fluid.get( (*i)[0], (*i)[1] ).f[6] - 1/3*(this->getVelocity()[0]+this->getVelocity()[1]));
                //swap(fluid.get( (*i)[0],(*i)[1] ).f[2], fluid.get( (*i)[0],(*i)[1] ).f[4]);
            }
        }
        //bottom-left boundary
        if(pow((*i)[2] - this->position[0],2.0) + pow((*i)[3] - this->position[1],2.0) >= pow(this->radius,2.0) ){
            if(find(LeftRightLink.begin(), LeftRightLink.end(), Array<int,4>((*i)[2]-1,(*i)[3],(*i)[2], (*i)[3])) != LeftRightLink.end() &&
               find(UpDownLink.begin(),    UpDownLink.end(),    Array<int,4>((*i)[2],(*i)[3],(*i)[2], (*i)[3]-1)) != UpDownLink.end() ){
                printf("swap 6--8\n");
                fluid.get( (*i)[0], (*i)[1] ).f[6] *= (1-zeta);
                fluid.get( (*i)[0], (*i)[1] ).f[6] += zeta*(fluid.get( (*i)[0], (*i)[1] ).f[4] + 1/3*(this->getVelocity()[0]+this->getVelocity()[1]));
                fluid.get( (*i)[0], (*i)[1] ).f[8] *= (1-zeta);
                fluid.get( (*i)[0], (*i)[1] ).f[8] += zeta*(fluid.get( (*i)[0], (*i)[1] ).f[2] + 1/3*(this->getVelocity()[0]+this->getVelocity()[1]));
                //swap(fluid.get( (*i)[2],(*i)[3] ).f[6], fluid.get( (*i)[2],(*i)[3] ).f[8]);
            }
        }
    }
//=================================================================================
    printf("antiDiagonalLink adjustment.\n");
    for(std::vector<Array<int,4> >::iterator i=antiDiagonalLink.begin();   i!=antiDiagonalLink.end(); ++i){
        //up-left boundary
        if(pow((*i)[0] - this->position[0],2.0) + pow((*i)[1] - this->position[1],2.0) >= pow(this->radius,2.0) ){
            if(find(LeftRightLink.begin(), LeftRightLink.end(), Array<int,4>((*i)[0]-1,(*i)[1],(*i)[0], (*i)[1])) != LeftRightLink.end() &&
               find(UpDownLink.begin(),    UpDownLink.end(),    Array<int,4>((*i)[0],(*i)[1]+1,(*i)[0], (*i)[1])) != UpDownLink.end() ){
                printf("swap 4--6\n");
                fluid.get( (*i)[0], (*i)[1] ).f[4] *= (1-zeta);
                fluid.get( (*i)[0], (*i)[1] ).f[4] += zeta*(fluid.get( (*i)[0], (*i)[1] ).f[2] + 1/3*(this->getVelocity()[0]-this->getVelocity()[1]));
                fluid.get( (*i)[0], (*i)[1] ).f[6] *= (1-zeta);
                fluid.get( (*i)[0], (*i)[1] ).f[6] += zeta*(fluid.get( (*i)[0], (*i)[1] ).f[8] + 1/3*(this->getVelocity()[0]-this->getVelocity()[1]));
                //swap(fluid.get( (*i)[0],(*i)[1] ).f[4], fluid.get( (*i)[0],(*i)[1] ).f[6]);
            }
        }
        //bottom-right boundary
        if(pow((*i)[2] - this->position[0],2.0) + pow((*i)[3] - this->position[1],2.0) >= pow(this->radius,2.0) ){
            if(find(LeftRightLink.begin(), LeftRightLink.end(), Array<int,4>((*i)[2],(*i)[3],(*i)[2]+1, (*i)[3])) != LeftRightLink.end() &&
               find(UpDownLink.begin(),    UpDownLink.end(),    Array<int,4>((*i)[2],(*i)[3],(*i)[2], (*i)[3]-1)) != UpDownLink.end() ){
                printf("swap 2--8\n");
                fluid.get( (*i)[0], (*i)[1] ).f[2] *= (1-zeta);
                fluid.get( (*i)[0], (*i)[1] ).f[2] += zeta*(fluid.get( (*i)[0], (*i)[1] ).f[4] - 1/3*(this->getVelocity()[0]-this->getVelocity()[1]));
                fluid.get( (*i)[0], (*i)[1] ).f[8] *= (1-zeta);
                fluid.get( (*i)[0], (*i)[1] ).f[8] += zeta*(fluid.get( (*i)[0], (*i)[1] ).f[2] - 1/3*(this->getVelocity()[0]-this->getVelocity()[1]));
                //swap(fluid.get( (*i)[2],(*i)[3] ).f[2], fluid.get( (*i)[2],(*i)[3] ).f[8]);
            }
        }
    }
    return;
}
*/

}
#endif
