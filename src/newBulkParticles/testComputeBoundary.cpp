#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#define TRUE 1
#define FALSE 0

using namespace std;

class Box2D{
public:
    Box2D(int x0_,int x1_, int y0_, int y1_):x0(x0_),x1(x1_),y0(y0_),y1(y1_){}

    int x0, x1;
    int y0, y1;
};

template<typename T>
class Array {
public:
    Array() { }
    Array(T x, T y) {
        data[0] = x;
        data[1] = y;
    }
    template<typename U>
    Array(Array<U> const& rhs) {
        data[0] = (T) rhs[0];
        data[1] = (T) rhs[1];
    }
    template<typename U>
    Array<T>& operator=(Array<U> const& rhs) {
        data[0] = (T) rhs[0];
        data[1] = (T) rhs[1];
        return *this;
    }
    T& operator[](int index) {
        return data[index];
    }
    
    bool operator==(Array<T> const& rhs) {
         if(data[0] == rhs.data[0] &&
            data[1] == rhs.data[1] )return 1;
         else return 0;
    }
private:
    T data[2];
};

template<typename T>
class BulkParticle2D{
public:
    BulkParticle2D(Array<T> const& position_, double AngularPosition_, double radius_)
    : position(position_),
      AngularPosition(AngularPosition_),
      radius(radius_)
    { }
    void computeBoundaries(Box2D domain);
    bool inside(T x, T y);
private:
    Array<T> position;
    double AngularPosition;
    double radius;
    vector<Array<int> > fluidBoundary;
    vector<Array<int> > solidBoundary;
    vector<Array<int> > internalSolid;
};

template<typename T>
bool BulkParticle2D<T>::inside(T x, T y){
    return ( (x-position[0])*(x-position[0]) + (y-position[1])*(y-position[1]) <= radius*radius );
}

template<typename T>    
void BulkParticle2D<T>::computeBoundaries(Box2D domain) {
        for (int iX=domain.x0; iX<=domain.x1; ++iX) {
            for (int iY=domain.y0; iY<=domain.y1; ++iY) {
                if( inside(iX, iY) && 
                                 (!inside(iX-1,iY) ||
                                  !inside(iX+1,iY) ||
                                  !inside(iX,iY-1) ||
                                  !inside(iX,iY+1) ||
                                  !inside(iX-1,iY-1) ||
                                  !inside(iX-1,iY+1) ||
                                  !inside(iX+1,iY-1) ||
                                  !inside(iX+1,iY+1)) ){
                    solidBoundary.push_back(Array<int>(iX, iY));
                }
                else if( inside(iX, iY) ) {
                    internalSolid.push_back(Array<int>(iX, iY));
                }
            }
        }
        cout<<"solidBoundary:"<<endl;
        for(vector<Array<int> >::iterator i=solidBoundary.begin(); i!=solidBoundary.end(); ++i)
            cout<<(*i)[0]<<" "<<(*i)[1]<<endl;
            
        for (int iX=domain.x0; iX<=domain.x1; ++iX) {
            for (int iY=domain.y0; iY<=domain.y1; ++iY) {
                if( !inside(iX, iY) &&
                                (find(solidBoundary.begin(),solidBoundary.end(),Array<int>(iX-1,iY))!=solidBoundary.end() ||
                                 find(solidBoundary.begin(),solidBoundary.end(),Array<int>(iX+1,iY))!=solidBoundary.end() ||
                                 find(solidBoundary.begin(),solidBoundary.end(),Array<int>(iX,iY-1))!=solidBoundary.end() ||
                                 find(solidBoundary.begin(),solidBoundary.end(),Array<int>(iX,iY+1))!=solidBoundary.end() ||
                                 find(solidBoundary.begin(),solidBoundary.end(),Array<int>(iX-1,iY-1))!=solidBoundary.end() ||
                                 find(solidBoundary.begin(),solidBoundary.end(),Array<int>(iX+1,iY-1))!=solidBoundary.end() ||
                                 find(solidBoundary.begin(),solidBoundary.end(),Array<int>(iX-1,iY+1))!=solidBoundary.end() ||
                                 find(solidBoundary.begin(),solidBoundary.end(),Array<int>(iX+1,iY+1))!=solidBoundary.end() ) ){                                
                    fluidBoundary.push_back(Array<int>(iX, iY));
                }
            }
        }
        cout<<"fluidBoundary:"<<endl;
        for(vector<Array<int> >::iterator i=fluidBoundary.begin(); i!=fluidBoundary.end(); ++i)
            cout<<(*i)[0]<<" "<<(*i)[1]<<endl;
}

int main()
{
    Box2D myDomain(0,10,0,10);
    BulkParticle2D<int> myBulkParticle(Array<int>(5,5),0, 2);
    myBulkParticle.computeBoundaries(myDomain);
    system("PAUSE");
    return 0;
}
    
