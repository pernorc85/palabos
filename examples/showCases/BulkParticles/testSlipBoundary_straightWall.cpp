#include "palabos2D.h"
#include "palabos2D.hh"
//#include <cstdlib>
//#include <iostream>

using namespace plb;
using namespace std;

// Use double-precision arithmetics
typedef double T;
const double PI = 3.1415926535897932;
// Use a grid which additionally to the f's stores two variables for
//   the external force term.
#define DESCRIPTOR descriptors::D2Q9Descriptor
FILE *dataFile;

int main(int argc, char *argv[])
{
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");
    //srand(global::mpi().getRank());
    plint nx = 640;
    plint ny = 600;
    const int maxIter  = 1600000;
    const int saveIter = 100;
    const int statIter = 10;
    double omega = 0.8;
    double rho1 = 1;
    double rho0 = 0;

    BlockLattice2D<T,DESCRIPTOR> myFluid (nx, ny, new BGKdynamics<T,DESCRIPTOR>(omega));

    OnLatticeBoundaryCondition2D<T,DESCRIPTOR>*
        boundaryCondition = createLocalBoundaryCondition2D<T,DESCRIPTOR>();
        boundaryCondition = createInterpBoundaryCondition2D<T,DESCRIPTOR>();
    Box2D inlet(2,nx-2, 0, 0);
    Box2D outlet(2,nx-2, ny-1, ny-1);
    boundaryCondition->setVelocityConditionOnBlockBoundaries(myFluid, inlet, boundary::outflow);
    boundaryCondition->setVelocityConditionOnBlockBoundaries(myFluid, outlet, boundary::outflow);
//    boundaryCondition.setPressureConditionOnBlockBoundaries (
//                lattice, Box2D(1,nx-2, 0, 0) );
//    boundaryCondition.setPressureConditionOnBlockBoundaries (
//                lattice, Box2D(1,nx-2, ny-1, ny-1) );

    initializeAtEquilibrium(myFluid, myFluid.getBoundingBox(), 1., Array<T,2>(0.,0.), 0. );
    //myFluid.periodicity().toggle(1,true);
    std::cout<<"ok0"<<std::endl;
    //setBoundaryVelocity(myFluid, Box2D(1, nx, ny-1, ny-1), Array<T,2>(10.,0.) );
    //initializeAtEquilibrium(myFluid, Box2D(1, nx, ny-1, ny-1), 1., Array<T,2>(10.,0.), 0. );
    //std::cout<<"ok1"<<std::endl;
    myFluid.initialize();


    T cos_theta = 600/sqrt(600*600+100*100);
    T sin_theta = 100/sqrt(600*600+100*100);
    movingInclinedWall2D<T,DESCRIPTOR> myWall1(1, 190, 0, 600, 600/100, Array<T,2> (-0.05*sin_theta,-0.05*cos_theta));
    movingInclinedWall2D<T,DESCRIPTOR> myWall2(2, 350, 0, 600, 600/100, Array<T,2> (-0.05*sin_theta,-0.05*cos_theta));

    vector<T> velocity_parallel;

    for (plint iT=0; iT<maxIter; ++iT) {
        std::cout<<"=====================Identify wall boundary links==================\n\n"<<std::endl;

        Box2D myBox(myWall1.get_x0() - 10., myWall1.get_x0() + (myWall1.get_y1()-myWall1.get_y0())/myWall1.get_k() + 10.,
                    myWall1.get_y0(),       myWall1.get_y1());
        myWall1.computeBoundaryLinks(myBox);
        Box2D myBox2(myWall2.get_x0()- 10., myWall2.get_x0() + (myWall2.get_y1()-myWall2.get_y0())/myWall2.get_k() + 10,
                     myWall2.get_y0(),      myWall2.get_y1());
        myWall2.computeBoundaryLinks(myBox2);


        std::cout<<"=============================Collide and stream==========================\n\n"<<std::endl;
        myFluid.collideAndStream();

        std::cout<<"=============movingObject and fluid interaction==========================\n\n"<<std::endl;
        myWall1.objectToFluid(myFluid);
        myWall2.objectToFluid(myFluid);
        std::cout<<"=============particle and movingObject advance=============================\n\n"<<std::endl;
        if(iT<=2000){
            myWall1.advance();
            myWall2.advance();
        }

        ImageWriter<T> imageWriter("leeloo.map");
        if(iT%10 == 0){
            //imageWriter.writeScaledGif(createFileName("fluid_", iT, 6),
            //                   *computeDensity(myFluid));
            imageWriter.writeGif(createFileName("fluid_", iT, 6),
                           //*computeVelocityNorm(myFluid), 0, 1);
                             *computeVelocityComponent(myFluid,1), -0.1, 0.1);
        }

        if(iT%100 == 0){
            velocity_parallel.clear();
            T v_parallel;
            for(plint x=0;x<nx;x++){
                plint y=ny/2;
                T v_y =(myFluid.get(x,y)[1] + myFluid.get(x,y)[8] + myFluid.get(x,y)[7])
                        -(myFluid.get(x,y)[3] + myFluid.get(x,y)[4] + myFluid.get(x,y)[5]);
                T v_x =(myFluid.get(x,y)[5] + myFluid.get(x,y)[6] + myFluid.get(x,y)[7])
                        -(myFluid.get(x,y)[1] + myFluid.get(x,y)[2] + myFluid.get(x,y)[3]);
                v_parallel = v_y*cos_theta + v_x*sin_theta;
                velocity_parallel.push_back(v_parallel);
            }
            dataFile = fopen("velocity_parallel.txt","a");
            for(plint i=0;i<velocity_parallel.size();i++){
                fprintf(dataFile, "%f\t", velocity_parallel[i]);
            }
            fprintf(dataFile, "\n");
            fclose(dataFile);
        }
    }
}
