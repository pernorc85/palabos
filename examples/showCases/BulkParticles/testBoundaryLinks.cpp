#include "palabos2D.h"
#include "palabos2D.hh"
//#include <cstdlib>
//#include <iostream>

using namespace plb;
using namespace std;

// Use double-precision arithmetics
typedef double T;

// Use a grid which additionally to the f's stores two variables for
//   the external force term.
#define DESCRIPTOR descriptors::D2Q9Descriptor

int main(int argc, char *argv[])
{
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");
    //srand(global::mpi().getRank());
    const double PI = 3.1415926535897932;
    plint nx = 1200;
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
    boundaryCondition->setVelocityConditionOnBlockBoundaries(myFluid);
    setBoundaryVelocity(myFluid, myFluid.getBoundingBox(), Array<T,2>(0.,0.) );
    initializeAtEquilibrium(myFluid, myFluid.getBoundingBox(), 1., Array<T,2>(0.,0.), 0. );
    //myFluid2.periodicity().toggle(0,true);
    std::cout<<"ok0"<<std::endl;
    //setBoundaryVelocity(myFluid, Box2D(1, nx, ny-1, ny-1), Array<T,2>(10.,0.) );
    //initializeAtEquilibrium(myFluid, Box2D(1, nx, ny-1, ny-1), 1., Array<T,2>(10.,0.), 0. );
    //std::cout<<"ok1"<<std::endl;
    myFluid.initialize();

    //MultiParticleField2D<BulkParticleField2D<T,DESCRIPTOR> >* particles    =0;

    BulkParticleField2D<T,DESCRIPTOR> myParticleField(nx, ny);

    BulkParticle2D<T,DESCRIPTOR> myParticle1(1, 1, Array<T,2>(300-100,300), 0,
                                                   Array<T,2>(0.0, 0), 0.0, 20, 5);
    BulkParticle2D<T,DESCRIPTOR> myParticle2(2, 1, Array<T,2>(300+100*cos(2*PI/3),300+100*sin(2*PI/3)), 0,
                                                   Array<T,2>(-0.0*cos(2*PI/3), -0.0*sin(2*PI/3)), 0, 20, 5);
    BulkParticle2D<T,DESCRIPTOR> myParticle3(3, 1, Array<T,2>(300+100*cos(PI/3),  300+100*sin(PI/3)), 0,
                                                   Array<T,2>(-0.0*cos(PI/3),   -0.0*sin(PI/3)),   0, 20, 5);
    BulkParticle2D<T,DESCRIPTOR> myParticle4(4, 1, Array<T,2>(300+100,300), 0,
                                                   Array<T,2>(-0.0, 0), 0.0, 30, 5);
    BulkParticle2D<T,DESCRIPTOR> myParticle5(5, 1, Array<T,2>(300+100*cos(5*PI/3),300+100*sin(5*PI/3)), 0,
                                                   Array<T,2>(-0.0*cos(5*PI/3), -0.0*sin(5*PI/3)),0, 20, 5);
    BulkParticle2D<T,DESCRIPTOR> myParticle6(6, 1, Array<T,2>(300+100*cos(4*PI/3),300+100*sin(4*PI/3)), 0,
                                                Array<T,2>(-0.0*cos(4*PI/3), -0.0*sin(4*PI/3)),0, 20, 5);


    movingLumen2D<T,DESCRIPTOR> myLumen1(7, 450, 750, 200, 250, Array<T,2> (0.,-0.05));
    movingLumen2D<T,DESCRIPTOR> myLumen2(8, 450, 750, 350, 400, Array<T,2> (0., 0.05));


    //BulkParticle2D<T,DESCRIPTOR> myParticle7(7, 2, Array<T,2>(300,50), 0,
    //                                               Array<T,2>(0,0.1),0, 30, 5);

    //BulkParticle2D(plint id_, plint tag_, Array<T,2> const& position_, T AngularPosition_,
    //                                      Array<T,2> const& velocity_, T AngularVelocity_, T radius_, T density_);


    myParticleField.addParticle(Box2D(0,nx-1, 0,ny-1), &myParticle1);
    myParticleField.addParticle(Box2D(0,nx-1, 0,ny-1), &myParticle2);
    myParticleField.addParticle(Box2D(0,nx-1, 0,ny-1), &myParticle3);
    myParticleField.addParticle(Box2D(0,nx-1, 0,ny-1), &myParticle4);
    myParticleField.addParticle(Box2D(0,nx-1, 0,ny-1), &myParticle5);
    myParticleField.addParticle(Box2D(0,nx-1, 0,ny-1), &myParticle6);
    //myParticleField.addParticle(Box2D(0,nx-1, 0,ny-1), &myParticle7);


    for (plint iT=0; iT<maxIter; ++iT) {
        std::cout<<"=====================Identify particles boundary links==================\n\n"<<std::endl;

        for(plint i=0;i<myParticleField.particles.size(); i++){
            Box2D myBox(myParticleField.particles[i]->getPosition()[0] - myParticleField.particles[i]->getRadius()*1.414,
                        myParticleField.particles[i]->getPosition()[0] + myParticleField.particles[i]->getRadius()*1.414,
                        myParticleField.particles[i]->getPosition()[1] - myParticleField.particles[i]->getRadius()*1.414,
                        myParticleField.particles[i]->getPosition()[1] + myParticleField.particles[i]->getRadius()*1.414);
            myParticleField.particles[i]->computeBoundaryLinks(myBox);
            //myParticleChunk.particles[i]->computeBoundaries(myBox);
            Array<T,2> velocity = myParticleField.particles[i]->getVelocity();
            Array<T,2> position = myParticleField.particles[i]->getPosition();
            std::cout<<"Position="<<position[0]<<","<<position[1]<<std::endl;
            std::cout<<"Velocity="<<velocity[0]<<","<<velocity[1]<<std::endl;
            //system("PAUSE");
        }
        myLumen1.computeBoundaryLinks(myFluid.getBoundingBox());
        myLumen2.computeBoundaryLinks(myFluid.getBoundingBox());

        std::cout<<"=============================Collide and stream==========================\n\n"<<std::endl;
        myFluid.collideAndStream();

        std::cout<<"=============particle, movingObject, wall and fluid interaction=============\n\n"<<std::endl;
        for(plint i=0;i<myParticleField.particles.size(); i++){
            myParticleField.particles[i]->particleToParticle(myParticleField);
        }

        for(plint i=0;i<myParticleField.particles.size(); i++){
            myParticleField.particles[i]->wallToParticle(myFluid);
        }

        for(plint i=0;i<myParticleField.particles.size(); i++){
            myLumen1.objectToParticle(myParticleField.particles[i]);
            myLumen2.objectToParticle(myParticleField.particles[i]);
        }

        myParticleField.particleToFluidCoupling(myFluid);
        myLumen1.objectToFluid(myFluid);
        myLumen2.objectToFluid(myFluid);
        std::cout<<"=============particle and movingObject advance=============================\n\n"<<std::endl;
        if(iT<=2000){
            myLumen1.advance();
            myLumen2.advance();
        }

        if(iT==5000){
            myLumen1.setVelocity(Array<T,2> (0., 0.07));
            myLumen2.setVelocity(Array<T,2> (0.,-0.07));
        }
        if(iT>=5000 && iT<=7000){
            myLumen1.advance();
            myLumen2.advance();
        }

        myParticleField.advanceParticles(myFluid.getBoundingBox(), -1.);
        ImageWriter<T> imageWriter("leeloo.map");
        if(iT%10 == 0){
            //imageWriter.writeScaledGif(createFileName("fluid_", iT, 6),
            //                   *computeDensity(myFluid));
            imageWriter.writeGif(createFileName("fluid_", iT, 6),
                           //*computeVelocityNorm(myFluid), 0, 1);
                             *computeVelocityComponent(myFluid,0), -0.1, 0.1);
        }
    }
}
