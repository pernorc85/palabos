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
    plint nx = 1200;
    plint ny = 400;
    const int maxIter  = 16000;
    const int saveIter = 100;
    const int statIter = 10;
    double omega = 0.7;
    double rho1 = 1;
    double rho0 = 0;

    FP_BGKdynamics<T,DESCRIPTOR> myDynamics(omega);
    solidBlockLattice2D<T,DESCRIPTOR> myFluid (nx, ny, new FP_BGKdynamics<T,DESCRIPTOR>(omega));
    //BlockLattice2D<T,DESCRIPTOR> myFluid2 (nx, ny, new BGKdynamics<T,DESCRIPTOR>(omega));

    //defineDynamics(myFluid, Box2D(0,nx-1, 0,0), new BounceBack<T, DESCRIPTOR>(rho0) );
    //defineDynamics(myFluid, Box2D(0,nx-1, ny-1,ny-1), new BounceBack<T, DESCRIPTOR>(rho0) );
    //defineDynamics(myFluid, Box2D(0,0, 0,ny-1), new BounceBack<T, DESCRIPTOR>(rho0) );
    //defineDynamics(myFluid, Box2D(nx-1,nx-1, 0,ny-1), new BounceBack<T, DESCRIPTOR>(rho0) );


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
/*
    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice);

    setBoundaryVelocity(lattice, lattice.getBoundingBox(), Array<T,2>(0.,0.) );
    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), 1., Array<T,2>(0.,0.) );

    T u = parameters.getLatticeU();
    setBoundaryVelocity(lattice, Box2D(1, nx, ny-1, ny-1), Array<T,2>(u,0.) );
    initializeAtEquilibrium(lattice, Box2D(1, nx, ny-1, ny-1), 1., Array<T,2>(u,0.) );

    lattice.initialize();
*/
    //MultiParticleField2D<BulkParticleField2D<T,DESCRIPTOR> >* particles    =0;

    BulkParticleField2D<T,DESCRIPTOR> myParticleField(nx, ny);

    BulkParticle2D<T,DESCRIPTOR> myParticle1(1, Array<T,2>(200,200), 0,
                                                Array<T,2>(0,0),  10, 50, 1);
    //BulkParticle2D<T,DESCRIPTOR> myParticle2(2, Array<T,2>(200,200), 0,
    //                                            Array<T,2>(0,0),  -10, 50, 1);
    //BulkParticle2D<T,DESCRIPTOR> myParticle3(3, Array<T,2>(0,1), 0,
    //                                            Array<T,2>(50,50), 2, 30, 5);
    //BulkParticle2D<T,DESCRIPTOR> myParticle4(4, Array<T,2>(500,1), 0,
    //                                            Array<T,2>(-50,50),2, 30, 5);

    //BulkParticle2D(plint tag_, Array<T,2> const& position_, T AngularPosition_,
    //                           Array<T,2> const& velocity_, T AngularVelocity_, T radius_, T density_);


    myParticleField.addParticle(Box2D(0,nx-1, 0,ny-1), &myParticle1);
    //myParticleField.addParticle(Box2D(0,nx-1, 0,ny-1), &myParticle2);
    //myParticleField.addParticle(Box2D(0,nx-1, 0,ny-1), &myParticle3);
    //myParticleField.addParticle(Box2D(0,nx-1, 0,ny-1), &myParticle4);


    for (plint iT=0; iT<maxIter; ++iT) {
        std::cout<<"========Update particles boundaries and attribute to fluid field========\n\n"<<std::endl;

        for(plint i=0;i<myParticleField.particles.size(); i++){
            Box2D myBox(myParticleField.particles[i]->getPosition()[0] - myParticleField.particles[i]->getRadius()*1.414,
                        myParticleField.particles[i]->getPosition()[0] + myParticleField.particles[i]->getRadius()*1.414,
                        myParticleField.particles[i]->getPosition()[1] - myParticleField.particles[i]->getRadius()*1.414,
                        myParticleField.particles[i]->getPosition()[1] + myParticleField.particles[i]->getRadius()*1.414);
            myParticleField.particles[i]->computeBoundaries(myBox);
            Array<T,2> velocity = myParticleField.particles[i]->getVelocity();
            Array<T,2> position = myParticleField.particles[i]->getPosition();
            std::cout<<"Position="<<position[0]<<","<<position[1]<<std::endl;
            std::cout<<"Velocity="<<velocity[0]<<","<<velocity[1]<<std::endl;
            system("PAUSE");
        }
        myFluid.attributeSolidFracAndDensity(myParticleField);
        myFluid.attributeSolidVelocity(myParticleField);
/*
        for(pluint iX = 0; iX < myFluid.getNx(); iX++){
            for(pluint iY = 0; iY < myFluid.getNy(); iY++){
                iniCellAtEquilibrium(myFluid.getSolidCell(iX, iY), 1., myFluid.getSolidCell(iX,iY).getSolidVelocity() );
            }
        }
*/
        std::cout<<"=============================Collide and stream==========================\n\n"<<std::endl;
        myFluid.collideAndStream();
        ImageWriter<T> imageWriter("leeloo.map");
        if(iT%1 == 0){
            //imageWriter.writeScaledGif(createFileName("fluid_", iT, 6),
            //                   *computeDensity(myFluid));
            imageWriter.writeGif(createFileName("fluid_", iT, 6),
                            *computeVelocityNorm(myFluid), 0, 10);
        }

        std::cout<<"=============particle to wall, fluid and particle interaction=============\n\n"<<std::endl;
        for(plint i=0;i<myParticleField.particles.size(); i++){
            myParticleField.particles[i]->particleToParticle(myParticleField);
        }
        for(plint i=0;i<myParticleField.particles.size(); i++){
            myParticleField.particles[i]->wallToParticle(myFluid);
        }

        myParticleField.fluidToParticleCoupling(myFluid.getBoundingBox(),myFluid);
        myParticleField.advanceParticles(myFluid.getBoundingBox(), -1.);

    }
}

/*
        for(plint i=0;i<myFluid.getNx();i++){
            for(plint j=0;j<myFluid.getNy();j++){
                Array<T,DESCRIPTOR<T>::q>& f = myFluid.getSolidCell(i,j).getRawPopulations();
                if(f[0]!=0){
                    std::cout<<"after collideAndStream"<<f[0]<<","<<f[1]<<","<<f[2]<<","<<f[3]<<","<<f[4]<<","<<f[5]<<","<<f[6]<<","<<f[7]<<","<<f[8]<<std::endl;
                    std::cout<<"density is"<<myFluid.getSolidCell(i,j).computeDensity()<<std::endl;
                }
            }
            std::cout<<endl;
        }
*/

/*
            ScalarField2D<T> myField = *computeDensity(myFluid);
            for(plint i=0;i<myField.getNx();i++){
                for(plint j=0;j<myField.getNy();j++){
                    if(myField.get(i,j)!=1)
                        std::cout<<myField.get(i,j)<<" ";
                }
                std::cout<<endl;
            }
*/
