/** \file
 * Definition of a LB cell -- generic implementation.
 */
#ifndef SOLID_CELL_HH
#define SOLID_CELL_HH

#include <algorithm>
#include "newBulkParticles/solidCell.h"
#include "core/util.h"
#include <cstring>

namespace plb {

////////////////////////// Class Cell /////////////////////////////

/** The possibility to default construct Cell objects facilitates
 * their use in various types of containers. However, they can not
 * be used directly after construction; the method attributeDynamics()
 * must be used first.
 */
template<typename T, template<typename U> class Descriptor>
solidCell<T,Descriptor>::solidCell()
    : Cell<T,Descriptor>()
{
    iniSolidFrac();
}

/** This constructor initializes the dynamics, but not the values
 * of the distribution functions. Remember that the dynamics is not
 * owned by the Cell object, the user must ensure its proper
 * destruction and a sufficient life time.
 */
template<typename T, template<typename U> class Descriptor>
solidCell<T,Descriptor>::solidCell(Dynamics<T,Descriptor>* dynamics_)
    : Cell<T,Descriptor>(dynamics_)
{
    iniSolidFrac();
}


template<typename T, template<typename U> class Descriptor>
void solidCell<T,Descriptor>::iniSolidFrac() {
    solidFrac = 0;
}

template<typename T, template<typename U> class Descriptor>
void solidCell<T,Descriptor>::solidAttributeDynamics(Dynamics<T,Descriptor>* dynamics_) {
    //solidDynamics = dynamics_;
    this->dynamics = dynamics_;
}

/*
template<typename T, template<typename U> class Descriptor>
void solidCell<T,Descriptor>::serialize(char* data) const {
    const plint numPop = Descriptor<T>::numPop;
    const plint numExt = Descriptor<T>::ExternalField::numScalars;
    memcpy((void*)data, (const void*)(&f[0]), numPop*sizeof(T));
    if (numExt>0) {
        memcpy((void*)(data+numPop*sizeof(T)), (const void*)(external.get(0)), numExt*sizeof(T));
    }
    //serialize for solidFrac and solidVelocity
    memcpy((void*)(data+numPop*sizeof(T)+numExt*sizeof(T)), (const void*)(&solidFrac), sizeof(T));
    memcpy((void*)(data+numPop*sizeof(T)+numExt*sizeof(T)+sizeof(T)), (const void*)(&solidFrac), 2*sizeof(T));
}

template<typename T, template<typename U> class Descriptor>
void solidCell<T,Descriptor>::unSerialize(char const* data) {
    const plint numPop = Descriptor<T>::numPop;
    const plint numExt = Descriptor<T>::ExternalField::numScalars;

    memcpy((void*)(&f[0]), (const void*)data, numPop*sizeof(T));
    if (numExt>0) {
        memcpy((void*)(external.get(0)), (const void*)(data+numPop*sizeof(T)), numExt*sizeof(T));
    }
    //unSerialize for solidFrac...
    memcpy((void*)(&solidFrac), (const void*)(data+numPop*sizeof(T)+numExt*sizeof(T)), sizeof(T));
    memcpy((void*)(&solidVelocity), (const void*)(data+numPop*sizeof(T)+numExt*sizeof(T)+2*sizeof(T), 2*sizeof(T));
}
*/

template<typename T, template<typename U> class Descriptor>
void iniCellAtEquilibrium(solidCell<T,Descriptor>& cell, T density, Array<T,Descriptor<T>::d> const& velocity) {
    //std::cout<<"iniCellAtEquilibrium for solidCell."<<std::endl;
    Array<T,Descriptor<T>::d> j;
    VectorTemplate<T,Descriptor>::multiplyByScalar(velocity, density, j);
    T jSqr = VectorTemplate<T,Descriptor>::normSqr(j);
    T rhoBar = Descriptor<T>::rhoBar(density);
    for (plint iPop=0; iPop<Descriptor<T>::numPop; ++iPop) {
        cell[iPop] = cell.computeEquilibrium(iPop, rhoBar, j, jSqr);
    }
}

template<typename T, template<typename U> class Descriptor>
void iniCellAtEquilibrium(solidCell<T,Descriptor>& cell, T density, Array<T,Descriptor<T>::d> const& velocity, T temperature) {
    Array<T,Descriptor<T>::d> j;
    VectorTemplate<T,Descriptor>::multiplyByScalar(velocity, density, j);
    T jSqr = VectorTemplate<T,Descriptor>::normSqr(j);
    T rhoBar = Descriptor<T>::rhoBar(density);
    for (plint iPop=0; iPop<Descriptor<T>::numPop; ++iPop) {
        cell[iPop] = cell.computeEquilibrium(iPop, rhoBar, j, jSqr, temperature-(T)1);
    }
}

}  // namespace plb

#endif  // SOLID_CELL_HH
