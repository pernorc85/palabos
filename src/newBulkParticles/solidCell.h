#ifndef SOLID_CELL_H
#define SOLID_CELL_H

#include "core/globalDefs.h"
#include "core/plbDebug.h"
#include "latticeBoltzmann/externalFields.h"
#include "core/dynamics.h"
#include "core/array.h"
#include "core/hierarchicSerializer.h"
#include "core/cell.h"

namespace plb {


template<typename T, template<typename U> class Descriptor>
class solidCell : public Cell<T,Descriptor> {
//public:
    /// Additional per-cell scalars for external fields, e.g. forces
//    typedef ExternalFieldArray<T, typename Descriptor<T>::ExternalField> External;
public:
    /// Default constructor.
    solidCell();
    /// Constructor, to be used whenever possible.
    solidCell(Dynamics<T,Descriptor>* dynamics_);
public:
//    T& operator[](plint iPop) {}
//    T const& operator[](plint iPop) const { }
//    Array<T,Descriptor<T>::numPop>& getRawPopulations() {}
//    Array<T,Descriptor<T>::numPop> const& getRawPopulations() const {}
      T getSolidFrac() { return solidFrac; }
      T getSolidDensity() { return solidDensity; }
      Array<T,2> getSolidVelocity() { return solidVelocity; }
      plint getBoundaryType() { return BoundaryType; }
//    Cell<T,Descriptor>& attributeF(Cell<T,Descriptor> const& rhs) {}
//    Cell<T,Descriptor>& attributeValues(Cell<T,Descriptor> const& rhs) {}
//    T* getExternal(plint offset) {}
//    T const* getExternal(plint offset) const {}
//    Dynamics<T,Descriptor> const& getDynamics() const;
//    Dynamics<T,Descriptor>& getDynamics();
//    bool takesStatistics() const {}
//    void specifyStatisticsStatus(bool status) {}
private:
    /// You can't use this method. Use BlockLattice::attributeDynamics instead.
    /** This is one of the rare cases of a method accepting a pointer but
     *  not managing memory itself. Memory of the dynamics is handled by
     *  the BlockLattice.
     */
    void solidAttributeDynamics(Dynamics<T,Descriptor>* dynamics_);
    // Declare the BlockLatticeXD as a friend, to enable access to attributeDynamics.
    template<typename T_, template<typename U_> class Descriptor_> friend class solidBlockLattice2D;
    template<typename T_, template<typename U_> class Descriptor_> friend class solidBlockLattice3D;
#ifdef PLB_MPI_PARALLEL
    template<typename T_, template<typename U_> class Descriptor_> friend class ParallelCellAccess2D;
    template<typename T_, template<typename U_> class Descriptor_> friend class ParallelCellAccess3D;
#endif

// The following helper functions forward the function call
// to the Dynamics object
public:
    void collide(BlockStatistics& statistics) {
        //std::cout<<"solidCell collide"<<std::endl;
        PLB_PRECONDITION( dynamics );
        this->dynamics->collide(*this, statistics, 1);
    }
//    T computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
//                         T jSqr, T thetaBar=T()) const{}
//    void regularize(T rhoBar, Array<T,Descriptor<T>::d> const& j, T jSqr,
//                    Array<T,SymmetricTensor<T,Descriptor>::n> const& PiNeq, T thetaBar=T() ){}
//    T computeDensity() const;
//    T computeDensity() const {
//        PLB_PRECONDITION( solidDynamics );
//        //std::cout<<"computeDensity in Cell."<<std::endl;
//        //dynamics->report();
//        return solidDynamics->computeDensity(*this);
//    }
//    T computePressure() const;
//    void computeVelocity(Array<T,Descriptor<T>::d>& u) const {}
//    T computeTemperature() const {}
//    void computePiNeq (Array<T,SymmetricTensor<T,Descriptor>::n>& PiNeq ) const {}
//    void computeShearStress (Array<T,SymmetricTensor<T,Descriptor>::n>& stress ) const {}
//    void computeHeatFlux(Array<T,Descriptor<T>::d>& q) const {}
//    void computeMoment(plint momentId, T* moment) {}
//    void getPopulations(Array<T,Descriptor<T>::numPop>& f) const{}
//    void getExternalField(plint pos, plint size, T* ext) const {}
      void defineSolidFrac(T solidFrac_) {
           solidFrac = solidFrac_;
      }

      void defineSolidDensity(T solidDensity_) {
           solidDensity = solidDensity_;
      }

      void defineSolidVelocity(Array<T,2> solidVelocity_) {
           solidVelocity = solidVelocity_;
      }

      void defineBoundaryType(plint BoundaryType_) {
           BoundaryType = BoundaryType_;
      }
//    void defineDensity(T rho) {}
//    void defineVelocity(Array<T,Descriptor<T>::d> const& u) {}
//    void defineTemperature(T temperature) {}
//    void defineHeatFlux(Array<T,Descriptor<T>::d> const& q) {}
//    void definePiNeq (Array<T,SymmetricTensor<T,Descriptor>::n> const& PiNeq ) {}
//    void defineMoment(plint momentId, T const* value) {}
//    void setPopulations(Array<T,Descriptor<T>::numPop> const& f) {}
//    void setExternalField(plint pos, plint size, const T* ext) {}
    /// Revert ("bounce-back") the distribution functions.
//    void revert();
//    void serialize(char* data) const;
//    void unSerialize(char const* data);
private:
    void iniSolidFrac();
protected:
//    FP_BGKdynamics<T,Descriptor>*        solidDynamics;  ///< local LB dynamics
    plint      BoundaryType;
    // 0---internalSolid;
    // 1---solidBoundary;
    // 2---fluidBoundary;
    T          solidFrac;
    Array<T,2> solidVelocity;
    T          solidDensity;
};


}  // namespace plb

#endif  // SOLID_CELL_H
