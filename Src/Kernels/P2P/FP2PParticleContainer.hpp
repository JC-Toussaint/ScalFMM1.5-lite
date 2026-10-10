// See LICENCE file at project root
#ifndef FP2PPARTICLECONTAINER_HPP
#define FP2PPARTICLECONTAINER_HPP

#include "Components/FBasicParticleContainer.hpp"

template<class FReal, int NRHS = 1, int NLHS = 1, int NVALS = 1>
class FP2PParticleContainer : public FBasicParticleContainer<FReal, NVALS*(NRHS+NLHS), FReal> {
    typedef FBasicParticleContainer<FReal, NVALS*(NRHS+NLHS), FReal> Parent;

public:
    static const int NbAttributes = NVALS*(NRHS+NLHS);
    typedef FReal AttributesClass;

    FReal* getPhysicalValues(const int idxVals = 0, const int idxRhs = 0){
      return Parent::getAttribute((0+idxRhs)*NVALS+idxVals);
    }

    const FReal* getPhysicalValues(const int idxVals = 0, const int idxRhs = 0) const {
        return Parent::getAttribute((0+idxRhs)*NVALS+idxVals);
    }




    FReal* getPotentials(const int idxVals = 0, const int idxLhs = 0){
        return Parent::getAttribute((NRHS+idxLhs)*NVALS+idxVals);
    }

    const FReal* getPotentials(const int idxVals = 0, const int idxLhs = 0) const {
        return Parent::getAttribute((NRHS+idxLhs)*NVALS+idxVals);
    }
















    int getNVALS() const {
        return NVALS;
    }

};

#endif // FP2PPARTICLECONTAINER_HPP
