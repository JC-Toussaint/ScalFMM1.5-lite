// See LICENCE file at project root
#ifndef FP2PPARTICLECONTAINERINDEXED_HPP
#define FP2PPARTICLECONTAINERINDEXED_HPP

#include <vector>

#include "FP2PParticleContainer.hpp"
#include "Components/FParticleType.hpp"

template<class FReal, int NRHS = 1, int NLHS = 1, int NVALS = 1>
class FP2PParticleContainerIndexed : public FP2PParticleContainer<FReal, NRHS,NLHS,NVALS> {
    typedef FP2PParticleContainer<FReal, NRHS,NLHS,NVALS> Parent;

    std::vector<FSize> indexes;

public:
    template<typename... Args>
    void push(const FPoint<FReal>& inParticlePosition, const FSize index, Args... args){
        Parent::push(inParticlePosition, args... );
        indexes.push_back(index);
    }

    template<typename... Args>
    void push(const FPoint<FReal>& inParticlePosition, const FParticleType particleType, const FSize index, Args... args){
        Parent::push(inParticlePosition, particleType, args... );
        indexes.push_back(index);
    }

    const std::vector<FSize>& getIndexes() const{
        return indexes;
    }

    void clear(){
        indexes.clear();
        Parent::clear();
    }
};

#endif // FP2PPARTICLECONTAINERINDEXED_HPP
