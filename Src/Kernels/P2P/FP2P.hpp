// See LICENCE file at project root
#ifndef FP2P_HPP
#define FP2P_HPP

namespace FP2P {

/**
   * @brief MutualParticles (generic version)
   * P2P mutual interaction,
   * this function computes the interaction for 2 particles.
   *
   * Formulas are:
   * \f[
   * F = - q_1 * q_2 * grad K{12}
   * P_1 = q_2 * K{12} ; P_2 = q_1 * K_{12}
   * \f]
   * In details for \f$K(x,y)=1/|x-y|=1/r\f$ :
   * \f$ P_1 = \frac{ q_2 }{ r } \f$
   * \f$ P_2 = \frac{ q_1 }{ r } \f$
   * \f$ F(x) = \frac{ \Delta_x * q_1 * q_2 }{ r^2 } \f$
   *
   * @param sourceX
   * @param sourceY
   * @param sourceZ
   * @param sourcePhysicalValue
   * @param targetX
   * @param targetY
   * @param targetZ
   * @param targetPhysicalValue
   * @param targetPotential
   * @param MatrixKernel pointer to an interaction kernel evaluator
   */

/**
   * @brief NonMutualParticles (generic version)
   * P2P mutual interaction,
   * this function computes the interaction for 2 particles.
   *
   * Formulas are:
   * \f[
   * F = - q_1 * q_2 * grad K{12}
   * P_1 = q_2 * K{12} ; P_2 = q_1 * K_{12}
   * \f]
   * In details for \f$K(x,y)=1/|x-y|=1/r\f$ :
   * \f$ P_1 = \frac{ q_2 }{ r } \f$
   * \f$ P_2 = \frac{ q_1 }{ r } \f$
   * \f$ F(x) = \frac{ \Delta_x * q_1 * q_2 }{ r^2 } \f$
   */

template <class FReal, class ContainerClass, class MatrixKernelClass, class ComputeClass, int NbFRealInComputeClass>
static void GenericFullMutual(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
                              const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){

    const FSize nbParticlesTargets = inTargets->getNbParticles();
    const FReal*const targetsPhysicalValues = inTargets->getPhysicalValues();
    const FReal*const targetsX = inTargets->getPositions()[0];
    const FReal*const targetsY = inTargets->getPositions()[1];
    const FReal*const targetsZ = inTargets->getPositions()[2];
    FReal*const targetsPotentials = inTargets->getPotentials();

    for(FSize idxNeighbors = 0 ; idxNeighbors < limiteNeighbors ; ++idxNeighbors){
        if( inNeighbors[idxNeighbors] ){
            const FSize nbParticlesSources = (inNeighbors[idxNeighbors]->getNbParticles()+NbFRealInComputeClass-1)/NbFRealInComputeClass;
            const ComputeClass*const sourcesPhysicalValues = (const ComputeClass*)inNeighbors[idxNeighbors]->getPhysicalValues();
            const ComputeClass*const sourcesX = (const ComputeClass*)inNeighbors[idxNeighbors]->getPositions()[0];
            const ComputeClass*const sourcesY = (const ComputeClass*)inNeighbors[idxNeighbors]->getPositions()[1];
            const ComputeClass*const sourcesZ = (const ComputeClass*)inNeighbors[idxNeighbors]->getPositions()[2];
            ComputeClass*const sourcesPotentials = (ComputeClass*)inNeighbors[idxNeighbors]->getPotentials();

            for(FSize idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
                const ComputeClass tx = FMath::ConvertTo<ComputeClass, const FReal*>(&targetsX[idxTarget]);
                const ComputeClass ty = FMath::ConvertTo<ComputeClass, const FReal*>(&targetsY[idxTarget]);
                const ComputeClass tz = FMath::ConvertTo<ComputeClass, const FReal*>(&targetsZ[idxTarget]);
                const ComputeClass tv = FMath::ConvertTo<ComputeClass, const FReal*>(&targetsPhysicalValues[idxTarget]);
                ComputeClass  tpo = FMath::Zero<ComputeClass>();

                for(FSize idxSource = 0 ; idxSource < nbParticlesSources ; ++idxSource){
                    ComputeClass Kxy[1];
                                        Kxy[0] = MatrixKernel->evaluate(tx,ty,tz,
                                                             sourcesX[idxSource],sourcesY[idxSource],sourcesZ[idxSource]);
                    const ComputeClass mutual_coeff = FMath::ConvertTo<ComputeClass, FReal>(MatrixKernel->getMutualCoefficient());; // 1 if symmetric; -1 if antisymmetric

                    tpo += Kxy[0] * sourcesPhysicalValues[idxSource];

                    sourcesPotentials[idxSource] += mutual_coeff * Kxy[0] * tv;
                }

                targetsPotentials[idxTarget] += FMath::ConvertTo<FReal, ComputeClass>(tpo);
            }
        }
    }
}

template <class FReal, class ContainerClass, class MatrixKernelClass, class ComputeClass, int NbFRealInComputeClass>
static void GenericInner(ContainerClass* const FRestrict inTargets, const MatrixKernelClass *const MatrixKernel){

    const FSize nbParticlesTargets = inTargets->getNbParticles();
    const FReal*const targetsPhysicalValues = inTargets->getPhysicalValues();
    const FReal*const targetsX = inTargets->getPositions()[0];
    const FReal*const targetsY = inTargets->getPositions()[1];
    const FReal*const targetsZ = inTargets->getPositions()[2];
    FReal*const targetsPotentials = inTargets->getPotentials();

    {//In this part, we compute (vectorially) the interaction
        //within the target leaf.

        const FSize nbParticlesSources = (nbParticlesTargets+NbFRealInComputeClass-1)/NbFRealInComputeClass;
        const ComputeClass*const sourcesPhysicalValues = (const ComputeClass*)targetsPhysicalValues;
        const ComputeClass*const sourcesX = (const ComputeClass*)targetsX;
        const ComputeClass*const sourcesY = (const ComputeClass*)targetsY;
        const ComputeClass*const sourcesZ = (const ComputeClass*)targetsZ;
        ComputeClass*const sourcesPotentials = (ComputeClass*)targetsPotentials;

        for(FSize idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
            const ComputeClass tx = FMath::ConvertTo<ComputeClass, const FReal*>(&targetsX[idxTarget]);
            const ComputeClass ty = FMath::ConvertTo<ComputeClass, const FReal*>(&targetsY[idxTarget]);
            const ComputeClass tz = FMath::ConvertTo<ComputeClass, const FReal*>(&targetsZ[idxTarget]);
            const ComputeClass tv = FMath::ConvertTo<ComputeClass, const FReal*>(&targetsPhysicalValues[idxTarget]);
            ComputeClass  tpo = FMath::Zero<ComputeClass>();

            for(FSize idxSource = (idxTarget+NbFRealInComputeClass)/NbFRealInComputeClass ; idxSource < nbParticlesSources ; ++idxSource){
                ComputeClass Kxy[1];
                                Kxy[0] = MatrixKernel->evaluate(tx,ty,tz,
                                                         sourcesX[idxSource],sourcesY[idxSource],sourcesZ[idxSource]);
                const ComputeClass mutual_coeff = FMath::ConvertTo<ComputeClass, FReal>(MatrixKernel->getMutualCoefficient()); // 1 if symmetric; -1 if antisymmetric

        		tpo = FMath::FMAdd(Kxy[0],sourcesPhysicalValues[idxSource],tpo);

        		sourcesPotentials[idxSource] = FMath::FMAdd(mutual_coeff * Kxy[0],tv,sourcesPotentials[idxSource]);
            }

            targetsPotentials[idxTarget] += FMath::ConvertTo<FReal, ComputeClass>(tpo);
        }
    }

    for(FSize idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
        const FSize limitForTarget = NbFRealInComputeClass-(idxTarget%NbFRealInComputeClass);
        for(FSize idxS = 1 ; idxS < limitForTarget ; ++idxS){
            const FSize idxSource = idxTarget + idxS;
            FReal Kxy[1];
            Kxy[0] = MatrixKernel->evaluate(targetsX[idxTarget],targetsY[idxTarget],targetsZ[idxTarget],
                                                     targetsX[idxSource],targetsY[idxSource],targetsZ[idxSource]);
            const FReal mutual_coeff = MatrixKernel->getMutualCoefficient(); // 1 if symmetric; -1 if antisymmetric

            targetsPotentials[idxTarget] += Kxy[0] * targetsPhysicalValues[idxSource];

            targetsPotentials[idxSource] += mutual_coeff * Kxy[0] * targetsPhysicalValues[idxTarget];
        }
    }
}

template <class FReal, class ContainerClass, class MatrixKernelClass, class ComputeClass, int NbFRealInComputeClass>
static void GenericFullRemote(ContainerClass* const FRestrict inTargets, const ContainerClass* const inNeighbors[],
                              const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){

    const FSize nbParticlesTargets = inTargets->getNbParticles();
    const FReal*const targetsX = inTargets->getPositions()[0];
    const FReal*const targetsY = inTargets->getPositions()[1];
    const FReal*const targetsZ = inTargets->getPositions()[2];
    FReal*const targetsPotentials = inTargets->getPotentials();

    for(FSize idxNeighbors = 0 ; idxNeighbors < limiteNeighbors ; ++idxNeighbors){
        if( inNeighbors[idxNeighbors] ){
            const FSize nbParticlesSources = (inNeighbors[idxNeighbors]->getNbParticles()+NbFRealInComputeClass-1)/NbFRealInComputeClass;
            const ComputeClass*const sourcesPhysicalValues = (const ComputeClass*)inNeighbors[idxNeighbors]->getPhysicalValues();
            const ComputeClass*const sourcesX = (const ComputeClass*)inNeighbors[idxNeighbors]->getPositions()[0];
            const ComputeClass*const sourcesY = (const ComputeClass*)inNeighbors[idxNeighbors]->getPositions()[1];
            const ComputeClass*const sourcesZ = (const ComputeClass*)inNeighbors[idxNeighbors]->getPositions()[2];

            for(FSize idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
                const ComputeClass tx = FMath::ConvertTo<ComputeClass, const FReal*>(&targetsX[idxTarget]);
                const ComputeClass ty = FMath::ConvertTo<ComputeClass, const FReal*>(&targetsY[idxTarget]);
                const ComputeClass tz = FMath::ConvertTo<ComputeClass, const FReal*>(&targetsZ[idxTarget]);
                ComputeClass  tpo = FMath::Zero<ComputeClass>();

                for(FSize idxSource = 0 ; idxSource < nbParticlesSources ; ++idxSource){
                    ComputeClass Kxy[1];
                                        Kxy[0] = MatrixKernel->evaluate(tx,ty,tz,
                                                             sourcesX[idxSource],sourcesY[idxSource],sourcesZ[idxSource]);

                    tpo += Kxy[0] * sourcesPhysicalValues[idxSource];
                }

                targetsPotentials[idxTarget] += FMath::ConvertTo<FReal, ComputeClass>(tpo);
            }
        }
    }
}

} // End namespace

template <class FReal>
struct FP2PT{
};

#if defined(SCALFMM_USE_AVX)
template <>
struct FP2PT<double>{
    template <class ContainerClass, class MatrixKernelClass>
    static void FullMutual(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
                           const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericFullMutual<double, ContainerClass, MatrixKernelClass, __m256d, 4>(inTargets, inNeighbors, limiteNeighbors, MatrixKernel);
    }

    template <class ContainerClass, class MatrixKernelClass>
    static void Inner(ContainerClass* const FRestrict inTargets, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericInner<double, ContainerClass, MatrixKernelClass, __m256d, 4>(inTargets, MatrixKernel);
    }

    template <class ContainerClass, class MatrixKernelClass>
    static void FullRemote(ContainerClass* const FRestrict inTargets, const ContainerClass* const inNeighbors[],
                           const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericFullRemote<double, ContainerClass, MatrixKernelClass, __m256d, 4>(inTargets, inNeighbors, limiteNeighbors, MatrixKernel);
    }
};

template <>
struct FP2PT<float>{
    template <class ContainerClass, class MatrixKernelClass>
    static void FullMutual(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
                           const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericFullMutual<float, ContainerClass, MatrixKernelClass, __m256, 8>(inTargets, inNeighbors, limiteNeighbors, MatrixKernel);
    }

    template <class ContainerClass, class MatrixKernelClass>
    static void Inner(ContainerClass* const FRestrict inTargets, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericInner<float, ContainerClass, MatrixKernelClass, __m256, 8>(inTargets, MatrixKernel);
    }

    template <class ContainerClass, class MatrixKernelClass>
    static void FullRemote(ContainerClass* const FRestrict inTargets, const ContainerClass* const inNeighbors[],
                           const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericFullRemote<float, ContainerClass, MatrixKernelClass, __m256, 8>(inTargets, inNeighbors, limiteNeighbors, MatrixKernel);
    }
};
#elif defined(SCALFMM_USE_AVX2)
template <>
struct FP2PT<double>{
    template <class ContainerClass, class MatrixKernelClass>
    static void FullMutual(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
                           const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericFullMutual<double, ContainerClass, MatrixKernelClass, __m512d, 8>(inTargets, inNeighbors, limiteNeighbors, MatrixKernel);
    }

    template <class ContainerClass, class MatrixKernelClass>
    static void Inner(ContainerClass* const FRestrict inTargets, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericInner<double, ContainerClass, MatrixKernelClass, __m512d, 8>(inTargets, MatrixKernel);
    }

    template <class ContainerClass, class MatrixKernelClass>
    static void FullRemote(ContainerClass* const FRestrict inTargets, const ContainerClass* const inNeighbors[],
                           const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericFullRemote<double, ContainerClass, MatrixKernelClass, __m512d, 8>(inTargets, inNeighbors, limiteNeighbors, MatrixKernel);
    }
};

template <>
struct FP2PT<float>{
    template <class ContainerClass, class MatrixKernelClass>
    static void FullMutual(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
                           const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericFullMutual<float, ContainerClass, MatrixKernelClass, __m512, 16>(inTargets, inNeighbors, limiteNeighbors, MatrixKernel);
    }

    template <class ContainerClass, class MatrixKernelClass>
    static void Inner(ContainerClass* const FRestrict inTargets, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericInner<float, ContainerClass, MatrixKernelClass, __m512, 16>(inTargets, MatrixKernel);
    }

    template <class ContainerClass, class MatrixKernelClass>
    static void FullRemote(ContainerClass* const FRestrict inTargets, const ContainerClass* const inNeighbors[],
                           const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericFullRemote<float, ContainerClass, MatrixKernelClass, __m512, 16>(inTargets, inNeighbors, limiteNeighbors, MatrixKernel);
    }
};
#elif defined(SCALFMM_USE_SSE)
template <>
struct FP2PT<double>{
    template <class ContainerClass, class MatrixKernelClass>
    static void FullMutual(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
                           const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericFullMutual<double, ContainerClass, MatrixKernelClass, __m128d, 2>(inTargets, inNeighbors, limiteNeighbors, MatrixKernel);
    }

    template <class ContainerClass, class MatrixKernelClass>
    static void Inner(ContainerClass* const FRestrict inTargets, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericInner<double, ContainerClass, MatrixKernelClass, __m128d, 2>(inTargets, MatrixKernel);
    }

    template <class ContainerClass, class MatrixKernelClass>
    static void FullRemote(ContainerClass* const FRestrict inTargets, const ContainerClass* const inNeighbors[],
                           const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericFullRemote<double, ContainerClass, MatrixKernelClass, __m128d, 2>(inTargets, inNeighbors, limiteNeighbors, MatrixKernel);
    }
};

template <>
struct FP2PT<float>{
    template <class ContainerClass, class MatrixKernelClass>
    static void FullMutual(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
                           const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericFullMutual<float, ContainerClass, MatrixKernelClass, __m128, 4>(inTargets, inNeighbors, limiteNeighbors, MatrixKernel);
    }

    template <class ContainerClass, class MatrixKernelClass>
    static void Inner(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
                           const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericInner<float, ContainerClass, MatrixKernelClass, __m128, 4>(inTargets, MatrixKernel);
    }

    template <class ContainerClass, class MatrixKernelClass>
    static void FullRemote(ContainerClass* const FRestrict inTargets, const ContainerClass* const inNeighbors[],
                           const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericFullRemote<float, ContainerClass, MatrixKernelClass, __m128, 4>(inTargets, inNeighbors, limiteNeighbors, MatrixKernel);
    }
};
#else
template <>
struct FP2PT<double>{
    template <class ContainerClass, class MatrixKernelClass>
    static void FullMutual(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
                           const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericFullMutual<double, ContainerClass, MatrixKernelClass, double, 1>(inTargets, inNeighbors, limiteNeighbors, MatrixKernel);
    }

    template <class ContainerClass, class MatrixKernelClass>
    static void Inner(ContainerClass* const FRestrict inTargets, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericInner<double, ContainerClass, MatrixKernelClass, double, 1>(inTargets, MatrixKernel);
    }

    template <class ContainerClass, class MatrixKernelClass>
    static void FullRemote(ContainerClass* const FRestrict inTargets, const ContainerClass* const inNeighbors[],
                           const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericFullRemote<double, ContainerClass, MatrixKernelClass, double, 1>(inTargets, inNeighbors, limiteNeighbors, MatrixKernel);
    }
};

template <>
struct FP2PT<float>{
    template <class ContainerClass, class MatrixKernelClass>
    static void FullMutual(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
                           const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericFullMutual<float, ContainerClass, MatrixKernelClass, float, 1>(inTargets, inNeighbors, limiteNeighbors, MatrixKernel);
    }

    template <class ContainerClass, class MatrixKernelClass>
    static void Inner(ContainerClass* const FRestrict inTargets, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericInner<float, ContainerClass, MatrixKernelClass, float, 1>(inTargets, MatrixKernel);
    }

    template <class ContainerClass, class MatrixKernelClass>
    static void FullRemote(ContainerClass* const FRestrict inTargets, const ContainerClass* const inNeighbors[],
                           const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericFullRemote<float, ContainerClass, MatrixKernelClass, float, 1>(inTargets, inNeighbors, limiteNeighbors, MatrixKernel);
    }
};
#endif

#endif // FP2P_HPP
