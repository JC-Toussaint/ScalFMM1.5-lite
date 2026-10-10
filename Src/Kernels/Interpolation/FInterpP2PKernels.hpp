#ifndef FINTERPP2PKERNELS_HPP
#define FINTERPP2PKERNELS_HPP


#include "../P2P/FP2P.hpp"
#include "../P2P/FP2PR.hpp"

///////////////////////////////////////////////////////
// P2P Wrappers
///////////////////////////////////////////////////////

// P2P of the interpolation kernels: scalar matrix kernel (NCMP = 1), single right-hand side
// (NVALS = 1), potentials only
template <class FReal, int NCMP, int NVALS>
struct DirectInteractionComputer
{
    static_assert(NCMP == 1 && NVALS == 1, "only scalar matrix kernels with a single right-hand side");

  template <typename ContainerClass, typename MatrixKernelClass>
  static void P2P( ContainerClass* const FRestrict TargetParticles,
                   ContainerClass* const NeighborSourceParticles[],
                   const int inSize,
                   const MatrixKernelClass *const MatrixKernel){
      FP2PT<FReal>::template FullMutual<ContainerClass,MatrixKernelClass> (TargetParticles,NeighborSourceParticles,inSize,MatrixKernel);
  }

  template <typename ContainerClass, typename MatrixKernelClass>
  static void P2PInner( ContainerClass* const FRestrict TargetParticles,
                   const MatrixKernelClass *const MatrixKernel){
      FP2PT<FReal>::template Inner<ContainerClass, MatrixKernelClass>(TargetParticles,MatrixKernel);
  }

  template <typename ContainerClass, typename MatrixKernelClass>
  static void P2PRemote( ContainerClass* const FRestrict inTargets,
                         const ContainerClass* const inNeighbors[],
                         const int inSize,
                         const MatrixKernelClass *const MatrixKernel){
      FP2PT<FReal>::template FullRemote<ContainerClass,MatrixKernelClass>(inTargets,inNeighbors,inSize,MatrixKernel);
  }
};

#endif // FINTERPP2PKERNELS_HPP
