// See LICENCE file at project root
//
// Chebyshev interpolation kernel (1/r, symmetric M2L) with typed leaves (targets / sources).

#include "Components/FTypedLeaf.hpp"
#include "Containers/FOctree.hpp"
#include "Kernels/Chebyshev/FChebCell.hpp"
#include "Kernels/Chebyshev/FChebSymKernel.hpp"
#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"

#include "TsmCheck.hpp"

int main()
{
    typedef double FReal;
    const unsigned int ORDER = 7;
    typedef FTypedChebCell<FReal, ORDER>                                                CellClass;
    typedef FP2PParticleContainerIndexed<FReal>                                         ContainerClass;
    typedef FTypedLeaf<FReal, ContainerClass>                                           LeafClass;
    typedef FOctree<FReal, CellClass, ContainerClass, LeafClass>                        OctreeClass;
    typedef FInterpMatrixKernelR<FReal>                                                 MatrixKernelClass;
    typedef FChebSymKernel<FReal, CellClass, ContainerClass, MatrixKernelClass, ORDER>  KernelClass;

    const int NbLevels = 5, SizeSubLevels = 3;
    const FReal boxWidth = 2.01;
    const FPoint<FReal> centerOfBox(0., 0., 0.);
    const MatrixKernelClass MatrixKernel;

    OctreeClass tree(NbLevels, SizeSubLevels, boxWidth, centerOfBox);
    KernelClass kernels(NbLevels, boxWidth, centerOfBox, &MatrixKernel);
    return checkTsm<FReal, OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass>("Chebyshev (ORDER=7)", tree, kernels, 1e-5);
}
