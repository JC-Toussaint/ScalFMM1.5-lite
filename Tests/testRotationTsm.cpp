// See LICENCE file at project root
//
// Rotation kernel (spherical harmonics, order P) with typed leaves (targets / sources),
// same configuration as in feeLLGood.

#include "Components/FTypedLeaf.hpp"
#include "Containers/FOctree.hpp"
#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"
#include "Kernels/Rotation/FRotationCell.hpp"
#include "Kernels/Rotation/FRotationKernel.hpp"

#include "TsmCheck.hpp"

int main()
{
    typedef double FReal;
    const int P = 9;
    typedef FTypedRotationCell<FReal, P>                                CellClass;
    typedef FP2PParticleContainerIndexed<FReal>                         ContainerClass;
    typedef FTypedLeaf<FReal, ContainerClass>                           LeafClass;
    typedef FOctree<FReal, CellClass, ContainerClass, LeafClass>        OctreeClass;
    typedef FRotationKernel<FReal, CellClass, ContainerClass, P>        KernelClass;

    const int NbLevels = 5, SizeSubLevels = 3;
    const FReal boxWidth = 2.01;
    const FPoint<FReal> centerOfBox(0., 0., 0.);

    OctreeClass tree(NbLevels, SizeSubLevels, boxWidth, centerOfBox);
    KernelClass kernels(NbLevels, boxWidth, centerOfBox);
    return checkTsm<FReal, OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass>("Rotation (P=9)", tree, kernels, 1e-3);
}
