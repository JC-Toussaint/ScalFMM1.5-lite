// See LICENCE file at project root
//
// Validation of a target/source (Tsm) FMM on fixed particles:
// potentials at the targets, created by charged sources, compared with the direct sum of q/r.

#ifndef TSMCHECK_HPP
#define TSMCHECK_HPP

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "Components/FParticleType.hpp"
#include "Containers/FVector.hpp"
#include "Core/FFmmAlgorithmThreadTsm.hpp"
#include "Utils/FPoint.hpp"
#include "Utils/FTic.hpp"

template <class FReal, class OctreeClass, class CellClass, class ContainerClass, class KernelClass, class LeafClass>
int checkTsm(const char* name, OctreeClass& tree, KernelClass& kernels, const FReal tolerance)
{
    const int NbTargets = 2000, NbSources = 5000;
    const FReal halfWidth = 0.999;

    struct Particle { FReal x, y, z, q; };
    std::vector<Particle> targets(NbTargets), sources(NbSources);
    srand48(12345);
    auto random = [&]() { return (2. * drand48() - 1.) * halfWidth; };

    // indexes: targets in [0, NbTargets), sources in [NbTargets, NbTargets+NbSources)
    for (int i = 0; i < NbTargets; ++i) {
        targets[i] = {random(), random(), random(), 0.};
        tree.insert(FPoint<FReal>(targets[i].x, targets[i].y, targets[i].z), FParticleType::FParticleTypeTarget, i, 0.);
    }
    for (int i = 0; i < NbSources; ++i) {
        sources[i] = {random(), random(), random(), random()};
        tree.insert(FPoint<FReal>(sources[i].x, sources[i].y, sources[i].z), FParticleType::FParticleTypeSource,
                    NbTargets + i, sources[i].q);
    }

    FTic timer;
    FFmmAlgorithmThreadTsm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass> algo(&tree, &kernels);
    algo.execute();
    const double tFmm = timer.tacAndElapsed();

    std::vector<FReal> potential(NbTargets, 0.);
    tree.forEachLeaf([&](LeafClass* leaf) {
        const FReal* const potentials = leaf->getTargets()->getPotentials();
        const FVector<FSize>& indexes = leaf->getTargets()->getIndexes();
        for (FSize idx = 0; idx < leaf->getTargets()->getNbParticles(); ++idx)
            potential[indexes[idx]] = potentials[idx];
    });

    double errL2 = 0., refL2 = 0.;
    for (int i = 0; i < NbTargets; ++i) {
        double direct = 0.;
        for (const Particle& s : sources)
            direct += s.q / std::sqrt((targets[i].x - s.x) * (targets[i].x - s.x) + (targets[i].y - s.y) * (targets[i].y - s.y)
                                      + (targets[i].z - s.z) * (targets[i].z - s.z));
        errL2 += (potential[i] - direct) * (potential[i] - direct);
        refL2 += direct * direct;
    }
    const double relErr = std::sqrt(errL2 / refL2);
    const bool ok = relErr < tolerance;
    std::printf("%s: %d targets, %d sources, FMM %.3f s, relative L2 error on potential %.3e (tolerance %.1e) -> %s\n",
                name, NbTargets, NbSources, tFmm, relErr, double(tolerance), ok ? "OK" : "FAILED");
    return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}

#endif // TSMCHECK_HPP
