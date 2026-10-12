// See LICENCE file at project root
//
// Validation of a target/source (Tsm) FMM on fixed particles:
// potentials at the targets, created by charged sources, compared with the direct sum of q/r.
// Command line of the tests: <nTargets> <nSources> <seed> (e.g. 2000 5000 12345, as ctest).
// The same seed gives the same targets, sources and charges (drand48).

#ifndef TSMCHECK_HPP
#define TSMCHECK_HPP

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "Components/FParticleType.hpp"
#include "Core/FFmmAlgorithmThreadTsm.hpp"
#include "Utils/FPoint.hpp"
#include "Utils/FTic.hpp"

/** Number of targets and sources and seed of the random generator from the command line:
    <nTargets> <nSources> <seed>, all required. Returns false (after reminding the arguments and
    their order) if they are missing or invalid. */
inline bool readTestArguments(const int argc, char* argv[], int& nTargets, int& nSources, long& seed)
{
    bool ok = (argc == 4);
    long values[3] = {0, 0, 0};
    for (int i = 1; ok && i < argc; ++i) {
        char* end = nullptr;
        values[i - 1] = std::strtol(argv[i], &end, 10);
        ok = (end != argv[i] && *end == '\0' && values[i - 1] >= (i == 3 ? 0 : 1)
              && (i == 3 || values[i - 1] <= 2000000000L));
    }
    if (!ok) {
        std::fprintf(stderr, "usage: %s <nTargets> <nSources> <seed>\n"
                             "  nTargets: number of targets (potential computed by FMM and by direct sum)\n"
                             "  nSources: number of charged sources\n"
                             "  seed:     seed of the random generator (non-negative integer); the same seed gives\n"
                             "            the same positions of targets and sources and the same charges\n"
                             "  example:  %s 2000 5000 12345\n", argv[0], argv[0]);
        return false;
    }
    nTargets = int(values[0]);
    nSources = int(values[1]);
    seed = values[2];
    return true;
}

template <class FReal, class OctreeClass, class CellClass, class ContainerClass, class KernelClass, class LeafClass>
int checkTsm(const char* name, OctreeClass& tree, KernelClass& kernels, const FReal tolerance,
             const int NbTargets, const int NbSources, const long seed)
{
    const FReal halfWidth = 0.999;

    struct Particle { FReal x, y, z, q; };
    std::vector<Particle> targets(NbTargets), sources(NbSources);
    srand48(seed);
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
        const auto& indexes = leaf->getTargets()->getIndexes();
        for (FSize idx = 0; idx < leaf->getTargets()->getNbParticles(); ++idx)
            potential[indexes[idx]] = potentials[idx];
    });

    // reference: direct sum (NbTargets x NbSources interactions, OpenMP)
    double errL2 = 0., refL2 = 0.;
    timer.tic();
#pragma omp parallel for schedule(dynamic, 64) reduction(+ : errL2, refL2)
    for (int i = 0; i < NbTargets; ++i) {
        double direct = 0.;
        for (const Particle& s : sources)
            direct += s.q / std::sqrt((targets[i].x - s.x) * (targets[i].x - s.x) + (targets[i].y - s.y) * (targets[i].y - s.y)
                                      + (targets[i].z - s.z) * (targets[i].z - s.z));
        errL2 += (potential[i] - direct) * (potential[i] - direct);
        refL2 += direct * direct;
    }
    const double tDirect = timer.tacAndElapsed();
    const double relErr = std::sqrt(errL2 / refL2);
    const bool ok = relErr < tolerance;
    std::printf("%s: %d targets, %d sources, seed %ld, FMM %.3f s (direct %.3f s), relative L2 error on potential %.3e (tolerance %.1e) -> %s\n",
                name, NbTargets, NbSources, seed, tFmm, tDirect, relErr, double(tolerance), ok ? "OK" : "FAILED");
    return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}

#endif // TSMCHECK_HPP
