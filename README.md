# ScalFMM 1.5 lite

Version épurée de [ScalFMM](https://gitlab.inria.fr/solverstack/ScalFMM) 1.5
(commit `22b9e4f6cf4ea721d71198a71e3f5d2c5ae5e7cc`), pour feeLLGood :

- particules **fixes** : octree `FOctree`, feuilles `FSimpleLeaf` / `FTypedLeaf`
  (pas d'arrangers ni d'algorithmes de déplacement de particules) ;
- algorithmes séquentiels et OpenMP : `FFmmAlgorithm`, `FFmmAlgorithmThread`,
  et leurs variantes cibles/sources `FFmmAlgorithmTsm`, `FFmmAlgorithmThreadTsm` ;
- deux noyaux seulement :
  - **Rotation** (harmoniques sphériques) : `Kernels/Rotation`,
  - **Chebyshev** (interpolation) : `Kernels/Chebyshev` (`FChebKernel`, `FChebSymKernel`,
    `FChebDenseKernel`) et `Kernels/Interpolation` ;
- P2P direct : `Kernels/P2P`.

Modifications par rapport à l'original (performances du noyau Chebyshev, résultats inchangés) :

- `Utils/FBlas.hpp` : les petits produits matriciels (M*N*K ≤ `SCALFMM_SMALL_BLAS_MAX`, 4096 par
  défaut : P2M, M2M, L2L, L2P jusqu'à l'ordre 8) sont faits par des boucles au lieu d'appels BLAS.
  Sur ces tailles l'appel BLAS coûte plus que le calcul, et OpenBLAS (version pthread) sérialise
  les appels simultanés des threads OpenMP : P2M et M2M étaient plus lents à 8 threads qu'à 1 ;
- `FBlas::setSingleThreaded()`, appelé par le constructeur des noyaux Chebyshev : OpenBLAS sur
  un seul thread, BLAS étant appelé depuis les threads OpenMP de l'algorithme (sans cela, M2L à
  l'ordre 8 est jusqu'à 10 fois plus lent).

Supprimés : MPI, StarPU, CUDA/OpenCL, GroupTree, périodicité, FFT et noyaux
Uniform/Taylor/Spherical, lecteurs de fichiers, Addons, documentation, tests unitaires.
Les autres sources conservées sont identiques à l'original (patch `memcpy` de
`FBasicParticleContainer.hpp` inclus).

## Dépendances

CMake ≥ 3.16, compilateur C++14, OpenMP, BLAS et LAPACK (requis par Chebyshev).

## Compilation et installation

    mkdir Build && cd Build
    cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local
    make
    ctest              # validation Rotation et Chebyshev contre le calcul direct
    sudo make install

Installe `include/<Core|Containers|Components|Kernels|Utils|Extensions>/...`,
`include/ScalFmmConfig.h` et `lib/libscalfmm.a` (même disposition que ScalFMM 1.5).

Options : `SCALFMM_USE_AVX` (ON), `SCALFMM_USE_SSE` (OFF), `SCALFMM_USE_NATIVE`
(`-march=native`, ON), `SCALFMM_USE_ASSERT` (ON), `SCALFMM_USE_LOG` (OFF),
`SCALFMM_USE_MEM_STATS` (OFF), `SCALFMM_BLAS_MANGLING` (`ADD_`), `SCALFMM_BUILD_TESTS` (ON).

## Licence

CeCILL-C, voir `LICENCE`.
