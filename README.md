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
- P2P direct : `Kernels/P2P` ;
- **potentiel seul** (feeLLGood n'utilise pas les forces) : le conteneur de particules ne garde que
  la charge et le potentiel (2 valeurs par particule au lieu de 5), P2P et L2P ne calculent plus
  les forces. Les opérations du potentiel sont inchangées : potentiels identiques au bit près à
  ceux de ScalFMM 1.5. Noyaux d'interpolation : matrice scalaire et un seul second membre.

Modifications par rapport à l'original (performances du noyau Chebyshev, résultats inchangés) :

- `Utils/FBlas.hpp` : les petits produits matriciels (M*N*K ≤ `SCALFMM_SMALL_BLAS_MAX`, 4096 par
  défaut : P2M, M2M, L2L, L2P jusqu'à l'ordre 8) sont faits par des boucles au lieu d'appels BLAS.
  Sur ces tailles l'appel BLAS coûte plus que le calcul, et OpenBLAS (version pthread) sérialise
  les appels simultanés des threads OpenMP : P2M et M2M étaient plus lents à 8 threads qu'à 1 ;
- `FBlas::setSingleThreaded()`, appelé par le constructeur des noyaux Chebyshev : OpenBLAS sur
  un seul thread, BLAS étant appelé depuis les threads OpenMP de l'algorithme (sans cela, M2L à
  l'ordre 8 est jusqu'à 10 fois plus lent).

Bibliothèque standard (C++17) à la place des équivalents maison, résultats identiques au bit près :

| Avant | Maintenant |
|---|---|
| `FVector` | `std::vector` (`FVector.hpp` n'est plus qu'un alias, pour le code existant) |
| `FSmartPointer` | `std::shared_ptr` (données précalculées partagées par les copies des noyaux) |
| `FNoCopyable` | constructeur de copie et affectation `= delete` |
| `FMemUtils` | `std::copy_n`, `std::transform` |
| `FAlignedMemory` | `operator new(taille, std::align_val_t)` |
| `FBasicBlockAllocator` / `FListBlockAllocator` | `new` / `delete` (plus de paramètre d'allocateur pour `FOctree`) |
| `FEnv` | `std::getenv` (`SCALFMM_ALGO_NUM_THREADS`) |
| `FTic` | réduit à `std::chrono` |
| `FMath::Sin`, `Cos`, `Sqrt`, `Max`, `Min`... | `std::sin`, `std::cos`, `std::sqrt`, `std::max`, `std::min`... |
| `FComplex` | `std::complex` (noyau Rotation) |

Conservés faute d'équivalent standard : l'abstraction SIMD de `FMath` (P2P vectorisé SSE/AVX),
`FMath::pow(x, n)` à exposant entier (n multiplications ; `std::pow` arrondit autrement),
`FPoint` (déjà fondé sur `std::array`), `FAssert`. Dans le noyau Rotation, le produit de deux
`std::complex` est écrit explicitement (`FRotationKernel::mul`) : sans `-ffast-math`, l'opérateur
`*=` de `std::complex` ajoute à chaque produit un test NaN et un appel possible à `__muldc3`.

Octree allégé (idées reprises du scalFMMlight de feeLLGood) : sans périodicité ni recherches de
voisins inutilisées, `FTreeCoordinate` autonome (ne dérive plus de `FPoint`).

Supprimés : MPI (dont la sérialisation des cellules et conteneurs), forces (P2P multi-seconds
membres et tensoriel, L2P du gradient), StarPU, CUDA/OpenCL,
GroupTree, périodicité, FFT et noyaux Uniform/Taylor/Spherical, lecteurs de fichiers, Addons,
journal de débogage (`FLog`) et statistiques mémoire (`FMemStats`), documentation, tests unitaires,
bibliothèque compilée `libscalfmm.a` (en-têtes seuls).

## Dépendances

CMake ≥ 3.16, compilateur C++17, OpenMP, BLAS et LAPACK (requis par Chebyshev).

BLAS est appelé depuis les threads OpenMP de l'algorithme : il doit être séquentiel.
Sans `-DBLA_VENDOR=...`, CMake prend OpenBLAS, sinon MKL séquentiel (`Intel10_64lp_seq`),
sinon le premier BLAS trouvé. Ne pas utiliser MKL multithread (`Intel10_64lp`) avec GCC :
son runtime OpenMP (Intel) est incompatible avec celui de GCC (plantage dans `dgeqrf`).

## Compilation et installation

La bibliothèque est faite uniquement d'en-têtes : rien à compiler, `make` ne construit que les
tests. CMake génère `ScalFmmConfig.h` (options ci-dessous).

    mkdir Build && cd Build
    cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local      # ou $HOME/local, sans sudo
    make
    ctest              # validation Rotation et Chebyshev contre le calcul direct
    sudo make install

Installe `include/<Core|Containers|Components|Kernels|Utils|Extensions>/...` et
`include/ScalFmmConfig.h` (même disposition que ScalFMM 1.5), ainsi que la cible CMake
`scalfmm::scalfmm` (`lib/cmake/ScalFMM`). Il suffit au code utilisateur de mettre
`<prefix>/include` dans ses chemins d'inclusion et de lier OpenMP, BLAS et LAPACK.

Options : `SCALFMM_USE_AVX` (ON), `SCALFMM_USE_SSE` (OFF), `SCALFMM_USE_NATIVE`
(`-march=native`, ON), `SCALFMM_USE_ASSERT` (ON),
`SCALFMM_BLAS_MANGLING` (`ADD_`), `SCALFMM_BUILD_TESTS` (ON).

## Licence

CeCILL-C, voir `LICENCE`.
