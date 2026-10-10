// See LICENCE file at project root
#ifndef FFMMALGORITHMTSM_HPP
#define FFMMALGORITHMTSM_HPP


#include <algorithm>
#include "../Utils/FAssert.hpp"

#include "../Utils/FTic.hpp"

#include "../Containers/FOctree.hpp"
#include "FCoreCommon.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FFmmAlgorithmTsm
* @brief
* Please read the license
*
* This class is a basic FMM algorithm
* It just iterates on a tree and call the kernels with good arguments.
*
* Of course this class does not deallocate pointer given in arguements.
*
* The differences with FmmAlgorithm is that it used target source model.
*/
template<class OctreeClass, class CellClass, class ContainerClass, class KernelClass, class LeafClass>
class FFmmAlgorithmTsm : public FAbstractAlgorithm{

    OctreeClass* const tree;                                                     //< The octree to work on
    KernelClass* const kernels;    //< The kernels

    const int OctreeHeight;

    const int leafLevelSeparationCriteria;


public:
    /** The constructor need the octree and the kernels used for computation
      * @param inTree the octree to work on
      * @param inKernels the kernels to call
      * An assert is launched if one of the arguments is null
      */
    FFmmAlgorithmTsm(OctreeClass* const inTree, KernelClass* const inKernels, const int inLeafLevelSeperationCriteria = 1)
        : tree(inTree) , kernels(inKernels) , OctreeHeight(tree->getHeight()), leafLevelSeparationCriteria(inLeafLevelSeperationCriteria){

        FAssertLF(tree, "tree cannot be null");
        FAssertLF(kernels, "kernels cannot be null");
        FAssertLF(leafLevelSeparationCriteria < 3, "Separation criteria should be < 3");

        FAbstractAlgorithm::setNbLevelsInTree(tree->getHeight());

    }

    /** Default destructor */
    virtual ~FFmmAlgorithmTsm(){
    }

protected:
    /**
      * To execute the fmm algorithm
      * Call this function to run the complete algorithm
      */
    void executeCore(const unsigned operationsToProceed) override {

        if(operationsToProceed & FFmmP2M) bottomPass();

        if(operationsToProceed & FFmmM2M) upwardPass();

        if(operationsToProceed & FFmmM2L) transferPass();

        if(operationsToProceed & FFmmL2L) downardPass();

        if((operationsToProceed & FFmmP2P) || (operationsToProceed & FFmmL2P)) directPass((operationsToProceed & FFmmP2P),(operationsToProceed & FFmmL2P));
    }

    /** P2M */
    void bottomPass(){

        typename OctreeClass::Iterator octreeIterator(tree);

        // Iterate on leafs
        octreeIterator.gotoBottomLeft();
        do{
            // We need the current cell that represent the leaf
            // and the list of particles
            ContainerClass* const sources = octreeIterator.getCurrentListSrc();
            if(sources->getNbParticles()){
                octreeIterator.getCurrentCell()->setSrcChildTrue();
                kernels->P2M( octreeIterator.getCurrentCell() , sources);
            }
            if(octreeIterator.getCurrentListTargets()->getNbParticles()){
                octreeIterator.getCurrentCell()->setTargetsChildTrue();
            }
        } while(octreeIterator.moveRight());


    }

    /** M2M */
    void upwardPass(){

        // Start from leal level - 1
        typename OctreeClass::Iterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        octreeIterator.moveUp();

        for(int idxLevel = OctreeHeight - 2 ; idxLevel > FAbstractAlgorithm::lowerWorkingLevel-1 ; --idxLevel){
            octreeIterator.moveUp();
        }

        typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

        // for each levels
        for(int idxLevel = std::min(OctreeHeight - 2, FAbstractAlgorithm::lowerWorkingLevel - 1) ; idxLevel >= FAbstractAlgorithm::upperWorkingLevel ; --idxLevel ){
            // for each cells
            do{
                // We need the current cell and the child
                // child is an array (of 8 child) that may be null

                CellClass* potentialChild[8];
                CellClass** const realChild = octreeIterator.getCurrentChild();
                CellClass* const currentCell = octreeIterator.getCurrentCell();
                for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
                    potentialChild[idxChild] = nullptr;
                    if(realChild[idxChild]){
                        if(realChild[idxChild]->hasSrcChild()){
                            currentCell->setSrcChildTrue();
                            potentialChild[idxChild] = realChild[idxChild];
                        }
                        if(realChild[idxChild]->hasTargetsChild()){
                            currentCell->setTargetsChildTrue();
                        }
                    }
                }
                kernels->M2M( currentCell , potentialChild, idxLevel);

            } while(octreeIterator.moveRight());

            avoidGotoLeftIterator.moveUp();
            octreeIterator = avoidGotoLeftIterator;// equal octreeIterator.moveUp(); octreeIterator.gotoLeft();
        }


    }

    /** M2L */
    void transferPass(){

        typename OctreeClass::Iterator octreeIterator(tree);
        octreeIterator.moveDown();

        for(int idxLevel = 2 ; idxLevel < FAbstractAlgorithm::upperWorkingLevel ; ++idxLevel){
            octreeIterator.moveDown();
        }

        typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

        const CellClass* neighbors[342];
        int neighborPositions[342];

        // for each levels
        for(int idxLevel = FAbstractAlgorithm::upperWorkingLevel ; idxLevel < FAbstractAlgorithm::lowerWorkingLevel ; ++idxLevel ){
            const int separationCriteria = (idxLevel != FAbstractAlgorithm::lowerWorkingLevel-1 ? 1 : leafLevelSeparationCriteria);
            // for each cells
            do{
                CellClass* const currentCell = octreeIterator.getCurrentCell();

                if(currentCell->hasTargetsChild()){
                    const int counter = tree->getInteractionNeighbors(neighbors, neighborPositions, octreeIterator.getCurrentGlobalCoordinate(),idxLevel, separationCriteria);
                    if( counter ){
                        int counterWithSrc = 0;
                        for(int idxRealNeighbors = 0 ; idxRealNeighbors < counter ; ++idxRealNeighbors ){
                            if(neighbors[idxRealNeighbors]->hasSrcChild()){
                                neighbors[counterWithSrc] = neighbors[idxRealNeighbors];
                                neighborPositions[counterWithSrc] = neighborPositions[idxRealNeighbors];
                                ++counterWithSrc;
                            }
                        }
                        if(counterWithSrc){
                            kernels->M2L( currentCell , neighbors, neighborPositions, counterWithSrc, idxLevel);
                        }
                    }
                }
            } while(octreeIterator.moveRight());

            kernels->finishedLevelM2L(idxLevel);

            avoidGotoLeftIterator.moveDown();
            octreeIterator = avoidGotoLeftIterator;
        }

    }

    /** L2L */
    void downardPass(){

        typename OctreeClass::Iterator octreeIterator(tree);
        octreeIterator.moveDown();

        for(int idxLevel = 2 ; idxLevel < FAbstractAlgorithm::upperWorkingLevel ; ++idxLevel){
            octreeIterator.moveDown();
        }

        typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

        const int heightMinusOne = FAbstractAlgorithm::lowerWorkingLevel - 1;
        // for each levels exepted leaf level
        for(int idxLevel = FAbstractAlgorithm::upperWorkingLevel ; idxLevel < heightMinusOne ; ++idxLevel ){
            // for each cells
            do{
                if( octreeIterator.getCurrentCell()->hasTargetsChild() ){
                    CellClass* potentialChild[8];
                    CellClass** const realChild = octreeIterator.getCurrentChild();
                    CellClass* const currentCell = octreeIterator.getCurrentCell();
                    for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
                        if(realChild[idxChild] && realChild[idxChild]->hasTargetsChild()){
                            potentialChild[idxChild] = realChild[idxChild];
                        }
                        else{
                            potentialChild[idxChild] = nullptr;
                        }
                    }
                    kernels->L2L( currentCell , potentialChild, idxLevel);
                }
            } while(octreeIterator.moveRight());

            avoidGotoLeftIterator.moveDown();
            octreeIterator = avoidGotoLeftIterator;
        }

    }



    /** P2P */
    void directPass(const bool p2pEnabled, const bool l2pEnabled){

        const int heightMinusOne = OctreeHeight - 1;

        typename OctreeClass::Iterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        // There is a maximum of 26 neighbors
        ContainerClass* neighbors[26];
        int neighborPositions[26];
        // for each leafs
        do{
            if( octreeIterator.getCurrentCell()->hasTargetsChild() ){
                if(l2pEnabled){
                    kernels->L2P(octreeIterator.getCurrentCell(), octreeIterator.getCurrentListTargets());
                }
                if(p2pEnabled){
                    if(octreeIterator.getCurrentCell()->hasSrcChild()){
                        kernels->P2P( octreeIterator.getCurrentGlobalCoordinate(), octreeIterator.getCurrentListTargets(),
                                      octreeIterator.getCurrentListSrc() , neighbors, neighborPositions, 0);
                    }
                    // need the current particles and neighbors particles
                    const int counter = tree->getLeafsNeighbors(neighbors, neighborPositions, octreeIterator.getCurrentGlobalCoordinate(), heightMinusOne);
                    kernels->P2PRemote( octreeIterator.getCurrentGlobalCoordinate(), octreeIterator.getCurrentListTargets(),
                              octreeIterator.getCurrentListSrc() , neighbors, neighborPositions, counter);
                }
            }
        } while(octreeIterator.moveRight());


    }

};


#endif //FFMMALGORITHMTSM_HPP


