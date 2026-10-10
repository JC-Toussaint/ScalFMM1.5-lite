// See LICENCE file at project root

#ifndef FCHEBCELL_HPP
#define FCHEBCELL_HPP
#include <cstring>
#include <iostream>

#include "../../Components/FBasicCell.hpp"

#include "./FChebTensor.hpp"
#include "../../Extensions/FExtendCellType.hpp"

/**
 * @author Matthias Messner (matthias.messner@inria.fr)
 * @class FChebCell
 * Please read the license
 *
 * This class defines a cell used in the Chebyshev based FMM.
 * @param NVALS is the number of right hand side.
 */
template <class FReal, int ORDER, int NRHS = 1, int NLHS = 1, int NVALS = 1>
class FChebCell : public FBasicCell
{
    // nnodes = ORDER^3
    // we multiply by 2 because we store the  Multipole expansion end the compressed one.
    static const int VectorSize = TensorTraits<ORDER>::nnodes * 2;

    FReal multipole_exp[NRHS * NVALS * VectorSize]; //< Multipole expansion
    FReal         local_exp[NLHS * NVALS * VectorSize]; //< Local expansion

public:
    FChebCell(){
        memset(multipole_exp, 0, sizeof(FReal) * NRHS * NVALS * VectorSize);
        memset(local_exp, 0, sizeof(FReal) * NLHS * NVALS * VectorSize);
    }

    ~FChebCell() {}

    /** Get Multipole */
    const FReal* getMultipole(const int inRhs) const
    {	return this->multipole_exp + inRhs*VectorSize;
    }
    /** Get Local */
    const FReal* getLocal(const int inRhs) const{
        return this->local_exp + inRhs*VectorSize;
    }

    /** Get Multipole */
    FReal* getMultipole(const int inRhs){
        return this->multipole_exp + inRhs*VectorSize;
    }
    /** Get Local */
    FReal* getLocal(const int inRhs){
        return this->local_exp + inRhs*VectorSize;
    }

    /** To get the leading dim of a vec */
    int getVectorSize() const{
        return VectorSize;
    }

    /** Make it like the begining */
    void resetToInitialState(){
        memset(multipole_exp, 0, sizeof(FReal) * NRHS * NVALS * VectorSize);
        memset(local_exp,         0, sizeof(FReal) * NLHS * NVALS * VectorSize);
    }


    //	template <class StreamClass>
    //	const void print(StreamClass& output) const{
    template <class StreamClass>
    friend StreamClass& operator<<(StreamClass& output, const FChebCell<FReal, ORDER, NRHS, NLHS, NVALS>&  cell){
        //	const void print() const{
        output <<"  Multipole exp NRHS " <<NRHS <<" NVALS "  <<NVALS << " VectorSize/2 "  << cell.getVectorSize() *0.5<< std::endl;
        for (int rhs= 0 ; rhs < NRHS ; ++rhs) {
            const FReal* pole = cell.getMultipole(rhs);
            for (int val= 0 ; val < NVALS ; ++val) {
                output<< "      val : " << val << " exp: " ;
                for (int i= 0 ; i < cell.getVectorSize()/2  ; ++i) {
                    output<< pole[i] << " ";
                }
                output << std::endl;
            }
        }
        return output;
    }

};

template <class FReal, int ORDER, int NRHS = 1, int NLHS = 1, int NVALS = 1>
class FTypedChebCell : public FChebCell<FReal, ORDER,NRHS,NLHS,NVALS>, public FExtendCellType {
public:
    void resetToInitialState(){
        FChebCell<FReal,ORDER,NRHS,NLHS,NVALS>::resetToInitialState();
        FExtendCellType::resetToInitialState();
    }


};
#endif //FCHEBCELL_HPP
