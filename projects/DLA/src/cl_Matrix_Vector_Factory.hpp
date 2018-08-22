/*
 * Sparse_Matrix_Factory.hpp
 *
 *  Created on: Mar 27, 2018
 *      Author: schmidt
 */
#ifndef SRC_DISTLINALG_SPARSE_MATRIX_FACTORY_HPP_
#define SRC_DISTLINALG_SPARSE_MATRIX_FACTORY_HPP_

#include <memory>
#include "cl_Map_Class.hpp"
//#include "cl_Solver_Input.hpp"

class Sparse_Matrix;

namespace moris
{
class Dist_Vector;
class Map_Class;
class Solver_Input;
class Matrix_Vector_Factory
{
private:
protected:

public:
    Matrix_Vector_Factory();

    Sparse_Matrix * create_matrix(       Solver_Input * aInput,
                                   const Map_Class    * aMap );

    Dist_Vector * create_vector(       Solver_Input    * aInput,
                                 const Map_Class       * aMap,
                                 const enum VectorType   aVectorType );

    Map_Class * create_map( const moris::uint        & aNumMyDofs,
                            const moris::Mat< int >  & aMyGlobalElements,
                            const moris::Mat< uint > & aMyConstraintDofs,
                            const moris::Mat< int >  & aOverlappingLocaltoGlobalMap );
};
}
#endif /* SRC_DISTLINALG_SPARSE_MATRIX_FACTORY_HPP_ */
