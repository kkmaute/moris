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

namespace moris
{
class Sparse_Matrix;
class Dist_Vector;
class Map_Class;
class Solver_Interface;
class Matrix_Vector_Factory
{
private:
protected:

public:
    Matrix_Vector_Factory();

    Sparse_Matrix * create_matrix(       Solver_Interface * aInput,
                                   const Map_Class        * aMap );

    Dist_Vector * create_vector(       Solver_Interface    * aInput,
                                 const Map_Class           * aMap,
                                 const enum VectorType       aVectorType );

    Map_Class * create_map( const moris::uint             & aNumMyDofs,
                            const moris::Matrix< DDSMat > & aMyGlobalElements,
                            const moris::Matrix< DDUMat > & aMyConstraintDofs,
                            const moris::Matrix< DDSMat > & aOverlappingLocaltoGlobalMap );
};
}
#endif /* SRC_DISTLINALG_SPARSE_MATRIX_FACTORY_HPP_ */
