/*
 * Sparse_Matrix_Factory.hpp
 *
 *  Created on: Mar 27, 2018
 *      Author: schmidt
 */
#ifndef SRC_DISTLINALG_SPARSE_MATRIX_FACTORY_HPP_
#define SRC_DISTLINALG_SPARSE_MATRIX_FACTORY_HPP_

#include <memory>
//#include "cl_DLA_Enums.hpp"
//#include "linalg_typedefs.hpp"
//#include "cl_Matrix.hpp"

#include "cl_DLA_Enums.hpp"

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
        enum MapType mMapType = MapType::Epetra;
    protected:

    public:
        Matrix_Vector_Factory( const enum MapType aMapType = MapType::Epetra );

        Sparse_Matrix * create_matrix(       Solver_Interface * aInput,
                                       const Map_Class        * aMap );

        Sparse_Matrix * create_matrix( const moris::uint aRows,
                                       const moris::uint aCols );

        Dist_Vector * create_vector(       Solver_Interface    * aInput,
                                     const Map_Class           * aMap,
                                     const enum VectorType       aVectorType );

        Dist_Vector * create_vector( );

        Map_Class * create_map( const moris::uint             & aNumMaxDofs,
                                const moris::Matrix< DDSMat > & aMyGlobalElements,
                                const moris::Matrix< DDUMat > & aMyConstraintDofs,
                                const moris::Matrix< DDSMat > & aOverlappingLocaltoGlobalMap );
    };
}
#endif /* SRC_DISTLINALG_SPARSE_MATRIX_FACTORY_HPP_ */
