/*
 * Sparse_Matrix_Factory.hpp
 *
 *  Created on: Mar 27, 2018
 *      Author: schmidt
 */
#ifndef SRC_DISTLINALG_SPARSE_MATRIX_FACTORY_HPP_
#define SRC_DISTLINALG_SPARSE_MATRIX_FACTORY_HPP_

#include <memory>
#include "cl_Matrix.hpp"
//#include "linalg_typedefs.hpp"
#include "typedefs.hpp"                       //MRS/COR/src

#include "cl_SOL_Enums.hpp"
//#include "cl_SOL_Dist_Map.hpp"

namespace moris
{
    class Dist_Matrix;
    class Dist_Vector;
    class Dist_Map;
    class Solver_Interface;
    class Matrix_Vector_Factory
    {
    private:
        enum sol::MapType mMapType = sol::MapType::Epetra;
    protected:

    public:
        Matrix_Vector_Factory( const enum sol::MapType aMapType = sol::MapType::Epetra );

        Dist_Matrix * create_matrix(       Solver_Interface * aInput,
                                       const Dist_Map        * aMap );

        Dist_Matrix * create_matrix( const moris::uint aRows,
                                       const moris::uint aCols );

        moris::Dist_Vector * create_vector(       moris::Solver_Interface * aInput,
                                                  moris::Dist_Map        * aMap,
                                            const sint                      aNumVectors = 1 );

        Dist_Vector * create_vector( );

        moris::Dist_Map * create_map( const moris::Matrix< DDSMat > & aMyGlobalIds,
                                      const moris::Matrix< DDUMat > & aMyConstraintIds );

        moris::Dist_Map * create_map( const moris::Matrix< DDSMat > & aMyGlobalIds );
    };
}
#endif /* SRC_DISTLINALG_SPARSE_MATRIX_FACTORY_HPP_ */
