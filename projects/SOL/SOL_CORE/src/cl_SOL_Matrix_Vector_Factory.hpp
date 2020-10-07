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
#include "typedefs.hpp"
#include "cl_SOL_Enums.hpp"

namespace moris
{
    namespace sol
    {
        class Dist_Matrix;
        class Dist_Vector;
        class Dist_Map;
    }
    class Solver_Interface;
    class Matrix_Vector_Factory
    {
    private:
        enum sol::MapType mMapType = sol::MapType::Epetra;
    protected:

    public:
        Matrix_Vector_Factory( const enum sol::MapType aMapType = sol::MapType::Epetra );

        sol::Dist_Matrix * create_matrix(       Solver_Interface * aInput,
                                          const sol::Dist_Map        * aMap );
										  
		sol::Dist_Matrix * create_matrix( const sol::Dist_Map * aRowMap,
                                          const sol::Dist_Map * aColMap );

        sol::Dist_Matrix * create_matrix( const moris::uint aRows,
                                          const moris::uint aCols );

        moris::sol::Dist_Vector * create_vector(       moris::Solver_Interface * aInput,
                                                       moris::sol::Dist_Map        * aMap,
                                                 const sint                      aNumVectors = 1 );

        moris::sol::Dist_Vector * create_vector(       moris::sol::Dist_Map        * aMap,
                                                 const sint                      aNumVectors = 1 );

        moris::sol::Dist_Map * create_map( const moris::Matrix< DDSMat > & aMyGlobalIds,
                                           const moris::Matrix< DDUMat > & aMyConstraintIds );

        moris::sol::Dist_Map * create_map( const moris::Matrix< DDSMat > & aMyGlobalIds );
    };
}
#endif /* SRC_DISTLINALG_SPARSE_MATRIX_FACTORY_HPP_ */
