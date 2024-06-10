/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SOL_Matrix_Vector_Factory.hpp
 *
 */

#ifndef SRC_DISTLINALG_SPARSE_MATRIX_FACTORY_HPP_
#define SRC_DISTLINALG_SPARSE_MATRIX_FACTORY_HPP_

#include <memory>
#include "cl_Matrix.hpp"
#include "moris_typedefs.hpp"
#include "cl_SOL_Enums.hpp"

namespace moris
{
    class Solver_Interface;

    namespace sol
    {
        class Dist_Matrix;
        class Dist_Vector;
        class Dist_Map;

        class Matrix_Vector_Factory
        {
          private:
            enum MapType            mMapType = MapType::Epetra;
            static Matrix< DDSMat > mDummyMatrix;

          public:
            Matrix_Vector_Factory( const enum MapType aMapType = MapType::Epetra );

            Dist_Matrix *create_matrix(
                    Solver_Interface *aInput,
                    Dist_Map         *aMap,
                    bool              aPointMap   = false,
                    bool              aBuildGraph = false );

            Dist_Matrix *create_matrix(
                    Dist_Map *aRowMap,
                    Dist_Map *aColMap );

            Dist_Matrix *create_matrix(
                    const moris::uint aRows,
                    const moris::uint aCols );

            Dist_Vector *create_vector(
                    moris::Solver_Interface *aInput,
                    Dist_Map                *aMap,
                    const sint               aNumVectors = 1,
                    bool                     aPointMap   = false,
                    bool                     aManageMap  = false );

            Dist_Vector *create_vector(
                    Dist_Map  *aMap,
                    const sint aNumVectors = 1,
                    bool       aPointMap   = false,
                    bool       aManageMap  = false );

            /**
             * Creates a distributed vector/matrix map
             *
             * @param aMyGlobalIds     Owned global IDs
             * @param aMyConstraintIds Constraint IDs
             * @return Distributed map
             */
            Dist_Map *create_map(
                    const moris::Matrix< DDSMat > &aMyGlobalIds,
                    const moris::Matrix< DDUMat > &aMyConstraintIds );

            Dist_Map *create_map(
                    const moris::Matrix< DDSMat > &aMyGlobalIds );

            //------------------------------------------------------------------------------------
            /**
             * Creates a distributed vector/matrix map
             *
             * @param aMyGlobalIds Owned global IDs
             * @return Distributed map
             */
            Dist_Map *create_map( const moris::Matrix< DDSMat > &aMyGlobalIds,
                    const moris::Matrix< DDSMat >               &aMyGlobalOwnedAndSharedIds );
            //------------------------------------------------------------------------------------
            /**
             * Creates a distributed vector/matrix map
             *
             * @param aMyGlobalIds Owned global IDs
             * @return Distributed map
             */
            Dist_Map* create_map( const Vector< sint >& aMyGlobalIds );

            /**
             * Creates a distributed vector map specifically for full vector
             *
             * @param aMyGlobalIds Owned global IDs
             * @param aMyGlobalIds Owned and Shared global IDs
             * @return Distributed map
             */
            Dist_Map *create_full_map(
                    const moris::Matrix< DDSMat > &aMyGlobalOwnedIds,
                    const moris::Matrix< DDSMat > &aMyGlobalOwnedAndSharedIds );
        };
    }    // namespace sol
}    // namespace moris
#endif /* SRC_DISTLINALG_SPARSE_MATRIX_FACTORY_HPP_ */
