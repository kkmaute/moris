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
    class Solver_Interface;

    namespace sol
    {
        class Dist_Matrix;
        class Dist_Vector;
        class Dist_Map;
        
        class Matrix_Vector_Factory
        {
        private:
            enum MapType mMapType = MapType::Epetra;
    
        public:
            Matrix_Vector_Factory(const enum MapType aMapType = MapType::Epetra);
    
            Dist_Matrix *create_matrix(
                    Solver_Interface *aInput,
                    Dist_Map* aMap,
                    bool aPointMap = false,
                    bool aBuildGraph = false);
    
            Dist_Matrix *create_matrix(
                    Dist_Map* aRowMap,
                    Dist_Map* aColMap);
    
            Dist_Matrix *create_matrix(
                    const moris::uint aRows,
                    const moris::uint aCols);
    
            Dist_Vector *create_vector(
                    moris::Solver_Interface *aInput,
                    Dist_Map* aMap,
                    const sint aNumVectors = 1,
                    bool aPointMap = false,
                    bool aManageMap = false);
    
            Dist_Vector *create_vector(
                    Dist_Map* aMap,
                    const sint aNumVectors = 1,
                    bool aPointMap = false,
                    bool aManageMap = false);

            /**
             * Creates a distributed vector/matrix map
             *
             * @param aMyGlobalIds Owned global IDs
             * @param aMyConstraintIds Constraint IDs
             * @return Distributed map
             */
            Dist_Map* create_map(
                    const moris::Matrix<DDSMat> &aMyGlobalIds,
                    const moris::Matrix<DDUMat> &aMyConstraintIds);

            /**
             * Creates a distributed vector/matrix map
             *
             * @param aMyGlobalIds Owned global IDs
             * @return Distributed map
             */
            Dist_Map* create_map(const moris::Matrix<DDSMat> &aMyGlobalIds);
        };
    }
}
#endif /* SRC_DISTLINALG_SPARSE_MATRIX_FACTORY_HPP_ */
