
/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MultiVector_PETSc.hpp
 *
 */

#pragma once

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

//#include "cl_MatrixPETSc.hpp"
#include "cl_SOL_Dist_Vector.hpp"
#include "cl_SOL_Dist_Map.hpp"

#include "cl_Map_PETSc.hpp"
#include "cl_SOL_Dist_Map_Custom.hpp"

namespace moris
{
    class MultiVector_PETSc : public moris::sol::Dist_Vector
    {
      private:
        Dist_Map_Custom* mMap;

        Mat mPetscVector = nullptr;

        moris::Matrix< DDUMat > mDirichletBCVec;

        uint mNumVectors = 1;

        void dirichlet_BC_vector(
                moris::Matrix< DDUMat >&       aDirichletBCVec,
                const moris::Matrix< DDUMat >& aMyConstraintDofs );

      protected:

      public:
        /** Default constructor */
        MultiVector_PETSc(
                moris::Solver_Interface* aInput,
                sol::Dist_Map*           aMap,
                const sint               aNumVectors,
                bool                     aManageMap = false );

        /** Destructor */
        ~MultiVector_PETSc();

        sol::Dist_Map*
        get_map()
        {
            return dynamic_cast< sol::Dist_Map* >( mMap );
        }

        /**
         * Gets a value in the distributed vector based on a given ID.
         *
         * @param aGlobalId Global ID
         * @return Value
         */
        real& operator()( sint aGlobalId, uint aVectorIndex = 0 );

        void sum_into_global_values(
                const moris::Matrix< DDSMat >& aGlobalIds,
                const moris::Matrix< DDRMat >& aValues,
                const uint&                    aVectorIndex = 0 );

        void sum_into_global_values(
                const Vector< sint >&          aGlobalIds,
                const moris::Matrix< DDRMat >& aValues,
                const uint&                    aVectorIndex = 0 );

        void replace_global_values(
                const moris::Matrix< DDSMat >& aGlobalIds,
                const moris::Matrix< DDRMat >& aValues,
                const uint&                    aVectorIndex = 0 );

        void replace_global_values(
                const Vector< sint >& aGlobalIds,
                const Vector< real >& aValues );

        void vector_global_assembly();

        void vec_plus_vec(
                const moris::real& aScaleA,
                sol::Dist_Vector&  aVecA,
                const moris::real& aScaleThis );

        void scale_vector(
                const moris::real& aValue,
                const moris::uint& aVecIndex = 0 );

        void import_local_to_global(
                sol::Dist_Vector& aSourceVec );
        
        //-----------------------------------------------------------------------------
        /**
         * @brief function to add owned source vector to
         * the full vector of the destination vector, used in slepc
         * 
         * @param aSourceVec 
         */
        void import_local_to_global( Vec aSourceVec, uint aVecIndex, Dist_Map_Custom* tSourceMap );

        //-----------------------------------------------------------------------------
        
        void vec_put_scalar( const moris::real& aValue );

        void random();

        moris::sint vec_local_length() const;

        moris::sint vec_global_length() const;

        Vector< moris::real > vec_norm2();

        void extract_copy( moris::Matrix< DDRMat >& LHSValues );

        void extract_copy( Vector< real >& aVector );

        void extract_my_values( const moris::uint&      aNumIndices,
                const moris::Matrix< DDSMat >&          aGlobalBlockRows,
                const moris::uint&                      aBlockRowOffsets,
                Vector< moris::Matrix< DDRMat > >& LHSValues );

        void print() const;

        void save_vector_to_matrix_market_file( const char* aFilename ){};

        void save_vector_to_matlab_file( const char* aFilename );

        void save_vector_to_HDF5( const char* aFilename );

        void read_vector_from_HDF5(
                const char* aFilename,
                std::string aGroupName   = "LHS",
                sint        aVectorindex = 0 );

        void check_vector();

        moris::real* get_values_pointer();

        //-----------------------------------------------------------------------------

        virtual Mat
        get_petsc_vector()
        {
            return mPetscVector;
        };

        virtual Mat
        get_petsc_vector() const
        {
            return mPetscVector;
        };
    };
}    // namespace moris
