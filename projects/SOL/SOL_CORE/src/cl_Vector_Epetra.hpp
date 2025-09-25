/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Vector_Epetra.hpp
 *
 */

#ifndef SRC_DISTLINALG_VECTOREPETRA_HPP_
#define SRC_DISTLINALG_VECTOREPETRA_HPP_

// MORIS header files.
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Communication_Tools.hpp"    // COM/src

// C++ system files
#include <cstdio>
#include <iostream>

// Project header files
#include "cl_Map_Epetra.hpp"
#include "cl_SOL_Dist_Vector.hpp"
#include "cl_SOL_Dist_Map.hpp"

namespace moris
{
    class Vector_Epetra : public sol::Dist_Vector
    {
      private:
        Map_Epetra* mMap = nullptr;

        Epetra_MultiVector* mEpetraVector = nullptr;

        // Pointer to MultiVector values
        moris::real* mValuesPtr;

        const bool mVecBuildWithPointMap;

      protected:

      public:
        Vector_Epetra(
                sol::Dist_Map* aMapClass,
                const sint     aNumVectors,
                bool           aPointMap  = false,
                bool           aManageMap = false );

        /** Destructor */
        ~Vector_Epetra() override;

        /**
         * Gets a value in the distributed vector based on a given ID.
         *
         * @param aGlobalId Global ID
         * @return Value
         */
        real& operator()( sint aGlobalId, uint aVectorIndex = 0 ) override;

        sol::Dist_Map*
        get_map() override
        {
            return dynamic_cast< sol::Dist_Map* >( mMap );
        }

        Epetra_MultiVector*
        get_epetra_vector()
        {
            return mEpetraVector;
        }
        Epetra_MultiVector*
        get_epetra_vector() const
        {
            return mEpetraVector;
        }

        void replace_global_values(
                const moris::Matrix< DDSMat >& aGlobalIds,
                const moris::Matrix< DDRMat >& aValues,
                const uint&                    aVectorIndex = 0 ) override;

        void replace_global_values(
                const Vector< sint >& aGlobalIds,
                const Vector< real >& aValues ) override;

        void sum_into_global_values(
                const moris::Matrix< DDSMat >& aGlobalIds,
                const moris::Matrix< DDRMat >& aValues,
                const uint&                    aVectorIndex = 0 ) override;

        void sum_into_global_values(
                const Vector< sint >&          aGlobalIds,
                const moris::Matrix< DDRMat >& aValues,
                const uint&                    aVectorIndex = 0 ) override;

        void vector_global_assembly() override;

        void vec_plus_vec(
                const moris::real& aScaleA,
                sol::Dist_Vector&  aVecA,
                const moris::real& aScaleThis ) override;

        void scale_vector(
                const moris::real& aValue,
                const moris::uint& aVecIndex = 0 ) override;

        void import_local_to_global( sol::Dist_Vector& aSourceVec ) override;

        void vec_put_scalar( const moris::real& aValue ) override;

        void random() override;

        moris::sint vec_local_length() const override;

        moris::sint vec_global_length() const override;

        Vector< moris::real > vec_norm2() override;

        void extract_copy( moris::Matrix< DDRMat >& LHSValues ) override;

        void extract_copy( Vector< real >& aVector ) override;

        void extract_my_values(
                const moris::uint&                      aNumIndices,
                const moris::Matrix< DDSMat >&          aGlobalRows,
                const moris::uint&                      aRowOffsets,
                Vector< moris::Matrix< DDRMat > >& LHSValues ) override;

        void extract_my_values(
        const moris::uint&                      aNumIndices,
        const Vector< sint >&                   aGlobalRows,
        const moris::uint&                      aRowOffsets,
        Vector< moris::Matrix< DDRMat > >& LHSValues ) override;


        void print() const override;

        void save_vector_to_matrix_market_file( const char* aFilename ) override;

        void save_vector_to_matlab_file( const char* aFilename ) override;

        void save_vector_to_HDF5( const char* aFilename ) override;

        void read_vector_from_HDF5(
                const char* aFilename,
                std::string aGroupName   = "LHS",
                sint        aVectorindex = 0 ) override;

        //----------------------------------------------------------------------------------------------

        moris::real*
        get_values_pointer() override
        {
            return mValuesPtr;
        };

        //----------------------------------------------------------------------------------------------

        void check_vector();
    };
}    // namespace moris
#endif /* SRC_DISTLINALG_VECTOREPETRA_HPP_ */
