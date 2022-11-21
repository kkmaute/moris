/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SOL_Dist_Vector.hpp
 *
 */

#ifndef SRC_DISTLINALG_CL_VECTOR_HPP_
#define SRC_DISTLINALG_CL_VECTOR_HPP_

// MORIS header files.
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

// TPL header files
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_FEVector.h"
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Import.h"

#include "EpetraExt_MultiVectorIn.h"
#include "EpetraExt_VectorIn.h"
#include "EpetraExt_VectorOut.h"
#include "EpetraExt_MultiVectorOut.h"
#include <EpetraExt_HDF5.h>

#ifdef MORIS_HAVE_PETSC
#include <petsc.h>
#include <petscis.h>
#include <petscao.h>
#include <petscsys.h>
#endif

#include "cl_SOL_Enums.hpp"

namespace moris
{
    namespace sol
    {
        class Dist_Map;

        class Dist_Vector
        {
          protected:
            Epetra_Import* mImporter = nullptr;
            Dist_Map*      mMap;
            bool           mManageMap;
            moris::sint    mNumVectors = 1;

          public:
            /**
             * Constructor
             *
             * @param aMapClass Distributed vector map
             */
            Dist_Vector( sol::Dist_Map* aMapClass, bool aManageMap );

            /**
             * Destructor (deletes map if desired)
             */
            virtual ~Dist_Vector();

            /**
             * Get distributed vector map
             *
             * @return Map
             */
            sol::Dist_Map* get_map();

            /**
             * Gets a value in the distributed vector based on a given ID.
             *
             * @param aGlobalId Global ID
             * @return Value
             */
            virtual real& operator()( sint aGlobalId, uint aVectorIndex = 0 ) = 0;

            /**
             * Returns the number of vectors.
             */
            moris::sint
            get_num_vectors()
            {
                return mNumVectors;
            };

            /**
             * Replace global vector entries.
             */
            virtual void replace_global_values(
                    const moris::Matrix< DDSMat >& aGlobalIds,
                    const moris::Matrix< DDRMat >& aValues,
                    const uint&                    aVectorIndex = 0 ) = 0;

            /**
             * Add global valued to the distributed vector.
             *
             * @param[in] aNumMyDof            Number of entries which will be inserted.
             * @param[in] aEleDofConectivity   Position where to place the entriess.
             * @param[in] aRHSVal              Array with values.
             *
             */
            virtual void sum_into_global_values(
                    const moris::Matrix< DDSMat >& aGlobalIds,
                    const moris::Matrix< DDRMat >& aValues,
                    const uint&                    aVectorIndex = 0 ) = 0;

            /**
             * Gather any overlapping/shared data into the non-overlapping partitioning defined by the Map.
             */
            virtual void vector_global_assembly() = 0;

            /**
             * Transfers data from a free vector to a corresponding full vector.
             *
             * @param[in] aSourceVec    Dist_Vector.
             */
            virtual void import_local_to_global( Dist_Vector& aSourceVec ) = 0;

            /**
             * Adds the scaled entries of the argument vector to the scaled entries of this vector.
             *
             * @param[in] aScaleA       Scaling value for argument vector.
             * @param[in] aSourceVec    Dist_Vector.
             * @param[in] aScaleThis    Scaling value for this vector.
             */
            virtual void vec_plus_vec(
                    const moris::real& aScaleA,
                    Dist_Vector&       aVecA,
                    const moris::real& aScaleThis ) = 0;

            /**
             * Scales this vector.
             *
             * @param[in] aScaleThis    Scaling value for this vector.
             */
            virtual void scale_vector(
                    const moris::real& aValue,
                    const moris::uint& aVecIndex = 0 ) = 0;

            /**
             * Inserts the argument value into all entries of this vector.
             *
             * @param[in] aValue    value to insert.
             */
            virtual void vec_put_scalar( const moris::real& aValue ) = 0;

            /**
             * Inserts random values into vactor.
             *
             */
            virtual void random() = 0;

            /**
             * Returns the local length of this vector.
             *
             * @return Local vector length.
             */
            virtual moris::sint vec_local_length() const = 0;

            /**
             * Returns the global length of this vector.
             *
             * @return Local vector length.
             */
            virtual moris::sint vec_global_length() const = 0;

            /**
             * Returns the euclidean norm of the vector.
             *
             * @return  Euclidean Vector norm.
             */
            virtual Cell< moris::real > vec_norm2() = 0;

            /**
             * Prints this vector.
             */
            virtual void print() const = 0;

            /**
             * Saves this vector to a Matrix Market file.
             *
             * @param aFilename File name to save to
             */
            virtual void save_vector_to_matrix_market_file( const char* aFilename ) = 0;

            /**
             * Saves this vector to a Matlab file.
             *
             * @param aFilename File name to save to
             */
            virtual void save_vector_to_matlab_file( const char* aFilename ) = 0;

            /**
             * Saves this vector to an HDF5 file.
             *
             * @param aFilename File name to save to
             */
            virtual void save_vector_to_HDF5( const char* aFilename ) = 0;

            /**
             * Reads a vector stored in an HDF5 file into this class.
             *
             * @param aFilename File name to read from
             */
            virtual void read_vector_from_HDF5( const char* aFilename ) = 0;

            /**
             * Extracts a full copy of this vector into a DDRMat format.
             *
             * @param LHSValues Matrix to extract into
             */
            virtual void extract_copy( moris::Matrix< DDRMat >& LHSValues ) = 0;

            /**
             * Extracts owned values in this vector into a DDRMat.
             *
             * @param aNumIndices Number of indices to extract
             * @param aGlobalBlockRows Global block rows
             * @param aBlockRowOffsets Global row offsets
             * @param LHSValues Matrix to extract into
             */
            virtual void extract_my_values(
                    const moris::uint&                      aNumIndices,
                    const moris::Matrix< DDSMat >&          aGlobalBlockRows,
                    const moris::uint&                      aBlockRowOffsets,
                    moris::Cell< moris::Matrix< DDRMat > >& LHSValues ) = 0;

            /**
             * Gets a pointer to the real values stored in this vector.
             *
             * @return real pointer
             */
            virtual moris::real* get_values_pointer() = 0;
        };
    }    // namespace sol
}    // namespace moris
#endif /* SRC_DISTLINALG_CL_VECTOR_HPP_ */
