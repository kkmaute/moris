/*
 * VectorEpetra.hpp
 *
 *  Created on: Jan 7, 2018
 *      Author: schmidt
 */
#ifndef SRC_DISTLINALG_VECTOREPETRA_HPP_
#define SRC_DISTLINALG_VECTOREPETRA_HPP_

// MORIS header files.
#ifdef MORIS_HAVE_PARALLEL
 #include <mpi.h>
#endif

// C++ system files
#include <cstdio>
#include <iostream>

// Project header files
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Communication_Tools.hpp" // COM/src

#include "cl_Map_Epetra.hpp"
#include "cl_Vector.hpp"

namespace moris
{
class Vector_Epetra : public Dist_Vector
{
private:
    // Pointer to MultiVector values
    moris::real * mValuesPtr;

protected:

public:
    Vector_Epetra(){};

    Vector_Epetra( const Map_Class       * aMapClass,
                   const enum VectorType   aVectorType );

    /** Destructor */
    ~Vector_Epetra();

    void replace_global_values();

    void sum_into_global_values(const moris::uint               & aNumMyDofs,
                                const moris::Matrix< DDSMat >         & aElementTopology,
                                const moris::Matrix< DDRMat > & aRHSVal);

    void vector_global_asembly();

    void vec_plus_vec( const moris::real & aScaleA,
                             Dist_Vector & aVecA,
                       const moris::real & aScaleThis );

    void scale_vector( const moris::real & aValue,
                       const moris::uint & aVecIndex = 0 );

    void import_local_to_global( Dist_Vector & aSourceVec );

    void vec_put_scalar( const moris::real & aValue );

    moris::sint vec_local_length() const;

    moris::sint vec_global_length() const;

    moris::real vec_norm2();

    void extract_copy( moris::Matrix< DDRMat > & LHSValues );

    void extract_my_values( const moris::uint               & aNumIndices,
                            const moris::Matrix< DDSMat > & aGlobalRows,
                            const moris::uint               & aRowOffsets,
                                  moris::Matrix< DDRMat > & LHSValues );

    void save_vector_to_matrix_market_file( const char* aFilename );

    void save_vector_to_HDF5( const char* aFilename );

    void read_vector_from_HDF5( const char* aFilename );
//----------------------------------------------------------------------------------------------

    void check_vector();
};
}
#endif /* SRC_DISTLINALG_VECTOREPETRA_HPP_ */
