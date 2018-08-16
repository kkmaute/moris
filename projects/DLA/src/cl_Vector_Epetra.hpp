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
#include "linalg.hpp"
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
    //Vector_Epetra( const Map_Class * aMap )
    Vector_Epetra( const Map_Class       * aMapClass,
                   const enum VectorType   aVectorType );

    /** Destructor */
    ~Vector_Epetra();

    void replace_global_values();

    void sum_into_global_values(const moris::uint               & aNumMyDofs,
                                const moris::Mat< int >         & aElementTopology,
                                const moris::Mat< moris::real > & aRHSVal);

    void vector_global_asembly();

    void vec_plus_vec( const moris::real & aScaleA,
                             Dist_Vector & aVecA,
                       const moris::real & aScaleThis );

    void scale_vector( const moris::real & aValue,
                       const moris::uint & aVecIndex = 0 );

    void import_local_to_global( const Dist_Vector & aSourceVec );

    void vec_put_scalar( const moris::real & aValue );

    moris::sint vec_local_length() const;

    moris::sint vec_global_length() const;

    moris::real vec_norm2();

    void extract_copy( moris::Mat< moris::real > & LHSValues );

    void extract_my_values( const moris::uint               & aNumIndices,
                            const moris::Mat< moris::sint > & aGlobalRows,
                            const moris::uint               & aRowOffsets,
                                  moris::Mat< moris::real > & LHSValues );

    void save_vector_to_matrix_market_file( const char* aFilename );
//----------------------------------------------------------------------------------------------

    void check_vector();
};
}
#endif /* SRC_DISTLINALG_VECTOREPETRA_HPP_ */
