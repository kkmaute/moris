/*
 * cl_SOL_Dist_Vector.hpp
 *
 *  Created on: Mar 27, 2018
 *      Author: schmidt
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

#include <petsc.h>
#include <petscis.h>
#include <petscao.h>
#include <petscsys.h>

#include "cl_SOL_Enums.hpp"

namespace moris
{
class Dist_Map;
class Dist_Vector
{
private:
protected:

          Epetra_Import   * mImporter;
          Dist_Map       * mMap;

          moris::sint       mNumVectors;
public:
     Dist_Vector(): mImporter( NULL ),
                    mMap( NULL )
          {};

    Dist_Vector( Dist_Map * aMapClass ): mImporter( NULL ),
                                          mMap( aMapClass )
    {
    };

    virtual  ~Dist_Vector()
    {
    };

    /**
     * @brief Replace global vector entries.
     */
    virtual void replace_global_values() = 0;

    /**
     * @brief Add global valued to the distributed vector.
     *
     * @param[in] aNumMyDof            Number of entries which will be inserted.
     * @param[in] aEleDofConectivity   Position where to place the entriess.
     * @param[in] aRHSVal              Array with values.
     *
     */
    virtual void sum_into_global_values( const moris::uint             & aNumMyDof,
                                         const moris::Matrix< DDSMat > & aEleDofConectivity,
                                         const moris::Matrix< DDRMat > & aRHSVal,
                                         const uint                    & aVectorIndex = 0 ) = 0;

    /**
     * @brief Gather any overlapping/shared data into the non-overlapping partitioning defined by the Map.
     *
     */
    virtual void vector_global_asembly() = 0;

    /**
     * @brief Transfers data from a free vector to a corresponding full vector.
     *
     * @param[in] aSourceVec    Dist_Vector.
     *
     */
    virtual void import_local_to_global( Dist_Vector & aSourceVec) = 0;

    /**
     * @brief Adds the scaled entries of the argument vector to the scaled entries of this vector.
     *
     * @param[in] aScaleA       Scaling value for argument vector.
     * @param[in] aSourceVec    Dist_Vector.
     * @param[in] aScaleThis    Scaling value for this vector.
     *
     */
    virtual void vec_plus_vec( const moris::real & aScaleA,
                                     Dist_Vector & aVecA,
                               const moris::real & aScaleThis ) = 0;

    /**
     * @brief Scales this vector.
     *
     * @param[in] aScaleThis    Scaling value for this vector.
     *
     */
    virtual void scale_vector( const moris::real & aValue,
                               const moris::uint & aVecIndex=0 ) = 0;

    /**
     * @brief Inserts the argument value into all entries of this vector.
     *
     * @param[in] aValue    value to insert.
     *
     */
    virtual void vec_put_scalar( const moris::real & aValue ) = 0;

    /**
     * @brief Returns the local length of this vector.
     *
     * @return Local vector length.
     */
    virtual moris::sint vec_local_length() const = 0;

    /**
     * @brief Returns the global length of this vector.
     *
     * @return Local vector length.
     */
    virtual moris::sint vec_global_length() const = 0;

    /**
     * @brief Returns the euclidean norm of the vector.
     *
     * @return  Euclidean Vector norm.
     */
    virtual moris::real vec_norm2() = 0;

    virtual void print() const = 0;

    virtual void save_vector_to_matrix_market_file( const char* aFilename ) = 0;

    virtual void save_vector_to_HDF5( const char* aFilename ) = 0;

    virtual void read_vector_from_HDF5( const char* aFilename ) = 0;

    virtual void extract_copy( moris::Matrix< DDRMat > & LHSValues ) = 0;

    virtual void extract_my_values( const moris::uint             & aNumIndices,
                                    const moris::Matrix< DDSMat > & aGlobalBlockRows,
                                    const moris::uint             & aBlockRowOffsets,
                                          moris::Matrix< DDRMat > & LHSValues ) = 0;

    virtual moris::real* get_values_pointer()
    {
        MORIS_ASSERT( false, "Not implemented");
        return nullptr;
    };

    virtual moris::real* get_values_pointer() const
    {
        MORIS_ASSERT( false, "Not implemented");
        return nullptr;
    };

    //------------------------------------------------------------------------------------------
    /**
     * @brief Get vector.
     *
     * @return  Vector of type Epetra_Vector or Vec
     */
    virtual Epetra_MultiVector * get_epetra_vector()
    {
    	MORIS_ERROR( false, "get_epetra_vector() function has no child implementation" );
    	return nullptr;
    };

    virtual Epetra_MultiVector * get_epetra_vector() const
    {
    	MORIS_ERROR( false, "get_epetra_vector() function has no child implementation" );
    	return nullptr;
    };

    virtual Vec get_petsc_vector()
    {
    	MORIS_ERROR( false, "get_petsc_vector() function has no child implementation" );
    	return nullptr;
    };

    virtual Vec get_petsc_vector() const
    {
    	MORIS_ERROR( false, "get_petsc_vector() function has no child implementation" );
    	return nullptr;
    };


    const Dist_Map * get_map() { return mMap; };



};
}
#endif /* SRC_DISTLINALG_CL_VECTOR_HPP_ */
