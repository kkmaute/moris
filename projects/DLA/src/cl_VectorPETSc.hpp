/*
 * cl_VectorPETSc.hpp
 *
 *  Created on: Mar 25, 2018
 *      Author: schmidt
 */
#ifndef SRC_DISTLINALG_CL_VECTORPETSC_HPP_
#define SRC_DISTLINALG_CL_VECTORPETSC_HPP_

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

//#include "cl_MatrixPETSc.hpp"
#include "cl_Vector.hpp"
#include "cl_Map_PETSc.hpp"
namespace moris
{

class Vector_PETSc : public moris::Dist_Vector
{
private:

    moris::Matrix< DDUMat >   DirichletBCVec;

    void dirichlet_BC_vector(       moris::Matrix< DDUMat > & aDirichletBCVec,
                              const moris::Matrix< DDUMat > & aMyConstraintDofs );

protected:

public:
    /** Default contructor */
    Vector_PETSc(       moris::Solver_Interface * aInput,
                  const moris::Map_Class        * aMap,
                  const enum moris::VectorType    aVectorType );

    /** Destructor */
    ~Vector_PETSc();

    void sum_into_global_values( const moris::uint             & aNumMyDof,
                                 const moris::Matrix< DDSMat > & aEleDofConectivity,
                                 const moris::Matrix< DDRMat > & aRHSVal );

    void replace_global_values(){};

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

    void extract_my_values( const moris::uint             & aNumIndices,
                            const moris::Matrix< DDSMat > & aGlobalBlockRows,
                            const moris::uint             & aBlockRowOffsets,
                                  moris::Matrix< DDRMat > & LHSValues );

    void save_vector_to_matrix_market_file( const char* aFilename ){};

    void save_vector_to_HDF5( const char* aFilename ){};

    void read_vector_from_HDF5( const char* aFilename ){};

    //-------------------------------------------------------------------------

    void check_vector();
};

}

#endif /* SRC_DISTLINALG_CL_VECTORPETSC_HPP_ */
