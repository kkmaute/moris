/*
 * cl_VectorPETSc.hpp
 *
 *  Created on: Mar 25, 2018
 *      Author: schmidt
 */
#ifndef SRC_DISTLINALG_CL_VECTORPETSC_HPP_
#define SRC_DISTLINALG_CL_VECTORPETSC_HPP_

#include "linalg.hpp"

//#include "cl_MatrixPETSc.hpp"
#include "cl_Vector.hpp"
#include "cl_Map_PETSc.hpp"

class Vector_PETSc : public moris::Dist_Vector
{
private:

    moris::Mat< moris::uint >   DirichletBCVec;

    void dirichlet_BC_vector(       moris::Mat< moris::uint > & aDirichletBCVec,
                              const moris::Mat< moris::uint > & aMyConstraintDofs );

protected:

public:

    Vector_PETSc(       moris::Solver_Input * aInput,
                  const moris::Map_Class    * aMap,
                  const enum moris::VectorType       aVectorType );


    ~Vector_PETSc();

    void sum_into_global_values( const moris::uint              & aNumMyDof,
                                 const moris::Mat< int >        & aEleDofConectivity,
                                 const moris::Mat< moris::real >& aRHSVal );

    void replace_global_values(){};

    void vector_global_asembly();

    void vec_plus_vec( const moris::real & aScaleA,
                             Dist_Vector & aVecA,
                       const moris::real & aScaleThis );

    void scale_vector( const moris::real & aValue,
                       const moris::uint & aVecIndex = 0 );

    void import_local_to_global( const Dist_Vector & aSourceVec ){};

    void vec_put_scalar( const moris::real & aValue );

    moris::sint vec_local_length() const;

    moris::sint vec_global_length() const;

    moris::real vec_norm2();

    void save_vector_to_matrix_market_file( const char* aFilename ){};

    //-------------------------------------------------------------------------

    void check_vector();
};

#endif /* SRC_DISTLINALG_CL_VECTORPETSC_HPP_ */
