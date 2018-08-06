/*
 * cl_Matrix_Vector_Factory.cpp
 *
 *  Created on: Jun 28, 2018
 *      Author: schmidt
 */
#include "cl_Matrix_Vector_Factory.hpp"

#include "cl_Sparse_Matrix_EpetraFECrs.hpp"
#include "cl_MatrixPETSc.hpp"
#include "cl_Vector_Epetra.hpp"
#include "cl_VectorPETSc.hpp"
#include "cl_Map_Epetra.hpp"
#include "cl_Map_PETSc.hpp"

moris::Matrix_Vector_Factory::Matrix_Vector_Factory()
{
}

Sparse_Matrix * moris::Matrix_Vector_Factory::create_matrix(       moris::Solver_Input * aInput,
                                                             const moris::Map_Class    * aMap )
{
    Sparse_Matrix * tSparseMatrix;

    switch(0)
    {
    case (0):
        tSparseMatrix = new Sparse_Matrix_EpetraFECrs( aInput, aMap );
        break;
//    case (1):
//        tSparseMatrix = new Matrix_PETSc( aInput, aMap );
//        break;
    default:
        MORIS_ASSERT( false, "No matrix type specified." );
        break;
    }
    return tSparseMatrix;
}

//-------------------------------------------------------------------------------------------------
moris::Dist_Vector * moris::Matrix_Vector_Factory::create_vector(       moris::Solver_Input * aInput,
                                                                  const moris::Map_Class    * aMap,
                                                                  const enum VectorType       aVectorType )
{
moris::Dist_Vector * tDistVector;

    switch(0)
    {
    case (0):
        tDistVector = new moris::Vector_Epetra( aMap, aVectorType );
        break;
//    case (1):
//        tDistVector = new Vector_PETSc( aInput, aMap, aVectorType );
//        break;
    default:
        MORIS_ASSERT( false, "No vector type specified." );
        break;
    }
    return tDistVector;
}

//-------------------------------------------------------------------------------------------------
moris::Map_Class * moris::Matrix_Vector_Factory::create_map( const moris::uint        & aNumMyDofs,
                                                             const moris::Mat< int >  & aMyGlobalElements,
                                                             const moris::Mat< uint > & aMyConstraintDofs )
{
    moris::Map_Class * tMap;

    switch(0)
        {
        case (0):
            tMap = new moris::Map_Epetra ( aNumMyDofs, aMyGlobalElements, aMyConstraintDofs );
            break;
//        case (1):
//            tMap = new Map_PETSc ( aNumMyDofs, aMyGlobalElements, aMyConstraintDofs );
//            break;
        default:
            MORIS_ASSERT( false, "No map type specified" );
            break;
        }
        return tMap;
}
