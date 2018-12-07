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
#include "cl_Map_Class.hpp"

using namespace moris;

moris::Matrix_Vector_Factory::Matrix_Vector_Factory( const enum MapType aMapType )
{
    mMapType = aMapType;
}

Sparse_Matrix * moris::Matrix_Vector_Factory::create_matrix(       moris::Solver_Interface * aInput,
                                                             const moris::Map_Class        * aMap )
{
    Sparse_Matrix * tSparseMatrix;

    switch( mMapType )
    {
    case (MapType::Epetra):
        tSparseMatrix = new Sparse_Matrix_EpetraFECrs( aInput, aMap );
        break;
    case (MapType::Petsc):
        tSparseMatrix = new Matrix_PETSc( aInput, aMap );
        break;
    default:
        MORIS_ASSERT( false, "No matrix type specified." );
        break;
    }
    return tSparseMatrix;
}

//-------------------------------------------------------------------------------------------------
moris::Dist_Vector * moris::Matrix_Vector_Factory::create_vector(       moris::Solver_Interface * aInput,
                                                                  const moris::Map_Class        * aMap,
                                                                  const enum VectorType           aVectorType )
{
moris::Dist_Vector * tDistVector;

    switch( mMapType )
    {
    case (MapType::Epetra):
        tDistVector = new moris::Vector_Epetra( aMap, aVectorType );
        break;
    case (MapType::Petsc):
        tDistVector = new Vector_PETSc( aInput, aMap, aVectorType );
        break;
    default:
        MORIS_ASSERT( false, "No vector type specified." );
        break;
    }
    return tDistVector;
}

moris::Dist_Vector * moris::Matrix_Vector_Factory::create_vector()
{
moris::Dist_Vector * tDistVector;

    switch( mMapType )
    {
    case (MapType::Epetra):
        tDistVector = new moris::Vector_Epetra();
        break;
//    case (MapType::Petsc):
//        tDistVector = new Vector_PETSc( aInput, aMap, aVectorType );
//        break;
    default:
        MORIS_ASSERT( false, "No vector type specified." );
        break;
    }
    return tDistVector;
}

//-------------------------------------------------------------------------------------------------
moris::Map_Class * moris::Matrix_Vector_Factory::create_map( const moris::uint             & aNumMyDofs,
                                                             const moris::Matrix< DDSMat > & aMyGlobalElements,
                                                             const moris::Matrix< DDUMat > & aMyConstraintDofs,
                                                             const moris::Matrix< DDSMat > & aOverlappingLocaltoGlobalMap )
{
    moris::Map_Class * tMap;

    switch( mMapType )
        {
        case (MapType::Epetra):
            tMap = new moris::Map_Epetra ( aNumMyDofs, aMyGlobalElements, aMyConstraintDofs, aOverlappingLocaltoGlobalMap );
            break;
        case (MapType::Petsc):
            tMap = new Map_PETSc ( aNumMyDofs, aMyGlobalElements, aMyConstraintDofs );
            break;
        default:
            MORIS_ASSERT( false, "No map type specified" );
            break;
        }
        return tMap;
}
