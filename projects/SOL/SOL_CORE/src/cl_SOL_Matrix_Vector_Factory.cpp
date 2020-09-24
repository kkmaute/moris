/*
 * cl_Matrix_Vector_Factory.cpp
 *
 *  Created on: Jun 28, 2018
 *      Author: schmidt
 */
#include "cl_SOL_Matrix_Vector_Factory.hpp"

#include "cl_Sparse_Matrix_EpetraFECrs.hpp"
#include "cl_MatrixPETSc.hpp"
#include "cl_Vector_Epetra.hpp"
#include "cl_VectorPETSc.hpp"
#include "cl_Map_Epetra.hpp"
#include "cl_Map_PETSc.hpp"
#include "cl_SOL_Dist_Map.hpp"

using namespace moris;

moris::Matrix_Vector_Factory::Matrix_Vector_Factory( const enum sol::MapType aMapType )
{
    mMapType = aMapType;
}

sol::Dist_Matrix * moris::Matrix_Vector_Factory::create_matrix(       moris::Solver_Interface * aInput,
                                                             const moris::sol::Dist_Map        * aMap )
{
	sol::Dist_Matrix * tSparseMatrix = nullptr;

    switch( mMapType )
    {
    case (sol::MapType::Epetra):
        tSparseMatrix = new Sparse_Matrix_EpetraFECrs( aInput, aMap );
        break;
    case (sol::MapType::Petsc):
        tSparseMatrix = new Matrix_PETSc( aInput, aMap );
        break;
    default:
        MORIS_ERROR( false, "No matrix type specified." );
        break;
    }
    return tSparseMatrix;
}

sol::Dist_Matrix * moris::Matrix_Vector_Factory::create_matrix( const moris::uint aRows,
                                                             const moris::uint aCols )
{
	sol::Dist_Matrix * tSparseMatrix = nullptr;

    switch( mMapType )
    {
    case (sol::MapType::Epetra):
        tSparseMatrix = new Sparse_Matrix_EpetraFECrs( aRows, aCols );
        break;
    case (sol::MapType::Petsc):
        tSparseMatrix = new Matrix_PETSc( aRows, aCols );
        break;
    default:
        MORIS_ERROR( false, "No matrix type specified." );
        break;
    }
    return tSparseMatrix;
}

//-------------------------------------------------------------------------------------------------
moris::sol::Dist_Vector * moris::Matrix_Vector_Factory::create_vector(       moris::Solver_Interface * aInput,
                                                                        moris::sol::Dist_Map        * aMap,
                                                                  const sint                      aNumVectors )
{
moris::sol::Dist_Vector * tDistVector = nullptr;

    switch( mMapType )
    {
    case (sol::MapType::Epetra):
        tDistVector = new moris::Vector_Epetra( aMap, aNumVectors );
        break;
    case (sol::MapType::Petsc):
        MORIS_ERROR( aNumVectors == 1, "Multivector not implemented for petsc");
        tDistVector = new Vector_PETSc( aInput, aMap, aNumVectors );
        break;
    default:
        MORIS_ERROR( false, "No vector type specified." );
        break;
    }
    return tDistVector;
}
//-------------------------------------------------------------------------------------------------
moris::sol::Dist_Vector * moris::Matrix_Vector_Factory::create_vector(       moris::sol::Dist_Map * aMap,
                                                                       const sint                   aNumVectors )
{
moris::sol::Dist_Vector * tDistVector = nullptr;

    switch( mMapType )
    {
    case (sol::MapType::Epetra):
        tDistVector = new moris::Vector_Epetra( aMap, aNumVectors );
        break;
//    case (sol::MapType::Petsc):
//        MORIS_ERROR( aNumVectors == 1, "Multivector not implemented for petsc");
//        tDistVector = new Vector_PETSc( aInput, aMap, aNumVectors );
//        break;
    default:
        MORIS_ERROR( false, "No vector type specified." );
        break;
    }
    return tDistVector;
}


moris::sol::Dist_Vector * moris::Matrix_Vector_Factory::create_vector()
{
moris::sol::Dist_Vector * tDistVector = nullptr;

    switch( mMapType )
    {
    case (sol::MapType::Epetra):
        tDistVector = new moris::Vector_Epetra();
        break;
//    case (sol::MapType::Petsc):
//        tDistVector = new Vector_PETSc( aInput, aMap, aVectorType );
//        break;
    default:
        MORIS_ERROR( false, "No vector type specified." );
        break;
    }
    return tDistVector;
}

//-------------------------------------------------------------------------------------------------
moris::sol::Dist_Map * moris::Matrix_Vector_Factory::create_map( const moris::Matrix< DDSMat > & aMyGlobalIds,
                                                                  const moris::Matrix< DDUMat > & aMyConstraintIds )
{
    moris::sol::Dist_Map * tMap = nullptr;

    switch( mMapType )
    {
        case (sol::MapType::Epetra):
            tMap = new moris::Map_Epetra ( aMyGlobalIds, aMyConstraintIds );
            break;
        case (sol::MapType::Petsc):
            tMap = new Map_PETSc ( aMyGlobalIds, aMyConstraintIds );
            break;
        default:
            MORIS_ERROR( false, "Matrix_Vector_Factory::create_map(), map type not specified" );
            break;
    }
        return tMap;
}

//-------------------------------------------------------------------------------------------------
moris::sol::Dist_Map * moris::Matrix_Vector_Factory::create_map( const moris::Matrix< DDSMat > & aMyGlobalIds )
{
    moris::sol::Dist_Map * tMap = nullptr;

    switch( mMapType )
        {
        case (sol::MapType::Epetra):
            tMap = new moris::Map_Epetra ( aMyGlobalIds );
            break;
        case (sol::MapType::Petsc):
            tMap = new Map_PETSc ( aMyGlobalIds );
            break;
        default:
            MORIS_ERROR( false, "No map type specified" );
            break;
        }
        return tMap;
}


