/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SOL_Matrix_Vector_Factory.cpp
 *
 */

#include "cl_SOL_Matrix_Vector_Factory.hpp"

#include "cl_Sparse_Matrix_EpetraFECrs.hpp"
#include "cl_Vector_Epetra.hpp"
#include "cl_Map_Epetra.hpp"
#include "cl_SOL_Dist_Map.hpp"

#ifdef MORIS_HAVE_PETSC
#include "cl_MatrixPETSc.hpp"
#include "cl_Vector_PETSc.hpp"
#include "cl_Vector_PETSc_Multi.hpp"
#include "cl_Map_PETSc.hpp"
#endif

namespace moris
{
    namespace sol
    {
        //-------------------------------------------------------------------------------------------

        sol::Matrix_Vector_Factory::Matrix_Vector_Factory( const enum MapType aMapType )
        {
            mMapType = aMapType;
        }

        //-------------------------------------------------------------------------------------------

        Dist_Matrix*
        sol::Matrix_Vector_Factory::create_matrix(
                Solver_Interface* aInput,
                Dist_Map*         aMap,
                bool              aPointMap,
                bool              aBuildGraph )
        {
            Dist_Matrix* tSparseMatrix = nullptr;

            switch ( mMapType )
            {
                case MapType::Epetra:
                {
                    tSparseMatrix = new Sparse_Matrix_EpetraFECrs( aInput, aMap, aPointMap, aBuildGraph );
                    break;
                }
                case MapType::Petsc:
                {
#ifdef MORIS_HAVE_PETSC
                    tSparseMatrix = new Matrix_PETSc( aInput, aMap );
#else
                    MORIS_ERROR( false, "MORIS is configured with out PETSC support." );
#endif
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "No matrix type specified." );
                    break;
                }
            }
            return tSparseMatrix;
        }

        //-------------------------------------------------------------------------------------------

        Dist_Matrix*
        sol::Matrix_Vector_Factory::create_matrix(
                Dist_Map* aRowMap,
                Dist_Map* aColMap )
        {
            Dist_Matrix* tSparseMatrix = nullptr;

            switch ( mMapType )
            {
                case MapType::Epetra:
                {
                    tSparseMatrix = new Sparse_Matrix_EpetraFECrs( aRowMap, aColMap );
                    break;
                }
                // case (MapType::Petsc):
                //     tSparseMatrix = new Matrix_PETSc( aInput, aMap );
                //     break;
                default:
                {
                    MORIS_ERROR( false, "No matrix type specified." );
                    break;
                }
            }
            return tSparseMatrix;
        }

        //-------------------------------------------------------------------------------------------

        Dist_Matrix*
        sol::Matrix_Vector_Factory::create_matrix(
                const uint aRows,
                const uint aCols )
        {
            Dist_Matrix* tSparseMatrix = nullptr;

            switch ( mMapType )
            {
                case MapType::Epetra:
                {
                    tSparseMatrix = new Sparse_Matrix_EpetraFECrs( aRows, aCols );
                    break;
                }
                case MapType::Petsc:
                {
#ifdef MORIS_HAVE_PETSC
                    tSparseMatrix = new Matrix_PETSc( aRows, aCols );
#else
                    MORIS_ERROR( false, "MORIS is configured with out PETSC support." );
#endif
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "No matrix type specified." );
                    break;
                }
            }
            return tSparseMatrix;
        }

        //-------------------------------------------------------------------------------------------------

        Dist_Vector*
        sol::Matrix_Vector_Factory::create_vector(
                Solver_Interface* aInput,
                Dist_Map*         aMap,
                const sint        aNumVectors,
                bool              aPointMap,
                bool              aManageMap )
        {
            Dist_Vector* tDistVector = nullptr;

            switch ( mMapType )
            {
                case MapType::Epetra:
                {
                    tDistVector = new Vector_Epetra( aMap, aNumVectors, aPointMap, aManageMap );
                    break;
                }
                case MapType::Petsc:
                {
#ifdef MORIS_HAVE_PETSC
                    tDistVector = new MultiVector_PETSc( aInput, aMap, aNumVectors, aManageMap );
#else
                    MORIS_ERROR( false, "MORIS is configured with out PETSC support." );
#endif
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "No vector type specified." );
                    break;
                }
            }
            return tDistVector;
        }

        //-------------------------------------------------------------------------------------------------

        Dist_Vector*
        sol::Matrix_Vector_Factory::create_vector(
                Dist_Map*  aMap,
                const sint aNumVectors,
                bool       aPointMap,
                bool       aManageMap )
        {
            Dist_Vector* tDistVector = nullptr;

            switch ( mMapType )
            {
                case MapType::Epetra:
                {
                    tDistVector = new Vector_Epetra( aMap, aNumVectors, aPointMap, aManageMap );
                    break;
                }
                //    case (MapType::Petsc):
                //        MORIS_ERROR( aNumVectors == 1, "Multivector not implemented for petsc");
                //        tDistVector = new Vector_PETSc( aInput, aMap, aNumVectors );
                //        break;
                default:
                {
                    MORIS_ERROR( false, "No vector type specified." );
                    break;
                }
            }
            return tDistVector;
        }

        //-------------------------------------------------------------------------------------------------
        Dist_Map*
        sol::Matrix_Vector_Factory::create_map(
                const Matrix< DDSMat >& aMyGlobalIds,
                const Matrix< DDUMat >& aMyConstraintIds )
        {
            Dist_Map* tMap = nullptr;

            switch ( mMapType )
            {
                case MapType::Epetra:
                {
                    tMap = new Map_Epetra( aMyGlobalIds, aMyConstraintIds );
                    break;
                }
                case MapType::Petsc:
                {
#ifdef MORIS_HAVE_PETSC
                    tMap = new Map_PETSc( aMyGlobalIds, aMyConstraintIds );
#else
                    MORIS_ERROR( false, "MORIS is configured with out PETSC support." );
#endif
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "Matrix_Vector_Factory::create_map(), map type not specified" );
                    break;
                }
            }
            return tMap;
        }

        //-------------------------------------------------------------------------------------------------

        Dist_Map*
        sol::Matrix_Vector_Factory::create_map( const Matrix< DDSMat >& aMyGlobalIds )
        {
            Dist_Map* tMap = nullptr;

            switch ( mMapType )
            {
                case MapType::Epetra:
                {
                    tMap = new Map_Epetra( aMyGlobalIds );
                    break;
                }
                case MapType::Petsc:
                {
#ifdef MORIS_HAVE_PETSC
                    tMap = new Map_PETSc( aMyGlobalIds );
#else
                    MORIS_ERROR( false, "MORIS is configured with out PETSC support." );
#endif
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "No map type specified" );
                    break;
                }
            }
            return tMap;
        }

        //-------------------------------------------------------------------------------------------------

        Dist_Map*
        sol::Matrix_Vector_Factory::create_map( const Vector< sint >& aMyGlobalIds )
        {
            Dist_Map* tMap = nullptr;

            switch ( mMapType )
            {
                case MapType::Epetra:
                {
                    tMap = new Map_Epetra( aMyGlobalIds );
                    break;
                }
                case MapType::Petsc:
                {
#ifdef MORIS_HAVE_PETSC
                    tMap = new Map_PETSc( aMyGlobalIds );
#else
                    MORIS_ERROR( false, "MORIS is configured with out PETSC support." );
#endif
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "No map type specified" );
                    break;
                }
            }
            return tMap;
        }

        //-------------------------------------------------------------------------------------------------

        Dist_Map*
        sol::Matrix_Vector_Factory::create_full_map(
                const moris::Matrix< DDSMat >& aMyGlobalOwnedIds,
                const moris::Matrix< DDSMat >& aMyGlobalOwnedAndSharedIds )
        {
            Dist_Map* tMap = nullptr;

            switch ( mMapType )
            {
                case MapType::Epetra:
                {
                    tMap = new Map_Epetra( aMyGlobalOwnedAndSharedIds );
                    break;
                }
                case MapType::Petsc:
                {
#ifdef MORIS_HAVE_PETSC
                    tMap = new Map_PETSc(
                            aMyGlobalOwnedIds,
                            aMyGlobalOwnedAndSharedIds,
                            true );
#else
                    MORIS_ERROR( false, "MORIS is configured with out PETSC support." );
#endif
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "No map type specified" );
                    break;
                }
            }
            return tMap;
        }
    }    // namespace sol
}    // namespace moris
