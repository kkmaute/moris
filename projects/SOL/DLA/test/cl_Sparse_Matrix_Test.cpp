/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Sparse_Matrix_Test.cpp
 *
 */

#include "catch.hpp"
#include "fn_equal_to.hpp"    // ALG/src
#include "typedefs.hpp"       // COR/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_Communication_Tools.hpp"          // COM/src/
#include "cl_SOL_Matrix_Vector_Factory.hpp"    // DLA/src/
#include "cl_Solver_Interface_Proxy.hpp"       // DLA/src/
#include "cl_SOL_Dist_Vector.hpp"              // DLA/src/
#include "cl_SOL_Dist_Matrix.hpp"              // DLA/src/

#ifdef MORIS_HAVE_PETSC
#include "cl_MatrixPETSc.hpp"    // DLA/src/
#endif

namespace moris
{
    namespace sol
    {
        TEST_CASE( "Sparse Mat", "[Sparse Mat],[DistLinAlg]" )
        {
            // Determine process rank
            size_t rank = par_rank();
            size_t size = par_size();

            if ( size == 4 )
            {
                // Build Input Class
                Solver_Interface* tSolverInput = new Solver_Interface_Proxy();

                // Build matrix factory
                Matrix_Vector_Factory tMatFactory;

                // Build map
                Dist_Map* tLocalMap = tMatFactory.create_map( tSolverInput->get_my_local_global_map(),
                        tSolverInput->get_constrained_Ids() );

                // Create pointer to sparse matrix
                sol::Dist_Matrix* tMat = tMatFactory.create_matrix( tSolverInput, tLocalMap );

                // Build sparse matrix graph
                for ( moris::uint Ii = 0; Ii < tSolverInput->get_num_my_elements(); Ii++ )
                {
                    Matrix< DDSMat > tElementTopology;
                    tSolverInput->get_element_topology( Ii, tElementTopology );

                    tMat->build_graph( tElementTopology.n_rows(), tElementTopology );
                }

                // Call Global Asemby to ship information between processes
                tMat->matrix_global_assembly();

                // Fill element matrices into global matrix
                for ( uint Ii = 0; Ii < tSolverInput->get_num_my_elements(); Ii++ )
                {
                    Matrix< DDSMat > tElementTopology;
                    tSolverInput->get_element_topology( Ii, tElementTopology );

                    Matrix< DDRMat > tElementMatrix;
                    tSolverInput->get_equation_object_operator( Ii, tElementMatrix );

                    tMat->fill_matrix( tElementTopology.n_rows(), tElementMatrix, tElementTopology );
                }

                // Call Global Asemby to ship information between processes
                tMat->matrix_global_assembly();

                // tMat->print();

                // Set up output matrix
                sint tGlobalRow  = 8;
                sint tLength     = 13;
                sint tNumEntries = 5;

                moris::Matrix< DDRMat > tValues( tLength, 1, 0.0 );

                // Get matrix values
                tMat->get_matrix()->ExtractGlobalRowCopy( tGlobalRow, tLength, tNumEntries, tValues.data() );

                // Compare to true values.
                if ( rank == 0 )
                {
                    CHECK( equal_to( tValues( 0, 0 ), 24 ) );
                    CHECK( equal_to( tValues( 4, 0 ), -6 ) );
                    CHECK( equal_to( tValues( 8, 0 ), -3 ) );
                }
                delete ( tSolverInput );
                delete ( tLocalMap );
                delete ( tMat );
            }
        }

        TEST_CASE( "Scale Sparse Mat", "[Scale Sparse Mat],[DistLinAlg]" )
        {
            // Determine process rank
            size_t rank = par_rank();
            size_t size = par_size();

            if ( size == 4 )
            {
                // Build Input Class
                Solver_Interface* tSolverInput = new Solver_Interface_Proxy();

                // Build matrix factory
                Matrix_Vector_Factory tMatFactory;

                // Build map
                Dist_Map* tMap = tMatFactory.create_map( tSolverInput->get_my_local_global_map(),
                        tSolverInput->get_constrained_Ids() );

                // build distributed vector
                sol::Dist_Vector* tVectorScale = tMatFactory.create_vector( tSolverInput, tMap, 1 );

                // Create pointer to sparse matrix
                sol::Dist_Matrix* tMat = tMatFactory.create_matrix( tSolverInput, tMap );

                // Build sparse matrix graph
                for ( moris::uint Ii = 0; Ii < tSolverInput->get_num_my_elements(); Ii++ )
                {
                    Matrix< DDSMat > tElementTopology;
                    tSolverInput->get_element_topology( Ii, tElementTopology );

                    tMat->build_graph( tElementTopology.n_rows(), tElementTopology );
                }

                // Call Global Asemby to ship information between processes
                tMat->matrix_global_assembly();

                // Fill element matrices into global matrix
                for ( uint Ii = 0; Ii < tSolverInput->get_num_my_elements(); Ii++ )
                {
                    Matrix< DDSMat > tElementTopology;
                    tSolverInput->get_element_topology( Ii, tElementTopology );

                    Matrix< DDRMat > tElementMatrix;
                    tSolverInput->get_equation_object_operator( Ii, tElementMatrix );

                    tMat->fill_matrix( tElementTopology.n_rows(), tElementMatrix, tElementTopology );
                }

                // Call Global Asemby to ship information between processes
                tMat->matrix_global_assembly();

                // tMat->print();

                // Put Values in scaling vector
                tVectorScale->vec_put_scalar( 2.25 );

                tMat->sparse_mat_left_scale( *tVectorScale );

                // Set up output matrix
                moris::sint             tGlobalRow  = 8;
                moris::sint             tLength     = 13;
                moris::sint             tNumEntries = 5;
                moris::Matrix< DDRMat > tValues( tLength, 1, 0.0 );

                // Get matrix values
                tMat->get_matrix()->ExtractGlobalRowCopy( tGlobalRow, tLength, tNumEntries, tValues.data() );

                // Compare to true values.
                if ( rank == 0 )
                {
                    CHECK( equal_to( tValues( 0, 0 ), 54 ) );
                    CHECK( equal_to( tValues( 4, 0 ), -13.5 ) );
                    CHECK( equal_to( tValues( 8, 0 ), -6.75 ) );
                }
                delete ( tSolverInput );
                delete ( tMap );
                delete ( tVectorScale );
                delete ( tMat );
            }
        }

        TEST_CASE( "Diagonal Sparse Mat", "[Diagonal Sparse Mat],[DistLinAlg]" )
        {
            // Determine process rank
            size_t rank = par_rank();
            size_t size = par_size();

            if ( size == 4 )
            {
                // Build Input Class
                Solver_Interface* tSolverInput = new Solver_Interface_Proxy();

                // Build matrix factory
                Matrix_Vector_Factory tMatFactory;

                // Build map
                Dist_Map* tMap = tMatFactory.create_map( tSolverInput->get_my_local_global_map(),
                        tSolverInput->get_constrained_Ids() );

                // build distributed vector
                sol::Dist_Vector* tVectorDiagonal = tMatFactory.create_vector( tSolverInput, tMap, 1 );

                // Create pointer to sparse matrix
                sol::Dist_Matrix* tMat = tMatFactory.create_matrix( tSolverInput, tMap );

                // Build sparse matrix graph
                for ( moris::uint Ii = 0; Ii < tSolverInput->get_num_my_elements(); Ii++ )
                {
                    Matrix< DDSMat > tElementTopology;
                    tSolverInput->get_element_topology( Ii, tElementTopology );

                    tMat->build_graph( tElementTopology.n_rows(), tElementTopology );
                }

                // Call Global Asemby to ship information between processes
                tMat->matrix_global_assembly();

                // Fill element matrices into global matrix
                for ( uint Ii = 0; Ii < tSolverInput->get_num_my_elements(); Ii++ )
                {
                    Matrix< DDSMat > tElementTopology;
                    tSolverInput->get_element_topology( Ii, tElementTopology );

                    Matrix< DDRMat > tElementMatrix;
                    tSolverInput->get_equation_object_operator( Ii, tElementMatrix );

                    tMat->fill_matrix( tElementTopology.n_rows(), tElementMatrix, tElementTopology );
                }

                // Call Global Asemby to ship information between processes
                tMat->matrix_global_assembly();

                // tMat->print();

                // extract diagonal value into Dist Vector
                tMat->get_diagonal( *tVectorDiagonal );

                moris::Matrix< DDRMat > tDiagonal( 15, 1, 0.0 );

                // needed as offset parameter for Epetra. =0
                sint tMyLDA = 0;

                // Get solution and output it in moris::Mat LHSValues
                dynamic_cast< Vector_Epetra* >( tVectorDiagonal )->get_epetra_vector()->ExtractCopy( tDiagonal.data(), tMyLDA );

                // Compare to true values.
                if ( rank == 0 )
                {
                    CHECK( equal_to( tDiagonal( 0, 0 ), 24 ) );
                    CHECK( equal_to( tDiagonal( 1, 0 ), 24 ) );
                    CHECK( equal_to( tDiagonal( 2, 0 ), 48 ) );
                    CHECK( equal_to( tDiagonal( 3, 0 ), 48 ) );
                    CHECK( equal_to( tDiagonal( 4, 0 ), 24 ) );
                    CHECK( equal_to( tDiagonal( 5, 0 ), 24 ) );
                }
                if ( rank == 3 )
                {
                    CHECK( equal_to( tDiagonal( 0, 0 ), 24 ) );
                    CHECK( equal_to( tDiagonal( 1, 0 ), 24 ) );
                    CHECK( equal_to( tDiagonal( 2, 0 ), 12 ) );
                    CHECK( equal_to( tDiagonal( 3, 0 ), 12 ) );
                }

                // Set all entries of vector to 100.33
                tVectorDiagonal->vec_put_scalar( 100.33 );

                // Replace Matrix diagonal values with tVectorDiagonal values
                tMat->replace_diagonal_values( *tVectorDiagonal );

                // Set up output matrix
                moris::sint tGlobalRow  = 8;
                moris::sint tLength     = 13;
                moris::sint tNumEntries = 5;

                moris::Matrix< DDRMat > tValues( tLength, 1, 0.0 );

                // Get matrix values
                tMat->get_matrix()->ExtractGlobalRowCopy( tGlobalRow, tLength, tNumEntries, tValues.data() );

                // Compare to true values.
                if ( rank == 0 )
                {
                    CHECK( equal_to( tValues( 0, 0 ), 100.33 ) );
                    CHECK( equal_to( tValues( 4, 0 ), -6 ) );
                    CHECK( equal_to( tValues( 8, 0 ), -3 ) );
                }
                delete ( tSolverInput );
                delete ( tMap );
                delete ( tVectorDiagonal );
                delete ( tMat );
            }
        }

#ifdef MORIS_HAVE_PETSC
        TEST_CASE( "Get Matrix Values", "[Get_Matrix_Values],[DistLinAlg]" )
        {
            // Determine process rank
            size_t size = par_size();

            if ( size == 1 )
            {
                PetscInitializeNoArguments();
                // Create pointer to sparse matrix
                sol::Dist_Matrix* tMat = new moris::Matrix_PETSc( 4, 4 );

                Matrix< DDSMat > tElementIds = { { 0, 1, 2, 3 } };

                // Fill element matrices into global matrix
                Matrix< DDRMat > tElementMatrix = { { 0, 1, 2, 3 },
                    { 4, 5, 6, 7 },
                    { 8, 9, 10, 11 },
                    { 12, 13, 14, 15 } };

                tMat->insert_values( tElementIds, tElementIds, tElementMatrix );

                // Call Global Asemby to ship information between processes
                tMat->matrix_global_assembly();

                // Get matrix values
                Matrix< DDSMat > tRequestedIds = { { 1, 2, 3 } };
                Matrix< DDRMat > tValues( 3, 3, 0.0 );

                tMat->get_matrix_values( tRequestedIds, tValues );

                CHECK( equal_to( tValues( 0, 0 ), 5 ) );
                CHECK( equal_to( tValues( 0, 1 ), 6 ) );
                CHECK( equal_to( tValues( 0, 2 ), 7 ) );
                CHECK( equal_to( tValues( 1, 0 ), 9 ) );
                CHECK( equal_to( tValues( 1, 1 ), 10 ) );
                CHECK( equal_to( tValues( 1, 2 ), 11 ) );
                CHECK( equal_to( tValues( 2, 0 ), 13 ) );
                CHECK( equal_to( tValues( 2, 1 ), 14 ) );
                CHECK( equal_to( tValues( 2, 2 ), 15 ) );

                delete ( tMat );

                PetscFinalize();
            }
        }
#endif

        TEST_CASE( "Non-square matrix", "[Non-square matrix],[DistLinAlg]" )
        {
            // Determine process rank
            size_t rank = par_rank();
            size_t size = par_size();

            if ( size == 2 )
            {
                // Build matrix factory
                sol::Matrix_Vector_Factory tMatFactory;

                moris::Matrix< DDSMat > tRowMapVal;
                moris::Matrix< DDSMat > tColMapVal;
                if ( rank == 0 )
                {
                    tRowMapVal = { { 0 } };
                    tColMapVal = { { 0 }, { 1 } };
                }
                else if ( rank == 1 )
                {
                    tRowMapVal = { { 1 } };
                    tColMapVal = { { 2 }, { 3 } };
                }

                // Build map
                sol::Dist_Map* tRowMap = tMatFactory.create_map( tRowMapVal );
                sol::Dist_Map* tColMap = tMatFactory.create_map( tColMapVal );

                // build distributed vector
                sol::Dist_Vector* tVector1 = tMatFactory.create_vector( tColMap );
                sol::Dist_Vector* tVector2 = tMatFactory.create_vector( tRowMap );

                tVector1->vec_put_scalar( 1.0 );

                // Create pointer to sparse matrix
                sol::Dist_Matrix* tMat = tMatFactory.create_matrix( tRowMap, tColMap );

                moris::Matrix< DDSMat > tRows;
                moris::Matrix< DDSMat > tCols;
                moris::Matrix< DDRMat > tVals;

                if ( rank == 0 )
                {
                    tRows = { { 0 }, { 1 } };
                    tCols = { { 0 }, { 1 }, { 2 }, { 3 } };
                    tVals = { { 1, 1, 0, 0 }, { 0, 0, 1, 1 } };

                    tMat->insert_values( tRows, tCols, tVals );
                }

                // Call Global Asemby to ship information between processes
                tMat->matrix_global_assembly();

                // tMat->print();

                tMat->mat_vec_product( *tVector1, *tVector2, false );

                // tVector2->print();

                moris::Matrix< DDRMat > tResult( 1, 1, 0.0 );

                // needed as offset parameter for Epetra. =0
                sint tMyLDA = 0;

                // Get solution and output it in moris::Mat LHSValues
                static_cast< Vector_Epetra* >( tVector2 )->get_epetra_vector()->ExtractCopy( tResult.data(), tMyLDA );

                // Compare to true values.
                if ( rank == 0 )
                {
                    CHECK( equal_to( tResult( 0, 0 ), 2.0 ) );
                }
                if ( rank == 1 )
                {
                    CHECK( equal_to( tResult( 0, 0 ), 2.0 ) );
                }

                delete ( tVector1 );
                delete ( tVector2 );
                delete ( tMat );
                delete tRowMap;
                delete tColMap;
            }
        }

        TEST_CASE( "Non-square matrix sum into values", "[Non-square matrix],[DistLinAlg],[sum into values]" )
        {
            // Determine process rank
            size_t rank = par_rank();
            size_t size = par_size();

            if ( size == 2 )
            {
                // Build matrix factory
                sol::Matrix_Vector_Factory tMatFactory;

                moris::Matrix< DDSMat > tRowMapVal;
                moris::Matrix< DDSMat > tColMapVal;
                if ( rank == 0 )
                {
                    tRowMapVal = { { 0 } };
                    tColMapVal = { { 0 }, { 1 } };
                }
                else if ( rank == 1 )
                {
                    tRowMapVal = { { 1 } };
                    tColMapVal = { { 2 }, { 3 } };
                }

                // Build map
                sol::Dist_Map* tRowMap = tMatFactory.create_map( tRowMapVal );
                sol::Dist_Map* tColMap = tMatFactory.create_map( tColMapVal );

                // build distributed vector
                sol::Dist_Vector* tVector1 = tMatFactory.create_vector( tColMap );
                sol::Dist_Vector* tVector2 = tMatFactory.create_vector( tRowMap );

                tVector1->vec_put_scalar( 1.0 );

                // Create pointer to sparse matrix
                sol::Dist_Matrix* tMat = tMatFactory.create_matrix( tRowMap, tColMap );

                moris::Matrix< DDSMat > tRows;
                moris::Matrix< DDSMat > tCols;
                moris::Matrix< DDRMat > tVals;

                if ( rank == 0 )
                {
                    tRows = { { 0 }, { 1 } };
                    tCols = { { 0 }, { 1 }, { 2 }, { 3 } };
                    tVals = { { 0, 0, 0, 0 }, { 0, 0, 0, 0 } };

                    tMat->insert_values( tRows, tCols, tVals );
                }

                // Call Global Asemby to ship information between processes
                tMat->matrix_global_assembly();

                if ( rank == 0 )
                {
                    tVals = { { 1, 1, 0, 0 }, { 0, 0, 1, 1 } };
                    tMat->sum_into_values( tRows, tCols, tVals );
                }

                tMat->matrix_global_assembly();

                tMat->print();

                tMat->mat_vec_product( *tVector1, *tVector2, false );

                tVector2->print();

                moris::Matrix< DDRMat > tResult( 1, 1, 0.0 );

                // needed as offset parameter for Epetra. =0
                sint tMyLDA = 0;

                // Get solution and output it in moris::Mat LHSValues
                static_cast< Vector_Epetra* >( tVector2 )->get_epetra_vector()->ExtractCopy( tResult.data(), tMyLDA );

                // Compare to true values.
                if ( rank == 0 )
                {
                    CHECK( equal_to( tResult( 0, 0 ), 2.0 ) );
                }
                if ( rank == 1 )
                {
                    CHECK( equal_to( tResult( 0, 0 ), 2.0 ) );
                }

                delete ( tVector1 );
                delete ( tVector2 );
                delete ( tMat );
                delete tRowMap;
                delete tColMap;
            }
        }
    }    // namespace sol
}    // namespace moris
