/* Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details. */

#include "catch.hpp"
#include "cl_MatrixPETSc.hpp"
#include "cl_Vector_PETSc_Multi.hpp"
#include "cl_SOL_Dist_Map_Custom.hpp"
#include "cl_DLA_Preconditioner_PETSc.hpp"
#include "cl_DLA_Linear_Problem.hpp"
#include "cl_Solver_Interface_Proxy.hpp"
#include "fn_PRM_SOL_Parameters.hpp"

namespace moris
{
    namespace
    {
        // Supply a prescribed DOF without constructing a full analysis model.
        class BDDC_Interface : public Solver_Interface_Proxy
        {
          public:
            sint mConstraint = -1;

            Matrix< DDUMat > get_constrained_Ids() override
            {
                if ( mConstraint < 0 )
                {
                    return {};
                }
                return { { static_cast< uint >( mConstraint ) } };
            }
        };

        // Expose the assembled matrix to the production preconditioner builder.
        class BDDC_Problem : public dla::Linear_Problem
        {
          public:
            explicit BDDC_Problem( sol::Dist_Matrix* aMatrix ) { mMat = aMatrix; }
            void get_solution( Matrix< DDRMat >& ) override {}
            sint solve_linear_system() override { return 0; }
        };
    }    // namespace

    TEST_CASE( "PETSc MATIS assembly and BDDC", "[.][PETSc_BDDC]" )
    {
        // Each dedicated CTest run owns one PETSc lifetime.
        PetscBool tInitialized;
        REQUIRE( PetscInitialized( &tInitialized ) == 0 );
        if ( !tInitialized )
        {
            REQUIRE( PetscInitializeNoArguments() == 0 );
        }

        // Exercise both graph paths, retained constraints, and asymmetric ghost maps.
        for ( sint tCase = 0; tCase < 8; ++tCase )
        {
            bool tPreallocate = tCase % 2;
            bool tConstrained = tCase % 4 >= 2;
            bool tUnusedGhost = tCase >= 4 && par_size() > 1 && par_rank() == 0;

            // Two elements per rank; the right endpoint is owned by the next rank.
            // Noncontiguous application IDs and reversed overlap order exercise both maps.
            sint             tRank       = par_rank();
            sint             tSize       = par_size();
            sint             tOwnedCount = tRank == tSize - 1 ? 3 : 2;
            Matrix< DDSMat > tOwned( tOwnedCount, 1 );
            for ( sint i = 0; i < tOwnedCount; ++i )
            {
                tOwned( i ) = 100 + 7 * ( 2 * tRank + i );
            }

            Matrix< DDSMat > tShared( tUnusedGhost ? 4 : 3, 1 );
            for ( sint i = 0; i < 3; ++i )
            {
                tShared( i ) = 100 + 7 * ( 2 * tRank + 2 - i );
            }

            if ( tUnusedGhost )
            {
                tShared( 3 ) = 100 + 7 * ( tRank == 0 ? 2 * tSize : 0 );
            }

            Dist_Map_Custom tMap( tOwned, tShared );
            BDDC_Interface  tInput;
            PetscInt        tConstraint = 2;
            if ( tConstrained )
            {
                tInput.mConstraint = 100 + 7 * tConstraint;
            }

            // Assemble the same problem with the default AIJ and optional MATIS types.
            Matrix_PETSc tAIJ( &tInput, &tMap );
            Matrix_PETSc tIS( &tInput, &tMap, "is" );
            Mat          tA = tAIJ.get_petsc_matrix();
            Mat          tB = tIS.get_petsc_matrix();
            PetscBool    tMatches;
            REQUIRE( PetscObjectTypeCompare( (PetscObject)tB, MATIS, &tMatches ) == 0 );
            REQUIRE( tMatches );
            REQUIRE( PetscObjectTypeCompare( (PetscObject)tA, tSize == 1 ? MATSEQAIJ : MATMPIAIJ, &tMatches ) == 0 );
            REQUIRE( tMatches );

            // Cover both sparsity preallocation and element-by-element graph creation.
            Vector< moris_id > tDiagonal( tOwnedCount, tOwnedCount );
            Vector< moris_id > tOffDiagonal( tOwnedCount, tSize > 1 ? 2 : 0 );
            if ( tPreallocate )
            {
                tAIJ.build_graph( tDiagonal, tOffDiagonal );
                tIS.build_graph( tDiagonal, tOffDiagonal );
            }
            else
            {
                for ( sint e = 0; e < 2; ++e )
                {
                    Matrix< DDSMat > tIds = { { 100 + 7 * ( 2 * tRank + e ) }, { 100 + 7 * ( 2 * tRank + e + 1 ) } };
                    tAIJ.build_graph( 2, tIds );
                    tIS.build_graph( 2, tIds );
                }
                tAIJ.initial_matrix_global_assembly();
                tIS.initial_matrix_global_assembly();

                // Match the column-major convention used by the optimized graph path.
                REQUIRE( MatSetOption( tA, MAT_ROW_ORIENTED, PETSC_FALSE ) == 0 );
            }

            // Use a known solution with zero correction at the prescribed DOF.
            Vec tX, tY, tZ, tSolution;
            REQUIRE( MatCreateVecs( tA, &tX, &tY ) == 0 );
            REQUIRE( VecDuplicate( tY, &tZ ) == 0 );
            REQUIRE( VecDuplicate( tX, &tSolution ) == 0 );

            PetscInt tStart, tEnd;
            REQUIRE( VecGetOwnershipRange( tX, &tStart, &tEnd ) == 0 );
            for ( PetscInt i = tStart; i < tEnd; ++i )
            {
                REQUIRE( VecSetValue( tX, i, ( tConstrained && i == tConstraint ) ? 0.0 : 1.0 + i, INSERT_VALUES ) == 0 );
            }
            REQUIRE( VecAssemblyBegin( tX ) == 0 );
            REQUIRE( VecAssemblyEnd( tX ) == 0 );

            // Repeat assembly to catch retained values and stale preconditioner state.
            for ( sint pass = 0; pass < 3; ++pass )
            {
                tAIJ.mat_put_scalar( 0.0 );
                tIS.mat_put_scalar( 0.0 );
                for ( sint e = 0; e < 2; ++e )
                {
                    Matrix< DDSMat > tIds     = { { 100 + 7 * ( 2 * tRank + e ) }, { 100 + 7 * ( 2 * tRank + e + 1 ) } };
                    Matrix< DDRMat > tElement = { { 2.0 + pass, -1.0 }, { pass == 2 ? -0.5 : -1.0, 2.0 + pass } };
                    tAIJ.fill_matrix( 2, tElement, tIds );
                    tIS.fill_matrix( 2, tElement, tIds );
                }

                // Match MATIS constraint elimination in the AIJ reference matrix.
                if ( tConstrained )
                {
                    REQUIRE( MatSetOption( tA, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE ) == 0 );
                    if ( tStart <= tConstraint && tConstraint < tEnd )
                    {
                        REQUIRE( MatSetValue( tA, tConstraint, tConstraint, 0.0, ADD_VALUES ) == 0 );
                    }
                }

                tAIJ.matrix_global_assembly();
                tIS.matrix_global_assembly();
                if ( tConstrained )
                {
                    REQUIRE( MatZeroRowsColumns( tA, 1, &tConstraint, 1.0, nullptr, nullptr ) == 0 );
                }

                // Compare the action of both assembled matrices on the known solution.
                REQUIRE( MatMult( tA, tX, tY ) == 0 );
                REQUIRE( MatMult( tB, tX, tZ ) == 0 );
                REQUIRE( VecAXPY( tZ, -1.0, tY ) == 0 );

                PetscReal tNorm;
                REQUIRE( VecNorm( tZ, NORM_INFINITY, &tNorm ) == 0 );
                CHECK( tNorm < 1.e-12 );

                // Check the MATIS multivector path with two independent RHS columns.
                if ( pass == 0 )
                {
                    MultiVector_PETSc tInputMulti( &tInput, &tMap, 2, false );
                    MultiVector_PETSc tOutputMulti( &tInput, &tMap, 2, false );
                    for ( PetscInt c = 0; c < 2; ++c )
                    {
                        Vec tColumn;
                        REQUIRE( MatDenseGetColumnVecWrite( tInputMulti.get_petsc_vector(), c, &tColumn ) == 0 );
                        REQUIRE( VecCopy( tX, tColumn ) == 0 );
                        REQUIRE( VecScale( tColumn, c + 1.0 ) == 0 );
                        REQUIRE( MatDenseRestoreColumnVecWrite( tInputMulti.get_petsc_vector(), c, &tColumn ) == 0 );
                    }

                    tIS.mat_vec_product( tInputMulti, tOutputMulti, false );

                    for ( PetscInt c = 0; c < 2; ++c )
                    {
                        Vec tColumn;
                        REQUIRE( MatDenseGetColumnVecRead( tOutputMulti.get_petsc_vector(), c, &tColumn ) == 0 );
                        REQUIRE( VecCopy( tColumn, tZ ) == 0 );
                        REQUIRE( VecAXPY( tZ, -( c + 1.0 ), tY ) == 0 );
                        REQUIRE( VecNorm( tZ, NORM_INFINITY, &tNorm ) == 0 );
                        CHECK( tNorm < 1.e-12 );
                        REQUIRE( MatDenseRestoreColumnVecRead( tOutputMulti.get_petsc_vector(), c, &tColumn ) == 0 );
                    }
                }

                // Verify against the analytically assembled chain, not just AIJ.
                const PetscScalar* tValues;
                REQUIRE( VecGetArrayRead( tY, &tValues ) == 0 );
                for ( PetscInt i = tStart; i < tEnd; ++i )
                {
                    double tExpected = ( ( i == 0 || i == 2 * tSize ) ? 1 : 2 ) * ( 2.0 + pass ) * ( i + 1 );
                    if ( i > 0 )
                    {
                        tExpected -= ( pass == 2 ? 0.5 : 1.0 ) * i;
                    }
                    if ( i < 2 * tSize )
                    {
                        tExpected -= i + 2;
                    }
                    if ( tConstrained )
                    {
                        if ( i == tConstraint )
                        {
                            tExpected = 0.0;
                        }
                        else if ( i == tConstraint - 1 )
                        {
                            tExpected += tConstraint + 1;
                        }
                        else if ( i == tConstraint + 1 )
                        {
                            tExpected += ( pass == 2 ? 0.5 : 1.0 ) * ( tConstraint + 1 );
                        }
                    }
                    CHECK( tValues[ i - tStart ] == Approx( tExpected ) );
                }
                REQUIRE( VecRestoreArrayRead( tY, &tValues ) == 0 );

                if ( pass == 2 )
                {
                    continue;    // Nonsymmetric assembly checks column orientation.
                }

                // Solve the symmetric cases with the production BDDC setup.
                BDDC_Problem   tProblem( &tIS );
                Parameter_List tParameters = prm::create_preconditioner_parameter_list( sol::PreconditionerType::PETSC );
                tParameters.set( "PCType", "bddc" );
                tParameters.set( "use_bddc_defaults", tCase % 2 != 0 );
                tParameters.set( "bddc_use_edges", tCase % 4 != 0 );
                tParameters.set( "bddc_use_faces", false );
                tParameters.set( "bddc_use_stiffness_scaling", true );

                // Exercise deluxe scaling without symmetry/positive-definiteness
                // assumptions, as required by coupled thermoelastic operators.
                tParameters.set( "bddc_use_deluxe_scaling", true );
                tParameters.set( "bddc_sub_schurs_symmetric", false );
                tParameters.set( "bddc_sub_schurs_posdef", false );

                dla::Preconditioner_PETSc tPreconditioner( tParameters );
                KSP                       tKSP;
                REQUIRE( KSPCreate( PETSC_COMM_WORLD, &tKSP ) == 0 );
                REQUIRE( KSPSetOptionsPrefix( tKSP, "moris_bddc_parameter_test_" ) == 0 );
                REQUIRE( KSPSetOperators( tKSP, tB, tB ) == 0 );
                REQUIRE( KSPSetType( tKSP, KSPGMRES ) == 0 );
                REQUIRE( KSPSetTolerances( tKSP, 1.e-12, PETSC_DEFAULT, PETSC_DEFAULT, 100 ) == 0 );

                tPreconditioner.build_preconditioner( &tProblem, tKSP );

                // Exercise managed switches with a nonempty options prefix. The
                // default path must leave native BDDC choices untouched, and the
                // temporary settings must not contaminate subsequent solvers.
                PetscOptions tOptions;
                REQUIRE( PetscOptionsCreate( &tOptions ) == 0 );
                tPreconditioner.set_parameter_options( tOptions );

                PetscBool tUseEdges = PETSC_FALSE;
                PetscBool tFound    = PETSC_FALSE;
                REQUIRE( PetscOptionsGetBool( tOptions, "moris_bddc_parameter_test_", "-pc_bddc_use_edges", &tUseEdges, &tFound ) == 0 );
                CHECK( bool( tFound ) == ( tCase % 2 == 0 ) );
                if ( tFound )
                {
                    CHECK( bool( tUseEdges ) == ( tCase % 4 != 0 ) );
                }

                // Sub-Schur options share the PC prefix but do not contain the
                // pc_bddc_ infix. Verify that defaults leave these choices unset.
                const char* tSchurOptions[] = { "-pc_bddc_use_deluxe_scaling", "-sub_schurs_symmetric", "-sub_schurs_posdef" };
                for ( uint i = 0; i < 3; ++i )
                {
                    PetscBool tValue = PETSC_FALSE;
                    REQUIRE( PetscOptionsGetBool( tOptions, "moris_bddc_parameter_test_", tSchurOptions[ i ], &tValue, &tFound ) == 0 );
                    CHECK( bool( tFound ) == ( tCase % 2 == 0 ) );
                    if ( tFound )
                    {
                        CHECK( bool( tValue ) == ( i == 0 ) );
                    }
                }

                REQUIRE( PetscOptionsPush( tOptions ) == 0 );
                REQUIRE( KSPSetFromOptions( tKSP ) == 0 );
                REQUIRE( KSPSetUp( tKSP ) == 0 );
                REQUIRE( PetscOptionsPop() == 0 );
                REQUIRE( PetscOptionsDestroy( &tOptions ) == 0 );

                REQUIRE( PetscOptionsHasName( nullptr, "moris_bddc_parameter_test_", "-pc_bddc_use_edges", &tFound ) == 0 );
                CHECK_FALSE( tFound );
                for ( const auto* tOption : tSchurOptions )
                {
                    REQUIRE( PetscOptionsHasName( nullptr, "moris_bddc_parameter_test_", tOption, &tFound ) == 0 );
                    CHECK_FALSE( tFound );
                }

                REQUIRE( KSPSolve( tKSP, tY, tSolution ) == 0 );

                KSPConvergedReason tReason;
                REQUIRE( KSPGetConvergedReason( tKSP, &tReason ) == 0 );
                CHECK( tReason > 0 );

                // Convergence must also recover the independently specified solution.
                REQUIRE( VecAXPY( tSolution, -1.0, tX ) == 0 );
                REQUIRE( VecNorm( tSolution, NORM_INFINITY, &tNorm ) == 0 );
                CHECK( tNorm < 1.e-9 );
                REQUIRE( KSPDestroy( &tKSP ) == 0 );
            }

            REQUIRE( VecDestroy( &tX ) == 0 );
            REQUIRE( VecDestroy( &tY ) == 0 );
            REQUIRE( VecDestroy( &tZ ) == 0 );
            REQUIRE( VecDestroy( &tSolution ) == 0 );
        }

        if ( !tInitialized )
        {
            REQUIRE( PetscFinalize() == 0 );
        }
    }
}    // namespace moris
