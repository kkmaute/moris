/* Copyright (c) 2026 University of Colorado. Licensed under the MIT license. */

#include "catch.hpp"

#ifdef MORIS_HAVE_PETSC
#include <memory>
#include <utility>
#include <petscsys.h>

#include "cl_Communication_Tools.hpp"
#include "cl_SOL_Dist_Map.hpp"
#include "cl_SOL_Dist_Vector.hpp"
#include "cl_SOL_Matrix_Vector_Factory.hpp"
#include "cl_Solver_Interface_Proxy.hpp"
#include "cl_Vector_PETSc_Multi.hpp"

namespace moris::sol
{
    namespace
    {
        class Carry_Input : public Solver_Interface_Proxy
        {
          public:
            Matrix< DDUMat > get_constrained_Ids() override
            {
                return {};
            }
        };

        // This hidden test runs in its own CTest process. Other DLA tests may
        // initialize/finalize PETSc, which cannot be initialized a second time.
        struct Petsc_Session
        {
            Petsc_Session()
            {
                REQUIRE( PetscInitializeNoArguments() == 0 );
            }

            ~Petsc_Session()
            {
                PetscFinalize();
            }
        };

        sint dof_id( sint aIndex )
        {
            return 17 + 13 * aIndex;
        }

        real value( sint aId, uint aColumn )
        {
            return 0.01 * aId + 7.0 * ( aColumn + 1 );
        }

        // Fill the raw dense columns so the reference does not depend on the
        // transfer or map-translation routines under test. PETSc row numbering
        // follows local ownership order, whereas expected values use MORIS IDs.
        void fill( Dist_Vector& aVector, const Matrix< DDSMat >& aIds, uint aColumns, real aScale )
        {
            Mat      tMatrix = dynamic_cast< MultiVector_PETSc& >( aVector ).get_petsc_vector();
            PetscInt tFirst, tEnd;
            REQUIRE( MatGetOwnershipRange( tMatrix, &tFirst, &tEnd ) == 0 );
            REQUIRE( tEnd - tFirst == (PetscInt)aIds.numel() );

            for ( uint c = 0; c < aColumns; ++c )
            {
                for ( uint i = 0; i < aIds.numel(); ++i )
                {
                    const PetscInt    tRow    = tFirst + i;
                    const PetscInt    tColumn = c;
                    const PetscScalar tValue  = aScale * value( aIds( i ), c );
                    REQUIRE( MatSetValues( tMatrix, 1, &tRow, 1, &tColumn, &tValue, INSERT_VALUES ) == 0 );
                }
            }
            aVector.vector_global_assembly();
        }

        void check( Dist_Vector& aVector, const Matrix< DDSMat >& aIds, uint aColumns, real aScale )
        {
            Matrix< DDRMat > tValues;
            aVector.extract_copy( tValues );
            REQUIRE( tValues.n_rows() == aIds.numel() );
            REQUIRE( tValues.n_cols() == aColumns );
            for ( uint c = 0; c < aColumns; ++c )
            {
                for ( uint i = 0; i < aIds.numel(); ++i )
                {
                    CHECK( tValues( i, c ) == Approx( aScale * value( aIds( i ), c ) ).margin( 1e-12 ) );
                }
            }
        }

        void check_transfer( const Matrix< DDSMat >& aSourceOwned,
                const Matrix< DDSMat >&              aSourceShared,
                const Matrix< DDSMat >&              aTargetOwned,
                const Matrix< DDSMat >&              aTargetShared,
                bool                                 aFullSource,
                bool                                 aFullTarget,
                uint                                 aColumns )
        {
            Carry_Input                    tInput;
            Matrix_Vector_Factory          tFactory( MapType::Petsc );
            Matrix< DDUMat >               tConstraints;
            std::unique_ptr< Dist_Map >    tSourceMap( aFullSource
                                                            ? tFactory.create_full_map( aSourceOwned, aSourceShared )
                                                            : tFactory.create_map( aSourceOwned, tConstraints ) );
            std::unique_ptr< Dist_Map >    tTargetMap( aFullTarget
                                                            ? tFactory.create_full_map( aTargetOwned, aTargetShared )
                                                            : tFactory.create_map( aTargetOwned, tConstraints ) );
            std::unique_ptr< Dist_Vector > tSource( tFactory.create_vector( &tInput, tSourceMap.get(), aColumns ) );
            std::unique_ptr< Dist_Vector > tTarget( tFactory.create_vector( &tInput, tTargetMap.get(), aColumns ) );
            std::unique_ptr< Dist_Vector > tWorkspace( tFactory.create_vector( &tInput, tTargetMap.get(), aColumns ) );
            const auto&                    tSourceIds = aFullSource ? aSourceShared : aSourceOwned;
            const auto&                    tTargetIds = aFullTarget ? aTargetShared : aTargetOwned;
            fill( *tSource, tSourceIds, aColumns, 1.0 );
            fill( *tTarget, tTargetIds, aColumns, 0.5 );
            fill( *tWorkspace, tTargetIds, aColumns, 99.0 );

            // Repeated residual evaluations add the immutable carry once each;
            // stale workspace values must never enter the result.
            for ( uint i = 0; i < 2; ++i )
            {
                tTarget->add_vector( *tSource, *tWorkspace );
                check( *tTarget, tTargetIds, aColumns, 1.5 + i );
                check( *tSource, tSourceIds, aColumns, 1.0 );
            }
            fill( *tSource, tSourceIds, aColumns, 3.0 );
            tTarget->add_vector( *tSource, *tWorkspace );
            check( *tTarget, tTargetIds, aColumns, 5.5 );
            check( *tSource, tSourceIds, aColumns, 3.0 );
        }
    }    // namespace

    // Full-state carries are restricted to reordered subsystem maps. Exercise
    // off-rank ownership, local ghost copies, all objective columns, repeated
    // additions and empty ranks against explicit ID-based numerical references.
    TEST_CASE( "PETSc adjoint carry follows MORIS IDs", "[.][PETScAdjointCarry]" )
    {
        Petsc_Session    tSession;
        const sint       tRank = par_rank();
        const sint       tSize = par_size();
        Matrix< DDSMat > tOwned( 4, 1 ), tShared( 4 * tSize, 1 );
        for ( sint i = 0; i < 4; ++i )
        {
            tOwned( i ) = dof_id( 4 * tRank + i );
        }
        for ( sint i = 0; i < 4 * tSize; ++i )
        {
            tShared( i ) = dof_id( 4 * tSize - i - 1 );
        }

        for ( bool tSubset : { false, true } )
        {
            const sint       tRows = tSubset ? 2 : 4;
            Matrix< DDSMat > tTargetOwned( tRows, 1 ), tTargetShared( tRows * tSize, 1 );
            for ( sint i = 0; i < tRows; ++i )
            {
                tTargetOwned( i ) = dof_id( 4 * ( ( tRank + 1 ) % tSize ) + 3 - ( tSubset ? 2 : 1 ) * i );
            }
            for ( sint r = 0; r < tSize; ++r )
            {
                for ( sint i = 0; i < tRows; ++i )
                {
                    tTargetShared( tRows * r + i ) = dof_id( 4 * r + 3 - ( tSubset ? 2 : 1 ) * i );
                }
            }
            for ( bool tFullSource : { false, true } )
            {
                for ( bool tFullTarget : { false, true } )
                {
                    for ( uint tColumns : { 1u, 3u } )
                    {
                        INFO( "subset=" << tSubset << ", full source=" << tFullSource
                                        << ", full target=" << tFullTarget << ", RHS=" << tColumns );
                        check_transfer( tOwned, tShared, tTargetOwned, tTargetShared, tFullSource, tFullTarget, tColumns );
                    }
                }
            }
        }

        // Only rank zero owns destination rows; all ranks must still enter the
        // scatter. At 2/4 ranks the requested source entries live on another rank.
        Matrix< DDSMat > tEmptyTarget( tRank == 0 ? 2 : 0, 1 );
        if ( tRank == 0 )
        {
            tEmptyTarget( 0 ) = dof_id( 4 * ( tSize - 1 ) + 3 );
            tEmptyTarget( 1 ) = dof_id( 4 * ( tSize - 1 ) + 1 );
        }
        check_transfer( tOwned, tShared, tEmptyTarget, tEmptyTarget, false, false, 3 );

        // Empty source ranks still participate when destination ownership moves
        // away from rank zero. This also covers an empty source/destination pair.
        Matrix< DDSMat > tOnlySource( tRank == 0 ? 4 : 0, 1 );
        Matrix< DDSMat > tMovedTarget( tRank == tSize - 1 ? 2 : 0, 1 );
        if ( tRank == 0 )
        {
            for ( sint i = 0; i < 4; ++i )
            {
                tOnlySource( i ) = dof_id( i );
            }
        }
        if ( tRank == tSize - 1 )
        {
            tMovedTarget( 0 ) = dof_id( 3 );
            tMovedTarget( 1 ) = dof_id( 1 );
        }
        check_transfer( tOnlySource, tOnlySource, tMovedTarget, tMovedTarget, false, false, 3 );

        // Missing IDs on one rank must fail collectively before changing RHS.
        // A second failure verifies that equal sizes cannot hide a reordered
        // workspace and that a bad transfer does not poison the next valid one.
        Carry_Input           tInput;
        Matrix_Vector_Factory tFactory( MapType::Petsc );
        Matrix< DDUMat >      tConstraints;
        auto                  tBadIds = tOwned;
        if ( tRank == 0 )
        {
            tBadIds( 0 ) = 999999;
        }
        std::unique_ptr< Dist_Map >    tSourceMap( tFactory.create_map( tOwned, tConstraints ) );
        std::unique_ptr< Dist_Map >    tBadMap( tFactory.create_map( tBadIds, tConstraints ) );
        std::unique_ptr< Dist_Vector > tSource( tFactory.create_vector( &tInput, tSourceMap.get(), 3 ) );
        std::unique_ptr< Dist_Vector > tTarget( tFactory.create_vector( &tInput, tBadMap.get(), 3 ) );
        std::unique_ptr< Dist_Vector > tWorkspace( tFactory.create_vector( &tInput, tBadMap.get(), 3 ) );
        std::unique_ptr< Dist_Vector > tWrongColumns( tFactory.create_vector( &tInput, tBadMap.get(), 1 ) );
        fill( *tSource, tOwned, 3, 1.0 );
        fill( *tTarget, tBadIds, 3, 0.5 );
        CHECK_THROWS( tTarget->add_vector( *tSource, *tWorkspace ) );
        CHECK_THROWS( tTarget->add_vector( *tWrongColumns, *tWorkspace ) );
        CHECK_THROWS( tTarget->add_vector( *tSource, *tSource ) );
        check( *tTarget, tBadIds, 3, 0.5 );
        check( *tSource, tOwned, 3, 1.0 );

        auto tReorderedIds = tOwned;
        std::swap( tReorderedIds( 0 ), tReorderedIds( 1 ) );
        std::unique_ptr< Dist_Map >    tReorderedMap( tFactory.create_map( tReorderedIds, tConstraints ) );
        std::unique_ptr< Dist_Vector > tReordered( tFactory.create_vector( &tInput, tReorderedMap.get(), 3 ) );
        fill( *tReordered, tReorderedIds, 3, 0.5 );
        CHECK_THROWS( tSource->vec_plus_vec( 1.0, *tReordered, 1.0 ) );
        CHECK_THROWS( tReordered->add_vector( *tSource, *tWorkspace ) );
        CHECK_THROWS( tReordered->add_vector( *tSource, *tReordered ) );
        std::unique_ptr< Dist_Vector > tValidWorkspace( tFactory.create_vector( &tInput, tReorderedMap.get(), 3 ) );
        tReordered->add_vector( *tSource, *tValidWorkspace );
        check( *tReordered, tReorderedIds, 3, 1.5 );
        check( *tSource, tOwned, 3, 1.0 );
    }
}    // namespace moris::sol
#endif
