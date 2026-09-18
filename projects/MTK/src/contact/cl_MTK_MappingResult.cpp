/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_MTK_MappingResult.cpp
 *
 */
#include "cl_MTK_MappingResult.hpp"
#include "cl_Json_Object.hpp"

#include "fn_assert.hpp"

// Assert parametric coordinate in [-1,1] for lines/quads
inline void assert_param_in_bounds_box( const moris::Matrix< moris::DDRMat >& aParam, const char* aContext )
{
    if ( aParam.n_rows() >= 1 )
    {
        MORIS_ASSERT(
                aParam( 0 ) >= -1.0 - 1e-12 && aParam( 0 ) <= 1.0 + 1e-12,
                "MTK parametric coordinate out of bounds (param[0]=%e, context=%s)",
                aParam( 0 ),
                aContext );
    }
    if ( aParam.n_rows() >= 2 )
    {
        MORIS_ASSERT(
                aParam( 1 ) >= -1.0 - 1e-12 && aParam( 1 ) <= 1.0 + 1e-12,
                "MTK parametric coordinate out of bounds (param[1]=%e, context=%s)",
                aParam( 1 ),
                aContext );
    }
}

// Assert parametric coordinate in [0,1] for triangles/simplex
inline void assert_param_in_bounds_simplex( const moris::Matrix< moris::DDRMat >& aParam, const char* aContext )
{
    if ( aParam.n_rows() >= 1 )
    {
        MORIS_ASSERT(
                aParam( 0 ) >= 0.0 - 1e-12 && aParam( 0 ) <= 1.0 + 1e-12,
                "MTK TRI parametric coordinate out of bounds (eta=%e, context=%s)",
                aParam( 0 ),
                aContext );
    }
    if ( aParam.n_rows() >= 2 )
    {
        MORIS_ASSERT(
                aParam( 1 ) >= 0.0 - 1e-12 && aParam( 1 ) <= 1.0 + 1e-12,
                "MTK TRI parametric coordinate out of bounds (zeta=%e, context=%s)",
                aParam( 1 ),
                aContext );
        MORIS_ASSERT(
                aParam( 0 ) + aParam( 1 ) <= 1.0 + 1e-12,
                "MTK TRI parametric coordinate out of bounds (eta+zeta=%e, context=%s)",
                aParam( 0 ) + aParam( 1 ),
                aContext );
    }
}

namespace moris::mtk
{
    MappingResult::MappingResult(
            moris_index aSourceMeshIndex,
            uint        aPhysicalDimension,
            uint        aNumberOfPoints )
            : mSourceMeshIndex( aSourceMeshIndex )
            , mSourcePhysicalCoordinate( aPhysicalDimension, aNumberOfPoints )
            , mSourceCellIndex( aNumberOfPoints, -1 )
            , mSourceClusterIndex( aNumberOfPoints, -1 )
            , mTargetClusterIndex( aNumberOfPoints, -1 )
            , mTargetCellIndices( aNumberOfPoints, -1 )
            , mTargetPhysicalCoordinate( aPhysicalDimension, aNumberOfPoints )
            , mTargetParametricCoordinate( aPhysicalDimension - 1, aNumberOfPoints )
            , mTargetSideSetIndices( aNumberOfPoints, -1 )
            , mNormals( aPhysicalDimension, aNumberOfPoints )
            , mReferenceNormals( aPhysicalDimension, aNumberOfPoints )
            , mSignedDistance( aNumberOfPoints )
            , mNormalsNonlinear( aPhysicalDimension, aNumberOfPoints )
            , mSourcePhysicalCoordinateNonlinear( aPhysicalDimension, aNumberOfPoints )
    {
        // Assert all parametric coordinates are in bounds at construction
        for ( uint i = 0; i < aNumberOfPoints; ++i )
        {
            if ( mTargetParametricCoordinate.n_rows() == 2 )
            {
                assert_param_in_bounds_box( mTargetParametricCoordinate.get_column( i ), "MappingResult::MappingResult (box)" );
            }
            else if ( mTargetParametricCoordinate.n_rows() == 3 )
            {
                assert_param_in_bounds_simplex( mTargetParametricCoordinate.get_column( i ), "MappingResult::MappingResult (simplex)" );
            }
        }
    }

    Json MappingResult::to_json()
    {
        Json tMappingResult;

        auto& tResults = tMappingResult.add_child( "results", Json() );
        for ( size_t iPoint = 0; iPoint < mTargetCellIndices.size(); iPoint++ )
        {
            Json tPoint;
            tPoint.add_child( "source_coordinate", moris::to_json( trans( mSourcePhysicalCoordinate.get_column( iPoint ) ) ) );
            tPoint.add_child( "target_coordinate", moris::to_json( trans( mTargetPhysicalCoordinate.get_column( iPoint ) ) ) );
            tPoint.add_child( "target_parametric_coordinate", moris::to_json( trans( mTargetParametricCoordinate.get_column( iPoint ) ) ) );
            tPoint.add_child( "normal", moris::to_json( trans( mNormals.get_column( iPoint ) ) ) );
            tPoint.add( "distance", mSignedDistance( iPoint ) );
            tPoint.add( "source_cell_index", iPoint );
            tPoint.add( "target_cell_index", mTargetCellIndices( iPoint ) );
            tPoint.add( "target_side_set_index", mTargetSideSetIndices( iPoint ) );

            tResults.push_back( { "", tPoint } );
        }

        return tMappingResult;
    }
}    // namespace moris::mtk
