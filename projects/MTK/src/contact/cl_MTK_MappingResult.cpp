//
// Created by frank on 12/8/23.
//

#include "cl_MTK_MappingResult.hpp"
#include "cl_Json_Object.hpp"

namespace moris::mtk
{
    MappingResult::MappingResult( moris_index aSourceMeshIndex, uint aPhysicalDimension, uint aNumberOfPoints )
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
            , mDistances( aNumberOfPoints )
    {
    }

    Json MappingResult::to_json()
    {
        Json tMappingResult;

        auto &tResults = tMappingResult.add_child( "results", Json() );
        for ( size_t iPoint = 0; iPoint < mTargetCellIndices.size(); iPoint++ )
        {
            Json tPoint;
            tPoint.add_child( "source_coordinate", moris::to_json( trans( mSourcePhysicalCoordinate.get_column( iPoint ) ) ) );
            tPoint.add_child( "target_coordinate", moris::to_json( trans( mTargetPhysicalCoordinate.get_column( iPoint ) ) ) );
            tPoint.add_child( "target_parametric_coordinate", moris::to_json( trans( mTargetParametricCoordinate.get_column( iPoint ) ) ) );
            tPoint.add_child( "normal", moris::to_json( trans( mNormals.get_column( iPoint ) ) ) );
            tPoint.add( "distance", mDistances( iPoint ) );
            tPoint.add( "source_cell_index", iPoint );
            tPoint.add( "target_cell_index", mTargetCellIndices( iPoint ) );
            tPoint.add( "target_side_set_index", mTargetSideSetIndices( iPoint ) );

            tResults.push_back( { "", tPoint } );
        }

        return tMappingResult;
    }
}    // namespace moris::mtk
