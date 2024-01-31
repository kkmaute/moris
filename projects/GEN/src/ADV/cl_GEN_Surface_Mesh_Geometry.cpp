/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Surface_Mesh_Geometry.cpp
 *
 */

#include <string>

#include "cl_GEN_Surface_Mesh_Geometry.hpp"
#include "cl_GEN_Intersection_Node_Surface_Mesh.hpp"

namespace moris::ge
{

    //--------------------------------------------------------------------------------------------------------------

    Surface_Mesh_Parameters::Surface_Mesh_Parameters( const ParameterList& aParameterList )
            : Field_Parameters( aParameterList )
            , Design_Parameters( aParameterList )
            , mFilePath( aParameterList.get< std::string >( "file_path" ) )
            , mIntersectionTolerance( aParameterList.get< real >( "intersection_tolerance" ) )
    {
        string_to_cell( aParameterList.get< std::string >( "offset" ), mOffsets );
        string_to_cell( aParameterList.get< std::string >( "scale" ), mScale );
    }

    //--------------------------------------------------------------------------------------------------------------

    Surface_Mesh_Geometry::Surface_Mesh_Geometry( Surface_Mesh_Parameters aParameters )
            : Geometry( aParameters, aParameters.mIntersectionTolerance )
            , Object( aParameters.mFilePath, aParameters.mIntersectionTolerance, aParameters.mOffsets, aParameters.mScale )
            , mParameters( aParameters )
    {
        mName = aParameters.mFilePath.substr( aParameters.mFilePath.find_last_of( "/" ) + 1, aParameters.mFilePath.find_last_of( "." ) );
    }

    //--------------------------------------------------------------------------------------------------------------

    Geometric_Region Surface_Mesh_Geometry::get_geometric_region(
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aNodeCoordinates )
    {
        // BRENDAN DELETE PRINTS
        PRINT( aNodeCoordinates );

        sdf::Object_Region tRegion = raycast_point( *this, aNodeCoordinates );

        std::cout << "geometric region: " << tRegion << std::endl;

        switch ( tRegion )
        {
            case sdf::Object_Region::INSIDE:
            {
                return Geometric_Region::NEGATIVE;
                break;
            }
            case sdf::Object_Region::OUTSIDE:
            {
                return Geometric_Region::POSITIVE;
                break;
            }
            case sdf::Object_Region::INTERFACE:
            {
                return Geometric_Region::INTERFACE;
                break;
            }
            default:
            {
                MORIS_ERROR( false, "Unexpected sdf::Object_Region of %d returned from raycast.", tRegion );
                return Geometric_Region::INTERFACE;
                break;
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Intersection_Node* Surface_Mesh_Geometry::create_intersection_node(
            uint                     aNodeIndex,
            const Cell< Node* >&     aBaseNodes,
            const Parent_Node&       aFirstParentNode,
            const Parent_Node&       aSecondParentNode,
            mtk::Geometry_Type       aBackgroundGeometryType,
            mtk::Interpolation_Order aBackgroundInterpolationOrder )
    {
        // Create linear intersection node
        return new Intersection_Node_Surface_Mesh(
                aNodeIndex,
                aBaseNodes,
                aFirstParentNode,
                aSecondParentNode,
                aBackgroundGeometryType,
                aBackgroundInterpolationOrder,
                shared_from_this() );
    }

    //--------------------------------------------------------------------------------------------------------------

    Cell< std::shared_ptr< mtk::Field > > Surface_Mesh_Geometry::get_mtk_fields()
    {
        // TODO BRENDAN: maybe?
        return {};
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Surface_Mesh_Geometry::import_advs( sol::Dist_Vector* aOwnedADVs )
    {
        // TODO BRENDAN
        return;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Surface_Mesh_Geometry::reset_nodal_data( mtk::Interpolation_Mesh* aInterpolationMesh )
    {
        // TODO BRENDAN
        return;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Surface_Mesh_Geometry::discretize(
            mtk::Mesh_Pair          aMeshPair,
            sol::Dist_Vector*       aOwnedADVs,
            const Matrix< DDSMat >& aSharedADVIds,
            uint                    aADVOffsetID )
    {
        // TODO BRENDAN
        return;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Surface_Mesh_Geometry::discretize(
            std::shared_ptr< mtk::Field > aMTKField,
            mtk::Mesh_Pair                aMeshPair,
            sol::Dist_Vector*             aOwnedADVs,
            const Matrix< DDSMat >&       aSharedADVIds,
            uint                          aADVOffsetID )
    {
        // TODO BRENDAN
        return;
    }

    void
    Surface_Mesh_Geometry::get_design_info(
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aCoordinates,
            Cell< real >&           aOutputDesignInfo )
    {
        // TODO BRENDAN
        aOutputDesignInfo.resize( 0 );
    }

    bool Surface_Mesh_Geometry::intended_discretization()
    {
        return ( mParameters.mDiscretizationIndex >= 0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    moris_index
    Surface_Mesh_Geometry::get_discretization_mesh_index() const
    {
        MORIS_ASSERT( mParameters.mDiscretizationIndex >= 0,
                "A discretization is not intended for this field. Check this with intended_discretization() first." );

        return mParameters.mDiscretizationIndex;
    }

    //--------------------------------------------------------------------------------------------------------------

    real
    Surface_Mesh_Geometry::get_discretization_lower_bound()
    {
        return mParameters.mDiscretizationLowerBound;
    }

    //--------------------------------------------------------------------------------------------------------------

    real
    Surface_Mesh_Geometry::get_discretization_upper_bound()
    {
        return mParameters.mDiscretizationUpperBound;
    }
}    // namespace moris::ge
