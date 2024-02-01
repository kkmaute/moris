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
#include "cl_GEN_Parent_Node.hpp"
#include "fn_cross.hpp"
#include "fn_eye.hpp"

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
        sdf::Object_Region tRegion = raycast_point( *this, aNodeCoordinates );

        switch ( tRegion )
        {
            case sdf::Object_Region::INSIDE:
            {
                return Geometric_Region::NEGATIVE;
            }
            case sdf::Object_Region::OUTSIDE:
            {
                return Geometric_Region::POSITIVE;
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
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Intersection_Node* Surface_Mesh_Geometry::create_intersection_node(
            uint                     aNodeIndex,
            const Cell< Node* >&     aBackgroundNodes,
            const Parent_Node&       aFirstParentNode,
            const Parent_Node&       aSecondParentNode,
            mtk::Geometry_Type       aBackgroundGeometryType,
            mtk::Interpolation_Order aBackgroundInterpolationOrder )
    {
        // Create linear intersection node
        return new Intersection_Node_Surface_Mesh(
                aNodeIndex,
                aBackgroundNodes,
                aFirstParentNode,
                aSecondParentNode,
                aBackgroundGeometryType,
                aBackgroundInterpolationOrder,
                *this );
    }
    
    //--------------------------------------------------------------------------------------------------------------
    
    real Surface_Mesh_Geometry::compute_intersection_local_coordinate(
            const Cell< Node* >& aBackgroundNodes,
            const Parent_Node&   aFirstParentNode,
            const Parent_Node&   aSecondParentNode )
    {
        // transform the interface geometry to local coordinates
        this->transform_surface_mesh_to_local_coordinate( aFirstParentNode, aSecondParentNode );
        
        // Compute the distance to the facets
        Matrix< DDRMat > tCastPoint( this->get_dimension(), 1 );
        tCastPoint.fill( 0.0 );
        Cell< real > tLocalCoordinate = sdf::compute_distance_to_facets( *this, tCastPoint, 0 );

        // shift local coordinate to be between -1 and 1
        for ( uint iIntersection = 0; iIntersection < tLocalCoordinate.size(); iIntersection++ )
        {
            tLocalCoordinate( iIntersection ) += -1.0;
        }

        // reset the object
        this->reset_coordinates();

        if ( tLocalCoordinate.size() == 0 )
        {
            return MORIS_REAL_MAX;
        }
        // FIXME: the case where 3 or more intersections occur between the two parents needs to be carefully considered
        else if ( tLocalCoordinate.size() > 2 )
        {
            MORIS_ERROR( tLocalCoordinate( 2 ) <= 1.0, "GEN - Intersection Node Surface Mesh: Parent nodes are in different geometric regions, and multiple intersections detected along parent edge." );
        }

        return tLocalCoordinate( 0 );
    }

    //--------------------------------------------------------------------------------------------------------------
    
    void Surface_Mesh_Geometry::transform_surface_mesh_to_local_coordinate(
            const Parent_Node& aFirstParentNode,
            const Parent_Node& aSecondParentNode )
    {
        // step 1: shift the object so the first parent is at the origin
        Matrix< DDRMat > tFirstParentNodeGlobalCoordinates = aFirstParentNode.get_global_coordinates();
        Cell< real >     tShift( this->get_dimension() );
        MORIS_ASSERT( tFirstParentNodeGlobalCoordinates.numel() == tShift.size(), "Intersection Node Surface Mesh::transform_mesh_to_local_coordinates() inconsistent parent node and interface geometry dimensions." );
        for ( uint iCoord = 0; iCoord < tShift.size(); iCoord++ )
        {
            tShift( iCoord ) = -1.0 * tFirstParentNodeGlobalCoordinates( iCoord );
        }
        this->shift( tShift );

        // step 2: rotate the object
        // get unit axis to rotate to
        Matrix< DDRMat > tParentVector = aSecondParentNode.get_global_coordinates() - aFirstParentNode.get_global_coordinates();

        // augment with zero if 2D
        if ( tParentVector.numel() == 2 )
        {
            tParentVector.reshape( 3, 1 );
            tParentVector( 2, 0 ) = 0.0;
        }

        real tParentVectorNorm = norm( tParentVector );

        tParentVector = tParentVector / tParentVectorNorm;

        // create vector orthogonal to parent vector and cast axis
        // in 2D, this vector is the z axis
        Matrix< DDRMat > tRotationMatrix( 3, 1 );
        Matrix< DDRMat > tCastAxis = { { 1.0 }, { 0.0 }, { 0.0 } };

        if ( norm( tParentVector + tCastAxis ) < this->get_intersection_tolerance() )
        {
            tRotationMatrix = { { -1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 } };
        }
        else
        {
            Matrix< DDRMat > tAntisymmetricCrossProduct = { { 0, tParentVector( 1 ), tParentVector( 2 ) },
                { -tParentVector( 1 ), 0.0, 0.0 },
                { -tParentVector( 2 ), 0.0, 0.0 } };

            Matrix< DDRMat > tAntisymmetricCrossProductSquared = { { -std::pow( tParentVector( 1 ), 2 ) - std::pow( tParentVector( 2 ), 2 ), 0.0, 0.0 },
                { 0.0, -std::pow( tParentVector( 1 ), 2 ), -tParentVector( 1 ) * tParentVector( 2 ) },
                { 0.0, -tParentVector( 1 ) * tParentVector( 2 ), -std::pow( tParentVector( 2 ), 2 ) } };

            tRotationMatrix = eye( 3, 3 ) + tAntisymmetricCrossProduct + ( 1 / ( 1 + tParentVector( 0 ) ) ) * tAntisymmetricCrossProductSquared;
        }

        // check that the rotation matrix is correct by ensuring the parent vector was rotated to the x axis
        MORIS_ASSERT( norm( tRotationMatrix * tParentVector - tCastAxis ) < this->get_intersection_tolerance(), "Rotation matrix should rotate the parent vector to the x axis." );

        // trim the transformation matrix if 2D
        if ( this->get_dimension() == 2 )
        {
            tRotationMatrix.resize( 2, 2 );
        }

        // rotate the object
        this->rotate( tRotationMatrix );

        // step 3: scale the object
        Cell< real > tScaling( this->get_dimension(), 2.0 / tParentVectorNorm );
        this->scale( tScaling );
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
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Surface_Mesh_Geometry::reset_nodal_data( mtk::Interpolation_Mesh* aInterpolationMesh )
    {
        // TODO BRENDAN
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
    }
    
    //--------------------------------------------------------------------------------------------------------------

    void
    Surface_Mesh_Geometry::get_design_info(
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aCoordinates,
            Cell< real >&           aOutputDesignInfo )
    {
        // TODO BRENDAN
        aOutputDesignInfo.resize( 0 );
    }
    
    //--------------------------------------------------------------------------------------------------------------

    bool Surface_Mesh_Geometry::intended_discretization()
    {
        return ( mParameters.mDiscretizationIndex >= 0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    moris_index
    Surface_Mesh_Geometry::get_discretization_mesh_index()
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
    
    //--------------------------------------------------------------------------------------------------------------
    
}    // namespace moris::ge
