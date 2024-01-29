/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Intersection_Node_Surface_Mesh.cpp
 *
 */

#include <limits>
#include "cl_GEN_Intersection_Node_Surface_Mesh.hpp"
#include "cl_GEN_Surface_Mesh_Geometry.hpp"
#include "cl_GEN_Parent_Node.hpp"

#include "fn_norm.hpp"
#include "fn_trans.hpp"
#include "fn_cross.hpp"
#include "fn_eye.hpp"

namespace moris::ge
{
    Intersection_Node_Surface_Mesh::Intersection_Node_Surface_Mesh(
            uint                                     aNodeIndex,
            const Cell< Node* >&                     aBaseNodes,
            const Parent_Node&                       aFirstParentNode,
            const Parent_Node&                       aSecondParentNode,
            mtk::Geometry_Type                       aBackgroundGeometryType,
            mtk::Interpolation_Order                 aBackgroundInterpolationOrder,
            std::shared_ptr< Surface_Mesh_Geometry > aInterfaceGeometry )
            : Intersection_Node(
                    aNodeIndex,
                    aBaseNodes,
                    aFirstParentNode,
                    aSecondParentNode,
                    Intersection_Node_Surface_Mesh::compute_local_coordinate( aFirstParentNode, aSecondParentNode, aInterfaceGeometry ),
                    aBackgroundGeometryType,
                    aBackgroundInterpolationOrder,
                    aInterfaceGeometry )
            , mInterfaceGeometry( aInterfaceGeometry )
    {
    }

    void
    Intersection_Node_Surface_Mesh::transform_surface_mesh_to_local_coordinate(
            const Parent_Node&                       aFirstParentNode,
            const Parent_Node&                       aSecondParentNode,
            std::shared_ptr< Surface_Mesh_Geometry > aInterfaceGeometry,
            uint&                                    aRotationAxis )
    {
        // step 1: shift the object so the first parent is at the origin
        Matrix< DDRMat > tFirstParentNodeGlobalCoordinates = aFirstParentNode.get_global_coordinates();
        Cell< real >     tShift( aInterfaceGeometry->get_dimension() );
        MORIS_ASSERT( tFirstParentNodeGlobalCoordinates.numel() == tShift.size(), "Intersection Node Surface Mesh::transform_mesh_to_local_coordinates() inconsistent parent node and interface geometry dimensions." );
        for ( uint iCoord = 0; iCoord < tShift.size(); iCoord++ )
        {
            tShift( iCoord ) = -1.0 * tFirstParentNodeGlobalCoordinates( iCoord );
        }
        aInterfaceGeometry->shift( tShift );

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

        if ( norm( tParentVector + tCastAxis ) < aInterfaceGeometry->get_intersection_tolerance() )
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
        MORIS_ASSERT( norm( tRotationMatrix * tParentVector - tCastAxis ) < aInterfaceGeometry->get_intersection_tolerance(), "Rotation matrix should rotate the parent vector to the x axis." );

        // trim the transformation matrix if 2D
        if ( aInterfaceGeometry->get_dimension() == 2 )
        {
            tRotationMatrix.resize( 2, 2 );
        }

        // rotate the object
        aInterfaceGeometry->rotate( tRotationMatrix );

        // step 3: scale the object
        Cell< real > tScaling( aInterfaceGeometry->get_dimension(), 2.0 / tParentVectorNorm );
        aInterfaceGeometry->scale( tScaling );
    }

    real Intersection_Node_Surface_Mesh::compute_local_coordinate(
            const Parent_Node&                       aFirstParentNode,
            const Parent_Node&                       aSecondParentNode,
            std::shared_ptr< Surface_Mesh_Geometry > aInterfaceGeometry )
    {
        // transform the interface geometry to local coordinates
        uint tRotatedAxis;
        transform_surface_mesh_to_local_coordinate( aFirstParentNode, aSecondParentNode, aInterfaceGeometry, tRotatedAxis );

        // Compute the distance from the first parent edge to each facet in the local coordinate frame
        Matrix< DDRMat > tCastPoint( aInterfaceGeometry->get_dimension(), 1 );
        tCastPoint.fill( 0.0 );
        Cell< real > tLocalCoordinate = sdf::compute_distance_to_facets( *aInterfaceGeometry, tCastPoint, 0 );

        // shift local coordinate to be between -1 and 1
        for ( uint iIntersection = 0; iIntersection < tLocalCoordinate.size(); iIntersection++ )
        {
            tLocalCoordinate( iIntersection ) += -1.0;
        }

        // reset the object
        aInterfaceGeometry->reset_coordinates();

        // return not intersected if raycast found no facets
        if ( tLocalCoordinate.size() == 0 )
        {
            return MORIS_REAL_MAX;
        }
        // FIXME: need logic for odd number of intersections between parent edges
        else if ( tLocalCoordinate.size() > 2
                  && tLocalCoordinate.size() % 2
                  && std::all_of( tLocalCoordinate.begin(), tLocalCoordinate.end(), []( real tCoord ) { return not( tCoord < -1 ) and not( tCoord > 1 ); } ) )
        {
            MORIS_ERROR( false,
                    "Intersection_Node_Surface_Mesh::compute_local_coordinate() Parent nodes in different regions, and 3 or more facet intersections detected." );
        }

        real tEdgeCoordinate = std::abs( tLocalCoordinate( 0 ) + 1 ) <  aInterfaceGeometry->get_intersection_tolerance() ? -1
                             : std::abs( tLocalCoordinate( 0 ) - 1 ) <  aInterfaceGeometry->get_intersection_tolerance() ? 1
                                                                                                                               : tLocalCoordinate( 0 );
        return tEdgeCoordinate;
    }

    void Intersection_Node_Surface_Mesh::append_dcoordinate_dadv( Matrix< DDRMat >& aCoordinateSensitivities, const Matrix< DDRMat >& aSensitivityFactor )
    {
        // TODO
        MORIS_ERROR( false, "Intersection_Node_Surface_Mesh::get_dcoordinate_dadv() not implemented yet." );

        return;
    };

    Matrix< DDSMat >
    Intersection_Node_Surface_Mesh::get_coordinate_determining_adv_ids()
    {
        // TODO
        MORIS_ERROR( false, "Intersection_Node_Surface_Mesh - get_coordinate_determining_adv_ids() not implemented yet." );
        return { { -1 } };
    }
}    // namespace moris::ge
