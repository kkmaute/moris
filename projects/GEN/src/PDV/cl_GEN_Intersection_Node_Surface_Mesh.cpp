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
#include "cl_GEN_Parent_Node.hpp"
#include "cl_GEN_Intersection_Node_Surface_Mesh.hpp"
#include "cl_GEN_Surface_Mesh_Geometry.hpp"

#include "fn_norm.hpp"
#include "fn_trans.hpp"
#include "fn_cross.hpp"

namespace moris::ge
{
    Intersection_Node_Surface_Mesh::Intersection_Node_Surface_Mesh(
            uint                                     aNodeIndex,
            const Cell< Node* >&                     aBaseNodes,
            const Parent_Node&                       aFirstParentNode,
            const Parent_Node&                       aSecondParentNode,
            mtk::Geometry_Type                       aBaseGeometryType,
            std::shared_ptr< Surface_Mesh_Geometry > aInterfaceGeometry )
            : Intersection_Node(
                    aNodeIndex,
                    aBaseNodes,
                    aFirstParentNode,
                    aSecondParentNode,
                    Intersection_Node_Surface_Mesh::compute_local_coordinate( aFirstParentNode, aSecondParentNode, aInterfaceGeometry ),
                    aBaseGeometryType,
                    aInterfaceGeometry )
            , mInterfaceGeometry( aInterfaceGeometry )
    {
    }

    bool
    Intersection_Node_Surface_Mesh::determine_is_intersected()
    {
        // lock interface geometry
        std::shared_ptr< Surface_Mesh_Geometry > tLockedInterfaceGeometry = mInterfaceGeometry.lock();

        // determine if the parents are on the interface
        bool tFirstParentOnInterface  = ( tLockedInterfaceGeometry->get_geometric_region( 0, mFirstParentNode.get_global_coordinates() ) == Geometric_Region::INTERFACE );
        bool tSecondParentOnInterface = ( tLockedInterfaceGeometry->get_geometric_region( 0, mSecondParentNode.get_global_coordinates() ) == Geometric_Region::INTERFACE );

        // FIXME: add assert statements
        if ( tFirstParentOnInterface or tSecondParentOnInterface )
        {
            return true;
        }
        else
        {
            return std::abs( this->get_local_coordinate() ) <= 1.0;
        }
    }

    void
    Intersection_Node_Surface_Mesh::transform_surface_mesh_to_local_coordinate(
            const Parent_Node&                       aFirstParentNode,
            const Parent_Node&                       aSecondParentNode,
            std::shared_ptr< Surface_Mesh_Geometry > aInterfaceGeometry,
            uint&                                    aRotationAxis )
    {
        // step 1: shift the object so the first parent is at the origin
        aInterfaceGeometry->shift_object( -1.0 * trans( aFirstParentNode.get_global_coordinates() ) );

        // BRENDAN
        std::cout << "shift done" << std::endl;

        // step 2: rotate the object
        // get unit axis to rotate to
        Matrix< DDRMat > tParentVector( 3, 1 );
        if ( aFirstParentNode.get_global_coordinates().numel() == 2 )
        {
            // BRENDAN [2ndParent 0] - [1stParent 0] is what should be here
            tParentVector =
        }
        else
        {
            tParentVector = aSecondParentNode.get_global_coordinates() - aFirstParentNode.get_global_coordinates();
        }
        tParentVector = tParentVector / norm( tParentVector );

        // create vector orthogonal to parent vector and coordinate axis
        Matrix< DDRMat > tFirstBasis = cross( tParentVector, { { 1.0, 0.0, 0.0 } } );
        aRotationAxis                = 0;
        if ( norm( tFirstBasis ) < MORIS_REAL_EPS )
        {
            tFirstBasis   = cross( tParentVector, { { 0.0, 1.0, 0.0 } } );
            aRotationAxis = 1;

            // rotate along z axis only if basis is 3D
            if ( norm( tFirstBasis ) < MORIS_REAL_EPS && aFirstParentNode.get_global_coordinates().numel() > 2 )
            {
                tFirstBasis   = cross( tParentVector, { { 0.0, 0.0, 1.0 } } );
                aRotationAxis = 2;
            }
        }

        // create second vector orthogonal to parent vector and first basis
        Matrix< DDRMat > tSecondBasis = cross( tParentVector, tFirstBasis );

        // BRENDAN create transformation matrix from set of vectors
        // basically: tRotationMatrix = [ tParentVector tSecondBasis tFirstBasis ]
        Matrix< DDRMat > tRotationMatrix( 3, 3 );

        // rotate the object
        aInterfaceGeometry->rotate_object( trans( tRotationMatrix ) );

        // step 3: scale the object
        Matrix< DDRMat > tScaling( aInterfaceGeometry->get_dimension(), 1 );
        tScaling.fill( 0.5 );
        aInterfaceGeometry->scale_object( tScaling );

        // BRENDAN
        std::cout << "scaling happened already" << std::endl;
    }

    real Intersection_Node_Surface_Mesh::compute_local_coordinate(
            const Parent_Node&                       aFirstParentNode,
            const Parent_Node&                       aSecondParentNode,
            std::shared_ptr< Surface_Mesh_Geometry > aInterfaceGeometry )
    {
        // transform the interface geometry to local coordinates
        uint tRotatedAxis;
        transform_surface_mesh_to_local_coordinate( aFirstParentNode, aSecondParentNode, aInterfaceGeometry, tRotatedAxis );

        // Compute the distance to the facets
        Matrix< DDRMat > tCastPoint( aInterfaceGeometry->get_dimension(), 1 );
        tCastPoint.fill( 0.0 );
        Cell< real > tLocalCoordinate = sdf::compute_distance_to_facets( *aInterfaceGeometry, tCastPoint, tRotatedAxis );

        // shift local coordinate to be between -1 and 1
        for ( uint iIntersection = 0; iIntersection < tLocalCoordinate.size(); iIntersection++ )
        {
            tLocalCoordinate( iIntersection ) += -1.0;
        }

        // reset the object
        aInterfaceGeometry->reset_object_coordinates();

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

    void Intersection_Node_Surface_Mesh::append_dcoordinate_dadv( Matrix< DDRMat >& aCoordinateSensitivities, const Matrix< DDRMat >& aSensitivityFactor )
    {
        // TODO
        MORIS_ERROR( false, "Intersection_Node_Surface_Mesh - get_dcoordinate_dadv() not implemented yet." );

        return;
    };

    Matrix< DDSMat >
    Intersection_Node_Surface_Mesh::get_coordinate_determining_adv_ids()
    {
        // TODO
        MORIS_ERROR( false, "Intersection_Node_Surface_Mesh - get_coordinate_determining_adv_ids() not implemented yet." );
        return { { -1 } };
    }

    real Intersection_Node_Surface_Mesh::get_dxi_dfield_from_ancestor( uint aAncestorIndex )
    {
        // TODO
        MORIS_ERROR( false, "Intersection_Node_Surface_Mesh - get_dxi_dfield_from_ancestor() not implemented yet." );

        return std::numeric_limits< double >::quiet_NaN();
    }

    Matrix< DDRMat >
    Intersection_Node_Surface_Mesh::get_dxi_dcoordinate_first_parent()
    {
        // TODO
        MORIS_ERROR( false, "Intersection_Node_Surface_Mesh - get_dxi_dcoordinate_first_parent() not implemented yet." );
        return { { std::numeric_limits< double >::quiet_NaN() } };
    }

    Matrix< DDRMat >
    Intersection_Node_Surface_Mesh::get_dxi_dcoordinate_second_parent()
    {
        // TODO
        MORIS_ERROR( false, "Intersection_Node_Surface_Mesh - get_dxi_dcoordinate_second_parent() not implemented yet." );
        return { { std::numeric_limits< double >::quiet_NaN() } };
    }
}    // namespace moris::ge
