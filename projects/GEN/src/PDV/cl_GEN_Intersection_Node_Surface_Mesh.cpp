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
        aInterfaceGeometry->shift( -1.0 * trans( aFirstParentNode.get_global_coordinates() ) );

        // step 2: rotate the object
        // get unit axis to rotate to
        Matrix< DDRMat > tTransformationMatrix( 3, 3 );

        Matrix< DDRMat > tParentVector = aSecondParentNode.get_global_coordinates() - aFirstParentNode.get_global_coordinates();

        // augment with zero if 2D
        if ( tParentVector.numel() == 2 )
        {
            tParentVector.reshape( 3, 1 );
            tParentVector( 2, 0 ) = 0.0;
        }

        real tParentVectorNorm = norm( tParentVector );

        tParentVector = tParentVector / tParentVectorNorm;

        // create vector orthogonal to parent vector and coordinate axis
        // in 2D, this vector is the z axis
        tTransformationMatrix.set_column( 2, cross( tParentVector, { { 1.0, 0.0, 0.0 } } ) );
        aRotationAxis = 0;
        if ( norm( tTransformationMatrix.get_column( 2 ) ) < MORIS_REAL_EPS )
        {
            tTransformationMatrix.set_column( 2, cross( tParentVector, { { 0.0, 1.0, 0.0 } } ) );
            aRotationAxis = 1;

            // rotate along z axis only if basis is 3D
            if ( norm( tTransformationMatrix.get_column( 2 ) ) < MORIS_REAL_EPS && aFirstParentNode.get_global_coordinates().numel() > 2 )
            {
                tTransformationMatrix.set_column( 2, cross( tParentVector, { { 0.0, 0.0, 1.0 } } ) );
                aRotationAxis = 2;
            }
        }
        tTransformationMatrix.set_column( 2, tTransformationMatrix.get_column( 2 ) / norm( tTransformationMatrix.get_column( 2 ) ) );

        // create a second vector orthogonal to parent vector and first basis
        tTransformationMatrix.set_column( 1, cross( tParentVector, tTransformationMatrix.get_column( 2 ) ) );

        // the third vector of the transformation matrix is the parent vector
        tTransformationMatrix.set_column( 0, tParentVector );

        // trim the transformation matrix if 2D
        if ( aInterfaceGeometry->get_dimension() == 2 )
        {
            tTransformationMatrix.resize( 2, 2 );
        }

        // rotate the object
        aInterfaceGeometry->rotate( tTransformationMatrix );

        // step 3: scale the object
        Matrix< DDRMat > tScaling( aInterfaceGeometry->get_dimension(), 1 );
        tScaling.fill( 2.0 / tParentVectorNorm );
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
        aInterfaceGeometry->reset_coordinates();

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
}    // namespace moris::ge
