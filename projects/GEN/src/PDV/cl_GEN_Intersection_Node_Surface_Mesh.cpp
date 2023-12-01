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

#include "fn_norm.hpp"
#include "fn_eye.hpp"

namespace moris::ge
{
    Intersection_Node_Surface_Mesh::Intersection_Node_Surface_Mesh(
            std::shared_ptr< Intersection_Node >     aFirstParentNode,
            std::shared_ptr< Intersection_Node >     aSecondParentNode,
            uint                                     aFirstParentNodeIndex,
            uint                                     aSecondParentNodeIndex,
            const Matrix< DDRMat >&                  aFirstParentNodeLocalCoordinates,
            const Matrix< DDRMat >&                  aSecondParentNodeLocalCoordinates,
            std::shared_ptr< Surface_Mesh_Geometry > aInterfaceGeometry )
            : Intersection_Node(
                    compute_local_coordinate(
                            aFirstParentNodeLocalCoordinates,
                            aSecondParentNodeLocalCoordinates ),
                    aFirstParentNode,
                    aSecondParentNode,
                    aFirstParentNodeIndex,
                    aSecondParentNodeIndex,
                    { { -1 } },
                    { { 1 } },
                    { { aFirstParentNodeIndex, aSecondParentNodeIndex } },
                    { { aFirstParentNodeLocalCoordinates }, { aSecondParentNodeLocalCoordinates } },
                    Element_Interpolation_Type::Linear_1D )
            , mInterfaceGeometry( aInterfaceGeometry )
    {
    }

    Matrix< DDRMat >
    Intersection_Node_Surface_Mesh::compute_global_coordinates()
    {
        // TODO
        MORIS_ERROR( false, "Intersection_Node_Surface_Mesh - compute_global_coordinates() not implemented yet." );
        return { { std::numeric_limits< double >::quiet_NaN() } };
    }

    bool
    Intersection_Node_Surface_Mesh::determine_is_intersected(
            const Element_Interpolation_Type aAncestorBasisFunction,
            const Matrix< DDRMat >&          aFirstParentNodeLocalCoordinates,
            const Matrix< DDRMat >&          aSecondParentNodeLocalCoordinates )
    {
        // lock interface geometry
         std::shared_ptr< Surface_Mesh_Geometry > tLockedInterfaceGeometry = mInterfaceGeometry.lock();

        // compute the global coordinates of the parent nodes
        Matrix< DDRMat > tFirstParentNodeGlobalCoordinates;
        Matrix< DDRMat > tSecondParentNodeGlobalCoordinates;

        // determine if the parents are on the interface
        bool tFirstParentOnInterface = ( tLockedInterfaceGeometry->get_geometric_region( 0, tFirstParentNodeGlobalCoordinates ) == Geometric_Region::INTERFACE );
        bool tSecondParentOnInterface = ( tLockedInterfaceGeometry->get_geometric_region( 0, tSecondParentNodeGlobalCoordinates ) == Geometric_Region::INTERFACE );

        // FIXME: add assert statements
        if( tFirstParentOnInterface or tSecondParentOnInterface )
        {
            return true;
        }
        else
        {
            return std:: abs( mLocalCoordinate ) <= 1.0;
        }
    }

    Matrix< DDRMat >
    Intersection_Node_Surface_Mesh::compute_raycast_rotation()
    {
        // FIXME: mParentVector interpolates the parent global coordinates with the ancestor coordinates
        // which may mean the parent vector is computed improperly for surface mesh geometries. Look here if problems arise
        Matrix< DDRMat > tUnitParentVector = mParentVector / norm( mParentVector );

        real tCosineAngle = tUnitParentVector( 2 );

        Matrix< DDRMat > tRotationMatrix( mParentVector.numel(), mParentVector.numel() );
        if ( abs( tCosineAngle ) - 1.0 < MORIS_REAL_EPS )
        {
            // if the parent vector points in the opposite direction of the z axis, reflect it
            tRotationMatrix = { { 1.0, 0, 0 }, { 0, 1.0, 0 }, { 0, 0, -1.0 } };
        }
        else
        {
            // else, the rotation matrix is defined by rodrigues' rotation formula
            Matrix< DDRMat > tAntiSymmetricCrossProduct = { { 0, 0, tUnitParentVector( 0 ) },
                { 0, 0, tUnitParentVector( 1 ) },
                { -tUnitParentVector( 0 ), -tUnitParentVector( 1 ), 0 } };

            tRotationMatrix = eye( mParentVector.numel(), mParentVector.numel() ) + tAntiSymmetricCrossProduct + tAntiSymmetricCrossProduct * tAntiSymmetricCrossProduct / ( 1 + tCosineAngle );
        }
        return tRotationMatrix;
    }

    real
    Intersection_Node_Surface_Mesh::compute_local_coordinate(
            const Matrix< DDRMat >& aFirstParentNodeCoordinates,
            const Matrix< DDRMat >& aSecondParentNodeCoordinates )
    {
        // compute the rotation matrix
        Matrix< DDRMat > tRotationMatrix = this->compute_raycast_rotation();

        // rotate the geometry and the cast point
        std::shared_ptr< Surface_Mesh_Geometry > tLockedInterfaceGeometry = mInterfaceGeometry.lock();
        tLockedInterfaceGeometry->rotate_object( tRotationMatrix );
        Matrix< DDRMat > tRotatedFirstParentCoordinates = tRotationMatrix * aFirstParentNodeCoordinates;

        // Compute the distance to the facets
        // FIXME: the interface geometry needs to be put in local coordinate frame
        Cell< real > tLocalCoordinate = sdf::compute_distance_to_facets( *tLockedInterfaceGeometry, tRotatedFirstParentCoordinates, 2 );

        // Find all intersection locations and return the closest one
        // FIXME: this ignores the possibility that there are multiple intersection locations between each facet
        return tLocalCoordinate( 0 );
    }

    void
    Intersection_Node_Surface_Mesh::append_dcoordinate_dadv( Matrix< DDRMat >& aCoordinateSensitivities, const Matrix< DDRMat >& aSensitivityFactor )
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

    real
    Intersection_Node_Surface_Mesh::get_dxi_dfield_from_ancestor( uint aAncestorIndex )
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
