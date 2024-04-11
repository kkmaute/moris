/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Intersection_Node_Surface_Mesh.cpp
 *
 */

#include "cl_GEN_Intersection_Node_Surface_Mesh.hpp"
#include "cl_GEN_Parent_Node.hpp"
#include "cl_GEN_Surface_Mesh_Geometry.hpp"

#include "fn_norm.hpp"
#include "fn_eye.hpp"
#include "fn_trans.hpp"

namespace moris::gen
{
    Intersection_Node_Surface_Mesh::Intersection_Node_Surface_Mesh(
            uint                              aNodeIndex,
            const Vector< Background_Node* >& aBackgroundNodes,
            const Parent_Node&                aFirstParentNode,
            const Parent_Node&                aSecondParentNode,
            real                              aLocalCoordinate,
            sdf::Facet*                       aParentFacet,
            mtk::Geometry_Type                aBackgroundGeometryType,
            mtk::Interpolation_Order          aBackgroundInterpolationOrder,
            Surface_Mesh_Geometry&            aInterfaceGeometry )
            : Intersection_Node(
                    aNodeIndex,
                    aBackgroundNodes,
                    aFirstParentNode,
                    aSecondParentNode,
                    aLocalCoordinate,
                    aBackgroundGeometryType,
                    aBackgroundInterpolationOrder )
            , mInterfaceGeometry( aInterfaceGeometry )
    {
    }

    void
    Intersection_Node_Surface_Mesh::get_rotation_matrix( Matrix< DDRMat >& aRotationMatrix ) const
    {
        // get parent vector
        Matrix< DDRMat > tParentVector = this->get_first_parent_node().get_global_coordinates() - this->get_second_parent_node().get_global_coordinates();

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

        if ( norm( tParentVector + tCastAxis ) < mInterfaceGeometry.get_intersection_tolerance() )
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
    }

    //--------------------------------------------------------------------------------------------------------------

    Geometry& Intersection_Node_Surface_Mesh::get_interface_geometry()
    {

        return mInterfaceGeometry;
    }

    //--------------------------------------------------------------------------------------------------------------

    const Geometry& Intersection_Node_Surface_Mesh::get_interface_geometry() const
    {
        return mInterfaceGeometry;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    Intersection_Node_Surface_Mesh::compute_dxi_dfacet() const
    {
        // Get the parent vector and its norm from the intersection node
        Matrix< DDRMat > tParentVector     = this->get_first_parent_node().get_global_coordinates() - this->get_second_parent_node().get_global_coordinates();
        real             tParentVectorNorm = norm( tParentVector );

        // Get the rotation matrix from the intersection node
        Matrix< DDRMat > tRotationMatrix;
        this->get_rotation_matrix( tRotationMatrix );

        // get the normal vector rotated in the local coordinate frame
        Matrix< DDRMat > tNormalPrime = tRotationMatrix * mParentFacet->get_normal();

        // get the center vector in the local coordinate frame
        Matrix< DDRMat > tCenterPrime = 2.0 / ( 3.0 * tParentVectorNorm ) * tRotationMatrix * mParentFacet->get_center();

        // compute vector between facet vertices, and unnormalized normal vector
        Matrix< DDRMat > tFacetVector = mParentFacet->get_vertex_coords().get_row( 1 ) - mParentFacet->get_vertex_coords().get_row( 0 );
        Matrix< DDRMat > tNormal      = { { mParentFacet->get_vertex_coord( 0, 1 ) - mParentFacet->get_vertex_coord( 1, 1 ), mParentFacet->get_vertex_coord( 1, 0 ) - mParentFacet->get_vertex_coord( 0, 0 ) } };

        // derivative of normal vector
        Matrix< DDRMat > tdNormaldVertex1 = { { 1.0, -1.0 } };
        Matrix< DDRMat > tdNormaldVertex2 = { { -1.0, 1.0 } };

        // Sensitivity of the normal vector to the vertices
        Matrix< DDRMat > tdNormalPrimedVertex1 = tRotationMatrix * ( 1.0 / norm( tFacetVector ) * tdNormaldVertex1 - 0.5 * tNormal * ( trans( tFacetVector ) * 1.0 * eye( 2, 2 ) - eye( 2, 2 ) * ( tFacetVector ) ) );
        Matrix< DDRMat > tdNormalPrimedVertex2 = tRotationMatrix * ( 1.0 / norm( tFacetVector ) * tdNormaldVertex2 - 0.5 * tNormal * ( trans( tFacetVector ) * 1.0 * eye( 2, 2 ) - eye( 2, 2 ) * ( tFacetVector ) ) );


        Matrix< DDRMat > tdCenterdVertices = 2.0 / ( 3.0 * tParentVectorNorm ) * tRotationMatrix;

        return join_cols( ( tdNormalPrimedVertex1 * tCenterPrime + tdCenterdVertices * tNormalPrime ) * mParentFacet->get_normal()( 0 ) - tdNormalPrimedVertex1( 0 ) * ( tNormalPrime * tCenterPrime ),
                ( tdNormalPrimedVertex2 * tCenterPrime + tdCenterdVertices * tNormalPrime ) * mParentFacet->get_normal()( 0 ) - tdNormalPrimedVertex2( 0 ) * ( tNormalPrime * tCenterPrime ) );
    }

    //--------------------------------------------------------------------------------------------------------------

    void Intersection_Node_Surface_Mesh::append_dcoordinate_dadv( Matrix< DDRMat >& aCoordinateSensitivities, const Matrix< DDRMat >& aSensitivityFactor ) const
    {
        // Get parent nodes
        const Basis_Node& tFirstParentNode  = this->get_first_parent_node();
        const Basis_Node& tSecondParentNode = this->get_second_parent_node();

        // Compute parent vector
        Matrix< DDRMat > tParentVector = trans( tSecondParentNode.get_global_coordinates() - tFirstParentNode.get_global_coordinates() );

        Matrix< DDRMat > tLocalCoordinateFacetVertexSensitivities = this->compute_dxi_dfacet();

        Matrix< DDRMat > tSensitivitiesToAdd = .5 * aSensitivityFactor * tParentVector * ( tLocalCoordinateFacetVertexSensitivities.get_column( 0 ) * mInterfaceGeometry.get_dvertex_dadv( mParentFacet->get_vertex_inds()( 0 ) ) + tLocalCoordinateFacetVertexSensitivities.get_column( 1 ) * mInterfaceGeometry.get_dvertex_dadv( mParentFacet->get_vertex_inds()( 1 ) ) );

        // Resize sensitivities
        uint tJoinedSensitivityLength = aCoordinateSensitivities.n_cols();
        aCoordinateSensitivities.resize( tSensitivitiesToAdd.n_rows(),
                tJoinedSensitivityLength + tSensitivitiesToAdd.n_cols() );

        // Join sensitivities
        for ( uint iCoordinateIndex = 0; iCoordinateIndex < tSensitivitiesToAdd.n_rows(); iCoordinateIndex++ )
        {
            for ( uint iAddedSensitivity = 0; iAddedSensitivity < tSensitivitiesToAdd.n_cols(); iAddedSensitivity++ )
            {
                aCoordinateSensitivities( iCoordinateIndex, tJoinedSensitivityLength + iAddedSensitivity ) =
                        tSensitivitiesToAdd( iCoordinateIndex, iAddedSensitivity );
            }
        }

        // Add first parent coordinate sensitivities
        if ( tFirstParentNode.depends_on_advs() )
        {
            Matrix< DDRMat > tLocCoord          = ( 1.0 - this->get_local_coordinate() ) * eye( tParentVector.n_rows(), tParentVector.n_rows() );
            Matrix< DDRMat > tSensitivityFactor = 0.5 * ( tLocCoord + tParentVector * this->get_dxi_dcoordinate_first_parent() );
            tFirstParentNode.append_dcoordinate_dadv( aCoordinateSensitivities, tSensitivityFactor );
        }

        // Add second parent coordinate sensitivities
        if ( tSecondParentNode.depends_on_advs() )
        {
            Matrix< DDRMat > tLocCoord          = ( 1.0 + this->get_local_coordinate() ) * eye( tParentVector.n_rows(), tParentVector.n_rows() );
            Matrix< DDRMat > tSensitivityFactor = 0.5 * ( tLocCoord + tParentVector * this->get_dxi_dcoordinate_second_parent() );
            tSecondParentNode.append_dcoordinate_dadv( aCoordinateSensitivities, tSensitivityFactor );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDSMat >
    Intersection_Node_Surface_Mesh::get_coordinate_determining_adv_ids() const
    {
        // Initialize ADV IDs
        Matrix< DDSMat > tCoordinateDeterminingADVIDs;

        // Get sensitivity values from parents
        const Vector< Basis_Node >& tFieldBasisNodes = this->get_locator_nodes();
        for ( uint iFieldBasisNode = 0; iFieldBasisNode < tFieldBasisNodes.size(); iFieldBasisNode++ )
        {
            // Get geometry field sensitivity with respect to ADVs
            const Matrix< DDSMat >& tAncestorADVIDs = mInterfaceGeometry.get_determining_adv_ids(
                    tFieldBasisNodes( iFieldBasisNode ).get_index(),
                    tFieldBasisNodes( iFieldBasisNode ).get_global_coordinates() );

            // Join IDs
            Intersection_Node::join_adv_ids( tCoordinateDeterminingADVIDs, tAncestorADVIDs );
        }

        // Get parent nodes
        const Basis_Node& tFirstParentNode = this->get_first_parent_node();
        const Basis_Node& tSecondParentNode = this->get_second_parent_node();

        // Add parent IDs
        if ( tFirstParentNode.depends_on_advs() )
        {
            Intersection_Node::join_adv_ids( tCoordinateDeterminingADVIDs, tFirstParentNode.get_coordinate_determining_adv_ids() );
        }
        if ( tSecondParentNode.depends_on_advs() )
        {
            Intersection_Node::join_adv_ids( tCoordinateDeterminingADVIDs, tSecondParentNode.get_coordinate_determining_adv_ids() );
        }

        // Return joined ADV IDs
        return tCoordinateDeterminingADVIDs;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    Intersection_Node_Surface_Mesh::get_dxi_dcoordinate_first_parent() const
    {
        MORIS_ERROR( false, "Intersection_Node_Surface_Mesh does not yet support intersections on intersections." );
        return { { 0 } };
    }

    //--------------------------------------------------------------------------------------------------------------


    Matrix< DDRMat >
    Intersection_Node_Surface_Mesh::get_dxi_dcoordinate_second_parent() const
    {
        MORIS_ERROR( false, "Intersection_Node_Surface_Mesh does not yet support intersections on intersections." );
        return { { 0 } };
    }

}    // namespace moris::gen
