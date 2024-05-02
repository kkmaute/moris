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
#include "fn_dot.hpp"

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
            , mParentFacet( aParentFacet )
            , mInterfaceGeometry( aInterfaceGeometry )
    {
    }

    Matrix< DDRMat >
    Intersection_Node_Surface_Mesh::get_rotation_matrix() const
    {
        // get parent vector
        Matrix< DDRMat > tParentVector = this->get_second_parent_node().get_global_coordinates() - this->get_first_parent_node().get_global_coordinates();

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

        // check that the rotation matrix is correct by ensuring the parent vector was rotated to the x axis
        MORIS_ASSERT( norm( tRotationMatrix * tParentVector - tCastAxis ) < mInterfaceGeometry.get_intersection_tolerance(), "Rotation matrix should rotate the parent vector to the x axis." );

        // trim the transformation matrix if 2D
        if ( mInterfaceGeometry.get_dimension() == 2 )
        {
            tRotationMatrix.resize( 2, 2 );
        }

        return tRotationMatrix;
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

    bool Intersection_Node_Surface_Mesh::depends_on_advs() const
    {
        // Check if the facet vertices depend on advs
        for ( auto iParentFacetVertexIndex : mParentFacet->get_vertex_inds() )
        {
            if ( mInterfaceGeometry.facet_vertex_depends_on_advs( iParentFacetVertexIndex ) )
            {
                return true;
            }
        }

        // Get parent nodes
        const Basis_Node& tFirstParentNode  = this->get_first_parent_node();
        const Basis_Node& tSecondParentNode = this->get_second_parent_node();

        return tFirstParentNode.depends_on_advs() or tSecondParentNode.depends_on_advs();
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    Intersection_Node_Surface_Mesh::compute_dxi_dfacet() const
    {
        // Error for 3D problems for now
        MORIS_ERROR( mInterfaceGeometry.get_dimension() < 3, "Intersection_Node_Surface_Mesh does not yet support 3D problems. " );

        // Get the parent vector and its norm from the intersection node
        Matrix< DDRMat > tParentVector     = this->get_second_parent_node().get_global_coordinates() - this->get_first_parent_node().get_global_coordinates();
        real             tParentVectorNorm = norm( tParentVector );

        // Get the rotation matrix from the intersection node
        Matrix< DDRMat > tRotationMatrix = this->get_rotation_matrix();

        // Get the facet vertices
        Matrix< DDRMat > tVertex1Coordinates = mParentFacet->get_vertex_coords().get_row( 0 );
        Matrix< DDRMat > tVertex2Coordinates = mParentFacet->get_vertex_coords().get_row( 1 );

        // get the normal vector rotated in the local coordinate frame
        Matrix< DDRMat > tNormalPrime = tRotationMatrix * mParentFacet->get_normal();

        // get the center vector in the local coordinate frame
        Matrix< DDRMat > tCenterPrime = 2.0 / ( (real)mInterfaceGeometry.get_dimension() * tParentVectorNorm ) * tRotationMatrix * ( mParentFacet->get_center() - trans( this->get_first_parent_node().get_global_coordinates() ) );

        // Sensitivity of the normal vector to the vertices

        Matrix< DDRMat > tdNormalPrimedVertex1 = { { ( tVertex2Coordinates( 1 ) - tVertex1Coordinates( 1 ) ) * ( tVertex1Coordinates( 0 ) - tVertex2Coordinates( 0 ) ), -1.0 * std::pow( tVertex2Coordinates( 1 ) - tVertex1Coordinates( 1 ), 2.0 ) }, { std::pow( tVertex2Coordinates( 0 ) - tVertex1Coordinates( 0 ), 2.0 ), ( tVertex1Coordinates( 0 ) - tVertex2Coordinates( 0 ) ) * ( tVertex1Coordinates( 1 ) - tVertex2Coordinates( 1 ) ) } };
        tdNormalPrimedVertex1                  = tRotationMatrix / std::pow( norm( tVertex2Coordinates - tVertex1Coordinates ), 3.0 ) * tdNormalPrimedVertex1;
        Matrix< DDRMat > tdNormalPrimedVertex2 = { { ( tVertex1Coordinates( 1 ) - tVertex2Coordinates( 1 ) ) * ( tVertex1Coordinates( 0 ) - tVertex2Coordinates( 0 ) ), std::pow( tVertex2Coordinates( 1 ) - tVertex1Coordinates( 1 ), 2.0 ) }, { -1.0 * std::pow( tVertex2Coordinates( 0 ) - tVertex1Coordinates( 0 ), 2.0 ), ( tVertex2Coordinates( 0 ) - tVertex1Coordinates( 0 ) ) * ( tVertex1Coordinates( 1 ) - tVertex2Coordinates( 1 ) ) } };
        tdNormalPrimedVertex2                  = tRotationMatrix / std::pow( norm( tVertex2Coordinates - tVertex1Coordinates ), 3.0 ) * tdNormalPrimedVertex2;

        Matrix< DDRMat > tdCenterdVertices = 2.0 / ( (real)mInterfaceGeometry.get_dimension() * tParentVectorNorm ) * tRotationMatrix;

        return join_rows( tdNormalPrimedVertex1 * tCenterPrime / tNormalPrime( 0 ) + tdCenterdVertices * tNormalPrime / tNormalPrime( 0 ) - dot( tCenterPrime, tNormalPrime ) * trans( tdNormalPrimedVertex1.get_row( 0 ) ) / std::pow( tNormalPrime( 0 ), 2.0 ),
                tdNormalPrimedVertex2 * tCenterPrime / tNormalPrime( 0 ) + tdCenterdVertices * tNormalPrime / tNormalPrime( 0 ) - dot( tCenterPrime, tNormalPrime ) * trans( tdNormalPrimedVertex2.get_row( 0 ) ) / std::pow( tNormalPrime( 0 ), 2.0 ) );
    }

    //--------------------------------------------------------------------------------------------------------------

    void Intersection_Node_Surface_Mesh::append_dcoordinate_dadv( Matrix< DDRMat >& aCoordinateSensitivities, const Matrix< DDRMat >& aSensitivityFactor ) const
    {
        // Get parent nodes
        const Basis_Node& tFirstParentNode  = this->get_first_parent_node();
        const Basis_Node& tSecondParentNode = this->get_second_parent_node();

        // BRENDAN
        if ( mParentFacet->get_index() == 0 and std::abs( tFirstParentNode.get_global_coordinates()( 0 ) - 1.0 ) < 1e-4 and std::abs( tSecondParentNode.get_global_coordinates()( 0 ) - 1.0 ) < 1e-4 )
        {
            std::cout << "this is the node we are testing\n";
        }

        // BRENDAN
        PRINT( tFirstParentNode.get_global_coordinates() );
        PRINT( tSecondParentNode.get_global_coordinates() );

        // Compute parent vector
        Matrix< DDRMat > tParentVector = trans( tSecondParentNode.get_global_coordinates() - tFirstParentNode.get_global_coordinates() );

        // Get the sensitivity of the local coordinate wrt the facet vertices
        Matrix< DDRMat > tLocalCoordinateFacetVertexSensitivities = this->compute_dxi_dfacet();

        Matrix< DDRMat > tDenseSensitivities;
        // Loop over the facet parents
        for ( uint iParentFacetVertexIndex : mParentFacet->get_vertex_inds() )
        {
            // Counter for local vertex index
            uint tLocalFacetVertexIndex = 0;
            if ( mInterfaceGeometry.facet_vertex_depends_on_advs( iParentFacetVertexIndex ) )
            {
                Matrix< DDRMat > tDenseSensitivities = .5 * aSensitivityFactor * tParentVector * ( trans( tLocalCoordinateFacetVertexSensitivities.get_column( tLocalFacetVertexIndex ) ) * mInterfaceGeometry.get_dvertex_dadv( iParentFacetVertexIndex ) );

                // Expand the sensitivity matrix to be per ADV ID
                Matrix< DDRMat > tSensitivitiesToAdd( tDenseSensitivities.n_rows(), tDenseSensitivities.n_cols() * tDenseSensitivities.n_rows() );
                for ( uint iDimensionIndex = 0; iDimensionIndex < tDenseSensitivities.n_rows(); iDimensionIndex++ )
                {
                    for ( uint iADVIndex = 0; iADVIndex < tDenseSensitivities.n_cols(); iADVIndex++ )
                    {
                        tSensitivitiesToAdd( iDimensionIndex, tDenseSensitivities.n_rows() * iADVIndex + iDimensionIndex ) = tDenseSensitivities( iDimensionIndex, iADVIndex );
                    }
                }

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
            }
            // Increment local vertex index
            tLocalFacetVertexIndex++;
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
        // BRENDAN
        if ( std::abs( this->get_local_coordinate() + 0.5 ) < 1e-4 or std::abs( this->get_local_coordinate() - 1.5 ) < 1e-4 )
        {
            std::cout << "this is the node we are testing\n";
        }

        // Initialize ADV IDs
        Matrix< DDSMat > tCoordinateDeterminingADVIDs;

        // Get ADV IDs from facet parents
        for ( uint tParentFacetVertex : mParentFacet->get_vertex_inds() )
        {
            if ( mInterfaceGeometry.facet_vertex_depends_on_advs( tParentFacetVertex ) )
            {    // Get the IDs for this vertex
                Matrix< DDSMat > tVertexADVIds = mInterfaceGeometry.get_vertex_adv_ids( tParentFacetVertex );

                // Join IDs
                Intersection_Node::join_adv_ids( tCoordinateDeterminingADVIDs, tVertexADVIds );
            }
        }

        // Get parent nodes
        const Basis_Node& tFirstParentNode  = this->get_first_parent_node();
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
