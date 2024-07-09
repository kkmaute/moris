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
        // debug checks to make sure the intersection node is being created in the right spot
        // for ( uint iDimension = 0; iDimension < mInterfaceGeometry.get_dimension(); iDimension++ )
        // {
        //     bool tIsBetween = ( this->get_global_coordinates()( iDimension ) - aFirstParentNode.get_global_coordinates()( iDimension ) ) * ( this->get_global_coordinates()( iDimension ) - aSecondParentNode.get_global_coordinates()( iDimension ) ) < 1e-6;
        //     if ( not tIsBetween )
        //     {
        //         PRINT( this->get_global_coordinates() );
        //         PRINT( aFirstParentNode.get_global_coordinates() );
        //         PRINT( aSecondParentNode.get_global_coordinates() );
        //         std::cout << "First parent node index: " << aFirstParentNode.get_index() << std::endl;
        //         std::cout << "Second parent node index: " << aSecondParentNode.get_index() << std::endl;
        //     }

        //     MORIS_ASSERT( tIsBetween, "Intersection node coordinate %d is not between the parent nodes' coordinates.", iDimension );
        // }
    }

    Matrix< DDRMat >
    Intersection_Node_Surface_Mesh::get_rotation_matrix() const
    {
        // get parent vector
        Matrix< DDRMat > tParentVector = trans( this->get_second_parent_node().get_global_coordinates() - this->get_first_parent_node().get_global_coordinates() );

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
        Matrix< IndexMat > tFacetVertexIndices = mParentFacet->get_vertex_inds();

        // Return true if either parent node depends on advs, or any of the parent facet's vertices depend on advs
        return this->get_first_parent_node().depends_on_advs()
            or this->get_second_parent_node().depends_on_advs()
            or std::any_of( tFacetVertexIndices.cbegin(),
                    tFacetVertexIndices.cend(),
                    [ & ]( uint iParentFacetVertexIndex ) { return mInterfaceGeometry.facet_vertex_depends_on_advs( iParentFacetVertexIndex ); } );
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    Intersection_Node_Surface_Mesh::compute_dxi_dfacet() const
    {
        // Get the parent vector and its norm from the intersection node
        Matrix< DDRMat > tParentVector     = this->get_second_parent_node().get_global_coordinates() - this->get_first_parent_node().get_global_coordinates();
        real             tParentVectorNorm = norm( tParentVector );

        // Get the rotation matrix from the intersection node
        Matrix< DDRMat > tRotationMatrix = this->get_rotation_matrix();

        // Get the facet vertices
        Matrix< DDRMat > tVertexCoordinates = mParentFacet->get_vertex_coords();

        // get the normal vector rotated in the local coordinate frame
        Matrix< DDRMat > tNormalPrime = tRotationMatrix * mParentFacet->get_normal();

        // get the center vector in the local coordinate frame
        Matrix< DDRMat > tCenterPrime = 2.0 / tParentVectorNorm * tRotationMatrix * ( mParentFacet->get_center() - trans( this->get_first_parent_node().get_global_coordinates() ) );

        // get the jacobian of the center vector wrt to the facet vertices (same for all vertices)
        Matrix< DDRMat > tdCenterdVertices = 2.0 / ( (real)mInterfaceGeometry.get_dimension() * tParentVectorNorm ) * tRotationMatrix;

        // get the jacobians of the normal vector wrt to the facet vertices
        Vector< Matrix< DDRMat > > tNormalVectorSensitivities( mInterfaceGeometry.get_dimension() );
        switch ( mInterfaceGeometry.get_dimension() )
        {
            case 2:    // 2D surface mesh
            {
                // magnitude of the normal vector
                real tNormalVectorNorm = norm( tVertexCoordinates.get_row( 1 ) - tVertexCoordinates.get_row( 0 ) );

                tNormalVectorSensitivities( 0 ) = { { ( tVertexCoordinates( 0, 0 ) - tVertexCoordinates( 1, 0 ) ) * ( tVertexCoordinates( 1, 1 ) - tVertexCoordinates( 0, 1 ) ), std::pow( tVertexCoordinates( 1, 0 ) - tVertexCoordinates( 0, 0 ), 2.0 ) }, { -1.0 * std::pow( tVertexCoordinates( 1, 1 ) - tVertexCoordinates( 0, 1 ), 2.0 ), ( tVertexCoordinates( 0, 0 ) - tVertexCoordinates( 1, 0 ) ) * ( tVertexCoordinates( 0, 1 ) - tVertexCoordinates( 1, 1 ) ) } };
                tNormalVectorSensitivities( 1 ) = { { ( tVertexCoordinates( 0, 1 ) - tVertexCoordinates( 1, 1 ) ) * ( tVertexCoordinates( 0, 0 ) - tVertexCoordinates( 1, 0 ) ), -1.0 * std::pow( tVertexCoordinates( 1, 0 ) - tVertexCoordinates( 0, 0 ), 2.0 ) }, { std::pow( tVertexCoordinates( 1, 1 ) - tVertexCoordinates( 0, 1 ), 2.0 ), ( tVertexCoordinates( 1, 0 ) - tVertexCoordinates( 0, 0 ) ) * ( tVertexCoordinates( 0, 1 ) - tVertexCoordinates( 1, 1 ) ) } };
                tNormalVectorSensitivities( 0 ) = tRotationMatrix / std::pow( tNormalVectorNorm, 3.0 ) * tNormalVectorSensitivities( 0 );
                tNormalVectorSensitivities( 1 ) = tRotationMatrix / std::pow( tNormalVectorNorm, 3.0 ) * tNormalVectorSensitivities( 1 );
                break;
            }
            case 3:    // 3D surface mesh
            {
                Vector< Matrix< DDRMat > > tNormalVectorNormSensitivity( mInterfaceGeometry.get_dimension() );

                // Compute the normal vector (not unit)
                Matrix< DDRMat > tNormal = cross( tVertexCoordinates.get_row( 1 ) - tVertexCoordinates.get_row( 0 ), tVertexCoordinates.get_row( 2 ) - tVertexCoordinates.get_row( 0 ) );

                // magnitude of the normal vector
                real tNormalVectorNorm = norm( tNormal );

                // jacobians of the normal vector wrt to the facet vertices
                tNormalVectorSensitivities( 0 ) = { { 0.0, tVertexCoordinates( 1, 2 ) - tVertexCoordinates( 2, 2 ), tVertexCoordinates( 2, 1 ) - tVertexCoordinates( 1, 1 ) }, { tVertexCoordinates( 2, 2 ) - tVertexCoordinates( 1, 2 ), 0.0, tVertexCoordinates( 1, 0 ) - tVertexCoordinates( 2, 0 ) }, { tVertexCoordinates( 1, 1 ) - tVertexCoordinates( 2, 1 ), tVertexCoordinates( 2, 0 ) - tVertexCoordinates( 1, 0 ), 0.0 } };
                tNormalVectorSensitivities( 1 ) = { { 0.0, tVertexCoordinates( 2, 2 ) - tVertexCoordinates( 0, 2 ), tVertexCoordinates( 0, 1 ) - tVertexCoordinates( 2, 1 ) }, { tVertexCoordinates( 0, 2 ) - tVertexCoordinates( 2, 2 ), 0.0, tVertexCoordinates( 2, 0 ) - tVertexCoordinates( 0, 0 ) }, { tVertexCoordinates( 2, 1 ) - tVertexCoordinates( 0, 1 ), tVertexCoordinates( 0, 0 ) - tVertexCoordinates( 2, 0 ), 0.0 } };
                tNormalVectorSensitivities( 2 ) = { { 0.0, tVertexCoordinates( 0, 2 ) - tVertexCoordinates( 1, 2 ), tVertexCoordinates( 1, 1 ) - tVertexCoordinates( 0, 1 ) }, { tVertexCoordinates( 1, 2 ) - tVertexCoordinates( 0, 2 ), 0.0, tVertexCoordinates( 0, 0 ) - tVertexCoordinates( 1, 0 ) }, { tVertexCoordinates( 0, 1 ) - tVertexCoordinates( 1, 1 ), tVertexCoordinates( 1, 0 ) - tVertexCoordinates( 0, 0 ), 0.0 } };

                tNormalVectorNormSensitivity( 0 ) = { { tNormal( 2 ) * ( tVertexCoordinates( 1, 1 ) - tVertexCoordinates( 2, 1 ) ) - tNormal( 1 ) * ( tVertexCoordinates( 1, 2 ) - tVertexCoordinates( 2, 2 ) ), -tNormal( 2 ) * ( tVertexCoordinates( 1, 0 ) - tVertexCoordinates( 2, 0 ) ) + tNormal( 0 ) * ( tVertexCoordinates( 1, 2 ) - tVertexCoordinates( 2, 2 ) ), tNormal( 1 ) * ( tVertexCoordinates( 1, 0 ) - tVertexCoordinates( 2, 0 ) ) - tNormal( 0 ) * ( tVertexCoordinates( 1, 1 ) - tVertexCoordinates( 2, 1 ) ) } };
                tNormalVectorNormSensitivity( 1 ) = { { -tNormal( 2 ) * ( tVertexCoordinates( 0, 1 ) - tVertexCoordinates( 2, 1 ) ) + tNormal( 1 ) * ( tVertexCoordinates( 0, 2 ) - tVertexCoordinates( 2, 2 ) ), tNormal( 2 ) * ( tVertexCoordinates( 0, 0 ) - tVertexCoordinates( 2, 0 ) ) - tNormal( 0 ) * ( tVertexCoordinates( 0, 2 ) - tVertexCoordinates( 2, 2 ) ), -tNormal( 1 ) * ( tVertexCoordinates( 0, 0 ) - tVertexCoordinates( 2, 0 ) ) + tNormal( 0 ) * ( tVertexCoordinates( 0, 1 ) - tVertexCoordinates( 2, 1 ) ) } };
                tNormalVectorNormSensitivity( 2 ) = { { tNormal( 2 ) * ( tVertexCoordinates( 0, 1 ) - tVertexCoordinates( 1, 1 ) ) - tNormal( 1 ) * ( tVertexCoordinates( 0, 1 ) - tVertexCoordinates( 1, 2 ) ), -tNormal( 2 ) * ( tVertexCoordinates( 0, 0 ) - tVertexCoordinates( 1, 0 ) ) + tNormal( 0 ) * ( tVertexCoordinates( 0, 2 ) - tVertexCoordinates( 1, 2 ) ), tNormal( 1 ) * ( tVertexCoordinates( 0, 0 ) - tVertexCoordinates( 1, 0 ) ) - tNormal( 0 ) * ( tVertexCoordinates( 0, 1 ) - tVertexCoordinates( 1, 1 ) ) } };
                tNormalVectorNormSensitivity( 0 ) = 1.0 / tNormalVectorNorm * tNormalVectorNormSensitivity( 0 );
                tNormalVectorNormSensitivity( 1 ) = 1.0 / tNormalVectorNorm * tNormalVectorNormSensitivity( 1 );
                tNormalVectorNormSensitivity( 2 ) = 1.0 / tNormalVectorNorm * tNormalVectorNormSensitivity( 2 );


                for ( uint iDimension = 0; iDimension < 3; iDimension++ )
                {
                    tNormalVectorSensitivities( iDimension ) = tRotationMatrix / tNormalVectorNorm * ( tNormalVectorSensitivities( iDimension ) - mParentFacet->get_normal() / tNormalVectorNorm * tNormalVectorNormSensitivity( iDimension ) );
                }

                break;
            }
        }

        // Compute the local coordinate sensitivity wrt the facet vertices
        Matrix< DDRMat > tdXidFacet( mInterfaceGeometry.get_dimension(), mInterfaceGeometry.get_dimension() );
        for ( uint iDimension = 0; iDimension < mInterfaceGeometry.get_dimension(); iDimension++ )
        {
            tdXidFacet.set_column( iDimension,
                    trans( tNormalVectorSensitivities( iDimension ) ) * tCenterPrime / tNormalPrime( 0 )
                            + trans( tdCenterdVertices ) * tNormalPrime / tNormalPrime( 0 )
                            - dot( tCenterPrime, tNormalPrime ) * trans( tNormalVectorSensitivities( iDimension ).get_row( 0 ) ) / std::pow( tNormalPrime( 0 ), 2.0 ) );
        }

        return tdXidFacet;
    }

    //--------------------------------------------------------------------------------------------------------------

    void Intersection_Node_Surface_Mesh::append_dcoordinate_dadv(
            Matrix< DDRMat >&       aCoordinateSensitivities,
            const Matrix< DDRMat >& aSensitivityFactor ) const
    {
        // Get parent nodes
        const Basis_Node& tFirstParentNode  = this->get_first_parent_node();
        const Basis_Node& tSecondParentNode = this->get_second_parent_node();

        // Compute parent vector
        Matrix< DDRMat > tParentVector = trans( tSecondParentNode.get_global_coordinates() - tFirstParentNode.get_global_coordinates() );

        // Get the sensitivity of the local coordinate wrt the facet vertices
        Matrix< DDRMat > tLocalCoordinateFacetVertexSensitivities = this->compute_dxi_dfacet();

        // Counter for local vertex index
        uint tLocalFacetVertexIndex = 0;

        // Loop over the facet parents
        for ( uint iParentFacetVertexIndex : mParentFacet->get_vertex_inds() )
        {
            if ( mInterfaceGeometry.facet_vertex_depends_on_advs( iParentFacetVertexIndex ) )
            {
                Matrix< DDRMat > tSensitivitiesToAdd = .5 * aSensitivityFactor * tParentVector *    //
                                                       ( trans( tLocalCoordinateFacetVertexSensitivities.get_column( tLocalFacetVertexIndex ) ) * mInterfaceGeometry.get_dvertex_dadv( iParentFacetVertexIndex ) );

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

    Vector< sint >
    Intersection_Node_Surface_Mesh::get_coordinate_determining_adv_ids() const
    {
        // Initialize ADV IDs
        Vector< sint > tCoordinateDeterminingADVIDs;

        // Get ADV IDs from facet parents
        for ( uint tParentFacetVertex : mParentFacet->get_vertex_inds() )
        {
            if ( mInterfaceGeometry.facet_vertex_depends_on_advs( tParentFacetVertex ) )
            {    // Get the IDs for this vertex
                Vector< sint > tVertexADVIds = mInterfaceGeometry.get_vertex_adv_ids( tParentFacetVertex );

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
