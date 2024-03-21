/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SDF_Line.cpp
 *
 */

#include "cl_SDF_Line.hpp"

#include "assert.hpp"
#include "fn_dot.hpp"
#include "fn_norm.hpp"
#include "fn_trans.hpp"
#include "op_times.hpp"
#include "SDF_Tools.hpp"
#include "fn_stringify_matrix.hpp"

namespace moris::sdf
{
    //-------------------------------------------------------------------------------

    Line::Line(
            moris_index                                aIndex,
            Vector< std::shared_ptr< Facet_Vertex > >& aVertices )
            : Facet( aIndex, aVertices, 2 )
    {
        this->update_data();
    }

    //-------------------------------------------------------------------------------

    void
    Line::update_data()
    {
        // step 1: compute center
        this->compute_center();

        // step 2: compute the bounding box
        this->compute_min_and_max_coordinates();

        // compute the normal and hesse
        this->calculate_hesse_normal_form();
    }

    void
    Line::calculate_hesse_normal_form()
    {
        // Normal vector defined by (-dy, dx)
        mNormal( 0 ) = -( mVertices( 1 )->get_coord( 1 ) - mVertices( 0 )->get_coord( 1 ) );
        mNormal( 1 ) = mVertices( 1 )->get_coord( 0 ) - mVertices( 0 )->get_coord( 0 );

        // Make normal a unit vector
        real tNorm = norm( mNormal );
        for ( moris_index iDim = 0; iDim < 2; iDim++ )
        {
            mNormal( iDim ) /= tNorm;
        }

        // compute Hesse distance
        mHesse = dot( mNormal, mCenter );
    }

    //-------------------------------------------------------------------------------
    // SDF functions
    //-------------------------------------------------------------------------------

    bool
    Line::check_edge(
            const uint              aEdge,
            const uint              aAxis,
            const Matrix< DDRMat >& aPoint )
    {
        MORIS_ERROR( false, "SDF_Line()::check_edge() not implemented for Line." );
        return false;
    }

    //-------------------------------------------------------------------------------

    real
    Line::get_distance_to_point(
            const Matrix< DDRMat >& aPoint )
    {
        MORIS_ERROR( false, "SDF_Line()::get_distance_to_point() not implemented for Line." );
        return -1.0;
    }

    //-------------------------------------------------------------------------------

    Matrix< DDRMat >
    Line::compute_dxi_dvertices( moris::gen::Intersection_Node_Surface_Mesh& aIntersectionNode )
    {
        // Get the parent vector and its norm from the intersection node
        Matrix< DDRMat > tParentVector     = aIntersectionNode.get_first_parent_node().get_global_coordinates() - aIntersectionNode.get_second_parent_node().get_global_coordinates();
        real             tParentVectorNorm = norm( tParentVector );

        // Get the rotation matrix from the intersection node
        Matrix< DDRMat > tRotationMatrix;
        aIntersectionNode.get_rotation_matrix( tRotationMatrix );

        // get the normal vector rotated in the local coordinate frame
        Matrix< DDRMat > tNormalPrime = tRotationMatrix * mNormal;

        // get the center vector in the local coordinate frame
        Matrix< DDRMat > tCenterPrime = 2.0 / ( 3.0 * tParentVectorNorm ) * tRotationMatrix * mCenter;

        // compute vector between facet vertices, and unnormalized normal vector
        Matrix< DDRMat > tFacetVector = mVertices( 1 )->get_coords() - mVertices( 0 )->get_coords();
        Matrix< DDRMat > tNormal      = { { mVertices( 0 )->get_coord( 1 ) - mVertices( 1 )->get_coord( 1 ), mVertices( 1 )->get_coord( 0 ) - mVertices( 0 )->get_coord( 0 ) } };

        // derivative of normal vector
        Matrix< DDRMat > tdNormaldVertex1 = { { 1.0, -1.0 } };
        Matrix< DDRMat > tdNormaldVertex2 = { { -1.0, 1.0 } };

        // Sensitivity of the normal vector to the vertices
        Matrix< DDRMat > tdNormalPrimedVertex1 = tRotationMatrix * ( 1.0 / norm( tFacetVector ) * tdNormaldVertex1 - 0.5 * tNormal * ( trans( tFacetVector ) * 1.0 * eye( 2, 2 ) - eye( 2, 2 ) * ( tFacetVector ) ) );
        Matrix< DDRMat > tdNormalPrimedVertex2 = tRotationMatrix * ( 1.0 / norm( tFacetVector ) * tdNormaldVertex2 - 0.5 * tNormal * ( trans( tFacetVector ) * 1.0 * eye( 2, 2 ) - eye( 2, 2 ) * ( tFacetVector ) ) );


        Matrix< DDRMat > tdCenterdVertices = 2.0 / ( 3.0 * tParentVectorNorm ) * tRotationMatrix;

        return join_cols( ( tdNormalPrimedVertex1 * tCenterPrime + tdCenterdVertices * tNormalPrime ) * mNormal( 0 ) - tdNormalPrimedVertex1( 0 ) * ( tNormalPrime * tCenterPrime ),
                ( tdNormalPrimedVertex2 * tCenterPrime + tdCenterdVertices * tNormalPrime ) * mNormal( 0 ) - tdNormalPrimedVertex2( 0 ) * ( tNormalPrime * tCenterPrime ) );
    }

} /* namespace moris::sdf */
