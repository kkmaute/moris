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
        MORIS_ERROR( false, "SDF_Line()::check_edge() not implemented for Lines yet." );
        return false;
    }

    //-------------------------------------------------------------------------------

    real
    Line::get_distance_to_point(
            const Matrix< DDRMat >& aPoint )
    {
        MORIS_ERROR( false, "SDF_Line()::get_distance_to_point() not implemented for Lines yet." );
        return -1.0;
    }

} /* namespace moris::sdf */
