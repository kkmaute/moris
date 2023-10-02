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
            moris_index                   aIndex,
            moris::Cell< Facet_Vertex* >& aVertices )
            : Facet( aIndex, aVertices, 2 )
    {
        this->update_data();
    }

    //-------------------------------------------------------------------------------

    void
    Line::update_data()
    {
        // step 1: copy node coordinates and determine center
        this->copy_node_coords_and_inds( mVertices, 2 );

        this->calculate_hesse_normal_form();
    }

    void
    Line::calculate_hesse_normal_form()
    {
        // Normal vector defined by (-dy, dx)
        mNormal( 0 ) = -( mNodeCoords( 1, 1 ) - mNodeCoords( 1, 0 ) );
        mNormal( 1 ) = mNodeCoords( 0, 1 ) - mNodeCoords( 0, 0 );

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
        // BRENDAN TODO - MAYBE?
        return false;
    }

    //-------------------------------------------------------------------------------

    real
    Line::get_distance_to_point(
            const Matrix< DDRMat >& aPoint )
    {
        // BRENDAN TODO - MAYBE?
        return -1.0;
    }

} /* namespace moris::sdf */
