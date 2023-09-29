/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SDF_Facet.cpp
 *
 */

#include "cl_SDF_Facet.hpp"

#include "cl_Cell.hpp"

#include "assert.hpp"
#include "fn_dot.hpp"
#include "fn_norm.hpp"
#include "op_times.hpp"
#include "SDF_Tools.hpp"
#include "fn_stringify_matrix.hpp"

namespace moris
{
    namespace sdf
    {
        //-------------------------------------------------------------------------------

        Facet::Facet(
                moris_index                   aIndex,
                moris::Cell< Facet_Vertex* >& aVertices,
                uint                          aDimension )
                : mIndex( aIndex )
                , mVertices( aVertices )
                , mNodeCoords( aDimension, aDimension )
                , mNodeIndices( aDimension, 1 )
                , mCenter( aDimension, 1 )
                , mNormal( aDimension, 1 )
                , mMinCoord( aDimension, 1 )
                , mMaxCoord( aDimension, 1 )
        {
        }

        //-------------------------------------------------------------------------------
        // SDF Functions
        //-------------------------------------------------------------------------------

        void
        Facet::intersect_with_coordinate_axis(
                const Matrix< DDRMat >& aPoint,
                const uint              aAxis,
                real&                   aCoordinate,
                bool&                   aError )
        {
            if ( std::abs( mNormal( aAxis ) ) < gSDFepsilon )
            {
                aCoordinate = 0;
                aError      = true;
            }
            else
            {
                aCoordinate = aPoint( aAxis ) + ( mHesse - dot( mNormal, aPoint ) ) / mNormal( aAxis );
                aError      = false;
            }
        }

        //-------------------------------------------------------------------------------
        // MTK Interface
        //-------------------------------------------------------------------------------
        uint
        Facet::get_number_of_vertices() const
        {
            return mVertices.size();
        }


        Cell< mtk::Vertex* >
        Facet::get_vertex_pointers() const
        {
            uint                        tDimension = get_number_of_vertices();
            moris::Cell< mtk::Vertex* > aVertices( tDimension, nullptr );

            for ( uint k = 0; k < tDimension; ++k )
            {
                aVertices( k ) = mVertices( k );
            }
            return aVertices;
        }

        //-------------------------------------------------------------------------------

        // TODO MESH-CLEANUP
        void
        Facet::remove_vertex_pointer( moris_index aIndex )
        {
            // std::cout<<"In SDF facet"<<std::endl;
        }

        //-------------------------------------------------------------------------------

        Matrix< IdMat >
        Facet::get_vertex_ids() const
        {
            uint            tDimension = get_number_of_vertices();
            Matrix< IdMat > aIDs( tDimension, 1 );
            for ( uint k = 0; k < tDimension; ++k )
            {
                aIDs( k ) = mVertices( k )->get_id();
            }

            return aIDs;
        }

        //-------------------------------------------------------------------------------

        Matrix< IndexMat >
        Facet::get_vertex_inds() const
        {
            uint               tDimenion = get_number_of_vertices();
            Matrix< IndexMat > aINDs( tDimenion, 1 );

            for ( uint k = 0; k < tDimenion; ++k )
            {
                aINDs( k ) = mVertices( k )->get_index();
            }
            return aINDs;
        }

        //-------------------------------------------------------------------------------

        Matrix< DDRMat >
        Facet::get_vertex_coords() const
        {
            uint             tDimension = get_number_of_vertices();
            Matrix< DDRMat > aCoords( tDimension, tDimension );

            for ( uint k = 0; k < tDimension; ++k )
            {
                aCoords.set_row( k, mVertices( k )->get_coords() );
            }

            return aCoords;
        }

        mtk::Geometry_Type
        Facet::get_geometry_type() const
        {
            switch ( get_number_of_vertices() )
            {
                case 2:
                {
                    return mtk::Geometry_Type::LINE;
                }
                case 3:
                {
                    return mtk::Geometry_Type::TRI;
                }
                default:
                {
                    MORIS_ERROR( false, "Geometry type not implemented for %d vertex facets.", get_number_of_vertices() );
                    return mtk::Geometry_Type::UNDEFINED;
                }
            }
        }

        void
        Facet::copy_node_coords_and_inds(
                moris::Cell< Facet_Vertex* >& aVertices,
                uint                          aDimension )
        {
            // make sure that the length is correct
            MORIS_ASSERT( aVertices.size() >= aDimension,
                    "Facet() - Dimension of %d and %lu vertices. Number of vertices needs to be at least equal to the number of dimensions.",
                    aDimension,
                    aVertices.size() );

            // step 1: copy node coordinates into member variables
            //         and calculate center

            // reset center
            mCenter.fill( 0 );

            // loop over all nodes
            for ( uint iVertexNumber = 0; iVertexNumber < aDimension; ++iVertexNumber )
            {
                // get vertex coordinates
                auto tNodeCoords = aVertices( iVertexNumber )->get_coords();

                // copy coordinates into member matrix
                for ( uint iAxis = 0; iAxis < aDimension; ++iAxis )
                {
                    mNodeCoords( iAxis, iVertexNumber ) = tNodeCoords( iAxis );
                    mCenter( iAxis ) += tNodeCoords( iAxis );
                }

                // remember node indices
                mNodeIndices( iVertexNumber ) = aVertices( iVertexNumber )->get_index();
            }

            // identify minimum and maximum coordinate
            for ( uint iAxis = 0; iAxis < aDimension; ++iAxis )
            {
                // FIXME: there's probably an easier way to determine the min and max without the use of a switch case
                switch ( aDimension )
                {
                    case 2:
                    {
                        mMinCoord( iAxis ) = min( mNodeCoords( iAxis, 0 ),
                                mNodeCoords( iAxis, 1 ) );

                        mMaxCoord( iAxis ) = max( mNodeCoords( iAxis, 0 ),
                                mNodeCoords( iAxis, 1 ) );
                        break;
                    }
                    case 3:
                    {
                        mMinCoord( iAxis ) = min( mNodeCoords( iAxis, 0 ),
                                mNodeCoords( iAxis, 1 ),
                                mNodeCoords( iAxis, 2 ) );

                        mMaxCoord( iAxis ) = max( mNodeCoords( iAxis, 0 ),
                                mNodeCoords( iAxis, 1 ),
                                mNodeCoords( iAxis, 2 ) );
                        break;
                    }
                        default:
                        {
                            MORIS_ASSERT( false, "SDF Facet() - mMinCoord and mMaxCoord not properly computed for dimension %d", aDimension );
                        }
                }
            }

            // divide center by dimension
            for ( uint iAxis = 0; iAxis < aDimension; ++iAxis )
            {
                mCenter( iAxis ) /= (real)aDimension;
            }
        }
    } /* namespace sdf */
} /* namespace moris */
