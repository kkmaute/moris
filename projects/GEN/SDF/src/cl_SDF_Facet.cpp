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

#include "cl_Vector.hpp"

#include "assert.hpp"
#include "fn_dot.hpp"
#include "fn_norm.hpp"
#include "op_times.hpp"
#include "fn_stringify_matrix.hpp"
#include "op_equal_equal.hpp"
#include "fn_all_true.hpp"

namespace moris
{
    namespace sdf
    {
        //-------------------------------------------------------------------------------

        Facet::Facet(
                moris_index                                aIndex,
                Vector< std::shared_ptr< Facet_Vertex > >& aVertices,
                uint                                       aDimension,
                real                                       aIntersectionTolerance )
                : mIndex( aIndex )
                , mVertices( aVertices )
                , mCenter( aDimension, 1 )
                , mNormal( aDimension, 1 )
                , mMinCoord( aDimension, 1 )
                , mMaxCoord( aDimension, 1 )
                , mIntersectionTolerance( aIntersectionTolerance )
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
            if ( std::abs( mNormal( aAxis ) ) < mIntersectionTolerance )
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

        real
        Facet::get_vertex_coord( uint aVertexIndex, uint aAxis )
        {
            return mVertices( aVertexIndex )->get_coord( aAxis );
        }

        //-------------------------------------------------------------------------------
        // MTK Interface
        //-------------------------------------------------------------------------------
        uint
        Facet::get_number_of_vertices() const
        {
            return mVertices.size();
        }


        Vector< mtk::Vertex* >
        Facet::get_vertex_pointers() const
        {
            uint                   tDimension = get_number_of_vertices();
            Vector< mtk::Vertex* > tVertices( tDimension, nullptr );

            for ( uint k = 0; k < tDimension; ++k )
            {
                tVertices( k ) = mVertices( k ).get();
            }
            return tVertices;
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
            uint             tDimension = this->get_number_of_vertices();
            Matrix< DDRMat > tCoords( tDimension, tDimension );

            for ( uint iVertex = 0; iVertex < tDimension; ++iVertex )
            {
                tCoords.set_row( iVertex, mVertices( iVertex )->get_coords() );
            }

            return tCoords;
        }

        //-------------------------------------------------------------------------------

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
        Facet::compute_center()
        {
            uint tDimension = mVertices.size();

            // reset center
            mCenter.fill( 0 );

            // loop over all nodes
            for ( uint iVertexNumber = 0; iVertexNumber < tDimension; ++iVertexNumber )
            {
                // add to the center total
                for ( uint iAxis = 0; iAxis < tDimension; ++iAxis )
                {
                    mCenter( iAxis ) += mVertices( iVertexNumber )->get_coord( iAxis );
                }
            }

            // divide center by dimension
            for ( uint iAxis = 0; iAxis < tDimension; ++iAxis )
            {
                mCenter( iAxis ) /= (real)tDimension;
            }
        }

        void
        Facet::compute_min_and_max_coordinates()
        {
            uint tDimension = mVertices.size();

            // identify minimum and maximum coordinate
            for ( uint iAxis = 0; iAxis < tDimension; ++iAxis )
            {
                // FIXME: there's probably an easier way to determine the min and max without the use of a switch case
                switch ( tDimension )
                {
                    case 2:
                    {
                        mMinCoord( iAxis ) = min( mVertices( 0 )->get_coord( iAxis ),
                                mVertices( 1 )->get_coord( iAxis ) );

                        mMaxCoord( iAxis ) = max( mVertices( 0 )->get_coord( iAxis ),
                                mVertices( 1 )->get_coord( iAxis ) );
                        break;
                    }
                    case 3:
                    {
                        mMinCoord( iAxis ) = min( mVertices( 0 )->get_coord( iAxis ),
                                mVertices( 1 )->get_coord( iAxis ),
                                mVertices( 2 )->get_coord( iAxis ) );

                        mMaxCoord( iAxis ) = max( mVertices( 0 )->get_coord( iAxis ),
                                mVertices( 1 )->get_coord( iAxis ),
                                mVertices( 2 )->get_coord( iAxis ) );
                        break;
                    }
                    default:
                    {
                        MORIS_ASSERT( false, "SDF Facet() - mMinCoord and mMaxCoord not properly computed for dimension %d", tDimension );
                    }
                }
            }
        }

        bool
        Facet::operator==( const Facet& aRHS ) const
        {
            return mVertices.size() != aRHS.get_number_of_vertices()
                && all_true( this->get_vertex_ids() == aRHS.get_vertex_ids() )
                && all_true( mNormal == aRHS.get_normal() )
                && std::abs( mHesse - aRHS.get_hesse() ) < mIntersectionTolerance;
        }
    } /* namespace sdf */
} /* namespace moris */
