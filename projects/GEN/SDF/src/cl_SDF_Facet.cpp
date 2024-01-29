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
#include "fn_stringify_matrix.hpp"
#include "op_equal_equal.hpp"
#include "fn_all_true.hpp"

namespace moris
{
    namespace sdf
    {
        //-------------------------------------------------------------------------------

        Facet::Facet(
                moris_index                                     aIndex,
                moris::Cell< std::shared_ptr< Facet_Vertex > >& aVertices,
                uint                                            aDimension,
                real                                            aIntersectionTolerance )
                : mIndex( aIndex )
                , mVertices( aVertices )
                , mCenter( aDimension, 1 )
                , mNormal( aDimension, 1 )
                , mMinCoord( aDimension, 1 )
                , mMaxCoord( aDimension, 1 )
                , mIntersectionTolerance( aIntersectionTolerance )
        {
        }

        void
        Facet::rotate( const Matrix< DDRMat >& aRotationMatrix )
        {
            for ( uint iVertex = 0; iVertex < mVertices.size(); iVertex++ )
            {
                if ( not mVertices( iVertex )->is_transformed() )
                {
                    mVertices( iVertex )->rotate_node_coords( aRotationMatrix );
                }
            }
        }

        //-------------------------------------------------------------------------------

        void
        Facet::scale( const moris::Cell< real >& aScaling )
        {
            for ( uint iVertex = 0; iVertex < mVertices.size(); iVertex++ )
            {
                if ( not mVertices( iVertex )->is_transformed() )
                {
                    mVertices( iVertex )->scale_node_coords( aScaling );
                }
            }
        }

        //-------------------------------------------------------------------------------

        void
        Facet::shift( const moris::Cell< real >& aShift )
        {
            for ( uint iVertex = 0; iVertex < mVertices.size(); iVertex++ )
            {
                if ( not mVertices( iVertex )->is_transformed() )
                {
                    mVertices( iVertex )->shift_node_coords( aShift );
                }
            }
        }

        //-------------------------------------------------------------------------------

        void
        Facet::reset_vertex_transformed_flags()
        {
            for ( uint iVertex = 0; iVertex < mVertices.size(); iVertex++ )
            {
                mVertices( iVertex )->reset_transformed_flag();
            }
        }

        //-------------------------------------------------------------------------------

        void
        Facet::reset_coordinates()
        {
            for ( uint iVertex = 0; iVertex < mVertices.size(); iVertex++ )
            {
                mVertices( iVertex )->reset_node_coords();
            }
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
        // MTK Interface
        //-------------------------------------------------------------------------------
        uint
        Facet::get_number_of_vertices() const
        {
            return mVertices.size();
        }


        moris::Cell< mtk::Vertex* >
        Facet::get_vertex_pointers() const
        {
            uint                        tDimension = get_number_of_vertices();
            moris::Cell< mtk::Vertex* > tVertices( tDimension, nullptr );

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
            uint             tDimension = get_number_of_vertices();
            Matrix< DDRMat > aCoords( tDimension, tDimension );

            for ( uint iVertex = 0; iVertex < tDimension; ++iVertex )
            {
                for ( uint iDimensionIndex = 0; iDimensionIndex < tDimension; iDimensionIndex++ )
                    aCoords( iVertex, iDimensionIndex ) = mVertices( iVertex )->get_coord( iDimensionIndex );
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
            if ( mVertices.size() != aRHS.get_number_of_vertices() )
            {
                return false;
            }

            return all_true( this->get_vertex_ids() == aRHS.get_vertex_ids() ) && all_true( mNormal == aRHS.get_normal() ) && std::abs( mHesse - aRHS.get_hesse() ) < mIntersectionTolerance;
        }
    } /* namespace sdf */
} /* namespace moris */
