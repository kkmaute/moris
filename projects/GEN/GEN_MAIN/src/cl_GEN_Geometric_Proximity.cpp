/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Geometric_Proximity.cpp
 *
 */

#include "cl_GEN_Geometric_Proximity.hpp"

namespace moris
{
    namespace ge
    {
        //-----------------------------------------------------------------------------------------

        Geometric_Proximity::Geometric_Proximity()
        {
            MORIS_ERROR( false, "ge::Geometric_Proximity::Geometric_Proximity() - Default constructor not implemented." );
        }

        //-----------------------------------------------------------------------------------------

        Geometric_Proximity::Geometric_Proximity(
                moris_index const &aNumGeometries )
                : mAssociatedVertexIndex( MORIS_INDEX_MAX )
                , mGeometricProximity( aNumGeometries, MORIS_INDEX_MAX )
        {
            // do nothing else
        }

        //-----------------------------------------------------------------------------------------

        Geometric_Proximity::~Geometric_Proximity()
        {
            // do nothing
        }

        //-----------------------------------------------------------------------------------------

        void
        Geometric_Proximity::set_num_geometries( const uint aNumGeometries )
        {
            // check that this function is only used if the number of geometries is not initialized yet
            MORIS_ASSERT( mGeometricProximity.size() == 0, "Geometric_Proximity::set_num_geometries() - Number of geometries already set." );

            // resize the vector
            mGeometricProximity.resize( aNumGeometries, MORIS_INDEX_MAX );
        }

        //-----------------------------------------------------------------------------------------

        void
        Geometric_Proximity::set_geometric_proximity(
                const moris_index aGeometricProximity,
                const moris_index aGeometryIndex )
        {
            MORIS_ASSERT( aGeometryIndex < (moris_index)mGeometricProximity.size(),
                    "GEN::Geometric_Proximity::set_geometric_proximity() - Geometry index out of bounds" );

            // make sure the proximity is not already set
            MORIS_ERROR( mGeometricProximity( aGeometryIndex ) == MORIS_INDEX_MAX,
                    "GEN::Geometric_Proximity::set_geometric_proximity() - Geometric proximity already set." );

            mGeometricProximity( aGeometryIndex ) = aGeometricProximity;
        }

        //-----------------------------------------------------------------------------------------

        bool
        Geometric_Proximity::is_geometric_proximity_set( const moris_index aGeometricProximity ) const
        {
            return ( mGeometricProximity( aGeometricProximity ) != MORIS_INDEX_MAX );
        }

        //-----------------------------------------------------------------------------------------

        moris_index
        Geometric_Proximity::get_geometric_proximity( const moris_index aGeometryIndex ) const
        {
            moris_index tProximity = mGeometricProximity( aGeometryIndex );

            MORIS_ERROR( 
                    tProximity != MORIS_INDEX_MAX,
                    "GEN::Geometric_Proximity::get_geometric_proximity() - "
                    "Geometric proximity for geometry #%i not set yet.", 
                    aGeometryIndex );

            return tProximity;
        }

        //-----------------------------------------------------------------------------------------

    }    // namespace ge
}    // namespace moris
