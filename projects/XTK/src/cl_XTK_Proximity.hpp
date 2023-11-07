/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_XTK_Proximity.hpp
 *
 */
#pragma once

#include "cl_Cell.hpp"
#include "enums.hpp"

using namespace moris;

namespace xtk
{

    /**
     * @brief Holds the geometric proximity of a given entity (e.g. cell, facet, vertex) wrt. to each geometry.
     * 
     * Policies to ensure robustness:
     * 1. Proximity values for a given entity must only be set once and are not allowed to be modified thereafter.
     * 2. UNDEFINED proximity values must never be exposed to the outside --> if we find one it always indicates and error.
     * 3. The occurrence of both an INSIDE and OUTSIDE value in a proximity decision always indicates an error.
     * 
     */
    class Proximity
    {
        // ----------------------------------------------------------------------------------

      private:
        Cell< Geometric_Proximity > mGeomProximity;

        // ----------------------------------------------------------------------------------

      public:
        // ----------------------------------------------------------------------------------

        Proximity(
                uint aNumGeometries )
                : mGeomProximity( aNumGeometries, xtk::Geometric_Proximity::UNDEFINED )
        {
            // check that proximity table is initialized with non-zero number of geometries
            MORIS_ASSERT(
                    aNumGeometries > 0,
                    "xtk::Proximity::Proximity() - "
                    "Needs to be constructed with fixed number of geometries greater than one." );
        }

        // ----------------------------------------------------------------------------------

        ~Proximity() {}

        // ----------------------------------------------------------------------------------

        void
        set_value(
                const moris_index              aGeometryIndex,
                const xtk::Geometric_Proximity aProximity );

        // ----------------------------------------------------------------------------------

        xtk::Geometric_Proximity
        get_value( const moris_index aGeometryIndex ) const;

        // ----------------------------------------------------------------------------------

        bool
        is_value_set( const moris_index aGeometryIndex ) const;

        // ----------------------------------------------------------------------------------

    };    // end class: xtk::Proximity

    // ----------------------------------------------------------------------------------
    // ----------------------------------------------------------------------------------

    bool
    do_proximity_values_contradict(
            xtk::Geometric_Proximity aProximityValue1,
            xtk::Geometric_Proximity aProximityValue2 );

    // ----------------------------------------------------------------------------------

    xtk::Geometric_Proximity
    add_two_proximity_values(
            xtk::Geometric_Proximity aProximityValue1,
            xtk::Geometric_Proximity aProximityValue2 );

    // ----------------------------------------------------------------------------------

    xtk::Geometric_Proximity
    decide_proximity_from_parent_proximities( Cell< xtk::Geometric_Proximity > const & aParentProximities );

    // ----------------------------------------------------------------------------------

    xtk::Geometric_Proximity
    decide_proximity_from_parent_proximities(
            Cell< const Proximity* > const & aParentProximities,
            const moris_index                aGeometryIndex );

    // ----------------------------------------------------------------------------------

}    // namespace xtk