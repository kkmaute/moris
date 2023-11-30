/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Geometric_Proximity.hpp
 *
 */

#ifndef MORIS_CL_Geometric_Proximity_HPP_
#define MORIS_CL_Geometric_Proximity_HPP_

#include "cl_Cell.hpp"
#include "typedefs.hpp"
namespace moris
{
    namespace ge
    {
        class Geometric_Proximity
        {
          public:
            //-----------------------------------------------------------------------------------------

            moris_index mAssociatedVertexIndex = MORIS_INDEX_MAX;

            // Keeps track of a vertex proximity to each geometry ( NumVerts x NumGeometries)
            // 0 - G(x) < threshold (outside)
            // 1 - G(x) == threshold (interface)
            // 2 - G(x) > threshold (inside)
            // MORIS_INDEX_MAX - proximity not set
            moris::Cell< moris_index > mGeometricProximity;    // input: vertex index || output: proximity value

            //-----------------------------------------------------------------------------------------
            //-----------------------------------------------------------------------------------------

            Geometric_Proximity();

            //-----------------------------------------------------------------------------------------

            Geometric_Proximity( moris_index const &aNumGeometries );

            //-----------------------------------------------------------------------------------------

            ~Geometric_Proximity();

            //-----------------------------------------------------------------------------------------

            void
            set_num_geometries( const uint aNumGeometries );

            //-----------------------------------------------------------------------------------------

            void
            set_geometric_proximity( 
                    const moris_index aGeometricProximity, 
                    const moris_index aGeometryIndex );

            //-----------------------------------------------------------------------------------------

            bool
            is_geometric_proximity_set( const moris_index aGeometricProximity ) const;

            //-----------------------------------------------------------------------------------------

            moris_index
            get_geometric_proximity( const moris_index aGeometryIndex ) const;

            //-----------------------------------------------------------------------------------------

        };    // class Geometric_Proximity

        //-----------------------------------------------------------------------------------------

    }    // namespace ge
}    // namespace moris

#endif    // MORIS_CL_Geometric_Proximity_HPP_
