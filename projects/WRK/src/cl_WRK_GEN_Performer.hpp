/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_WRK_GEN_Performer.hpp
 *
 */

#pragma once

#include "cl_WRK_Performer.hpp"
#include "cl_GEN_Geometry_Engine.hpp"
#include "cl_Matrix.hpp"
#include <memory>

namespace moris::wrk
{
    // Performer interface to the geometry engine
    class Gen_Performer : public Performer
    {
        //------------------------------------------------------------------------------------------------------------------

      private:
        std::shared_ptr< moris::gen::Geometry_Engine > mGeometryEngine;

        //------------------------------------------------------------------------------------------------------------------

      public:
        //------------------------------------------------------------------------------------------------------------------

        Gen_Performer( std::shared_ptr< moris::gen::Geometry_Engine > aGeometryEngine );

        //------------------------------------------------------------------------------------------------------------------

        uint
        get_num_refinement_fields() override;

        /**
         * Gets the geometric region of a node with respect to a given geometry.
         *
         * @param aGeometryIndex Geometry index
         * @param aNodeIndex Node index
         * @param aCoordinates Node coordinates
         * @return Geometric region
         */
        gen::Geometric_Region get_geometric_region(
                uint                    aFieldIndex,
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates ) override;

        //------------------------------------------------------------------------------------------------------------------

        const Cell< uint >&
        get_num_refinements( uint aFieldIndex ) override;

        //------------------------------------------------------------------------------------------------------------------

        const Cell< uint >&
        get_refinement_mesh_indices( uint aFieldIndex ) override;

        //------------------------------------------------------------------------------------------------------------------

        sint
        get_refinement_function_index(
                uint aFieldIndex,
                uint aRefinementIndex ) override;

        //------------------------------------------------------------------------------------------------------------------
    };
}
