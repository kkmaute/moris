/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_WRK_GEN_Performer.hpp
 *
 */

#ifndef MORIS_CL_WRK_GEN_PERFORMER_HPP
#define MORIS_CL_WRK_GEN_PERFORMER_HPP
#include "cl_WRK_Performer.hpp"
#include "cl_GEN_Geometry_Engine.hpp"
#include "cl_Matrix.hpp"
#include <memory>

namespace moris
{
    namespace wrk
    {
        // Performer interface to the geometry engine
        class Gen_Performer : public Performer
        {
            //------------------------------------------------------------------------------------------------------------------

          private:
            std::shared_ptr< moris::ge::Geometry_Engine > mGeometryEngine;

            //------------------------------------------------------------------------------------------------------------------

          public:
            //------------------------------------------------------------------------------------------------------------------

            Gen_Performer( std::shared_ptr< moris::ge::Geometry_Engine > aGeometryEngine );

            //------------------------------------------------------------------------------------------------------------------

            uint
            get_num_refinement_fields();

            //------------------------------------------------------------------------------------------------------------------

            real
            get_field_value(
                    uint                    aFieldIndex,
                    uint                    aNodeIndex,
                    const Matrix< DDRMat >& aCoordinates );

            //------------------------------------------------------------------------------------------------------------------

            const Cell< uint >&
            get_num_refinements( uint aFieldIndex );
            //------------------------------------------------------------------------------------------------------------------

            const Cell< uint >&
            get_refinement_mesh_indices( uint aFieldIndex );

            //------------------------------------------------------------------------------------------------------------------

            sint
            get_refinement_function_index(
                    uint aFieldIndex,
                    uint aRefinementIndex );

            //------------------------------------------------------------------------------------------------------------------
        };
    }    // namespace wrk
}    // namespace moris

#endif
