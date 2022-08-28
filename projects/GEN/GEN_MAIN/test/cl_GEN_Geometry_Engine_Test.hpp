/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Geometry_Engine_Test.hpp
 *
 */

#ifndef MORIS_CL_GEN_GEOMETRY_ENGINE_TEST_HPP
#define MORIS_CL_GEN_GEOMETRY_ENGINE_TEST_HPP

#include "cl_GEN_Geometry_Engine.hpp"

namespace moris
{
    namespace ge
    {
        /**
         * Alternate geometry engine class for testing (provides access to some needed protected members)
         */
        class Geometry_Engine_Test : public Geometry_Engine
        {
        public:
            /**
             * Constructor
             *
             * @param aMesh Mesh for getting B-spline information
             * @param aParameters Optional geometry engine parameters
             */
            Geometry_Engine_Test(
                    mtk::Interpolation_Mesh*   aMesh = nullptr,
                    Geometry_Engine_Parameters aParameters = {});

            /**
             * Gets a geometry
             *
             * @param aGeometryIndex Geometry index
             * @return Geometry
             */
            std::shared_ptr<Geometry> get_geometry(uint aGeometryIndex);

            /**
             * Gets a property
             *
             * @param aPropertyIndex Geometry index
             * @return Property
             */
            std::shared_ptr<Property> get_property(uint aPropertyIndex);
        };
    }
}

#endif //MORIS_CL_GEN_GEOMETRY_ENGINE_TEST_HPP

