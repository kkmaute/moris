/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_GEN_create_properties.hpp
 *
 */

#ifndef MORIS_FN_GEN_CREATE_PROPERTIES_HPP
#define MORIS_FN_GEN_CREATE_PROPERTIES_HPP

#include "cl_GEN_Property.hpp"
#include "cl_GEN_Geometry.hpp"
#include "cl_Param_List.hpp"
#include "cl_Library_IO.hpp"

namespace moris
{
    namespace mtk
    {
        class Mesh;
    }

    namespace ge
    {
        /**
         * Higher-level call for creating a cell of properties, which ensures that all property dependencies are
         * resolved correctly.
         *
         * @tparam Vector_Type Type of vector where ADVs are stored
         * @param aPropertyParameterLists Parameter lists for creating Property classes
         * @param aADVs ADV vector
         * @param aGeometries Already-created fields for finding dependencies
         * @param aLibrary Pointer to library for loading user-defined functions
         * @return Pointer to specific Property class
         */
        template <typename Vector_Type>
        Vector<std::shared_ptr<Property>> create_properties(
                Vector<ParameterList>             aPropertyParameterLists,
                Vector_Type&                    aADVs,
                Vector<std::shared_ptr<Geometry>> aGeometries = {},
                std::shared_ptr<Library_IO>     aLibrary = nullptr);

        /**
         * Creates an instance of the specified Property class and returns a shared pointer to it.
         *
         * @tparam Vector_Type Type of vector where ADVs are stored
         * @param aPropertyParameterList Parameter list for creating a Property class
         * @param aADVs ADV vector
         * @param aLibrary Pointer to library for loading user-defined functions
         * @return Pointer to specific Property class
         */
        template <typename Vector_Type>
        std::shared_ptr<Property> create_property(
                ParameterList                aPropertyParameterList,
                Vector_Type&                 aADVs,
                Vector<std::shared_ptr<Field>> aFieldDependencies = {},
                std::shared_ptr<Library_IO>  aLibrary = nullptr);
    }
}

#endif //MORIS_FN_GEN_CREATE_PROPERTIES_HPP

