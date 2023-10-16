/*
* Copyright (c) 2022 University of Colorado
* Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
*
*------------------------------------------------------------------------------------
*
* cl_GEN_Design_Factory.hpp
*
*/

#pragma once

#include "cl_Cell.hpp"
#include "cl_GEN_Level_Set_Geometry.hpp"
#include "cl_GEN_Property.hpp"
#include "cl_Param_List.hpp"

// Forward declarations
namespace moris
{
    class Library_IO;
    namespace mtk
    {
        class Mesh;
    }
}

namespace moris::ge
{
    class Design_Factory
    {
      private:
        Cell< std::shared_ptr< Design_Field > > mDesigns;
        Cell< std::shared_ptr< Level_Set_Geometry > > mGeometries;
        Cell< std::shared_ptr< Property > > mProperties;

      public:
        /**
         * Design factory constructor
         *
         * @param aParameterLists Parameter lists for creating designs and related objects such as fields.
         * @param aADVs ADV vector
         * @param aLibrary Pointer to library for loading user-defined functions
         * @param aMesh MTK mesh used for defining some fields
         */
        Design_Factory(
                Cell< ParameterList >         aParameterLists,
                Matrix< DDRMat >&             aADVs,
                std::shared_ptr< Library_IO > aLibrary = nullptr,
                mtk::Mesh*                    aMesh = nullptr );

        /**
         * Gets the geometries that the factory has created.
         *
         * @return vector of geometries
         */
        Cell< std::shared_ptr< Level_Set_Geometry > > get_geometries();

        /**
         * Gets the properties that the factory has created.
         *
         * @return vector of properties
         */
        Cell< std::shared_ptr< Property > > get_properties();
    };
}
