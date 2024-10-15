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

#include "cl_Vector.hpp"
#include "cl_GEN_Level_Set_Geometry.hpp"
#include "cl_GEN_Surface_Mesh_Geometry.hpp"
#include "cl_GEN_Geometry.hpp"
#include "cl_GEN_Property.hpp"
#include "cl_Submodule_Parameter_Lists.hpp"

// Forward declarations
namespace moris
{
    class Library_IO;
    namespace gen
    {
        class ADV_Manager;
    }
    namespace mtk
    {
        class Mesh;
    }
}

namespace moris::gen
{
    class Design_Factory
    {
      private:
        Vector< std::shared_ptr< Field > > mFields;
        Vector< std::shared_ptr< Geometry > > mGeometries;
        Vector< std::shared_ptr< Property > > mProperties;

      public:
        /**
         * Design factory constructor
         *
         * @param aParameterLists Parameter lists for creating designs and related objects such as fields.
         * @param aADVManager ADV manager, constaining an ADV vector and bounds
         * @param aLibrary Pointer to library for loading user-defined functions
         * @param aMesh MTK mesh used for defining some fields
         * @param aNodeManager Node manager from the geometry engine, if applicable
         */
        Design_Factory(
                Submodule_Parameter_Lists            aParameterLists,
                ADV_Manager&                         aADVManager,
                const std::shared_ptr< Library_IO >& aLibrary     = nullptr,
                mtk::Mesh*                           aMesh        = nullptr,
                Node_Manager&                        aNodeManager = Node_Manager::get_trivial_instance() );

        /**
         * Gets the geometries that the factory has created.
         *
         * @return vector of geometries
         */
        Vector< std::shared_ptr< Geometry > > get_geometries();

        /**
         * Gets the properties that the factory has created.
         *
         * @return vector of properties
         */
        Vector< std::shared_ptr< Property > > get_properties();
    };
}
