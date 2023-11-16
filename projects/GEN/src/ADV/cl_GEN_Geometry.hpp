/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Geometry.hpp
 *
 */

#pragma once

#include <memory>

namespace moris
{
    template< typename T >
    class Cell;
    namespace mtk
    {
        class Field;
    }
}

namespace moris::ge
{

    class Geometry
    {
      public:

        /**
         * Constructor
         */
        Geometry();

        /**
         * Gets an MTK field, if this geometry uses one that needs to be remapped to a new mesh
         *
         * @return MTK field
         */
        virtual Cell< std::shared_ptr< mtk::Field > > get_mtk_fields() = 0;
    };
}
