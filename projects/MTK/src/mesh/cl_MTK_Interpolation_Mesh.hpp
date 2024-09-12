/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Interpolation_Mesh.hpp
 *
 */

#pragma once

#include "cl_MTK_Mesh_Core.hpp"

#include "assert.hpp"
#include "cl_Matrix.hpp"

namespace moris::mtk
{
    class Interpolation_Mesh : public virtual Mesh
    {
      public:
        Interpolation_Mesh(){};

        ~Interpolation_Mesh() override{};

        /**
         * Gets a background mesh of this interpolation mesh (could be this mesh itself)
         *
         * @return Background mesh pointer
         */
        virtual const Interpolation_Mesh& get_background_mesh()
        {
            // There is currently no use of this base class definition, so for now we give a warning.
            // In the future if this functionality is intended and implemented for all child classes where it makes sense, this warning can be removed.
            MORIS_LOG_WARNING( "A background mesh was requested from an interpolation mesh which doesn't yet return one. Assuming it is already a background mesh." );
            return *this;
        };
    };
}    // namespace moris::mtk
