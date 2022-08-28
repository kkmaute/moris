/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Field_Proxy.cpp
 *
 */

#include "../MTK_Test_Proxy/cl_MTK_Field_Proxy.hpp"

#include "fn_dot.hpp"

#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src

namespace moris
{
    namespace mtk
    {
        Field_Proxy::Field_Proxy(
                mtk::Mesh_Pair   aMeshPairs,
                uint const     & aDiscretizationMeshIndex )
        : Field( aMeshPairs,
                 aDiscretizationMeshIndex )
        {
        }

        // ----------------------------------------------------------------------------------------------

        Field_Proxy::~Field_Proxy()
        {

        }

        // ----------------------------------------------------------------------------------------------

    }
}

