/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Set_User_Info.cpp
 *
 */
#include "cl_FEM_Set_User_Info.hpp"
#include "cl_FEM_IWG.hpp"
#include "cl_FEM_IQI.hpp"

namespace moris::fem
{
    void Set_User_Info::print_names()
    {
        // print the mesh set name
        std::cout << "Mesh set name: " << mMeshSetName << std::endl;

        // print the bool for time sideset
        std::cout << "Bool for time sideset: " << mTimeContinuity << std::endl;

        // print IWG names
        for ( uint iIWG = 0; iIWG < mIWGs.size(); iIWG++ )
        {
            std::cout << "IWG name: " << mIWGs( iIWG )->get_name() << std::endl;
        }

        // print IQI names
        for ( uint iIQI = 0; iIQI < mIQIs.size(); iIQI++ )
        {
            std::cout << "IQI name: " << mIQIs( iIQI )->get_name() << std::endl;
        }
    }
}    // namespace moris::fem
