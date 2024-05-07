/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * st_GEN_Geometry_Engine_Parameters.hpp
 *
 */

#pragma once

#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_GEN_ADV_Manager.hpp"
#include "cl_GEN_Level_Set_Geometry.hpp"
#include "cl_GEN_Surface_Mesh_Geometry.hpp"
#include "cl_GEN_Property.hpp"
#include "GEN_Data_Types.hpp"

namespace moris::gen
{
    struct Geometry_Engine_Parameters
    {
        /**
         * @var mADVs ADVs which the created geometries are based on
         * @var mGeometries Geometries to be used
         * @var mProperties Properties to be used
         * @var mBulkPhases Bulk phases for creating a custom phase table
         * @var mRequestedIQIs Requested IQI names
         * @var mGeometryFieldFile File name for writing geometry field values
         * @var mOutputMeshFile File name for writing an exodus mesh
         * @var mTimeOffset Time offset for writing sequential meshes
         */
        ADV_Manager                           mADVManager;
        Vector< std::shared_ptr< Geometry > > mGeometries        = {};
        Vector< std::shared_ptr< Property > > mProperties        = {};
        Matrix< DDUMat >                      mBulkPhases        = {{}};
        Vector< std::string >                 mRequestedIQIs     = {};
        std::string                           mGeometryFieldFile = "";
        std::string                           mOutputMeshFile    = "";
        real                                  mTimeOffset        = 0.0;
    };
}
