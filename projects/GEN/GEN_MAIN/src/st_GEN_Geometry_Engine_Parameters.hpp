/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * st_GEN_Geometry_Engine_Parameters.hpp
 *
 */

#ifndef MORIS_ST_GEN_GEOMETRY_ENGINE_PARAMETERS_HPP
#define MORIS_ST_GEN_GEOMETRY_ENGINE_PARAMETERS_HPP

#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_GEN_Geometry.hpp"
#include "cl_GEN_Property.hpp"
#include "cl_GEN_Pdv_Enums.hpp"

namespace moris
{
    namespace ge
    {
        struct Geometry_Engine_Parameters
        {
            /**
             * @var mADVs ADVs which the created geometries are based on
             * @var mGeometries Geometries to be used
             * @var mProperties Properties to be used
             * @var mBulkPhases Bulk phases for creating a custom phase table
             * @var mIntersectionMode Intersection mode
             * @var mRequestedIQIs Requested IQI names
             * @var mGeometryFieldFile File name for writing geometry field values
             * @var mOutputMeshFile File name for writing an exodus mesh
             * @var mTimeOffset Time offset for writing sequential meshes
             */
            Matrix<DDRMat>                  mADVs = {{}};
            Vector<std::shared_ptr<Geometry>> mGeometries = {};
            Vector<std::shared_ptr<Property>> mProperties = {};
            Matrix<DDUMat>                  mBulkPhases = {{}};
            Vector<std::string>               mRequestedIQIs = {};
            std::string                     mGeometryFieldFile = "";
            std::string                     mOutputMeshFile = "";
            real                            mTimeOffset = 0.0;
        };
    }
}

#endif //MORIS_ST_GEN_GEOMETRY_ENGINE_PARAMETERS_HPP

