/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Multigeometry.cpp
 *
 */

#include "cl_GEN_Multigeometry.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Multigeometry::Multigeometry(Vector<std::shared_ptr<Geometry>> aGeometries)
                : Field(aGeometries(0))
                , Geometry(aGeometries(0))
                , mGeometries(aGeometries)
        {
            MORIS_ERROR(mGeometries.size() > 0, "A GEN Multigeometry must be created with at least one geometry.");
        }

        //--------------------------------------------------------------------------------------------------------------

        real Multigeometry::get_field_value(
                uint                  aNodeIndex,
                const Matrix<DDRMat>& aCoordinates)
        {
            real tResult = mGeometries(0)->get_field_value(aNodeIndex, aCoordinates);
            for (uint tGeometryIndex = 1; tGeometryIndex < mGeometries.size(); tGeometryIndex++)
            {
                tResult = std::min(tResult, mGeometries(tGeometryIndex)->get_field_value(aNodeIndex, aCoordinates));
            }
            return tResult;
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix<DDRMat>& Multigeometry::get_dfield_dadvs(
                uint                  aNodeIndex,
                const Matrix<DDRMat>& aCoordinates)
        {
            // Find which geometry is the minimum
            real tMin = mGeometries(0)->get_field_value(aNodeIndex, aCoordinates);
            uint tMinGeometryIndex = 0;
            for (uint tGeometryIndex = 1; tGeometryIndex < mGeometries.size(); tGeometryIndex++)
            {
                real tResult = mGeometries(tGeometryIndex)->get_field_value(aNodeIndex, aCoordinates);
                if (tResult < tMin)
                {
                    tMin = tResult;
                    tMinGeometryIndex = tGeometryIndex;
                }
            }

            // Return relevant sensitivity
            return mGeometries(tMinGeometryIndex)->get_dfield_dadvs(aNodeIndex, aCoordinates);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Multigeometry::get_dfield_dcoordinates(
                uint                  aNodeIndex,
                const Matrix<DDRMat>& aCoordinates,
                Matrix<DDRMat>&       aSensitivities)
        {
            // Find which geometry is the minimum
            real tMin = mGeometries(0)->get_field_value(aNodeIndex, aCoordinates);
            uint tMinGeometryIndex = 0;
            for (uint tGeometryIndex = 1; tGeometryIndex < mGeometries.size(); tGeometryIndex++)
            {
                real tResult = mGeometries(tGeometryIndex)->get_field_value(aNodeIndex, aCoordinates);
                if (tResult < tMin)
                {
                    tMin = tResult;
                    tMinGeometryIndex = tGeometryIndex;
                }
            }

            // Get relevant sensitivity
            mGeometries(tMinGeometryIndex)->get_dfield_dcoordinates(aNodeIndex, aCoordinates, aSensitivities);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Multigeometry::add_geometry(std::shared_ptr<Geometry> aGeometry)
        {
            mGeometries.push_back(aGeometry);
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}

