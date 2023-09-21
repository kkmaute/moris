/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Level_Set_Geometry.hpp
 *
 */

#pragma once

#include "cl_GEN_Design_Field.hpp"
#include "cl_GEN_Geometry.hpp"
#include "GEN_Data_Types.hpp"

namespace moris::ge
{
    enum class Int_Interpolation
    {
        LINEAR,
        MULTILINEAR
    };

    /**
     * This is a struct used to simplify \ref moris::ge::Level_Set_Geometry constructors. It contains all field and level-set parameters.
     */
    struct Level_Set_Parameters : public Field_Parameters
    {
        Int_Interpolation mIntersectionInterpolation = Int_Interpolation::LINEAR; // The type of interpolation used to determine intersection location
        real          mIsocontourThreshold = 0.0;                         // Level set isocontour level
        real          mIsocontourTolerance = 1e-12;                       // Interface tolerance based on geometry value
        real          mIntersectionTolerance = 1e-12;                     // Interface tolerance based on intersecction distance
    };

    class Level_Set_Geometry : public Design_Field, public Geometry
    {
    private:
        Level_Set_Parameters mParameters;

    public:

        /**
         * Constructor taking in a field pointer and a set of parameters.
         *
         * @param aField Field for computing nodal values
         * @param aParameters Field parameters
         */
        Level_Set_Geometry(
              std::shared_ptr< Field > aField,
              Level_Set_Parameters     aParameters = {} );

        /**
         * Gets the intersection interpolation type for this geometry.
         *
         * @return Intersection interpolation
         */
        Int_Interpolation get_intersection_interpolation();

        /**
         * Gets the mode of intersection used for this geometry
         *
         * @return Intersection_Mode enum
         */
        Intersection_Mode get_intersection_mode();

        /**
         * Accesses the isocontour level that determines the interface for this geometry
         *
         * @return the isocontour level that determines the geometry interface
         */
        real get_isocontour_threshold();

        /**
         * Acccesses the isocontour tolerance for this geometry
         *
         * @return isocontour tolerance
         */
        real get_isocontour_tolerance();

        /**
         * Accesses the intersection tolerance for this geometry
         *
         * @return The real value of the intersection tolerance
         */
        real get_intersection_tolerance();

    };
}
