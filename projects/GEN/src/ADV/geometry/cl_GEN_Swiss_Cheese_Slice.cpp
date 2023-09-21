/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Swiss_Cheese_Slice.cpp
 *
 */

#include "cl_GEN_Swiss_Cheese_Slice.hpp"
#include "cl_GEN_Superellipse.hpp"

namespace moris::ge
{
    //--------------------------------------------------------------------------------------------------------------

    Swiss_Cheese_Slice::Swiss_Cheese_Slice(
            real             aLeftBound,
            real             aRightBound,
            real             aBottomBound,
            real             aTopBound,
            uint             aNumXHoles,
            uint             aNumYHoles,
            real             aXSemidiameter,
            real             aYSemidiameter,
            real             aExponent,
            real             aScaling,
            real             aRegularization,
            real             aShift,
            real             aOffset,
            Level_Set_Parameters aParameters)
            : Multigeometry(Cell<std::shared_ptr< Field > >(1, std::make_shared<Superellipse>(
                      aLeftBound,
                      aBottomBound,
                      aXSemidiameter,
                      aYSemidiameter,
                      aExponent,
                      aScaling,
                      aRegularization,
                      aShift,
                      aParameters)))
    {
        this->create_holes(
                aLeftBound,
                aBottomBound,
                aNumXHoles,
                aNumYHoles,
                (aRightBound - aLeftBound) / (aNumXHoles - 1),
                (aTopBound - aBottomBound) / (aNumYHoles - 1),
                aOffset,
                aXSemidiameter,
                aYSemidiameter,
                aExponent,
                aScaling,
                aRegularization,
                aShift,
                aParameters);
    }

    //--------------------------------------------------------------------------------------------------------------

    Swiss_Cheese_Slice::Swiss_Cheese_Slice(real aLeftBound,
            real             aRightBound,
            real             aBottomBound,
            real             aTopBound,
            real             aTargetXSpacing,
            real             aTargetYSpacing,
            real             aXSemidiameter,
            real             aYSemidiameter,
            real             aExponent,
            real             aScaling,
            real             aRegularization,
            real             aShift,
            real             aOffset,
            bool             aAllowLessThanTargetSpacing,
            Level_Set_Parameters                 aParameters)
            : Multigeometry(Cell<std::shared_ptr< Field > >(1, std::make_shared<Superellipse>(
                      aLeftBound,
                      aBottomBound,
                      aXSemidiameter,
                      aYSemidiameter,
                      aExponent,
                      aScaling,
                      aRegularization,
                      aShift,
                      aParameters)))
    {
        real tNumXHoles = ((aRightBound - aLeftBound) / aTargetXSpacing) + 1.0;
        real tNumYHoles = ((aTopBound - aBottomBound) / aTargetYSpacing) + 1.0;
        if (aAllowLessThanTargetSpacing)
        {
            tNumXHoles += 0.5;
            tNumYHoles += 0.5;
        }
        this->create_holes(
                aLeftBound,
                aBottomBound,
                (uint)tNumXHoles,
                (uint)tNumYHoles,
                (aRightBound - aLeftBound) / ((uint)tNumXHoles - 1.0),
                (aTopBound - aBottomBound) / ((uint)tNumYHoles - 1.0),
                aOffset,
                aXSemidiameter,
                aYSemidiameter,
                aExponent,
                aScaling,
                aRegularization,
                aShift,
                aParameters);
    }

    //--------------------------------------------------------------------------------------------------------------

    void Swiss_Cheese_Slice::create_holes(
            real             aXOrigin,
            real             aYOrigin,
            uint             aNumXHoles,
            uint             aNumYHoles,
            real             aXDelta,
            real             aYDelta,
            real             aOffset,
            real             aXSemidiameter,
            real             aYSemidiameter,
            real             aExponent,
            real             aScaling,
            real             aRegularization,
            real             aShift,
            Level_Set_Parameters aParameters)
    {
        for (uint tYHoleIndex = 0; tYHoleIndex < aNumYHoles; tYHoleIndex++)
        {
            for (uint tXHoleIndex = 0; tXHoleIndex < aNumXHoles; tXHoleIndex++)
            {
//                if (tXHoleIndex > 0 or tYHoleIndex > 0) // This hole was already created
//                {
//                    this->add_geometry(std::make_shared<Superellipse>(
//                            aXOrigin + (tXHoleIndex * aXDelta) + std::fmod((tYHoleIndex * aOffset), aXDelta),
//                            aYOrigin + (tYHoleIndex * aYDelta),
//                            aXSemidiameter,
//                            aYSemidiameter,
//                            aExponent,
//                            aScaling,
//                            aRegularization,
//                            aShift,
//                            aParameters));
//                }
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

}
