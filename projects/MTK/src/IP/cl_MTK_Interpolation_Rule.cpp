/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Interpolation_Rule.cpp
 *
 */

#include "cl_MTK_Interpolation_Rule.hpp"                //MTK/src
#include "cl_MTK_Interpolation_Function_Base.hpp"       //MTK/src
#include "cl_MTK_Interpolation_Function_Factory.hpp"    //MTK/src

namespace moris
{
    namespace mtk
    {
        //------------------------------------------------------------------------------

        Interpolation_Rule::Interpolation_Rule(
                const Geometry_Type&       aGeometryType,
                const Interpolation_Type&  aSpaceInterpolationType,
                const Interpolation_Order& aSpaceInterpolationOrder,
                const Interpolation_Type&  aTimeInterpolationType,
                const Interpolation_Order& aTimeInterpolationOrder )
                : mGeometryType( aGeometryType )
                , mSpaceInterpolationType( aSpaceInterpolationType )
                , mSpaceInterpolationOrder( aSpaceInterpolationOrder )
                , mTimeGeometryType( Geometry_Type::LINE )
                , mTimeInterpolationType( aTimeInterpolationType )
                , mTimeInterpolationOrder( aTimeInterpolationOrder )
        {
        }

        Interpolation_Rule::Interpolation_Rule(
                const Geometry_Type&       aGeometryType,
                const Interpolation_Type&  aSpaceInterpolationType,
                const Interpolation_Order& aSpaceInterpolationOrder,
                const Geometry_Type&       aTimeGeometryType,
                const Interpolation_Type&  aTimeInterpolationType,
                const Interpolation_Order& aTimeInterpolationOrder )
                : mGeometryType( aGeometryType )
                , mSpaceInterpolationType( aSpaceInterpolationType )
                , mSpaceInterpolationOrder( aSpaceInterpolationOrder )
                , mTimeGeometryType( aTimeGeometryType )
                , mTimeInterpolationType( aTimeInterpolationType )
                , mTimeInterpolationOrder( aTimeInterpolationOrder )
        {
        }

        //------------------------------------------------------------------------------

        Interpolation_Function_Base*
        Interpolation_Rule::create_space_interpolation_function() const
        {
            // create the factory
            Interpolation_Function_Factory tFactory;

            // return new interpolation function pointer
            return tFactory.create_interpolation_function(
                    mGeometryType,
                    mSpaceInterpolationType,
                    mSpaceInterpolationOrder );
        }

        //------------------------------------------------------------------------------

        Interpolation_Function_Base*
        Interpolation_Rule::create_time_interpolation_function() const
        {
            // create the factory
            Interpolation_Function_Factory tFactory;

            // return new interpolation function pointer
            return tFactory.create_interpolation_function(
                    mTimeGeometryType,
                    mTimeInterpolationType,
                    mTimeInterpolationOrder );
        }

        //----------------------------------------------------------------------------

        uint
        Interpolation_Rule::get_number_of_param_dimensions() const
        {
            switch ( mGeometryType )
            {
                case Geometry_Type::POINT:
                {
                    return 1;
                }
                case Geometry_Type::LINE:
                {
                    return 1;
                }
                case Geometry_Type::QUAD:
                {
                    return 2;
                }
                case Geometry_Type::TRI:
                {
                    return 2;
                }
                case Geometry_Type::HEX:
                {
                    return 3;
                }
                case Geometry_Type::TET:
                {
                    return 4;
                }
                default:
                {
                    MORIS_ERROR( false,
                            "Interpolation_Rule::get_number_of_param_dimensions - Geometry type not defined here" );
                    return 0;
                }
            }
        }

        //----------------------------------------------------------------------------

    } /* namespace mtk */
} /* namespace moris */

