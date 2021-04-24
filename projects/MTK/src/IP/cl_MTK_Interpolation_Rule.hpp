/*
 * cl_MTK_Interpolation_Rule.hpp
 *
 *  Created on: Jul 18, 2018
 *      Author: messe
 */

#ifndef SRC_MTK_CL_MTK_INTERPOLATION_RULE_HPP_
#define SRC_MTK_CL_MTK_INTERPOLATION_RULE_HPP_

#include "cl_Matrix.hpp"   //LINALG/src

#include "cl_MTK_Enums.hpp" //MTK/src

#include "cl_MTK_Enums.hpp"                          //MTK/src
#include "cl_MTK_Interpolation_Function_Base.hpp"    //MTK/src
#include "cl_MTK_Interpolation_Function_Factory.hpp" //MTK/src

namespace moris
{
    namespace mtk
    {
        //------------------------------------------------------------------------------

        // forward declaration for shape function
        class Interpolation_Function_Base;

        //------------------------------------------------------------------------------

        /**
         *  /brief a container that defines an interpolation function
         */
        class Interpolation_Rule
        {
                const Geometry_Type       mGeometryType;
                const Interpolation_Type  mSpaceInterpolationType;
                const Interpolation_Order mSpaceInterpolationOrder;
                const Geometry_Type       mTimeGeometryType;
                const Interpolation_Type  mTimeInterpolationType;
                const Interpolation_Order mTimeInterpolationOrder;

                //------------------------------------------------------------------------------
            public:
                //------------------------------------------------------------------------------

                /**
                 * constructs an interpolation rule
                 **
                 * @param[ in ] aGeometryType             eg. QUAD, HEX ...
                 * @param[ in ] aSpaceInterpolationType   eg. Constant, Lagrange, Bezier ...
                 * @param[ in ] aSpaceInterpolationOrder  eg. LINEAR, SERENDIPITY, QUADRATIC
                 * @param[ in ] aTimeInterpolationType    eg. Constant, Lagrange, Bezier ...
                 * @param[ in ] aTimeInterpolationOrder   eg. CONSTANT, LINEAR, SERENDIPITY, QUADRATIC
                 *
                 */
                Interpolation_Rule(
                        const Geometry_Type        & aGeometryType,
                        const Interpolation_Type   & aSpaceInterpolationType,
                        const Interpolation_Order  & aSpaceInterpolationOrder,
                        const Interpolation_Type   & aTimeInterpolationType,
                        const Interpolation_Order  & aTimeInterpolationOrder);

                /**
                 * constructs an interpolation rule
                 **
                 * @param[ in ] aGeometryType             eg. QUAD, HEX ...
                 * @param[ in ] aSpaceInterpolationType   eg. Constant, Lagrange, Bezier ...
                 * @param[ in ] aSpaceInterpolationOrder  eg. LINEAR, SERENDIPITY, QUADRATIC
                 * @param[ in ] aTimeInterpolationType    eg. Constant, Lagrange, Bezier ...
                 * @param[ in ] aTimeInterpolationOrder   eg. CONSTANT, LINEAR, SERENDIPITY, QUADRATIC
                 *
                 */
                Interpolation_Rule(
                        const Geometry_Type        & aGeometryType,
                        const Interpolation_Type   & aSpaceInterpolationType,
                        const Interpolation_Order  & aSpaceInterpolationOrder,
                        const Geometry_Type        & aTimeGeometryType,
                        const Interpolation_Type   & aTimeInterpolationType,
                        const Interpolation_Order  & aTimeInterpolationOrder);

                //------------------------------------------------------------------------------

                /**
                 * trivial destructor
                 */
                ~Interpolation_Rule(){};

                //------------------------------------------------------------------------------

                /**
                 * creates a function out of the defined rule
                 */
                Interpolation_Function_Base * create_space_interpolation_function() const;

                //------------------------------------------------------------------------------

                /**
                 * creates a function out of the defined rule
                 */
                Interpolation_Function_Base * create_time_interpolation_function() const;

                //------------------------------------------------------------------------------

                /**
                 * returns the interpolation primitive
                 */
                Geometry_Type get_geometry_type() const
                {
                    return mGeometryType;
                }

                //------------------------------------------------------------------------------

                /**
                 * returns the interpolation primitive
                 */
                Geometry_Type get_time_geometry_type() const
                {
                    return mTimeGeometryType;
                }

                //------------------------------------------------------------------------------
                /**
                 * returns the interpolation type in space
                 */
                Interpolation_Type get_space_interpolation_type() const
                {
                    return mSpaceInterpolationType;
                }

                //------------------------------------------------------------------------------
                /**
                 * returns the interpolation order in space
                 */
                Interpolation_Order get_space_interpolation_order() const
                {
                    return mSpaceInterpolationOrder;
                }

                //------------------------------------------------------------------------------
                /**
                 * returns the interpolation type in time
                 */
                Interpolation_Type get_time_interpolation_type() const
                {
                    return mTimeInterpolationType;
                }

                //------------------------------------------------------------------------------
                /**
                 * returns the interpolation order in time
                 */
                Interpolation_Order get_time_interpolation_order() const
                {
                    return mTimeInterpolationOrder;
                }

                //------------------------------------------------------------------------------
                /**
                 * returns the number of param dimensions without having to create
                 * an interpolation function base
                 */
                uint get_number_of_param_dimensions() const;

        };

        //------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */

#endif /* SRC_MTK_CL_MTK_INTERPOLATION_RULE_HPP_ */
