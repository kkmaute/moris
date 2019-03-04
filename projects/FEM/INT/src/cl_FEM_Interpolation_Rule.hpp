/*
 * cl_FEM_Interpolation_Rule.cpp
 *
 *  Created on: Jul 18, 2018
 *      Author: messe
 */

#ifndef SRC_FEM_CL_FEM_INTERPOLATION_RULE_HPP_
#define SRC_FEM_CL_FEM_INTERPOLATION_RULE_HPP_

#include "cl_Matrix.hpp"   //LINALG/src

#include "cl_MTK_Enums.hpp" //MTK/src

#include "cl_FEM_Enums.hpp"                          //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Base.hpp"    //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Factory.hpp" //FEM/INT/src

namespace moris
{
    namespace fem
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
            const mtk::Geometry_Type       mGeometryType;
            const Interpolation_Type       mSpaceInterpolationType;
            const mtk::Interpolation_Order mSpaceInterpolationOrder;
            const Interpolation_Type       mTimeInterpolationType;
            const mtk::Interpolation_Order mTimeInterpolationOrder;

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
        Interpolation_Rule( const mtk::Geometry_Type       & aGeometryType,
                            const Interpolation_Type       & aSpaceInterpolationType,
                            const mtk::Interpolation_Order & aSpaceInterpolationOrder,
                            const Interpolation_Type       & aTimeInterpolationType,
                            const mtk::Interpolation_Order & aTimeInterpolationOrder);

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
        mtk::Geometry_Type get_geometry_type() const
        {
            return mGeometryType;
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
        mtk::Interpolation_Order get_space_interpolation_order() const
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
        mtk::Interpolation_Order get_time_interpolation_order() const
        {
            return mTimeInterpolationOrder;
        }

    };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_INTERPOLATION_RULE_HPP_ */
