/*
 * cl_FEM_Interpolation_Rule.cpp
 *
 *  Created on: Jul 18, 2018
 *      Author: messe
 */

#ifndef SRC_FEM_CL_FEM_INTERPOLATION_RULE_HPP_
#define SRC_FEM_CL_FEM_INTERPOLATION_RULE_HPP_

#include "cl_MTK_Enums.hpp" //MTK/src
#include "cl_FEM_Enums.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Base.hpp" //FEM/INT/src
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
            const Interpolation_Type       mSpaceTimeInterpolationType;
            const mtk::Interpolation_Order mSpaceTimeInterpolationOrder;

            //! flag telling if integration rule is a combination of two
            const bool                     mHasTwoRulesFlag;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

        /**
         * constructs an interpolation rule
         *
         ** @param[ in ] aGeometryType            eg. QUAD, HEX ...
         * @param[ in ] aInterpolationType        eg. Lagrange, Bezier ...
         * @param[ in ] aInterpolationOrder       eg. LINEAR, SERENDIPITY, QUADRATIC
         *
         */
        Interpolation_Rule(
                const mtk::Geometry_Type      		& aGeometryType,
                const Interpolation_Type      		& aSpaceTimeInterpolationType,
                const mtk::Interpolation_Order     	& aSpaceTimeInterpolationOrder );
//------------------------------------------------------------------------------

        Interpolation_Rule(
        		const mtk::Geometry_Type      		& aGeometryType,
                const Interpolation_Type      		& aSpaceInterpolationType,
                const mtk::Interpolation_Order     	& aSpaceInterpolationOrder,
                const Interpolation_Type      		& aTimeInterpolationType,
                const mtk::Interpolation_Order     	& aTimeInterpolationOrder);

//------------------------------------------------------------------------------

        /**
         * trivial destructor
         */
        ~Interpolation_Rule(){};

//------------------------------------------------------------------------------

        auto
        has_two_rules() const -> decltype( mHasTwoRulesFlag )
        {
            return mHasTwoRulesFlag;
        }

//------------------------------------------------------------------------------

        /**
         * creates a function out of the defined rule
         */
        Interpolation_Function_Base * create_space_time_interpolation_function() const;

//------------------------------------------------------------------------------

        /**
         * creates a function out of the defined rule
         */
        Interpolation_Function_Base *
        create_space_interpolation_function() const;

//------------------------------------------------------------------------------

        /**
         * creates a function out of the defined rule
         */
        Interpolation_Function_Base *
        create_time_interpolation_function() const;

//------------------------------------------------------------------------------

        /**
         * returns the interpolation primitive
         */
        auto
        get_geometry_type() const -> decltype ( mGeometryType )
        {
            return mGeometryType;
        }

//------------------------------------------------------------------------------

        Interpolation_Type
        get_type_in_space() const;

//------------------------------------------------------------------------------

        mtk::Interpolation_Order
        get_order_in_space() const;
//------------------------------------------------------------------------------

    };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_INTERPOLATION_RULE_HPP_ */
