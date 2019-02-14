/*
 * cl_FEM_Interpolation_Rule_Bis.hpp
 *
 *  Created on: Jan 28, 2019
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_INTERPOLATION_RULE_BIS_HPP_
#define SRC_FEM_CL_FEM_INTERPOLATION_RULE_BIS_HPP_

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
        class Interpolation_Rule_Bis
        {
            const mtk::Geometry_Type      		mGeometryType;
            const Interpolation_Type      		mSpaceInterpolationType;
            const mtk::Interpolation_Order     	mSpaceInterpolationOrder;
            const Interpolation_Type      		mTimeInterpolationType;
            const mtk::Interpolation_Order     	mTimeInterpolationOrder;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------
        /**
         * constructs an interpolation rule
         **
         * @param[ in ] aGeometryType             eg. LINE, QUAD, HEX, ...
         * @param[ in ] aInterpolationType        eg. LAGRANGE, BEZIER, ...
         * @param[ in ] aInterpolationOrder       eg. CONSTANT, LINEAR, QUADRATIC, ...
         */
         Interpolation_Rule_Bis(
                const mtk::Geometry_Type      		& aGeometryType,
                const Interpolation_Type      		& aSpaceInterpolationType,
                const mtk::Interpolation_Order     	& aSpaceInterpolationOrder );

//------------------------------------------------------------------------------
         Interpolation_Rule_Bis(
        		const mtk::Geometry_Type      		& aGeometryType,
                const Interpolation_Type      		& aSpaceInterpolationType,
                const mtk::Interpolation_Order     	& aSpaceInterpolationOrder,
                const Interpolation_Type      		& aTimeInterpolationType,
                const mtk::Interpolation_Order     	& aTimeInterpolationOrder);

//------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~Interpolation_Rule_Bis(){};

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
         * creates a function out of the defined rule
         */
        Interpolation_Function_Base *
        create_space_time_interpolation_function() const;

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
        /**
        * returns the space interpolation type
        */
        Interpolation_Type
        get_space_interpolation_type() const;

//------------------------------------------------------------------------------
        /**
         * returns the space interpolation order
         */
        mtk::Interpolation_Order
		get_space_interpolation_order() const;

//------------------------------------------------------------------------------
        /**
         * returns the time interpolation type
         */
         Interpolation_Type
         get_time_interpolation_type() const;

//------------------------------------------------------------------------------
         /**
          * returns the time interpolation order
          */
         mtk::Interpolation_Order
         get_time_interpolation_order() const;
    };


    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_INTERPOLATION_RULE_BIS_HPP_ */
