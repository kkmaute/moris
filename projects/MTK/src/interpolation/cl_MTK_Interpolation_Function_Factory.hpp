/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Interpolation_Function_Factory.hpp
 *
 */

#ifndef SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_FACTORY_HPP_
#define SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_FACTORY_HPP_

#include "moris_typedefs.hpp"    //MRS/COR/src
#include "cl_MTK_Enums.hpp"      //MTK/src

namespace moris::mtk
{
    //------------------------------------------------------------------------------

    // forward declarations for shape function rule
    class Interpolation_Function_Base;

    // forward declarations for shape function base class
    class Interpolation_Rule;
    //------------------------------------------------------------------------------

    class Interpolation_Function_Factory
    {
        //------------------------------------------------------------------------------

      public:
        //------------------------------------------------------------------------------

        /**
         * trivial constructor
         */
        Interpolation_Function_Factory(){};

        //------------------------------------------------------------------------------

        /**
         * trivial destructor
         */
        ~Interpolation_Function_Factory(){};

        //------------------------------------------------------------------------------

        /**
         * creates an interpolation function according to a given rule
         */
        Interpolation_Function_Base *create_interpolation_function(
                const Geometry_Type       &aGeometryType,
                const Interpolation_Type  &aInterpolationType,
                const Interpolation_Order &aInterpolationOrder );

        //------------------------------------------------------------------------------

        //------------------------------------------------------------------------------

      private:
        //------------------------------------------------------------------------------

        Interpolation_Function_Base *
        create_lagrange_quad( const Interpolation_Order &aInterpolationOrder );

        //------------------------------------------------------------------------------

        Interpolation_Function_Base *
        create_lagrange_hex( const Interpolation_Order &aInterpolationOrder );

        //------------------------------------------------------------------------------

        Interpolation_Function_Base *
        create_lagrange_bar( const Interpolation_Order &aInterpolationOrder );

        //------------------------------------------------------------------------------

        Interpolation_Function_Base *
        create_lagrange_tri( const Interpolation_Order &aInterpolationOrder );

        //------------------------------------------------------------------------------

        Interpolation_Function_Base *
        create_lagrange_tet( const Interpolation_Order &aInterpolationOrder );

        //------------------------------------------------------------------------------

        Interpolation_Function_Base *
        create_constant_bar( const Interpolation_Order &aInterpolationOrder );

        //------------------------------------------------------------------------------

        Interpolation_Function_Base *
        create_constant_point( const Interpolation_Order &aInterpolationOrder );

        //------------------------------------------------------------------------------
    };
    //------------------------------------------------------------------------------
}    // namespace moris::mtk

#endif /* SRC_MTK_CL_MTK_INTERPOLATION_FUNCTION_FACTORY_HPP_ */
