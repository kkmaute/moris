/*
 * cl_FEM_Interpolation_Function_Factory.hpp
 *
 *  Created on: Jul 9, 2018
 *      Author: messe
 */


#ifndef SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_FACTORY_HPP_
#define SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_FACTORY_HPP_


#include "typedefs.hpp" //MRS/COR/src
#include "cl_MTK_Enums.hpp" //MTK/src
#include "cl_FEM_Enums.hpp" //FEM/INT/src

namespace moris
{
    namespace fem
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
            Interpolation_Function_Base *  create_interpolation_function(
                    const mtk::Geometry_Type       & aGeometryType,
                    const Interpolation_Type       & aInterpolationType,
                    const mtk::Interpolation_Order & aInterpolationOrder );

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
       private:
//------------------------------------------------------------------------------

            Interpolation_Function_Base *
            create_lagrange_quad( const mtk::Interpolation_Order   & aInterpolationOrder );

//------------------------------------------------------------------------------

            Interpolation_Function_Base *
            create_lagrange_hex( const mtk::Interpolation_Order   & aInterpolationOrder );

//------------------------------------------------------------------------------

            Interpolation_Function_Base *
            create_lagrange_bar( const mtk::Interpolation_Order   & aInterpolationOrder );

//------------------------------------------------------------------------------

            Interpolation_Function_Base *
            create_lagrange_tri( const mtk::Interpolation_Order   & aInterpolationOrder );

//------------------------------------------------------------------------------

            Interpolation_Function_Base *
            create_lagrange_tet( const mtk::Interpolation_Order   & aInterpolationOrder );

//------------------------------------------------------------------------------

            Interpolation_Function_Base *
            create_constant_bar( const mtk::Interpolation_Order   & aInterpolationOrder );

//------------------------------------------------------------------------------

            Interpolation_Function_Base *
            create_constant_point( const mtk::Interpolation_Order   & aInterpolationOrder );

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */


#endif /* SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_FACTORY_HPP_ */
