/*
 * cl_FEM_Interpolation_Function_Base_SpaceTime.hpp
 *
 *  Created on: Jan 29, 2019
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_BASE_SPACETIME_HPP_
#define SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_BASE_SPACETIME_HPP_

#include "cl_FEM_Interpolation_Matrix.hpp"
#include "cl_MTK_Enums.hpp" //MTK/src
#include "cl_FEM_Enums.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Base.hpp" //FEM/INT/src

namespace moris
{
    namespace fem
    {

//------------------------------------------------------------------------------

        /**
         * shape function base class for space time
         */
        class Interpolation_Function_Base_SpaceTime
        {

        	Interpolation_Function_Base* mSpaceInterpolation;
        	Interpolation_Function_Base* mTimeInterpolation;


//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             * constructor using two Interpolation_Function_Base
             */
        	Interpolation_Function_Base_SpaceTime(
        			Interpolation_Function_Base & aSpaceInterpolation,
					Interpolation_Function_Base & aTimeInterpolation){};

//------------------------------------------------------------------------------

            /**
             * trivial destructor
             */
            virtual ~Interpolation_Function_Base_SpaceTime(){};

//------------------------------------------------------------------------------
            /**
             * evaluates the shape function at a given point
             *
             * @param[ out ] aN  shape function as
             *                   ( 1 x <number of nodes> )
             *
             * @param[ in ]  aXi parameter coordinates
             *                   ( <number of dimensions>  x 1 )
             * @param[ in ]  aTau parameter coordinates
             *                   ( <number of dimensions>  x 1 )
             */
            virtual void
            eval_N(  Interpolation_Matrix  	& aN,
                     const Matrix< DDRMat > & aXi,
					 const Matrix< DDRMat > & aTau) const = 0;

        };
//------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */



#endif /* SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_BASE_SPACETIME_HPP_ */
