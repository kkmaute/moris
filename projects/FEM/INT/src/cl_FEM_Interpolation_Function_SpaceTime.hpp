/*
 * cl_FEM_Interpolation_Function_SpaceTime.hpp
 *
 *  Created on: Jan 29, 2019
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_SPACETIME_HPP_
#define SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_SPACETIME_HPP_

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
        class Interpolation_Function_SpaceTime : public Interpolation_Function_Base
        {



//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             * constructor using two Interpolation_Function_Base
             */
        	Interpolation_Function_SpaceTime(){};

//------------------------------------------------------------------------------

            /**
             * trivial destructor
             */
            virtual ~Interpolation_Function_SpaceTime(){};

//------------------------------------------------------------------------------
            void
            eval_N(       Interpolation_Matrix  & aN,
                    const Matrix< DDRMat > 		& aXi ) const;

//------------------------------------------------------------------------------

            void
            eval_dNdXi(       Interpolation_Matrix & adNdXi,
                        const Matrix< DDRMat >    & aXi ) const;

//------------------------------------------------------------------------------

            void
            eval_d2NdXi2 (		 Interpolation_Matrix & ad2NdXi2,
                           const Matrix< DDRMat >     & aXi ) const;

//------------------------------------------------------------------------------

            void
            get_param_coords( Matrix< DDRMat > & aXihat ) const;

//------------------------------------------------------------------------------

            uint
            get_number_of_basis() const
            {
             return B;
             }

//------------------------------------------------------------------------------

        	uint
            get_number_of_dimensions() const
            {
        		return N;
             }

//------------------------------------------------------------------------------

        	mtk::Interpolation_Order
            get_interpolation_order() const;

//------------------------------------------------------------------------------

        	Interpolation_Type
            get_interpolation_type() const
            {
        		return T;
            }


//------------------------------------------------------------------------------

        	Interpolation_Matrix
            create_matrix( const uint & aNumberOfFields,
                           const uint & aDerivativeInSpace,
            			   const uint & aDerivativeInTime ) const;

//------------------------------------------------------------------------------

        	Interpolation_Matrix *
            create_matrix_pointer( const uint & aNumberOfFields,
                                   const uint & aDerivativeInSpace,
            					   const uint & aDerivativeInTime ) const;

//------------------------------------------------------------------------------

        };
//------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */



#endif /* SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_SPACETIME_HPP_ */
