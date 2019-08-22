/*
 * cl_FEM_IWG_Olsson_CLS_Interface.hpp
 *
 *  Created on: Feb 13, 2019
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_IWG_OLSSON_CLS_INTERFACE_HPP_
#define SRC_FEM_CL_FEM_IWG_OLSSON_CLS_INTERFACE_HPP_

#include "typedefs.hpp"                     //MRS/COR/src
#include "cl_Cell.hpp"                      //MRS/CON/src

#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src

#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_IWG.hpp"                   //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        class IWG_Olsson_CLS_Interface : public IWG
        {
            // field upper and lower bound
            real mPhiUB;
            real mPhiLB;

            // Olsson CLS epsilon parameter
            real mEpsilon;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------
            /*
             *  constructor
             */
//            IWG_Olsson_CLS_Interface( const real aFieldUpperBound,
//                                      const real aFieldLowerBound,
//                                      const real aEpsilonParameter );
            IWG_Olsson_CLS_Interface();

//------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~IWG_Olsson_CLS_Interface(){};

//------------------------------------------------------------------------------
            /**
             * compute the residual - Olsson et al. (2007)
             * r = Nt * ( ( phi - phi_L) * ( phi_U - phi ) * n_phi ) * n
             *   - Nt * epsilon * ( ( gradx( phi ) * n_phi ) * n_phi ) * n
             *
             * @param[ in ] aResidual            residual vector to fill
             * @param[ in ] aFieldInterpolators  list of active field interpolators
             */
            void compute_residual( Matrix< DDRMat >                   & aResidual );

//------------------------------------------------------------------------------
            /**
             * compute the jacobian - Olsson et al. (2007)
             * j = Nt * ( ( phi_L + phi_U - 2 * phi ) * N * n_phi ) * n
             *   - Nt * epsilon * ( ( Bx * n_phi ) * n_phi ) * n
             *
             * @param[ in ] aJacobians           list of jacobian matrices to fill
             * @param[ in ] aFieldInterpolators  list of active field interpolators
             *
             */
            void compute_jacobian( moris::Cell< Matrix< DDRMat > >    & aJacobians );

//------------------------------------------------------------------------------
            /**
             * compute the residual and the jacobian
             *
             * @param[ in ] aJacobians           list of jacobian matrices to fill
             * @param[ in ] aResidual            residual vector to fill
             * @param[ in ] aFieldInterpolators  list of active field interpolators
             */
            void compute_jacobian_and_residual( moris::Cell< Matrix< DDRMat > >    & aJacobians,
                                                Matrix< DDRMat >                   & aResidual );

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_IWG_OLSSON_CLS_INTERFACE_HPP_ */
