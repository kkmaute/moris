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
            // pointer to interpolator
            Field_Interpolator * mFieldInterpolator = nullptr;

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
            IWG_Olsson_CLS_Interface(       Field_Interpolator * aFieldInterpolator,
                                      const real                 aFieldUpperBound,
                                      const real                 aFieldLowerBound,
                                      const real                 aEpsilonParameter );

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
             * @param[ in ] aResidual        residual vector to fill
             * @param[ in ] aFieldNormal     field normal at evaluation point
             * @param[ in ] aInterfaceNormal interface normal at evaluation point
             */
            void compute_residual( Matrix< DDRMat > & aResidual,
                                   Matrix< DDRMat >   aFieldNormal,
                                   Matrix< DDRMat >   aInterfaceNormal );

//------------------------------------------------------------------------------
            /**
             * compute the jacobian - Olsson et al. (2007)
             * j = Nt * ( ( phi_L + phi_U - 2 * phi ) * N * n_phi ) * n
             *   - Nt * epsilon * ( ( Bx * n_phi ) * n_phi ) * n
             *
             * @param[ in ] aJacobian        jacobian matrix to fill
             * @param[ in ] aFieldNormal     field normal at evaluation point
             * @param[ in ] aInterfaceNormal interface normal at evaluation point
             *
             */
            void compute_jacobian( Matrix< DDRMat > & aJacobian,
                                   Matrix< DDRMat >   aFieldNormal,
                                   Matrix< DDRMat >   aInterfaceNormal );

//------------------------------------------------------------------------------
            /**
             * compute the residual and the jacobian
             *
             * @param[ in ] aJacobian        jacobian matrix to fill
             * @param[ in ] aResidual        residual vector to fill
             * @param[ in ] aFieldNormal     field normal at evaluation point
             * @param[ in ] aInterfaceNormal interface normal at evaluation point
             *
             */
            void compute_jacobian_and_residual( Matrix< DDRMat > & aJacobian,
                                                Matrix< DDRMat > & aResidual,
                                                Matrix< DDRMat >   aFieldNormal,
                                                Matrix< DDRMat >   aInterfaceNormal );

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_IWG_OLSSON_CLS_INTERFACE_HPP_ */
