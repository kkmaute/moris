/*
 * cl_FEM_IWG_Hamilton_Jacobi_Bulk2.hpp
 *
 *  Created on: Feb 27, 2019
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_IWG_HAMILTON_JACOBI_BULK2_HPP_
#define SRC_FEM_CL_FEM_IWG_HAMILTON_JACOBI_BULK2_HPP_

#include "typedefs.hpp"                     //MRS/COR/src
#include "cl_Cell.hpp"                      //MRS/CON/src

#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src

#include "cl_FEM_IWG.hpp"                   //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        class IWG_Hamilton_Jacobi_Bulk2 : public IWG
        {

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------
            /*
             *  constructor
             */
            IWG_Hamilton_Jacobi_Bulk2();

//------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~IWG_Hamilton_Jacobi_Bulk2(){};

//------------------------------------------------------------------------------
            /**
             * compute the residual
             * r = N_phit  * gradt( phi ) + N_phit * vN * gradx( phi )
             *   = N_phit  * ( gradt( phi ) + vN * gradx( phi ) )
             *
             * @param[ in ] aResidual            residual vector to fill
             * @param[ in ] aFieldInterpolators  list of active field interpolators
             */
            void compute_residual( Matrix< DDRMat >            & aResidual,
                                   Cell< Field_Interpolator* > & aFieldInterpolators );

//------------------------------------------------------------------------------
            /**
             * compute the jacobian
             * j_phiphi = N_phit * Bt_phi + N_phit * vN * Bx_phi = N_phit * ( Bt_phi + vN * Bx_phi )
             * j_phivN  = N_phit * N_vN * gradx( phi )
             *
             * @param[ in ] aJacobians           list of jacobian matrices to fill
             * @param[ in ] aFieldInterpolators  list of active field interpolators
             */
            void compute_jacobian( Cell< Matrix< DDRMat > >    & aJacobians,
                                   Cell< Field_Interpolator* > & aFieldInterpolators );

//------------------------------------------------------------------------------
            /**
             * compute the residual and the jacobian
             *
             * @param[ in ] aJacobians           list of jacobian matrices to fill
             * @param[ in ] aResidual            residual vector to fill
             * @param[ in ] aFieldInterpolators  list of active field interpolators
             *
             */
            void compute_jacobian_and_residual( Cell< Matrix< DDRMat > >    & aJacobians,
                                                Matrix< DDRMat >            & aResidual,
                                                Cell< Field_Interpolator* > & aFieldInterpolators );

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_IWG_HAMILTON_JACOBI_BULK_HPP_ */
