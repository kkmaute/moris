/*
 * cl_FEM_IWG_Hamilton_Jacobi_Bulk.hpp
 *
 *  Created on: Feb 13, 2019
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_IWG_HAMILTON_JACOBI_BULK_HPP_
#define SRC_FEM_CL_FEM_IWG_HAMILTON_JACOBI_BULK_HPP_

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

        class IWG_Hamilton_Jacobi_Bulk : public IWG
        {
            // pointer to interpolator
            Field_Interpolator * mFieldInterpolator = nullptr;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------
            /*
             *  constructor
             */
            IWG_Hamilton_Jacobi_Bulk( Field_Interpolator * aFieldInterpolator );

//------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~IWG_Hamilton_Jacobi_Bulk(){};

//------------------------------------------------------------------------------
            /**
             * compute the residual
             * r = Nt  * gradt( phi ) + Nt * v * gradx( phi )
             *   = Nt  * ( gradt( phi ) + v * gradx( phi ) )
             *
             * @param[ in ] aResidual      residual vector to fill
             * @param[ in ] aVelocityField velocity field at evaluation point
             */
            void compute_residual( Matrix< DDRMat > & aResidual,
                                   Matrix< DDRMat >   aVelocityField );

//------------------------------------------------------------------------------
            /**
             * compute the jacobian
             * j = Nt * Bt + Nt * v * Bx
             *   = Nt * (Bt + v * Bx)
             *
             * @param[ in ] aJacobian jacobian matrix to fill
             * @param[ in ] aVelocityField velocity field at evaluation point
             */
            void compute_jacobian( Matrix< DDRMat > & aJacobian,
                                   Matrix< DDRMat >   aVelocityField );

//------------------------------------------------------------------------------
            /**
             * compute the residual and the jacobian
             *
             * @param[ in ] aJacobian jacobian matrix to fill
             * @param[ in ] aResidual residual vector to fill
             * @param[ in ] aVelocityField velocity field at evaluation point
             *
             */
            void compute_jacobian_and_residual( Matrix< DDRMat > & aJacobian,
                                                Matrix< DDRMat > & aResidual,
                                                Matrix< DDRMat >   aVelocityField );

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_IWG_HAMILTON_JACOBI_BULK_HPP_ */
