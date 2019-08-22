/*
 * cl_FEM_IWG_Hamilton_Jacobi_Bulk_Test.hpp
 *
 *  Created on: Mar 15, 2019
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_IWG_HAMILTON_JACOBI_BULK_TEST_HPP_
#define SRC_FEM_CL_FEM_IWG_HAMILTON_JACOBI_BULK_TEST_HPP_

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

        class IWG_Hamilton_Jacobi_Bulk_Test : public IWG
        {

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------
            /*
             *  constructor
             */
            IWG_Hamilton_Jacobi_Bulk_Test()
        {
            // set the residual dof type
            mResidualDofType = { MSI::Dof_Type::LS1 };

            // set the active dof type
            mActiveDofTypes = {{ MSI::Dof_Type::LS1 }};
        };

//------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~IWG_Hamilton_Jacobi_Bulk_Test(){};

//------------------------------------------------------------------------------
            /**
             * compute the residual
             * r = N_phit  * gradt( phi ) + N_phit * vN * gradx( phi )
             *   = N_phit  * ( gradt( phi ) + vN * gradx( phi ) )
             *
             * @param[ in ] aResidual            residual vector to fill
             * @param[ in ] aFieldInterpolators  list of active field interpolators
             */
            void compute_residual( Matrix< DDRMat > & aResidual )
            {
                // set field interpolators
                Field_Interpolator* phi = mMasterFI( 0 );

                // velocity field value
                Matrix< DDRMat > aVN( 1, 3, 1.0 );

               //compute the residual
               aResidual = trans( phi->N() ) * ( phi->gradt( 1 ) + aVN * phi->gradx( 1 ) );
            }

//------------------------------------------------------------------------------
            /**
             * compute the jacobian
             * j_phiphi = N_phit * Bt_phi + N_phit * vN * Bx_phi = N_phit * ( Bt_phi + vN * Bx_phi )
             * j_phivN  = N_phit * N_vN * gradx( phi )
             *
             * @param[ in ] aJacobians           list of jacobian matrices to fill
             * @param[ in ] aFieldInterpolators  list of active field interpolators
             */
            void compute_jacobian( moris::Cell< Matrix< DDRMat > > & aJacobians )
            {
                // set field interpolators
                Field_Interpolator* phi = mMasterFI( 0 );

                // velocity field value
                Matrix< DDRMat > aVN( 1, 3, 1.0 );

                // set the jacobian size
                aJacobians.resize( 1 );

                // compute the jacobian Jphiphi
                aJacobians( 0 ) = trans( phi->N() ) * ( phi->Bt() + aVN * phi->Bx() );

            };

//------------------------------------------------------------------------------
            /**
             * compute the residual and the jacobian
             *
             * @param[ in ] aJacobians           list of jacobian matrices to fill
             * @param[ in ] aResidual            residual vector to fill
             * @param[ in ] aFieldInterpolators  list of active field interpolators
             *
             */
            void compute_jacobian_and_residual( moris::Cell< Matrix< DDRMat > >    & aJacobians,
                                                Matrix< DDRMat >                   & aResidual )
            {
                // set field interpolators
                Field_Interpolator* phi = mMasterFI( 0 );

                // velocity field value
                Matrix< DDRMat > aVN( 1, 3, 1.0 );

                //compute the residual
                aResidual = trans( phi->N() ) * ( phi->gradt( 1 ) + aVN * phi->gradx( 1 ) );

                // set the jacobian size
                aJacobians.resize( 1 );

                // compute the jacobian Jphiphi
                aJacobians( 0 ) = trans( phi->N() ) * ( phi->Bt() + aVN * phi->Bx() );
            }

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_IWG_HAMILTON_JACOBI_BULK_TEST_HPP_ */
