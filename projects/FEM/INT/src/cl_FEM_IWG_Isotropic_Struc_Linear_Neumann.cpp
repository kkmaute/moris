
#include "cl_FEM_IWG_Isotropic_Struc_Linear_Neumann.hpp"

#include "fn_trans.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Neumann::compute_residual( moris::Cell< Matrix< DDRMat > > & aResidual )
        {
            // check master field interpolators, properties and constitutive models
            this->check_dof_field_interpolators();
            this->check_dv_field_interpolators();
            this->check_properties();
            this->check_constitutive_models();

            // set residual size
            this->set_residual( aResidual );

            // compute the residual r_U
            aResidual( 0 ) = - trans( mMasterFI( 0 )->N() ) * trans( mMasterProp( 0 )->val() );
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Neumann::compute_jacobian( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians )
        {
            // check master field interpolators, properties and constitutive models
            this->check_dof_field_interpolators();
            this->check_dv_field_interpolators();
            this->check_properties();
            this->check_constitutive_models();

            // set jacobian size
            this->set_jacobian( aJacobians );

            // compute the jacobian for direct IWG dof dependencies
            // None

            // compute the jacobian for indirect IWG dof dependencies through properties
            for( uint iDOF = 0; iDOF < mMasterGlobalDofTypes.size(); iDOF++ )
            {
                // get dof type
                Cell< MSI::Dof_Type > tDofType = mMasterGlobalDofTypes( iDOF );

                // if dependency in the dof type
                if ( mMasterProp( 0 )->check_dof_dependency( tDofType ) )
                {
                    aJacobians( 0 )( iDOF ).matrix_data()
                    += - trans( mMasterFI( 0 )->N() ) * mMasterProp( 0 )->dPropdDOF( tDofType );
                }
            }
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Neumann::compute_jacobian_and_residual( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                                                                     moris::Cell< Matrix< DDRMat > >                & aResidual )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Neumann::compute_jacobian_and_residual - This function does nothing.");
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
