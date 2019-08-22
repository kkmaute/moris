
#include "cl_FEM_IWG_Isotropic_Spatial_Diffusion_Dirichlet.hpp"

#include "fn_trans.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        IWG_Isotropic_Spatial_Diffusion_Dirichlet::IWG_Isotropic_Spatial_Diffusion_Dirichlet()
        {
            // set the residual dof type
            mResidualDofType = { MSI::Dof_Type::TEMP };

            // set the active dof types
            mActiveDofTypes = { { MSI::Dof_Type::TEMP } };

            // set the active property types
            mActivePropertyTypes = { fem::Property_Type::CONDUCTIVITY,
                                     fem::Property_Type::TEMP_DIRICHLET };

            // FIXME set a penalty
            mGamma = 1.0;

            //FIXME forced diffusion parameter
            //      forced dimensions for 3D
            eye( mSpaceDim, mSpaceDim, mKappa );
            mKappa = 1.0 * mKappa;
        }

//------------------------------------------------------------------------------

        void IWG_Isotropic_Spatial_Diffusion_Dirichlet::compute_residual( Matrix< DDRMat > & aResidual )
        {
            // check master field interpolators
            this->check_field_interpolators();

            // fixme interpolated TBar from nodal values
            Matrix< DDRMat > tTBar = mMasterFI( 0 )->N() * mNodalWeakBCs;

            // compute the residual r_T
            aResidual = - trans( mMasterFI( 0 )->N() ) * dot( mKappa * mMasterFI( 0 )->gradx( 1 ), mNormal )
                        + trans( mKappa * mMasterFI( 0 )->Bx() ) * mNormal * ( mMasterFI( 0 )->val()( 0 ) - tTBar( 0 ) )
                        + mGamma * trans( mMasterFI( 0 )->N() ) * ( mMasterFI( 0 )->val()( 0 ) - tTBar( 0 ) );
        }

//------------------------------------------------------------------------------

        void IWG_Isotropic_Spatial_Diffusion_Dirichlet::compute_jacobian( moris::Cell< Matrix< DDRMat > > & aJacobians )
        {
            // check master field interpolators
            this->check_field_interpolators();

            // set the jacobian size
            aJacobians.resize( 1 );

            // compute the jacobian j_T_T
            aJacobians( 0 ) = - trans( mMasterFI( 0 )->N() ) * trans( mNormal ) * mKappa * mMasterFI( 0 )->Bx()
                              + trans( mKappa * mMasterFI( 0 )->Bx() ) * mNormal * mMasterFI( 0 )->N()
                              + mGamma * trans( mMasterFI( 0 )->N() ) * mMasterFI( 0 )->N();

        }

//------------------------------------------------------------------------------

        void IWG_Isotropic_Spatial_Diffusion_Dirichlet::compute_jacobian_and_residual( moris::Cell< Matrix< DDRMat > > & aJacobians,
                                                                                       Matrix< DDRMat >                & aResidual )
        {
            // check master field interpolators
            this->check_field_interpolators();

            // interpolated TBar from nodal values
            Matrix< DDRMat > tTBar = mMasterFI( 0 )->N() * mNodalWeakBCs;

            // compute the residual r_T
            aResidual = - trans( mMasterFI( 0 )->N() ) * dot( mKappa * mMasterFI( 0 )->gradx( 1 ), mNormal )
                        + trans( mKappa * mMasterFI( 0 )->Bx() ) * mNormal * ( mMasterFI( 0 )->val()( 0 ) - tTBar( 0 ) )
                        + mGamma * trans( mMasterFI( 0 )->N() ) * ( mMasterFI( 0 )->val()( 0 ) - tTBar( 0 ) );

            // set the jacobian size
            aJacobians.resize( 1 );

            // compute the jacobian j_T_T
            aJacobians( 0 ) = - trans( mMasterFI( 0 )->N() ) * trans( mNormal ) * mKappa * mMasterFI( 0 )->Bx()
                              + trans( mKappa * mMasterFI( 0 )->Bx() ) * mNormal * mMasterFI( 0 )->N()
                              + mGamma * trans( mMasterFI( 0 )->N() ) * mMasterFI( 0 )->N();
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
