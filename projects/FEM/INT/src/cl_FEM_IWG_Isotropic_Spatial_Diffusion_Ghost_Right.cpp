
#include "cl_FEM_IWG_Isotropic_Spatial_Diffusion_Ghost_Right.hpp"

#include "fn_trans.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        IWG_Isotropic_Spatial_Diffusion_Ghost_Right::IWG_Isotropic_Spatial_Diffusion_Ghost_Right()
        {
            // set the residual dof type
            mResidualDofType = { MSI::Dof_Type::TEMP };

            // set the active dof type
            mActiveDofTypes = { { MSI::Dof_Type::TEMP } };

            // FIXME set a penalty
            mGammaGhost = 0.0001;

            // FIXME set mesh parameter
            // Mesh parameter or parameters which it depends on must be fed to IWG from outside
            mMeshParameter = 1;

            // FIXME set order of Shape functions
            // Order must be fed to IWG from outside
            mOrder = 1;

            //FIXME forced diffusion parameter
            //      forced dimensions for 3D
            eye( mSpaceDim, mSpaceDim, mKappa );
            mKappa = 1.0 * mKappa;
        }

//------------------------------------------------------------------------------

        void IWG_Isotropic_Spatial_Diffusion_Ghost_Right::compute_residual
            ( Matrix< DDRMat >            & aResidual,
              Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            // check, if order of shape functions is supported
        	MORIS_ERROR( mOrder < 2, " Ghost stabilization for this order of shape fncts. not supported yet. " );

        	// set field interpolator
            Field_Interpolator* tTemp = aFieldInterpolators( 0 );

            // compute the residual r_T

            uint tResSize = tTemp->get_number_of_space_time_coefficients();
            aResidual.set_size( tResSize, 1, 0.0);
            aResidual =   aResidual
            		    + mGammaGhost * mMeshParameter // * mKappa                            // scaling parameters
            		    * trans(tTemp->Bx()) * mNormal * trans(mNormal) * tTemp->gradx( 1 );  // matrices & vectors
        }

//------------------------------------------------------------------------------

        void IWG_Isotropic_Spatial_Diffusion_Ghost_Right::compute_jacobian
            ( Cell< Matrix< DDRMat > >    & aJacobians,
              Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            // check, if order of shape functions is supported
        	MORIS_ERROR( mOrder < 2, " Ghost stabilization for this order of shape fncts. not supported yet. " );

        	// set field interpolator
            Field_Interpolator* tTemp  = aFieldInterpolators( 0 );

            // set the jacobian size
            aJacobians.resize( 1 );

            // compute the jacobian j_T_T

            uint tNumOfBases = tTemp->get_number_of_space_time_bases();
            aJacobians( 0 ).set_size( tNumOfBases, tNumOfBases, 0.0 );

            aJacobians( 0 ) =   aJacobians( 0 )
            		          + mGammaGhost * mMeshParameter //* mKappa                       // scaling parameters
                              * trans(tTemp->Bx()) * mNormal * trans(mNormal) * tTemp->Bx();  // matrices & vectors
        }

//------------------------------------------------------------------------------

        void IWG_Isotropic_Spatial_Diffusion_Ghost_Right::compute_jacobian_and_residual
            ( Cell< Matrix< DDRMat > >    & aJacobians,
              Matrix< DDRMat >            & aResidual,
              Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            // check, if order of shape functions is supported
        	MORIS_ERROR( mOrder < 2, " Ghost stabilization for this order of shape fncts. not supported yet. " );

        	// set field interpolator
            Field_Interpolator* tTemp = aFieldInterpolators( 0 );

            // compute the residual r_T

            uint tResSize = tTemp->get_number_of_space_time_coefficients();
            aResidual.set_size( tResSize, 1, 0.0);
            aResidual =   aResidual
            		    + mGammaGhost * mMeshParameter // * mKappa                            // scaling parameters
            		    * trans(tTemp->Bx()) * mNormal * trans(mNormal) * tTemp->gradx( 1 );  // matrices & vectors

            // set the jacobian size
                       aJacobians.resize( 1 );

            // compute the jacobian j_T_T

            uint tNumOfBases = tTemp->get_number_of_space_time_bases();
            aJacobians( 0 ).set_size( tNumOfBases, tNumOfBases, 0.0 );

            aJacobians( 0 ) =   aJacobians( 0 )
                       		  + mGammaGhost * mMeshParameter //* mKappa                       // scaling parameters
                              * trans(tTemp->Bx()) * mNormal * trans(mNormal) * tTemp->Bx();  // matrices & vectors
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
