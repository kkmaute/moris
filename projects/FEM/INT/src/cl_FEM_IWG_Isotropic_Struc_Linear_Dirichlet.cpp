
#include "cl_FEM_IWG_Isotropic_Struc_Linear_Dirichlet.hpp"

#include "fn_trans.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {
        // default selection matrix value function
        moris::Matrix< moris::DDRMat > tMValFunction( moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                                                      moris::Cell< moris::fem::Field_Interpolator* > & aDofFI,
                                                      moris::Cell< moris::fem::Field_Interpolator* > & aDvFI,
                                                      moris::fem::Geometry_Interpolator              * aGeometryInterpolator )
        {
            // get spatial dimension
            uint tSpaceDim = aParameters( 0 ).numel();

            // build selection matrix
            Matrix< DDRMat > tM( tSpaceDim, tSpaceDim, 0.0 );

            // populate selection matrix
            for ( uint i = 0; i < tSpaceDim; i++ )
            {
                tM( i, i ) = aParameters( 0 )( i );
            }

            return tM;
        }

        IWG_Isotropic_Struc_Linear_Dirichlet::IWG_Isotropic_Struc_Linear_Dirichlet()
        {
            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Dirichlet" ] = IWG_Property_Type::DIRICHLET;
            mPropertyMap[ "Select" ] = IWG_Property_Type::SELECT;

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "ElastLinIso" ] = IWG_Constitutive_Type::ELAST_LIN_ISO;

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "DirichletNitsche" ] = IWG_Stabilization_Type::DIRICHLET_NITSCHE;

        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Dirichlet::compute_residual( moris::Cell< Matrix< DDRMat > > & aResidual )
        {
            // check master field interpolators, properties and constitutive models
            this->check_dof_field_interpolators();
            this->check_dv_field_interpolators();

            // selection matrix
            Matrix< DDRMat > tM;

            // set a default selection matrix
            if ( mMasterProp( static_cast< uint >( IWG_Property_Type::SELECT ) ) == nullptr )
            {
                // get spatial dimension
                uint tSpaceDim = mMasterFI( 0 )->get_dof_type().size();

                // set selection matrix as identity
                eye( tSpaceDim, tSpaceDim, tM );
            }
            else
            {
                tM = mMasterProp( static_cast< uint >( IWG_Property_Type::SELECT ) )->val();
            }

            // compute jump
            Matrix< DDRMat > tJump = mMasterFI( 0 )->val() - mMasterProp( static_cast< uint >( IWG_Property_Type::DIRICHLET ) )->val();

            // set residual size
            this->set_residual( aResidual );

            // compute the residual
            aResidual( 0 ) = - trans( mMasterFI( 0 )->N() ) * tM * mMasterCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) )->traction( mNormal )
                             + mMasterCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) )->testTraction( mNormal ) * tM * tJump
                             + mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::DIRICHLET_NITSCHE ) )->val()( 0 ) * trans( mMasterFI( 0 )->N() ) * tM * tJump;
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Dirichlet::compute_jacobian( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians )
        {
            // check master field interpolators, properties and constitutive models
            this->check_dof_field_interpolators();
            this->check_dv_field_interpolators();

//            // set a default selection matrix
//            if ( mMasterProp( static_cast< uint >( IWG_Property_Type::SELECT ) ) == nullptr )
//            {
//                // default selection matrix
//                std::shared_ptr< fem::Property > tPropSelect = std::make_shared< fem::Property >();
//                Matrix< DDRMat > tMDiag( 1, mMasterFI( 0 )->get_dof_type().size(), 1.0 );
//                tPropSelect->set_parameters( { tMDiag } );
//                tPropSelect->set_val_function( tMValFunction );
//                mMasterProp( static_cast< uint >( IWG_Property_Type::SELECT ) ) = tPropSelect;
//            }

            // selection matrix
            Matrix< DDRMat > tM;

            // set a default selection matrix
            if ( mMasterProp( static_cast< uint >( IWG_Property_Type::SELECT ) ) == nullptr )
            {
                // get spatial dimension
                uint tSpaceDim = mMasterFI( 0 )->get_dof_type().size();

                // set selection matrix as identity
                eye( tSpaceDim, tSpaceDim, tM );
            }
            else
            {
                tM = mMasterProp( static_cast< uint >( IWG_Property_Type::SELECT ) )->val();
            }

            // compute jump
            Matrix< DDRMat > tJump = mMasterFI( 0 )->val() - mMasterProp( static_cast< uint >( IWG_Property_Type::DIRICHLET ) )->val();

            // set the jacobian size
            this->set_jacobian( aJacobians );

            // compute the jacobian for direct dof dependencies
            aJacobians( 0 )( 0 ) = mMasterCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) )->testTraction( mNormal ) * tM * mMasterFI( 0 )->N()
                                 + mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::DIRICHLET_NITSCHE ) )->val()( 0 ) * trans( mMasterFI( 0 )->N() ) * tM * mMasterFI( 0 )->N();

            // compute the jacobian for indirect dof dependencies through properties
            uint tNumDofDependencies = mMasterGlobalDofTypes.size();
            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the dof type
                Cell< MSI::Dof_Type > tDofType = mMasterGlobalDofTypes( iDOF );

                // if dependency on the dof type
                if ( mMasterProp( static_cast< uint >( IWG_Property_Type::DIRICHLET ) )->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    aJacobians( 0 )( iDOF ).matrix_data()
                    += -1.0 * mMasterCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) )->testTraction( mNormal ) * tM * mMasterProp( static_cast< uint >( IWG_Property_Type::DIRICHLET ) )->dPropdDOF( tDofType )
                       - mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::DIRICHLET_NITSCHE ) )->val()( 0 ) * trans( mMasterFI( 0 )->N() ) * tM * mMasterProp( static_cast< uint >( IWG_Property_Type::DIRICHLET ) )->dPropdDOF( tDofType );
                }

                // if dependency on the dof type
                if ( mMasterCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) )->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    aJacobians( 0 )( iDOF ).matrix_data()
                    += - trans( mMasterFI( 0 )->N() ) *  tM * mMasterCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) )->dTractiondDOF( tDofType, mNormal ) ;
                    // + mMasterCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) )->dTestTractiondDOF( tDofType, mNormal ) * tM * tJump;
                }

                // if dependency on the dof type
                if ( mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::DIRICHLET_NITSCHE ) )->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    aJacobians( 0 )( iDOF ).matrix_data()
                    += trans( mMasterFI( 0 )->N() ) * tM * tJump * mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::DIRICHLET_NITSCHE ) )->dSPdMasterDOF( tDofType );
                }
            }
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Dirichlet::compute_jacobian_and_residual( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                                                                  moris::Cell< Matrix< DDRMat > >                & aResidual )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Dirichlet::compute_jacobian_and_residual - This function does nothing.");
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
