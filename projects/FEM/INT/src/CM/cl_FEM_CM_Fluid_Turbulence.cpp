//FEM/INt/src
#include "cl_FEM_CM_Fluid_Turbulence.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"
//LINALG/src
#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"
#include "op_minus.hpp"

namespace moris
{
    namespace fem
    {

        //--------------------------------------------------------------------------------------------------------------

        CM_Fluid_Turbulence::CM_Fluid_Turbulence()
        {
            // set the property pointer cell size
            mProperties.resize( static_cast< uint >( CM_Property_Type::MAX_ENUM ), nullptr );

            // populate the map
            mPropertyMap[ "Density" ]   = static_cast< uint >( CM_Property_Type::DENSITY );
            mPropertyMap[ "Viscosity" ] = static_cast< uint >( CM_Property_Type::VISCOSITY );
            mPropertyMap[ "KinViscosity" ] = static_cast< uint >( CM_Property_Type::KIN_VISCOSITY );
        }

        //------------------------------------------------------------------------------

        void CM_Fluid_Turbulence::set_dof_type_list(
                moris::Cell< moris::Cell< MSI::Dof_Type > > aDofTypes,
                moris::Cell< std::string >                  aDofStrings )
        {
            // set dof type list
            Constitutive_Model::set_dof_type_list( aDofTypes );

            // loop over the provided dof type
            for( uint iDof = 0; iDof < aDofTypes.size(); iDof++ )
            {
                // get dof type string
                std::string tDofString = aDofStrings( iDof );

                // get dof type
                MSI::Dof_Type tDofType = aDofTypes( iDof )( 0 );

                // switch on dof type string
                if( tDofString == "Velocity" )
                {
                    mDofVelocity = tDofType;
                }
                else if( tDofString == "Pressure" )
                {
                    mDofPressure = tDofType;
                }
                else if( tDofString == "Viscosity" )
                {
                    mDofViscosity = tDofType;
                }
                else
                {
                    // error unknown dof string
                    MORIS_ERROR( false ,
                            "CM_Fluid_Turbulence::set_dof_type_list - Unknown aDofString : %s \n",
                            tDofString.c_str() );
                }
            }
        }

        //------------------------------------------------------------------------------

        void CM_Fluid_Turbulence::set_local_properties()
        {
            // set the density property
            mPropDensity = get_property( "Density" );

            // set the dynamic viscosity property
            mPropViscosity = get_property( "Viscosity" );

            // set the kinematic viscosity property
            mPropKinViscosity = get_property( "KinViscosity" );
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Fluid_Turbulence::eval_flux()
        {
            // get viscous contribution to flux
            CM_Fluid_Incompressible::eval_flux();

            // get the turbulence viscosity
            real tViscosityT = this->compute_turbulence_viscosity();

            // compute flux
            mFlux += 2.0 * tViscosityT * this->strain();
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Fluid_Turbulence::eval_dFluxdDOF(
                const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get viscous contribution to dFluxdDOF
            CM_Fluid_Incompressible::eval_dFluxdDOF( aDofTypes );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // if velocity dof
            if( aDofTypes( 0 ) == mDofVelocity )
            {
                // get the turbulence viscosity
                real tViscosityT = this->compute_turbulence_viscosity();

                // build dfluxdv
                mdFluxdDof( tDofIndex ) +=
                        2.0 * tViscosityT * this->dStraindDOF( aDofTypes );
            }

            // evaluate derivative of the turbulence viscosity
            Matrix< DDRMat > tdviscositytdu;
            this->compute_dviscositytdu( aDofTypes, tdviscositytdu );

            // add contribution from turbulence viscosity
            mdFluxdDof( tDofIndex ) +=
                    2.0 * this->strain() * tdviscositytdu;
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Fluid_Turbulence::eval_divflux()
        {
            // get viscous contribution to div flux
            CM_Fluid_Incompressible::eval_divflux();

            // get the turbulence viscosity
            real tViscosityT = this->compute_turbulence_viscosity();

            // get the gradx for the turbulence viscosity
            Matrix< DDRMat > tdViscosityTdx;
            this->compute_dviscositytdx( tdViscosityTdx );

            // flatten dviscositytdx
            Matrix< DDRMat > tdViscosityTdxFlat;
            this->flatten_normal( tdViscosityTdx, tdViscosityTdxFlat );

            // compute flux
            mDivFlux +=
                    2.0 * tViscosityT * this->divstrain() +
                    2.0 * tdViscosityTdxFlat * this->strain();
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Fluid_Turbulence::eval_ddivfluxdu( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get viscous contribution to ddivfluxdu
            CM_Fluid_Incompressible::eval_ddivfluxdu( aDofTypes );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // if velocity dof
            if( aDofTypes( 0 ) == mDofVelocity )
            {
                // get the turbulence viscosity
                real tViscosityT = this->compute_turbulence_viscosity();

                // get the gradx for the turbulence viscosity
                Matrix< DDRMat > tdViscosityTdx;
                this->compute_dviscositytdx( tdViscosityTdx );

                // flatten dviscositytdx
                Matrix< DDRMat > tdViscosityTdxFlat;
                this->flatten_normal( tdViscosityTdx, tdViscosityTdxFlat );

                // add contribution to ddivstrain/dv
                mddivfluxdu( tDofIndex ) +=
                        2.0 * tViscosityT * this->ddivstraindu( aDofTypes ) +
                        2.0 * tdViscosityTdxFlat * this->dStraindDOF( aDofTypes );
            }

            // get the turbulence viscosity
            Matrix< DDRMat > tdviscositytdu;
            this->compute_dviscositytdu( aDofTypes, tdviscositytdu );

            // add contribution to ddivstrain/du
            mddivfluxdu( tDofIndex ) +=
                    2.0 * this->divstrain() * tdviscositytdu;

            // get the derivative of gradx for the turbulence viscosity wrt dof
            Matrix< DDRMat > tdViscosityTdxdu;
            this->compute_dviscositytdxdu( aDofTypes, tdViscosityTdxdu );

            const Matrix< DDRMat > & tStrain = this->strain();
            Matrix< DDRMat > tStrainFull;
            switch ( mSpaceDim )
            {
                case 2:
                {
                    tStrainFull = {
                            { tStrain( 0 ), tStrain( 2 ) },
                            { tStrain( 2 ), tStrain( 1 ) } };
                    break;
                }

                case 3:
                {
                    tStrainFull = {
                            { tStrain( 0 ), tStrain( 5 ), tStrain( 4 ) },
                            { tStrain( 5 ), tStrain( 1 ), tStrain( 3 ) },
                            { tStrain( 4 ), tStrain( 3 ), tStrain( 2 ) }};
                    break;
                }

                default:
                    MORIS_ERROR( false, "CM_Fluid_Turbulence::eval_ddivfluxdu - only 2 or 3D" );
            }

            //
            mddivfluxdu( tDofIndex ) += 2.0 * tStrainFull * tdViscosityTdxdu;
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Fluid_Turbulence::eval_dfluxdx( uint aOrder )
        {
            // get viscous contribution dfluxdx
            CM_Fluid_Incompressible::eval_dfluxdx( aOrder );

            // only 1st order supported
            MORIS_ERROR( aOrder == 1, "CM_Fluid_Incompressible::eval_dfluxdx - only 1st order supported." );

            // get the turbulence viscosity
            real tViscosityT = this->compute_turbulence_viscosity();

            // evaluate dfluxdx
            mdFluxdx( aOrder - 1 ) -=
                    2.0 * tViscosityT * this->dstraindx( aOrder );

            // FIXME need d/dx part related to tViscosityT
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Fluid_Turbulence::eval_traction( const Matrix< DDRMat > & aNormal )
        {
            // flatten the normal
            Matrix< DDRMat > tFlatNormal;
            this->flatten_normal( aNormal, tFlatNormal );

            // compute the traction
            mTraction = tFlatNormal * this->flux();
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Fluid_Turbulence::eval_testTraction(
                const Matrix< DDRMat >             & aNormal,
                const moris::Cell< MSI::Dof_Type > & aTestDofTypes )
        {
            // get viscous contribution to test traction
            CM_Fluid_Incompressible::eval_testTraction( aNormal, aTestDofTypes );

            // get test dof type index
            uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

            // flatten the normal
            Matrix< DDRMat > tFlatNormal;
            this->flatten_normal( aNormal, tFlatNormal );

            // if test traction wrt velocity
            if( aTestDofTypes( 0 ) == mDofVelocity )
            {
                // get the turbulence viscosity
                real tViscosityT = this->compute_turbulence_viscosity();

                // compute test traction wrt velocity
                mTestTraction( tTestDofIndex ) +=
                        2.0 * tViscosityT * tFlatNormal * this->testStrain();
            }

            // evaluate test turbulence viscosity
            Matrix< DDRMat > tdviscositytdutest;
            this->compute_dviscositytdu( aTestDofTypes, tdviscositytdutest );

            // add contribution from turbulence viscosity
            mTestTraction( tTestDofIndex ) +=
                    2.0 * tFlatNormal * this->strain() * tdviscositytdutest;
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Fluid_Turbulence::eval_dTractiondDOF(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                const Matrix< DDRMat >             & aNormal )
        {
            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // flatten normal
            Matrix< DDRMat > tFlatNormal;
            this->flatten_normal( aNormal, tFlatNormal );

            // compute dtractiondu
            mdTractiondDof( tDofIndex ) = tFlatNormal * this->dFluxdDOF( aDofTypes );
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Fluid_Turbulence::eval_dTestTractiondDOF(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                const Matrix< DDRMat >             & aNormal,
                const Matrix< DDRMat >             & aJump,
                const moris::Cell< MSI::Dof_Type > & aTestDofTypes )
        {
            // get viscous contribution
            CM_Fluid_Incompressible::eval_dTestTractiondDOF(
                    aDofTypes,
                    aNormal,
                    aJump,
                    aTestDofTypes );

            // get test dof type index
            uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

            // get the test dof FI
            Field_Interpolator * tFITest =
                    mFIManager->get_field_interpolators_for_type( aTestDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the derivative dof FI
            Field_Interpolator * tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // evaluate derivative of the turbulence viscosity
            Matrix< DDRMat > tdviscositytduder;
            this->compute_dviscositytdu( aDofTypes, tdviscositytduder );

            // evaluate test of the turbulence viscosity
            Matrix< DDRMat > tdviscositytdutest;
            this->compute_dviscositytdu( aTestDofTypes, tdviscositytdutest );

            // flatten normal
            Matrix< DDRMat > tFlatNormal;
            this->flatten_normal( aNormal, tFlatNormal );

            // compute contribution to dTestTractiondDof
            mdTestTractiondDof( tTestDofIndex )( tDofIndex ) +=
                    2.0 * trans( tFlatNormal * this->dStraindDOF( aTestDofTypes ) ) * aJump * tdviscositytduder +
                    2.0 * trans( tdviscositytdutest ) * trans( aJump ) * tFlatNormal * this->dStraindDOF( aDofTypes );

            // if viscosity is the test dof
            if( aTestDofTypes( 0 ) == mDofViscosity && aDofTypes( 0 ) == mDofViscosity )
            {
                MORIS_LOG_INFO( "Missing second order derivative - FD for now" );

                mdTestTractiondDof( tTestDofIndex )( tDofIndex ).set_size(
                        tFITest->get_number_of_space_time_coefficients(),
                        tFIDer->get_number_of_space_time_coefficients() );

                Constitutive_Model::eval_dtesttractiondu_FD(
                        aDofTypes,
                        aTestDofTypes,
                        mdTestTractiondDof( tTestDofIndex )( tDofIndex ),
                        1e-6,
                        aNormal,
                        aJump,
                        fem::FDScheme_Type::POINT_1_FORWARD );

                //MORIS_ERROR( false, "CM_Fluid_Turbulence::eval_dTestTractiondDOF - Case not implemented so far, require d2viscositytdviscosity2" );
            }
        }

        //------------------------------------------------------------------------------

        real CM_Fluid_Turbulence::compute_turbulence_viscosity()
        {
            // init the turbulence viscosity
            real tViscosityT = 0.0;

            // get the viscosity dof type FI
            Field_Interpolator * tFIModViscosity =
                    mFIManager->get_field_interpolators_for_type( mDofViscosity );

            // get the modified viscosity value
            real tModViscosity = tFIModViscosity->val()( 0 );

            // if modified viscosity is positive
            if( tModViscosity >= 0.0 )
            {
                // compute fv1
                real tFv1 = compute_fv1(
                        { mDofViscosity },
                        mFIManager,
                        mPropKinViscosity );

                // compute turbulent viscosity
                tViscosityT = mPropDensity->val()( 0 ) * tModViscosity * tFv1;
            }

            // return the turbulence viscosity
            return tViscosityT;
        }

        //------------------------------------------------------------------------------

        void CM_Fluid_Turbulence::compute_dviscositytdu(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adviscositytdu )
        {
            // get the dof type FI
            Field_Interpolator * tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set matrix size
            adviscositytdu.set_size(
                    1,
                    tFIDer->get_number_of_space_time_coefficients() );

            // get the viscosity dof type FI
            Field_Interpolator * tFIModViscosity =
                    mFIManager->get_field_interpolators_for_type( mDofViscosity );

            // get the modified viscosity value
            real tModViscosity = tFIModViscosity->val()( 0 );

            // if modified viscosity is positive
            if( tModViscosity >= 0.0 )
            {
                // compute fv1
                real tFv1 = compute_fv1(
                        { mDofViscosity },
                        mFIManager,
                        mPropKinViscosity );

                // compute dfv1du
                Matrix< DDRMat > tdfv1du;
                compute_dfv1du(
                        { mDofViscosity },
                        mFIManager,
                        mPropKinViscosity,
                        aDofTypes,
                        tdfv1du );

                // add contribution from dfv1du
                adviscositytdu = mPropDensity->val()( 0 ) * tFIModViscosity->val() * tdfv1du;

                // if dof type is viscosity
                if( aDofTypes( 0 ) == mDofViscosity )
                {
                    // add contribution to dSPdu
                    adviscositytdu += mPropDensity->val()( 0 ) * tFv1 * tFIModViscosity->N();
                }

                // if density depends on dof
                if( mPropDensity->check_dof_dependency( aDofTypes ) )
                {
                    // add contribution from drhodu
                    adviscositytdu += tFIModViscosity->val() * tFv1 * mPropDensity->dPropdDOF( aDofTypes );
                }
            }
            else
            {
                adviscositytdu.fill( 0.0 );
            }
        }

        //------------------------------------------------------------------------------

        void CM_Fluid_Turbulence::compute_dviscositytdx(
                Matrix< DDRMat > & adviscositytdx )
        {
            // set matrix size
            adviscositytdx.set_size( mSpaceDim, 1 );

            // get the viscosity dof type FI
            Field_Interpolator * tFIModViscosity =
                    mFIManager->get_field_interpolators_for_type( mDofViscosity );

            // get the modified viscosity value
            real tModViscosity = tFIModViscosity->val()( 0 );

            // if modified viscosity is positive
            if( tModViscosity >= 0.0 )
            {
                // compute fv1
                real tFv1 = compute_fv1(
                        { mDofViscosity },
                        mFIManager,
                        mPropKinViscosity );

                // compute dfv1dx
                Matrix< DDRMat > tdfv1dx;
                compute_dfv1dx(
                        { mDofViscosity },
                        mFIManager,
                        mPropKinViscosity,
                        tdfv1dx );

                // compute dviscositytdx
                adviscositytdx =
                        mPropDensity->val()( 0 ) * tFIModViscosity->gradx( 1 ) * tFv1 +
                        mPropDensity->val()( 0 ) * tdfv1dx * tModViscosity;
                        // FIXME tFv1 * tModViscosity * mPropDensity->dPropdx();
            }
            else
            {
                adviscositytdx.fill( 0.0 );
            }
        }

        //------------------------------------------------------------------------------

        void CM_Fluid_Turbulence::compute_dviscositytdxdu(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adviscositytdxdu )
        {
            // get the dof type FI
            Field_Interpolator * tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set matrix size
            adviscositytdxdu.set_size(
                    mSpaceDim,
                    tFIDer->get_number_of_space_time_coefficients() );

            // get the viscosity dof type FI
            Field_Interpolator * tFIModViscosity =
                    mFIManager->get_field_interpolators_for_type( mDofViscosity );

            // get the modified viscosity value
            real tModViscosity = tFIModViscosity->val()( 0 );

            // if modified viscosity is positive
            if( tModViscosity >= 0.0 )
            {
                // compute fv1
                real tFv1 = compute_fv1(
                        { mDofViscosity },
                        mFIManager,
                        mPropKinViscosity );

                // compute dfv1dx
                Matrix< DDRMat > tdfv1dx;
                compute_dfv1dx(
                        { mDofViscosity },
                        mFIManager,
                        mPropKinViscosity,
                        tdfv1dx );

                // compute dfv1du
                Matrix< DDRMat > tdfv1du;
                compute_dfv1du(
                        { mDofViscosity },
                        mFIManager,
                        mPropKinViscosity,
                        aDofTypes,
                        tdfv1du );

                // compute dfv1dxdu
                Matrix< DDRMat > tdfv1dxdu;
                compute_dfv1dxdu(
                        { mDofViscosity },
                        mFIManager,
                        mPropKinViscosity,
                        aDofTypes,
                        tdfv1dxdu );

                // add contribution from dfv1du
                adviscositytdxdu =
                        mPropDensity->val()( 0 ) * tFIModViscosity->gradx( 1 ) * tdfv1du +
                        mPropDensity->val()( 0 ) * tFIModViscosity->val()( 0 ) * tdfv1dxdu;
                // FIXME mPropDensity->dPropdx() * tFIModViscosity->val()( 0 ) * tdfv1du

                // if dof type is viscosity
                if( aDofTypes( 0 ) == mDofViscosity )
                {
                    // add contribution to dviscositytdxdu
                    adviscositytdxdu +=
                            mPropDensity->val()( 0 ) * tFv1 * tFIModViscosity->dnNdxn( 1 ) +
                            mPropDensity->val()( 0 ) * tdfv1dx * tFIModViscosity->N();
                    // FIXME tFv1 * mPropDensity->dPropdx() * tFIModViscosity->N()
                 }

                // if density depends on dof
                if( mPropDensity->check_dof_dependency( aDofTypes ) )
                {
                    adviscositytdxdu +=
                            tFIModViscosity->gradx( 1 ) * tFv1 * mPropDensity->dPropdDOF( aDofTypes ) +
                             tdfv1dx * tFIModViscosity->val()( 0 ) * mPropDensity->dPropdDOF( aDofTypes );
                    // FIXME + tFIModViscosity->val()( 0 ) * tFv1 * mPropDensity->dPropdxdDOF( aDofTypes )
                }
            }
            else
            {
                adviscositytdxdu.fill( 0.0 );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
