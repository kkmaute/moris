
#include "cl_FEM_CM_Spalart_Allmaras_Turbulence.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_dot.hpp"
#include "fn_eye.hpp"
#include "fn_clip_value.hpp"
#include "op_minus.hpp"

namespace moris
{
    namespace fem
    {
        //--------------------------------------------------------------------------------------------------------------

        CM_Spalart_Allmaras_Turbulence::CM_Spalart_Allmaras_Turbulence()
        {
            // set the property pointer cell size
            mProperties.resize( static_cast< uint >( CM_Property_Type::MAX_ENUM ), nullptr );

            // populate the map
            mPropertyMap[ "KinViscosity" ] = static_cast< uint >( CM_Property_Type::KINVISCOSITY );
            mPropertyMap[ "WallDistance" ] = static_cast< uint >( CM_Property_Type::WALL_DISTANCE );

            // FIXME for now only 1st order allowed
            uint tOrder = 1;

            // init storage for evaluation
            mdChidx.resize( tOrder );
            mdFndx.resize( tOrder );
            mdDiffusionCoeffdx.resize( tOrder );

            // init flag for evaluation
            mdChidxEval.set_size( tOrder, 1, true );
            mdFndxEval.set_size( tOrder, 1, true );
            mdDiffusionCoeffdxEval.set_size( tOrder, 1, true );
        }

        //------------------------------------------------------------------------------

        void
        CM_Spalart_Allmaras_Turbulence::set_parameters(
                moris::Cell< Matrix< DDRMat > > aParameters )
        {
            // FIXME not necessary
            // set mParameters
            mParameters = aParameters;

            // get number of parameters
            uint tParamSize = aParameters.size();

            // check for proper size of constant function parameters
            MORIS_ERROR( tParamSize <= 2,
                    "CM_Spalart_Allmaras_Turbulence::set_parameters - max 2 constant parameters can to be set." );

            // flag turning on/off ft2
            if ( tParamSize > 0 )
            {
                // check for proper parameter type; here just a scalar
                MORIS_ERROR( aParameters( 0 ).numel() == 1,
                        "CM_Spalart_Allmaras_Turbulence::set_parameters - 1st parameter is not a scalar but a vector." );

                // check for proper parameter value
                MORIS_ERROR( aParameters( 0 )( 0 ) == 0.0 || aParameters( 0 )( 0 ) == 1.0,
                        "CM_Spalart_Allmaras_Turbulence::set_parameters - invalid mUseFt2 parameter \n" );

                // consider or not ft2 value
                mUseFt2 = aParameters( 0 )( 0 ) > 0 ? true : false;
            }

            // if alpha specification
            if ( tParamSize > 1 )
            {
                // check for proper parameter type; here just a scalar
                MORIS_ERROR( aParameters( 1 ).numel() == 1,
                        "CM_Spalart_Allmaras_Turbulence::set_parameters - 2nd parameter is not a scalar but a vector." );

                // check for proper parameter value
                MORIS_ERROR( aParameters( 1 )( 0 ) >= 1.0,
                        "CM_Spalart_Allmaras_Turbulence::set_parameters - invalid mAlpha parameter." );

                // set alpha value
                mAlpha = aParameters( 1 )( 0 );
            }
        }

        //------------------------------------------------------------------------------

        void
        CM_Spalart_Allmaras_Turbulence::set_local_properties()
        {
            // set the viscosity property
            mPropKinViscosity = get_property( "KinViscosity" );

            // set the wall distance property
            mPropWallDistance = get_property( "WallDistance" );
        }

        //------------------------------------------------------------------------------

        void
        CM_Spalart_Allmaras_Turbulence::set_dof_type_list(
                moris::Cell< moris::Cell< MSI::Dof_Type > > aDofTypes,
                moris::Cell< std::string >                  aDofStrings )
        {
            // set dof type list
            Constitutive_Model::set_dof_type_list( aDofTypes );

            // loop over the provided dof type
            for ( uint iDof = 0; iDof < aDofTypes.size(); iDof++ )
            {
                // get dof type string
                const std::string& tDofString = aDofStrings( iDof );

                // get dof type
                MSI::Dof_Type tDofType = aDofTypes( iDof )( 0 );

                // switch on dof type string
                if ( tDofString == "Velocity" )
                {
                    mDofVelocity = tDofType;
                }
                else if ( tDofString == "Viscosity" )
                {
                    mDofViscosity = tDofType;
                }
                else
                {
                    // error unknown dof string
                    MORIS_ERROR( false,
                            "CM_Spalart_Allmaras_Turbulence::set_dof_type_list - Unknown aDofString : %s \n",
                            tDofString.c_str() );
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        CM_Spalart_Allmaras_Turbulence::build_global_dof_type_list()
        {
            // call parent implementation
            Constitutive_Model::build_global_dof_type_list();

            // get number of dof types
            uint tNumGlobalDofTypes = mGlobalDofTypes.size();

            // init child specific eval flags
            mdProductionCoeffduEval.set_size( tNumGlobalDofTypes, 1, true );
            mdProductionTermduEval.set_size( tNumGlobalDofTypes, 1, true );

            mdWallDestructionCoeffduEval.set_size( tNumGlobalDofTypes, 1, true );
            mdWallDestructionTermduEval.set_size( tNumGlobalDofTypes, 1, true );

            mdDiffusionCoeffduEval.set_size( tNumGlobalDofTypes, 1, true );

            mdModVelocityduEval.set_size( tNumGlobalDofTypes, 1, true );

            mdModVelocityLinearizedduEval.set_size( tNumGlobalDofTypes, 1, true );

            mdChiduEval.set_size( tNumGlobalDofTypes, 1, true );

            mdFt2duEval.set_size( tNumGlobalDofTypes, 1, true );

            mdWduEval.set_size( tNumGlobalDofTypes, 1, true );
            mdSduEval.set_size( tNumGlobalDofTypes, 1, true );
            mdFv1duEval.set_size( tNumGlobalDofTypes, 1, true );
            mdFv2duEval.set_size( tNumGlobalDofTypes, 1, true );
            mdSBarduEval.set_size( tNumGlobalDofTypes, 1, true );
            mdSModduEval.set_size( tNumGlobalDofTypes, 1, true );
            mdSTildeduEval.set_size( tNumGlobalDofTypes, 1, true );

            mdRduEval.set_size( tNumGlobalDofTypes, 1, true );
            mdGduEval.set_size( tNumGlobalDofTypes, 1, true );
            mdFwduEval.set_size( tNumGlobalDofTypes, 1, true );
            mdFnduEval.set_size( tNumGlobalDofTypes, 1, true );

            // FIXME for now only 1st order allowed
            uint tOrder = 1;
            mdChidxduEval.set_size( tOrder, tNumGlobalDofTypes, true );
            mdFndxduEval.set_size( tOrder, tNumGlobalDofTypes, true );
            mdDiffusionCoeffdxduEval.set_size( tOrder, tNumGlobalDofTypes, true );

            // init child specific storage
            mdProductionCoeffdu.resize( tNumGlobalDofTypes );
            mdProductionTermdu.resize( tNumGlobalDofTypes );

            mdWallDestructionCoeffdu.resize( tNumGlobalDofTypes );
            mdWallDestructionTermdu.resize( tNumGlobalDofTypes );

            mdDiffusionCoeffdu.resize( tNumGlobalDofTypes );

            mdModVelocitydu.resize( tNumGlobalDofTypes );

            mdModVelocityLinearizeddu.resize( tNumGlobalDofTypes );

            mdChidu.resize( tNumGlobalDofTypes );

            mdFt2du.resize( tNumGlobalDofTypes );

            mdWdu.resize( tNumGlobalDofTypes );
            mdSdu.resize( tNumGlobalDofTypes );
            mdFv1du.resize( tNumGlobalDofTypes );
            mdFv2du.resize( tNumGlobalDofTypes );
            mdSBardu.resize( tNumGlobalDofTypes );
            mdSModdu.resize( tNumGlobalDofTypes );
            mdSTildedu.resize( tNumGlobalDofTypes );

            mdRdu.resize( tNumGlobalDofTypes );
            mdGdu.resize( tNumGlobalDofTypes );
            mdFwdu.resize( tNumGlobalDofTypes );

            mdFndu.resize( tNumGlobalDofTypes );

            mdFndxdu.resize( tOrder );
            mdDiffusionCoeffdxdu.resize( tOrder );
            mdChidxdu.resize( tOrder );
            for ( uint iOrder = 0; iOrder < tOrder; iOrder++ )
            {
                mdChidxdu( iOrder ).resize( tNumGlobalDofTypes );
                mdDiffusionCoeffdxdu( iOrder ).resize( tNumGlobalDofTypes );
                mdFndxdu( iOrder ).resize( tNumGlobalDofTypes );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Spalart_Allmaras_Turbulence::reset_eval_flags()
        {
            // call parent implementation
            Constitutive_Model::reset_eval_flags();

            // reset child specific eval flags for SA model
            mProductionCoeffEval = true;
            mdProductionCoeffduEval.fill( true );
            mProductionTermEval = true;
            mdProductionTermduEval.fill( true );

            mWallDestructionCoeffEval = true;
            mdWallDestructionCoeffduEval.fill( true );
            mWallDestructionTermEval = true;
            mdWallDestructionTermduEval.fill( true );

            mDiffusionCoeffEval = true;
            mdDiffusionCoeffduEval.fill( true );
            mdDiffusionCoeffdxEval.fill( true );
            mdDiffusionCoeffdxduEval.fill( true );

            mModVelocityEval = true;
            mdModVelocityduEval.fill( true );

            mModVelocityLinearizedEval = true;
            mdModVelocityLinearizedduEval.fill( true );

            mChiEval = true;
            mdChiduEval.fill( true );
            mdChidxEval.fill( true );
            mdChidxduEval.fill( true );

            mFt2Eval = true;
            mdFt2duEval.fill( true );

            mWEval = true;
            mdWduEval.fill( true );
            mSEval = true;
            mdSduEval.fill( true );
            mFv1Eval = true;
            mdFv1duEval.fill( true );
            mFv2Eval = true;
            mdFv2duEval.fill( true );
            mSBarEval = true;
            mdSBarduEval.fill( true );
            mSModEval = true;
            mdSModduEval.fill( true );
            mSTildeEval = true;
            mdSTildeduEval.fill( true );

            mREval = true;
            mdRduEval.fill( true );
            mGEval = true;
            mdGduEval.fill( true );
            mFwEval = true;
            mdFwduEval.fill( true );

            mFnEval = true;
            mdFnduEval.fill( true );
            mdFndxEval.fill( true );
            mdFndxduEval.fill( true );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Spalart_Allmaras_Turbulence::eval_divflux()
        {
            // get the residual viscosity FI
            Field_Interpolator* tFIViscosity =
                    mFIManager->get_field_interpolators_for_type( mDofViscosity );

            // compute the divergence of the flux
            mDivFlux =
                    trans( this->ddiffusioncoeffdx( 1 ) ) * tFIViscosity->gradx( 1 )
                    + this->diffusion_coefficient() * sum( tFIViscosity->gradx( 2 )( { 0, mSpaceDim - 1 }, { 0, 0 } ) );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Spalart_Allmaras_Turbulence::eval_ddivfluxdu(
                const moris::Cell< MSI::Dof_Type >& aDofTypes )
        {
            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the corresponding FI
            Field_Interpolator* tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set size for ddivflux/du
            mddivfluxdu( tDofIndex ).set_size( mSpaceDim, tFIDer->get_number_of_space_time_coefficients() );

            // get the residual viscosity FI
            Field_Interpolator* tFIViscosity =
                    mFIManager->get_field_interpolators_for_type( mDofViscosity );

            // add contribution from diffusion coefficient to derivative
            mddivfluxdu( tDofIndex ) =
                    trans( tFIViscosity->gradx( 1 ) ) * this->ddiffusioncoeffdxdu( aDofTypes, 1 )
                    + sum( tFIViscosity->gradx( 2 )( { 0, mSpaceDim - 1 }, { 0, 0 } ) ) * this->ddiffusioncoeffdu( aDofTypes );

            // if derivative wrt to residual dof type (here viscosity)
            if ( aDofTypes( 0 ) == mDofViscosity )
            {
                // get number of space time bases
                uint tNumBases = tFIViscosity->get_number_of_space_time_bases();

                // get second order shape function derivatives
                Matrix< DDRMat > td2Ndx2( 1, tNumBases, 0.0 );
                for ( uint iSpace = 0; iSpace < mSpaceDim; iSpace++ )
                {
                    td2Ndx2 += tFIViscosity->dnNdxn( 2 ).get_row( iSpace );
                }

                // add contribution to derivative
                mddivfluxdu( tDofIndex ) +=
                        trans( this->ddiffusioncoeffdx( 1 ) ) * tFIViscosity->dnNdxn( 1 ) + this->diffusion_coefficient() * td2Ndx2;
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Spalart_Allmaras_Turbulence::eval_traction(
                const Matrix< DDRMat >& aNormal )
        {
            // get the viscosity dof FI
            Field_Interpolator* tFIViscosity =
                    mFIManager->get_field_interpolators_for_type( mDofViscosity );

            // compute the traction
            mTraction = this->diffusion_coefficient() * trans( tFIViscosity->gradx( 1 ) ) * aNormal;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Spalart_Allmaras_Turbulence::eval_testTraction(
                const Matrix< DDRMat >&             aNormal,
                const moris::Cell< MSI::Dof_Type >& aTestDofTypes )
        {
            // get the test FI
            Field_Interpolator* tFITest =
                    mFIManager->get_field_interpolators_for_type( aTestDofTypes( 0 ) );

            // get test dof type index
            uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

            // init the derivative of the test traction
            mTestTraction( tTestDofIndex ).set_size( tFITest->get_number_of_space_time_coefficients(), 1 );

            // get the residual dof FI (here viscosity)
            Field_Interpolator* tFIViscosity =
                    mFIManager->get_field_interpolators_for_type( mDofViscosity );

            // add contribution of diffusion to test traction
            mTestTraction( tTestDofIndex ) =
                    trans( this->ddiffusioncoeffdu( aTestDofTypes ) ) * trans( tFIViscosity->gradx( 1 ) ) * aNormal;

            // if derivative dof type is viscosity dof
            if ( aTestDofTypes( 0 ) == mDofViscosity )
            {
                // add contribution to test traction
                mTestTraction( tTestDofIndex ) +=
                        trans( tFIViscosity->dnNdxn( 1 ) ) * aNormal * this->diffusion_coefficient();
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Spalart_Allmaras_Turbulence::eval_dTractiondDOF(
                const moris::Cell< MSI::Dof_Type >& aDofTypes,
                const Matrix< DDRMat >&             aNormal )
        {
            // get the dof type as a uint
            const uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            const uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the der FI
            Field_Interpolator* tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init the derivative of the traction wrt to dof
            mdTractiondDof( tDofIndex ).set_size( 1, tFIDer->get_number_of_space_time_coefficients() );

            // get the residual dof FI (here viscosity)
            Field_Interpolator* tFIViscosity =
                    mFIManager->get_field_interpolators_for_type( mDofViscosity );

            // add contribution to dtractiondu
            mdTractiondDof( tDofIndex ) =
                    trans( aNormal ) * tFIViscosity->gradx( 1 ) * this->ddiffusioncoeffdu( aDofTypes );

            // if derivative dof type is viscosity dof type
            if ( aDofTypes( 0 ) == mDofViscosity )
            {
                // add contribution to dtractiondu
                mdTractiondDof( tDofIndex ) +=
                        this->diffusion_coefficient() * trans( aNormal ) * tFIViscosity->dnNdxn( 1 );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Spalart_Allmaras_Turbulence::eval_dTestTractiondDOF(
                const moris::Cell< MSI::Dof_Type >& aDofTypes,
                const Matrix< DDRMat >&             aNormal,
                const moris::Cell< MSI::Dof_Type >& aTestDofTypes )
        {
            // get test dof type index
            const uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

            // get the dof type index
            const uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the test dof FI
            Field_Interpolator* tFITest =
                    mFIManager->get_field_interpolators_for_type( aTestDofTypes( 0 ) );

            // get the derivative dof FI
            Field_Interpolator* tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init the dTestTractiondDof
            mdTestTractiondDof( tTestDofIndex )( tDofIndex ).set_size( tFITest->get_number_of_space_time_coefficients(), tFIDer->get_number_of_space_time_coefficients() );

            // get the modified viscosity dof FI
            Field_Interpolator* tFIModViscosity =
                    mFIManager->get_field_interpolators_for_type( mDofViscosity );

            // get modified viscosity value
            real tModViscosity = tFIModViscosity->val()( 0 );

            if ( ( tModViscosity >= 0.0 ) && ( aTestDofTypes( 0 ) == mDofViscosity ) )
            {
                // add contribution
                mdTestTractiondDof( tTestDofIndex )( tDofIndex ) =
                        trans( tFIModViscosity->dnNdxn( 1 ) ) * aNormal * this->ddiffusioncoeffdu( aDofTypes );

                // if derivative dof type is viscosity dof type
                if ( aDofTypes( 0 ) == mDofViscosity )
                {
                    // add contribution
                    mdTestTractiondDof( tTestDofIndex )( tDofIndex ) +=
                            trans( this->ddiffusioncoeffdu( aTestDofTypes ) ) * trans( aNormal ) * tFIModViscosity->dnNdxn( 1 );
                }

                // FIXME assumed that second order derivative of diffusion coeff is zero
            }
            else
            {
                // FIXME compute the test traction derivative by FD
                this->eval_dtesttractiondu_FD(
                        aDofTypes,
                        aTestDofTypes,
                        mdTestTractiondDof( tTestDofIndex )( tDofIndex ),
                        1e-6,
                        aNormal,
                        fem::FDScheme_Type::POINT_3_CENTRAL );
            }
        }

        //------------------------------------------------------------------------------

        void
        CM_Spalart_Allmaras_Turbulence::eval_production_coefficient()
        {
            // get the viscosity FI
            Field_Interpolator* tFIModViscosity =
                    mFIManager->get_field_interpolators_for_type( mDofViscosity );

            // get the viscosity value
            real tModViscosity = tFIModViscosity->val()( 0 );

            // if viscosity is positive or zero
            if ( tModViscosity >= 0.0 )
            {
                // compute production coefficient
                mProductionCoeff = { { mCb1 * ( 1.0 - this->ft2() ) * this->stilde() } };
            }
            // if viscosity is negative
            else
            {
                // compute production coefficient
                mProductionCoeff = { { mCb1 * ( 1.0 - mCt3 ) * this->s() } };
            }
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Spalart_Allmaras_Turbulence::production_coefficient(
                enum CM_Function_Type aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Spalart_Allmaras_Turbulence::production_coefficient - Only DEFAULT CM function type known in base class." );

            // if the production coefficient was not evaluated
            if ( mProductionCoeffEval )
            {
                // evaluate the turbulent dynamic viscosity
                this->eval_production_coefficient();

                // set bool for evaluation
                mProductionCoeffEval = false;
            }

            // return the production coefficient value
            return mProductionCoeff;
        }

        //------------------------------------------------------------------------------

        void
        CM_Spalart_Allmaras_Turbulence::eval_dproductioncoeffdu(
                const moris::Cell< MSI::Dof_Type >& aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get derivative dof type FI
            Field_Interpolator* tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // initialize the matrix for dEffConddu
            mdProductionCoeffdu( tDofIndex ).set_size( 1, tFIDer->get_number_of_space_time_coefficients() );

            // get the viscosity FI
            Field_Interpolator* tFIModViscosity =
                    mFIManager->get_field_interpolators_for_type( mDofViscosity );

            // get the viscosity value
            real tModViscosity = tFIModViscosity->val()( 0 );

            // if viscosity is positive or zero
            if ( tModViscosity >= 0.0 )
            {
                // compute dproductiondu
                mdProductionCoeffdu( tDofIndex ) = mCb1 * ( 1.0 - this->ft2() ) * this->dstildedu( aDofTypes );

                // if contribution from ft2
                if ( mUseFt2 )
                {
                    mdProductionCoeffdu( tDofIndex ) -= mCb1 * this->stilde() * this->dft2du( aDofTypes );
                }
            }
            // if viscosity is negative
            else
            {
                // compute dproductiondu
                mdProductionCoeffdu( tDofIndex ) = mCb1 * ( 1.0 - mCt3 ) * this->dsdu( aDofTypes );
            }

            MORIS_ASSERT( isfinite( mdProductionCoeffdu( tDofIndex ) ),
                    "CM_Spalart_Allmaras_Turbulence::eval_dproductioncoeffdu - mdProductionCoeffdu contains NAN or INF, exiting!" );
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Spalart_Allmaras_Turbulence::dproductioncoeffdu(
                const moris::Cell< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type               aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Spalart_Allmaras_Turbulence::dproductioncoeffdu - Only DEFAULT CM function type known in base class." );

            // if aDofType is not an active dof type for the CM
            MORIS_ERROR(
                    this->check_dof_dependency( aDofType ),
                    "CM_Spalart_Allmaras_Turbulence::dproductioncoeffdu - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if ( mdProductionCoeffduEval( tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_dproductioncoeffdu( aDofType );

                // set bool for evaluation
                mdProductionCoeffduEval( tDofIndex ) = false;
            }

            // return the derivative
            return mdProductionCoeffdu( tDofIndex );
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Spalart_Allmaras_Turbulence::production_term(
                enum CM_Function_Type aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Spalart_Allmaras_Turbulence::production_term - Only DEFAULT CM function type known in base class." );

            // if the production term was not evaluated
            if ( mProductionTermEval )
            {
                // evaluate the production term
                this->eval_production_term();

                // set bool for evaluation
                mProductionTermEval = false;
            }
            // return the production term
            return mProductionTerm;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Spalart_Allmaras_Turbulence::eval_production_term()
        {
            // get the viscosity FI
            Field_Interpolator* tFIModViscosity =
                    mFIManager->get_field_interpolators_for_type( mDofViscosity );

            // get the viscosity value
            real tModViscosity = tFIModViscosity->val()( 0 );

            // compute production term
            mProductionTerm = this->production_coefficient() * tModViscosity;
        }

        //------------------------------------------------------------------------------

        void
        CM_Spalart_Allmaras_Turbulence::eval_dproductiontermdu(
                const moris::Cell< MSI::Dof_Type >& aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get derivative dof type FI
            Field_Interpolator* tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // initialize the matrix for dProductionTermdu
            mdProductionTermdu( tDofIndex ).set_size( 1, tFIDer->get_number_of_space_time_coefficients() );

            // get the viscosity FI
            Field_Interpolator* tFIModViscosity =
                    mFIManager->get_field_interpolators_for_type( mDofViscosity );

            // get the viscosity value
            real tModViscosity = tFIModViscosity->val()( 0 );

            // compute the derivative of the production coefficient
            mdProductionTermdu( tDofIndex ) =
                    this->dproductioncoeffdu( aDofTypes ) * tModViscosity;

            // if derivative dof type is viscosity
            if ( aDofTypes( 0 ) == mDofViscosity )
            {
                // add contribution to derivative
                mdProductionTermdu( tDofIndex ) +=
                        this->production_coefficient()( 0 ) * tFIModViscosity->N();
            }

            MORIS_ASSERT( isfinite( mdProductionTermdu( tDofIndex ) ),
                    "CM_Spalart_Allmaras_Turbulence::eval_dproductiontermdu - mdProductionTermdu contains NAN or INF, exiting!" );
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Spalart_Allmaras_Turbulence::dproductiontermdu(
                const moris::Cell< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type               aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Spalart_Allmaras_Turbulence::dproductiontermdu - Only DEFAULT CM function type known in base class." );

            // if aDofType is not an active dof type for the CM
            MORIS_ERROR(
                    this->check_dof_dependency( aDofType ),
                    "CM_Spalart_Allmaras_Turbulence::dproductiontermdu - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if ( mdProductionTermduEval( tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_dproductiontermdu( aDofType );

                // set bool for evaluation
                mdProductionTermduEval( tDofIndex ) = false;
            }

            // return the derivative
            return mdProductionTermdu( tDofIndex );
        }

        //------------------------------------------------------------------------------

        void
        CM_Spalart_Allmaras_Turbulence::eval_wall_destruction_coefficient()
        {
            // get the viscosity FI
            Field_Interpolator* tFIModViscosity =
                    mFIManager->get_field_interpolators_for_type( mDofViscosity );

            // get the viscosity value
            real tModViscosity = tFIModViscosity->val()( 0 );

            // get the wall distance value
            real tWallDistance = mPropWallDistance->val()( 0 );

            // threshold wall distance
            real tWallDistance2 = clip_value( std::pow( tWallDistance, 2.0 ), mEpsilon );

            // if viscosity is positive or zero
            if ( tModViscosity >= 0.0 )
            {
                // compute wall destruction coefficient
                mWallDestructionCoeff = { { ( mCw1 * this->fw() - mCb1 * this->ft2() / std::pow( mKappa, 2.0 ) ) * tModViscosity / tWallDistance2 } };

                // xxx to be removed
                clip_value( std::pow( tWallDistance, 2.0 ), mEpsilon, ( mCw1 * this->fw() - mCb1 * this->ft2() / std::pow( mKappa, 2.0 ) ) * tModViscosity );
            }
            // if viscosity is negative
            else
            {
                // compute wall destruction coefficient
                mWallDestructionCoeff = { { -mCw1 * mAlpha * tModViscosity / tWallDistance2 } };

                // xxx to be removed
                clip_value( std::pow( tWallDistance, 2.0 ), mEpsilon, -mCw1 * mAlpha * tModViscosity );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Spalart_Allmaras_Turbulence::wall_destruction_coefficient(
                enum CM_Function_Type aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Spalart_Allmaras_Turbulence::wall_destruction_coefficient - Only DEFAULT CM function type known in base class." );

            // if the wall destruction coefficient was not evaluated
            if ( mWallDestructionCoeffEval )
            {
                // evaluate the wall destruction coefficient
                this->eval_wall_destruction_coefficient();

                // set bool for evaluation
                mWallDestructionCoeffEval = false;
            }
            // return the wall destruction coefficient
            return mWallDestructionCoeff;
        }

        //------------------------------------------------------------------------------

        void
        CM_Spalart_Allmaras_Turbulence::eval_dwalldestructioncoeffdu(
                const moris::Cell< MSI::Dof_Type >& aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get derivative dof type FI
            Field_Interpolator* tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // initialize the matrix for dwalldestructioncoeffdu
            mdWallDestructionCoeffdu( tDofIndex ).set_size( 1, tFIDer->get_number_of_space_time_coefficients() );

            // get the viscosity FI
            Field_Interpolator* tFIModViscosity =
                    mFIManager->get_field_interpolators_for_type( mDofViscosity );

            // get the viscosity value
            real tModViscosity = tFIModViscosity->val()( 0 );

            // get the wall distance value
            real tWallDistance = mPropWallDistance->val()( 0 );

            // threshold wall distance
            real tWallDistance2 = clip_value( std::pow( tWallDistance, 2.0 ), mEpsilon );

            // if viscosity is positive or zero
            if ( tModViscosity >= 0.0 )
            {
                // add contribution to dwalldestructiondu
                mdWallDestructionCoeffdu( tDofIndex ) = mCw1 * this->dfwdu( aDofTypes ) * tModViscosity / tWallDistance2;

                // xxx to be removed
                clip_value( std::pow( tWallDistance, 2.0 ), mEpsilon, mCw1 * tModViscosity );

                MORIS_ASSERT( isfinite( mdWallDestructionCoeffdu( tDofIndex ) ),
                        "CM_Spalart_Allmaras_Turbulence::eval_dwalldestructioncoeffdu - mdWallDestructionCoeffdu contains NAN or INF (1)" );

                // if contribution from ft2
                if ( mUseFt2 )
                {
                    mdWallDestructionCoeffdu( tDofIndex ) -= mCb1 * this->dft2du( aDofTypes ) * tModViscosity / std::pow( mKappa, 2.0 ) / tWallDistance2;

                    // xxx to be removed
                    clip_value( std::pow( tWallDistance, 2.0 ), mEpsilon, mCb1 * tModViscosity / std::pow( mKappa, 2.0 ) );

                    MORIS_ASSERT( isfinite( mdWallDestructionCoeffdu( tDofIndex ) ),
                            "CM_Spalart_Allmaras_Turbulence::eval_dwalldestructioncoeffdu - mdWallDestructionCoeffdu contains NAN or INF (2)" );
                }

                // if derivative dof type is viscosity
                if ( aDofTypes( 0 ) == mDofViscosity )
                {
                    // add contribution to dwalldestructiondu
                    mdWallDestructionCoeffdu( tDofIndex ) += mCw1 * this->fw() * tFIModViscosity->N() / tWallDistance2;

                    // xxx to be removed
                    clip_value( std::pow( tWallDistance, 2.0 ), mEpsilon, mCw1 * this->fw() );

                    MORIS_ASSERT( isfinite( mdWallDestructionCoeffdu( tDofIndex ) ),
                            "CM_Spalart_Allmaras_Turbulence::eval_dwalldestructioncoeffdu - mdWallDestructionCoeffdu contains NAN or INF (3)" );

                    // if contribution from ft2
                    if ( mUseFt2 )
                    {
                        mdWallDestructionCoeffdu( tDofIndex ) -= mCb1 * this->ft2() * tFIModViscosity->N() / std::pow( mKappa, 2.0 ) / tWallDistance2;

                        // xxx to be removed
                        clip_value( std::pow( tWallDistance, 2.0 ), mEpsilon, mCb1 * this->ft2() / std::pow( mKappa, 2.0 ) );

                        MORIS_ASSERT( isfinite( mdWallDestructionCoeffdu( tDofIndex ) ),
                                "CM_Spalart_Allmaras_Turbulence::eval_dwalldestructioncoeffdu - mdWallDestructionCoeffdu contains NAN or INF (4)" );
                    }
                }

                // if wall distance depends on derivative dof type
                if ( ( mPropWallDistance->check_dof_dependency( aDofTypes ) ) && ( tWallDistance2 > mEpsilon ) )
                {
                    // threshold denominator
                    real tWallDistance3 = clip_value( std::pow( tWallDistance, 3.0 ), mEpsilonDeriv );

                    // add contribution to dwalldestructiondu
                    mdWallDestructionCoeffdu( tDofIndex ) -= 2.0 * mCw1 * this->fw() * tModViscosity * mPropWallDistance->dPropdDOF( aDofTypes ) / tWallDistance3;

                    // xxx to be removed
                    clip_value( std::pow( tWallDistance, 3.0 ), mEpsilonDeriv, 2.0 * mCw1 * this->fw() * tModViscosity );

                    MORIS_ASSERT( isfinite( mdWallDestructionCoeffdu( tDofIndex ) ),
                            "CM_Spalart_Allmaras_Turbulence::eval_dwalldestructioncoeffdu - mdWallDestructionCoeffdu contains NAN or INF (5)" );

                    // if contribution from ft2
                    if ( mUseFt2 )
                    {
                        mdWallDestructionCoeffdu( tDofIndex ) +=
                                2.0 * mCb1 * this->ft2() / std::pow( mKappa, 2.0 ) * tModViscosity * mPropWallDistance->dPropdDOF( aDofTypes ) / tWallDistance3;

                        MORIS_ASSERT( isfinite( mdWallDestructionCoeffdu( tDofIndex ) ),
                                "CM_Spalart_Allmaras_Turbulence::eval_dwalldestructioncoeffdu - mdWallDestructionCoeffdu contains NAN or INF (6)" );
                    }
                }
            }
            // if viscosity is negative
            else
            {
                // if derivative dof type is viscosity
                if ( aDofTypes( 0 ) == mDofViscosity )
                {
                    // add contribution to dwalldestructiondu
                    mdWallDestructionCoeffdu( tDofIndex ) =
                            -mCw1 * mAlpha * tFIModViscosity->N() / tWallDistance2;
                }
                else
                {
                    // fill with zeros
                    mdWallDestructionCoeffdu( tDofIndex ).fill( 0.0 );
                }

                // if wall distance depends on derivative dof type
                if ( ( mPropWallDistance->check_dof_dependency( aDofTypes ) ) && ( tWallDistance2 > mEpsilon ) )
                {
                    // threshold denominator
                    real tWallDistance3 = clip_value( std::pow( tWallDistance, 3.0 ), mEpsilonDeriv );

                    // add contribution to dwalldestructiondu
                    mdWallDestructionCoeffdu( tDofIndex ) +=
                            2.0 * mCw1 * mAlpha * tModViscosity * mPropWallDistance->dPropdDOF( aDofTypes ) / tWallDistance3;

                    // xxx to be removed
                    clip_value( std::pow( tWallDistance, 3.0 ), mEpsilonDeriv, 2.0 * mCw1 * mAlpha * tModViscosity );

                    MORIS_ASSERT( isfinite( mdWallDestructionCoeffdu( tDofIndex ) ),
                            "CM_Spalart_Allmaras_Turbulence::eval_dwalldestructioncoeffdu - mdWallDestructionCoeffdu contains NAN or INF (7)" );
                }
            }

            MORIS_ASSERT( isfinite( mdWallDestructionCoeffdu( tDofIndex ) ),
                    "CM_Spalart_Allmaras_Turbulence::eval_dwalldestructioncoeffdu - mdWallDestructionCoeffdu contains NAN or INF, exiting!" );
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Spalart_Allmaras_Turbulence::dwalldestructioncoeffdu(
                const moris::Cell< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type               aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Spalart_Allmaras_Turbulence::dwalldestructioncoeffdu - Only DEFAULT CM function type known in base class." );

            // if aDofType is not an active dof type for the CM
            MORIS_ERROR(
                    this->check_dof_dependency( aDofType ),
                    "CM_Spalart_Allmaras_Turbulence::dwalldestructioncoeffdu - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if ( mdWallDestructionCoeffduEval( tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_dwalldestructioncoeffdu( aDofType );

                // set bool for evaluation
                mdWallDestructionCoeffduEval( tDofIndex ) = false;
            }

            // return the derivative
            return mdWallDestructionCoeffdu( tDofIndex );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Spalart_Allmaras_Turbulence::eval_wall_destruction_term()
        {
            // get the viscosity FI
            Field_Interpolator* tFIModViscosity =
                    mFIManager->get_field_interpolators_for_type( mDofViscosity );

            // get the viscosity value
            real tModViscosity = tFIModViscosity->val()( 0 );

            // compute wall destruction term
            mWallDestructionTerm = this->wall_destruction_coefficient() * tModViscosity;
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Spalart_Allmaras_Turbulence::wall_destruction_term(
                enum CM_Function_Type aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Spalart_Allmaras_Turbulence::wall_destruction_term - Only DEFAULT CM function type known in base class." );

            // if the wall destruction term was not evaluated
            if ( mWallDestructionTermEval )
            {
                // evaluate the wall destruction term
                this->eval_wall_destruction_term();

                // set bool for evaluation
                mWallDestructionTermEval = false;
            }
            // return the wall destruction term
            return mWallDestructionTerm;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Spalart_Allmaras_Turbulence::eval_dwalldestructiontermdu(
                const moris::Cell< MSI::Dof_Type >& aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get derivative dof type FI
            Field_Interpolator* tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // initialize the matrix for dWallDestructionTermdu
            mdWallDestructionTermdu( tDofIndex ).set_size( 1, tFIDer->get_number_of_space_time_coefficients() );

            // get the viscosity FI
            Field_Interpolator* tFIModViscosity =
                    mFIManager->get_field_interpolators_for_type( mDofViscosity );

            // get the viscosity value
            real tModViscosity = tFIModViscosity->val()( 0 );

            // compute the derivative of the wall destruction coefficient
            mdWallDestructionTermdu( tDofIndex ) =
                    this->dwalldestructioncoeffdu( aDofTypes ) * tModViscosity;

            // if derivative dof type is viscosity
            if ( aDofTypes( 0 ) == mDofViscosity )
            {
                // add contribution to derivative
                mdWallDestructionTermdu( tDofIndex ) +=
                        this->wall_destruction_coefficient()( 0 ) * tFIModViscosity->N();
            }

            MORIS_ASSERT( isfinite( mdWallDestructionTermdu( tDofIndex ) ),
                    "CM_Spalart_Allmaras_Turbulence::eval_dwalldestructiontermdu - mdWallDestructionTermdu contains NAN or INF, exiting!" );
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Spalart_Allmaras_Turbulence::dwalldestructiontermdu(
                const moris::Cell< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type               aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Spalart_Allmaras_Turbulence::dwalldestructiontermdu - Only DEFAULT CM function type known in base class." );

            // if aDofType is not an active dof type for the CM
            MORIS_ERROR(
                    this->check_dof_dependency( aDofType ),
                    "CM_Spalart_Allmaras_Turbulence::dwalldestructiontermdu - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if ( mdWallDestructionTermduEval( tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_dwalldestructiontermdu( aDofType );

                // set bool for evaluation
                mdWallDestructionTermduEval( tDofIndex ) = false;
            }

            // return the derivative
            return mdWallDestructionTermdu( tDofIndex );
        }

        //------------------------------------------------------------------------------

        void
        CM_Spalart_Allmaras_Turbulence::eval_diffusion_coefficient()
        {
            // get the viscosity FI
            Field_Interpolator* tFIModViscosity =
                    mFIManager->get_field_interpolators_for_type( mDofViscosity );

            // get the viscosity value
            real tModViscosity = tFIModViscosity->val()( 0 );

            // get the fluid kinematic viscosity value
            real tKinViscosity = mPropKinViscosity->val()( 0 );

            // if viscosity is positive or zero
            if ( tModViscosity >= 0.0 )
            {
                // compute diffusion term
                mDiffusionCoeff = { { ( tKinViscosity + tModViscosity ) / mSigma } };
            }
            // if viscosity is negative
            else
            {
                // compute diffusion term
                mDiffusionCoeff = { { ( tKinViscosity + tModViscosity * this->fn() ) / mSigma } };
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Spalart_Allmaras_Turbulence::diffusion_coefficient(
                enum CM_Function_Type aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Spalart_Allmaras_Turbulence::diffusion_coefficient - Only DEFAULT CM function type known in base class." );

            // if the diffusion coefficient was not evaluated
            if ( mDiffusionCoeffEval )
            {
                // evaluate the diffusion coefficient
                this->eval_diffusion_coefficient();

                // set bool for evaluation
                mDiffusionCoeffEval = false;
            }
            // return the diffusion coefficient
            return mDiffusionCoeff;
        }

        //------------------------------------------------------------------------------

        void
        CM_Spalart_Allmaras_Turbulence::eval_ddiffusioncoeffdu(
                const moris::Cell< MSI::Dof_Type >& aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get derivative dof type FI
            Field_Interpolator* tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // initialize the matrix for dEffConddu
            mdDiffusionCoeffdu( tDofIndex ).set_size( 1, tFIDer->get_number_of_space_time_coefficients() );

            // get the viscosity FI
            Field_Interpolator* tFIModViscosity =
                    mFIManager->get_field_interpolators_for_type( mDofViscosity );

            // get the viscosity value
            real tModViscosity = tFIModViscosity->val()( 0 );

            // if viscosity is positive or zero
            if ( tModViscosity >= 0.0 )
            {
                // if derivative dof type is viscosity
                if ( aDofTypes( 0 ) == mDofViscosity )
                {
                    // add contribution to ddiffusiondu
                    mdDiffusionCoeffdu( tDofIndex ) = tFIModViscosity->N() / mSigma;
                }
                else
                {
                    mdDiffusionCoeffdu( tDofIndex ).fill( 0.0 );
                }

                // if kinematic viscosity depends on derivative dof type
                if ( mPropKinViscosity->check_dof_dependency( aDofTypes ) )
                {
                    // add contribution to ddiffusiondu
                    mdDiffusionCoeffdu( tDofIndex ) +=
                            mPropKinViscosity->dPropdDOF( aDofTypes ) / mSigma;
                }
            }
            // if viscosity is negative
            else
            {
                // add contribution from fn to ddiffusiondu
                mdDiffusionCoeffdu( tDofIndex ) = tModViscosity * this->dfndu( aDofTypes ) / mSigma;

                // if derivative dof type is viscosity
                if ( aDofTypes( 0 ) == mDofViscosity )
                {
                    // add contribution to ddiffusiondu
                    mdDiffusionCoeffdu( tDofIndex ) += this->fn() * tFIModViscosity->N() / mSigma;
                }

                // if kinematic viscosity depends on derivative dof type
                if ( mPropKinViscosity->check_dof_dependency( aDofTypes ) )
                {
                    // add contribution to ddiffusiondu
                    mdDiffusionCoeffdu( tDofIndex ) += mPropKinViscosity->dPropdDOF( aDofTypes ) / mSigma;
                }
            }

            MORIS_ASSERT( isfinite( mdDiffusionCoeffdu( tDofIndex ) ),
                    "CM_Spalart_Allmaras_Turbulence::eval_ddiffusioncoeffdu - mdDiffusionCoeffdu contains NAN or INF, exiting!" );
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Spalart_Allmaras_Turbulence::ddiffusioncoeffdu(
                const moris::Cell< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type               aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Spalart_Allmaras_Turbulence::ddiffusioncoeffdu - Only DEFAULT CM function type known in base class." );

            // if aDofType is not an active dof type for the CM
            MORIS_ERROR(
                    this->check_dof_dependency( aDofType ),
                    "CM_Spalart_Allmaras_Turbulence::ddiffusioncoeffdu - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if ( mdDiffusionCoeffduEval( tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_ddiffusioncoeffdu( aDofType );

                // set bool for evaluation
                mdDiffusionCoeffduEval( tDofIndex ) = false;
            }

            // return the derivative
            return mdDiffusionCoeffdu( tDofIndex );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Spalart_Allmaras_Turbulence::eval_ddiffusioncoeffdx( uint aOrder )
        {
            // FIXME work only for 1st order
            MORIS_ERROR( aOrder == 1,
                    "CM_Spalart_Allmaras_Turbulence::eval_ddiffusioncoeffdx - Works only for 1st order derivative for now." );

            // set matrix size
            mdDiffusionCoeffdx( aOrder - 1 ).set_size( mSpaceDim, 1 );

            // get the viscosity FI
            Field_Interpolator* tFIModViscosity =
                    mFIManager->get_field_interpolators_for_type( mDofViscosity );

            // get the viscosity value
            real tModViscosity = tFIModViscosity->val()( 0 );

            // if viscosity is positive or zero
            if ( tModViscosity >= 0.0 )
            {
                // compute diffusion term
                mdDiffusionCoeffdx( aOrder - 1 ) = tFIModViscosity->gradx( 1 ) / mSigma;

                // if kinematic viscosity depends on space
                if ( mPropKinViscosity->check_space_dependency( 1 ) )
                {
                    mdDiffusionCoeffdx( aOrder - 1 ) += mPropKinViscosity->dnPropdxn( 1 ) / mSigma;
                }
            }
            // if viscosity is negative
            else
            {
                // compute diffusion term
                mdDiffusionCoeffdx( aOrder - 1 ) =
                        ( tFIModViscosity->gradx( 1 ) * this->fn() + tFIModViscosity->val()( 0 ) * this->dfndx( 1 ) ) / mSigma;

                // if kinematic viscosity depends on space
                if ( mPropKinViscosity->check_space_dependency( 1 ) )
                {
                    mdDiffusionCoeffdx( aOrder - 1 ) += mPropKinViscosity->dnPropdxn( 1 ) / mSigma;
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Spalart_Allmaras_Turbulence::ddiffusioncoeffdx(
                uint                  aOrder,
                enum CM_Function_Type aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Spalart_Allmaras_Turbulence::ddiffusioncoeffdx - Only DEFAULT CM function type known in base class." );

            MORIS_ERROR( aOrder == 1,
                    "CM_Spalart_Allmaras_Turbulence::ddiffusioncoeffdx - Works only for 1st order derivative for now." );

            // if the derivative has not been evaluated yet
            if ( mdDiffusionCoeffdxEval( aOrder - 1 ) )
            {
                // evaluate the derivative
                this->eval_ddiffusioncoeffdx( aOrder );

                // set bool for evaluation
                mdDiffusionCoeffdxEval( aOrder - 1 ) = false;
            }

            // return the derivative
            return mdDiffusionCoeffdx( aOrder - 1 );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Spalart_Allmaras_Turbulence::eval_ddiffusioncoeffdxdu(
                const moris::Cell< MSI::Dof_Type >& aDofTypes,
                uint                                aOrder )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get derivative dof type FI
            Field_Interpolator* tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set matrix size
            mdDiffusionCoeffdxdu( aOrder - 1 )( tDofIndex ).set_size( mSpaceDim, tFIDer->get_number_of_space_time_coefficients() );

            // get the viscosity FI
            Field_Interpolator* tFIModViscosity =
                    mFIManager->get_field_interpolators_for_type( mDofViscosity );

            // get the viscosity value
            real tModViscosity = tFIModViscosity->val()( 0 );

            // if viscosity is positive or zero
            if ( tModViscosity >= 0.0 )
            {
                // if derivative dof type is viscosity
                if ( aDofTypes( 0 ) == mDofViscosity )
                {
                    // compute diffusion term
                    mdDiffusionCoeffdxdu( aOrder - 1 )( tDofIndex ) =
                            tFIModViscosity->dnNdxn( 1 ) / mSigma;

                    //                    //FIXME missing term aPropKinViscosity->dPropdxdu()
                    //                    // if kinematic viscosity depends on space
                    //                    if( mPropKinViscosity->check_space_dependency( 1 ) )
                    //                    {
                    //                         mdDiffusionCoeffdxdu( aOrder - 1 )( tDofIndex ) += mPropKinViscosity->dPropdxdu / mSigma;
                    //                    }
                }
                else
                {
                    // fill with zeros
                    mdDiffusionCoeffdxdu( aOrder - 1 )( tDofIndex ).fill( 0.0 );
                }
            }
            // if viscosity is negative
            else
            {
                // compute diffusion term
                mdDiffusionCoeffdxdu( aOrder - 1 )( tDofIndex ) =
                        ( tFIModViscosity->gradx( 1 ) * this->dfndu( aDofTypes ) + tFIModViscosity->val()( 0 ) * this->dfndxdu( aDofTypes, 1 ) ) / mSigma;

                // FIXME missing term aPropKinViscosity->dPropdxdu()
                //                     // if kinematic viscosity depends on space
                //                     if( aPropKinViscosity->check_space_dependency( 1 ) )
                //                     {
                //                         mdDiffusionCoeffdxdu( aOrder - 1 )( tDofIndex ) +=
                //                                 mPropKinViscosity->dPropdxdu() / mSigma;
                //                     }

                // if derivative dof type is viscosity
                if ( aDofTypes( 0 ) == mDofViscosity )
                {
                    mdDiffusionCoeffdxdu( aOrder - 1 )( tDofIndex ) +=
                            ( tFIModViscosity->dnNdxn( 1 ) * this->fn() + this->dfndx( 1 ) * tFIModViscosity->N() ) / mSigma;
                }
            }

            MORIS_ASSERT( isfinite( mdDiffusionCoeffdxdu( aOrder - 1 )( tDofIndex ) ),
                    "CM_Spalart_Allmaras_Turbulence::eval_ddiffusioncoeffdxdu - mdDiffusionCoeffdxdu contains NAN or INF, exiting!" );
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Spalart_Allmaras_Turbulence::ddiffusioncoeffdxdu(
                const moris::Cell< MSI::Dof_Type >& aDofType,
                uint                                aOrder,
                enum CM_Function_Type               aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Spalart_Allmaras_Turbulence::ddiffusioncoeffdxdu - Only DEFAULT CM function type known in base class." );

            MORIS_ERROR( aOrder == 1,
                    "CM_Spalart_Allmaras_Turbulence::ddiffusioncoeffdxdu - Works only for 1st order derivative for now." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if ( mdDiffusionCoeffdxduEval( aOrder - 1, tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_ddiffusioncoeffdxdu( aDofType, aOrder );

                // set bool for evaluation
                mdDiffusionCoeffdxduEval( aOrder - 1, tDofIndex ) = false;
            }

            // return the derivative
            return mdDiffusionCoeffdxdu( aOrder - 1 )( tDofIndex );
        }

        //------------------------------------------------------------------------------

        void
        CM_Spalart_Allmaras_Turbulence::eval_modified_velocity()
        {
            // get the viscosity FI
            Field_Interpolator* tFIModViscosity =
                    mFIManager->get_field_interpolators_for_type( mDofViscosity );

            // get the velocity FI
            Field_Interpolator* tFIVelocity =
                    mFIManager->get_field_interpolators_for_type( mDofVelocity );

            // compute modified velocity
            mModVelocity = tFIVelocity->val() - mCb2 * tFIModViscosity->gradx( 1 ) / mSigma;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Spalart_Allmaras_Turbulence::modified_velocity(
                enum CM_Function_Type aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Spalart_Allmaras_Turbulence::modified_velocity - Only DEFAULT CM function type known in base class." );

            // if the modified velocity was not evaluated
            if ( mModVelocityEval )
            {
                // evaluate the modified velocity
                this->eval_modified_velocity();

                // set bool for evaluation
                mModVelocityEval = false;
            }
            // return the modified velocity value
            return mModVelocity;
        }

        //------------------------------------------------------------------------------

        void
        CM_Spalart_Allmaras_Turbulence::eval_dmodvelocitydu(
                const moris::Cell< MSI::Dof_Type >& aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get derivative dof type FI
            Field_Interpolator* tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // initialize the matrix for dEffConddu
            mdModVelocitydu( tDofIndex ).set_size( mSpaceDim, tFIDer->get_number_of_space_time_coefficients() );

            // get the velocity FI
            Field_Interpolator* tFIVelocity =
                    mFIManager->get_field_interpolators_for_type( mDofVelocity );

            // get the viscosity FI
            Field_Interpolator* tFIModViscosity =
                    mFIManager->get_field_interpolators_for_type( mDofViscosity );

            // if dof type is velocity
            if ( aDofTypes( 0 ) == mDofVelocity )
            {
                // add contribution to mdPPdMasterDof
                mdModVelocitydu( tDofIndex ) = tFIVelocity->N();
            }
            // if dof type is viscosity
            else if ( aDofTypes( 0 ) == mDofViscosity )
            {
                // add contribution to mdPPdMasterDof
                mdModVelocitydu( tDofIndex ) = -mCb2 * tFIModViscosity->dnNdxn( 1 ) / mSigma;
            }
            else
            {
                mdModVelocitydu( tDofIndex ).fill( 0.0 );
            }

            MORIS_ASSERT( isfinite( mdModVelocitydu( tDofIndex ) ),
                    "CM_Spalart_Allmaras_Turbulence::eval_dmodvelocitydu - mdModVelocitydu contains NAN or INF, exiting!" );
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Spalart_Allmaras_Turbulence::dmodvelocitydu(
                const moris::Cell< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type               aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Spalart_Allmaras_Turbulence::dmodvelocitydu - Only DEFAULT CM function type known in base class." );

            // if aDofType is not an active dof type for the CM
            MORIS_ERROR(
                    this->check_dof_dependency( aDofType ),
                    "CM_Spalart_Allmaras_Turbulence::dmodvelocitydu - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if ( mdModVelocityduEval( tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_dmodvelocitydu( aDofType );

                // set bool for evaluation
                mdModVelocityduEval( tDofIndex ) = false;
            }

            // return the derivative
            return mdModVelocitydu( tDofIndex );
        }

        //------------------------------------------------------------------------------

        void
        CM_Spalart_Allmaras_Turbulence::eval_modified_velocity_linearized()
        {
            // get the viscosity FI
            Field_Interpolator* tFIModViscosity =
                    mFIManager->get_field_interpolators_for_type( mDofViscosity );

            // get the velocity FI
            Field_Interpolator* tFIVelocity =
                    mFIManager->get_field_interpolators_for_type( mDofVelocity );

            // compute modified velocity
            mModVelocityLinearized = tFIVelocity->val() - 2.0 * mCb2 * tFIModViscosity->gradx( 1 ) / mSigma;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Spalart_Allmaras_Turbulence::modified_velocity_linearized(
                enum CM_Function_Type aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Spalart_Allmaras_Turbulence::modified_velocity_linearized - Only DEFAULT CM function type known in base class." );

            // if the modified velocity was not evaluated
            if ( mModVelocityLinearizedEval )
            {
                // evaluate the modified velocity
                this->eval_modified_velocity_linearized();

                // set bool for evaluation
                mModVelocityLinearizedEval = false;
            }
            // return the modified velocity value
            return mModVelocityLinearized;
        }

        //------------------------------------------------------------------------------

        void
        CM_Spalart_Allmaras_Turbulence::eval_dmodvelocitylinearizeddu(
                const moris::Cell< MSI::Dof_Type >& aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get derivative dof type FI
            Field_Interpolator* tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // initialize the matrix for dEffConddu
            mdModVelocityLinearizeddu( tDofIndex ).set_size( mSpaceDim, tFIDer->get_number_of_space_time_coefficients() );

            // get the velocity FI
            Field_Interpolator* tFIVelocity =
                    mFIManager->get_field_interpolators_for_type( mDofVelocity );

            // get the viscosity FI
            Field_Interpolator* tFIModViscosity =
                    mFIManager->get_field_interpolators_for_type( mDofViscosity );

            // if dof type is velocity
            if ( aDofTypes( 0 ) == mDofVelocity )
            {
                // add contribution to mdPPdMasterDof
                mdModVelocityLinearizeddu( tDofIndex ) = tFIVelocity->N();
            }
            // if dof type is viscosity
            else if ( aDofTypes( 0 ) == mDofViscosity )
            {
                // add contribution to mdPPdMasterDof
                mdModVelocityLinearizeddu( tDofIndex ) = -2.0 * mCb2 * tFIModViscosity->dnNdxn( 1 ) / mSigma;
            }
            else
            {
                mdModVelocityLinearizeddu( tDofIndex ).fill( 0.0 );
            }

            MORIS_ASSERT( isfinite( mdModVelocitydu( tDofIndex ) ),
                    "CM_Spalart_Allmaras_Turbulence::eval_dmodvelocitylinearizeddu - mdModVelocitydu contains NAN or INF, exiting!" );
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Spalart_Allmaras_Turbulence::dmodvelocitylinearizeddu(
                const moris::Cell< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type               aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Spalart_Allmaras_Turbulence::dmodvelocitylinearizeddu - Only DEFAULT CM function type known in base class." );

            // if aDofType is not an active dof type for the CM
            MORIS_ERROR(
                    this->check_dof_dependency( aDofType ),
                    "CM_Spalart_Allmaras_Turbulence::dmodvelocitylinearizeddu - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if ( mdModVelocityLinearizedduEval( tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_dmodvelocitylinearizeddu( aDofType );

                // set bool for evaluation
                mdModVelocityLinearizedduEval( tDofIndex ) = false;
            }

            // return the derivative
            return mdModVelocityLinearizeddu( tDofIndex );
        }

        //--------------------------------------------------------------------------------------------------------------

        real
        CM_Spalart_Allmaras_Turbulence::chi(
                enum CM_Function_Type aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Spalart_Allmaras_Turbulence::chi - Only DEFAULT CM function type known in base class." );

            // if the diffusion coefficient was not evaluated
            if ( mChiEval )
            {
                // evaluate chi
                mChi = compute_chi(
                        { mDofViscosity },
                        mFIManager,
                        mPropKinViscosity );

                // set bool for evaluation
                mChiEval = false;
            }
            // return the diffusion coefficient
            return mChi;
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Spalart_Allmaras_Turbulence::dchidu(
                const moris::Cell< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type               aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Fluid_Turbulence::dchidu - Only DEFAULT CM function type known in base class." );

            // if aDofType is not an active dof type for the CM
            MORIS_ERROR(
                    this->check_dof_dependency( aDofType ),
                    "CM_Fluid_Turbulence::dchidu - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if ( mdChiduEval( tDofIndex ) )
            {
                // evaluate the derivative
                compute_dchidu(
                        { mDofViscosity },
                        mFIManager,
                        mPropKinViscosity,
                        aDofType,
                        mdChidu( tDofIndex ) );

                // set bool for evaluation
                mdChiduEval( tDofIndex ) = false;
            }

            // return the derivative
            return mdChidu( tDofIndex );
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Spalart_Allmaras_Turbulence::dchidx(
                uint                  aOrder,
                enum CM_Function_Type aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Fluid_Turbulence::dchidx - Only DEFAULT CM function type known in base class." );

            MORIS_ERROR( aOrder == 1,
                    "CM_Fluid_Turbulence::dchidx - Works only for 1st order derivative for now." );

            // if the derivative has not been evaluated yet
            if ( mdChidxEval( aOrder - 1 ) )
            {
                // evaluate the derivative
                compute_dchidx(
                        { mDofViscosity },
                        mFIManager,
                        mPropKinViscosity,
                        mdChidx( aOrder - 1 ) );

                // set bool for evaluation
                mdChidxEval( aOrder - 1 ) = false;
            }

            // return the derivative
            return mdChidx( aOrder - 1 );
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Spalart_Allmaras_Turbulence::dchidxdu(
                const moris::Cell< MSI::Dof_Type >& aDofType,
                uint                                aOrder,
                enum CM_Function_Type               aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Fluid_Turbulence::dchidxdu - Only DEFAULT CM function type known in base class." );

            MORIS_ERROR( aOrder == 1,
                    "CM_Fluid_Turbulence::dchidxdu - Works only for 1st order derivative for now." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if ( mdChidxduEval( aOrder - 1, tDofIndex ) )
            {
                // evaluate the derivative
                compute_dchidxdu(
                        { mDofViscosity },
                        mFIManager,
                        mPropKinViscosity,
                        aDofType,
                        mdChidxdu( aOrder - 1 )( tDofIndex ) );

                // set bool for evaluation
                mdChidxduEval( aOrder - 1, tDofIndex ) = false;
            }

            // return the derivative
            return mdChidxdu( aOrder - 1 )( tDofIndex );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Spalart_Allmaras_Turbulence::eval_ft2()
        {
            // compute ft2
            if ( mUseFt2 )
            {
                mFt2 = mCt3 * std::exp( -mCt4 * std::pow( this->chi(), 2.0 ) );
            }
            else
            {
                mFt2 = 0.0;
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        real
        CM_Spalart_Allmaras_Turbulence::ft2(
                enum CM_Function_Type aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Spalart_Allmaras_Turbulence::chi - Only DEFAULT CM function type known in base class." );

            // if the diffusion coefficient was not evaluated
            if ( mFt2Eval )
            {
                // evaluate chi
                this->eval_ft2();

                // set bool for evaluation
                mFt2Eval = false;
            }
            // return the diffusion coefficient
            return mFt2;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Spalart_Allmaras_Turbulence::eval_dft2du(
                const moris::Cell< MSI::Dof_Type >& aDofTypes )
        {
            // check that function is only called if ft2 is considered
            MORIS_ERROR( mUseFt2,
                    "CM_Spalart_Allmaras_Turbulence::eval_dft2du - function is called although mUseFt2 is false." );

            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get derivative dof type FI
            Field_Interpolator* tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // initialize the matrix for dEffConddu
            mdFt2du( tDofIndex ).set_size( 1, tFIDer->get_number_of_space_time_coefficients() );

            // compute dft2du
            mdFt2du( tDofIndex ) = -mCt4 * this->ft2() * 2.0 * this->chi() * this->dchidu( aDofTypes );

            MORIS_ASSERT( isfinite( mdFt2du( tDofIndex ) ),
                    "CM_Spalart_Allmaras_Turbulence::eval_dft2du - mdFt2du contains NAN or INF, exiting!" );
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Spalart_Allmaras_Turbulence::dft2du(
                const moris::Cell< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type               aCMFunctionType )
        {
            // check that function is only called if ft2 is considered
            MORIS_ERROR( mUseFt2,
                    "CM_Spalart_Allmaras_Turbulence::dft2du - function is called although mUseFt2 is false." );

            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Spalart_Allmaras_Turbulence::dft2du - Only DEFAULT CM function type known in base class." );

            // if aDofType is not an active dof type for the CM
            MORIS_ERROR(
                    this->check_dof_dependency( aDofType ),
                    "CM_Spalart_Allmaras_Turbulence::dft2du - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if ( mdFt2duEval( tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_dft2du( aDofType );

                // set bool for evaluation
                mdFt2duEval( tDofIndex ) = false;
            }

            // return the derivative
            return mdFt2du( tDofIndex );
        }

        //--------------------------------------------------------------------------------------------------------------

        // FIXME function pointer for 2D and 3D
        void
        CM_Spalart_Allmaras_Turbulence::eval_w()
        {
            // get the velocity FI
            Field_Interpolator* tFIVelocity =
                    mFIManager->get_field_interpolators_for_type( mDofVelocity );

            // get gradient of velocity
            const Matrix< DDRMat >& tGradVelocity = tFIVelocity->gradx( 1 );

            // switch on space dim
            uint tSpaceDim = tFIVelocity->val().numel();
            switch ( tSpaceDim )
            {
                case 2:
                {
                    // init aWij = [ w11 w12 w21 w22]
                    mW.set_size( 4, 1, 0.0 );

                    // compute Wij
                    mW( 1 ) = 0.5 * ( tGradVelocity( 1, 0 ) - tGradVelocity( 0, 1 ) );
                    mW( 2 ) = 0.5 * ( tGradVelocity( 0, 1 ) - tGradVelocity( 1, 0 ) );
                    break;
                }
                case 3:
                {
                    // init aWij = [ w11 w12 w13 w21 w22 w23 w31 w32 w33 ]
                    mW.set_size( 9, 1, 0.0 );

                    // compute Wij
                    mW( 1 ) = 0.5 * ( tGradVelocity( 1, 0 ) - tGradVelocity( 0, 1 ) );
                    mW( 2 ) = 0.5 * ( tGradVelocity( 2, 0 ) - tGradVelocity( 0, 2 ) );
                    mW( 3 ) = 0.5 * ( tGradVelocity( 0, 1 ) - tGradVelocity( 1, 0 ) );
                    mW( 5 ) = 0.5 * ( tGradVelocity( 2, 1 ) - tGradVelocity( 1, 2 ) );
                    mW( 6 ) = 0.5 * ( tGradVelocity( 0, 2 ) - tGradVelocity( 2, 0 ) );
                    mW( 7 ) = 0.5 * ( tGradVelocity( 1, 2 ) - tGradVelocity( 2, 1 ) );
                    break;
                }
                default:
                    MORIS_ERROR( false, "CM_Spalart_Allmaras_Turbulence::eval_w - space dim can only be 2 or 3" );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Spalart_Allmaras_Turbulence::w(
                enum CM_Function_Type aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Spalart_Allmaras_Turbulence::w - Only DEFAULT CM function type known in base class." );

            // if the diffusion coefficient was not evaluated
            if ( mWEval )
            {
                // evaluate the diffusion coefficient
                this->eval_w();

                // set bool for evaluation
                mWEval = false;
            }
            // return the diffusion coefficient
            return mW;
        }

        //--------------------------------------------------------------------------------------------------------------

        // FIXME function pointer for 2D, 3D
        void
        CM_Spalart_Allmaras_Turbulence::eval_dwdu(
                const moris::Cell< MSI::Dof_Type >& aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the der FI
            Field_Interpolator* tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // get the velocity FI
            Field_Interpolator* tFIVelocity =
                    mFIManager->get_field_interpolators_for_type( mDofVelocity );

            // switch on space dim
            uint tSpaceDim = tFIVelocity->val().numel();
            switch ( tSpaceDim )
            {
                case 2:
                {
                    // init aWij = [ w11 w12 w21 w22]
                    mdWdu( tDofIndex ).set_size( 4, tFIDer->get_number_of_space_time_coefficients(), 0.0 );

                    if ( aDofTypes( 0 ) == mDofVelocity )
                    {
                        // get gradient of velocity
                        const Matrix< DDRMat >& tdNdxVelocity = tFIVelocity->dnNdxn( 1 );

                        // get number of bases for displacement
                        uint tNumBases = tFIVelocity->get_number_of_space_time_bases();

                        // compute adwijdu
                        mdWdu( tDofIndex )( { 1, 1 }, { 0, tNumBases - 1 } )             = 0.5 * tdNdxVelocity.get_row( 1 );
                        mdWdu( tDofIndex )( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } ) = -0.5 * tdNdxVelocity.get_row( 0 );
                        mdWdu( tDofIndex )( { 2, 2 }, { 0, tNumBases - 1 } )             = -0.5 * tdNdxVelocity.get_row( 1 );
                        mdWdu( tDofIndex )( { 2, 2 }, { tNumBases, 2 * tNumBases - 1 } ) = 0.5 * tdNdxVelocity.get_row( 0 );
                    }
                    break;
                }
                case 3:
                {
                    // init aWij = [ w11 w12 w13 w21 w22 w23 w31 w32 w33 ]
                    mdWdu( tDofIndex ).set_size( 9, tFIDer->get_number_of_space_time_coefficients(), 0.0 );

                    if ( aDofTypes( 0 ) == mDofVelocity )
                    {
                        // get gradient of velocity
                        const Matrix< DDRMat >& tdNdxVelocity = tFIVelocity->dnNdxn( 1 );

                        // get number of bases for displacement
                        uint tNumBases = tFIVelocity->get_number_of_space_time_bases();

                        // compute adwijdu
                        mdWdu( tDofIndex )( { 1, 1 }, { 0, tNumBases - 1 } )             = 0.5 * tdNdxVelocity.get_row( 1 );
                        mdWdu( tDofIndex )( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } ) = -0.5 * tdNdxVelocity.get_row( 0 );

                        mdWdu( tDofIndex )( { 2, 2 }, { 0, tNumBases - 1 } )                 = 0.5 * tdNdxVelocity.get_row( 2 );
                        mdWdu( tDofIndex )( { 2, 2 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = -0.5 * tdNdxVelocity.get_row( 0 );

                        mdWdu( tDofIndex )( { 3, 3 }, { 0, tNumBases - 1 } )             = -0.5 * tdNdxVelocity.get_row( 1 );
                        mdWdu( tDofIndex )( { 3, 3 }, { tNumBases, 2 * tNumBases - 1 } ) = 0.5 * tdNdxVelocity.get_row( 0 );

                        mdWdu( tDofIndex )( { 5, 5 }, { tNumBases, 2 * tNumBases - 1 } )     = 0.5 * tdNdxVelocity.get_row( 2 );
                        mdWdu( tDofIndex )( { 5, 5 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = -0.5 * tdNdxVelocity.get_row( 1 );

                        mdWdu( tDofIndex )( { 6, 6 }, { 0, tNumBases - 1 } )                 = -0.5 * tdNdxVelocity.get_row( 2 );
                        mdWdu( tDofIndex )( { 6, 6 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = 0.5 * tdNdxVelocity.get_row( 0 );

                        mdWdu( tDofIndex )( { 7, 7 }, { tNumBases, 2 * tNumBases - 1 } )     = -0.5 * tdNdxVelocity.get_row( 2 );
                        mdWdu( tDofIndex )( { 7, 7 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = 0.5 * tdNdxVelocity.get_row( 1 );
                    }
                    break;
                }
                default:
                    MORIS_ERROR( false, "CM_Spalart_Allmaras_Turbulence::dWdu - space dim can only be 2 or 3" );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Spalart_Allmaras_Turbulence::dwdu(
                const moris::Cell< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type               aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Spalart_Allmaras_Turbulence::dwdu - Only DEFAULT CM function type known in base class." );

            // if aDofType is not an active dof type for the CM
            MORIS_ERROR(
                    this->check_dof_dependency( aDofType ),
                    "CM_Spalart_Allmaras_Turbulence::dwdu - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if ( mdWduEval( tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_dwdu( aDofType );

                // set bool for evaluation
                mdWduEval( tDofIndex ) = false;
            }

            // return the derivative
            return mdWdu( tDofIndex );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Spalart_Allmaras_Turbulence::eval_s()
        {
            // compute vorticity in vector form and its Frobenius norm
            real tWW = dot( this->w(), this->w() );

            mS = std::max( std::sqrt( 2.0 * tWW ), mEpsilon );
        }

        //--------------------------------------------------------------------------------------------------------------

        real
        CM_Spalart_Allmaras_Turbulence::s(
                enum CM_Function_Type aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Spalart_Allmaras_Turbulence::s - Only DEFAULT CM function type known in base class." );

            // if not evaluated
            if ( mSEval )
            {
                // evaluate
                this->eval_s();

                // set bool for evaluation
                mSEval = false;
            }
            // return
            return mS;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Spalart_Allmaras_Turbulence::eval_dsdu(
                const moris::Cell< MSI::Dof_Type >& aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get derivative dof type FI
            Field_Interpolator* tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init dsdu
            mdSdu( tDofIndex ).set_size( 1, tFIDer->get_number_of_space_time_coefficients() );

            // if s is greater than threshold
            if ( this->s() > mEpsilon )
            {
                // compute dsdu
                mdSdu( tDofIndex ) = 2.0 * trans( this->w() ) * this->dwdu( aDofTypes ) / this->s();
            }
            else
            {
                mdSdu( tDofIndex ).fill( 0.0 );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Spalart_Allmaras_Turbulence::dsdu(
                const moris::Cell< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type               aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Spalart_Allmaras_Turbulence::dsdu - Only DEFAULT CM function type known in base class." );

            // if aDofType is not an active dof type for the CM
            MORIS_ERROR(
                    this->check_dof_dependency( aDofType ),
                    "CM_Spalart_Allmaras_Turbulence::dsdu - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if ( mdSduEval( tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_dsdu( aDofType );

                // set bool for evaluation
                mdSduEval( tDofIndex ) = false;
            }

            // return the derivative
            return mdSdu( tDofIndex );
        }

        //--------------------------------------------------------------------------------------------------------------

        real
        CM_Spalart_Allmaras_Turbulence::fv1(
                enum CM_Function_Type aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Spalart_Allmaras_Turbulence::fv1 - Only DEFAULT CM function type known in base class." );

            // if not evaluated
            if ( mFv1Eval )
            {
                // evaluate
                mFv1 = compute_fv1(
                        { mDofViscosity },
                        mFIManager,
                        mPropKinViscosity );

                // set bool for evaluation
                mFv1Eval = false;
            }
            // return
            return mFv1;
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Spalart_Allmaras_Turbulence::dfv1du(
                const moris::Cell< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type               aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Spalart_Allmaras_Turbulence::dfv1du - Only DEFAULT CM function type known in base class." );

            // if aDofType is not an active dof type for the CM
            MORIS_ERROR(
                    this->check_dof_dependency( aDofType ),
                    "CM_Spalart_Allmaras_Turbulence::dfv1du - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if ( mdFv1duEval( tDofIndex ) )
            {
                // evaluate the derivative
                compute_dfv1du(
                        { mDofViscosity },
                        mFIManager,
                        mPropKinViscosity,
                        aDofType,
                        mdFv1du( tDofIndex ) );

                // set bool for evaluation
                mdFv1duEval( tDofIndex ) = false;
            }

            // return the derivative
            return mdFv1du( tDofIndex );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Spalart_Allmaras_Turbulence::eval_fv2()
        {
            // threshold deno
            real tDeno = clip_value( 1.0 + this->chi() * this->fv1(), mEpsilon );

            // compute fv2
            mFv2 = 1.0 - this->chi() / tDeno;

            // xxx to be removed
            clip_value( 1.0 + this->chi() * this->fv1(), mEpsilon, this->chi() );
        }

        //--------------------------------------------------------------------------------------------------------------

        real
        CM_Spalart_Allmaras_Turbulence::fv2(
                enum CM_Function_Type aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Spalart_Allmaras_Turbulence::fv2 - Only DEFAULT CM function type known in base class." );

            // if not evaluated
            if ( mFv2Eval )
            {
                // evaluate
                this->eval_fv2();

                // set bool for evaluation
                mFv2Eval = false;
            }
            // return
            return mFv2;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Spalart_Allmaras_Turbulence::eval_dfv2du(
                const moris::Cell< MSI::Dof_Type >& aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get derivative dof type FI
            Field_Interpolator* tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // initialize the matrix for dEffConddu
            mdFv2du( tDofIndex ).set_size( 1, tFIDer->get_number_of_space_time_coefficients() );

            // threshold deno
            real tDeno = clip_value( 1.0 + this->chi() * this->fv1(), mEpsilon );

            if ( std::abs( tDeno ) > mEpsilon )
            {
                // threshold denominator
                real tDeno2 = clip_value( std::pow( tDeno, 2.0 ), mEpsilonDeriv );

                // compute adfv2du
                mdFv2du( tDofIndex ) =
                        ( std::pow( this->chi(), 2.0 ) * this->dfv1du( aDofTypes ) - this->dchidu( aDofTypes ) ) / tDeno2;
            }
            else
            {
                // compute adfv2du
                mdFv2du( tDofIndex ) = -this->dchidu( aDofTypes ) / tDeno;
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Spalart_Allmaras_Turbulence::dfv2du(
                const moris::Cell< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type               aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Spalart_Allmaras_Turbulence::dfv2du - Only DEFAULT CM function type known in base class." );

            // if aDofType is not an active dof type for the CM
            MORIS_ERROR(
                    this->check_dof_dependency( aDofType ),
                    "CM_Spalart_Allmaras_Turbulence::dfv2du - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if ( mdFv2duEval( tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_dfv2du( aDofType );

                // set bool for evaluation
                mdFv2duEval( tDofIndex ) = false;
            }

            // return the derivative
            return mdFv2du( tDofIndex );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Spalart_Allmaras_Turbulence::eval_sbar()
        {
            // get the viscosity FI
            Field_Interpolator* tFIViscosity =
                    mFIManager->get_field_interpolators_for_type( mDofViscosity );

            // get wall distance
            real tWallDistance = mPropWallDistance->val()( 0 );

            // threshold denominator
            real tDeno = clip_value( std::pow( mKappa * tWallDistance, 2.0 ), mEpsilon );

            // compute sbar
            mSBar = this->fv2() * tFIViscosity->val()( 0 ) / tDeno;

            // xxx to be removed
            clip_value( std::pow( mKappa * tWallDistance, 2.0 ), mEpsilon, this->fv2() * tFIViscosity->val()( 0 ) );
        }

        //--------------------------------------------------------------------------------------------------------------

        real
        CM_Spalart_Allmaras_Turbulence::sbar(
                enum CM_Function_Type aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Spalart_Allmaras_Turbulence::sbar - Only DEFAULT CM function type known in base class." );

            // if not evaluated
            if ( mSBarEval )
            {
                // evaluate
                this->eval_sbar();

                // set bool for evaluation
                mSBarEval = false;
            }
            // return
            return mSBar;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Spalart_Allmaras_Turbulence::eval_dsbardu(
                const moris::Cell< MSI::Dof_Type >& aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get derivative dof type FI
            Field_Interpolator* tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init dsbardu
            mdSBardu( tDofIndex ).set_size( 1, tFIDer->get_number_of_space_time_coefficients() );

            // get the viscosity FI
            Field_Interpolator* tFIViscosity =
                    mFIManager->get_field_interpolators_for_type( mDofViscosity );

            // get the wall distance value
            real tWallDistance = mPropWallDistance->val()( 0 );

            // threshold denominator
            real tDeno = clip_value( std::pow( mKappa * tWallDistance, 2.0 ), mEpsilon );

            // compute dsbardu
            mdSBardu( tDofIndex ) =
                    tFIViscosity->val() * this->dfv2du( aDofTypes ) / tDeno;

            // xxx to be removed
            clip_value( std::pow( mKappa * tWallDistance, 2.0 ), mEpsilon, tFIViscosity->val()( 0 ) );

            // if derivative dof type is viscosity dof type
            if ( aDofTypes( 0 ) == mDofViscosity )
            {
                // add contribution
                mdSBardu( tDofIndex ) += this->fv2() * tFIViscosity->N() / tDeno;
            }

            // if wall distance depends on derivative dof type
            if ( ( mPropWallDistance->check_dof_dependency( aDofTypes ) ) && ( tDeno > mEpsilon ) )
            {
                // threshold denominator
                real tKappa2WallDistance3 = clip_value( std::pow( mKappa, 2.0 ) * std::pow( tWallDistance, 3.0 ), mEpsilonDeriv );

                // add contribution to dsbardu
                mdSBardu( tDofIndex ) -=
                        2.0 * this->fv2() * tFIViscosity->val()( 0 ) * mPropWallDistance->dPropdDOF( aDofTypes ) / tKappa2WallDistance3;

                // xxx to be removed
                clip_value( std::pow( mKappa, 2.0 ) * std::pow( tWallDistance, 3.0 ), mEpsilonDeriv, 2.0 * this->fv2() * tFIViscosity->val()( 0 ) );
            }

            MORIS_ASSERT( isfinite( mdSBardu( tDofIndex ) ),
                    "CM_Spalart_Allmaras_Turbulence::eval_dsbardu - mdSBardu contains NAN or INF, exiting!" );
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Spalart_Allmaras_Turbulence::dsbardu(
                const moris::Cell< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type               aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Spalart_Allmaras_Turbulence::dsbardu - Only DEFAULT CM function type known in base class." );

            // if aDofType is not an active dof type for the CM
            MORIS_ERROR(
                    this->check_dof_dependency( aDofType ),
                    "CM_Spalart_Allmaras_Turbulence::dsbardu - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if ( mdSBarduEval( tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_dsbardu( aDofType );

                // set bool for evaluation
                mdSBarduEval( tDofIndex ) = false;
            }

            // return the derivative
            return mdSBardu( tDofIndex );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Spalart_Allmaras_Turbulence::eval_smod()
        {
            // threshold deno
            real tDeno = clip_value( ( mCv3 - 2.0 * mCv2 ) * this->s() - this->sbar(), mEpsilon );

            // compute smod
            mSMod = this->s() * ( std::pow( mCv2, 2 ) * this->s() + mCv3 * this->sbar() ) / tDeno;

            // xxx to be removed
            clip_value( ( mCv3 - 2.0 * mCv2 ) * this->s() - this->sbar(), mEpsilon, this->s() * ( std::pow( mCv2, 2 ) * this->s() + mCv3 * this->sbar() ) );
        }

        //--------------------------------------------------------------------------------------------------------------

        real
        CM_Spalart_Allmaras_Turbulence::smod(
                enum CM_Function_Type aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Spalart_Allmaras_Turbulence::smod - Only DEFAULT CM function type known in base class." );

            // if not evaluated
            if ( mSModEval )
            {
                // evaluate
                this->eval_smod();

                // set bool for evaluation
                mSModEval = false;
            }
            // return
            return mSMod;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Spalart_Allmaras_Turbulence::eval_dsmoddu(
                const moris::Cell< MSI::Dof_Type >& aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get derivative dof type FI
            Field_Interpolator* tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init dSModdu
            mdSModdu( tDofIndex ).set_size( 1, tFIDer->get_number_of_space_time_coefficients() );

            // compute smod num
            real tSModNum = this->s() * ( std::pow( mCv2, 2 ) * this->s() + mCv3 * this->sbar() );

            // compute smod deno and threshold deno
            real tSModDeno = clip_value( ( mCv3 - 2.0 * mCv2 ) * this->s() - this->sbar(), mEpsilon );

            if ( std::abs( tSModDeno ) > mEpsilon )
            {
                // threshold denominator
                real tSModDeno2 = clip_value( std::pow( tSModDeno, 2.0 ), mEpsilonDeriv );

                // compute dsmoddu
                mdSModdu( tDofIndex ) =
                        ( ( this->dsdu( aDofTypes ) * ( std::pow( mCv2, 2 ) * this->s() + mCv3 * this->sbar() ) + this->s() * ( std::pow( mCv2, 2 ) * this->dsdu( aDofTypes ) + mCv3 * this->dsbardu( aDofTypes ) ) ) * tSModDeno - tSModNum * ( ( mCv3 - 2.0 * mCv2 ) * this->dsdu( aDofTypes ) - this->dsbardu( aDofTypes ) ) ) / tSModDeno2;
            }
            else
            {
                // compute dsmoddu
                mdSModdu( tDofIndex ) =
                        ( this->dsdu( aDofTypes ) * ( std::pow( mCv2, 2 ) * this->s() + mCv3 * this->sbar() ) + this->s() * ( std::pow( mCv2, 2 ) * this->dsdu( aDofTypes ) + mCv3 * this->dsbardu( aDofTypes ) ) ) / tSModDeno;
            }

            MORIS_ASSERT( isfinite( mdSModdu( tDofIndex ) ),
                    "CM_Spalart_Allmaras_Turbulence::eval_dsmoddu - mdSModdu contains NAN or INF, exiting!" );
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Spalart_Allmaras_Turbulence::dsmoddu(
                const moris::Cell< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type               aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Spalart_Allmaras_Turbulence::dsmoddu - Only DEFAULT CM function type known in base class." );

            // if aDofType is not an active dof type for the CM
            MORIS_ERROR(
                    this->check_dof_dependency( aDofType ),
                    "CM_Spalart_Allmaras_Turbulence::dsmoddu - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if ( mdSModduEval( tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_dsmoddu( aDofType );

                // set bool for evaluation
                mdSModduEval( tDofIndex ) = false;
            }

            // return the derivative
            return mdSModdu( tDofIndex );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Spalart_Allmaras_Turbulence::eval_stilde()
        {
            // init mSTilde
            mSTilde = this->s();

            // compute STilde
            if ( this->sbar() >= -mCv2 * this->s() )
            {
                mSTilde += this->sbar();
            }
            else
            {
                // compute sMod
                mSTilde += this->smod();
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        real
        CM_Spalart_Allmaras_Turbulence::stilde(
                enum CM_Function_Type aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Spalart_Allmaras_Turbulence::stilde - Only DEFAULT CM function type known in base class." );

            // if not evaluated
            if ( mSTildeEval )
            {
                // evaluate
                this->eval_stilde();

                // set bool for evaluation
                mSTildeEval = false;
            }
            // return
            return mSTilde;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Spalart_Allmaras_Turbulence::eval_dstildedu(
                const moris::Cell< MSI::Dof_Type >& aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get derivative dof type FI
            Field_Interpolator* tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set size dSTildedu
            mdSTildedu( tDofIndex ).set_size( 1, tFIDer->get_number_of_space_time_coefficients() );

            // init dstildedu
            mdSTildedu( tDofIndex ) = this->dsdu( aDofTypes );

            // compute dstildedu
            if ( this->sbar() >= -mCv2 * this->s() )
            {
                // add dSbardu
                mdSTildedu( tDofIndex ) += this->dsbardu( aDofTypes );
            }
            else
            {
                // add dsmoddu
                mdSTildedu( tDofIndex ) += this->dsmoddu( aDofTypes );
            }

            MORIS_ASSERT( isfinite( mdSTildedu( tDofIndex ) ),
                    "CM_Spalart_Allmaras_Turbulence::eval_dstildedu - mdSTildedu contains NAN or INF, exiting!" );
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Spalart_Allmaras_Turbulence::dstildedu(
                const moris::Cell< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type               aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Spalart_Allmaras_Turbulence::dstildedu - Only DEFAULT CM function type known in base class." );

            // if aDofType is not an active dof type for the CM
            MORIS_ERROR(
                    this->check_dof_dependency( aDofType ),
                    "CM_Spalart_Allmaras_Turbulence::dstildedu - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if ( mdSTildeduEval( tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_dstildedu( aDofType );

                // set bool for evaluation
                mdSTildeduEval( tDofIndex ) = false;
            }

            // return the derivative
            return mdSTildedu( tDofIndex );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Spalart_Allmaras_Turbulence::eval_r()
        {
            // get the residual dof FI (here viscosity)
            Field_Interpolator* tFIViscosity =
                    mFIManager->get_field_interpolators_for_type( mDofViscosity );

            // get the wall distance value
            real tWallDistance = mPropWallDistance->val()( 0 );

            // threshold deno
            real tRDeno = clip_value( this->stilde() * std::pow( mKappa * tWallDistance, 2.0 ), mEpsilon );

            // compute viscosity / ( stilde * kappa² * d² )
            mR = tFIViscosity->val()( 0 ) / tRDeno;

            // select minimum
            // mR = std::min( mR, mRLim );

            // xxx to be removed and line above commented out
            if ( mR < mRLim )
            {

                clip_value( this->stilde() * std::pow( mKappa * tWallDistance, 2.0 ), mEpsilon, tFIViscosity->val()( 0 ) );
            }
            else
            {
                mR = mRLim;
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        real
        CM_Spalart_Allmaras_Turbulence::r(
                enum CM_Function_Type aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Spalart_Allmaras_Turbulence::r - Only DEFAULT CM function type known in base class." );

            // if not evaluated
            if ( mREval )
            {
                // evaluate
                this->eval_r();

                // set bool for evaluation
                mREval = false;
            }
            // return
            return mR;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Spalart_Allmaras_Turbulence::eval_drdu(
                const moris::Cell< MSI::Dof_Type >& aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get derivative dof type FI
            Field_Interpolator* tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set size dRdu
            mdRdu( tDofIndex ).set_size( 1, tFIDer->get_number_of_space_time_coefficients() );

            // if r < rlim
            if ( this->r() < mRLim )
            {
                // get the viscosity dof FI
                Field_Interpolator* tFIViscosity =
                        mFIManager->get_field_interpolators_for_type( mDofViscosity );

                // get the wall distance value
                real tWallDistance = mPropWallDistance->val()( 0 );

                // threshold deno
                real tRDeno = clip_value( this->stilde() * std::pow( mKappa * tWallDistance, 2.0 ), mEpsilon );

                // if dof type is viscosity
                if ( aDofTypes( 0 ) == mDofViscosity )
                {
                    // add contribution from viscosity
                    mdRdu( tDofIndex ) = tFIDer->N() / tRDeno;
                }
                else
                {
                    mdRdu( tDofIndex ).fill( 0.0 );
                }

                // check deno
                if ( std::abs( tRDeno ) > mEpsilon )
                {
                    // threshold denominator
                    real tdRduDeno1 = clip_value( std::pow( this->stilde() * mKappa * tWallDistance, 2.0 ), mEpsilonDeriv );

                    // add contribution from dStildedu
                    mdRdu( tDofIndex ) -= tFIViscosity->val() * this->dstildedu( aDofTypes ) / tdRduDeno1;

                    // xxx to be removed
                    clip_value( std::pow( this->stilde() * mKappa * tWallDistance, 2.0 ), mEpsilonDeriv, tFIViscosity->val()( 0 ) );

                    // if wall distance depends on derivative dof type
                    if ( mPropWallDistance->check_dof_dependency( aDofTypes ) )
                    {
                        // threshold denominator
                        real tdRduDeno2 = clip_value(
                                this->stilde() * std::pow( mKappa, 2.0 ) * std::pow( tWallDistance, 3.0 ),
                                mEpsilonDeriv );

                        // add contribution from wall distance
                        mdRdu( tDofIndex ) -=
                                2.0 * tFIViscosity->val()( 0 ) * mPropWallDistance->dPropdDOF( aDofTypes ) / tdRduDeno2;

                        // xxx to be removed
                        clip_value(
                                this->stilde() * std::pow( mKappa, 2.0 ) * std::pow( tWallDistance, 3.0 ),
                                mEpsilonDeriv,
                                2.0 * tFIViscosity->val()( 0 ) );
                    }
                }
            }
            else
            {
                mdRdu( tDofIndex ).fill( 0.0 );
            }

            MORIS_ASSERT( isfinite( mdRdu( tDofIndex ) ),
                    "CM_Spalart_Allmaras_Turbulence::eval_drdu - mdRdu contains NAN or INF, exiting!" );
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Spalart_Allmaras_Turbulence::drdu(
                const moris::Cell< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type               aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Spalart_Allmaras_Turbulence::drdu - Only DEFAULT CM function type known in base class." );

            // if aDofType is not an active dof type for the CM
            MORIS_ERROR(
                    this->check_dof_dependency( aDofType ),
                    "CM_Spalart_Allmaras_Turbulence::drdu - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if ( mdRduEval( tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_drdu( aDofType );

                // set bool for evaluation
                mdRduEval( tDofIndex ) = false;
            }

            // return the derivative
            return mdRdu( tDofIndex );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Spalart_Allmaras_Turbulence::eval_g()
        {
            // compute g
            mG = this->r() + mCw2 * ( std::pow( this->r(), 6.0 ) - this->r() );
        }

        //--------------------------------------------------------------------------------------------------------------

        real
        CM_Spalart_Allmaras_Turbulence::g(
                enum CM_Function_Type aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Spalart_Allmaras_Turbulence::g - Only DEFAULT CM function type known in base class." );

            // if not evaluated
            if ( mGEval )
            {
                // evaluate
                this->eval_g();

                // set bool for evaluation
                mGEval = false;
            }
            // return
            return mG;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Spalart_Allmaras_Turbulence::eval_dgdu(
                const moris::Cell< MSI::Dof_Type >& aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get derivative dof type FI
            Field_Interpolator* tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set size dGdu
            mdGdu( tDofIndex ).set_size( 1, tFIDer->get_number_of_space_time_coefficients() );

            // compute dGdu
            mdGdu( tDofIndex ) = ( 1.0 + mCw2 * ( 6.0 * std::pow( this->r(), 5.0 ) - 1.0 ) ) * this->drdu( aDofTypes );

            MORIS_ASSERT( isfinite( mdGdu( tDofIndex ) ),
                    "CM_Spalart_Allmaras_Turbulence::eval_dgdu - mdGdu contains NAN or INF, exiting!" );
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Spalart_Allmaras_Turbulence::dgdu(
                const moris::Cell< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type               aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Spalart_Allmaras_Turbulence::dgdu - Only DEFAULT CM function type known in base class." );

            // if aDofType is not an active dof type for the CM
            MORIS_ERROR(
                    this->check_dof_dependency( aDofType ),
                    "CM_Spalart_Allmaras_Turbulence::dgdu - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if ( mdGduEval( tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_dgdu( aDofType );

                // set bool for evaluation
                mdGduEval( tDofIndex ) = false;
            }

            // return the derivative
            return mdGdu( tDofIndex );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Spalart_Allmaras_Turbulence::eval_fw()
        {
            // threshold deno
            real tFwDeno = clip_value( std::pow( this->g(), 6.0 ) + std::pow( mCw3, 6.0 ), mEpsilon );

            // compute fw
            mFw = ( 1.0 + std::pow( mCw3, 6.0 ) ) / tFwDeno;

            // xxx to be removed
            clip_value( std::pow( this->g(), 6.0 ) + std::pow( mCw3, 6.0 ), mEpsilon, ( 1.0 + std::pow( mCw3, 6.0 ) ) );

            mFw = this->g() * std::pow( mFw, 1.0 / 6.0 );
        }

        //--------------------------------------------------------------------------------------------------------------

        real
        CM_Spalart_Allmaras_Turbulence::fw(
                enum CM_Function_Type aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Spalart_Allmaras_Turbulence::fw - Only DEFAULT CM function type known in base class." );

            // if not evaluated
            if ( mFwEval )
            {
                // evaluate
                this->eval_fw();

                // set bool for evaluation
                mFwEval = false;
            }
            // return
            return mFw;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Spalart_Allmaras_Turbulence::eval_dfwdu(
                const moris::Cell< MSI::Dof_Type >& aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get derivative dof type FI
            Field_Interpolator* tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set size dFwdu
            mdFwdu( tDofIndex ).set_size( 1, tFIDer->get_number_of_space_time_coefficients() );

            // threshold deno
            real tFwDeno = clip_value( std::pow( this->g(), 6.0 ) + std::pow( mCw3, 6.0 ), mEpsilon );

            if ( tFwDeno > mEpsilon )
            {
                // threshold denominator
                real tdFwduDeno = clip_value( this->g() * tFwDeno, mEpsilonDeriv );

                // compute dfwdu
                mdFwdu( tDofIndex ) = ( this->fw() * std::pow( mCw3, 6.0 ) * this->dgdu( aDofTypes ) ) / tdFwduDeno;

                // xxx to be removed
                clip_value( this->g() * tFwDeno, mEpsilonDeriv, this->fw() * std::pow( mCw3, 6.0 ) );

                MORIS_ASSERT( isfinite( mdFwdu( tDofIndex ) ),
                        "CM_Spalart_Allmaras_Turbulence::eval_dfwdu - mdFwdu( tDofIndex ) has NAN or INF (1)" );
            }
            else
            {
                // compute dfwdu
                mdFwdu( tDofIndex ) = this->dgdu( aDofTypes ) * std::pow( ( 1.0 + std::pow( mCw3, 6.0 ) ) / tFwDeno, 1.0 / 6.0 );

                MORIS_ASSERT( isfinite( mdFwdu( tDofIndex ) ),
                        "CM_Spalart_Allmaras_Turbulence::eval_dfwdu - mdFwdu( tDofIndex ) has NAN or INF (2)" );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Spalart_Allmaras_Turbulence::dfwdu(
                const moris::Cell< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type               aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Spalart_Allmaras_Turbulence::dfwdu - Only DEFAULT CM function type known in base class." );

            // if aDofType is not an active dof type for the CM
            MORIS_ERROR(
                    this->check_dof_dependency( aDofType ),
                    "CM_Spalart_Allmaras_Turbulence::dfwdu - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if ( mdFwduEval( tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_dfwdu( aDofType );

                // set bool for evaluation
                mdFwduEval( tDofIndex ) = false;
            }

            // return the derivative
            return mdFwdu( tDofIndex );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Spalart_Allmaras_Turbulence::eval_fn()
        {
            // threshold deno
            real tFnDeno = clip_value( mCn1 - std::pow( this->chi(), 3.0 ), mEpsilon );

            // xxx to be removed
            clip_value( mCn1 - std::pow( this->chi(), 3.0 ), mEpsilon, ( mCn1 + std::pow( this->chi(), 3.0 ) ) );

            // compute fn
            mFn = ( mCn1 + std::pow( this->chi(), 3.0 ) ) / tFnDeno;
        }

        //--------------------------------------------------------------------------------------------------------------

        real
        CM_Spalart_Allmaras_Turbulence::fn(
                enum CM_Function_Type aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Spalart_Allmaras_Turbulence::fn - Only DEFAULT CM function type known in base class." );

            // if not evaluated
            if ( mFnEval )
            {
                // evaluate
                this->eval_fn();

                // set bool for evaluation
                mFnEval = false;
            }
            // return
            return mFn;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Spalart_Allmaras_Turbulence::eval_dfndu(
                const moris::Cell< MSI::Dof_Type >& aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get derivative dof type FI
            Field_Interpolator* tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set size dFndu
            mdFndu( tDofIndex ).set_size( 1, tFIDer->get_number_of_space_time_coefficients() );

            // threshold deno
            real tFnDeno = clip_value( mCn1 - std::pow( this->chi(), 3.0 ), mEpsilon );

            if ( std::abs( tFnDeno ) > mEpsilon )
            {
                // threshold denominator
                real tdFnduDeno = clip_value( std::pow( tFnDeno, 2.0 ), mEpsilonDeriv );

                // compute dfndu
                mdFndu( tDofIndex ) =
                        6.0 * mCn1 * std::pow( this->chi(), 2 ) * this->dchidu( aDofTypes ) / tdFnduDeno;

                // xxx to be removed
                clip_value( std::pow( tFnDeno, 2.0 ), mEpsilonDeriv, 6.0 * mCn1 * std::pow( this->chi(), 2 ) );
            }
            else
            {
                // compute dfndu
                mdFndu( tDofIndex ) =
                        -3.0 * mCn1 * std::pow( this->chi(), 2.0 ) * this->dchidu( aDofTypes ) / tFnDeno;

                // xxx to be removed
                clip_value( mCn1 - std::pow( this->chi(), 3.0 ), mEpsilon, -3.0 * mCn1 * std::pow( this->chi(), 2.0 ) );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Spalart_Allmaras_Turbulence::dfndu(
                const moris::Cell< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type               aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Spalart_Allmaras_Turbulence::dfndu - Only DEFAULT CM function type known in base class." );

            // if aDofType is not an active dof type for the CM
            MORIS_ERROR(
                    this->check_dof_dependency( aDofType ),
                    "CM_Spalart_Allmaras_Turbulence::dfndu - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if ( mdFnduEval( tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_dfndu( aDofType );

                // set bool for evaluation
                mdFnduEval( tDofIndex ) = false;
            }

            // return the derivative
            return mdFndu( tDofIndex );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Spalart_Allmaras_Turbulence::eval_dfndx( uint aOrder )
        {
            // FIXME work only for 1st order
            MORIS_ERROR( aOrder == 1,
                    "CM_Spalart_Allmaras_Turbulence::eval_dfndx - Works only for 1st order derivative for now." );

            // set matrix size
            mdFndx( aOrder - 1 ).set_size( mSpaceDim, 1 );

            // threshold deno
            real tFnDeno = clip_value( mCn1 - std::pow( this->chi(), 3.0 ), mEpsilon );

            if ( std::abs( tFnDeno ) > mEpsilon )
            {
                // threshold denominator
                real tdFnduDeno = clip_value( std::pow( tFnDeno, 2.0 ), mEpsilonDeriv );

                // compute dfndx
                mdFndx( aOrder - 1 ) =
                        6.0 * mCn1 * std::pow( this->chi(), 2 ) * this->dchidx( 1 ) / tdFnduDeno;

                // xxx to be removed
                clip_value( std::pow( tFnDeno, 2.0 ), mEpsilonDeriv, 6.0 * mCn1 * std::pow( this->chi(), 2 ) );
            }
            else
            {
                // compute dfndx
                mdFndx( aOrder - 1 ) =
                        -3.0 * mCn1 * std::pow( this->chi(), 2.0 ) * this->dchidx( 1 ) / tFnDeno;

                // xxx to be removed
                clip_value( mCn1 - std::pow( this->chi(), 3.0 ), mEpsilon, -3.0 * mCn1 * std::pow( this->chi(), 2.0 ) );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Spalart_Allmaras_Turbulence::dfndx(
                uint                  aOrder,
                enum CM_Function_Type aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Spalart_Allmaras_Turbulence::dfndx - Only DEFAULT CM function type known in base class." );

            MORIS_ERROR( aOrder == 1,
                    "CM_Spalart_Allmaras_Turbulence::dfndx - Works only for 1st order derivative for now." );

            // if the derivative has not been evaluated yet
            if ( mdFndxEval( aOrder - 1 ) )
            {
                // evaluate the derivative
                this->eval_dfndx( aOrder );

                // set bool for evaluation
                mdFndxEval( aOrder - 1 ) = false;
            }

            // return the derivative
            return mdFndx( aOrder - 1 );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Spalart_Allmaras_Turbulence::eval_dfndxdu(
                const moris::Cell< MSI::Dof_Type >& aDofTypes,
                uint                                aOrder )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get derivative dof type FI
            Field_Interpolator* tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set matrix size
            mdFndxdu( aOrder - 1 )( tDofIndex ).set_size( mSpaceDim, tFIDer->get_number_of_space_time_coefficients() );

            // threshold deno
            real tFnDeno = clip_value( mCn1 - std::pow( this->chi(), 3.0 ), mEpsilon );

            if ( std::abs( tFnDeno ) > mEpsilon )
            {
                // threshold denominator
                real tdFndxduDeno = clip_value( std::pow( tFnDeno, 3.0 ), mEpsilonDeriv );

                // compute dfndxdu
                mdFndxdu( aOrder - 1 )( tDofIndex ) =
                        6.0 * mCn1 * ( 2.0 * this->chi() * ( mCn1 + 2.0 * std::pow( this->chi(), 3.0 ) ) * this->dchidx( 1 ) * this->dchidu( aDofTypes ) + tFnDeno * std::pow( this->chi(), 2.0 ) * this->dchidxdu( aDofTypes, 1 ) ) / tdFndxduDeno;
            }
            else
            {
                mdFndxdu( aOrder - 1 )( tDofIndex ) = -3.0 * mCn1 * ( 2.0 * this->chi() * this->dchidx( 1 ) * this->dchidu( aDofTypes ) + std::pow( this->chi(), 2.0 ) * this->dchidxdu( aDofTypes, 1 ) ) / tFnDeno;
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Spalart_Allmaras_Turbulence::dfndxdu(
                const moris::Cell< MSI::Dof_Type >& aDofType,
                uint                                aOrder,
                enum CM_Function_Type               aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Spalart_Allmaras_Turbulence::dfndxdu - Only DEFAULT CM function type known in base class." );

            MORIS_ERROR( aOrder == 1,
                    "CM_Spalart_Allmaras_Turbulence::dfndxdu - Works only for 1st order derivative for now." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if ( mdFndxduEval( aOrder - 1, tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_dfndxdu( aDofType, aOrder );

                // set bool for evaluation
                mdFndxduEval( aOrder - 1, tDofIndex ) = false;
            }

            // return the derivative
            return mdFndxdu( aOrder - 1 )( tDofIndex );
        }

        //--------------------------------------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
