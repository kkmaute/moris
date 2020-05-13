#include "cl_FEM_CM_Diffusion_Linear_Isotropic_Phase_Change.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include <iostream>

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"
#include "fn_sum.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------
        CM_Diffusion_Linear_Isotropic_Phase_Change::CM_Diffusion_Linear_Isotropic_Phase_Change()
        {
            // set the property pointer cell size
            mProperties.resize( static_cast< uint >( Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Conductivity" ]          = Property_Type::CONDUCTIVITY;
            mPropertyMap[ "Density" ]               = Property_Type::DENSITY;
            mPropertyMap[ "Heat_Capacity" ]         = Property_Type::HEAT_CAPACITY;
            mPropertyMap[ "Latent_Heat" ]           = Property_Type::LATENT_HEAT;
            mPropertyMap[ "Lower_PC_Temp" ]         = Property_Type::LOWER_PC_TEMP;
            mPropertyMap[ "Upper_PC_Temp" ]         = Property_Type::UPPER_PC_TEMP;
            mPropertyMap[ "Phase_State_Function" ]  = Property_Type::PHASE_STATE_FUNCTION;
            mPropertyMap[ "Phase_Change_Const" ]    = Property_Type::PHASE_CHANGE_CONST;

            // populate dof map
            mDofMap[ "Temp" ] = MSI::Dof_Type::TEMP;
        }

        //------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic_Phase_Change::eval_gradHdot()
        {
            moris::real tDensity = mProperties( static_cast< uint >( Property_Type::DENSITY ) )->val()( 0 );
            moris::real tHeatCap = mProperties( static_cast< uint >( Property_Type::HEAT_CAPACITY ) )->val()( 0 );
            moris::real tLatHeat = mProperties( static_cast< uint >( Property_Type::LATENT_HEAT ) )->val()( 0 );

            // compute derivative of Phase State Function
            real tdfdT = this->eval_dFdTemp();

            // compute gradient of
            mGradHdot = tDensity * ( tHeatCap + tLatHeat * tdfdT ) *
                    mFIManager->get_field_interpolators_for_type( mDofMap[ "Temp" ] )->gradxt();
        }

        //------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic_Phase_Change::eval_graddivflux()
        {
            moris::real tK = mProperties( static_cast< uint >( Property_Type::CONDUCTIVITY ) )->val()( 0 );

            // matrix for purely isotropic case
            Matrix< DDRMat > tKijIsotropic;
            switch( mSpaceDim )
            {
                case 2:
                {
                    tKijIsotropic = {{tK,  0,  0, tK},
                            { 0, tK, tK,  0}};
                    break;
                }
                case 3:
                {
                    tKijIsotropic = {{tK, 0, 0, 0, 0,tK, 0,tK, 0, 0},
                            { 0,tK, 0,tK, 0, 0, 0, 0,tK, 0},
                            { 0, 0,tK, 0,tK, 0,tK, 0, 0, 0}};
                    break;
                }
                default:
                    MORIS_ASSERT(false, "CM_Diffusion_Linear_Isotropic_Phase_Change::eval_graddivflux: Number of spatial dimensions must be 2 or 3");
                    break;
            }

            mGradDivFlux = tKijIsotropic * mFIManager->get_field_interpolators_for_type( mDofMap[ "Temp" ] )->gradx(3);
        }

        //------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic_Phase_Change::eval_flux()
        {
            // get the conductivity property
            std::shared_ptr< Property > tPropConductivity =
                    mProperties( static_cast< uint >( Property_Type::CONDUCTIVITY ) );

            // compute flux
            mFlux = tPropConductivity->val()( 0 ) *
                    mFIManager->get_field_interpolators_for_type( mDofMap[ "Temp" ] )->gradx( 1 );
        }

        //------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic_Phase_Change::eval_divflux()
        {
            // get the conductivity property
            std::shared_ptr< Property > tPropConductivity
            = mProperties( static_cast< uint >( Property_Type::CONDUCTIVITY ) );

            // compute the divergence of the flux
            mDivFlux = tPropConductivity->val() * this->divstrain();
        }

        //------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic_Phase_Change::eval_traction( const Matrix< DDRMat > & aNormal )
        {
            // compute traction
            mTraction = trans( this->flux() ) * aNormal;
        }

        //------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic_Phase_Change::eval_testTraction( const Matrix< DDRMat > & aNormal,
                const moris::Cell< MSI::Dof_Type > & aTestDofTypes )
        {
            // get test dof type index
            uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

            // get the conductivity property
            std::shared_ptr< Property > tPropConductivity =
                    mProperties( static_cast< uint >( Property_Type::CONDUCTIVITY ) );

            // compute test traction
            mTestTraction( tTestDofIndex ) =
                    trans( mFIManager->get_field_interpolators_for_type( mDofMap[ "Temp" ] )->dnNdxn( 1 ) ) *
                    tPropConductivity->val()( 0 ) * aNormal;
        }

        //------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic_Phase_Change::eval_strain()
        {
            // compute strain
            mStrain = mFIManager->get_field_interpolators_for_type( mDofMap[ "Temp" ] )->gradx( 1 );
        }

        //------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic_Phase_Change::eval_divstrain()
        {
            // get the temperature gradient
            Matrix< DDRMat > tTempGrad
            = mFIManager->get_field_interpolators_for_type( mDofMap[ "Temp" ] )->gradx( 2 );

            // evaluate the divergence of the strain
            mDivStrain = sum( tTempGrad( { 0, mSpaceDim - 1 }, { 0, 0 } ) );
        }

        //------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic_Phase_Change::eval_testStrain()
        {
            // compute test strain
            mTestStrain = mFIManager->get_field_interpolators_for_type( mDofMap[ "Temp" ] )->dnNdxn( 1 );
        }

        //------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic_Phase_Change::eval_const()
        {
            // build an identity matrix
            Matrix< DDRMat > I;
            eye( mSpaceDim, mSpaceDim, I );

            // compute conductivity matrix
            mConst = mProperties( static_cast< uint >( Property_Type::CONDUCTIVITY ) )->val()( 0 ) * I;
        }

        //------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic_Phase_Change::eval_dFluxdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the corresponding FI
            Field_Interpolator * tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            //init the matrix
            mdFluxdDof( tDofIndex ).set_size( mSpaceDim, tFI->get_number_of_space_time_coefficients(), 0.0 );

            // get the conductivity property
            std::shared_ptr< Property > tPropConductivity
            = mProperties( static_cast< uint >( Property_Type::CONDUCTIVITY ) );

            // if direct dependency on the dof type
            if( aDofTypes( 0 ) == mDofMap[ "Temp" ] )
            {
                // compute derivative with direct dependency
                mdFluxdDof( tDofIndex ).matrix_data() +=
                        tPropConductivity->val()( 0 ) * mFIManager->get_field_interpolators_for_type( mDofMap[ "Temp" ] )->dnNdxn( 1 );
            }

            //            // if indirect dependency on the dof type
            //            if ( mProperties( static_cast< uint >( Property_Type::CONDUCTIVITY ) )->check_dof_dependency( aDofTypes ) )
            //            {
            //                // compute derivative with indirect dependency through properties
            //                mdFluxdDof( tDofIndex ).matrix_data()
            //                += mFIManager->get_field_interpolators_for_type( mDofMap[ "Temp" ] )->gradx( 1 )
            //                 * tPropConductivity->dPropdDOF( aDofTypes );
            //            }
        }

        //------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic_Phase_Change::eval_dHdotdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get properties
            moris::real tDensity = mProperties( static_cast< uint >( Property_Type::DENSITY ) )->val()( 0 );
            moris::real tHeatCap = mProperties( static_cast< uint >( Property_Type::HEAT_CAPACITY ) )->val()( 0 );
            moris::real tLatHeat = mProperties( static_cast< uint >( Property_Type::LATENT_HEAT ) )->val()( 0 );

            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the corresponding FI
            Field_Interpolator * tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init the matrix
            mHdotDof( tDofIndex ).set_size( 1, tFI->get_number_of_space_time_coefficients(), 0.0 );

            // compute derivative of Phase State Function
            real tdfdT = this->eval_dFdTemp();

            // if direct dependency on the dof type
            if( aDofTypes( 0 ) == mDofMap[ "Temp" ] )
            {
                // compute derivative with direct dependency
                mHdotDof( tDofIndex ).matrix_data() +=
                        tDensity * ( tHeatCap + tLatHeat * tdfdT ) *
                        mFIManager->get_field_interpolators_for_type( mDofMap[ "Temp" ] )->dnNdtn(1);
            }
        }

        //--------------------------------------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic_Phase_Change::eval_dGradHdotdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get properties
            moris::real tDensity = mProperties( static_cast< uint >( Property_Type::DENSITY ) )->val()( 0 );
            moris::real tHeatCap = mProperties( static_cast< uint >( Property_Type::HEAT_CAPACITY ) )->val()( 0 );
            moris::real tLatHeat = mProperties( static_cast< uint >( Property_Type::LATENT_HEAT ) )->val()( 0 );

            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the corresponding FI
            Field_Interpolator * tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init the matrix
            mGradHdotDof( tDofIndex ).set_size( mSpaceDim, tFI->get_number_of_space_time_coefficients(), 0.0 );

            // compute derivative of Phase State Function
            real tdfdT = this->eval_dFdTemp();

            // if direct dependency on the dof type
            if( aDofTypes( 0 ) == mDofMap[ "Temp" ] )
            {
                // compute derivative with direct dependency
                mGradHdotDof( tDofIndex ).matrix_data() +=
                        tDensity * ( tHeatCap + tLatHeat * tdfdT ) *
                        mFIManager->get_field_interpolators_for_type( mDofMap[ "Temp" ] )->d2Ndxt();
            }
        }

        //--------------------------------------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic_Phase_Change::eval_dGradDivFluxdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // gets added later
            mGradDivFlux = {{0}};

            moris::real tK = mProperties( static_cast< uint >( Property_Type::CONDUCTIVITY ) )->val()( 0 );

            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the corresponding FI
            Field_Interpolator * tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init the matrix
            mGradDivFluxDof( tDofIndex ).set_size( mSpaceDim, tFI->get_number_of_space_time_coefficients(), 0.0 );

            // matrix for purely isotropic case
            Matrix< DDRMat > tKijIsotropic;
            switch( mSpaceDim )
            {
                case 2:
                {
                    tKijIsotropic = {{tK,  0,  0, tK},
                            { 0, tK, tK,  0}};
                    break;
                }
                case 3:
                {
                    tKijIsotropic = {{tK, 0, 0, 0, 0,tK, 0,tK, 0, 0},
                            { 0,tK, 0,tK, 0, 0, 0, 0,tK, 0},
                            { 0, 0,tK, 0,tK, 0,tK, 0, 0, 0}};
                    break;
                }
                default:
                    MORIS_ASSERT(false, "CM_Diffusion_Linear_Isotropic_Phase_Change::eval_dGradDivFluxdDOF: Number of spatial dimensions must be 2 or 3");
                    break;
            }

            // if direct dependency on the dof type
            if( aDofTypes( 0 ) == mDofMap[ "Temp" ] )
            {
                // compute derivative with direct dependency
                mGradDivFluxDof( tDofIndex ).matrix_data() +=
                        tKijIsotropic * mFIManager->get_field_interpolators_for_type( mDofMap[ "Temp" ] )->dnNdxn( 3 ).matrix_data();
            }
        }


        //--------------------------------------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic_Phase_Change::eval_ddivfluxdu( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the corresponding FI
            Field_Interpolator * tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set size for ddivflux/du
            mddivfluxdu( tDofIndex ).set_size( 1, tFI->get_number_of_space_time_coefficients(), 0.0 );

            // get the conductivity property
            std::shared_ptr< Property > tPropConductivity
            = mProperties( static_cast< uint >( Property_Type::CONDUCTIVITY ) );

            // if temperature dof
            if( aDofTypes( 0 ) == mDofMap[ "Temp" ] )
            {
                // fill ddivstrain/dv
                mddivfluxdu( tDofIndex ).matrix_data() += tPropConductivity->val()( 0 ) * this->ddivstraindu( aDofTypes );
            }

            if( tPropConductivity->check_dof_dependency( aDofTypes ) )
            {
                // fill ddivstrain/du
                mddivfluxdu( tDofIndex ).matrix_data() += this->divstrain() * tPropConductivity->dPropdDOF( aDofTypes );
            }
        }

        //--------------------------------------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic_Phase_Change::eval_ddivstraindu( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the dof type FI
            Field_Interpolator * tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set size for ddivstrain/du
            mddivstraindu( tDofIndex ).set_size( 1, tFI->get_number_of_space_time_coefficients(), 0.0 );

            if( aDofTypes( 0 ) == mDofMap[ "Temp" ] )
            {
                // get the 2nd order derivative of the shape functions d2Ndx2
                Matrix< DDRMat > tTempd2Ndx2 = tFI->dnNdxn( 2 );

                // fill ddivstrain/du
                mddivstraindu( tDofIndex ) = tTempd2Ndx2.get_row( 0 ) + tTempd2Ndx2.get_row( 1 );

                if( tTempd2Ndx2.n_rows() == 6 )
                {
                    mddivstraindu( tDofIndex ).matrix_data() += tTempd2Ndx2.get_row( 2 );
                }
            }
        }

        //------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic_Phase_Change::eval_dTractiondDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes,
                const Matrix< DDRMat >             & aNormal )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // compute derivative
            mdTractiondDof( tDofIndex ) = trans( aNormal ) * this->dFluxdDOF( aDofTypes );
        }

        //------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic_Phase_Change::eval_dTestTractiondDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes,
                const Matrix< DDRMat >             & aNormal,
                const moris::Cell< MSI::Dof_Type > & aTestDofTypes )
        {
            // get test dof type index
            uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // init the dTestTractiondDof
            mdTestTractiondDof( tTestDofIndex )( tDofIndex ).set_size( mFIManager->get_field_interpolators_for_type( aTestDofTypes( 0 ) )->get_number_of_space_time_coefficients(),
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) )->get_number_of_space_time_coefficients(), 0.0 );

            // get the conductivity property
            std::shared_ptr< Property > tPropConductivity
            = mProperties( static_cast< uint >( Property_Type::CONDUCTIVITY ) );

            // if conductivity depends on dof type
            if( tPropConductivity->check_dof_dependency( aDofTypes ) )
            {
                // add contribution
                mdTestTractiondDof( tTestDofIndex )( tDofIndex ).matrix_data()
                                                                += trans( mFIManager->get_field_interpolators_for_type( mDofMap[ "Temp" ] )->dnNdxn( 1 ) )
                                                                * aNormal * tPropConductivity->dPropdDOF( aDofTypes );
            }
        }

        //------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic_Phase_Change::eval_dTestTractiondDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes,
                const Matrix< DDRMat >             & aNormal,
                const Matrix< DDRMat >             & aJump,
                const moris::Cell< MSI::Dof_Type > & aTestDofTypes )
        {
            // get test dof type index
            uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // init the dTestTractiondDof
            mdTestTractiondDof( tTestDofIndex )( tDofIndex ).set_size( mFIManager->get_field_interpolators_for_type( aTestDofTypes( 0 ) )->get_number_of_space_time_coefficients(),
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) )->get_number_of_space_time_coefficients(), 0.0 );

            // get the conductivity property
            std::shared_ptr< Property > tPropConductivity
            = mProperties( static_cast< uint >( Property_Type::CONDUCTIVITY ) );

            // if conductivity depends on dof type
            if( tPropConductivity->check_dof_dependency( aDofTypes ) )
            {
                // add contribution
                mdTestTractiondDof( tTestDofIndex )( tDofIndex ).matrix_data() +=
                        trans( mFIManager->get_field_interpolators_for_type( mDofMap[ "Temp" ] )->dnNdxn( 1 ) ) *
                        aNormal * tPropConductivity->dPropdDOF( aDofTypes );
            }
        }

        //------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic_Phase_Change::eval_dStraindDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // if direct dependency on the dof type
            if( aDofTypes( 0 ) == mDofMap[ "Temp" ] )
            {
                // compute derivative with direct dependency
                mdStraindDof( tDofIndex ) = mFIManager->get_field_interpolators_for_type( mDofMap[ "Temp" ] )->dnNdxn( 1 );
            }
            else
            {
                // reset the matrix
                mdStraindDof( tDofIndex ).set_size( mSpaceDim, mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) )->get_number_of_space_time_coefficients(), 0.0 );
            }
        }

        //------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic_Phase_Change::eval_dConstdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // reset the matrix
            mdConstdDof( tDofIndex ).set_size( 1, mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) )->get_number_of_space_time_coefficients(), 0.0 );

            // if indirect dependency on the dof type
            if ( mProperties( static_cast< uint >( Property_Type::CONDUCTIVITY ) )->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mdConstdDof( tDofIndex ) = mProperties( static_cast< uint >( Property_Type::CONDUCTIVITY ) )->dPropdDOF( aDofTypes );
            }
        }

        //------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic_Phase_Change::eval_dFluxdDV( const moris::Cell< PDV_Type > & aDvTypes )
        {
            MORIS_ASSERT( false, " CM_Diffusion_Linear_Isotropic::eval_dFluxdDV - This function is not implemented.");
        }

        //------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic_Phase_Change::eval_dStraindDV( const moris::Cell< PDV_Type > & aDvTypes )
        {
            MORIS_ASSERT( false, " CM_Diffusion_Linear_Isotropic::eval_dStraindDV - This function is not implemented.");
        }

        //------------------------------------------------------------------------------
        moris::real CM_Diffusion_Linear_Isotropic_Phase_Change::eval_dFdTemp()
        {
            moris::real tTlower  = mProperties( static_cast< uint >( Property_Type::LOWER_PC_TEMP ) )->val()( 0 );
            moris::real tTupper  = mProperties( static_cast< uint >( Property_Type::UPPER_PC_TEMP ) )->val()( 0 );
            moris::uint tPCfunc  = mProperties( static_cast< uint >( Property_Type::PHASE_STATE_FUNCTION ) )->val()( 0 );

            // get temperature
            real tTemp = mFIManager->get_field_interpolators_for_type( mDofMap[ "Temp" ] )->val()( 0 );

            moris::real tdfdT;

            switch ( tPCfunc )
            {
                // linear phase change model
                case 1:
                {
                    if ( (tTemp > tTupper) || (tTemp < tTlower) )
                        tdfdT = 0;
                    else
                        tdfdT = 1 / (tTupper - tTlower);
                    break;
                }

                // cubic phase change model
                case 2:
                {
                    // cubic function: f(T)=( (tTemp - tTlower)^2 * (2*tTemp + tTlower - 3*tTupper))/(tTlower - tTupper)^3
                    if ( (tTemp > tTupper) || (tTemp < tTlower) )
                        tdfdT = 0;
                    else
                        tdfdT = 6.0*(tTemp - tTlower)*(tTemp - tTupper)/std::pow(tTlower - tTupper,3.0);
                    break;
                }

                // logistic function
                case 3:
                {
                    moris::real tPCconst = mProperties( static_cast< uint >( Property_Type::PHASE_CHANGE_CONST ) )->val()( 0 );

                    // logistic function parameter k
                    real tLogisticParam = ( 2 * std::log(1/tPCconst - 1) ) / ( tTupper - 3 * tTlower );

                    real tExp = tLogisticParam * ( tTemp - (tTupper + tTlower)/2 );
                    tdfdT = tLogisticParam  / ( std::exp(-tExp) + 2 + std::exp(tExp) );
                    break;
                }
                default:
                    MORIS_ERROR(false,"wrong option for phase change function.");
            }

            return tdfdT;
        }

        //------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic_Phase_Change::eval_Hdot()
        {
            // get properties
            moris::real tDensity = mProperties( static_cast< uint >( Property_Type::DENSITY ) )->val()( 0 );
            moris::real tHeatCap = mProperties( static_cast< uint >( Property_Type::HEAT_CAPACITY ) )->val()( 0 );
            moris::real tLatHeat = mProperties( static_cast< uint >( Property_Type::LATENT_HEAT ) )->val()( 0 );

            // compute derivative of Phase State Function
            real tdfdT = this->eval_dFdTemp();

            // compute derivative of enthalpy
            mHdot = tDensity * ( tHeatCap + tLatHeat * tdfdT )
                             * mFIManager->get_field_interpolators_for_type( mDofMap[ "Temp" ] )->gradt( 1 );
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
