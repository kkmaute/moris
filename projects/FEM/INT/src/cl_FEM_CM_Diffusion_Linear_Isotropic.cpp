
#include "cl_FEM_CM_Diffusion_Linear_Isotropic.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

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

        CM_Diffusion_Linear_Isotropic::CM_Diffusion_Linear_Isotropic()
        {
            // set the property pointer cell size
            mProperties.resize( static_cast< uint >( Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Conductivity" ] = Property_Type::CONDUCTIVITY;
            mPropertyMap[ "Density" ]      = Property_Type::DENSITY;
            mPropertyMap[ "HeatCapacity" ] = Property_Type::HEAT_CAPACITY;
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic::set_dof_type_list(
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

                // if temperature dof type string
                if( tDofString == "Temperature" )
                {
                    mTempDof = tDofType;
                }
                else
                {
                    // create error string
                    std::string tErrMsg =
                            std::string( "CM_Diffusion_Linear_Isotropic::set_dof_type_list - Unknown aDofString : " ) +
                            tDofString;

                    // error
                    MORIS_ERROR( false , tErrMsg.c_str() );
                }
            }
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic::set_property(
                std::shared_ptr< fem::Property > aProperty,
                std::string                      aPropertyString )
        {
            // check that aPropertyString makes sense
            if ( mPropertyMap.find( aPropertyString ) == mPropertyMap.end() )
            {
                std::string tErrMsg =
                        std::string( "CM_Diffusion_Linear_Isotropic::set_property - Unknown aPropertyString : ") +
                        aPropertyString;

                MORIS_ERROR( false , tErrMsg.c_str() );
            }

            // set the property in the property cell
            mProperties( static_cast< uint >( mPropertyMap[ aPropertyString ] ) ) = aProperty;
        }

        //------------------------------------------------------------------------------

        std::shared_ptr< Property > CM_Diffusion_Linear_Isotropic::get_property(
                std::string aPropertyString )
        {
            // check that aPropertyString makes sense
            if ( mPropertyMap.find( aPropertyString ) == mPropertyMap.end() )
            {
                std::string tErrMsg =
                        std::string( "CM_Diffusion_Linear_Isotropic::get_property - Unknown aPropertyString : ") +
                        aPropertyString;

                MORIS_ERROR( false , tErrMsg.c_str() );
            }

            // get the property in the property cell
            return  mProperties( static_cast< uint >( mPropertyMap[ aPropertyString ] ) );
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic::eval_flux()
        {
            // get the conductivity property
            std::shared_ptr< Property > tPropConductivity =
                    mProperties( static_cast< uint >( Property_Type::CONDUCTIVITY ) );

            // compute flux
            mFlux = tPropConductivity->val()( 0 ) *
                    mFIManager->get_field_interpolators_for_type( mTempDof )->gradx( 1 );
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic::eval_EnergyDot()
        {
            // get properties
            std::shared_ptr< Property > tPropDensity = mProperties( static_cast< uint >( Property_Type::DENSITY ) );
            std::shared_ptr< Property > tPropHeatCap = mProperties( static_cast< uint >( Property_Type::HEAT_CAPACITY ) );

            if (tPropDensity != nullptr && tPropHeatCap != nullptr)
            {
                // compute rate of enthalpy
                mEnergyDot = tPropDensity->val()( 0 ) *  tPropHeatCap->val()( 0 ) *
                        mFIManager->get_field_interpolators_for_type( mTempDof )->gradt( 1 );
            }
            else
            {
                // if no capacity or density is given, set EnergyDot to zero
                mEnergyDot = 0.0 * mFIManager->get_field_interpolators_for_type( mTempDof )->gradt( 1 );
            }
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic::eval_gradEnergy()
        {
            // get properties
            std::shared_ptr< Property > tPropDensity = mProperties( static_cast< uint >( Property_Type::DENSITY ) );
            std::shared_ptr< Property > tPropHeatCap = mProperties( static_cast< uint >( Property_Type::HEAT_CAPACITY ) );

            if (tPropDensity != nullptr && tPropHeatCap != nullptr)
            {
                // compute rate of gradient of enthalpy
                mGradEnergy = tPropDensity->val()( 0 ) *  tPropHeatCap->val()( 0 ) *
                        mFIManager->get_field_interpolators_for_type( mTempDof )->gradx( 1 );
            }
            else
            {
                // if no capacity or density is given, set gradEnergy to zero
                mGradEnergy = 0.0 * mFIManager->get_field_interpolators_for_type( mTempDof )->gradx( 1 );
            }
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic::eval_gradEnergyDot()
        {
            // get properties
            std::shared_ptr< Property > tPropDensity = mProperties( static_cast< uint >( Property_Type::DENSITY ) );
            std::shared_ptr< Property > tPropHeatCap = mProperties( static_cast< uint >( Property_Type::HEAT_CAPACITY ) );

            if (tPropDensity != nullptr && tPropHeatCap != nullptr)
            {
                // compute rate of gradient of enthalpy
                mGradEnergyDot = tPropDensity->val()( 0 ) *  tPropHeatCap->val()( 0 ) *
                        mFIManager->get_field_interpolators_for_type( mTempDof )->gradxt();
            }
            else
            {
                // if no capacity or density is given, set gradEnergyDot to zero
                mGradEnergyDot = 0.0 * mFIManager->get_field_interpolators_for_type( mTempDof )->gradxt();;
            }
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic::eval_divflux()
        {
            // get the conductivity property
            std::shared_ptr< Property > tPropConductivity =
                    mProperties( static_cast< uint >( Property_Type::CONDUCTIVITY ) );

            // compute the divergence of the flux
            mDivFlux = tPropConductivity->val() * this->divstrain();
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic::eval_graddivflux()
        {
            moris::real tK = mProperties( static_cast< uint >( Property_Type::CONDUCTIVITY ) )->val()( 0 );

            // get spatial interpolation order
            mtk::Interpolation_Order tInterpOrder =
                    mFIManager->get_field_interpolators_for_type( mTempDof )-> get_space_interpolation_order();

            // only compute if interpolation is not 1 or 0, in that case simply keep the vector of zeros
            if (tInterpOrder == mtk::Interpolation_Order::LINEAR || tInterpOrder == mtk::Interpolation_Order::CONSTANT)
            {
                // initialize results with zeros
                mGradDivFlux.set_size( mSpaceDim, 1, 0.0 );
            }
            else
            {
                // matrix for purely isotropic case
                // FIXME: this implementation is slow and needs to be improved
                Matrix< DDRMat > tKijIsotropic;
                switch( mSpaceDim )
                {
                    case 2:
                    {
                        tKijIsotropic = {
                                {tK,  0,  0, tK},
                                { 0, tK, tK,  0}};
                        break;
                    }
                    case 3:
                    {
                        tKijIsotropic = {
                                {tK, 0, 0, 0, 0,tK, 0,tK, 0, 0},
                                { 0,tK, 0,tK, 0, 0, 0, 0,tK, 0},
                                { 0, 0,tK, 0,tK, 0,tK, 0, 0, 0}};
                        break;
                    }
                    default:
                        MORIS_ASSERT(false, "CM_Diffusion_Linear_Isotropic::eval_graddivflux: Number of spatial dimensions must be 2 or 3");
                        break;
                }

                // compute grad div flux
                mGradDivFlux = tKijIsotropic * mFIManager->get_field_interpolators_for_type( mTempDof )->gradx(3);
            }

        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic::eval_traction( const Matrix< DDRMat > & aNormal )
        {
            // compute traction
            mTraction = trans( this->flux() ) * aNormal;
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic::eval_testTraction(
                const Matrix< DDRMat > & aNormal,
                const moris::Cell< MSI::Dof_Type > & aTestDofTypes )
        {
            // get test dof type index
            uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

            // get the conductivity property
            std::shared_ptr< Property > tPropConductivity =
                    mProperties( static_cast< uint >( Property_Type::CONDUCTIVITY ) );

            // compute test traction
            mTestTraction( tTestDofIndex ) =
                    trans( mFIManager->get_field_interpolators_for_type( mTempDof )->dnNdxn( 1 ) ) *
                    tPropConductivity->val()( 0 ) * aNormal;
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic::eval_strain()
        {
            // compute strain
            mStrain = mFIManager->get_field_interpolators_for_type( mTempDof )->gradx( 1 );
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic::eval_divstrain()
        {
            // get the temperature gradient
            Matrix< DDRMat > tTempGrad =
                    mFIManager->get_field_interpolators_for_type( mTempDof )->gradx( 2 );

            // evaluate the divergence of the strain
            mDivStrain = sum( tTempGrad( { 0, mSpaceDim - 1 }, { 0, 0 } ) );
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic::eval_testStrain()
        {
            // compute test strain
            mTestStrain = mFIManager->get_field_interpolators_for_type( mTempDof )->dnNdxn( 1 );
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic::eval_const()
        {
            // build an identity matrix
            Matrix< DDRMat > I;
            eye( mSpaceDim, mSpaceDim, I );

            // compute conductivity matrix
            mConst = mProperties( static_cast< uint >( Property_Type::CONDUCTIVITY ) )->val()( 0 ) * I;
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic::eval_dFluxdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the conductivity property
            std::shared_ptr< Property > tPropConductivity =
                    mProperties( static_cast< uint >( Property_Type::CONDUCTIVITY ) );

            // initialize the matrix
            mdFluxdDof( tDofIndex ).set_size( mSpaceDim, mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) )->
                    get_number_of_space_time_coefficients(), 0.0 );

            // get temperature FI
            Field_Interpolator * tFITemp =
                    mFIManager->get_field_interpolators_for_type( mTempDof );

            // if direct dependency on the dof type
            if( aDofTypes( 0 ) == mTempDof )
            {
                // compute derivative with direct dependency
                mdFluxdDof( tDofIndex ) +=
                        tPropConductivity->val()( 0 ) * tFITemp->dnNdxn( 1 );
            }

            // if indirect dependency on the dof type
            if ( mProperties( static_cast< uint >( Property_Type::CONDUCTIVITY ) )->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mdFluxdDof( tDofIndex ) +=
                        tFITemp->gradx( 1 ) * tPropConductivity->dPropdDOF( aDofTypes );
            }
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic::eval_dEnergyDotdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get properties
            std::shared_ptr< Property > tPropDensity = mProperties( static_cast< uint >( Property_Type::DENSITY ) );
            std::shared_ptr< Property > tPropHeatCap = mProperties( static_cast< uint >( Property_Type::HEAT_CAPACITY ) );

            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the corresponding FI
            Field_Interpolator * tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // initialize the matrix
            mEnergyDotDof( tDofIndex ).set_size( 1, tFI->get_number_of_space_time_coefficients(), 0.0 );

            // check if density and heat capacity are set
            if (tPropDensity == nullptr || tPropHeatCap == nullptr)
            {
                return;
            }

            // get the temperature FI
            Field_Interpolator * tFITemp =
                    mFIManager->get_field_interpolators_for_type( mTempDof );

            // if direct dependency on the dof type
            if( aDofTypes( 0 ) == mTempDof )
            {
                // compute derivative with direct dependency
                mEnergyDotDof( tDofIndex ) +=
                        tPropDensity->val()( 0 ) *
                        tPropHeatCap->val()( 0 ) *
                        tFITemp->dnNdtn(1);
            }

            // if indirect dependency of density on the dof type
            if ( tPropDensity->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mEnergyDotDof( tDofIndex ) +=
                        tPropHeatCap->val()( 0 ) *
                        tFITemp->gradt(1) *
                        tPropDensity->dPropdDOF( aDofTypes );
            }

            // if indirect dependency of heat capacity on the dof type
            if ( tPropHeatCap->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mEnergyDotDof( tDofIndex ) +=
                        tPropDensity->val()( 0 ) *
                        tFITemp->gradt(1) *
                        tPropHeatCap->dPropdDOF( aDofTypes );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic::eval_dGradEnergydDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get properties
            std::shared_ptr< Property > tPropDensity = mProperties( static_cast< uint >( Property_Type::DENSITY ) );
            std::shared_ptr< Property > tPropHeatCap = mProperties( static_cast< uint >( Property_Type::HEAT_CAPACITY ) );

            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the corresponding FI
            Field_Interpolator * tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // initialize the matrix
            mGradEnergyDof( tDofIndex ).set_size( mSpaceDim, tFI->get_number_of_space_time_coefficients(), 0.0 );

            // check if density and heat capacity are set
            if ( tPropDensity == nullptr || tPropHeatCap == nullptr )
            {
                return;
            }

            // temperature dof type
            Field_Interpolator * tFITemp =
                    mFIManager->get_field_interpolators_for_type( mTempDof );

            // if direct dependency on the dof type
            if( aDofTypes( 0 ) == mTempDof )
            {
                // compute derivative with direct dependency
                mGradEnergyDof( tDofIndex ) +=
                        tPropDensity->val()( 0 ) *
                        tPropHeatCap->val()( 0 ) *
                        tFITemp->dnNdxn(1);
            }

            // if indirect dependency of density on the dof type
            if ( tPropDensity->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mGradEnergyDof( tDofIndex ) +=
                        tPropHeatCap->val()( 0 ) *
                        tFITemp->gradx(1) *
                        tPropDensity->dPropdDOF( aDofTypes );
            }

            // if indirect dependency of heat capacity on the dof type
            if ( tPropHeatCap->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mGradEnergyDof( tDofIndex ) +=
                        tPropDensity->val()( 0 ) *
                        tFITemp->gradx(1) *
                        tPropHeatCap->dPropdDOF( aDofTypes );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic::eval_dGradEnergyDotdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get properties
            std::shared_ptr< Property > tPropDensity = mProperties( static_cast< uint >( Property_Type::DENSITY ) );
            std::shared_ptr< Property > tPropHeatCap = mProperties( static_cast< uint >( Property_Type::HEAT_CAPACITY ) );

            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the corresponding FI
            Field_Interpolator * tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // initialize the matrix
            mGradEnergyDotDof( tDofIndex ).set_size( mSpaceDim, tFI->get_number_of_space_time_coefficients(), 0.0 );

            // check if density and heat capacity are set
            if ( tPropDensity == nullptr || tPropHeatCap == nullptr )
            {
                return;
            }

            // temperature dof type
            Field_Interpolator * tFITemp =
                    mFIManager->get_field_interpolators_for_type( mTempDof );

            // if direct dependency on the dof type
            if( aDofTypes( 0 ) == mTempDof )
            {
                // compute derivative with direct dependency
                mGradEnergyDotDof( tDofIndex ) +=
                        tPropDensity->val()( 0 ) *
                        tPropHeatCap->val()( 0 ) *
                        tFITemp->d2Ndxt();
            }

            // if indirect dependency of density on the dof type
            if ( tPropDensity->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mGradEnergyDotDof( tDofIndex ) +=
                        tPropHeatCap->val()( 0 ) *
                        tFITemp->gradxt() *
                        tPropDensity->dPropdDOF( aDofTypes );
            }

            // if indirect dependency of heat capacity on the dof type
            if ( tPropHeatCap->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mGradEnergyDotDof( tDofIndex ) +=
                        tPropDensity->val()( 0 ) *
                        tFITemp->gradxt() *
                        tPropHeatCap->dPropdDOF( aDofTypes );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic::eval_dGradDivFluxdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {

            moris::real tK = mProperties( static_cast< uint >( Property_Type::CONDUCTIVITY ) )->val()( 0 );

            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the corresponding FI
            Field_Interpolator * tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // initialize results with zeros
            mGradDivFluxDof( tDofIndex ).set_size( mSpaceDim, tFI->get_number_of_space_time_coefficients(), 0.0 );

            // get spatial interpolation order
            mtk::Interpolation_Order tInterpOrder = tFI->get_space_interpolation_order();

            // only compute if interpolation is not 1 or 0, in that case simply keep the vector of zeros
            if (tInterpOrder != mtk::Interpolation_Order::LINEAR && tInterpOrder != mtk::Interpolation_Order::CONSTANT)
            {
                // matrix for purely isotropic case
                // FIXME: this implementation is slow and doesn't support indirect DoF dependency of conductivity
                Matrix< DDRMat > tKijIsotropic;
                switch( mSpaceDim )
                {
                    case 2:
                    {
                        tKijIsotropic = {
                                {tK,  0,  0, tK},
                                { 0, tK, tK,  0}};
                        break;
                    }
                    case 3:
                    {
                        tKijIsotropic = {
                                {tK, 0, 0, 0, 0,tK, 0,tK, 0, 0},
                                { 0,tK, 0,tK, 0, 0, 0, 0,tK, 0},
                                { 0, 0,tK, 0,tK, 0,tK, 0, 0, 0}};
                        break;
                    }
                    default:
                        MORIS_ASSERT(false, "CM_Diffusion_Linear_Isotropic::eval_dGradDivFluxdDOF: Number of spatial dimensions must be 2 or 3");
                        break;
                }

                // FIXME: indirect dependencies missing, spatial derivatives of properties needed
                // if direct dependency on the dof type
                if( aDofTypes( 0 ) == mTempDof )
                {
                    // compute derivative with direct dependency
                    mGradDivFluxDof( tDofIndex ) +=
                            tKijIsotropic *
                            mFIManager->get_field_interpolators_for_type( mTempDof )->dnNdxn( 3 );
                }
            }
        }
        //--------------------------------------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic::eval_ddivfluxdu( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the corresponding FI
            Field_Interpolator * tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set size for ddivflux/du
            mddivfluxdu( tDofIndex ).set_size( 1, tFI->get_number_of_space_time_coefficients(), 0.0 );

            // get the conductivity property
            std::shared_ptr< Property > tPropConductivity =
                    mProperties( static_cast< uint >( Property_Type::CONDUCTIVITY ) );

            // if temperature dof
            if( aDofTypes( 0 ) == mTempDof )
            {
                // fill ddivstrain/dv
                mddivfluxdu( tDofIndex ) += tPropConductivity->val()( 0 ) * this->ddivstraindu( aDofTypes );
            }

            if( tPropConductivity->check_dof_dependency( aDofTypes ) )
            {
                // fill ddivstrain/du
                mddivfluxdu( tDofIndex ) += this->divstrain() * tPropConductivity->dPropdDOF( aDofTypes );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic::eval_ddivstraindu( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the dof type FI
            Field_Interpolator * tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set size for ddivstrain/du
            mddivstraindu( tDofIndex ).set_size( 1, tFI->get_number_of_space_time_coefficients(), 0.0 );

            // if temperature dof type
            if( aDofTypes( 0 ) == mTempDof )
            {
                // get the 2nd order derivative of the shape functions d2Ndx2
                Matrix< DDRMat > tTempd2Ndx2 = tFI->dnNdxn( 2 );

                // fill ddivstrain/du
                mddivstraindu( tDofIndex ) = tTempd2Ndx2.get_row( 0 ) + tTempd2Ndx2.get_row( 1 );

                if( tTempd2Ndx2.n_rows() == 6 )
                {
                    mddivstraindu( tDofIndex ) += tTempd2Ndx2.get_row( 2 );
                }
            }
        }

        //------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic::eval_dTractiondDOF(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
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

        void CM_Diffusion_Linear_Isotropic::eval_dTestTractiondDOF(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                const Matrix< DDRMat >             & aNormal,
                const moris::Cell< MSI::Dof_Type > & aTestDofTypes )
        {
            // get test dof type index
            uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // initialize the dTestTractiondDof
            mdTestTractiondDof( tTestDofIndex )( tDofIndex ).set_size(
                    mFIManager->get_field_interpolators_for_type( aTestDofTypes( 0 ) )->get_number_of_space_time_coefficients(),
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) )->get_number_of_space_time_coefficients(), 0.0 );

            // get the conductivity property
            std::shared_ptr< Property > tPropConductivity =
                    mProperties( static_cast< uint >( Property_Type::CONDUCTIVITY ) );


            // if conductivity depends on dof type
            if( tPropConductivity->check_dof_dependency( aDofTypes ) )
            {
                // add contribution
                mdTestTractiondDof( tTestDofIndex )( tDofIndex ) +=
                        trans( mFIManager->get_field_interpolators_for_type( mTempDof )->dnNdxn( 1 ) ) *
                        aNormal * tPropConductivity->dPropdDOF( aDofTypes );
            }
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic::eval_dTestTractiondDOF(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
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

            // initialize the dTestTractiondDof
            mdTestTractiondDof( tTestDofIndex )( tDofIndex ).set_size(
                    mFIManager->get_field_interpolators_for_type( aTestDofTypes( 0 ) )->get_number_of_space_time_coefficients(),
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) )->get_number_of_space_time_coefficients(), 0.0 );

            // get the conductivity property
            std::shared_ptr< Property > tPropConductivity = mProperties( static_cast< uint >( Property_Type::CONDUCTIVITY ) );

            // if conductivity depends on dof type
            if( tPropConductivity->check_dof_dependency( aDofTypes ) )
            {
                // add contribution
                mdTestTractiondDof( tTestDofIndex )( tDofIndex ) +=
                        trans( mFIManager->get_field_interpolators_for_type( mTempDof )->dnNdxn( 1 ) ) *
                        aNormal * tPropConductivity->dPropdDOF( aDofTypes );
            }
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic::eval_dStraindDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the derivative dof type FI
            Field_Interpolator * tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init the matrix for dStraindDof
            mdStraindDof( tDofIndex ).set_size(
                    mSpaceDim,
                    tFIDer->get_number_of_space_time_coefficients(),
                    0.0 );

            // if direct dependency on the dof type
            if( aDofTypes( 0 ) == mTempDof )
            {
                // compute derivative with direct dependency
                mdStraindDof( tDofIndex ) = tFIDer->dnNdxn( 1 );
            }
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic::eval_dConstdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the derivative dof type FI
            Field_Interpolator * tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // reset the matrix
            mdConstdDof( tDofIndex ).set_size(
                    1,
                    tFIDer->get_number_of_space_time_coefficients(),
                    0.0 );

            // if indirect dependency on the dof type
            if ( mProperties( static_cast< uint >( Property_Type::CONDUCTIVITY ) )->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mdConstdDof( tDofIndex ) = mProperties( static_cast< uint >( Property_Type::CONDUCTIVITY ) )->dPropdDOF( aDofTypes );
            }
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic::eval_dFluxdDV( const moris::Cell< PDV_Type > & aDvTypes )
        {
            MORIS_ASSERT( false, " CM_Diffusion_Linear_Isotropic::eval_dFluxdDV - This function is not implemented.");
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic::eval_dStraindDV( const moris::Cell< PDV_Type > & aDvTypes )
        {
            MORIS_ASSERT( false, " CM_Diffusion_Linear_Isotropic::eval_dStraindDV - This function is not implemented.");
        }

        //------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
