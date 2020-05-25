#include "cl_FEM_CM_Diffusion_Linear_Isotropic_Phase_Change.hpp"
#include "cl_FEM_CM_Diffusion_Linear_Isotropic.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include <iostream>

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"
#include "fn_sum.hpp"
#include "fn_FEM_CM_Phase_State_Functions.hpp"

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
            mPropertyMap[ "PC_Temp" ]               = Property_Type::PC_TEMP;
            mPropertyMap[ "Phase_State_Function" ]  = Property_Type::PHASE_STATE_FUNCTION;
            mPropertyMap[ "Phase_Change_Const" ]    = Property_Type::PHASE_CHANGE_CONST;

            // populate dof map
            mDofMap[ "Temp" ] = MSI::Dof_Type::TEMP;
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic_Phase_Change::set_property(
                 std::shared_ptr< fem::Property > aProperty,
                 std::string                      aPropertyString )
         {
             // check that aPropertyString makes sense
             if ( mPropertyMap.find( aPropertyString ) == mPropertyMap.end() )
             {
                 std::string tErrMsg =
                         std::string("CM_Diffusion_Linear_Isotropic_Phase_Change::set_property - Unknown aPropertyString : ") +
                         aPropertyString;

                 MORIS_ERROR( false , tErrMsg.c_str() );
             }

             // set the property in the property cell
             mProperties( static_cast< uint >( mPropertyMap[ aPropertyString ] ) ) = aProperty;
         }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic_Phase_Change::eval_Hdot()
        {
            // get properties
            moris::real tDensity = mProperties( static_cast< uint >( Property_Type::DENSITY ) )->val()( 0 );
            moris::real tHeatCap = mProperties( static_cast< uint >( Property_Type::HEAT_CAPACITY ) )->val()( 0 );
            moris::real tLatHeat = mProperties( static_cast< uint >( Property_Type::LATENT_HEAT ) )->val()( 0 );

            // temperature dof type
            MSI::Dof_Type tTempDof = mDofMap[ "Temp" ];

            // compute derivative of Phase State Function
            // real tdfdT = this->eval_dFdTemp();
            real tdfdT = eval_dFdTemp(
                    mProperties( static_cast< uint >( Property_Type::PC_TEMP ) )->val()( 0 ),
                    mProperties( static_cast< uint >( Property_Type::PHASE_CHANGE_CONST ) )->val()( 0 ),
                    mProperties( static_cast< uint >( Property_Type::PHASE_STATE_FUNCTION ) )->val()( 0 ),
                    mFIManager->get_field_interpolators_for_type( tTempDof ) );

            // compute derivative of enthalpy
            mHdot = tDensity * ( tHeatCap + tLatHeat * tdfdT ) *
                    mFIManager->get_field_interpolators_for_type( tTempDof )->gradt( 1 );
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic_Phase_Change::eval_gradHdot()
        {
            moris::real tDensity = mProperties( static_cast< uint >( Property_Type::DENSITY ) )->val()( 0 );
            moris::real tHeatCap = mProperties( static_cast< uint >( Property_Type::HEAT_CAPACITY ) )->val()( 0 );
            moris::real tLatHeat = mProperties( static_cast< uint >( Property_Type::LATENT_HEAT ) )->val()( 0 );

            // temperature dof type
            MSI::Dof_Type tTempDof = mDofMap[ "Temp" ];

            // compute derivative of Phase State Function
            // real tdfdT = this->eval_dFdTemp();
            real tdfdT = eval_dFdTemp(
                    mProperties( static_cast< uint >( Property_Type::PC_TEMP ) )->val()( 0 ),
                    mProperties( static_cast< uint >( Property_Type::PHASE_CHANGE_CONST ) )->val()( 0 ),
                    mProperties( static_cast< uint >( Property_Type::PHASE_STATE_FUNCTION ) )->val()( 0 ),
                    mFIManager->get_field_interpolators_for_type( tTempDof ) );

            // compute gradient of
            mGradHdot = tDensity * ( tHeatCap + tLatHeat * tdfdT ) *
                    mFIManager->get_field_interpolators_for_type( tTempDof )->gradxt();
        }


        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic_Phase_Change::eval_dHdotdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get properties
            std::shared_ptr< Property > tPropDensity = mProperties( static_cast< uint >( Property_Type::DENSITY ) );
            std::shared_ptr< Property > tPropHeatCap = mProperties( static_cast< uint >( Property_Type::HEAT_CAPACITY ) );
            std::shared_ptr< Property > tPropLatHeat = mProperties( static_cast< uint >( Property_Type::LATENT_HEAT ) );
            std::shared_ptr< Property > tPropPCtemp  = mProperties( static_cast< uint >( Property_Type::PC_TEMP ) );
            std::shared_ptr< Property > tPropPCconst = mProperties( static_cast< uint >( Property_Type::PHASE_CHANGE_CONST ) );
            std::shared_ptr< Property > tPropPSfunct = mProperties( static_cast< uint >( Property_Type::PHASE_STATE_FUNCTION ) );

            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the corresponding FI
            Field_Interpolator * tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // initialize the matrix
            mHdotDof( tDofIndex ).set_size( 1, tFI->get_number_of_space_time_coefficients(), 0.0 );

            // temperature dof type
            MSI::Dof_Type tTempDof = mDofMap[ "Temp" ];

            // compute derivative of Phase State Function
            // real tdfdT = this->eval_dFdTemp();
            real tdfdT = eval_dFdTemp(
                    tPropPCtemp->val()( 0 ),
                    tPropPCconst->val()( 0 ),
                    tPropPSfunct->val()( 0 ),
                    mFIManager->get_field_interpolators_for_type( tTempDof ) );

            // if direct dependency on the dof type
            if( aDofTypes( 0 ) == tTempDof )
            {
                // compute derivative with direct dependency
                mHdotDof( tDofIndex ).matrix_data() +=
                        tPropDensity->val()(0) * ( tPropHeatCap->val()(0) + tPropLatHeat->val()(0) * tdfdT ) *
                        mFIManager->get_field_interpolators_for_type( tTempDof )->dnNdtn(1);
            }

            // if indirect dependency of density on the dof type
            if ( tPropDensity->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mHdotDof( tDofIndex ).matrix_data() +=
                        ( tPropHeatCap->val()(0) + tPropLatHeat->val()(0) * tdfdT ) *
                        mFIManager->get_field_interpolators_for_type( tTempDof )->gradt(1) *
                        tPropDensity->dPropdDOF( aDofTypes );
            }

            // if indirect dependency of heat capacity on the dof type
            if ( tPropHeatCap->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mHdotDof( tDofIndex ).matrix_data() +=
                        tPropDensity->val()( 0 ) *
                        mFIManager->get_field_interpolators_for_type( tTempDof )->gradt(1) *
                        tPropHeatCap->dPropdDOF( aDofTypes );
            }
        }

        //--------------------------------------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic_Phase_Change::eval_dGradHdotdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get properties
            std::shared_ptr< Property > tPropDensity = mProperties( static_cast< uint >( Property_Type::DENSITY ) );
            std::shared_ptr< Property > tPropHeatCap = mProperties( static_cast< uint >( Property_Type::HEAT_CAPACITY ) );
            std::shared_ptr< Property > tPropLatHeat = mProperties( static_cast< uint >( Property_Type::LATENT_HEAT ) );
            std::shared_ptr< Property > tPropPCtemp  = mProperties( static_cast< uint >( Property_Type::PC_TEMP ) );
            std::shared_ptr< Property > tPropPCconst = mProperties( static_cast< uint >( Property_Type::PHASE_CHANGE_CONST ) );
            std::shared_ptr< Property > tPropPSfunct = mProperties( static_cast< uint >( Property_Type::PHASE_STATE_FUNCTION ) );

            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the corresponding FI
            Field_Interpolator * tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // initialize the matrix
            mGradHdotDof( tDofIndex ).set_size( mSpaceDim, tFI->get_number_of_space_time_coefficients(), 0.0 );

            // temperature dof type
            MSI::Dof_Type tTempDof = mDofMap[ "Temp" ];

            // compute derivative of Phase State Function
            // real tdfdT = this->eval_dFdTemp();
            real tdfdT = eval_dFdTemp(
                    tPropPCtemp->val()( 0 ),
                    tPropPCconst->val()( 0 ),
                    tPropPSfunct->val()( 0 ),
                    mFIManager->get_field_interpolators_for_type( tTempDof ) );

            // if direct dependency on the dof type
            if( aDofTypes( 0 ) == tTempDof )
            {
                // compute derivative with direct dependency
                mGradHdotDof( tDofIndex ).matrix_data() +=
                        tPropDensity->val()(0) * ( tPropHeatCap->val()(0) + tPropLatHeat->val()(0) * tdfdT ) *
                        mFIManager->get_field_interpolators_for_type( tTempDof )->d2Ndxt();
            }

            // if indirect dependency of density on the dof type
            if ( tPropDensity->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mGradHdotDof( tDofIndex ).matrix_data() +=
                        ( tPropHeatCap->val()(0) + tPropLatHeat->val()(0) * tdfdT ) *
                        mFIManager->get_field_interpolators_for_type( tTempDof )->gradxt() *
                        tPropDensity->dPropdDOF( aDofTypes );
            }

            // if indirect dependency of heat capacity on the dof type
            if ( tPropHeatCap->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mGradHdotDof( tDofIndex ).matrix_data() +=
                        tPropDensity->val()( 0 ) *
                        mFIManager->get_field_interpolators_for_type( tTempDof )->gradxt() *
                        tPropHeatCap->dPropdDOF( aDofTypes );
            }
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
