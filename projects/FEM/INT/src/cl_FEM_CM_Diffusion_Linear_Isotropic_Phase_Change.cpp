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
        void CM_Diffusion_Linear_Isotropic_Phase_Change::eval_Hdot()
        {
            // get properties
            moris::real tDensity = mProperties( static_cast< uint >( Property_Type::DENSITY ) )->val()( 0 );
            moris::real tHeatCap = mProperties( static_cast< uint >( Property_Type::HEAT_CAPACITY ) )->val()( 0 );
            moris::real tLatHeat = mProperties( static_cast< uint >( Property_Type::LATENT_HEAT ) )->val()( 0 );

            // compute derivative of Phase State Function
            // real tdfdT = this->eval_dFdTemp();
            real tdfdT = eval_dFdTemp(
                    mProperties( static_cast< uint >( Property_Type::PC_TEMP ) )->val()( 0 ),
                    mProperties( static_cast< uint >( Property_Type::PHASE_CHANGE_CONST ) )->val()( 0 ),
                    mProperties( static_cast< uint >( Property_Type::PHASE_STATE_FUNCTION ) )->val()( 0 ),
                    mFIManager->get_field_interpolators_for_type( mDofMap[ "Temp" ] ) );

            // compute derivative of enthalpy
            mHdot = tDensity * ( tHeatCap + tLatHeat * tdfdT )
                             * mFIManager->get_field_interpolators_for_type( mDofMap[ "Temp" ] )->gradt( 1 );
        }

        //------------------------------------------------------------------------------
        void CM_Diffusion_Linear_Isotropic_Phase_Change::eval_gradHdot()
        {
            moris::real tDensity = mProperties( static_cast< uint >( Property_Type::DENSITY ) )->val()( 0 );
            moris::real tHeatCap = mProperties( static_cast< uint >( Property_Type::HEAT_CAPACITY ) )->val()( 0 );
            moris::real tLatHeat = mProperties( static_cast< uint >( Property_Type::LATENT_HEAT ) )->val()( 0 );

            // compute derivative of Phase State Function
            // real tdfdT = this->eval_dFdTemp();
            real tdfdT = eval_dFdTemp(
                    mProperties( static_cast< uint >( Property_Type::PC_TEMP ) )->val()( 0 ),
                    mProperties( static_cast< uint >( Property_Type::PHASE_CHANGE_CONST ) )->val()( 0 ),
                    mProperties( static_cast< uint >( Property_Type::PHASE_STATE_FUNCTION ) )->val()( 0 ),
                    mFIManager->get_field_interpolators_for_type( mDofMap[ "Temp" ] ) );

            // compute gradient of
            mGradHdot = tDensity * ( tHeatCap + tLatHeat * tdfdT ) *
                    mFIManager->get_field_interpolators_for_type( mDofMap[ "Temp" ] )->gradxt();
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
            // real tdfdT = this->eval_dFdTemp();
            real tdfdT = eval_dFdTemp(
                    mProperties( static_cast< uint >( Property_Type::PC_TEMP ) )->val()( 0 ),
                    mProperties( static_cast< uint >( Property_Type::PHASE_CHANGE_CONST ) )->val()( 0 ),
                    mProperties( static_cast< uint >( Property_Type::PHASE_STATE_FUNCTION ) )->val()( 0 ),
                    mFIManager->get_field_interpolators_for_type( mDofMap[ "Temp" ] ) );

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
            // real tdfdT = this->eval_dFdTemp();
            real tdfdT = eval_dFdTemp(
                    mProperties( static_cast< uint >( Property_Type::PC_TEMP ) )->val()( 0 ),
                    mProperties( static_cast< uint >( Property_Type::PHASE_CHANGE_CONST ) )->val()( 0 ),
                    mProperties( static_cast< uint >( Property_Type::PHASE_STATE_FUNCTION ) )->val()( 0 ),
                    mFIManager->get_field_interpolators_for_type( mDofMap[ "Temp" ] ) );

            // if direct dependency on the dof type
            if( aDofTypes( 0 ) == mDofMap[ "Temp" ] )
            {
                // compute derivative with direct dependency
                mGradHdotDof( tDofIndex ).matrix_data() +=
                        tDensity * ( tHeatCap + tLatHeat * tdfdT ) *
                        mFIManager->get_field_interpolators_for_type( mDofMap[ "Temp" ] )->d2Ndxt();
            }
        }



        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
