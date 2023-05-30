/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_2D_Drag_Lift_Coefficient.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_2D_Drag_Lift_Coefficient.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IQI_Drag_Lift_Coefficient::IQI_Drag_Lift_Coefficient( sint aBeta )
        {
            // set the property pointer cell size
            mLeaderProp.resize( static_cast< uint >( Property_Type::MAX_ENUM ), nullptr );

            // populate the map
            // NOTE: the user can either directly give a reference pressure value,
            // NOTE: or give the density, velocity, and a reference lengthscale to normalize against the dynamic pressure
            mPropertyMap[ "RefDensity" ]     = static_cast< uint >( Property_Type::REF_DENSITY );
            mPropertyMap[ "RefVelocity" ]    = static_cast< uint >( Property_Type::REF_VELOCITY );
            mPropertyMap[ "LengthScale" ]    = static_cast< uint >( Property_Type::LENGTHSCALE );
            mPropertyMap[ "RefPressure" ]    = static_cast< uint >( Property_Type::REF_PRESSURE );

            // set size for the constitutive model pointer cell
            mLeaderCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "Fluid" ] =
                    static_cast< uint >( IQI_Constitutive_Type::FLUID );

            // set mBeta
            mBeta = aBeta;
        }

        //------------------------------------------------------------------------------

        void IQI_Drag_Lift_Coefficient::compute_QI( Matrix< DDRMat > & aQI )
        {
            // get the reference pressure property
            std::shared_ptr< fem::Property > tPropRefPressure = mLeaderProp( static_cast< uint >( Property_Type::REF_PRESSURE ) );

            // initialize variable storing the reference pressure
            real tRefPressure = 0.0;

            // check if user has defined the reference pressure directly
            if( tPropRefPressure != nullptr )
            {
                tRefPressure = tPropRefPressure->val()( 0 );
            }
            // otherwise compute dynamic pressure as reference
            else
            {
                // get the density property value
                std::shared_ptr< fem::Property > tPropRefDensity = mLeaderProp( static_cast< uint >( Property_Type::REF_DENSITY ) );

                // get the maximum velocity property value
                std::shared_ptr< fem::Property > tPropRefVelocity = mLeaderProp( static_cast< uint >( Property_Type::REF_VELOCITY ) );

                // get the diameter property value
                std::shared_ptr< fem::Property > tPropLengthScale = mLeaderProp( static_cast< uint >( Property_Type::LENGTHSCALE ) );

                // check that the properties are actually populated
                MORIS_ASSERT( tPropRefDensity != nullptr && tPropRefVelocity != nullptr && tPropLengthScale != nullptr,
                    "IQI_Drag_Lift_Coefficient::compute_QI() - Properties not fully populated. "
                    "User must specify either reference pressure, or give reference values for the density, velocity and length scale." );

                // compute the dynamic pressure as reference
                tRefPressure = 0.5 * tPropRefDensity->val()( 0 ) * tPropLengthScale->val()( 0 ) * tPropRefVelocity->val()( 0 ) * tPropRefVelocity->val()( 0 );
            }

            // check for zero in denominator
            MORIS_ERROR( std::abs( tRefPressure ) > MORIS_REAL_EPS,
                "IQI_Drag_Lift_Coefficient::compute_QI() - Zero reference pressure detected. Check properties. (must not be zero)" );

            // get the diffusion CM
            const std::shared_ptr< Constitutive_Model > & tFluidCM =
                    mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::FLUID ) );

            // compute QI
            if( mBeta == 1 ) // drag coefficient
            {
                aQI = {{ tFluidCM->traction( mNormal )( 0 ) / tRefPressure }};
            }
            else if( mBeta == -1 ) // lift coefficient
            {
                aQI = {{ tFluidCM->traction( mNormal )( 1 ) / tRefPressure }};
            }
            else
            {
                MORIS_ERROR( false, "IQI_Drag_Lift_Coefficient::compute_QI() - unknown mBeta." );
            }
        }

        //------------------------------------------------------------------------------

        void IQI_Drag_Lift_Coefficient::compute_QI( real aWStar )
        {
            // get the reference pressure property
            std::shared_ptr< fem::Property > tPropRefPressure = mLeaderProp( static_cast< uint >( Property_Type::REF_PRESSURE ) );

            // initialize variable storing the reference pressure
            real tRefPressure = 0.0;

            // check if user has defined the reference pressure directly
            if( tPropRefPressure != nullptr )
            {
                tRefPressure = tPropRefPressure->val()( 0 );
            }
            // otherwise compute dynamic pressure as reference
            else
            {
                // get the density property value
                std::shared_ptr< fem::Property > tPropRefDensity = mLeaderProp( static_cast< uint >( Property_Type::REF_DENSITY ) );

                // get the maximum velocity property value
                std::shared_ptr< fem::Property > tPropRefVelocity = mLeaderProp( static_cast< uint >( Property_Type::REF_VELOCITY ) );

                // get the diameter property value
                std::shared_ptr< fem::Property > tPropLengthScale = mLeaderProp( static_cast< uint >( Property_Type::LENGTHSCALE ) );

                // check that the properties are actually populated
                MORIS_ASSERT( tPropRefDensity != nullptr && tPropRefVelocity != nullptr && tPropLengthScale != nullptr,
                    "IQI_Drag_Lift_Coefficient::compute_QI() - Properties not fully populated. "
                    "User must specify either reference pressure, or give reference values for the density, velocity and length scale." );

                // compute the dynamic pressure as reference
                tRefPressure = 0.5 * tPropRefDensity->val()( 0 ) * tPropLengthScale->val()( 0 ) * tPropRefVelocity->val()( 0 ) * tPropRefVelocity->val()( 0 );
            }

            // check for zero in denominator
            MORIS_ERROR( std::abs( tRefPressure ) > MORIS_REAL_EPS,
                "IQI_Drag_Lift_Coefficient::compute_QI() - Zero reference pressure detected. Check properties. (must not be zero)" );

            // get the diffusion CM
            const std::shared_ptr< Constitutive_Model > & tFluidCM =
                    mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::FLUID ) );

            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            if( mBeta == 1 ) // new drag coefficient
            {
                mSet->get_QI()( tQIIndex ) += aWStar * ( tFluidCM->traction( mNormal )( 0 ) / tRefPressure );
            }
            else if( mBeta == -1 ) // new lift coefficient
            {
                mSet->get_QI()( tQIIndex ) += aWStar * ( tFluidCM->traction( mNormal )( 1 ) / tRefPressure );
            }
            else
            {
                MORIS_ERROR( false, "IQI_Drag_Lift_Coefficient::compute_QI() - unknown mBeta." );
            }
        }

        //------------------------------------------------------------------------------
    }/* end namespace fem */
}/* end namespace moris */

