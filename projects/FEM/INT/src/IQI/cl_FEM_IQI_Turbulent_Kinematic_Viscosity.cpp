/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Turbulent_Kinematic_Viscosity.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Turbulent_Kinematic_Viscosity.hpp"
#include "fn_norm.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IQI_Turbulent_Kinematic_Viscosity::IQI_Turbulent_Kinematic_Viscosity()
        {
            // set fem IQI type
            mFEMIQIType = fem::IQI_Type::TURBULENT_KINEMATIC_VISCOSITY;

            // set the property pointer cell size
            mMasterProp.resize( static_cast< uint >( Property_Type::MAX_ENUM ), nullptr );

            // populate the map
            mPropertyMap[ "Density" ]          = static_cast< uint >( Property_Type::DENSITY );
            mPropertyMap[ "DynamicViscosity" ] = static_cast< uint >( Property_Type::KINEMATIC_VISCOSITY );
        }

        //------------------------------------------------------------------------------

        void IQI_Turbulent_Kinematic_Viscosity::set_dof_type_list(
                moris::Cell< moris::Cell< MSI::Dof_Type > > & aDofTypes,
                moris::Cell< std::string >                  & aDofStrings,
                mtk::Master_Slave                             aIsMaster )
        {
            // switch on master slave
            switch ( aIsMaster )
            {
                case mtk::Master_Slave::MASTER :
                {
                    // set dof type list
                    mMasterDofTypes = aDofTypes;

                    // loop on dof type
                    for( uint iDof = 0; iDof < aDofTypes.size(); iDof++ )
                    {
                        // get dof string
                        std::string tDofString = aDofStrings( iDof );

                        // get dof type
                        MSI::Dof_Type tDofType = aDofTypes( iDof )( 0 );

                        // if viscosity
                        if( tDofString == "Viscosity" )
                        {
                            mMasterDofViscosity = tDofType;
                        }
                        else
                        {
                            // create error message
                            std::string tErrMsg =
                                    std::string( "IQI_Turbulent_Kinematic_Viscosity::set_dof_type_list - Unknown aDofString : ") +
                                    tDofString;
                            MORIS_ERROR( false , tErrMsg.c_str() );
                        }
                    }
                    break;
                }

                case mtk::Master_Slave::SLAVE :
                {
                    // set dof type list
                    mSlaveDofTypes = aDofTypes;
                    break;
                }

                default:
                    MORIS_ERROR( false, "IQI_Turbulent_Kinematic_Viscosity::set_dof_type_list - unknown master slave type." );
            }
        }

        //------------------------------------------------------------------------------

        void IQI_Turbulent_Kinematic_Viscosity::compute_QI( Matrix< DDRMat > & aQI )
        {
            // get field interpolator for modified viscosity
            Field_Interpolator * tFIModViscosity =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofViscosity );

            // get the density property
            const std::shared_ptr< Property > & tPropDensity =
                    mMasterProp( static_cast< uint >( Property_Type::DENSITY ) );

            // get the kinematic viscosity property
            const std::shared_ptr< Property > & tPropKinViscosity =
                    mMasterProp( static_cast< uint >( Property_Type::KINEMATIC_VISCOSITY ) );

            // compute fv1
            real tFv1 = compute_fv1(
                    { mMasterDofViscosity },
                    mMasterFIManager,
                    tPropKinViscosity );

            // compute turbulent kinematic viscosity
            aQI = tPropDensity->val() * tFIModViscosity->val() * tFv1;
        }

        //------------------------------------------------------------------------------

        void IQI_Turbulent_Kinematic_Viscosity::compute_QI( real aWStar )
        {
            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // get field interpolator for modified viscosity
            Field_Interpolator * tFIModViscosity =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofViscosity );

            // get the density property
            const std::shared_ptr< Property > & tPropDensity =
                    mMasterProp( static_cast< uint >( Property_Type::DENSITY ) );

            // get the kinematic viscosity property
            const std::shared_ptr< Property > & tPropKinViscosity =
                    mMasterProp( static_cast< uint >( Property_Type::KINEMATIC_VISCOSITY ) );

            // compute fv1
            real tFv1 = compute_fv1(
                    { mMasterDofViscosity },
                    mMasterFIManager,
                    tPropKinViscosity );

            // compute turbulent kinematic viscosity
            mSet->get_QI()( tQIIndex ) += aWStar * ( tPropDensity->val() * tFIModViscosity->val() * tFv1 );
        }
        //------------------------------------------------------------------------------
    }/* end_namespace_fem */
}/* end_namespace_moris */

