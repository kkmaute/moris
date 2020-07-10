/*
 * cl_FEM_IQI_Latent_Heat_Absorption.cpp
 *
 *  Created on: Jul 10, 2020
 *      Author: wunsch
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Latent_Heat_Absorption.hpp"
#include "fn_FEM_CM_Phase_State_Functions.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        IQI_Latent_Heat_Absorption::IQI_Latent_Heat_Absorption()
        {
            // set IQI type
            mIQIType = vis::Output_Type::LATENT_HEAT_ABSORPTION;

            // set fem IQI type
            mFEMIQIType = fem::IQI_Type::LATENT_HEAT_ABSORPTION;

            // set the property pointer cell size
            mMasterProp.resize( static_cast< uint >( IQI_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Density" ]             = IQI_Property_Type::DENSITY;
            mPropertyMap[ "LatentHeat" ]          = IQI_Property_Type::LATENT_HEAT;
            mPropertyMap[ "PCTemp" ]              = IQI_Property_Type::PC_TEMP;
            mPropertyMap[ "PhaseStateFunction" ]  = IQI_Property_Type::PHASE_STATE_FUNCTION;
            mPropertyMap[ "PhaseChangeConst" ]    = IQI_Property_Type::PHASE_CHANGE_CONST;
        }

        //------------------------------------------------------------------------------

        void IQI_Latent_Heat_Absorption::set_property(
                std::shared_ptr< Property > aProperty,
                std::string                 aPropertyString,
                mtk::Master_Slave           aIsMaster)
        {
            // check that aPropertyString makes sense
            if ( mPropertyMap.find( aPropertyString ) == mPropertyMap.end() )
            {
                std::string tErrMsg =
                        std::string( "CM_Diffusion_Linear_Isotropic_Phase_Change::set_property - Unknown aPropertyString : ") +
                        aPropertyString;

                MORIS_ERROR( false , tErrMsg.c_str() );
            }

            // set the property in the property cell
            mMasterProp( static_cast< uint >( mPropertyMap[ aPropertyString ] ) ) = aProperty;
        }

        //------------------------------------------------------------------------------

        void IQI_Latent_Heat_Absorption::compute_QI( Matrix< DDRMat > & aQI )
        {
            moris::real tDensity = mMasterProp( static_cast< uint >( IQI_Property_Type::DENSITY ) )->val()( 0 );
            moris::real tLatHeat = mMasterProp( static_cast< uint >( IQI_Property_Type::LATENT_HEAT ) )->val()( 0 );

            // get the temperature FI
            Field_Interpolator * tFITemp =
                    mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::TEMP );

            // compute derivative of Phase State Function
            // real tdfdT = this->eval_dFdTemp();
            real tdfdT = eval_dFdTemp(
                    mMasterProp( static_cast< uint >( IQI_Property_Type::PC_TEMP ) )->val()( 0 ),
                    mMasterProp( static_cast< uint >( IQI_Property_Type::PHASE_CHANGE_CONST ) )->val()( 0 ),
                    mMasterProp( static_cast< uint >( IQI_Property_Type::PHASE_STATE_FUNCTION ) )->val()( 0 ),
                    tFITemp );

            // evaluate the QI
            aQI = tDensity * tLatHeat * tdfdT * tFITemp->gradt( 1 );
        }

        //------------------------------------------------------------------------------

        void IQI_Latent_Heat_Absorption::compute_QI( moris::real aWStar )
        {
            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // get properties
            moris::real tDensity = mMasterProp( static_cast< uint >( IQI_Property_Type::DENSITY ) )->val()( 0 );
            moris::real tLatHeat = mMasterProp( static_cast< uint >( IQI_Property_Type::LATENT_HEAT ) )->val()( 0 );

            // get the temperature FI
            Field_Interpolator * tFITemp =
                    mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::TEMP );

            // compute derivative of Phase State Function
            // real tdfdT = this->eval_dFdTemp();
            real tdfdT = eval_dFdTemp(
                    mMasterProp( static_cast< uint >( IQI_Property_Type::PC_TEMP ) )->val()( 0 ),
                    mMasterProp( static_cast< uint >( IQI_Property_Type::PHASE_CHANGE_CONST ) )->val()( 0 ),
                    mMasterProp( static_cast< uint >( IQI_Property_Type::PHASE_STATE_FUNCTION ) )->val()( 0 ),
                    tFITemp );

            // evaluate the QI
            mSet->get_QI()( tQIIndex ).matrix_data() += aWStar * tDensity * tLatHeat * tdfdT * tFITemp->gradt( 1 );
        }

        //------------------------------------------------------------------------------

        // FIXME: this functionality has not been tested, yet
        void IQI_Latent_Heat_Absorption::compute_dQIdu( real aWStar )
        {

            // get the column index to assemble in residual
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // get the requested dof types
            moris::Cell < enum MSI::Dof_Type > tRequestedDofTypes =
                    this->get_requested_dof_types();

            // get properties
            std::shared_ptr< Property > tPropDensity = mMasterProp( static_cast< uint >( IQI_Property_Type::DENSITY ) );
            std::shared_ptr< Property > tPropLatHeat = mMasterProp( static_cast< uint >( IQI_Property_Type::LATENT_HEAT ) );
            std::shared_ptr< Property > tPropPCtemp  = mMasterProp( static_cast< uint >( IQI_Property_Type::PC_TEMP ) );
            std::shared_ptr< Property > tPropPCconst = mMasterProp( static_cast< uint >( IQI_Property_Type::PHASE_CHANGE_CONST ) );
            std::shared_ptr< Property > tPropPSfunct = mMasterProp( static_cast< uint >( IQI_Property_Type::PHASE_STATE_FUNCTION ) );

            // get the temperature FI
            Field_Interpolator * tFITemp =
                    mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::TEMP );

            // compute derivative of Phase State Function
            // real tdfdT = this->eval_dFdTemp();
            real tdfdT = eval_dFdTemp(
                    tPropPCtemp->val()( 0 ),
                    tPropPCconst->val()( 0 ),
                    tPropPSfunct->val()( 0 ),
                    tFITemp );

            // compute dQIdDof for indirect dof dependencies
            for( uint iDof = 0; iDof < tRequestedDofTypes.size(); iDof++ )
            {
                // get treated dof type
                MSI::Dof_Type tDofType = tRequestedDofTypes( iDof );

                // get the set index for dof type
                sint tDofIndex = mSet->get_dof_index_for_type( tDofType, mtk::Master_Slave::MASTER );

                // get start and end indices for assembly
                uint tStartRow = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 );
                uint tEndRow   = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 );

                // if direct dependency on the dof type
                if( tDofType == MSI::Dof_Type::TEMP )
                {
                    // compute Dof derivative of phase state function
                    const moris::Matrix<DDRMat> dfdDof = eval_dFdTempdDOF(
                            tPropPCtemp->val()(0),
                            tPropPCconst->val()(0),
                            tPropPSfunct->val()(0),
                            tFITemp);

                    // compute derivative with direct dependency
                    mSet->get_residual()( tQIIndex )( { tStartRow, tEndRow }, { 0, 0 } ) +=
                            tPropDensity->val()(0) *  tPropLatHeat->val()(0) * tdfdT * tFITemp->dnNdtn(1)
                            + tPropDensity->val()(0) * tPropLatHeat->val()(0) * tFITemp->gradt(1) * dfdDof;
                }

                // if indirect dependency of density on the dof type
                if ( tPropDensity->check_dof_dependency( {MSI::Dof_Type::TEMP} ) )
                {
                    // compute derivative with indirect dependency through properties
                    mSet->get_residual()( tQIIndex )( { tStartRow, tEndRow }, { 0, 0 } ) +=
                            tPropLatHeat->val()(0) * tdfdT * tFITemp->gradt(1) * tPropDensity->dPropdDOF( {MSI::Dof_Type::TEMP} );
                }

            } /* end: for each requested DoF type */

        } /* end: method compute_dQIdu */

        //------------------------------------------------------------------------------
    }/* end_namespace_fem */
}/* end_namespace_moris */



