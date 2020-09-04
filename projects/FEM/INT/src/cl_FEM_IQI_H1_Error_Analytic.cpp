/*
 * cl_FEM_IQI_H1_Error_Analytic.cpp
 *
 *  Created on: Feb 2, 2020
 *      Author: noel
 */
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_H1_Error_Analytic.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        IQI_H1_Error_Analytic::IQI_H1_Error_Analytic()
        {
            // set IQI type
            mIQIType = vis::Output_Type::H1_ERROR_ANALYTIC;

            // set FEM IQI type
            mFEMIQIType = fem::IQI_Type::H1_ERROR_ANALYTIC;

            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IQI_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "H1Check" ] = IQI_Property_Type::H1_CHECK;
        }

        //------------------------------------------------------------------------------

        void IQI_H1_Error_Analytic::set_property(
                std::shared_ptr< Property > aProperty,
                std::string                 aPropertyString,
                mtk::Master_Slave           aIsMaster )
        {
            // check that aPropertyString makes sense
            std::string tErrMsg =
                    std::string( "IQI_H1_Error_Analytic::set_property - Unknown aPropertyString: " ) +
                    aPropertyString;
            MORIS_ERROR( mPropertyMap.find( aPropertyString ) != mPropertyMap.end(), tErrMsg.c_str() );

            // check no slave allowed
            MORIS_ERROR( aIsMaster == mtk::Master_Slave::MASTER,
                    "IQI_H1_Error_Analytic::set_property - No slave allowed." );

            // set the property in the property cell
            this->get_properties( aIsMaster )( static_cast< uint >( mPropertyMap[ aPropertyString ] ) ) = aProperty;
        }

        //------------------------------------------------------------------------------

        void IQI_H1_Error_Analytic::compute_QI( Matrix< DDRMat > & aQI )
        {
            // get field interpolator
            Field_Interpolator * tFI =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofTypes( 0 )( 0 ) );

            // get analytical solution property
            std::shared_ptr< Property > tPropH1Check =
                    mMasterProp( static_cast< uint >( IQI_Property_Type::H1_CHECK ) );

            // get jump between value and analytic
            Matrix< DDRMat > tJump =
                    reshape( tFI->gradx( 1 ) - tPropH1Check->val(), tFI->gradx( 1 ).numel(), 1 );

            // evaluate the QI
            aQI = trans( tJump ) * tJump ;
        }

        //------------------------------------------------------------------------------

        void IQI_H1_Error_Analytic::compute_dQIdu(
                moris::Cell< MSI::Dof_Type > & aDofType,
                Matrix< DDRMat >             & adQIdu )
        {
            MORIS_ERROR( false, "IQI_H1_Error_Analytic::compute_dQIdu - Not implemented." );
        }

        //------------------------------------------------------------------------------
    }/* end_namespace_fem */
}/* end_namespace_moris */



