/*
 * cl_FEM_IQI_J_Integral.cpp
 *
 *  Created on: Feb 11, 2020
 *      Author: sonne
 */
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_J_Integral.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        void IQI_J_Integral::set_constitutive_model(
                std::shared_ptr< Constitutive_Model > aConstitutiveModel,
                std::string                           aConstitutiveString,
                mtk::Master_Slave                     aIsMaster )
        {
            // can only be master
            MORIS_ERROR( aIsMaster == mtk::Master_Slave::MASTER,
                    "IQI::set_constitutive model - can only be master." );

            // FIXME check that constitutive string makes sense?

            // set the constitutive model in the constitutive model cell
            this->get_constitutive_models( aIsMaster )( static_cast< uint >( mConstitutiveMap[ aConstitutiveString ] ) ) = aConstitutiveModel;
        }

        //------------------------------------------------------------------------------

        IQI_J_Integral::IQI_J_Integral()
        {
            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "ElastLinIso" ] = IQI_Constitutive_Type::ELAST_LIN_ISO;
        }

        //------------------------------------------------------------------------------

        void IQI_J_Integral::compute_QI( Matrix< DDRMat > & aQI )
        {
            /*
             * TODO: implement switch case for 2D or 3D
             */
            Matrix< DDRMat > tN = { {1.0},
                    {0.0} };

            // get indices for properties, CM and SP
            uint tElastLinIsoIndex = static_cast< uint >( IQI_Constitutive_Type::ELAST_LIN_ISO );

            /* evaluate the QI:
             * when crack grows straight ahead, G = J = int_{\Gamma} W n_1 - t_i \frac{\partial u_i}{\partial x_1}
             */

            // 2D
            aQI = ( trans( mMasterCM( tElastLinIsoIndex )->flux() ) * mMasterCM( tElastLinIsoIndex )->strain() )*tN(0)
                            - mMasterCM( tElastLinIsoIndex )->traction(tN)(0)*mMasterCM( tElastLinIsoIndex )->strain()(0)
                            - mMasterCM( tElastLinIsoIndex )->traction(tN)(1)*mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX )->gradx(1)(1,0);  //FIXME: need to ask for displacement dof

        }

        //------------------------------------------------------------------------------

        void IQI_J_Integral::compute_dQIdu(
                moris::Cell< MSI::Dof_Type > & aDofType,
                Matrix< DDRMat >             & adQIdu  )
        {
            MORIS_ERROR( false, "IQI_J_Integral::compute_dQIdu() - this function does nothing for the J-Integral" );
        }

        //------------------------------------------------------------------------------

    }   // end fem namespace
}       // end moris namespace

