/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_J_Integral.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_J_Integral.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IQI_J_Integral::IQI_J_Integral()
        {
            // set size for the constitutive model pointer cell
            mLeaderCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "ElastLinIso" ] = static_cast< uint >( IQI_Constitutive_Type::ELAST_LIN_ISO );
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
            aQI = ( trans( mLeaderCM( tElastLinIsoIndex )->flux() ) * mLeaderCM( tElastLinIsoIndex )->strain() )*tN(0)
                            - mLeaderCM( tElastLinIsoIndex )->traction(tN)(0)*mLeaderCM( tElastLinIsoIndex )->strain()(0)
                            - mLeaderCM( tElastLinIsoIndex )->traction(tN)(1)*mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX )->gradx(1)(1,0);  //FIXME: need to ask for displacement dof

        }

        //------------------------------------------------------------------------------

        void IQI_J_Integral::compute_QI( real aWStar )
        {
            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            /*
             * TODO: implement switch case for 2D or 3D
             */
            Matrix< DDRMat > tN = { {1.0}, {0.0} };

            // get indices for properties, CM and SP
            uint tElastLinIsoIndex = static_cast< uint >( IQI_Constitutive_Type::ELAST_LIN_ISO );

            /* evaluate the QI:
             * when crack grows straight ahead, G = J = int_{\Gamma} W n_1 - t_i \frac{\partial u_i}{\partial x_1}
             */

            // 2D
            // evaluate the QI
            mSet->get_QI()( tQIIndex ) += aWStar * (
                    trans( mLeaderCM( tElastLinIsoIndex )->flux() ) * mLeaderCM( tElastLinIsoIndex )->strain() )*tN(0)
                    - mLeaderCM( tElastLinIsoIndex )->traction(tN)(0)*mLeaderCM( tElastLinIsoIndex )->strain()(0)
                    - mLeaderCM( tElastLinIsoIndex )->traction(tN)(1)*mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX )->gradx(1)(1,0);  //FIXME: need to ask for displacement dof

        }

        //------------------------------------------------------------------------------

    }   // end fem namespace
}       // end moris namespace

