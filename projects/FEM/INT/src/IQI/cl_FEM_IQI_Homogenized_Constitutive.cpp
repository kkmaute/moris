/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Homogenized_Constitutive.cpp
 *
 */

#include "cl_FEM_IQI_Homogenized_Constitutive.hpp"

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IQI_Homogenized_Constitutive::IQI_Homogenized_Constitutive()
        {
            // set fem IQI type
            mFEMIQIType = fem::IQI_Type::HOMOGENIZED_CONSTITUTIVE;

            // set size for the constitutive model pointer cell
            mLeaderCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "Elast" ] = static_cast< uint >( IQI_Constitutive_Type::ELAST_LIN_ISO );

            // set size for the property pointer cell
            mLeaderProp.resize( static_cast< uint >( IQI_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "EigenStrain" ] = static_cast< uint >( IQI_Property_Type::EIGEN_STRAIN );
        }

        //------------------------------------------------------------------------------

        void
        IQI_Homogenized_Constitutive::compute_QI( Matrix< DDRMat >& aQI )
        {
            // get the const. model pointer element corresponding to elasticity from the cell
            std::shared_ptr< fem::Constitutive_Model >& tCMElasticity =
                    mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::ELAST_LIN_ISO ) );

            // Jacob Fish 2014, Practical Multiscaling, 2.2.1 Tow-Scale Formulation
            // note: strain = elastic strain

            // assign value of QI;
            aQI = ( trans( tCMElasticity->strain() ) * tCMElasticity->constitutive() * tCMElasticity->strain() );
        }

        //------------------------------------------------------------------------------

        void
        IQI_Homogenized_Constitutive::compute_QI( real aWStar )
        {
            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // get the const. model pointer element corresponding to elasticity from the cell
            std::shared_ptr< fem::Constitutive_Model >& tCMElasticity =
                    mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::ELAST_LIN_ISO ) );

            // Jacob Fish 2014, Practical Multiscaling, 2.2.1 Tow-Scale Formulation
            // note: strain = elastic strain

            // evaluate the QI
            mSet->get_QI()( tQIIndex ) +=
                    aWStar * (    //
                            trans( tCMElasticity->strain() ) * tCMElasticity->constitutive() * tCMElasticity->strain() );
        }
    }    // namespace fem
}    // namespace moris
