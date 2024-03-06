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
            init_constitutive_model( "Elast", IQI_Constitutive_Type::ELAST_LIN_ISO );
            init_property( "EigenStrain", IQI_Property_Type::EIGEN_STRAIN );
        }

        //------------------------------------------------------------------------------

        void
        IQI_Homogenized_Constitutive::compute_QI( Matrix< DDRMat >& aQI )
        {
            // get the const. model pointer element corresponding to elasticity from the cell
            std::shared_ptr< fem::Constitutive_Model > const & tCMElasticity = get_leader_constitutive_model( IQI_Constitutive_Type::ELAST_LIN_ISO );

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
            sint tQIIndex = mSet->get_QI_assembly_index( get_name() );

            // get the const. model pointer element corresponding to elasticity from the cell
            std::shared_ptr< fem::Constitutive_Model > const & tCMElasticity = get_leader_constitutive_model( IQI_Constitutive_Type::ELAST_LIN_ISO );

            // Jacob Fish 2014, Practical Multiscaling, 2.2.1 Tow-Scale Formulation
            // note: strain = elastic strain

            // evaluate the QI
            mSet->get_QI()( tQIIndex ) +=
                    aWStar * (    //
                            trans( tCMElasticity->strain() ) * tCMElasticity->constitutive() * tCMElasticity->strain() );
        }
    }    // namespace fem
}    // namespace moris
