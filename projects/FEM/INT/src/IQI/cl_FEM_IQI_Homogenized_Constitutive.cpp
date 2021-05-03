/*
 * cl_FEM_IQI_Homogebized_Constitutive.cpp
 *
 *  Created on: Apr 10, 2021
 *      Author: momo
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
            mFEMIQIType = fem::IQI_Type::HOMOGENIZED_CONSTITUTIVE ;

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "Elast" ] = static_cast< uint >( IQI_Constitutive_Type::ELAST_LIN_ISO );

            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IQI_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "EigenStrain" ] = static_cast< uint >( IQI_Property_Type::EIGEN_STRAIN );
                }

        //------------------------------------------------------------------------------

        void IQI_Homogenized_Constitutive::compute_QI( Matrix< DDRMat > & aQI )
        {
            //get the const. model pointer element corresponding to elasticity from the cell
            std::shared_ptr< fem::Constitutive_Model > & tCMElasticity =
                    mMasterCM( static_cast< uint >( IQI_Constitutive_Type::ELAST_LIN_ISO ) );

            //get the property pointer element corresponding to eigen-strain from the cell
            std::shared_ptr< fem::Property > & tPropEigenStrain =
                    mMasterProp( static_cast< uint >( IQI_Property_Type::EIGEN_STRAIN ) );

            //calclaue the elastic strain influence function
            //Jacob Fish 2014, Practical Multiscaling, 2.2.1 Tow-Scale Formulation
            Matrix< DDRMat > tStrainInfluence = tPropEigenStrain->val() + tCMElasticity->strain();

            //assign value of QI
            aQI = (trans( tStrainInfluence) * tCMElasticity->constitutive() * tStrainInfluence);
        }

        //------------------------------------------------------------------------------

        void IQI_Homogenized_Constitutive::compute_QI( real aWStar )
        {
            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            //get the const. model pointer element corresponding to elasticity from the cell
            std::shared_ptr< fem::Constitutive_Model > & tCMElasticity =
                    mMasterCM( static_cast< uint >( IQI_Constitutive_Type::ELAST_LIN_ISO ) );

            //get the property pointer element corresponding to eigen-strain from the cell
            std::shared_ptr< fem::Property > & tPropEigenStrain =
                    mMasterProp( static_cast< uint >( IQI_Property_Type::EIGEN_STRAIN ) );

            //calclaue the elastic strain influence function
            //Jacob Fish 2014, Practical Multiscaling, 2.2.1 Tow-Scale Formulation
            Matrix< DDRMat > tStrainInfluence = tPropEigenStrain->val() + tCMElasticity->strain();

            // evaluate the QI
            mSet->get_QI()( tQIIndex ) += aWStar * ( trans(tStrainInfluence) * tCMElasticity->constitutive() * tStrainInfluence );
        }
    }/* end_namespace_fem */
}/* end_namespace_moris */


