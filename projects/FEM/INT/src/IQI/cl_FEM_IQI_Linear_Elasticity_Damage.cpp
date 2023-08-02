/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Linear_Elasticity_Damage.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Linear_Elasticity_Damage.hpp"

#include "cl_FEM_CM_Struc_Linear_Isotropic_Damage.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IQI_Linear_Elasticity_Damage::IQI_Linear_Elasticity_Damage()
        {
            // set fem IQI type
            mFEMIQIType = fem::IQI_Type::LINEAR_ELASTICITY_DAMAGE;

            // set size for the constitutive model pointer cell
            mLeaderCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "ElasticDamage" ] =
                    static_cast< uint >( IQI_Constitutive_Type::ELASTIC_DAMAGE );
        }

        //------------------------------------------------------------------------------

        void
        IQI_Linear_Elasticity_Damage::compute_QI( Matrix< DDRMat >& aQI )
        {
            // get the elasticity with damage CM
            const std::shared_ptr< Constitutive_Model >& tCMElasticityDamage =
                    mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::ELASTIC_DAMAGE ) );

            // cast constitutive model base class pointer to elasticity damage constitutive model
            CM_Struc_Linear_Isotropic_Damage* tCMElasticityDamagePtr =
                    dynamic_cast< CM_Struc_Linear_Isotropic_Damage* >( tCMElasticityDamage.get() );

            // switch on type index
            switch ( mIQITypeIndex )
            {
                // equivalent strain
                case 0:
                {
                    // compute equivalent strain
                    aQI = tCMElasticityDamagePtr->equivalent_strain();
                    break;
                }
                // damage
                case 1:
                {
                    // compute damage
                    aQI = tCMElasticityDamagePtr->damage();
                    break;
                }
                    // smooth damage
                case 2:
                {
                    aQI = tCMElasticityDamagePtr->smooth_damage();
                    break;
                }
                case 3:
                {
                    MORIS_ERROR( mIQITypeIndex == -1,
                            "IQI_Linear_Elasticity_Damage::compute_QI - no specified dof, exiting!" );

                    // FIXME
                    aQI = tCMElasticityDamagePtr->traction( mNormal )( 0 );
                    break;
                }
                    // if none of the above
                default:
                {
                    MORIS_ERROR( false, "IQI_Linear_Elasticity_Damage::compute_QI - wrong mIQITypeIndex type" );
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        IQI_Linear_Elasticity_Damage::compute_QI( real aWStar )
        {
            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // init QI value matrix
            Matrix< DDRMat > tQI( 1, 1 );
            this->compute_QI( tQI );

            // evaluate the QI
            mSet->get_QI()( tQIIndex ) += aWStar * tQI( 0 );
        }

        //

        //------------------------------------------------------------------------------
    }    // namespace fem
}    // namespace moris

