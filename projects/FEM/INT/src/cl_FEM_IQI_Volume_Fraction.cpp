/*
 * cl_FEM_IQI_Volume_Fraction.cpp
 *
 *  Created on: Feb 28, 2020
 *      Author: schmidt
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Volume_Fraction.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        IQI_Volume_Fraction::IQI_Volume_Fraction()
        {
            // set IQI type
            mIQIType = vis::Output_Type::VOLUME_FRACTION;

            // set fem IQI type
            mFEMIQIType = fem::IQI_Type::VOLUME_FRACTION;

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IQI_Stabilization_Type::MAX_ENUM ), nullptr );

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IQI_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mStabilizationMap[ "Reciprocal_total_vol" ] = IQI_Stabilization_Type::RECIPROCAL_TOTAL_VOLUME;
        }

        //------------------------------------------------------------------------------

        void IQI_Volume_Fraction::compute_QI( Matrix< DDRMat > & aQI )
        {
            // set the size for the QI
            aQI.set_size( 1, 1 );

            // evaluate the QI
            aQI( 0 ) = mStabilizationParam( 0 )->val()( 0 );
        }

        //------------------------------------------------------------------------------

        void IQI_Volume_Fraction::compute_dQIdu( MSI::Dof_Type aDofType, Matrix< DDRMat > & adQIdu )
        {
            MORIS_ERROR(false, "compute_dQIdu not implemented for a volume fraction IQI.");
        }

        //------------------------------------------------------------------------------

    }/* end namespace fem */
}/* end namespace moris */


