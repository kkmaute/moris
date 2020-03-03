/*
 * cl_FEM_IQI_Volume_Fraction.cpp
 *
 *  Created on: Feb 28, 2020
 *      Author: schmidt
 */

#include "cl_FEM_IQI_Volume_Fraction.hpp"

namespace moris
{
    namespace fem
    {

//------------------------------------------------------------------------------
        void IQI_Volume_Fraction::compute_QI( Matrix< DDRMat > & aQI )
        {
            // set the size for the QI
            aQI.set_size( 1, 1 );

            // evaluate the QI
            aQI( 0 ) = 1 * mStabilizationParam( 0 )->val()( 0 );
        }

        void IQI_Volume_Fraction::compute_QI( moris::real aWStar )
        {
            // get indices for properties, CM and SP
            //uint tElastIndex = static_cast< uint >( IQI_Constitutive_Type::ELAST );

            // get index for QI
            //sint tQIIndex = mSet->get_QI_assembly_map()( static_cast< uint >( mIQIMatType ) )( static_cast< uint >( mFEMIQIType ) );

            //print( mSet->get_QI(), "mSet->get_QI()");
            // evaluate the QI
            //mSet->get_QI()( tQIIndex ).matrix_data() += aWStar * trans( mMasterCM( tElastIndex )->flux() ) * mMasterCM( tElastIndex )->strain();
        }

//------------------------------------------------------------------------------
    }/* end namespace fem */
}/* end namespace moris */


