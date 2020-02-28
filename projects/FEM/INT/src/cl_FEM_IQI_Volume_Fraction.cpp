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

//------------------------------------------------------------------------------
    }/* end namespace fem */
}/* end namespace moris */


