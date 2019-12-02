/*
 * cl_FEM_IQI_Volume.cpp
 *
 *  Created on: Nov 20, 2019
 *      Author: noel
 */

#include "cl_FEM_IQI_Volume.hpp"

namespace moris
{
    namespace fem
    {

//------------------------------------------------------------------------------
        void IQI_Volume::compute_QI( Matrix< DDRMat > & aQI )
        {
            // set the size for the QI
            aQI.set_size( 1, 1 );

            // evaluate the QI
            aQI = {{ 1.0 }};
        }

//------------------------------------------------------------------------------
    }/* end namespace fem */
}/* end namespace moris */


