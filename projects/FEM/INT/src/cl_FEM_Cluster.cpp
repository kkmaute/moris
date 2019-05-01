/*
 * cl_FEM_Cluster.hpp
 *
 *  Created on: Apr 11, 2019
 *      Author: schmidt
 */
#include <iostream>

#include "cl_FEM_Element.hpp" //FEM/INT/src
#include "cl_FEM_Cluster.hpp"   //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

        void Cluster::compute_jacobian()
        {
            this->get_my_pdof_values();

            mElementBlock->initialize_mJacobianElement();

            this->set_field_interpolators_coefficients();

            for ( uint Ik = 0; Ik < mInterpElements.size(); Ik++ )
            {
                mInterpElements( Ik )->compute_jacobian();
            }
        }

//------------------------------------------------------------------------------

        void Cluster::compute_residual()
        {
            this->get_my_pdof_values();

            mElementBlock->initialize_mResidualElement();

            this->set_field_interpolators_coefficients();

            for ( uint Ik = 0; Ik < mInterpElements.size(); Ik++ )
            {
                mInterpElements( Ik )->compute_residual();
            }
        }

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
