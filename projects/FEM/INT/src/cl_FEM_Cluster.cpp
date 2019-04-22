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

            mInterpElements( 0 )->compute_jacobian();
        }

//------------------------------------------------------------------------------

        void Cluster::compute_residual()
        {
            this->get_my_pdof_values();

            mInterpElements( 0 )->compute_residual();
        }

//------------------------------------------------------------------------------

//        void Cluster::compute_jacobian_and_residual()
//        {
//            MORIS_ERROR( false, " Element::compute_jacobian_and_residual - not implemented. ");
//        }

//------------------------------------------------------------------------------



//------------------------------------------------------------------------------



//------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
