#include "assert.hpp"
#include "cl_FEM_Element_HJ.hpp" //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        Element_HJ::Element_HJ( Cell< IWG * > aListOfIWGs,
                                Cell< MSI::Dof_Type > aListOfNodePDOFsTypes,
                                Cell< MSI::Dof_Type > aListOfElementPDOFsTypes )
        {
        }

//------------------------------------------------------------------------------

        void Element_HJ::compute_jacobian()
        {
            MORIS_ASSERT( false, " Element_HJ - compute_jacobian - not implemented. " );
        }

//------------------------------------------------------------------------------

        void Element_HJ::compute_residual()
        {
            MORIS_ASSERT( false, " Element_HJ - compute_residual - not implemented. " );
        }

//------------------------------------------------------------------------------

        void Element_HJ::compute_jacobian_and_residual()
        {
            MORIS_ASSERT( false, " Element_HJ - compute_jacobian_and_residual - not implemented. " );
        }
//------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
