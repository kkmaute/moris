
#include "cl_FEM_SP_Incompressible_Flow.hpp" //FEM/INT/src
#include "cl_FEM_Cluster.hpp"              //FEM/INT/src
#include "cl_FEM_Field_Interpolator_Manager.hpp"              //FEM/INT/src

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"
#include "op_div.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
        void SP_Incompressible_Flow::eval_SP()
        {
            MORIS_ERROR( false, "SP_Incompressible_Flow::eval_SP - Not implemented." );
        }

//------------------------------------------------------------------------------
        void SP_Incompressible_Flow::eval_dSPdMasterDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            MORIS_ERROR( false, "SP_Incompressible_Flow::eval_dSPdMasterDOF - Not implemented." );
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */


