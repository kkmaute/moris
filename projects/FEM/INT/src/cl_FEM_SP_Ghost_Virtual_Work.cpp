
#include "cl_FEM_SP_Ghost_Virtual_Work.hpp" //FEM/INT/src
#include "cl_FEM_Cluster.hpp"              //FEM/INT/src

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"
#include "op_div.hpp"

namespace moris
{
    namespace fem
    {
        SP_Ghost_Virtual_Work::SP_Ghost_Virtual_Work()
        {
            // set the list of cluster measures
            mClusterMeasures = { fem::Cluster_Measure::ELEMENT_SIZE };
        }
//------------------------------------------------------------------------------
        void SP_Ghost_Virtual_Work::eval_SP()
        {
            // compute stabilization parameter value
            mPPVal = mParameters( 0 ) * std::pow( mElementSize, 2 * ( mParameters( 1 )( 0 ) - 1 ) + 1 );
        }

//------------------------------------------------------------------------------
        void SP_Ghost_Virtual_Work::eval_dSPdMasterDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mMasterGlobalDofTypeMap( tDofType );

            // reset the matrix
            mdPPdMasterDof( tDofIndex ).set_size( 1, mMasterDofFI( tDofIndex )->get_number_of_space_time_coefficients(), 0.0 );
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */


