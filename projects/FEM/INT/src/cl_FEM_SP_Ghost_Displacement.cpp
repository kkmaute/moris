
#include "cl_FEM_SP_Ghost_Displacement.hpp" //FEM/INT/src
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
//------------------------------------------------------------------------------
        void SP_Ghost_Displacement::eval_SP()
        {
            // compute stabilization parameter value
            mPPVal = mParameters( 0 ) * std::pow( mElementSize, 2 * ( mParameters( 1 )( 0 ) - 1 ) + 1 ) * mMasterProp( 0 )->val()( 0 );
        }

//------------------------------------------------------------------------------
        void SP_Ghost_Displacement::eval_dSPdMasterDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mMasterGlobalDofTypeMap( tDofType );

            // reset the matrix
            mdPPdMasterDof( tDofIndex ).set_size( 1, mMasterDofFI( tDofIndex )->get_number_of_space_time_coefficients(), 0.0 );

            // if indirect dependency on the dof type
            for( uint iProp = 0; iProp < mMasterProp.size(); iProp++ )
            {
                if ( mMasterProp( iProp )->check_dof_dependency( aDofTypes ) )
                {
                    // compute derivative with indirect dependency through properties
                    mdPPdMasterDof( tDofIndex ).matrix_data()
                    += this->val()( 0 ) * mMasterProp( iProp )->dPropdDOF( aDofTypes ) / mMasterProp( iProp )->val()( 0 );
                }
            }
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */


