
#include "cl_FEM_SP_Dirichlet_Nitsche.hpp" //FEM/INT/src
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
        SP_Dirichlet_Nitsche::SP_Dirichlet_Nitsche()
        {
            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( SP_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Material" ] = SP_Property_Type::MATERIAL;
        }

//------------------------------------------------------------------------------
        void SP_Dirichlet_Nitsche::eval_SP()
        {
            // compute stabilization parameter value
            mPPVal = mParameters( 0 ) * mMasterProp( 0 )->val() / mElementSize;
        }

//------------------------------------------------------------------------------
        void SP_Dirichlet_Nitsche::eval_dSPdMasterDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
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
                    mdPPdMasterDof( tDofIndex ).matrix_data() += this->val()( 0 ) / mMasterProp( iProp )->val()( 0 ) * mMasterProp( iProp )->dPropdDOF( aDofTypes ) ;
                }
            }
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
