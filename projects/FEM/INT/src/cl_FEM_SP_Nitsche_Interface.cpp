#include "cl_FEM_SP_Nitsche_Interface.hpp" //FEM/INT/src
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
        void SP_Nitsche_Interface::eval_SP()
        {
            // compute stabilization parameter value
            mPPVal = 2.0 * mParameters( 0 ) * mInterfaceSurface / ( mMasterVolume / mMasterProp( 0 )->val()( 0 ) + mSlaveVolume / mSlaveProp( 0 )->val()( 0 ) );
        }

//------------------------------------------------------------------------------
        void SP_Nitsche_Interface::eval_dSPdMasterDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
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
                    += this->val()( 0 ) * mMasterVolume * mMasterProp( iProp )->dPropdDOF( aDofTypes ) / ( std::pow( mMasterProp( iProp )->val()( 0 ), 2 ) * ( mMasterVolume / mMasterProp( 0 )->val()( 0 ) + mSlaveVolume / mSlaveProp( 0 )->val()( 0 ) ) );
                }
            }
        }

//------------------------------------------------------------------------------
        void SP_Nitsche_Interface::eval_dSPdSlaveDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mSlaveGlobalDofTypeMap( tDofType );

            // reset the matrix
            mdPPdSlaveDof( tDofIndex ).set_size( 1, mSlaveDofFI( tDofIndex )->get_number_of_space_time_coefficients(), 0.0 );

            // if indirect dependency on the dof type
            for( uint iProp = 0; iProp < mSlaveProp.size(); iProp++ )
            {
                if ( mSlaveProp( iProp )->check_dof_dependency( aDofTypes ) )
                {
                    // compute derivative with indirect dependency through properties
                    mdPPdSlaveDof( tDofIndex ).matrix_data()
                    += this->val()( 0 ) * mSlaveVolume * mSlaveProp( iProp )->dPropdDOF( aDofTypes ) / ( std::pow( mSlaveProp( iProp )->val()( 0 ), 2 ) * ( mMasterVolume / mMasterProp( 0 )->val()( 0 ) + mSlaveVolume / mSlaveProp( 0 )->val()( 0 ) ) );
                }
            }
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

