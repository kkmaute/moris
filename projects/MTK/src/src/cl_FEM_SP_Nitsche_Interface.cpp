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
        SP_Nitsche_Interface::SP_Nitsche_Interface()
        {
            // set size for the property pointer cells
            mMasterProp.resize( static_cast< uint >( SP_Property_Type::MAX_ENUM ), nullptr );
            mSlaveProp.resize( static_cast< uint >( SP_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Material" ] = SP_Property_Type::MATERIAL;

            // set the list of cluster measures
            mClusterMeasures = { fem::Cluster_Measure::INTERFACE_SURFACE,
                                 fem::Cluster_Measure::MASTER_VOLUME,
                                 fem::Cluster_Measure::SLAVE_VOLUME };
        }

//------------------------------------------------------------------------------
        void SP_Nitsche_Interface::eval_SP()
        {
            // compute stabilization parameter value
            mPPVal = 2.0 * mParameters( 0 ) * mInterfaceSurface / ( mMasterVolume / mMasterProp( static_cast< uint >( SP_Property_Type::MATERIAL ) )->val()( 0 ) + mSlaveVolume / mSlaveProp( static_cast< uint >( SP_Property_Type::MATERIAL ) )->val()( 0 ) );
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
            if ( mMasterProp( static_cast< uint >( SP_Property_Type::MATERIAL ) )->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mdPPdMasterDof( tDofIndex ).matrix_data()
                += this->val()( 0 ) * mMasterVolume * mMasterProp( static_cast< uint >( SP_Property_Type::MATERIAL ) )->dPropdDOF( aDofTypes ) / ( std::pow( mMasterProp( static_cast< uint >( SP_Property_Type::MATERIAL ) )->val()( 0 ), 2 ) * ( mMasterVolume / mMasterProp( static_cast< uint >( SP_Property_Type::MATERIAL ) )->val()( 0 ) + mSlaveVolume / mSlaveProp( static_cast< uint >( SP_Property_Type::MATERIAL ) )->val()( 0 ) ) );
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
            if ( mSlaveProp( static_cast< uint >( SP_Property_Type::MATERIAL ) )->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mdPPdSlaveDof( tDofIndex ).matrix_data()
                += this->val()( 0 ) * mSlaveVolume * mSlaveProp( static_cast< uint >( SP_Property_Type::MATERIAL ) )->dPropdDOF( aDofTypes ) / ( std::pow( mSlaveProp( static_cast< uint >( SP_Property_Type::MATERIAL ) )->val()( 0 ), 2 ) * ( mMasterVolume / mMasterProp( static_cast< uint >( SP_Property_Type::MATERIAL ) )->val()( 0 ) + mSlaveVolume / mSlaveProp( static_cast< uint >( SP_Property_Type::MATERIAL ) )->val()( 0 ) ) );
            }
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

