#include "cl_FEM_SP_Nitsche_Interface.hpp" //FEM/INT/src
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
        SP_Nitsche_Interface::SP_Nitsche_Interface()
        {
            // set size for the property pointer cells
            mMasterProp.resize( static_cast< uint >( SP_Property_Type::MAX_ENUM ), nullptr );
            mSlaveProp.resize( static_cast< uint >( SP_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Material" ] = SP_Property_Type::MATERIAL;
        }

        //------------------------------------------------------------------------------
        void SP_Nitsche_Interface::reset_cluster_measures()
        {
            // evaluate cluster measures from the cluster
            mMasterVolume = mCluster->compute_cluster_cell_measure(
                    mtk::Primary_Void::INTERP,
                    mtk::Master_Slave::MASTER );
            mSlaveVolume = mCluster->compute_cluster_cell_measure(
                    mtk::Primary_Void::INTERP,
                    mtk::Master_Slave::SLAVE );
            mInterfaceSurface = mCluster->compute_cluster_cell_side_measure(
                    mtk::Primary_Void::PRIMARY,
                    mtk::Master_Slave::MASTER );
        }

        //------------------------------------------------------------------------------
        void SP_Nitsche_Interface::set_property(
                std::shared_ptr< Property > aProperty,
                std::string                 aPropertyString,
                mtk::Master_Slave           aIsMaster )
        {
            // check that aPropertyString makes sense
            MORIS_ERROR( mPropertyMap.find( aPropertyString ) != mPropertyMap.end(),
                    "SP_Nitsche_Interface::set_property - Unknown aPropertyString." );

            // set the property in the property cell
            this->get_properties( aIsMaster )( static_cast< uint >( mPropertyMap[ aPropertyString ] ) ) = aProperty;
        }

        //------------------------------------------------------------------------------
        void SP_Nitsche_Interface::eval_SP()
        {
            // get the master/slave material property
            std::shared_ptr< Property > tMasterPropMaterial =
                    mMasterProp( static_cast< uint >( SP_Property_Type::MATERIAL ) );
            std::shared_ptr< Property > tSlavePropMaterial =
                    mSlaveProp( static_cast< uint >( SP_Property_Type::MATERIAL ) );

            // compute stabilization parameter value
            mPPVal = 2.0 * mParameters( 0 ) * mInterfaceSurface /
                    ( mMasterVolume / tMasterPropMaterial->val()( 0 ) + mSlaveVolume / tSlavePropMaterial->val()( 0 ) );
        }

        //------------------------------------------------------------------------------
        void SP_Nitsche_Interface::eval_dSPdMasterDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mMasterGlobalDofTypeMap( tDofType );

            // get the derivative dof type FI
            Field_Interpolator * tFIDer =
                    mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // reset the matrix
            mdPPdMasterDof( tDofIndex ).set_size(
                    1,
                    tFIDer->get_number_of_space_time_coefficients(),
                    0.0 );

            // get the master/slave material property
            std::shared_ptr< Property > tMasterPropMaterial =
                    mMasterProp( static_cast< uint >( SP_Property_Type::MATERIAL ) );
            std::shared_ptr< Property > tSlavePropMaterial =
                    mSlaveProp( static_cast< uint >( SP_Property_Type::MATERIAL ) );

            // if indirect dependency on the dof type
            if ( tMasterPropMaterial->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mdPPdMasterDof( tDofIndex ).matrix_data() +=
                        this->val()( 0 ) * mMasterVolume * tMasterPropMaterial->dPropdDOF( aDofTypes ) /
                        ( std::pow( tMasterPropMaterial->val()( 0 ), 2 ) * ( mMasterVolume / tMasterPropMaterial->val()( 0 ) + mSlaveVolume / tSlavePropMaterial->val()( 0 ) ) );
            }
        }

        //------------------------------------------------------------------------------
        void SP_Nitsche_Interface::eval_dSPdSlaveDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mSlaveGlobalDofTypeMap( tDofType );

            // get the derivative dof type FI
            Field_Interpolator * tFIDer =
                    mSlaveFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // reset the matrix
            mdPPdSlaveDof( tDofIndex ).set_size(
                    1,
                    tFIDer->get_number_of_space_time_coefficients(),
                    0.0 );

            // get the master/slave material property
            std::shared_ptr< Property > tMasterPropMaterial =
                    mMasterProp( static_cast< uint >( SP_Property_Type::MATERIAL ) );
            std::shared_ptr< Property > tSlavePropMaterial =
                    mSlaveProp( static_cast< uint >( SP_Property_Type::MATERIAL ) );

            // if indirect dependency on the dof type
            if ( tSlavePropMaterial->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mdPPdSlaveDof( tDofIndex ).matrix_data() +=
                        this->val()( 0 ) * mSlaveVolume * tSlavePropMaterial->dPropdDOF( aDofTypes ) /
                        ( std::pow( tSlavePropMaterial->val()( 0 ), 2 ) * ( mMasterVolume / tMasterPropMaterial->val()( 0 ) + mSlaveVolume / tSlavePropMaterial->val()( 0 ) ) );
            }
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

