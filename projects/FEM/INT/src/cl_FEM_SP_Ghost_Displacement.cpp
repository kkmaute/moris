//FEM/INT/src
#include "cl_FEM_SP_Ghost_Displacement.hpp"
#include "cl_FEM_Cluster.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

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
        SP_Ghost_Displacement::SP_Ghost_Displacement()
        {
            // set size for the property pointer cells
            mMasterProp.resize( static_cast< uint >( SP_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Material" ] = SP_Property_Type::MATERIAL;
        }

        //------------------------------------------------------------------------------
        void SP_Ghost_Displacement::reset_cluster_measures()
        {
            // evaluate element size from the cluster
            mElementSize = mCluster->compute_cluster_cell_length_measure(
                    mtk::Primary_Void::PRIMARY,
                    mtk::Master_Slave::MASTER );
        }

        //------------------------------------------------------------------------------
        void SP_Ghost_Displacement::set_property(
                std::shared_ptr< Property > aProperty,
                std::string                 aPropertyString,
                mtk::Master_Slave           aIsMaster )
        {
            // check that aPropertyString makes sense
            std::string tErrMsg =
                    std::string( "SP_Ghost_Displacement::set_property - Unknown aPropertyString : " ) +
                    aPropertyString;
            MORIS_ERROR( mPropertyMap.find( aPropertyString ) != mPropertyMap.end() , tErrMsg.c_str() );

            // set the property in the property cell
            this->get_properties( aIsMaster )( static_cast< uint >( mPropertyMap[ aPropertyString ] ) ) = aProperty;
        }

        //------------------------------------------------------------------------------
        void SP_Ghost_Displacement::eval_SP()
        {
            // compute stabilization parameter value
            mPPVal = mParameters( 0 ) *
                    std::pow( mElementSize, 2 * ( mOrder - 1 ) + 1 ) *
                    mMasterProp( static_cast< uint >( SP_Property_Type::MATERIAL ) )->val()( 0 );
        }

        //------------------------------------------------------------------------------
        void SP_Ghost_Displacement::eval_dSPdMasterDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mMasterGlobalDofTypeMap( tDofType );

            // get FI for derivative dof type
            Field_Interpolator * tFIDer =
                    mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // reset the matrix
            mdPPdMasterDof( tDofIndex ).set_size(
                    1,
                    tFIDer->get_number_of_space_time_coefficients(),
                    0.0 );

            // get the material property
            std::shared_ptr< Property > tPropMaterial =
                    mMasterProp( static_cast< uint >( SP_Property_Type::MATERIAL ) );

            if ( tPropMaterial->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mdPPdMasterDof( tDofIndex ).matrix_data() +=
                        mParameters( 0 ) *
                        std::pow( mElementSize, 2 * ( mOrder - 1 ) + 1 ) *
                        tPropMaterial->dPropdDOF( aDofTypes );
            }
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */


