//FEM/INT/src
#include "cl_FEM_SP_Dirichlet_Nitsche.hpp"
#include "cl_FEM_Cluster.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
//LINALG/src
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

        void SP_Dirichlet_Nitsche::reset_cluster_measures()
        {
            // evaluate element size from the cluster
            mElementSize = mCluster->compute_cluster_cell_length_measure(
                    mtk::Primary_Void::PRIMARY,
                    mtk::Master_Slave::MASTER );
        }

        //------------------------------------------------------------------------------

        void SP_Dirichlet_Nitsche::set_property(
                std::shared_ptr< Property > aProperty,
                std::string                 aPropertyString,
                mtk::Master_Slave           aIsMaster )
        {
            // check that aPropertyString makes sense
            std::string tErrMsg =
                    std::string( "SP_Dirichlet_Nitsche::set_property - Unknown aPropertyString : " ) +
                    aPropertyString;
            MORIS_ERROR( mPropertyMap.find( aPropertyString ) != mPropertyMap.end() , tErrMsg.c_str() );

            // set the property in the property cell
            this->get_properties( aIsMaster )( static_cast< uint >( mPropertyMap[ aPropertyString ] ) ) = aProperty;
        }

        //------------------------------------------------------------------------------

        void SP_Dirichlet_Nitsche::eval_SP()
        {
            // get the material property
            std::shared_ptr< Property > tPropMaterial =
                    mMasterProp( static_cast< uint >( SP_Property_Type::MATERIAL ) );

            // compute stabilization parameter value
            mPPVal = mParameters( 0 ) * tPropMaterial->val() / mElementSize;
        }

        //------------------------------------------------------------------------------

        void SP_Dirichlet_Nitsche::eval_dSPdMasterDOF(
                const moris::Cell< MSI::Dof_Type > & aDofTypes )
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

            // if material property depends on the dof type
            if ( tPropMaterial->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mdPPdMasterDof( tDofIndex ) +=
                        mParameters( 0 ) * tPropMaterial->dPropdDOF( aDofTypes ) / mElementSize;
            }
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
