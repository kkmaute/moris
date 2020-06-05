//FEM/INT/src
#include "cl_FEM_SP_Ghost_Turbulence.hpp"
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
        SP_Ghost_Turbulence::SP_Ghost_Turbulence()
        {
            // set size for the property pointer cells
            mMasterProp.resize( static_cast< uint >( SP_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Viscosity" ] = SP_Property_Type::VISCOSITY;
        }

        //------------------------------------------------------------------------------
        void SP_Ghost_Turbulence::reset_cluster_measures()
        {
            // evaluate element size from the cluster
            mElementSize = mCluster->compute_cluster_cell_length_measure(
                    mtk::Primary_Void::PRIMARY,
                    mtk::Master_Slave::MASTER );
        }

        //------------------------------------------------------------------------------
        void SP_Ghost_Turbulence::set_dof_type_list(
                moris::Cell< moris::Cell< MSI::Dof_Type > > & aDofTypes,
                moris::Cell< std::string >                  & aDofStrings,
                mtk::Master_Slave                             aIsMaster )
        {
            // switch on master slave
            switch ( aIsMaster )
            {
                case mtk::Master_Slave::MASTER :
                {
                    // set dof type list
                    mMasterDofTypes = aDofTypes;

                    // loop on dof type
                    for( uint iDof = 0; iDof < aDofTypes.size(); iDof++ )
                    {
                        // get dof string
                        std::string tDofString = aDofStrings( iDof );

                        // get dof type
                        MSI::Dof_Type tDofType = aDofTypes( iDof )( 0 );

                        // if velocity
                        if( tDofString == "Viscosity" )
                        {
                            mMasterDofViscosity = tDofType;
                        }
                        else
                        {
                            // create error message
                            std::string tErrMsg =
                                    std::string( "SP_Ghost_Turbulence::set_dof_type_list - Unknown aDofString : ") +
                                    tDofString;
                            MORIS_ERROR( false , tErrMsg.c_str() );
                        }
                    }
                    break;
                }

                case mtk::Master_Slave::SLAVE :
                {
                    // set dof type list
                    mSlaveDofTypes = aDofTypes;
                    break;
                }

                default:
                    MORIS_ERROR( false, "SP_Ghost_Turbulence::set_dof_type_list - unknown master slave type." );
                    break;
            }
        }


        //------------------------------------------------------------------------------
        void SP_Ghost_Turbulence::set_property(
                std::shared_ptr< Property > aProperty,
                std::string                 aPropertyString,
                mtk::Master_Slave           aIsMaster )
        {
            // check that aPropertyString makes sense
            std::string tErrMsg =
                    std::string( "SP_Ghost_Turbulence::set_property - Unknown aPropertyString : " ) +
                    aPropertyString;
            MORIS_ERROR( mPropertyMap.find( aPropertyString ) != mPropertyMap.end() , tErrMsg.c_str() );

            // set the property in the property cell
            this->get_properties( aIsMaster )( static_cast< uint >( mPropertyMap[ aPropertyString ] ) ) = aProperty;
        }

        //------------------------------------------------------------------------------
        void SP_Ghost_Turbulence::eval_SP()
        {
            // get the velocity FI
            Field_Interpolator * tFIViscosity =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofViscosity );

            // get the material property
            std::shared_ptr< Property > tPropViscosity =
                    mMasterProp( static_cast< uint >( SP_Property_Type::VISCOSITY ) );

            // compute stabilization parameter value
            mPPVal = mParameters( 0 ) *
                    std::pow( mElementSize, 2 * ( mOrder - 1 ) + 1 ) *
                    ( tFIViscosity->val() + tPropViscosity->val() ) / mSigma;
        }

        //------------------------------------------------------------------------------
        void SP_Ghost_Turbulence::eval_dSPdMasterDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
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
            std::shared_ptr< Property > tPropViscosity =
                    mMasterProp( static_cast< uint >( SP_Property_Type::VISCOSITY ) );

            // if derivative dof type is viscosity dof type
            if( aDofTypes( 0 ) == mMasterDofViscosity )
            {
                mdPPdMasterDof( tDofIndex ).matrix_data() +=
                        mParameters( 0 ) * std::pow( mElementSize, 2 * ( mOrder - 1 ) + 1 ) *
                        tFIDer->N() /  mSigma;
            }

            // if material property depends on the dof type
            if ( tPropViscosity->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mdPPdMasterDof( tDofIndex ).matrix_data() +=
                        mParameters( 0 ) * std::pow( mElementSize, 2 * ( mOrder - 1 ) + 1 ) *
                        tPropViscosity->dPropdDOF( aDofTypes ) / mSigma;
            }
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */


