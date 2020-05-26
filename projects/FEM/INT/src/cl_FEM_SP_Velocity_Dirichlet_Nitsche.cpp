
//FEM/INT/src
#include "cl_FEM_SP_Velocity_Dirichlet_Nitsche.hpp"
#include "cl_FEM_Cluster.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------
        SP_Velocity_Dirichlet_Nitsche::SP_Velocity_Dirichlet_Nitsche()
        {
            // set the property pointer cell size
            mMasterProp.resize( static_cast< uint >( Property_Type::MAX_ENUM ), nullptr );

            // populate the map
            mPropertyMap[ "Viscosity" ] = Property_Type::VISCOSITY;
            mPropertyMap[ "Density" ]   = Property_Type::DENSITY;
        }

        //------------------------------------------------------------------------------
        void SP_Velocity_Dirichlet_Nitsche::reset_cluster_measures()
        {
            // evaluate element size from the cluster
            mElementSize = mCluster->compute_cluster_cell_length_measure(
                    mtk::Primary_Void::PRIMARY,
                    mtk::Master_Slave::MASTER );
        }

        //------------------------------------------------------------------------------
        void SP_Velocity_Dirichlet_Nitsche::set_dof_type_list(
                moris::Cell< moris::Cell< MSI::Dof_Type > > & aDofTypes,
                moris::Cell< std::string >                  & aDofStrings,
                mtk::Master_Slave                             aIsMaster )
        {
            // set dof type list
            mMasterDofTypes = aDofTypes;

            // switch on master slave
            switch ( aIsMaster )
            {
                case mtk::Master_Slave::MASTER :
                {
                    // loop on dof type
                    for( uint iDof = 0; iDof < aDofTypes.size(); iDof++ )
                    {
                        // get dof string
                        std::string tDofString = aDofStrings( iDof );

                        // get dof type
                        MSI::Dof_Type tDofType = aDofTypes( iDof )( 0 );

                        // if viscosity
                        if( "Velocity" )
                        {
                            mMasterDofVelocity = tDofType;
                        }
                        else
                        {
                            // create error message
                            std::string tErrMsg =
                                    std::string( "SP_Velocity_Dirichlet_Nitsche::set_dof_type_list - Unknown aDofString : ") +
                                    tDofString;
                            MORIS_ERROR( false , tErrMsg.c_str() );
                        }
                    }
                    break;
                }

                default:
                    MORIS_ERROR( false, "SP_Velocity_Dirichlet_Nitsche::set_dof_type_list - unknown or incorrect master slave type." );
                    break;
            }
        }

        //------------------------------------------------------------------------------
        void SP_Velocity_Dirichlet_Nitsche::set_property(
                std::shared_ptr< Property > aProperty,
                std::string                 aPropertyString,
                mtk::Master_Slave           aIsMaster )
        {
            // check that aPropertyString makes sense
            MORIS_ERROR( mPropertyMap.find( aPropertyString ) != mPropertyMap.end(),
                    "SP_Velocity_Dirichlet_Nitsche::set_property - Unknown aPropertyString." );

            // set the property in the property cell
            this->get_properties( aIsMaster )( static_cast< uint >( mPropertyMap[ aPropertyString ] ) ) = aProperty;
        }

        //------------------------------------------------------------------------------
        void SP_Velocity_Dirichlet_Nitsche::eval_SP()
        {
            // get the viscosity and density property
            std::shared_ptr< Property > tPropViscosity = mMasterProp( static_cast< uint >( Property_Type::VISCOSITY ) );
            std::shared_ptr< Property > tPropDensity   = mMasterProp( static_cast< uint >( Property_Type::DENSITY ) );

            // get the velocity FI
            Field_Interpolator * tVelocityFI = mMasterFIManager->get_field_interpolators_for_type( mMasterDofVelocity );

            // compute infinity norm of u
            real tInfinityNorm = std::abs( tVelocityFI->val()( 0 ) );
            for( uint iDim = 0; iDim < tVelocityFI->val().numel(); iDim++ )
            {
                real tAbsVelocity = std::abs( tVelocityFI->val()( iDim ) );
                if ( tInfinityNorm < tAbsVelocity )
                {
                    tInfinityNorm = tAbsVelocity;
                }
            }

            // compute deltaT
            real tDeltaT = mMasterFIManager->get_IP_geometry_interpolator()->get_time_step();

            // compute stabilization parameter value
            mPPVal = mParameters( 0 ) * (
                    tPropViscosity->val()( 0 ) / mElementSize +
                    tPropDensity->val()( 0 ) * tInfinityNorm / 6.0 +
                    tPropDensity->val()( 0 ) * mElementSize / ( 12.0 * mParameters( 1 )( 0 ) * tDeltaT ) );
        }

        //------------------------------------------------------------------------------
        void SP_Velocity_Dirichlet_Nitsche::eval_dSPdMasterDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type index
            uint tDofIndex = mMasterGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the dof type FI
            Field_Interpolator * tFI = mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set size for dSPdMasterDof
            mdPPdMasterDof( tDofIndex ).set_size( 1, tFI->get_number_of_space_time_coefficients(), 0.0 );

            // get the viscosity and density property
            std::shared_ptr< Property > tPropViscosity = mMasterProp( static_cast< uint >( Property_Type::VISCOSITY ) );
            std::shared_ptr< Property > tPropDensity   = mMasterProp( static_cast< uint >( Property_Type::DENSITY ) );

            // compute infinity norm
            Field_Interpolator * tVelocityFI = mMasterFIManager->get_field_interpolators_for_type( mMasterDofVelocity );
            uint tInfinityNormIndex = 0;
            real tInfinityNorm = std::abs( tVelocityFI->val()( 0 ) );
            for( uint iDim = 0; iDim < tVelocityFI->val().numel(); iDim++ )
            {
                real tAbsVelocity = std::abs( tVelocityFI->val()( iDim ) );
                if ( tInfinityNorm < tAbsVelocity )
                {
                    tInfinityNormIndex = iDim;
                    tInfinityNorm = tAbsVelocity;
                }
            }

            // FIXME infinity norm
            // if dof type == velocity
            if( aDofTypes( 0 ) == mMasterDofVelocity )
            {
                Matrix< DDRMat > tdInfinityNormdu;

                // compute derivative of infinity norm
                if( tInfinityNorm == 0.0 )
                {
                    tdInfinityNormdu = tVelocityFI->N().get_row( tInfinityNormIndex );
                }
                else
                {
                    tdInfinityNormdu = tVelocityFI->N().get_row( tInfinityNormIndex ) * tVelocityFI->val()( tInfinityNormIndex ) / tInfinityNorm;
                }

                // compute contribution from velocity
                mdPPdMasterDof( tDofIndex ).matrix_data()
                        += mParameters( 0 )( 0 ) * tPropDensity->val()( 0 ) * tdInfinityNormdu / 6.0;
            }

            // if viscosity depends on dof type
            if( tPropViscosity->check_dof_dependency( aDofTypes ) )
            {
                // compute contribution from viscosity
                mdPPdMasterDof( tDofIndex ).matrix_data()
                        += mParameters( 0 )( 0 ) * tPropViscosity->dPropdDOF( aDofTypes ) / mElementSize;
            }

            // if density depends on dof type
            if( tPropDensity->check_dof_dependency( aDofTypes ) )
            {
                // compute deltaT
                real tDeltaT = mMasterFIManager->get_IP_geometry_interpolator()->get_time_step();

                // compute contribution from density
                mdPPdMasterDof( tDofIndex ).matrix_data()
                        += mParameters( 0 )( 0 ) * tPropDensity->dPropdDOF( aDofTypes )
                        * ( tInfinityNorm / 6.0 + mElementSize / ( 12.0 * mParameters( 1 )( 0 ) * tDeltaT ) );
            }
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */


