
//FEM/INT/src
#include "cl_FEM_SP_Time_Velocity_Ghost.hpp"
#include "cl_FEM_Cluster.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        SP_Time_Velocity_Ghost::SP_Time_Velocity_Ghost()
        {
            // set the property pointer cell size
            mMasterProp.resize( static_cast< uint >( Property_Type::MAX_ENUM ), nullptr );

            // populate the map
            mPropertyMap[ "Density" ] = static_cast< uint >( Property_Type::DENSITY );
        }

        //------------------------------------------------------------------------------

        moris::Cell< std::tuple<
        fem::Measure_Type,
        mtk::Primary_Void,
        mtk::Master_Slave > > SP_Time_Velocity_Ghost::get_cluster_measure_tuple_list()
        {
            return { mElementSizeTuple };
        }

        //------------------------------------------------------------------------------

        void SP_Time_Velocity_Ghost::eval_SP()
        {
            // get element size cluster measure value
            real tElementSize = mCluster->get_cluster_measure(
                    std::get<0>( mElementSizeTuple ),
                    std::get<1>( mElementSizeTuple ),
                    std::get<2>( mElementSizeTuple ) )->val()( 0 );

            // get the density property
            const std::shared_ptr< Property > & tDensityProp =
                    mMasterProp( static_cast< uint >( Property_Type::DENSITY ) );

            // compute time increment deltat
            real tDeltaT = mMasterFIManager->get_IP_geometry_interpolator()->get_time_step();

            // compute stabilization parameter value
            mPPVal =
                    mParameters( 0 ) * std::pow( tElementSize, 2.0 * ( mOrder - 1.0 ) + 3.0 ) *
                    tDensityProp->val()( 0 )  / ( mParameters( 1 )( 0 ) * tDeltaT );
        }

        //------------------------------------------------------------------------------

        void SP_Time_Velocity_Ghost::eval_dSPdMasterDOF(
                const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get element size cluster measure value
            real tElementSize = mCluster->get_cluster_measure(
                    std::get<0>( mElementSizeTuple ),
                    std::get<1>( mElementSizeTuple ),
                    std::get<2>( mElementSizeTuple ) )->val()( 0 );

            // get the dof type index
            uint tDofIndex = mMasterGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the dof type FI
            Field_Interpolator * tFI = mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set size for dSPdMasterDof
            mdPPdMasterDof( tDofIndex ).set_size( 1, tFI->get_number_of_space_time_coefficients() );

            // get the density property
            const std::shared_ptr< Property > & tDensityProp =
                    mMasterProp( static_cast< uint >( Property_Type::DENSITY ) );

            // compute time increment deltat
            real tDeltaT = mMasterFIManager->get_IP_geometry_interpolator()->get_time_step();

            // if density depends on dof type
            if( tDensityProp->check_dof_dependency( aDofTypes ) )
            {
                // compute contribution from density
                mdPPdMasterDof( tDofIndex ) =
                        mParameters( 0 ) * std::pow( tElementSize, 2.0 * ( mOrder - 1.0 ) + 3.0 ) *
                        tDensityProp->dPropdDOF( aDofTypes ) / ( mParameters( 1 )( 0 ) * tDeltaT );
            }
            else
            {
                mdPPdMasterDof( tDofIndex ).fill( 0.0 );
            }
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */


