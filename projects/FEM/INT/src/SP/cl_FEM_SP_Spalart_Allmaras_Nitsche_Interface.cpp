//FEM/INT/src
#include "cl_FEM_SP_Spalart_Allmaras_Nitsche_Interface.hpp"
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

        SP_Spalart_Allmaras_Nitsche_Interface::SP_Spalart_Allmaras_Nitsche_Interface()
        {
            // set has slave flag to true
            mHasSlave = true;

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );
            mSlaveCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "SpalartAllmarasTurbulence" ] = static_cast< uint >( IWG_Constitutive_Type::SPALART_ALLMARAS_TURBULENCE );
        }

        //------------------------------------------------------------------------------

        void SP_Spalart_Allmaras_Nitsche_Interface::set_dof_type_list(
                moris::Cell< moris::Cell< MSI::Dof_Type > > & aDofTypes,
                moris::Cell< std::string >                  & aDofStrings,
                mtk::Master_Slave                             aIsMaster )
        {
            Stabilization_Parameter::set_dof_type_list( aDofTypes, aIsMaster );
        }

        //------------------------------------------------------------------------------

        moris::Cell< std::tuple<
        fem::Measure_Type,
        mtk::Primary_Void,
        mtk::Master_Slave > > SP_Spalart_Allmaras_Nitsche_Interface::get_cluster_measure_tuple_list()
        {
            return { mMasterVolumeTuple, mSlaveVolumeTuple, mInterfaceSurfaceTuple };
        }

        //------------------------------------------------------------------------------

        void SP_Spalart_Allmaras_Nitsche_Interface::eval_SP()
        {
            // get master volume cluster measure value
            real tMasterVolume = mCluster->get_cluster_measure(
                    std::get<0>( mMasterVolumeTuple ),
                    std::get<1>( mMasterVolumeTuple ),
                    std::get<2>( mMasterVolumeTuple ) )->val()( 0 );

            // get slave volume cluster measure value
            real tSlaveVolume = mCluster->get_cluster_measure(
                    std::get<0>( mSlaveVolumeTuple ),
                    std::get<1>( mSlaveVolumeTuple ),
                    std::get<2>( mSlaveVolumeTuple ) )->val()( 0 );

            // get interface surface cluster measure value
            real tInterfaceSurface = mCluster->get_cluster_measure(
                    std::get<0>( mInterfaceSurfaceTuple ),
                    std::get<1>( mInterfaceSurfaceTuple ),
                    std::get<2>( mInterfaceSurfaceTuple ) )->val()( 0 );

            // init SP values
            mPPVal.set_size( 3, 1, 0.0 );

            // get the SA turbulence CM
            const std::shared_ptr< Constitutive_Model > & tCMMasterSATurbulence =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::SPALART_ALLMARAS_TURBULENCE ) );
            const std::shared_ptr< Constitutive_Model > & tCMSlaveSATurbulence =
                    mSlaveCM( static_cast< uint >( IWG_Constitutive_Type::SPALART_ALLMARAS_TURBULENCE ) );

            // compute weighted property
            real tDeno =
                    tMasterVolume / tCMMasterSATurbulence->diffusion_coefficient()( 0 ) +
                    tSlaveVolume / tCMSlaveSATurbulence->diffusion_coefficient()( 0 );

            // compute Nitsche parameter value
            mPPVal( 0 ) = 2.0 * mParameters( 0 )( 0 ) * tInterfaceSurface / tDeno;

            // compute master weight value
            mPPVal( 1 ) = ( tMasterVolume / tCMMasterSATurbulence->diffusion_coefficient()( 0 ) ) / tDeno;

            // compute slave weight value
            mPPVal( 2 ) = ( tSlaveVolume / tCMSlaveSATurbulence->diffusion_coefficient()( 0 ) ) / tDeno;
        }

        //------------------------------------------------------------------------------

        void SP_Spalart_Allmaras_Nitsche_Interface::eval_dSPdMasterDOF(
                const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get master volume cluster measure value
            real tMasterVolume = mCluster->get_cluster_measure(
                    std::get<0>( mMasterVolumeTuple ),
                    std::get<1>( mMasterVolumeTuple ),
                    std::get<2>( mMasterVolumeTuple ) )->val()( 0 );

            // get slave volume cluster measure value
            real tSlaveVolume = mCluster->get_cluster_measure(
                    std::get<0>( mSlaveVolumeTuple ),
                    std::get<1>( mSlaveVolumeTuple ),
                    std::get<2>( mSlaveVolumeTuple ) )->val()( 0 );

            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mMasterGlobalDofTypeMap( tDofType );

            // get the derivative dof type FI
            Field_Interpolator * tFIDer =
                    mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // reset the matrix
            mdPPdMasterDof( tDofIndex ).set_size(
                    3,
                    tFIDer->get_number_of_space_time_coefficients() );

            // get the SA turbulence CM
            const std::shared_ptr< Constitutive_Model > & tCMMasterSATurbulence =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::SPALART_ALLMARAS_TURBULENCE ) );
            const std::shared_ptr< Constitutive_Model > & tCMSlaveSATurbulence =
                    mSlaveCM( static_cast< uint >( IWG_Constitutive_Type::SPALART_ALLMARAS_TURBULENCE ) );

            // if master turbulence CM depends on dof type
            if ( tCMMasterSATurbulence->check_dof_dependency( aDofTypes ) )
            {
                // compute weighted property
                real tDeno =
                        tMasterVolume / tCMMasterSATurbulence->diffusion_coefficient()( 0 ) +
                        tSlaveVolume / tCMSlaveSATurbulence->diffusion_coefficient()( 0 );

                // compute derivative of Nitsche SP FIXME ????
                mdPPdMasterDof( tDofIndex ).get_row( 0 ) =
                        this->val()( 0 ) * tMasterVolume * tCMMasterSATurbulence->ddiffusioncoeffdu( aDofTypes ) /
                        ( std::pow( tCMMasterSATurbulence->diffusion_coefficient()( 0 ), 2 ) * tDeno );

                // compute derivative of master weight SP
                mdPPdMasterDof( tDofIndex ).get_row( 1 ) =
                        - this->val()( 1 ) * tSlaveVolume * tCMMasterSATurbulence->ddiffusioncoeffdu( aDofTypes ) /
                        ( tCMMasterSATurbulence->diffusion_coefficient()( 0 ) * tCMSlaveSATurbulence->diffusion_coefficient()( 0 ) * tDeno );

                // compute derivative of the slave weight SP
                mdPPdMasterDof( tDofIndex ).get_row( 2 ) =
                        this->val()( 2 ) * tMasterVolume * tCMMasterSATurbulence->ddiffusioncoeffdu( aDofTypes ) /
                        ( std::pow( tCMMasterSATurbulence->diffusion_coefficient()( 0 ), 2 ) * tDeno );
            }
            else
            {
                mdPPdMasterDof( tDofIndex ).fill( 0.0 );
            }
        }

        //------------------------------------------------------------------------------

        void SP_Spalart_Allmaras_Nitsche_Interface::eval_dSPdSlaveDOF(
                const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get master volume cluster measure value
            real tMasterVolume = mCluster->get_cluster_measure(
                    std::get<0>( mMasterVolumeTuple ),
                    std::get<1>( mMasterVolumeTuple ),
                    std::get<2>( mMasterVolumeTuple ) )->val()( 0 );

            // get slave volume cluster measure value
            real tSlaveVolume = mCluster->get_cluster_measure(
                    std::get<0>( mSlaveVolumeTuple ),
                    std::get<1>( mSlaveVolumeTuple ),
                    std::get<2>( mSlaveVolumeTuple ) )->val()( 0 );

            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mSlaveGlobalDofTypeMap( tDofType );

            // get the derivative dof type FI
            Field_Interpolator * tFIDer =
                    mSlaveFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // reset the matrix
            mdPPdSlaveDof( tDofIndex ).set_size( 3, tFIDer->get_number_of_space_time_coefficients(), 0.0 );

            // get the SA turbulence CM
            const std::shared_ptr< Constitutive_Model > & tCMMasterSATurbulence =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::SPALART_ALLMARAS_TURBULENCE ) );
            const std::shared_ptr< Constitutive_Model > & tCMSlaveSATurbulence =
                    mSlaveCM( static_cast< uint >( IWG_Constitutive_Type::SPALART_ALLMARAS_TURBULENCE ) );

            // if slave turbulence CM depends on dof type
            if ( tCMSlaveSATurbulence->check_dof_dependency( aDofTypes ) )
            {
                // compute weighted property
                real tDeno =
                        tMasterVolume / tCMMasterSATurbulence->diffusion_coefficient()( 0 ) +
                        tSlaveVolume / tCMSlaveSATurbulence->diffusion_coefficient()( 0 );

                // compute derivative of Nitsche SP FIXME ????
                mdPPdSlaveDof( tDofIndex ).get_row( 0 ) =
                        this->val()( 0 ) * tSlaveVolume * tCMSlaveSATurbulence->ddiffusioncoeffdu( aDofTypes ) /
                        ( std::pow( tCMSlaveSATurbulence->diffusion_coefficient()( 0 ), 2 ) * tDeno );

                // compute derivative of master weight SP
                mdPPdSlaveDof( tDofIndex ).get_row( 1 ) =
                        this->val()( 1 ) * tSlaveVolume * tCMSlaveSATurbulence->ddiffusioncoeffdu( aDofTypes ) /
                        ( std::pow( tCMSlaveSATurbulence->diffusion_coefficient()( 0 ), 2 ) * tDeno );

                // compute derivative of the slave weight SP
                mdPPdSlaveDof( tDofIndex ).get_row( 2 ) =
                        - this->val()( 2 ) * tMasterVolume * tCMSlaveSATurbulence->ddiffusioncoeffdu( aDofTypes ) /
                        ( tCMMasterSATurbulence->diffusion_coefficient()( 0 ) * tCMSlaveSATurbulence->diffusion_coefficient()( 0 ) * tDeno );
            }
            else
            {
                mdPPdSlaveDof( tDofIndex ).fill( 0.0 );
            }
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

