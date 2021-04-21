
#include "cl_FEM_IWG_FS_Struc_Interface.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"

#include "fn_trans.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IWG_FS_Struc_Interface::IWG_FS_Struc_Interface(){}

        //------------------------------------------------------------------------------

        void IWG_FS_Struc_Interface::compute_residual( real aWStar )
        {
#ifdef DEBUG
            // check master and slave field interpolators
            this->check_field_interpolators( mtk::Master_Slave::MASTER );
            this->check_field_interpolators( mtk::Master_Slave::SLAVE );
#endif

            // get slave index for residual dof type, indices for assembly
            uint tSlaveDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::SLAVE );
            uint tSlaveResStartIndex = mSet->get_res_dof_assembly_map()( tSlaveDofIndex )( 0, 0 );
            uint tSlaveResStopIndex  = mSet->get_res_dof_assembly_map()( tSlaveDofIndex )( 0, 1 );

            // get master field interpolator for the residual dof type
            // FIXME protect dof type
            Field_Interpolator * tFIFluidPressure =
                    mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::P );

            // get master field interpolator for the residual dof type
            Field_Interpolator * tFISolidDispl =
                    mSlaveFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // compute slave residual
            mSet->get_residual()( 0 )(
                    { tSlaveResStartIndex, tSlaveResStopIndex },
                    { 0, 0 } ) -= aWStar * (
                            trans( tFISolidDispl->N() ) * tFIFluidPressure->val() );

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_FS_Struc_Interface::compute_residual - Residual contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_FS_Struc_Interface::compute_jacobian( real aWStar )
        {
#ifdef DEBUG
            // check master and slave field interpolators
            this->check_field_interpolators( mtk::Master_Slave::MASTER );
            this->check_field_interpolators( mtk::Master_Slave::SLAVE );
#endif

            // get master index for residual dof type, indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );

            // get slave index for residual dof type, indices for assembly
            uint tSlaveDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::SLAVE );
            uint tSlaveResStartIndex = mSet->get_res_dof_assembly_map()( tSlaveDofIndex )( 0, 0 );
            uint tSlaveResStopIndex  = mSet->get_res_dof_assembly_map()( tSlaveDofIndex )( 0, 1 );

            // get master field interpolator for the residual dof type
            // FIXME protect dof type
            Field_Interpolator * tFIFluidPressure =
                    mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::P );

            // get master field interpolator for the residual dof type
            Field_Interpolator * tFISolidDispl =
                    mSlaveFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get number of master dof dependencies
            uint tMasterNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

            // compute the jacobian for indirect dof dependencies through master constitutive models
            for( uint iDOF = 0; iDOF < tMasterNumDofDependencies; iDOF++ )
            {
                // get the dof type
                Cell< MSI::Dof_Type > & tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                // get the index for the dof type
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tMasterDepStartIndex = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 0 );
                uint tMasterDepStopIndex  = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 1 );

                // compute jacobian direct dependencies
                if ( tDofType( 0 ) == MSI::Dof_Type::P )
                {
                    mSet->get_jacobian()(
                            { tSlaveResStartIndex,  tSlaveResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) -= aWStar * (
                                    trans( tFISolidDispl->N() ) * tFIFluidPressure->N() );
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ) ,
                    "IWG_FS_Struc_Interface::compute_jacobian - Jacobian contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_FS_Struc_Interface::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_FS_Struc_Interface::compute_jacobian_and_residual - This function does nothing.");
        }

        //------------------------------------------------------------------------------

        void IWG_FS_Struc_Interface::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_FS_Struc_Interface::compute_dRdp - This function does nothing.");
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
