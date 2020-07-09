//FEM/INT/src
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG_Incompressible_NS_Convective_Normal_Velocity_Ghost.hpp"
//LINALG/src
#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        IWG_Incompressible_NS_Convective_Normal_Velocity_Ghost::IWG_Incompressible_NS_Convective_Normal_Velocity_Ghost()
        {
            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "ConvectiveGhost" ] = IWG_Stabilization_Type::CONVECTIVE_GHOST;
        }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Convective_Normal_Velocity_Ghost::set_stabilization_parameter(
                std::shared_ptr< Stabilization_Parameter > aStabilizationParameter,
                std::string                                aStabilizationString )
        {
            // check that aStabilizationString makes sense
            std::string tErrMsg =
                    std::string( "IWG_Incompressible_NS_Convective_Normal_Velocity_Ghost::set_stabilization_parameter - Unknown aStabilizationString: " ) +
                    aStabilizationString;
            MORIS_ERROR( mStabilizationMap.find( aStabilizationString ) != mStabilizationMap.end(), tErrMsg.c_str() );

            // set the stabilization parameter in the stabilization parameter cell
            this->get_stabilization_parameters()( static_cast< uint >( mStabilizationMap[ aStabilizationString ] ) ) = aStabilizationParameter;
        }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Convective_Normal_Velocity_Ghost::compute_residual( real aWStar )
        {
            // check master field interpolators
#ifdef DEBUG
            this->check_field_interpolators();
#endif

            // get master index for residual dof type, indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get slave index for residual dof type, indices for assembly
            uint tSlaveDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::SLAVE );
            uint tSlaveResStartIndex = mSet->get_res_dof_assembly_map()( tSlaveDofIndex )( 0, 0 );
            uint tSlaveResStopIndex  = mSet->get_res_dof_assembly_map()( tSlaveDofIndex )( 0, 1 );

            // get the master field interpolator for the residual dof type
            Field_Interpolator * tFIMaster = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get the slave field interpolator for the residual dof type
            Field_Interpolator * tFISlave  = mSlaveFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get the convective stabilization parameter
            std::shared_ptr< Stabilization_Parameter > tSPConvective =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::CONVECTIVE_GHOST ) );

            // get flattened directional derivatives for master and slave
            Matrix< DDRMat > tMasterdNdxFlat;
            this->compute_flat_dnNdxn( tMasterdNdxFlat, mtk::Master_Slave::MASTER );
            Matrix< DDRMat > tSlavedNdxFlat;
            this->compute_flat_dnNdxn( tSlavedNdxFlat, mtk::Master_Slave::SLAVE );

            // premultiply common terms
            Matrix< DDRMat > tConvectivePreMultiply =
                    tSPConvective->val()( 0 ) * trans( mNormal ) *
                    ( tFIMaster->gradx( 1 ) - tFISlave->gradx( 1 ) ) ;
            tConvectivePreMultiply = reshape( tConvectivePreMultiply, tConvectivePreMultiply.numel(), 1 );

            // compute master residual
            mSet->get_residual()( 0 )( { tMasterResStartIndex, tMasterResStopIndex }, { 0, 0 } ) += aWStar * (
                    tMasterdNdxFlat * tConvectivePreMultiply );

            // compute slave residual
            mSet->get_residual()( 0 )( { tSlaveResStartIndex, tSlaveResStopIndex }, { 0, 0 } ) -= aWStar * (
                    tSlavedNdxFlat * tConvectivePreMultiply );
        }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Convective_Normal_Velocity_Ghost::compute_jacobian( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            // get master index for residual dof type, indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get slave index for residual dof type, indices for assembly
            uint tSlaveDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::SLAVE );
            uint tSlaveResStartIndex = mSet->get_res_dof_assembly_map()( tSlaveDofIndex )( 0, 0 );
            uint tSlaveResStopIndex  = mSet->get_res_dof_assembly_map()( tSlaveDofIndex )( 0, 1 );

            // get the master field interpolator for residual dof type
            Field_Interpolator * tFIMaster =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get the slave field interpolator for residual dof type
            Field_Interpolator * tFISlave  =
                    mSlaveFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get the convective stabilization parameter
            std::shared_ptr< Stabilization_Parameter > tSPConvective =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::CONVECTIVE_GHOST ) );

            // get flattened directional derivatives for master and slave
            Matrix< DDRMat > tMasterdNdxFlat;
            this->compute_flat_dnNdxn( tMasterdNdxFlat, mtk::Master_Slave::MASTER );
            Matrix< DDRMat > tSlavedNdxFlat;
            this->compute_flat_dnNdxn( tSlavedNdxFlat, mtk::Master_Slave::SLAVE );

            // get number of master dependencies
            uint tMasterNumDofDependencies = mRequestedMasterGlobalDofTypes.size();
            uint tSlaveNumDofDependencies  = mRequestedSlaveGlobalDofTypes.size();

            // compute the jacobian for indirect dof dependencies through master
            for( uint iDOF = 0; iDOF < tMasterNumDofDependencies; iDOF++ )
            {
                // get the dof type
                Cell< MSI::Dof_Type > tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                // get the index for the dof type
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tMasterDepStartIndex = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 0 );
                uint tMasterDepStopIndex  = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 1 );

                // compute jacobian direct dependencies
                if ( tDofType( 0 ) == mResidualDofType( 0 ) )
                {
                    // dRM/dM
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    tMasterdNdxFlat * tSPConvective->val()( 0 ) * trans( tMasterdNdxFlat ) );

                    // dRS/dM
                    mSet->get_jacobian()(
                            { tSlaveResStartIndex,  tSlaveResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) -= aWStar * (
                                    tSlavedNdxFlat * tSPConvective->val()( 0 ) * trans( tMasterdNdxFlat ) );
                }

                // if convective stabilization parameter depends on dof type
                if ( tSPConvective->check_dof_dependency( tDofType ) )
                {
                    // premultiply common terms
                    Matrix< DDRMat > tPreMultiply = trans( mNormal ) * ( tFIMaster->gradx( 1 ) - tFISlave->gradx( 1 ) );
                    tPreMultiply = reshape( tPreMultiply , tPreMultiply.numel(), 1 );
                    tPreMultiply = tPreMultiply * tSPConvective->dSPdMasterDOF( tDofType );

                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    tMasterdNdxFlat * tPreMultiply );

                    mSet->get_jacobian()(
                            { tSlaveResStartIndex,  tSlaveResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) -= aWStar * (
                                    tSlavedNdxFlat  * tPreMultiply );
                }
            }

            // compute the jacobian for indirect dof dependencies through slave
            for( uint iDOF = 0; iDOF < tSlaveNumDofDependencies; iDOF++ )
            {
                // get the dof type
                Cell< MSI::Dof_Type > tDofType = mRequestedSlaveGlobalDofTypes( iDOF );

                // get the index for the dof type
                sint tDofDepIndex        = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::SLAVE );
                uint tSlaveDepStartIndex = mSet->get_jac_dof_assembly_map()( tSlaveDofIndex )( tDofDepIndex, 0 );
                uint tSlaveDepStopIndex  = mSet->get_jac_dof_assembly_map()( tSlaveDofIndex )( tDofDepIndex, 1 );

                // compute jacobian direct dependencies
                if ( tDofType( 0 ) == mResidualDofType( 0 ) )
                {
                    // dRM/dS
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tSlaveDepStartIndex,  tSlaveDepStopIndex } ) -= aWStar * (
                                    tMasterdNdxFlat * tSPConvective->val()( 0 ) * trans( tSlavedNdxFlat ) );

                    // dRS/dS
                    mSet->get_jacobian()(
                            { tSlaveResStartIndex, tSlaveResStopIndex },
                            { tSlaveDepStartIndex, tSlaveDepStopIndex } ) += aWStar * (
                                    tSlavedNdxFlat * tSPConvective->val()( 0 ) * trans( tSlavedNdxFlat ) );
                }
            }
        }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Convective_Normal_Velocity_Ghost::compute_jacobian_and_residual( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Incompressible_NS_Convective_Normal_Velocity_Ghost::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Convective_Normal_Velocity_Ghost::compute_dRdp( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Incompressible_NS_Convective_Normal_Velocity_Ghost::compute_dRdp - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Convective_Normal_Velocity_Ghost::compute_flat_dnNdxn(
                Matrix< DDRMat >  & aFlatdnNdxn,
                mtk::Master_Slave   aIsMaster )
        {
            // get the residual dof type FI (here velocity)
            Field_Interpolator * tVelocityFI = this->get_field_interpolator_manager( aIsMaster )->
                    get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get number of fields
            uint tNumFields = tVelocityFI->get_number_of_fields();

            // get flat dnNdxn (dnNdxn . normal)
            Matrix< DDRMat > tdnNdxn = trans( tVelocityFI->dnNdxn( 1 ) ) * mNormal;
            uint tNumRows = tdnNdxn.n_rows();
            uint tNumCols = tdnNdxn.n_cols();

            // set size for block flat dnNdxn (dnNdxn . normal)
            aFlatdnNdxn.set_size( tNumFields * tNumRows, tNumFields * tNumCols, 0.0 );

            // loop over the number of fields
            for( uint iField = 0; iField < tNumFields; iField++ )
            {
                // fill block flat dnNdxn
                aFlatdnNdxn(
                        { iField * tNumRows, ( iField + 1 ) * tNumRows - 1 },
                        { iField * tNumCols, ( iField + 1 ) * tNumCols - 1 } ) =
                                tdnNdxn.matrix_data();
            }
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
