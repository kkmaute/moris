//FEM/INT/src
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG_Incompressible_NS_Convective_Velocity_Ghost.hpp"
//LINALG/src
#include "fn_trans.hpp"
#include "fn_norm.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        IWG_Incompressible_NS_Convective_Velocity_Ghost::IWG_Incompressible_NS_Convective_Velocity_Ghost()
        {
            // set ghost flag
            mIsGhost = true;

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "ConvectiveGhost" ] = static_cast< uint >( IWG_Stabilization_Type::CONVECTIVE_GHOST );
        }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Convective_Velocity_Ghost::compute_residual( real aWStar )
        {
            // check master field interpolators
#ifdef DEBUG
            this->check_field_interpolators();
#endif

            // get master index for residual dof type, indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get slave index for residual dof type, indices for assembly
            uint tSlaveDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::SLAVE );
            uint tSlaveResStartIndex = mSet->get_res_dof_assembly_map()( tSlaveDofIndex )( 0, 0 );
            uint tSlaveResStopIndex  = mSet->get_res_dof_assembly_map()( tSlaveDofIndex )( 0, 1 );

            // get the master field interpolator for the residual dof type
            Field_Interpolator * tFIMaster = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get the slave field interpolator for the residual dof type
            Field_Interpolator * tFISlave  = mSlaveFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get the convective stabilization parameter
            std::shared_ptr< Stabilization_Parameter > & tSPConvective =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::CONVECTIVE_GHOST ) );

            // get flattened derivatives dnNdxn for master and slave
            Matrix< DDRMat > tMasterdNdx;
            this->compute_dnNdxn( tMasterdNdx, mtk::Master_Slave::MASTER );
            Matrix< DDRMat > tSlavedNdx;
            this->compute_dnNdxn( tSlavedNdx, mtk::Master_Slave::SLAVE );

            // premultiply common terms
            Matrix< DDRMat > tConvectivePreMultiply = vectorize(
                    tSPConvective->val()( 0 ) * ( tFIMaster->gradx( 1 ) - tFISlave->gradx( 1 ) ) );

            // compute master residual
            mSet->get_residual()( 0 )(
                    { tMasterResStartIndex, tMasterResStopIndex },
                    { 0, 0 } ) += aWStar * (
                            trans( tMasterdNdx ) * tConvectivePreMultiply );

            // compute slave residual
            mSet->get_residual()( 0 )(
                    { tSlaveResStartIndex, tSlaveResStopIndex },
                    { 0, 0 } ) -= aWStar * (
                            trans( tSlavedNdx ) * tConvectivePreMultiply );

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Incompressible_NS_Convective_Velocity_Ghost::compute_residual - Residual contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Convective_Velocity_Ghost::compute_jacobian( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            // get master index for residual dof type, indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get slave index for residual dof type, indices for assembly
            uint tSlaveDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::SLAVE );
            uint tSlaveResStartIndex = mSet->get_res_dof_assembly_map()( tSlaveDofIndex )( 0, 0 );
            uint tSlaveResStopIndex  = mSet->get_res_dof_assembly_map()( tSlaveDofIndex )( 0, 1 );

            // get the master field interpolator for residual dof type
            Field_Interpolator * tFIMaster = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get the slave field interpolator for residual dof type
            Field_Interpolator * tFISlave  = mSlaveFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get the convective stabilization parameter
            std::shared_ptr< Stabilization_Parameter > & tSPConvective =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::CONVECTIVE_GHOST ) );

            // get flattened derivatives dnNdxn for master and slave
            Matrix< DDRMat > tMasterdNdx;
            this->compute_dnNdxn( tMasterdNdx, mtk::Master_Slave::MASTER );
            Matrix< DDRMat > tSlavedNdx;
            this->compute_dnNdxn( tSlavedNdx, mtk::Master_Slave::SLAVE );

            // get number of master dependencies
            uint tMasterNumDofDependencies = mRequestedMasterGlobalDofTypes.size();
            uint tSlaveNumDofDependencies  = mRequestedSlaveGlobalDofTypes.size();

            // compute the jacobian for indirect dof dependencies through master
            for( uint iDOF = 0; iDOF < tMasterNumDofDependencies; iDOF++ )
            {
                // get the dof type
                Cell< MSI::Dof_Type > & tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                // get the index for the dof type
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tMasterDepStartIndex = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 0 );
                uint tMasterDepStopIndex  = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 1 );

                // compute jacobian direct dependencies
                if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    // dRM/dM
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    tSPConvective->val()( 0 ) * trans( tMasterdNdx ) * tMasterdNdx );

                    // dRS/dM
                    mSet->get_jacobian()(
                            { tSlaveResStartIndex,  tSlaveResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) -= aWStar * (
                                    tSPConvective->val()( 0 ) * trans( tSlavedNdx ) * tMasterdNdx );
                }

                // if stabilization parameter dependency on the dof type
                if ( tSPConvective->check_dof_dependency( tDofType ) )
                {
                    // premultiply common terms
                    Matrix< DDRMat > tConvectivePreMultiply = vectorize(
                            tFIMaster->gradx( 1 ) - tFISlave->gradx( 1 ) ) ;

                    tConvectivePreMultiply = tConvectivePreMultiply * tSPConvective->dSPdMasterDOF( tDofType );

                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    trans( tMasterdNdx ) * tConvectivePreMultiply );

                    mSet->get_jacobian()(
                            { tSlaveResStartIndex,  tSlaveResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) -= aWStar * (
                                    trans( tSlavedNdx ) * tConvectivePreMultiply );
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
                if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    // dRM/dS
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tSlaveDepStartIndex,  tSlaveDepStopIndex } ) -= aWStar * (
                                    tSPConvective->val()( 0 ) * trans( tMasterdNdx ) * tSlavedNdx );

                    // dRS/dS
                    mSet->get_jacobian()(
                            { tSlaveResStartIndex, tSlaveResStopIndex },
                            { tSlaveDepStartIndex, tSlaveDepStopIndex } ) += aWStar * (
                                    tSPConvective->val()( 0 ) * trans( tSlavedNdx ) * tSlavedNdx );
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ) ,
                    "IWG_Ghost_Normal_Field::compute_jacobian - Jacobian contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Convective_Velocity_Ghost::compute_jacobian_and_residual( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Incompressible_NS_Convective_Velocity_Ghost::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Convective_Velocity_Ghost::compute_dRdp( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Incompressible_NS_Convective_Velocity_Ghost::compute_dRdp - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Convective_Velocity_Ghost::compute_dnNdxn(
                Matrix< DDRMat >  & adNdx,
                mtk::Master_Slave   aIsMaster )
        {
            // get the field interpolator for residual dof type
            Field_Interpolator * tFIVelocity =
                    this->get_field_interpolator_manager( aIsMaster )->
                    get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // init size for dnNdtn
            uint tNumField = tFIVelocity->get_number_of_fields();
            uint tNumRow = tFIVelocity->dnNdxn( 1 ).n_rows();
            uint tNumCol = tFIVelocity->dnNdxn( 1 ).n_cols();
            adNdx.set_size( tNumField * tNumRow, tNumField * tNumCol , 0.0 );

            // loop over the fields
            for( uint iField = 0; iField < tNumField; iField++ )
            {
                // fill the matrix for each dimension
                adNdx(
                        { iField * tNumRow, ( iField + 1 ) * tNumRow - 1 },
                        { iField * tNumCol, ( iField + 1 ) * tNumCol - 1 } ) =
                                tFIVelocity->dnNdxn( 1 ).matrix_data();
            }
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
