//FEM/INT/src
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG_Incompressible_NS_Pressure_Ghost.hpp"
//LINALG/src
#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------
        IWG_Incompressible_NS_Pressure_Ghost::IWG_Incompressible_NS_Pressure_Ghost()
        {
            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "PressureGhost" ] = IWG_Stabilization_Type::PRESSURE_GHOST;
        }

        //------------------------------------------------------------------------------
        void IWG_Incompressible_NS_Pressure_Ghost::set_stabilization_parameter(
                std::shared_ptr< Stabilization_Parameter > aStabilizationParameter,
                std::string                                aStabilizationString )
        {
            // check that aStabilizationString makes sense
            std::string tErrMsg =
                    std::string( "IWG_Incompressible_NS_Pressure_Ghost::set_stabilization_parameter - Unknown aStabilizationString: " ) +
                    aStabilizationString;
            MORIS_ERROR( mStabilizationMap.find( aStabilizationString ) != mStabilizationMap.end(), tErrMsg.c_str() );

            // set the stabilization parameter in the stabilization parameter cell
            this->get_stabilization_parameters()( static_cast< uint >( mStabilizationMap[ aStabilizationString ] ) ) = aStabilizationParameter;
        }

        //------------------------------------------------------------------------------
        void IWG_Incompressible_NS_Pressure_Ghost::compute_residual( real aWStar )
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

            // FIXME the order should be set differently
            switch ( tFIMaster->get_space_interpolation_order() )
            {
                case mtk::Interpolation_Order::LINEAR :
                {
                    mOrder = 1;
                    break;
                }
                case mtk::Interpolation_Order::QUADRATIC :
                {
                    mOrder = 2;
                    break;
                }
                case mtk::Interpolation_Order::CUBIC :
                {
                    mOrder = 3;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "IWG_Incompressible_NS_Pressure_Ghost::compute_residual - order not supported");
                    break;
                }
            }

            // loop over the interpolation order
            for ( uint iOrder = 1; iOrder <= mOrder; iOrder++ )
            {
                // get the stabilization parameter
                std::shared_ptr< Stabilization_Parameter > tSPPressure
                = mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::PRESSURE_GHOST ) );

                // set the order for stabilization parameters
                tSPPressure->set_interpolation_order( iOrder );

                // get normal matrix
                Matrix< DDRMat > tFlatNormal;
                this->get_normal_matrix( tFlatNormal, iOrder );

                // premultiply common terms
                Matrix< DDRMat > tPreMultiply
                = tSPPressure->val()( 0 ) * tFlatNormal * ( tFIMaster->gradx( iOrder ) - tFISlave->gradx( iOrder ) ) ;
                tPreMultiply = reshape( tPreMultiply , tPreMultiply.numel(), 1 );

                // get flattened directional derivatives for master and slave
                Matrix< DDRMat > tMasterdNdxFlat;
                this->compute_flat_dnNdxn( tMasterdNdxFlat, iOrder, mtk::Master_Slave::MASTER );
                Matrix< DDRMat > tSlavedNdxFlat;
                this->compute_flat_dnNdxn( tSlavedNdxFlat, iOrder, mtk::Master_Slave::SLAVE );

                // compute master residual
                mSet->get_residual()( 0 )( { tMasterResStartIndex, tMasterResStopIndex }, { 0, 0 } )
                                                += aWStar * ( tMasterdNdxFlat * tPreMultiply );

                // compute slave residual
                mSet->get_residual()( 0 )( { tSlaveResStartIndex, tSlaveResStopIndex }, { 0, 0 } )
                                                -= aWStar *( tSlavedNdxFlat * tPreMultiply );
            }
        }

        //------------------------------------------------------------------------------
        void IWG_Incompressible_NS_Pressure_Ghost::compute_jacobian( real aWStar )
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
            Field_Interpolator * tFIMaster = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get the slave field interpolator for residual dof type
            Field_Interpolator * tFISlave  = mSlaveFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get the interpolation order
            switch ( tFIMaster->get_space_interpolation_order() )
            {
                case mtk::Interpolation_Order::LINEAR :
                {
                    mOrder = 1;
                    break;
                }
                case mtk::Interpolation_Order::QUADRATIC :
                {
                    mOrder = 2;
                    break;
                }
                case mtk::Interpolation_Order::CUBIC :
                {
                    mOrder = 3;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "IWG_Incompressible_NS_Pressure_Ghost::compute_residual - order not supported");
                    break;
                }
            }

            // loop over the interpolation orders
            for ( uint iOrder = 1; iOrder <= mOrder; iOrder++ )
            {
                // get the stabilization parameter
                std::shared_ptr< Stabilization_Parameter > tSPPressure
                = mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::PRESSURE_GHOST ) );

                // set the order for the stabilization parameter
                tSPPressure->set_interpolation_order( iOrder );

                // get flattened normal
                Matrix< DDRMat > tFlatNormal;
                this->get_normal_matrix( tFlatNormal, iOrder );

                // get flattened directional derivatives for master and slave
                Matrix< DDRMat > tMasterdNdxFlat;
                this->compute_flat_dnNdxn( tMasterdNdxFlat, iOrder, mtk::Master_Slave::MASTER );
                Matrix< DDRMat > tSlavedNdxFlat;
                this->compute_flat_dnNdxn( tSlavedNdxFlat, iOrder, mtk::Master_Slave::SLAVE );

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
                        mSet->get_jacobian()( { tMasterResStartIndex, tMasterResStopIndex },
                                { tMasterDepStartIndex, tMasterDepStopIndex } )
                                += aWStar * ( tMasterdNdxFlat * tSPPressure->val()( 0 ) * trans( tMasterdNdxFlat ) );

                        // dRS/dM
                        mSet->get_jacobian()( { tSlaveResStartIndex,  tSlaveResStopIndex },
                                { tMasterDepStartIndex, tMasterDepStopIndex } )
                                -= aWStar * ( tSlavedNdxFlat * tSPPressure->val()( 0 ) * trans( tMasterdNdxFlat ) );
                    }

                    // if stabilization parameter dependency on the dof type
                    if ( tSPPressure->check_dof_dependency( tDofType ) )
                    {
                        // premultiply common terms
                        Matrix< DDRMat > tPreMultiply
                        = tFlatNormal * ( tFIMaster->gradx( iOrder ) - tFISlave->gradx( iOrder ) );
                        tPreMultiply = reshape( tPreMultiply , tPreMultiply.numel(), 1 );
                        tPreMultiply = tPreMultiply * tSPPressure->dSPdMasterDOF( tDofType );

                        // add contribution to jacobian
                        mSet->get_jacobian()( { tMasterResStartIndex, tMasterResStopIndex },
                                { tMasterDepStartIndex, tMasterDepStopIndex } )
                                += aWStar * ( tMasterdNdxFlat * tPreMultiply );

                        mSet->get_jacobian()( { tSlaveResStartIndex,  tSlaveResStopIndex },
                                { tMasterDepStartIndex, tMasterDepStopIndex } )
                                -= aWStar * ( tSlavedNdxFlat  * tPreMultiply );
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
                        mSet->get_jacobian()( { tMasterResStartIndex, tMasterResStopIndex },
                                { tSlaveDepStartIndex,  tSlaveDepStopIndex } )
                                -= aWStar * ( tMasterdNdxFlat * tSPPressure->val()( 0 ) * trans( tSlavedNdxFlat ) );

                        // dRS/dS
                        mSet->get_jacobian()( { tSlaveResStartIndex, tSlaveResStopIndex },
                                { tSlaveDepStartIndex, tSlaveDepStopIndex } )
                                += aWStar * ( tSlavedNdxFlat * tSPPressure->val()( 0 ) * trans( tSlavedNdxFlat ) );
                    }
                }
            }
        }

        //------------------------------------------------------------------------------
        void IWG_Incompressible_NS_Pressure_Ghost::compute_jacobian_and_residual( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Incompressible_NS_Pressure_Ghost::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------
        void IWG_Incompressible_NS_Pressure_Ghost::compute_dRdp( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Incompressible_NS_Pressure_Ghost::compute_dRdp - Not implemented." );
        }

        //------------------------------------------------------------------------------
        void IWG_Incompressible_NS_Pressure_Ghost::compute_flat_dnNdxn(
                Matrix< DDRMat >  & aFlatdnNdxn,
                uint                aOrder,
                mtk::Master_Slave   aIsMaster )
        {
            // get flattened normal
            Matrix< DDRMat > tFlatNormal;
            this->get_normal_matrix( tFlatNormal, aOrder );

            // get the residual dof type FI (here pressure)
            Field_Interpolator * tPressureFI
            = this->get_field_interpolator_manager( aIsMaster )
            ->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get number of fields
            uint tNumFields = tPressureFI->get_number_of_fields();

            // get flat dnNdxn (dnNdxn . normal)
            Matrix< DDRMat > tdnNdxn = trans( tPressureFI->dnNdxn( aOrder ) ) * trans( tFlatNormal );
            uint tNumRows = tdnNdxn.n_rows();
            uint tNumCols = tdnNdxn.n_cols();

            // set size for block flat dnNdxn (dnNdxn . normal)
            aFlatdnNdxn.set_size( tNumFields * tNumRows, tNumFields * tNumCols, 0.0 );

            // loop over the number of fields
            for( uint iField = 0; iField < tNumFields; iField++ )
            {
                // fill block flat dnNdxn
                aFlatdnNdxn({ iField * tNumRows, ( iField + 1 ) * tNumRows - 1 },
                        { iField * tNumCols, ( iField + 1 ) * tNumCols - 1 } )
                                                = tdnNdxn.matrix_data();
            }
        }

        //------------------------------------------------------------------------------
        void IWG_Incompressible_NS_Pressure_Ghost::get_normal_matrix(
                Matrix< DDRMat > & aFlatNormal,
                uint               aOrder )
        {
            // get spatial dimensions
            uint tSpaceDim = mNormal.numel();

            // switch on the ghost order
            switch( aOrder )
            {
                case 1 :
                {
                    switch ( tSpaceDim )
                    {
                        case 2 :
                        {
                            aFlatNormal = trans( mNormal );
                            break;
                        }
                        case 3 :
                        {
                            aFlatNormal = trans( mNormal );
                            break;
                        }
                        default:
                        {
                            MORIS_ERROR( false, "IWG_Incompressible_NS_Pressure_Ghost::get_normal_matrix - Spatial dimensions can only be 2, 3." );
                            break;
                        }
                    }
                    break;
                }

                case 2 :
                {
                    switch ( tSpaceDim )
                    {
                        case 2 :
                        {
                            // set the normal matrix size
                            aFlatNormal.set_size( 2, 3, 0.0 );

                            // fill the normal matrix
                            aFlatNormal( 0, 0 ) = mNormal( 0 );
                            aFlatNormal( 1, 1 ) = mNormal( 1 );

                            aFlatNormal( 0, 2 ) = mNormal( 1 );
                            aFlatNormal( 1, 2 ) = mNormal( 0 );

                            break;
                        }
                        case 3 :
                        {
                            // set the normal matrix size
                            aFlatNormal.set_size( 3, 6, 0.0 );

                            // fill the normal matrix
                            aFlatNormal( 0, 0 ) = mNormal( 0 );
                            aFlatNormal( 1, 1 ) = mNormal( 1 );
                            aFlatNormal( 2, 2 ) = mNormal( 2 );

                            aFlatNormal( 1, 3 ) = mNormal( 2 );
                            aFlatNormal( 2, 3 ) = mNormal( 1 );

                            aFlatNormal( 0, 4 ) = mNormal( 2 );
                            aFlatNormal( 2, 4 ) = mNormal( 0 );

                            aFlatNormal( 0, 5 ) = mNormal( 1 );
                            aFlatNormal( 1, 5 ) = mNormal( 0 );

                            break;
                        }
                        default:
                        {
                            MORIS_ERROR( false, "IWG_Incompressible_NS_Pressure_Ghost::get_normal_matrix - Spatial dimensions can only be 2, 3." );
                            break;
                        }
                    }
                    break;
                }

                case 3 :
                {
                    switch ( tSpaceDim )
                    {
                        case 2 :
                        {
                            // set the normal matrix size
                            aFlatNormal.set_size( 3, 4, 0.0 );

                            aFlatNormal( 0, 0 ) = mNormal( 0 );
                            aFlatNormal( 1, 1 ) = mNormal( 1 );

                            aFlatNormal( 0, 2 ) = mNormal( 1 );
                            aFlatNormal( 1, 3 ) = mNormal( 0 );

                            real tSqrtOf2 = std::sqrt( 2 );

                            aFlatNormal( 2, 2 ) = tSqrtOf2 * mNormal( 0 );
                            aFlatNormal( 2, 3 ) = tSqrtOf2 * mNormal( 1 );
                            break;
                        }
                        case 3 :
                        {
                            // set the normal matrix size
                            aFlatNormal.set_size( 6, 10, 0.0 );

                            aFlatNormal( 0, 0 ) = mNormal( 0 );
                            aFlatNormal( 1, 1 ) = mNormal( 1 );
                            aFlatNormal( 2, 2 ) = mNormal( 2 );

                            aFlatNormal( 0, 3 ) = mNormal( 1 );
                            aFlatNormal( 0, 4 ) = mNormal( 2 );

                            aFlatNormal( 1, 5 ) = mNormal( 0 );
                            aFlatNormal( 1, 6 ) = mNormal( 2 );

                            aFlatNormal( 2, 7 ) = mNormal( 0 );
                            aFlatNormal( 2, 8 ) = mNormal( 1 );

                            real tSqrtOf2 = std::sqrt( 2 );

                            aFlatNormal( 3, 3 ) = tSqrtOf2 * mNormal( 0 );
                            aFlatNormal( 3, 5 ) = tSqrtOf2 * mNormal( 1 );
                            aFlatNormal( 3, 9 ) = tSqrtOf2 * mNormal( 2 );

                            aFlatNormal( 4, 6 ) = tSqrtOf2 * mNormal( 1 );
                            aFlatNormal( 4, 8 ) = tSqrtOf2 * mNormal( 2 );
                            aFlatNormal( 4, 9 ) = tSqrtOf2 * mNormal( 0 );

                            aFlatNormal( 5, 4 ) = tSqrtOf2 * mNormal( 0 );
                            aFlatNormal( 5, 7 ) = tSqrtOf2 * mNormal( 2 );
                            aFlatNormal( 5, 9 ) = tSqrtOf2 * mNormal( 1 );
                            break;
                        }
                        default:
                        {
                            MORIS_ERROR( false, "IWG_Incompressible_NS_Pressure_Ghost::get_normal_matrix - Spatial dimensions can only be 2, 3." );
                            break;
                        }
                    }
                    break;
                }

                default:
                {
                    MORIS_ERROR( false, "IWG_Incompressible_NS_Pressure_Ghost::get_normal_matrix - order not supported." );
                    break;
                }
            }
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
