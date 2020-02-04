
#include "cl_FEM_IWG_Isotropic_Spatial_Diffusion_Ghost.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_trans.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
        IWG_Isotropic_Spatial_Diffusion_Ghost::IWG_Isotropic_Spatial_Diffusion_Ghost()
        {
            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "GhostDisplOrder1" ] = IWG_Stabilization_Type::GHOST_DISPL_1;
            mStabilizationMap[ "GhostDisplOrder2" ] = IWG_Stabilization_Type::GHOST_DISPL_2;
            mStabilizationMap[ "GhostDisplOrder3" ] = IWG_Stabilization_Type::GHOST_DISPL_3;
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Spatial_Diffusion_Ghost::compute_residual( real aWStar )
        {
#ifdef DEBUG
            // check master and slave field interpolators
            this->check_field_interpolators( mtk::Master_Slave::MASTER );
            this->check_field_interpolators( mtk::Master_Slave::SLAVE );
#endif

            // get the master field interpolator for the residual dof type
            Field_Interpolator * tFIMaster = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get the slave field interpolator for the residual dof type
            Field_Interpolator * tFISlave  = mSlaveFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // FIXME the order should be set differently
            switch ( tFIMaster->get_space_interpolation_order() )
            {
                case( mtk::Interpolation_Order::LINEAR ):
                {
                    mOrder = 1;
                    break;
                }
                case( mtk::Interpolation_Order::QUADRATIC ):
                {
                    mOrder = 2;
                    break;
                }
                case( mtk::Interpolation_Order::CUBIC ):
                {
                    mOrder = 3;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "IWG_Isotropic_Spatial_Diffusion_Virtual_Work_Ghost::compute_residual - order not supported");
                    break;
                }
            }

            // get master index for residual dof type
            uint tDofIndexMaster = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );

            // get slave index for residual dof type
            uint tDofIndexSlave  = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::SLAVE );

            // loop over the interpolation order
            for ( uint iOrder = 1; iOrder <= mOrder; iOrder++ )
            {
                 // get normal matrix
                 Matrix< DDRMat > tNormalMatrix = this->get_normal_matrix( iOrder );

                 // premultiply common terms
                 Matrix< DDRMat > tPreMultiply = mStabilizationParam( iOrder - 1 )->val()( 0 )                                                         // penalty
                                               * trans( tNormalMatrix ) * tNormalMatrix                                // normals
                                               * ( tFIMaster->gradx( iOrder ) - tFISlave->gradx( iOrder ) ); // jump in iOrder order spatial gradient

                 // compute master residual
                 mSet->get_residual()( { mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 1 ) }, { 0, 0 } )
                 +=   trans( tFIMaster->dnNdxn( iOrder ) ) * tPreMultiply * aWStar;

                 // compute slave residual
                 mSet->get_residual()( { mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 1 ) }, { 0, 0 } )
                 += - trans( tFISlave->dnNdxn( iOrder ) )  * tPreMultiply* aWStar;
            }
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Spatial_Diffusion_Ghost::compute_jacobian( real aWStar )
        {

#ifdef DEBUG
            // check master and slave field interpolators
            this->check_field_interpolators( mtk::Master_Slave::MASTER );
            this->check_field_interpolators( mtk::Master_Slave::SLAVE );
#endif

            // get the master field interpolator for residual dof type
            Field_Interpolator * tFIMaster = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get the slave field interpolator for residual dof type
            Field_Interpolator * tFISlave  = mSlaveFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // FIXME the order should be set differently
            switch ( tFIMaster->get_space_interpolation_order() )
            {
                case( mtk::Interpolation_Order::LINEAR ):
                {
                    mOrder = 1;
                    break;
                }
                case( mtk::Interpolation_Order::QUADRATIC ):
                {
                    mOrder = 2;
                    break;
                }
                case( mtk::Interpolation_Order::CUBIC ):
                {
                    mOrder = 3;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "IWG_Isotropic_Spatial_Diffusion_Virtual_Work_Ghost::compute_residual - order not supported");
                    break;
                }
            }

            // get master index for residual dof type
            uint tDofIndexMaster = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );

            // get slave dof for residual dof type
            uint tDofIndexSlave  = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::SLAVE );

            // get number of master dependencies
            uint tMasterNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

            // loop over the interpolation orders
            for ( uint iOrder = 1; iOrder <= mOrder; iOrder++ )
            {
                // get normal matrix
                Matrix< DDRMat > tNormalMatrix = this->get_normal_matrix( iOrder );

                // premultiply common terms
                Matrix< DDRMat > tPreMultiply = mStabilizationParam( iOrder - 1 )->val()( 0 ) * trans( tNormalMatrix ) * tNormalMatrix;

                // compute Jacobian direct dependencies
                if ( mResidualDofTypeRequested )
                {
                    // dRM/dM
                    mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 1 ) },
                                          { mSet->get_jac_dof_assembly_map()( tDofIndexMaster )( tDofIndexMaster, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndexMaster )( tDofIndexMaster, 1 ) } )
                    +=   trans( tFIMaster->dnNdxn( iOrder ) ) * tPreMultiply * tFIMaster->dnNdxn( iOrder ) * aWStar;

                    // dRM/dS
                    mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 1 ) },
                                          { mSet->get_jac_dof_assembly_map()( tDofIndexMaster )( tDofIndexSlave, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndexMaster )( tDofIndexSlave, 1 ) } )
                    += - trans( tFIMaster->dnNdxn( iOrder ) ) * tPreMultiply * tFISlave->dnNdxn( iOrder ) * aWStar;

                    // dRS/dM
                    mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 1 ) },
                                          { mSet->get_jac_dof_assembly_map()( tDofIndexSlave )( tDofIndexMaster, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndexSlave )( tDofIndexMaster, 1 ) } )
                    += - trans( tFISlave->dnNdxn( iOrder ) )  * tPreMultiply * tFIMaster->dnNdxn( iOrder ) * aWStar;

                    // dRS/dS
                    mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 1 ) },
                                          { mSet->get_jac_dof_assembly_map()( tDofIndexSlave )( tDofIndexSlave, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndexSlave )( tDofIndexSlave, 1 ) } )
                    +=   trans( tFISlave->dnNdxn( iOrder ) )  * tPreMultiply * tFISlave->dnNdxn( iOrder ) * aWStar;
                }

                // compute the jacobian for indirect dof dependencies through master properties
                for( uint iDOF = 0; iDOF < tMasterNumDofDependencies; iDOF++ )
                {
                    // get the dof type
                    Cell< MSI::Dof_Type > tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                    // get the index for the dof type
                    sint tIndexDep = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );

                    // if dependency on the dof type
                    if ( mStabilizationParam( iOrder - 1 )->check_dof_dependency( tDofType ) )
                    {
                        // premultiply common terms
                        Matrix< DDRMat > tPreMultiply2 = trans( tNormalMatrix ) * tNormalMatrix
                                                         * ( tFIMaster->gradx( iOrder ) - tFISlave->gradx( iOrder ) )
                                                         * mStabilizationParam( iOrder - 1 )->dSPdMasterDOF( tDofType );

                        // add contribution to jacobian
                        mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 1 ) },
                                              { mSet->get_jac_dof_assembly_map()( tDofIndexMaster )( tIndexDep, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndexMaster )( tIndexDep, 1 ) } )
                        +=   trans( tFIMaster->dnNdxn( iOrder ) )  * tPreMultiply2 * aWStar;

                        mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 1 ) },
                                              { mSet->get_jac_dof_assembly_map()( tDofIndexSlave )( tIndexDep, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndexSlave )( tIndexDep, 1 ) } )
                        += - trans( tFISlave->dnNdxn( iOrder ) ) * tPreMultiply2 * aWStar;
                    }
                }
            }
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Spatial_Diffusion_Ghost::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Spatial_Diffusion_Ghost::compute_jacobian_and_residual - Not implemented." );
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Spatial_Diffusion_Ghost::compute_drdpdv( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Spatial_Diffusion_Ghost::compute_drdpdv - This function does nothing.");
        }

//------------------------------------------------------------------------------
        Matrix< DDRMat > IWG_Isotropic_Spatial_Diffusion_Ghost::get_normal_matrix ( uint aOrderGhost )
        {
            // init the normal matrix
            Matrix< DDRMat > tNormalMatrix;

            // get spatial dimensions
            uint tSpaceDim = mNormal.numel();

            // switch on the ghost order
            switch( aOrderGhost )
            {
                case ( 1 ):
                {
                    switch ( tSpaceDim )
                    {
                        case ( 2 ):
                        {
                            tNormalMatrix = trans( mNormal );
                            break;
                        }
                        case ( 3 ):
                        {
                            tNormalMatrix = trans( mNormal );
                            break;
                        }
                        default:
                        {
                            MORIS_ERROR( false, "IWG_Isotropic_Spatial_Diffusion_Ghost::get_normal_matrix - Spatial dimensions can only be 2, 3." );
                            break;
                        }
                    }
                    break;
                }

                case ( 2 ):
                {
                    switch ( tSpaceDim )
                    {
                        case ( 2 ):
                        {
                            // set the normal matrix size
                            tNormalMatrix.set_size( 2, 3, 0.0 );

                            // fill the normal matrix
                            tNormalMatrix( 0, 0 ) = mNormal( 0 );
                            tNormalMatrix( 1, 1 ) = mNormal( 1 );

                            tNormalMatrix( 0, 2 ) = mNormal( 1 );
                            tNormalMatrix( 1, 2 ) = mNormal( 0 );

                            break;
                        }
                        case ( 3 ):
                        {
                            // set the normal matrix size
                            tNormalMatrix.set_size( 3, 6, 0.0 );

                            // fill the normal matrix
                            tNormalMatrix( 0, 0 ) = mNormal( 0 );
                            tNormalMatrix( 1, 1 ) = mNormal( 1 );
                            tNormalMatrix( 2, 2 ) = mNormal( 2 );

                            tNormalMatrix( 1, 3 ) = mNormal( 2 );
                            tNormalMatrix( 2, 3 ) = mNormal( 1 );

                            tNormalMatrix( 0, 4 ) = mNormal( 2 );
                            tNormalMatrix( 2, 4 ) = mNormal( 0 );

                            tNormalMatrix( 0, 5 ) = mNormal( 1 );
                            tNormalMatrix( 1, 5 ) = mNormal( 0 );

                            break;
                        }
                        default:
                        {
                            MORIS_ERROR( false, "IWG_Isotropic_Spatial_Diffusion_Ghost::get_normal_matrix - Spatial dimensions can only be 2, 3." );
                            break;
                        }
                    }
                    break;
                }

                case ( 3 ):
                {
                    switch ( tSpaceDim )
                    {
                        case ( 2 ):
                        {
                            // set the normal matrix size
                            tNormalMatrix.set_size( 3, 4, 0.0 );

                            tNormalMatrix( 0, 0 ) = mNormal( 0 );
                            tNormalMatrix( 1, 1 ) = mNormal( 1 );

                            tNormalMatrix( 0, 2 ) = mNormal( 1 );
                            tNormalMatrix( 1, 3 ) = mNormal( 0 );

                            real tSqrtOf2 = std::sqrt( 2 );

                            tNormalMatrix( 2, 2 ) = tSqrtOf2 * mNormal( 0 );
                            tNormalMatrix( 2, 3 ) = tSqrtOf2 * mNormal( 1 );

                            break;
                        }
                        case ( 3 ):
                        {
                            // set the normal matrix size
                            tNormalMatrix.set_size( 6, 10, 0.0 );

                            tNormalMatrix( 0, 0 ) = mNormal( 0 );
                            tNormalMatrix( 1, 1 ) = mNormal( 1 );
                            tNormalMatrix( 2, 2 ) = mNormal( 2 );

                            tNormalMatrix( 0, 3 ) = mNormal( 1 );
                            tNormalMatrix( 0, 4 ) = mNormal( 2 );

                            tNormalMatrix( 1, 5 ) = mNormal( 0 );
                            tNormalMatrix( 1, 6 ) = mNormal( 2 );

                            tNormalMatrix( 2, 7 ) = mNormal( 0 );
                            tNormalMatrix( 2, 8 ) = mNormal( 1 );

                            real tSqrtOf2 = std::sqrt( 2 );

                            tNormalMatrix( 3, 3 ) = tSqrtOf2 * mNormal( 0 );
                            tNormalMatrix( 3, 5 ) = tSqrtOf2 * mNormal( 1 );
                            tNormalMatrix( 3, 9 ) = tSqrtOf2 * mNormal( 2 );

                            tNormalMatrix( 4, 6 ) = tSqrtOf2 * mNormal( 1 );
                            tNormalMatrix( 4, 8 ) = tSqrtOf2 * mNormal( 2 );
                            tNormalMatrix( 4, 9 ) = tSqrtOf2 * mNormal( 0 );

                            tNormalMatrix( 5, 4 ) = tSqrtOf2 * mNormal( 0 );
                            tNormalMatrix( 5, 7 ) = tSqrtOf2 * mNormal( 2 );
                            tNormalMatrix( 5, 9 ) = tSqrtOf2 * mNormal( 1 );

                            break;
                        }
                        default:
                        {
                            MORIS_ERROR( false, "IWG_Isotropic_Spatial_Diffusion_Ghost::get_normal_matrix - Spatial dimensions can only be 2, 3." );
                            break;
                        }
                    }
                    break;
                }

                default:
                {
                    MORIS_ERROR( false, "IWG_Isotropic_Spatial_Diffusion_Ghost::get_normal_matrix - order not supported." );
                    break;
                }
            }

            return tNormalMatrix;
        }

//------------------------------------------------------------------------------
//        void IWG_Isotropic_Spatial_Diffusion_Ghost::compute_residual( moris::Cell< Matrix< DDRMat > > & aResidual )
//        {
//            // check, if order of shape functions is supported
//            MORIS_ERROR( mOrder < 4, " Ghost stabilization for this order of shape fncts. not supported yet. " );
//
//            // check master and slave field interpolators
//            this->check_field_interpolators( mtk::Master_Slave::MASTER );
//            this->check_field_interpolators( mtk::Master_Slave::SLAVE );
//
//            // set field interpolators for master and slave
//            Field_Interpolator* tTemp_Left  = mMasterFI( 0 );
//            Field_Interpolator* tTemp_Right = mSlaveFI( 0 );
//
//            // compute the residual
//            uint tResSize_Left  = tTemp_Left ->get_number_of_space_time_coefficients();
//            uint tResSize_Right = tTemp_Right->get_number_of_space_time_coefficients();
//
//            // set residual size
//            this->set_residual( aResidual );
//
//            // need a residual that has the size of the two IP elements combined tResSize = tResSize_Left + tResSize_Right
//            Matrix< DDRMat > tResidual;
//            tResidual.set_size( tResSize_Left + tResSize_Right, 1, 0.0);
//
//            // --------- Ghost for linear shape functions ----------------------------------------------------------------
//
//            // assemble matrix containing both left and right shape functions and their derivatives wrt. (x,y,z)
//
//            // get derivatives from left and right shape functions
//            Matrix <DDRMat> tBx_Left = tTemp_Left->Bx();
//            uint tNumOfShapeFunctions_Left = tBx_Left.n_cols();
//
//            Matrix <DDRMat> tBx_Right = tTemp_Right->Bx();
//            uint tNumOfShapeFunctions_Right = tBx_Right.n_cols();
//
//            // put left and right B-Matrices next to each other in Block matrix Bx_comined = [Bx_L, Bx_R]
//            Matrix <DDRMat> tBx_combined( 3, tNumOfShapeFunctions_Left + tNumOfShapeFunctions_Right, 0.0 );
//            for(uint i = 0; i < tNumOfShapeFunctions_Left; i++ )
//            {
//                tBx_combined( 0, i ) = tBx_Left( 0, i );
//                tBx_combined( 1, i ) = tBx_Left( 1, i );
//                tBx_combined( 2, i ) = tBx_Left( 2, i );
//            }
//            for(uint i = 0; i < tNumOfShapeFunctions_Right; i++)
//            {
//                tBx_combined( 0, tNumOfShapeFunctions_Left + i ) = - tBx_Right( 0, i );
//                tBx_combined( 1, tNumOfShapeFunctions_Left + i ) = - tBx_Right( 1, i );
//                tBx_combined( 2, tNumOfShapeFunctions_Left + i ) = - tBx_Right( 2, i );
//            }
//
//            // create vector containing jump in gradient
//            Matrix <DDRMat> tGradX_combined = tTemp_Left->gradx( 1 ) - tTemp_Right->gradx( 1 );
//
//            // compute residual ---------------
//            tResidual.matrix_data() += mGammaGhost * mMeshParameter // * mKappa                                // scaling parameters
//                                     * trans( tBx_combined ) * mNormal * trans( mNormal ) * tGradX_combined;   // matrices & vectors
//
//            // --------- Ghost for quadratic shape functions -------------------------------------------------------------
//            if ( mOrder >= 2 )
//            {
//                // assemble matrix containing both left and right shape functions and their 2nd derivatives wrt. (x,y,z)
//
//                // get derivatives from left and right shape functions
//                Matrix <DDRMat> tB2x_Left = tTemp_Left->d2Ndx2();
//                uint tNumOfShapeFunctions_Left = tB2x_Left.n_cols();
//
//                Matrix <DDRMat> tB2x_Right = tTemp_Right->d2Ndx2();
//                uint tNumOfShapeFunctions_Right = tB2x_Right.n_cols();
//
//                // create combined B-Matrix, put left and right B2x*n next to each other B2x_comined = [B2x_L, -B2x_R]
//                Matrix <DDRMat> tB2x_combined( 6, tNumOfShapeFunctions_Left + tNumOfShapeFunctions_Right, 0.0 );
//
//                // assemble into combined matrix
//                for(uint i = 0; i < tNumOfShapeFunctions_Left; i++ )
//                {
//                    tB2x_combined( 0, i ) = tB2x_Left( 0, i );
//                    tB2x_combined( 1, i ) = tB2x_Left( 1, i );
//                    tB2x_combined( 2, i ) = tB2x_Left( 2, i );
//                    tB2x_combined( 3, i ) = tB2x_Left( 3, i );
//                    tB2x_combined( 4, i ) = tB2x_Left( 4, i );
//                    tB2x_combined( 5, i ) = tB2x_Left( 5, i );
//                }
//
//                // assemble into combined matrix
//                for(uint i = 0; i < tNumOfShapeFunctions_Right; i++ )
//                {
//                    tB2x_combined( 0, tNumOfShapeFunctions_Left + i ) = - tB2x_Right( 0, i );
//                    tB2x_combined( 1, tNumOfShapeFunctions_Left + i ) = - tB2x_Right( 1, i );
//                    tB2x_combined( 2, tNumOfShapeFunctions_Left + i ) = - tB2x_Right( 2, i );
//                    tB2x_combined( 3, tNumOfShapeFunctions_Left + i ) = - tB2x_Right( 3, i );
//                    tB2x_combined( 4, tNumOfShapeFunctions_Left + i ) = - tB2x_Right( 4, i );
//                    tB2x_combined( 5, tNumOfShapeFunctions_Left + i ) = - tB2x_Right( 5, i );
//                }
//
//                // assemble gradient vector
//                Matrix <DDRMat> tGrad2X_combined = tTemp_Left->gradx( 2 ) - tTemp_Right->gradx( 2 );
//
//                // get normal matrix
//                Matrix<DDRMat> tNormalMatrix_2 = this->get_normal_matrix(2);
//
//                // compute residual ---------------
//                tResidual.matrix_data() += mGammaGhost * ( std::pow( mMeshParameter,3 ) ) // * mKappa                              // scaling parameters (p=2)
//                                         * trans( tB2x_combined ) * trans( tNormalMatrix_2 ) * tNormalMatrix_2 * tGrad2X_combined; // matrices & vectors (p=2)
//            }
//
//            // --------- Ghost for cubic shape functions -------------------------------------------------------------
//            if ( mOrder >= 3 )
//            {
//                // assemble matrix containing both left and right shape functions and their 3nd derivatives wrt. (x,y,z)
//
//                // get derivatives from left and right shape functions
//                Matrix <DDRMat> tB3x_Left = tTemp_Left->d3Ndx3();
//                uint tNumOfShapeFunctions_Left = tB3x_Left.n_cols();
//
//                Matrix <DDRMat> tB3x_Right = tTemp_Right->d3Ndx3();
//                uint tNumOfShapeFunctions_Right = tB3x_Right.n_cols();
//
//                // create combined B-Matrix, put left and right B2x*n next to each other B3x_comined = [B3x_L, -B3x_R]
//                Matrix <DDRMat> tB3x_combined( 10, tNumOfShapeFunctions_Left + tNumOfShapeFunctions_Right, 0.0 );
//
//                // assemble into combined matrix
//                for( uint i = 0; i < tNumOfShapeFunctions_Left; i++ )
//                {
//                    tB3x_combined( 0, i ) = tB3x_Left( 0, i );
//                    tB3x_combined( 1, i ) = tB3x_Left( 1, i );
//                    tB3x_combined( 2, i ) = tB3x_Left( 2, i );
//                    tB3x_combined( 3, i ) = tB3x_Left( 3, i );
//                    tB3x_combined( 4, i ) = tB3x_Left( 4, i );
//                    tB3x_combined( 5, i ) = tB3x_Left( 5, i );
//                    tB3x_combined( 6, i ) = tB3x_Left( 6, i );
//                    tB3x_combined( 7, i ) = tB3x_Left( 7, i );
//                    tB3x_combined( 8, i ) = tB3x_Left( 8, i );
//                    tB3x_combined( 9, i ) = tB3x_Left( 9, i );
//                }
//
//                // assemble into combined matrix
//                for( uint i = 0; i < tNumOfShapeFunctions_Right; i++ )
//                {
//                    tB3x_combined( 0, tNumOfShapeFunctions_Left + i ) = - tB3x_Right( 0, i );
//                    tB3x_combined( 1, tNumOfShapeFunctions_Left + i ) = - tB3x_Right( 1, i );
//                    tB3x_combined( 2, tNumOfShapeFunctions_Left + i ) = - tB3x_Right( 2, i );
//                    tB3x_combined( 3, tNumOfShapeFunctions_Left + i ) = - tB3x_Right( 3, i );
//                    tB3x_combined( 4, tNumOfShapeFunctions_Left + i ) = - tB3x_Right( 4, i );
//                    tB3x_combined( 5, tNumOfShapeFunctions_Left + i ) = - tB3x_Right( 5, i );
//                    tB3x_combined( 6, tNumOfShapeFunctions_Left + i ) = - tB3x_Right( 6, i );
//                    tB3x_combined( 7, tNumOfShapeFunctions_Left + i ) = - tB3x_Right( 7, i );
//                    tB3x_combined( 8, tNumOfShapeFunctions_Left + i ) = - tB3x_Right( 8, i );
//                    tB3x_combined( 9, tNumOfShapeFunctions_Left + i ) = - tB3x_Right( 9, i );
//                }
//
//                // assemble gradient vector
//                 Matrix <DDRMat> tGrad3X_combined = tTemp_Left->gradx( 3 ) - tTemp_Right->gradx( 3 );
//
//                // get normal matrix
//                Matrix<DDRMat> tNormalMatrix_3 = this->get_normal_matrix( 3 );
//
//                // compute residual ---------------
//                tResidual.matrix_data() += mGammaGhost * ( std::pow( mMeshParameter, 5 ) ) // * mKappa                             // scaling parameters (p = 3)
//                                         * trans( tB3x_combined ) * trans( tNormalMatrix_3 ) * tNormalMatrix_3 * tGrad3X_combined; // matrices & vectors (p = 3)
//                print( tResidual, "tResidual" );
//            }
//            aResidual( 0 ) = tResidual({ 0, tResSize_Left-1 },{ 0, 0 });
//            aResidual( 1 ) = tResidual({ tResSize_Left, tResSize_Left + tResSize_Right - 1 },{ 0, 0 });
//        }

//        void IWG_Isotropic_Spatial_Diffusion_Ghost::compute_jacobian( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians )
//        {
//            // check, if order of shape functions is supported
//            MORIS_ERROR( mOrder < 4, " IWG_Isotropic_Spatial_Diffusion_Ghost::compute_jacobian - Ghost stabilization for this order of shape fncts. not supported yet. " );
//
//            // check master and slave field interpolators
//            this->check_field_interpolators( mtk::Master_Slave::MASTER );
//            this->check_field_interpolators( mtk::Master_Slave::SLAVE );
//
//            // set field interpolators for left and right
//            Field_Interpolator* tTemp_Left  = mMasterFI( 0 );
//            Field_Interpolator* tTemp_Right = mSlaveFI( 0 );
//
//            // set the jacobian size
//            this->set_jacobian( aJacobians );
//
//            // compute the jacobian j_T_T
//            uint tNumOfBases_Left = tTemp_Left->get_number_of_space_time_bases();
//            uint tNumOfBases_Right = tTemp_Right->get_number_of_space_time_bases();
//
//            Matrix< DDRMat > tJacobians( tNumOfBases_Left + tNumOfBases_Right, tNumOfBases_Left + tNumOfBases_Right, 0.0 );
//
//            // --------- Ghost for linear shape functions ----------------------------------------------------------------
//
//            // assemble matrix containing both left and right shape functions and their derivatives wrt. (x,y,z)
//
//            // get derivatives from left and right shape functions
//            Matrix <DDRMat> tBx_Left = tTemp_Left->Bx();
//            uint tNumOfShapeFunctions_Left = tBx_Left.n_cols();
//
//            Matrix <DDRMat> tBx_Right = tTemp_Right->Bx();
//            uint tNumOfShapeFunctions_Right = tBx_Right.n_cols();
//
//            // put left and right B-Matrices next to each other in Block matrix Bx_comined = [Bx_L, -Bx_R]
//            Matrix <DDRMat> tBx_combined( 3, tNumOfShapeFunctions_Left + tNumOfShapeFunctions_Right, 0.0 );
//            for( uint i = 0; i < tNumOfShapeFunctions_Left; i++ )
//            {
//                tBx_combined( 0, i ) = tBx_Left( 0, i );
//                tBx_combined( 1, i ) = tBx_Left( 1, i );
//                tBx_combined( 2, i ) = tBx_Left( 2, i );
//            }
//            for( uint i = 0; i < tNumOfShapeFunctions_Right; i++ )
//            {
//                tBx_combined( 0, tNumOfShapeFunctions_Left + i ) = - tBx_Right( 0, i );
//                tBx_combined( 1, tNumOfShapeFunctions_Left + i ) = - tBx_Right( 1, i );
//                tBx_combined( 2, tNumOfShapeFunctions_Left + i ) = - tBx_Right( 2, i );
//            }
//
//            // compute Jacobian ---------------
//            tJacobians.matrix_data() += mGammaGhost * mMeshParameter // * mKappa                           // scaling parameters
//                                     * trans( tBx_combined ) * mNormal * trans( mNormal ) * tBx_combined; // matrices & vectors
//
//            // --------- Ghost for quadratic shape functions -------------------------------------------------------------
//            if (mOrder >= 2)
//            {
//                // assemble matrix containing both left and right shape functions and their 2nd derivatives wrt. (x,y,z)
//
//                // get derivatives from left and right shape functions
//                Matrix <DDRMat> tB2x_Left = tTemp_Left->d2Ndx2();
//                uint tNumOfShapeFunctions_Left = tB2x_Left.n_cols();
//
//                Matrix <DDRMat> tB2x_Right = tTemp_Right->d2Ndx2();
//                uint tNumOfShapeFunctions_Right = tB2x_Right.n_cols();
//
//                // create combined B-Matrix, put left and right B2x*n next to each other B2x_comined = [B2x_L, -B2x_R]
//                Matrix <DDRMat> tB2x_combined( 6, tNumOfShapeFunctions_Left + tNumOfShapeFunctions_Right, 0.0 );
//
//                // assemble into combined matrix
//                for( uint i = 0; i < tNumOfShapeFunctions_Left; i++ )
//                {
//                    tB2x_combined( 0, i ) = tB2x_Left( 0, i );
//                    tB2x_combined( 1, i ) = tB2x_Left( 1, i );
//                    tB2x_combined( 2, i ) = tB2x_Left( 2, i );
//                    tB2x_combined( 3, i ) = tB2x_Left( 3, i );
//                    tB2x_combined( 4, i ) = tB2x_Left( 4, i );
//                    tB2x_combined( 5, i ) = tB2x_Left( 5, i );
//                }
//
//                //  assemble into combined matrix
//                for( uint i = 0; i < tNumOfShapeFunctions_Right; i++ )
//                {
//                    tB2x_combined( 0, tNumOfShapeFunctions_Left + i ) = - tB2x_Right( 0, i );
//                    tB2x_combined( 1, tNumOfShapeFunctions_Left + i ) = - tB2x_Right( 1, i );
//                    tB2x_combined( 2, tNumOfShapeFunctions_Left + i ) = - tB2x_Right( 2, i );
//                    tB2x_combined( 3, tNumOfShapeFunctions_Left + i ) = - tB2x_Right( 3, i );
//                    tB2x_combined( 4, tNumOfShapeFunctions_Left + i ) = - tB2x_Right( 4, i );
//                    tB2x_combined( 5, tNumOfShapeFunctions_Left + i ) = - tB2x_Right( 5, i );
//                }
//
//                // get normal matrix
//                Matrix<DDRMat> tNormalMatrix_2 = this->get_normal_matrix(2);
//
//                // compute Jacobian ---------------
//                tJacobians.matrix_data() += mGammaGhost * ( std::pow( mMeshParameter, 3 ) ) // * mKappa                              // scaling parameters (p=2)
//                                         * trans( tB2x_combined ) * trans( tNormalMatrix_2 ) * tNormalMatrix_2 * tB2x_combined; // matrices & vectors (p=2)
//            }
//
//            // --------- Ghost for cubic shape functions -------------------------------------------------------------
//            if ( mOrder >= 3 )
//            {
//                // assemble matrix containing both left and right shape functions and their 2nd derivatives wrt. (x,y,z)
//
//                // get derivatives from left and right shape functions
//                Matrix <DDRMat> tB3x_Left = tTemp_Left->d3Ndx3();
//                uint tNumOfShapeFunctions_Left = tB3x_Left.n_cols();
//
//                Matrix <DDRMat> tB3x_Right = tTemp_Right->d3Ndx3();
//                uint tNumOfShapeFunctions_Right = tB3x_Right.n_cols();
//
//                // create combined B-Matrix, put left and right B2x*n next to each other B3x_comined = [B3x_L, -B3x_R]
//                Matrix <DDRMat> tB3x_combined( 10, tNumOfShapeFunctions_Left + tNumOfShapeFunctions_Right, 0.0 );
//
//                // assemble into combined matrix
//                for( uint i = 0; i < tNumOfShapeFunctions_Left; i++ )
//                {
//                    tB3x_combined( 0, i ) = tB3x_Left( 0, i );
//                    tB3x_combined( 1, i ) = tB3x_Left( 1, i );
//                    tB3x_combined( 2, i ) = tB3x_Left( 2, i );
//                    tB3x_combined( 3, i ) = tB3x_Left( 3, i );
//                    tB3x_combined( 4, i ) = tB3x_Left( 4, i );
//                    tB3x_combined( 5, i ) = tB3x_Left( 5, i );
//                    tB3x_combined( 6, i ) = tB3x_Left( 6, i );
//                    tB3x_combined( 7, i ) = tB3x_Left( 7, i );
//                    tB3x_combined( 8, i ) = tB3x_Left( 8, i );
//                    tB3x_combined( 9, i ) = tB3x_Left( 9, i );
//                }
//
//                // assemble into combined matrix
//                for( uint i = 0; i < tNumOfShapeFunctions_Right; i++ )
//                {
//                    tB3x_combined( 0, tNumOfShapeFunctions_Left + i ) = - tB3x_Right( 0, i );
//                    tB3x_combined( 1, tNumOfShapeFunctions_Left + i ) = - tB3x_Right( 1, i );
//                    tB3x_combined( 2, tNumOfShapeFunctions_Left + i ) = - tB3x_Right( 2, i );
//                    tB3x_combined( 3, tNumOfShapeFunctions_Left + i ) = - tB3x_Right( 3, i );
//                    tB3x_combined( 4, tNumOfShapeFunctions_Left + i ) = - tB3x_Right( 4, i );
//                    tB3x_combined( 5, tNumOfShapeFunctions_Left + i ) = - tB3x_Right( 5, i );
//                    tB3x_combined( 6, tNumOfShapeFunctions_Left + i ) = - tB3x_Right( 6, i );
//                    tB3x_combined( 7, tNumOfShapeFunctions_Left + i ) = - tB3x_Right( 7, i );
//                    tB3x_combined( 8, tNumOfShapeFunctions_Left + i ) = - tB3x_Right( 8, i );
//                    tB3x_combined( 9, tNumOfShapeFunctions_Left + i ) = - tB3x_Right( 9, i );
//                }
//
//                // get normal matrix
//                Matrix<DDRMat> tNormalMatrix_3 = this->get_normal_matrix(3);
//
//                // compute Jacobian ---------------
//                tJacobians.matrix_data() += mGammaGhost * ( std::pow( mMeshParameter, 5 ) ) // * mKappa                         // scaling parameters (p=2)
//                                         * trans( tB3x_combined ) * trans( tNormalMatrix_3 ) * tNormalMatrix_3 * tB3x_combined; // matrices & vectors (p=2)
//            }
//
//            aJacobians( 0 )( 0 ) = tJacobians( { 0, tNumOfBases_Left - 1 },{ 0, tNumOfBases_Left - 1 });
//            aJacobians( 0 )( 1 ) = tJacobians( { 0, tNumOfBases_Left - 1 },{ tNumOfBases_Left, tNumOfBases_Left + tNumOfBases_Right - 1 });
//
//            aJacobians( 1 )( 0 ) = tJacobians( { tNumOfBases_Left, tNumOfBases_Left + tNumOfBases_Right - 1 },{ 0, tNumOfBases_Left - 1 });
//            aJacobians( 1 )( 1 ) = tJacobians( { tNumOfBases_Left, tNumOfBases_Left + tNumOfBases_Right - 1 },{ tNumOfBases_Left, tNumOfBases_Left + tNumOfBases_Right - 1 });
//        }


//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
