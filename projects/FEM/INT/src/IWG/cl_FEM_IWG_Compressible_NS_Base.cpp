/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Compressible_NS_Base.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG_Compressible_NS_Base.hpp"
#include "fn_FEM_IWG_Compressible_NS.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"

// debug - output to hdf5
#include "paths.hpp"
#include "HDF5_Tools.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------

        void
        IWG_Compressible_NS_Base::reset_spec_eval_flags()
        {
            // reset eval flags
            mYEval      = true;
            mdYdtEval   = true;
            mdYdxEval   = true;
            md2Ydx2Eval = true;

            mWEval      = true;
            mdWdtEval   = true;
            mdWdxEval   = true;
            md2Wdx2Eval = true;

            mWtransEval    = true;
            mdWtransdtEval = true;
            mdWtransdxEval = true;

            mAEval    = true;
            mKEval    = true;
            mKijiEval = true;

            mCEval    = true;
            mdCdYEval = true;

            // reset flags for child
            this->reset_child_eval_flags();
        }

        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------

        void
        IWG_Compressible_NS_Base::assemble_residual( const Matrix< DDRMat >& aStdRes )
        {
            // check that the size of the passed in residual makes sense
            MORIS_ASSERT( aStdRes.n_rows() == ( ( this->num_space_dims() + 2 ) * this->num_bases() ) and ( aStdRes.n_cols() == 1 ),
                    "IWG_Compressible_NS_Base::assemble_residual() - Size of residual vector passed in is incorrect." );

            // check residual dof types
            MORIS_ASSERT( check_residual_dof_types( mResidualDofType ),
                    "IWG_Compressible_NS_Base::assemble_residual() - Only pressure or density primitive variables supported for residual assembly." );

            // loop over residual dof types
            for ( uint iResDof = 0; iResDof < mResidualDofType.size(); iResDof++ )
            {
                // get index for residual dof types
                uint tLeaderDofIndex = mSet->get_dof_index_for_type( mResidualDofType( iResDof )( 0 ), mtk::Leader_Follower::LEADER );

                // get residual entry indices for assembly
                uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
                uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

                // get indices, where corresponding dof entries sit in standardized residual
                uint tStdResStartIndex = this->get_assembly_indices( mResidualDofType( iResDof )( 0 ) )( 0 );
                uint tStdResStopIndex  = this->get_assembly_indices( mResidualDofType( iResDof )( 0 ) )( 1 );

                // assemble into set residual
                mSet->get_residual()( 0 )( { tLeaderResStartIndex, tLeaderResStopIndex }, { 0, 0 } ) +=
                        aStdRes( { tStdResStartIndex, tStdResStopIndex }, { 0, 0 } );
            }
        }

        //------------------------------------------------------------------------------

        void
        IWG_Compressible_NS_Base::assemble_jacobian( const Matrix< DDRMat >& aStdJac )
        {
            // check that the size of the passed in residual makes sense
            MORIS_ASSERT( ( aStdJac.n_rows() == ( this->num_space_dims() + 2 ) * this->num_bases() ) and ( aStdJac.n_cols() == ( this->num_space_dims() + 2 ) * this->num_bases() ),
                    "IWG_Compressible_NS_Base::assemble_jacobian() - Size of Jacobian passed in is incorrect." );

            // check residual dof types
            MORIS_ASSERT( check_residual_dof_types( mResidualDofType ),
                    "IWG_Compressible_NS_Base::assemble_jacobian() - Only pressure or density primitive variables supported for jacobian assembly." );

            // check DoF dependencies
            MORIS_ASSERT( check_dof_dependencies( mSet, mResidualDofType, mRequestedLeaderGlobalDofTypes ),
                    "IWG_Compressible_NS_Base::assemble_jacobian() - Set of DoF dependencies not suppported." );

            // loop over residual dof types
            for ( uint iResDof = 0; iResDof < mResidualDofType.size(); iResDof++ )
            {
                // get index for residual dof types
                uint tLeaderDofIndex = mSet->get_dof_index_for_type( mResidualDofType( iResDof )( 0 ), mtk::Leader_Follower::LEADER );

                // get residual entry indices for assembly
                uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
                uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

                // get indices, where corresponding dof entries sit in standardized residual
                uint tStdResStartIndex = this->get_assembly_indices( mResidualDofType( iResDof )( 0 ) )( 0 );
                uint tStdResStopIndex  = this->get_assembly_indices( mResidualDofType( iResDof )( 0 ) )( 1 );

                // loop over dependent dof types
                for ( uint iDepDof = 0; iDepDof < mRequestedLeaderGlobalDofTypes.size(); iDepDof++ )
                {
                    // get index for dependent dof types
                    sint tDepDofIndex =
                            mSet->get_dof_index_for_type( mRequestedLeaderGlobalDofTypes( iDepDof )( 0 ), mtk::Leader_Follower::LEADER );

                    // get dependent variable indices for assembly
                    uint tLeaderDepStartIndex = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDepDofIndex, 0 );
                    uint tLeaderDepStopIndex  = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDepDofIndex, 1 );

                    // get indices, where corresponding dof entries sit in standardized residual
                    uint tStdDepStartIndex = this->get_assembly_indices( mRequestedLeaderGlobalDofTypes( iDepDof )( 0 ) )( 0 );
                    uint tStdDepStopIndex  = this->get_assembly_indices( mRequestedLeaderGlobalDofTypes( iDepDof )( 0 ) )( 1 );

                    // assemble into set jacobian
                    mSet->get_jacobian()( { tLeaderResStartIndex, tLeaderResStopIndex }, { tLeaderDepStartIndex, tLeaderDepStopIndex } ) +=
                            aStdJac( { tStdResStartIndex, tStdResStopIndex }, { tStdDepStartIndex, tStdDepStopIndex } );
                }
            }
        }

        //------------------------------------------------------------------------------

        Matrix< DDSMat >&
        IWG_Compressible_NS_Base::get_assembly_indices( const MSI::Dof_Type aDofType )
        {
            // initialize index vector
            mAssemblyIndices = { { -1 }, { -1 } };

            // get number of bases per DoF type
            uint tNumBases = this->num_bases();

            // get number of spatial dimensions
            uint tNumSpaceDims = this->num_space_dims();

            // check which Dof Type it is, and get corresponding indices
            if ( aDofType == mFirstDof )    // P
            {
                mAssemblyIndices( 0 ) = 0;
                mAssemblyIndices( 1 ) = tNumBases - 1;
            }
            else if ( aDofType == mVectorDof )    // VX
            {
                mAssemblyIndices( 0 ) = tNumBases;
                mAssemblyIndices( 1 ) = ( tNumSpaceDims + 1 ) * tNumBases - 1;
            }
            else if ( aDofType == mLastDof )    // TEMP
            {
                mAssemblyIndices( 0 ) = ( tNumSpaceDims + 1 ) * tNumBases;
                mAssemblyIndices( 1 ) = ( tNumSpaceDims + 2 ) * tNumBases - 1;
            }
            else
            {
                MORIS_ERROR( false,
                        "IWG_Compressible_NS_Base::get_assembly_indices - requested DoF type not part of variable set." );
            }

            // pass on vector with assembly indices
            return mAssemblyIndices;
        }

        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------

        uint
        IWG_Compressible_NS_Base::num_space_dims()
        {
            // get number of spatial dimensions from velocity, momentum, etc. field interpolator
            return mLeaderFIManager->get_field_interpolators_for_type( mVectorDof )->get_number_of_fields();
        }

        //------------------------------------------------------------------------------

        uint
        IWG_Compressible_NS_Base::num_bases()
        {
            // get number of bases from first state var FI
            return mLeaderFIManager->get_field_interpolators_for_type( mFirstDof )->get_number_of_space_time_bases();
        }

        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        IWG_Compressible_NS_Base::Y()
        {
            // check if the variable vectors have already been assembled
            if ( !mYEval )
            {
                return mY;
            }
            mYEval = false;

            // check residual DoF types
            MORIS_ASSERT( check_residual_dof_types( mResidualDofType ),
                    "IWG_Compressible_NS_Base::Y() - check for residual DoF types failed. See Error message above for more info." );

            // get field interpolators
            Field_Interpolator* tFI1 = mLeaderFIManager->get_field_interpolators_for_type( mFirstDof );
            Field_Interpolator* tFI2 = mLeaderFIManager->get_field_interpolators_for_type( mVectorDof );
            Field_Interpolator* tFI3 = mLeaderFIManager->get_field_interpolators_for_type( mLastDof );

            // clang-format off
            // construct Y - vector based on the number of space dims
            switch ( this->num_space_dims() )
            {
                // 2D
                case 2 :
                {
                    mY = {  { tFI1->val()( 0 ) },
                            { tFI2->val()( 0 ) },
                            { tFI2->val()( 1 ) },
                            { tFI3->val()( 0 ) } };
                    break;
                }

                // 3D
                case 3 :
                {
                    mY = {  { tFI1->val()( 0 ) },
                            { tFI2->val()( 0 ) },
                            { tFI2->val()( 1 ) },
                            { tFI2->val()( 2 ) },
                            { tFI3->val()( 0 ) } };
                    break;
                }

                // invalid number of spatial dimensions
                default :
                {
                    MORIS_ERROR( false, "IWG_Compressible_NS_Base::Y() - Number of spatial dimensions must be 2 or 3;" );
                    break;
                }
            }
            // clang-format on

            // check that thermodynamic variables are not zero
            MORIS_ASSERT( tFI1->val()( 0 ) != 0.0 and tFI3->val()( 0 ) != 0.0,
                    "cl_FEM_IWG_Compressible_NS_Base::Y() - Pressure, Density, or Temperature is zero, exiting." );

            // return Y
            return mY;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        IWG_Compressible_NS_Base::dYdt()
        {
            // check if the variable vectors have already been assembled
            if ( !mdYdtEval )
            {
                return mdYdt;
            }
            mdYdtEval = false;

            // check residual DoF types
            MORIS_ASSERT( check_residual_dof_types( mResidualDofType ),
                    "IWG_Compressible_NS_Base::dYdt() - check for residual DoF types failed. See Error message above for more info." );

            // get field interpolators
            Field_Interpolator* tFI1 = mLeaderFIManager->get_field_interpolators_for_type( mFirstDof );
            Field_Interpolator* tFI2 = mLeaderFIManager->get_field_interpolators_for_type( mVectorDof );
            Field_Interpolator* tFI3 = mLeaderFIManager->get_field_interpolators_for_type( mLastDof );

            // clang-format off
            // construct Y - vector based on the number of space dims
            switch ( this->num_space_dims() )
            {
                // 2D
                case 2 :
                {
                    mdYdt = {
                            { tFI1->gradt( 1 )( 0 ) },
                            { tFI2->gradt( 1 )( 0 ) },
                            { tFI2->gradt( 1 )( 1 ) },
                            { tFI3->gradt( 1 )( 0 ) } };
                    break;
                }

                // 3D
                case 3 :
                {
                    mdYdt = {
                            { tFI1->gradt( 1 )( 0 ) },
                            { tFI2->gradt( 1 )( 0 ) },
                            { tFI2->gradt( 1 )( 1 ) },
                            { tFI2->gradt( 1 )( 2 ) },
                            { tFI3->gradt( 1 )( 0 ) } };
                    break;
                }

                // invalid number of spatial dimensions
                default :
                {
                    MORIS_ERROR( false, "IWG_Compressible_NS_Base::dYdt() - Number of spatial dimensions must be 2 or 3;" );
                    break;
                }
            }
            // clang-format on

            // return dYdt
            return mdYdt;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        IWG_Compressible_NS_Base::dYdx( const uint aSpatialDirection )
        {
            // check if the variable vector has already been assembled
            if ( !mdYdxEval )
            {
                return mdYdx( aSpatialDirection );
            }
            mdYdxEval = false;

            // check residual DoF types
            MORIS_ASSERT( check_residual_dof_types( mResidualDofType ),
                    "IWG_Compressible_NS_Base::dYdx() - check for residual DoF types failed. See Error message above for more info." );

            // get field interpolators
            Field_Interpolator* tFI1 = mLeaderFIManager->get_field_interpolators_for_type( mFirstDof );
            Field_Interpolator* tFI2 = mLeaderFIManager->get_field_interpolators_for_type( mVectorDof );
            Field_Interpolator* tFI3 = mLeaderFIManager->get_field_interpolators_for_type( mLastDof );

            // clang-format off
            // construct Y - vector based on the number of space dims
            switch ( this->num_space_dims() )
            {
                // 2D
                case 2 :
                {
                    // set size of cell
                    mdYdx.resize( 2 );

                    // x-derivative
                    mdYdx( 0 ) = {
                            { tFI1->gradx( 1 )( 0 )    },
                            { tFI2->gradx( 1 )( 0, 0 ) },
                            { tFI2->gradx( 1 )( 0, 1 ) },
                            { tFI3->gradx( 1 )( 0 )    } };

                    // y-derivative
                    mdYdx( 1 ) = {
                            { tFI1->gradx( 1 )( 1 )    },
                            { tFI2->gradx( 1 )( 1, 0 ) },
                            { tFI2->gradx( 1 )( 1, 1 ) },
                            { tFI3->gradx( 1 )( 1 )    } };

                    break;
                }

                // 3D
                case 3 :
                {
                    // set size of cell
                    mdYdx.resize( 3 );

                    // x-derivative
                    mdYdx( 0 ) = {
                            { tFI1->gradx( 1 )( 0 )    },
                            { tFI2->gradx( 1 )( 0, 0 ) },
                            { tFI2->gradx( 1 )( 0, 1 ) },
                            { tFI2->gradx( 1 )( 0, 2 ) },
                            { tFI3->gradx( 1 )( 0 )    } };

                    // y-derivative
                    mdYdx( 1 ) = {
                            { tFI1->gradx( 1 )( 1 )    },
                            { tFI2->gradx( 1 )( 1, 0 ) },
                            { tFI2->gradx( 1 )( 1, 1 ) },
                            { tFI2->gradx( 1 )( 1, 2 ) },
                            { tFI3->gradx( 1 )( 1 )    } };

                    // z-derivative
                    mdYdx( 2 ) = {
                            { tFI1->gradx( 1 )( 2 )    },
                            { tFI2->gradx( 1 )( 2, 0 ) },
                            { tFI2->gradx( 1 )( 2, 1 ) },
                            { tFI2->gradx( 1 )( 2, 2 ) },
                            { tFI3->gradx( 1 )( 2 )    } };

                    break;
                }

                // invalid number of spatial dimensions
                default :
                {
                    MORIS_ERROR( false, "IWG_Compressible_NS_Base::dYdx() - Number of spatial dimensions must be 2 or 3;" );
                    break;
                }
            }
            // clang-format on

            // return dYdt
            return mdYdx( aSpatialDirection );
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        IWG_Compressible_NS_Base::d2Ydx2( const uint aI, const uint aJ )
        {
            // convert the two indices into one for condensed tensor
            uint tFlatIndex = convert_index_pair_to_flat( aI, aJ, this->num_space_dims() );

            // check if the variable vectors have already been assembled
            if ( !md2Ydx2Eval )
            {
                return md2Ydx2( tFlatIndex );
            }
            md2Ydx2Eval = false;

            // check residual DoF types
            MORIS_ASSERT( check_residual_dof_types( mResidualDofType ),
                    "IWG_Compressible_NS_Base::d2Ydx2() - check for residual DoF types failed. See Error message above for more info." );

            // get field interpolators
            Field_Interpolator* tFI1 = mLeaderFIManager->get_field_interpolators_for_type( mFirstDof );
            Field_Interpolator* tFI2 = mLeaderFIManager->get_field_interpolators_for_type( mVectorDof );
            Field_Interpolator* tFI3 = mLeaderFIManager->get_field_interpolators_for_type( mLastDof );

            // clang-format off
            // construct Y - vector based on the number of space dims
            switch ( this->num_space_dims() )
            {
                // 2D
                case 2 :
                {
                    // set size of cell
                    md2Ydx2.resize( 3 );

                    // loop over six derivative combinations
                    for ( uint iDeriv = 0; iDeriv < 3; iDeriv++ )
                    {
                        md2Ydx2( iDeriv ) = {
                                { tFI1->gradx( 2 )( iDeriv )    },
                                { tFI2->gradx( 2 )( iDeriv, 0 ) },
                                { tFI2->gradx( 2 )( iDeriv, 1 ) },
                                { tFI3->gradx( 2 )( iDeriv )    } };
                    }
                    break;
                }

                // 3D
                case 3 :
                {
                    // set size of cell
                    md2Ydx2.resize( 6 );

                    // loop over six derivative combinations
                    for ( uint iDeriv = 0; iDeriv < 6; iDeriv++ )
                    {
                        md2Ydx2( iDeriv ) = {
                                { tFI1->gradx( 2 )( iDeriv )    },
                                { tFI2->gradx( 2 )( iDeriv, 0 ) },
                                { tFI2->gradx( 2 )( iDeriv, 1 ) },
                                { tFI2->gradx( 2 )( iDeriv, 2 ) },
                                { tFI3->gradx( 2 )( iDeriv )    } };
                    }
                    break;
                }

                // invalid number of spatial dimensions
                default :
                {
                    MORIS_ERROR( false, "IWG_Compressible_NS_Base::d2Ydx2() - Number of spatial dimensions must be 2 or 3;" );
                    break;
                }
            }
            // clang-format off

            // return value
            return md2Ydx2( tFlatIndex );
        }

        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & IWG_Compressible_NS_Base::W()
        {
            // check if the variable vectors have already been assembled
            if( !mWEval )
            {
                return mW;
            }
            mWEval = false;

            // check residual DoF types
            MORIS_ASSERT( check_residual_dof_types( mResidualDofType  ),
                    "IWG_Compressible_NS_Bulk::W() - check for residual DoF types failed." );

            // get number of state variable fields
            uint tNumStateVars = this->num_space_dims() + 2;

            // initialize W matrix
            mW.set_size( tNumStateVars, tNumStateVars * this->num_bases(), 0.0 );

            // get representative values for the different basis function vectors
            // NOTE: only works under the assumption that all state variable fields are interpolated on the same mesh
            Field_Interpolator * tFI =  mLeaderFIManager->get_field_interpolators_for_type( mFirstDof );
            Matrix< DDRMat > tN = tFI->N()( { 0, 0 }, { 0, this->num_bases() - 1 } );

            // go through residual dof types and assemble the test function matrix
            for ( uint iVar = 0; iVar < tNumStateVars; iVar++ )
            {
                // put test functions into matrices
                mW( { iVar, iVar  }, { iVar * this->num_bases(), ( iVar + 1 ) * this->num_bases() - 1 } ) =
                        tN.matrix_data();
            }

            // return value
            return mW;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & IWG_Compressible_NS_Base::W_trans()
        {
            // check if the variable vectors have already been assembled
            if( !mWtransEval )
            {
                return mWtrans;
            }
            mWtransEval = false;

            // check residual DoF types
            MORIS_ASSERT( check_residual_dof_types( mResidualDofType  ),
                    "IWG_Compressible_NS_Bulk::W() - check for residual DoF types failed." );

            // evaluate transpose
            mWtrans = trans( this->W() );

            // return value
            return mWtrans;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & IWG_Compressible_NS_Base::dWdt()
        {
            // check if the variable vectors have already been assembled
            if( !mdWdtEval )
            {
                return mdWdt;
            }
            mdWdtEval = false;

            // check residual DoF types
            MORIS_ASSERT( check_residual_dof_types( mResidualDofType  ),
                    "IWG_Compressible_NS_Bulk::dWdt() - check for residual DoF types failed." );

            // get number of state variable fields
            uint tNumStateVars = this->num_space_dims() + 2;

            // initialize W matrix
            mdWdt.set_size( tNumStateVars, tNumStateVars * this->num_bases(), 0.0 );

            // get representative values for the different basis function vectors
            // NOTE: only works under the assumption that all state variable fields are interpolated on the same mesh
            Field_Interpolator * tFI =  mLeaderFIManager->get_field_interpolators_for_type( mFirstDof );
            Matrix< DDRMat > tdNdt = tFI->dnNdtn( 1 )( { 0, 0 }, { 0, this->num_bases() - 1 } );

            // go through residual dof types and assemble the test function matrix
            for ( uint iVar = 0; iVar < tNumStateVars; iVar++ )
            {
                // put test functions into matrices
                mdWdt( { iVar, iVar }, { iVar * this->num_bases(), ( iVar + 1 ) * this->num_bases() - 1 } ) =
                        tdNdt.matrix_data();
            }

            // return value
            return mdWdt;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & IWG_Compressible_NS_Base::dWdt_trans()
        {
            // check if the variable vectors have already been assembled
            if( !mdWtransdtEval )
            {
                return mdWtransdt;
            }
            mdWtransdtEval = false;

            // check residual DoF types
            MORIS_ASSERT( check_residual_dof_types( mResidualDofType  ),
                    "IWG_Compressible_NS_Bulk::W() - check for residual DoF types failed." );

            // evaluate transpose
            mdWtransdt = trans( this->dWdt() );

            // return value
            return mdWtransdt;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & IWG_Compressible_NS_Base::dWdx( const uint aSpatialDirection )
        {
            // check if the variable vectors have already been assembled
            if( !mdWdxEval )
            {
                return mdWdx( aSpatialDirection );
            }
            mdWdxEval = false;

            // check residual DoF types
            MORIS_ASSERT( check_residual_dof_types( mResidualDofType  ),
                    "IWG_Compressible_NS_Bulk::dWdx() - check for residual DoF types failed." );

            // get number of state variable fields
            uint tNumStateVars = this->num_space_dims() + 2;

            // create zero mat
            Matrix< DDRMat > tZeroMat( tNumStateVars, tNumStateVars * this->num_bases(), 0.0 );

            // initialize W matrix
            mdWdx.assign( this->num_space_dims(), tZeroMat );

            // get representative values for the different basis function vectors
            // NOTE: only works under the assumption that all state variable fields are interpolated on the same mesh
            Matrix< DDRMat > tdNdx = mLeaderFIManager->get_field_interpolators_for_type( mFirstDof )->dnNdxn( 1 );

            // go through residual dof types and assemble the test function matrix
            for ( uint iVar = 0; iVar < tNumStateVars; iVar++ )
            {
                for ( uint iDim = 0; iDim < this->num_space_dims(); iDim++ )
                {
                    mdWdx( iDim )( { iVar, iVar }, { iVar * this->num_bases(), ( iVar + 1 ) * this->num_bases() - 1 } ) =
                            tdNdx( { iDim, iDim }, { 0, this->num_bases() - 1 } );
                }
            }

            // return value
            return mdWdx( aSpatialDirection );
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & IWG_Compressible_NS_Base::dWdx_trans( const uint aSpatialDirection )
        {
            // check if the variable vectors have already been assembled
            if( !mdWtransdxEval )
            {
                return mdWtransdx( aSpatialDirection );
            }
            mdWtransdxEval = false;

            // check residual DoF types
            MORIS_ASSERT( check_residual_dof_types( mResidualDofType  ),
                    "IWG_Compressible_NS_Bulk::W() - check for residual DoF types failed." );

            // initialize cell
            mdWtransdx.resize( this->num_space_dims() );

            // evaluate transpose
            for ( uint iDim = 0; iDim < this->num_space_dims(); iDim++ )
            {
                mdWtransdx( iDim ) = trans( this->dWdx( iDim ) );
            }

            // return value
            return mdWtransdx( aSpatialDirection );
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & IWG_Compressible_NS_Base::d2Wdx2( const uint aI, const uint aJ )
        {
            // convert the two indices into one for condensed tensor
            uint tFlatIndex = convert_index_pair_to_flat( aI, aJ, this->num_space_dims() );

            // check if the variable vectors have already been assembled
            if( !md2Wdx2Eval )
            {
                return md2Wdx2( tFlatIndex );
            }
            md2Wdx2Eval = false;

            // check residual DoF types
            MORIS_ASSERT( check_residual_dof_types( mResidualDofType  ),
                    "IWG_Compressible_NS_Bulk::d2Wdx2() - check for residual DoF types failed." );

            // get number of state variable fields
            uint tNumStateVars = this->num_space_dims() + 2;
            uint tNumSecondSpaceDerivs = 3 * this->num_space_dims() - 3;

            // create zero mat
            Matrix< DDRMat > tZeroMat( tNumStateVars, tNumStateVars * this->num_bases(), 0.0 );

            // initialize W matrix
            md2Wdx2.assign( tNumSecondSpaceDerivs, tZeroMat );

            // get representative values for the different basis function vectors
            // NOTE: only works under the assumption that all state variable fields are interpolated on the same mesh
            Matrix< DDRMat > td2Ndx2 = mLeaderFIManager->get_field_interpolators_for_type( mFirstDof )->dnNdxn( 2 );

            // go through residual dof types and assemble the test function matrix
            for ( uint iVar = 0; iVar < tNumStateVars; iVar++ )
            {
                for ( uint iRow = 0; iRow < tNumSecondSpaceDerivs; iRow++ )
                {
                    md2Wdx2( iRow )( { iVar, iVar }, { iVar * this->num_bases(), ( iVar + 1 ) * this->num_bases() - 1 } ) =
                            td2Ndx2( { iRow, iRow }, { 0, this->num_bases() - 1 } );
                }
            }

            // return value
            return md2Wdx2( tFlatIndex );
        }

        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & IWG_Compressible_NS_Base::A( const uint aK )
        {
            // check that indices are not out of bounds
            MORIS_ASSERT( aK <= this->num_space_dims(),
                    "IWG_Compressible_NS_Base::A() - index  out of bounds." );

            // check if the variable vectors have already been assembled
            if( !mAEval )
            {
                return mA( aK );
            }

            // set the eval flag
            mAEval = false;

            // get the material and constitutive models
            std::shared_ptr< Material_Model > tMM = mLeaderMM( static_cast< uint >( IWG_Material_Type::FLUID_MM ) );
            std::shared_ptr< Constitutive_Model > tCM = mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_CM ) );

            // evaluate A-matrices and store them
            eval_A( tMM, tCM, mLeaderFIManager, mResidualDofType, mA );

            // return requested value
            return mA( aK );
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & IWG_Compressible_NS_Base::K(
                const uint aI,
                const uint aJ )
        {
            // check that indices are not out of bounds
            MORIS_ASSERT( ( aI < this->num_space_dims() ) and ( aJ < this->num_space_dims() ),
                    "IWG_Compressible_NS_Base::K() - indices out of bounds." );

            // check if K matrices have already been evaluated
            if ( !mKEval )
            {
                return mK( aI )( aJ );
            }

            // set the eval flag
            mKEval = false;

            // get the viscosity
            std::shared_ptr< Property > tPropDynamicViscosity = mLeaderProp( static_cast< uint >( IWG_Property_Type::DYNAMIC_VISCOSITY ) );
            std::shared_ptr< Property > tPropThermalConductivity = mLeaderProp( static_cast< uint >( IWG_Property_Type::THERMAL_CONDUCTIVITY ) );

            // eval K matrices and store them
            eval_K( tPropDynamicViscosity, tPropThermalConductivity, mLeaderFIManager, mK );

            // return requested K matrix
            return mK( aI )( aJ );
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & IWG_Compressible_NS_Base::Kiji( const uint aJ )
        {
            // check that indices are not out of bounds
            MORIS_ASSERT( aJ < this->num_space_dims(),
                    "IWG_Compressible_NS_Base::Kiji() - index out of bounds." );

            // check if Kiji matrices have already been evaluated
            if ( !mKijiEval )
            {
                return mKiji( aJ );
            }

            // set the eval flag
            mKijiEval = false;

            // get the viscosity
            std::shared_ptr< Property > tPropDynamicViscosity = mLeaderProp( static_cast< uint >( IWG_Property_Type::DYNAMIC_VISCOSITY ) );
            std::shared_ptr< Property > tPropThermalConductivity = mLeaderProp( static_cast< uint >( IWG_Property_Type::THERMAL_CONDUCTIVITY ) );

            // eval spatial derivatives of K matrices and store them
            eval_dKijdxi( tPropDynamicViscosity, tPropThermalConductivity, mLeaderFIManager, mKiji );

            // return requested Kiji matrix
            return mKiji( aJ );
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & IWG_Compressible_NS_Base::C()
        {
            // check if matrix is already evaluated
            if ( !mCEval )
            {
                return mC;
            }

            // set the eval flag
            mCEval = false;

            // get the body heat load
            std::shared_ptr< Property > tPropBodyHeatLoad = mLeaderProp( static_cast< uint >( IWG_Property_Type::BODY_HEAT_LOAD ) );
            real tQ = 0.0;
            if ( tPropBodyHeatLoad != nullptr )
            {
                tQ = tPropBodyHeatLoad->val()( 0 );
            }

            // FIXME: body force not considered yet
            // std::shared_ptr< Property > tPropBodyForce = mLeaderProp( static_cast< uint >( IWG_Property_Type::BODY_FORCE ) );

            // get the material model
            std::shared_ptr< Material_Model > tMM = mLeaderMM( static_cast< uint >( IWG_Material_Type::FLUID_MM ) );

            // get number of state variables
            uint tNumStateVars = this->num_space_dims() + 2;

            // define coefficient matrix
            mC.set_size( tNumStateVars, tNumStateVars, 0.0 );
            mC( tNumStateVars - 1, tNumStateVars - 1 ) = -1.0 * tQ / tMM->temperature()( 0 );

            // return coefficient matrix
            return mC;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & IWG_Compressible_NS_Base::dCdY_VR( const Matrix< DDRMat > aVR )
        {
            // get number of state variables
            uint tNumStateVars = this->num_space_dims() + 2;

            // check input vector
            MORIS_ASSERT( aVR.length() == tNumStateVars,
                    "IWG_Compressible_NS_Base::dCdY_VR() - pre-multiplication vector of incorrect size." );

            // get the body heat load
            std::shared_ptr< Property > tPropBodyHeatLoad = mLeaderProp( static_cast< uint >( IWG_Property_Type::BODY_HEAT_LOAD ) );
            real tQ = 0.0;
            if ( tPropBodyHeatLoad != nullptr )
            {
                tQ = tPropBodyHeatLoad->val()( 0 );
            }

            // FIXME: body force not considered yet
            // std::shared_ptr< Property > tPropBodyForce = mLeaderProp( static_cast< uint >( IWG_Property_Type::BODY_FORCE ) );

            // get the material model temperature
            std::shared_ptr< Material_Model > tMM = mLeaderMM( static_cast< uint >( IWG_Material_Type::FLUID_MM ) );
            real tTemp = tMM->temperature()( 0 );

            // get last entry of Vr
            mdCdYVR.set_size( tNumStateVars, tNumStateVars, 0.0 );
            mdCdYVR( tNumStateVars - 1, tNumStateVars - 1 ) = aVR( tNumStateVars - 1 ) * tQ / tTemp / tTemp;

            // return
            return mdCdYVR;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & IWG_Compressible_NS_Base::dCdY( const uint aYind )
        {
            // check that index are not out of bounds
            MORIS_ASSERT( aYind < this->num_space_dims() + 2,
                    "IWG_Compressible_NS_Base::dCdY() - state variable index out of bounds." );

            // check if matrix is already evaluated
            if ( !mdCdYEval )
            {
                return mdCdY( aYind );
            }

            // set the eval flag
            mdCdYEval = false;

            // get the body heat load
            std::shared_ptr< Property > tPropBodyHeatLoad = mLeaderProp( static_cast< uint >( IWG_Property_Type::BODY_HEAT_LOAD ) );
            real tQ = 0.0;
            if ( tPropBodyHeatLoad != nullptr )
            {
                tQ = tPropBodyHeatLoad->val()( 0 );
            }

            // FIXME: body force not considered yet
            // std::shared_ptr< Property > tPropBodyForce = mLeaderProp( static_cast< uint >( IWG_Property_Type::BODY_FORCE ) );

            // get the material model temperature
            std::shared_ptr< Material_Model > tMM = mLeaderMM( static_cast< uint >( IWG_Material_Type::FLUID_MM ) );
            real tTemp = tMM->temperature()( 0 );

            // get number of state variables
            uint tNumStateVars = this->num_space_dims() + 2;

            // define coefficient matrix
            Matrix< DDRMat > tZeroMat( tNumStateVars, tNumStateVars, 0.0 );
            mdCdY.resize( tNumStateVars, tZeroMat );
            mdCdY( tNumStateVars - 1 )( tNumStateVars - 1, tNumStateVars - 1 ) = tQ / tTemp / tTemp;

            // return coefficient matrix state var derivative
            return mdCdY( aYind );
        }

        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------

        const Matrix< DDRMat > IWG_Compressible_NS_Base::get_elemental_index_vector()
        {
            // get total number of bases
            real tTotNumBases = ( this->num_space_dims() + 2 ) * this->num_bases();

            // initialize vector
            Matrix< DDRMat > tIndexVec( tTotNumBases, 1, 0.0 );

            // loop over entries and fill
            for ( uint iRow = 0; iRow < tTotNumBases; iRow++ )
            {
                tIndexVec( iRow ) = (real) iRow;
            }

            // return index vector
            return tIndexVec;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > IWG_Compressible_NS_Base::get_elemental_index_matrix()
        {
            // get total number of bases
            real tTotNumBases = ( this->num_space_dims() + 2 ) * this->num_bases();

            // initialize vector
            Matrix< DDRMat > tIndexMat( tTotNumBases, tTotNumBases, 0.0 );

            // loop over entries and fill
            for ( uint iRow = 0; iRow < tTotNumBases; iRow++ )
            {
                for ( uint iCol = 0; iCol < tTotNumBases; iCol++ )
                {
                    tIndexMat( iRow, iCol ) = (real) iRow * 1000.0 + (real) iCol;

                    // for string to double conversion, use std::stod
                }
            }

            // return index matrix
            return tIndexMat;
        }

        //------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */

