/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Compressible_NS_Advective_Momentum_Flux_Boundary.cpp
 *
 */

#include "cl_FEM_IWG_Compressible_NS_Advective_Momentum_Flux_Boundary.hpp"

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_norm.hpp"
#include "fn_eye.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        IWG_Compressible_NS_Advective_Momentum_Flux_Boundary::IWG_Compressible_NS_Advective_Momentum_Flux_Boundary()
        {
            // set size for the constitutive model pointer cell
            mLeaderCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "Fluid" ] = static_cast< uint >( IWG_Constitutive_Type::FLUID );
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Advective_Momentum_Flux_Boundary::compute_residual( real aWStar )
        {
            // check leader field interpolators
#ifdef MORIS_HAVE_DEBUG
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type (here velocity), indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get the field interpolators
            Field_Interpolator * tDensityFI = mLeaderFIManager->get_field_interpolators_for_type( mDofDensity );
            Field_Interpolator * tVelocityFI = mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // build dyadic product of velocity vectors
            Matrix< DDRMat > tUiUj;
            this->compute_uiuj( tUiUj );

            // build flattened normal
            Matrix< DDRMat > tNormalMatrix;
            this->compute_normal_matrix( tNormalMatrix );

            // compute the residual weak form
            mSet->get_residual()( 0 )( { tLeaderResStartIndex, tLeaderResStopIndex }, { 0, 0 } ) += aWStar * (
                    tDensityFI->val()( 0 ) * tVelocityFI->N_trans() * tNormalMatrix * tUiUj );

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Compressible_NS_Advective_Momentum_Flux_Boundary::compute_residual - Residual contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------
        void IWG_Compressible_NS_Advective_Momentum_Flux_Boundary::compute_jacobian( real aWStar )
        {
            // check leader field interpolators
#ifdef MORIS_HAVE_DEBUG
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type (here velocity), indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get the FIs
            Field_Interpolator * tFIDensity =  mLeaderFIManager->get_field_interpolators_for_type( mDofDensity );
            Field_Interpolator * tFIVelocity =  mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // compute the jacobian for dof dependencies
            uint tNumDofDependencies = mRequestedLeaderGlobalDofTypes.size();
            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the treated dof type
                Vector< MSI::Dof_Type > & tDofType = mRequestedLeaderGlobalDofTypes( iDOF );

                // get the index for dof type, indices for assembly
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                uint tLeaderDepStartIndex = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 0 );
                uint tLeaderDepStopIndex  = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 1 );

                // if dof type is velocity, add diagonal term (velocity-velocity DoF types)
                if( tDofType( 0 ) == mDofDensity )
                {
                    // build dyadic product of velocity vectors
                    Matrix< DDRMat > tUiUj;
                    this->compute_uiuj( tUiUj );

                    // build flattened normal
                    Matrix< DDRMat > tNormalMatrix;
                    this->compute_normal_matrix( tNormalMatrix );

                    // add contribution
                    mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) += aWStar * (
                                    tFIVelocity->N_trans() * tNormalMatrix * tUiUj * tFIDensity->N() );
                }

                // if dof type is velocity, add diagonal term (velocity-velocity DoF types)
                if( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    // build derivative of dyadic product of velocity vectors wrt to dofs
                    Matrix< DDRMat > tdUiUjdDOF;
                    this->compute_duiujdDOF( tdUiUjdDOF );

                    // build flattened normal
                    Matrix< DDRMat > tNormalMatrix;
                    this->compute_normal_matrix( tNormalMatrix );

                    // add contribution
                    mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) += aWStar * (
                                    tFIDensity->val()( 0 ) * tFIVelocity->N_trans() * tNormalMatrix * tdUiUjdDOF );
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ) ,
                    "IWG_Compressible_NS_Advective_Momentum_Flux_Boundary::compute_jacobian - Jacobian contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Advective_Momentum_Flux_Boundary::compute_jacobian_and_residual( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Compressible_NS_Advective_Momentum_Flux_Boundary::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Advective_Momentum_Flux_Boundary::compute_dRdp( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Compressible_NS_Advective_Momentum_Flux_Boundary::compute_dRdp - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Advective_Momentum_Flux_Boundary::compute_residual_strong_form(
                Matrix< DDRMat > & aRM,
                real             & aRC )
        {
            MORIS_ERROR( false, "IWG_Compressible_NS_Advective_Momentum_Flux_Boundary::compute_residual_strong_form - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Advective_Momentum_Flux_Boundary::compute_jacobian_strong_form(
                Vector< MSI::Dof_Type >   aDofTypes,
                Matrix< DDRMat >             & aJM,
                Matrix< DDRMat >             & aJC )
        {
            MORIS_ERROR( false, "IWG_Compressible_NS_Advective_Mass_Flux_Boundary::compute_jacobian_strong_form - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Advective_Momentum_Flux_Boundary::compute_uiuj(Matrix< DDRMat > & auiuj)
        {
            // get the velocity vector
            Field_Interpolator * tFIVelocity =  mLeaderFIManager->get_field_interpolators_for_type( mDofVelocity );

            // get the velocity vector
            Matrix< DDRMat > tVelocityVec = tFIVelocity->val();

            // assembly into flattened tensor
            // for 2D
            if( tFIVelocity->get_number_of_fields() == 2 )
            {
                auiuj.set_size( 3, 1, 0.0 );

                auiuj( 0 ) = std::pow( tVelocityVec( 0 ), 2.0 );
                auiuj( 1 ) = std::pow( tVelocityVec( 1 ), 2.0 );
                auiuj( 2 ) = tVelocityVec( 0 ) * tVelocityVec( 1 );

            }
            // for 3D
            else
            {
                auiuj.set_size( 6, 1, 0.0 );

                auiuj( 0 ) = std::pow( tVelocityVec( 0 ), 2.0 );
                auiuj( 1 ) = std::pow( tVelocityVec( 1 ), 2.0 );
                auiuj( 2 ) = std::pow( tVelocityVec( 2 ), 2.0 );
                auiuj( 3 ) = tVelocityVec( 1 ) * tVelocityVec( 2 );
                auiuj( 4 ) = tVelocityVec( 0 ) * tVelocityVec( 2 );
                auiuj( 5 ) = tVelocityVec( 0 ) * tVelocityVec( 1 );
            }
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Advective_Momentum_Flux_Boundary::compute_duiujdDOF(Matrix< DDRMat > & aduiujdDOF)
        {
            // get the velocity vector
            Field_Interpolator * tFIVelocity =  mLeaderFIManager->get_field_interpolators_for_type( mDofVelocity );

            // get the velocity vector and shape functions
            Matrix< DDRMat > tUvec = tFIVelocity->val();
            Matrix< DDRMat > tNmat = tFIVelocity->N();

            // get number of bases
            uint tNumBases = mLeaderFIManager->get_field_interpolators_for_type( mDofVelocity )->get_number_of_space_time_bases();

            // assembly
            // for 2D
            if( tFIVelocity->get_number_of_fields() == 2 )
            {
                // initialize
                aduiujdDOF.set_size( 3, tNumBases * 2, 0.0 );

                // fill
                aduiujdDOF( { 0, 0 }, { 0, tNumBases - 1 } )             = 2.0 * tUvec( 0 ) * tNmat( { 0, 0 }, { 0, tNumBases - 1 } );
                aduiujdDOF( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } ) = 2.0 * tUvec( 1 ) * tNmat( { 0, 0 }, { 0, tNumBases - 1 } );
                aduiujdDOF( { 2, 2 }, { 0, tNumBases - 1 } )             = tUvec( 1 ) * tNmat( { 0, 0 }, { 0, tNumBases - 1 } );
                aduiujdDOF( { 2, 2 }, { tNumBases, 2 * tNumBases - 1 } ) = tUvec( 0 ) * tNmat( { 0, 0 }, { 0, tNumBases - 1 } );
            }
            // for 3D
            else
            {
                // initialize
                aduiujdDOF.set_size( 6, tNumBases * 3, 0.0 );

                // fill
                aduiujdDOF( { 0, 0 }, { 0, tNumBases - 1 } )                 = 2.0 * tUvec( 0 ) * tNmat( { 0, 0 }, { 0, tNumBases - 1 } );
                aduiujdDOF( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } )     = 2.0 * tUvec( 1 ) * tNmat( { 0, 0 }, { 0, tNumBases - 1 } );
                aduiujdDOF( { 2, 2 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = 2.0 * tUvec( 2 ) * tNmat( { 0, 0 }, { 0, tNumBases - 1 } );

                aduiujdDOF( { 3, 3 }, { tNumBases, 2 * tNumBases - 1 } )     = tUvec( 2 ) * tNmat( { 0, 0 }, { 0, tNumBases - 1 } );
                aduiujdDOF( { 3, 3 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = tUvec( 1 ) * tNmat( { 0, 0 }, { 0, tNumBases - 1 } );

                aduiujdDOF( { 4, 4 }, { 0, tNumBases - 1 } )                 = tUvec( 2 ) * tNmat( { 0, 0 }, { 0, tNumBases - 1 } );
                aduiujdDOF( { 4, 4 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = tUvec( 0 ) * tNmat( { 0, 0 }, { 0, tNumBases - 1 } );

                aduiujdDOF( { 5, 5 }, { 0, tNumBases - 1 } )                 = tUvec( 1 ) * tNmat( { 0, 0 }, { 0, tNumBases - 1 } );
                aduiujdDOF( { 5, 5 }, { tNumBases, 2 * tNumBases - 1 } )     = tUvec( 0 ) * tNmat( { 0, 0 }, { 0, tNumBases - 1 } );
            }
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Advective_Momentum_Flux_Boundary::compute_normal_matrix( Matrix< DDRMat > & aNormalMatrix )
        {
            // clang-format off
            // assembly into matrix
            // for 2D
            if( mLeaderFIManager->get_field_interpolators_for_type( mDofVelocity )->get_number_of_fields() == 2 )
            {
                aNormalMatrix = {
                        { mNormal( 0 ),         0.0 , mNormal( 1 ) },
                        {         0.0 , mNormal( 1 ), mNormal( 0 ) } };
            }
            // for 3D
            else
            {
                aNormalMatrix = {
                        { mNormal( 0 ),         0.0 ,         0.0 ,         0.0 , mNormal( 2 ), mNormal( 1 ) },
                        {         0.0 , mNormal( 1 ),         0.0 , mNormal( 2 ),         0.0 , mNormal( 0 ) },
                        {         0.0 ,         0.0 , mNormal( 2 ), mNormal( 1 ), mNormal( 0 ),         0.0  } };
            }
            // clang-format on
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

