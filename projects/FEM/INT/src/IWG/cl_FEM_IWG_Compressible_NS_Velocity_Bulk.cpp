/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Compressible_NS_Velocity_Bulk.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG_Compressible_NS_Velocity_Bulk.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        IWG_Compressible_NS_Velocity_Bulk::IWG_Compressible_NS_Velocity_Bulk()
        {
            init_property( "BodyForce", IWG_Property_Type::BODY_FORCE );
            init_constitutive_model( "Fluid", IWG_Constitutive_Type::FLUID );
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Velocity_Bulk::compute_residual( real aWStar )
        {
            // check leader field interpolators
#ifdef MORIS_HAVE_DEBUG
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type (here velocity), indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get the velocity FI
            Field_Interpolator *tFIVelocity = get_leader_fi_manager()->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get Density FI
            Field_Interpolator *tFIDensity = get_leader_fi_manager()->get_field_interpolators_for_type( mDofDensity );

            // get the mass body force property
            std::shared_ptr< Property > tPropBodyForce = get_leader_property(IWG_Property_Type::BODY_FORCE);

            // get the compressible fluid constitutive model
            std::shared_ptr< Constitutive_Model > tCMFluid = get_leader_constitutive_model(IWG_Constitutive_Type::FLUID);

            // build dyadic product of velocity vectors
            Matrix< DDRMat > tUiUj;
            this->compute_uiuj( tUiUj );

            // compute the residual
            mSet->get_residual()( 0 )( { tLeaderResStartIndex, tLeaderResStopIndex }, { 0, 0 } ) += aWStar * ( tFIDensity->gradt( 1 )( 0 ) * trans( tFIVelocity->N() ) * tFIVelocity->val() + tFIDensity->val()( 0 ) * trans( tFIVelocity->N() ) * trans( tFIVelocity->gradt( 1 ) ) + ( -1.0 ) * tFIDensity->val()( 0 ) * trans( tCMFluid->testStrain() ) * this->MultipMat() * tUiUj + trans( tCMFluid->testStrain() ) * this->MultipMat() * tCMFluid->flux( CM_Function_Type::MECHANICAL ) );

            // if there is a body force
            if ( tPropBodyForce != nullptr )
            {
                // add gravity to residual weak form
                mSet->get_residual()( 0 )( { tLeaderResStartIndex, tLeaderResStopIndex }, { 0, 0 } ) -= aWStar * ( tFIDensity->val()( 0 ) * trans( tFIVelocity->N() ) * tPropBodyForce->val() );
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Compressible_NS_Velocity_Bulk::compute_residual - Residual contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Velocity_Bulk::compute_jacobian( real aWStar )
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
            Field_Interpolator *tFIVelocity = get_leader_fi_manager()->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );
            Field_Interpolator *tFIDensity  = get_leader_fi_manager()->get_field_interpolators_for_type( mDofDensity );

            // get the mass body force property
            std::shared_ptr< Property > tPropBodyForce = get_leader_property(IWG_Property_Type::BODY_FORCE);

            // get the compressible fluid constitutive model
            std::shared_ptr< Constitutive_Model > tCMFluid = get_leader_constitutive_model(IWG_Constitutive_Type::FLUID);

            // compute the jacobian for dof dependencies
            uint tNumDofDependencies = get_requested_leader_dof_types().size();
            for ( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the treated dof type
                Vector< MSI::Dof_Type > const &tDofType = get_requested_leader_dof_types()( iDOF );

                // get the index for dof type, indices for assembly
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                uint tLeaderDepStartIndex = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 0 );
                uint tLeaderDepStopIndex  = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 1 );

                // if fluid CM depends on dof type
                if ( tCMFluid->check_dof_dependency( tDofType ) )
                {
                    // add contribution
                    mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) += aWStar * ( trans( tCMFluid->testStrain() ) * this->MultipMat() * tCMFluid->dFluxdDOF( tDofType, CM_Function_Type::MECHANICAL ) );
                }

                // if dof type is velocity, add diagonal term (velocity-velocity DoF types)
                if ( tDofType( 0 ) == mDofDensity )
                {
                    // build dyadic product of velocity vectors
                    Matrix< DDRMat > tUiUj;
                    this->compute_uiuj( tUiUj );

                    // add contribution
                    mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) += aWStar * ( trans( tFIVelocity->N() ) * trans( tFIVelocity->gradt( 1 ) ) * tFIDensity->N() + trans( tFIVelocity->N() ) * tFIVelocity->val() * tFIDensity->dnNdtn( 1 ) + ( -1.0 ) * trans( tCMFluid->testStrain() ) * this->MultipMat() * tUiUj * tFIDensity->N() );

                    // if a body force is present
                    if ( tPropBodyForce != nullptr )
                    {
                        // compute the jacobian contribution
                        mSet->get_jacobian()(
                                { tLeaderResStartIndex, tLeaderResStopIndex },
                                { tLeaderDepStartIndex, tLeaderDepStopIndex } ) += aWStar * ( ( -1.0 ) * trans( tFIVelocity->N() ) * tPropBodyForce->val() * tFIDensity->N() );
                    }
                }

                // if dof type is velocity, add diagonal term (velocity-velocity DoF types)
                if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    // build derivative of dyadic product of velocity vectors wrt to dofs
                    Matrix< DDRMat > tdUiUjdDOF;
                    this->compute_duiujdDOF( tdUiUjdDOF );

                    // build multi-D shape function matrix of time derivatives for velocity
                    Matrix< DDRMat > tdnNveldtn;
                    this->compute_dnNdtn( tdnNveldtn );

                    // add contribution
                    mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) += aWStar * ( tFIDensity->gradt( 1 )( 0 ) * trans( tFIVelocity->N() ) * tFIVelocity->N() + tFIDensity->val()( 0 ) * trans( tFIVelocity->N() ) * tdnNveldtn + ( -1.0 ) * tFIDensity->val()( 0 ) * trans( tCMFluid->testStrain() ) * this->MultipMat() * tdUiUjdDOF );
                }

                // if a body force is present
                if ( tPropBodyForce != nullptr )
                {
                    // if the body force depends on the dof type -> indirect dependency
                    if ( tPropBodyForce->check_dof_dependency( tDofType ) )
                    {
                        // compute the jacobian contribution
                        mSet->get_jacobian()(
                                { tLeaderResStartIndex, tLeaderResStopIndex },
                                { tLeaderDepStartIndex, tLeaderDepStopIndex } ) += aWStar * ( ( -1.0 ) * tFIDensity->val()( 0 ) * trans( tFIVelocity->N() ) * tPropBodyForce->dPropdDOF( tDofType ) );
                    }
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ),
                    "IWG_Compressible_NS_Velocity_Bulk::compute_jacobian - Jacobian contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Velocity_Bulk::compute_jacobian_and_residual( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Compressible_NS_Velocity_Bulk::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Velocity_Bulk::compute_dRdp( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Compressible_NS_Velocity_Bulk::compute_dRdp - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Velocity_Bulk::compute_residual_strong_form(
                Matrix< DDRMat > &aRM,
                real             &aRC )
        {
            MORIS_ERROR( false, "IWG_Compressible_NS_Velocity_Bulk::compute_residual_strong_form - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Velocity_Bulk::compute_jacobian_strong_form(
                Vector< MSI::Dof_Type > aDofTypes,
                Matrix< DDRMat >       &aJM,
                Matrix< DDRMat >       &aJC )
        {
            MORIS_ERROR( false, "IWG_Compressible_NS_Velocity_Bulk::compute_jacobian_strong_form - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Velocity_Bulk::compute_uiuj( Matrix< DDRMat > &auiuj )
        {
            // get the velocity vector
            Field_Interpolator *tFIVelocity = get_leader_fi_manager()->get_field_interpolators_for_type( mDofVelocity );

            // get the velocity vector
            Matrix< DDRMat > tVelocityVec = tFIVelocity->val();

            // assembly into flattened tensor
            // for 2D
            if ( tFIVelocity->get_number_of_fields() == 2 )
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

        void IWG_Compressible_NS_Velocity_Bulk::compute_duiujdDOF( Matrix< DDRMat > &aduiujdDOF )
        {
            // get the velocity vector
            Field_Interpolator *tFIVelocity = get_leader_fi_manager()->get_field_interpolators_for_type( mDofVelocity );

            // get the velocity vector and shape functions
            Matrix< DDRMat > tUvec = tFIVelocity->val();
            Matrix< DDRMat > tNmat = tFIVelocity->N();

            // get number of bases
            uint tNumBases = get_leader_fi_manager()->get_field_interpolators_for_type( mDofVelocity )->get_number_of_space_time_bases();

            // assembly
            // for 2D
            if ( tFIVelocity->get_number_of_fields() == 2 )
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

                aduiujdDOF( { 5, 5 }, { 0, tNumBases - 1 } )             = tUvec( 1 ) * tNmat( { 0, 0 }, { 0, tNumBases - 1 } );
                aduiujdDOF( { 5, 5 }, { tNumBases, 2 * tNumBases - 1 } ) = tUvec( 0 ) * tNmat( { 0, 0 }, { 0, tNumBases - 1 } );
            }
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > &IWG_Compressible_NS_Velocity_Bulk::MultipMat()
        {
            // check if multiplication matrix exists
            if ( mMultipMatIsBuild )
            {
                return mMultipMat;
            }

            // build multiplication matrix otherwise
            else
            {
                // build multiplication matrix
                // for 2D
                if ( get_leader_fi_manager()->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) )->get_number_of_fields() == 2 )
                {
                    mMultipMat = {
                        { 1.0, 0.0, 0.0 },
                        { 0.0, 1.0, 0.0 },
                        { 0.0, 0.0, 2.0 }
                    };
                }
                // for 3D
                else
                {
                    mMultipMat = {
                        { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                        { 0.0, 1.0, 0.0, 0.0, 0.0, 0.0 },
                        { 0.0, 0.0, 1.0, 0.0, 0.0, 0.0 },
                        { 0.0, 0.0, 0.0, 2.0, 0.0, 0.0 },
                        { 0.0, 0.0, 0.0, 0.0, 2.0, 0.0 },
                        { 0.0, 0.0, 0.0, 0.0, 0.0, 2.0 }
                    };
                }

                // FIXME: for unit testing mMultipMatIsBuild needs to be reset
                // set evaluation flag
                // mMultipMatIsBuild = true;

                // return matrix
                return mMultipMat;
            }
        }

        //------------------------------------------------------------------------------

        // FIXME provided directly by the field interpolator?
        void IWG_Compressible_NS_Velocity_Bulk::compute_dnNdtn(
                Matrix< DDRMat > &adnNdtn )
        {
            // get the residual dof type FI (here velocity)
            Field_Interpolator *tVelocityFI = get_leader_fi_manager()->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // init size for dnNdtn
            uint tNumRowt = tVelocityFI->get_number_of_fields();
            uint tNumColt = tVelocityFI->dnNdtn( 1 ).n_cols();
            adnNdtn.set_size( tNumRowt, tNumRowt * tNumColt, 0.0 );

            // loop over the fields
            for ( uint iField = 0; iField < tNumRowt; iField++ )
            {
                // fill the matrix for each dimension
                adnNdtn( { iField, iField }, { iField * tNumColt, ( iField + 1 ) * tNumColt - 1 } ) =
                        tVelocityFI->dnNdtn( 1 ).matrix_data();
            }
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
