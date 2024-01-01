/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Compressible_NS_Bulk.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG_Compressible_NS_Bulk.hpp"
#include "fn_FEM_IWG_Compressible_NS.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"

// debug - output to hdf5
#include "paths.hpp"
#include "HDF5_Tools.hpp"
#include "fn_FEM_Check.hpp"

#include "fn_sqrtmat.hpp"
//#include "../test/FEM_Test_Proxy/fn_FEM_Convert_Dimensions.cpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        IWG_Compressible_NS_Bulk::IWG_Compressible_NS_Bulk()
        {
            // set size for the property pointer cell
            mLeaderProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "DynamicViscosity" ]     = static_cast< uint >( IWG_Property_Type::DYNAMIC_VISCOSITY );
            mPropertyMap[ "ThermalConductivity" ]  = static_cast< uint >( IWG_Property_Type::THERMAL_CONDUCTIVITY );
            mPropertyMap[ "BodyForce" ]            = static_cast< uint >( IWG_Property_Type::BODY_FORCE );
            mPropertyMap[ "BodyHeatLoad" ]         = static_cast< uint >( IWG_Property_Type::BODY_HEAT_LOAD );

            // set size for the material model pointer cell
            mLeaderMM.resize( static_cast< uint >( IWG_Material_Type::MAX_ENUM ), nullptr );

            // populate the material map
            mMaterialMap[ "FluidMM" ] = static_cast< uint >( IWG_Material_Type::FLUID_MM );

            // set size for the constitutive model pointer cell
            mLeaderCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "FluidCM" ] = static_cast< uint >( IWG_Constitutive_Type::FLUID_CM );

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the stabilization parameter map
            mStabilizationMap[ "GLS" ] = static_cast< uint >( IWG_Stabilization_Type::GLS );
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Bulk::reset_child_eval_flags()
        {
            // reset eval flags specific to this child IWG
            mLYEval = true;
            mLWEval = true;
            mLDofYEval = true;

            mGLSTestFuncEval = true;

            mA0invEval = true;
            mdA0invdYEval = true;

            mGEval = true;

            mMEval = true;
            mMinvEval = true;
            mSqrtMinvEval = true;
            mdMdYEval = true;

            mTauEval = true;
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Bulk::compute_residual( real aWStar )
        {
            // check leader field interpolators
#ifdef MORIS_HAVE_DEBUG
            this->check_field_interpolators();
#endif
            // check residual dof types
            MORIS_ASSERT( check_residual_dof_types( mResidualDofType  ),
                    "IWG_Compressible_NS_Bulk::compute_jacobian() - Only pressure or density primitive variables supported for now." );

            // get number of space dimensions
            uint tNumSpaceDims = this->num_space_dims();

            // get total number of DoFs on Comp Flow Element
            uint tNumTotalBases = ( tNumSpaceDims + 2 ) * this->num_bases();

            // construct temporary Vector for residual
            Matrix< DDRMat > tTempRes( tNumTotalBases, 1, 0.0 );
            auto tRes = tTempRes( { 0, tNumTotalBases - 1 }, { 0, 0 } );

            // A0 matrix contribution
            tRes += aWStar * this->W_trans() * this->A( 0 ) * this->dYdt();

            // loop over A-Matrices
            for ( uint iA = 1; iA < tNumSpaceDims + 1; iA++ )
            {
                // compute residual
                tRes += aWStar * this->W_trans() * this->A( iA ) * this->dYdx( iA - 1 );
            }

            // loop over K-matrices
            for ( uint iDim = 0; iDim < tNumSpaceDims; iDim++ )
            {
                for ( uint jDim = 0; jDim < tNumSpaceDims; jDim++ )
                {
                    // compute residual
                    tRes += aWStar * this->dWdx_trans( iDim ) * this->K( iDim, jDim ) * this->dYdx( jDim );
                }
            }

            // contribution from body loads
            tRes += aWStar * this->W_trans() * this->C() * this->Y();

            // get the Stabilization Parameter
            const std::shared_ptr< Stabilization_Parameter > & tSP = mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::GLS ) );

            // add contribution of stabilization term if stabilization parameter has been set
            if ( tSP != nullptr )
            {
                // GLS stabilization term
                tRes += aWStar * tSP->val()( 0 ) * this->GLSTestFunc() * this->Tau() * this->LY();
            }

            // assemble into set residual
            this->assemble_residual( tTempRes );

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Compressible_NS_Bulk::compute_residual - Residual contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Bulk::compute_jacobian( real aWStar )
        {
            // check leader field interpolators
#ifdef MORIS_HAVE_DEBUG
            this->check_field_interpolators();
#endif
            // check residual dof types
            MORIS_ASSERT( check_residual_dof_types( mResidualDofType  ),
                    "IWG_Compressible_NS_Bulk::compute_jacobian() - Only pressure or density primitive variables supported for now." );

            // check DoF dependencies
            MORIS_ASSERT( check_dof_dependencies( mSet, mResidualDofType, mRequestedLeaderGlobalDofTypes ),
                    "IWG_Compressible_NS_Bulk::compute_jacobian - Set of DoF dependencies not suppported." );

            // get number of space dimensions
            uint tNumSpaceDims = this->num_space_dims();

            // get total number of DoFs on Comp Flow Element
            uint tNumTotalBases = ( tNumSpaceDims + 2 ) * this->num_bases();

            // construct temporary Matrix for elemental Jacobian in standardized format
            Matrix< DDRMat > tTempJac( tNumTotalBases, tNumTotalBases, 0.0 );
            auto tJac = tTempJac( { 0, tNumTotalBases - 1 }, { 0, tNumTotalBases - 1 } );

            // get the material and constitutive models
            std::shared_ptr< Material_Model > tMM = mLeaderMM( static_cast< uint >( IWG_Material_Type::FLUID_MM ) );
            std::shared_ptr< Constitutive_Model > tCM = mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_CM ) );

            // add contribution from d(A0)/dDof * Y,t
            Matrix< DDRMat > tdAdY;
            eval_dAdY_VR( tMM, tCM, mLeaderFIManager, mResidualDofType, this->dYdt(), 0, tdAdY );
            tJac += aWStar * this->W_trans() * ( tdAdY * this->W() + this->A( 0 ) * this->dWdt() );

            // loop over contributions from A-matrices
            for ( uint iDim = 0; iDim < tNumSpaceDims; iDim++ )
            {
                // evaluate d(Ai)/dDof * Y,i
                eval_dAdY_VR( tMM, tCM, mLeaderFIManager, mResidualDofType, this->dYdx( iDim ), iDim + 1, tdAdY );

                // add contribution
                tJac += aWStar * this->W_trans() * ( tdAdY * this->W() + this->A( iDim + 1 ) * this->dWdx( iDim ) );

            }

            // get properties for K-matrices
            std::shared_ptr< Property > tPropMu    = mLeaderProp( static_cast< uint >( IWG_Property_Type::DYNAMIC_VISCOSITY ) );
            std::shared_ptr< Property > tPropKappa = mLeaderProp( static_cast< uint >( IWG_Property_Type::THERMAL_CONDUCTIVITY ) );

            // loop over contributions from K-matrices
            Matrix< DDRMat > dKdY;
            for ( uint iDim = 0; iDim < tNumSpaceDims; iDim++ )
            {
                for ( uint jDim = 0; jDim < tNumSpaceDims; jDim++ )
                {
                    // get dKij/dY * Y,ij
                    eval_dKdY_VR( tPropMu, tPropKappa, mLeaderFIManager, this->dYdx( jDim ), iDim, jDim, dKdY );

                    // add contributions from K-matrices
                    tJac += aWStar * this->dWdx_trans( iDim ) * ( dKdY * this->W() + K( iDim, jDim ) * this->dWdx( jDim ) );
                }
            }

            // contribution from body loads
            tJac += aWStar * this->W_trans() * ( this->C() + this->dCdY_VR( this->Y() ) ) * this->W();

            // get the Stabilization Parameter
            const std::shared_ptr< Stabilization_Parameter > & tSP = mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::GLS ) );

            // add contribution of stabilization term if stabilization parameter has been set
            if ( tSP != nullptr )
            {
                // GLS stabilization term
                tJac += aWStar * tSP->val()( 0 ) * (
                        this->GLSTestFunc() * ( this->dTaudY( this->LY() ) * this->W() + this->Tau() * ( this->LW() + this->dLdDofY() ) ) +
                        this->dGLSTestFuncdDof( this->Tau() * this->LY() ) );
            }

            // assemble jacobian into set jacobian
            this->assemble_jacobian( tTempJac );

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ) ,
                    "IWG_Compressible_NS_Bulk::compute_jacobian - Jacobian contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Bulk::compute_jacobian_and_residual( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Compressible_NS_Bulk::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Bulk::compute_dRdp( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Compressible_NS_Bulk::compute_dRdp - Not implemented." );
        }

        //------------------------------------------------------------------------------

        Matrix< DDRMat > IWG_Compressible_NS_Bulk::dSqrtMinvdu_FD( const uint aColInd, const real aPerturbation )
        {
            // get residual dof type index in set, start and end indices for residual dof type
            uint tLeaderFirstDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );

            // get number of state variables
            uint tNumStateVars = this->num_space_dims() + 2;

            // get number of leader and follower rows
            uint tNumRows = tNumStateVars;

            // get number of cols for jacobian
            uint tNumCols = tNumStateVars * this->num_bases();;

            // set size for FD jacobian
            Matrix< DDRMat > tFDdMdu(tNumRows, tNumCols, 0.0 );

            // compute jacobian by FD
            // this->compute_jacobian_FD( aWStar, aPerturbation );
            {
                // get the FD scheme info
                moris::Vector< moris::Vector< real > > tFDScheme;
                fd_scheme( fem::FDScheme_Type::POINT_5, tFDScheme );
                uint tNumFDPoints = tFDScheme( 0 ).size();

                // get leader number of dof types
                uint tLeaderNumDofTypes = mRequestedLeaderGlobalDofTypes.size();

                // loop over the IWG dof types
                for( uint iFI = 0; iFI < tLeaderNumDofTypes; iFI++ )
                {
                    // init dof counter
                    uint tDofCounter = 0;

                    // get the dof type
                    Vector< MSI::Dof_Type > & tDofType = mRequestedLeaderGlobalDofTypes( iFI );

                    // get the index for the dof type
                    sint tLeaderDepDofIndex   = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                    uint tLeaderDepStartIndex = mSet->get_jac_dof_assembly_map()( tLeaderFirstDofIndex )( tLeaderDepDofIndex, 0 );

                    // get field interpolator for dependency dof type
                    Field_Interpolator * tFI =
                            mLeaderFIManager->get_field_interpolators_for_type( tDofType( 0 ) );

                    // get number of leader FI bases and fields
                    uint tDerNumBases  = tFI->get_number_of_space_time_bases();
                    uint tDerNumFields = tFI->get_number_of_fields();

                    // coefficients for dof type wrt which derivative is computed
                    Matrix< DDRMat > tCoeff = tFI->get_coeff();

                    // loop over the coefficient column
                    for( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
                    {
                        // loop over the coefficient row
                        for( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++  )
                        {
                            // compute the perturbation absolute value
                            real tDeltaH = aPerturbation ;//* tCoeff( iCoeffRow, iCoeffCol );

                            // check that perturbation is not zero
                            if( std::abs( tDeltaH ) < 1e-12 )
                            {
                                tDeltaH = aPerturbation;
                            }

                            // set starting point for FD
                            uint tStartPoint = 0;

                            // loop over the points for FD
                            for( uint iPoint = tStartPoint; iPoint < tNumFDPoints; iPoint++ )
                            {
                                // reset the perturbed coefficients
                                Matrix< DDRMat > tCoeffPert = tCoeff;

                                // perturb the coefficient
                                tCoeffPert( iCoeffRow, iCoeffCol ) += tFDScheme( 0 )( iPoint ) * tDeltaH;

                                // set the perturbed coefficients to FI
                                tFI->set_coeff( tCoeffPert );
                                tFI->reset_eval_flags(); // not useful

                                // reset properties, CM and SP for IWG
                                this->reset_eval_flags();

                                // index of dof being FDed
                                uint tDofIndex = tLeaderDepStartIndex + tDofCounter;

                                // evaluate column of M, Minv, or sqrtMinv
                                Matrix< DDRMat > Mcol = this->SqrtMinv()( { 0, tNumStateVars - 1 }, { aColInd, aColInd } );

                                // assemble the Jacobian
                                tFDdMdu( { 0, tNumStateVars - 1 }, {tDofIndex, tDofIndex } ) +=
                                                tFDScheme( 1 )( iPoint ) *
                                                Mcol / ( tFDScheme( 2 )( 0 ) * tDeltaH );
                            }
                            // update dof counter
                            tDofCounter++;
                        }
                    }
                    // reset the coefficients values
                    tFI->set_coeff( tCoeff );
                }
            }

            // return
            return tFDdMdu;
        }

        //------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */

