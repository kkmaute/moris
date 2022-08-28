/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_SP_YZBeta_Advection.cpp
 *
 */

#include "cl_FEM_SP_YZBeta_Advection.hpp"
#include "cl_FEM_Cluster.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
// LINALG/src
#include "fn_norm.hpp"
#include "fn_dot.hpp"

#include "fn_FEM_CM_Phase_State_Functions.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        SP_YZBeta_Advection::SP_YZBeta_Advection()
        {
            // set the property pointer cell size
            mMasterProp.resize( static_cast< uint >( Property_Type::MAX_ENUM ), nullptr );

            // populate the map
            mPropertyMap[ "Beta" ]           = static_cast< uint >( Property_Type::BETA_CONSTANT );
            mPropertyMap[ "ReferenceState" ] = static_cast< uint >( Property_Type::REFERENCE_STATE );

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "Diffusion" ] = static_cast< uint >( IWG_Constitutive_Type::DIFFUSION );
        }

        //------------------------------------------------------------------------------

        void
        SP_YZBeta_Advection::set_dof_type_list(
                moris::Cell< moris::Cell< MSI::Dof_Type > >& aDofTypes,
                moris::Cell< std::string >&                  aDofStrings,
                mtk::Master_Slave                            aIsMaster )
        {
            // switch on master slave
            switch ( aIsMaster )
            {
                case mtk::Master_Slave::MASTER:
                {
                    // set dof type list
                    mMasterDofTypes = aDofTypes;

                    // loop on dof type
                    for ( uint iDof = 0; iDof < aDofTypes.size(); iDof++ )
                    {
                        // get dof string
                        const std::string& tDofString = aDofStrings( iDof );

                        // get dof type
                        const MSI::Dof_Type tDofType = aDofTypes( iDof )( 0 );

                        // if scalar field
                        if ( tDofString == "ScalarField" )
                        {
                            mMasterDofScalarField = tDofType;
                        }
                        else
                        {
                            // error unknown dof string
                            MORIS_ERROR( false,
                                    "SP_YZBeta_Advection::set_dof_type_list - Unknown aDofString : %s \n",
                                    tDofString.c_str() );
                        }
                    }
                    break;
                }
                case mtk::Master_Slave::SLAVE:
                {
                    // set dof type list
                    mSlaveDofTypes = aDofTypes;
                    break;
                }
                default:
                    MORIS_ERROR( false, "SP_YZBeta_Advection::set_dof_type_list - unknown master slave type." );
            }
        }

        //------------------------------------------------------------------------------

        void
        SP_YZBeta_Advection::eval_SP()
        {
            // get the FI for scalar field
            Field_Interpolator* tScalarFI =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofScalarField );

            // get the diffusion CM
            const std::shared_ptr< Constitutive_Model >& tCMDiffusion =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::DIFFUSION ) );

            MORIS_ASSERT( tCMDiffusion != nullptr,
                    "SP_YZBeta_Advection::eval_SP - constitutive model not defined\n" );

            // get the beta property
            const std::shared_ptr< Property >& tPropBeta =
                    mMasterProp( static_cast< uint >( Property_Type::BETA_CONSTANT ) );

            MORIS_ASSERT( tPropBeta != nullptr,
                    "SP_YZBeta_Advection::eval_SP - beta parameter not defined\n" );

            const real tBeta = tPropBeta->val()( 0 );

            // get the reference state property
            const std::shared_ptr< Property >& tPropYref =
                    mMasterProp( static_cast< uint >( Property_Type::REFERENCE_STATE ) );

            MORIS_ASSERT( tPropBeta != nullptr,
                    "SP_YZBeta_Advection::eval_SP - reference state not defined\n" );

            const real tYref = tPropYref->val()( 0 );

            // compute norm of spatial gradient of scalar field; apply threshold
            const real tGradNorm = std::max( norm( tCMDiffusion->gradEnergy() ), mEpsilon );

            // get unit vector in direction of spatial gradient of scalar field
            Matrix< DDRMat > tJ = tCMDiffusion->gradEnergy() / tGradNorm;

            // compute hdc ( Eqn 14 in Bazilevs (2007) )
            const uint tNumNodes = tScalarFI->dnNdxn( 1 ).n_cols();

            // auxiliary variables to compute denominator in Eqn. 14
            real tHdcAux = 0.0;

            for ( uint iNode = 0; iNode < tNumNodes; iNode++ )
            {
                tHdcAux += std::abs( dot( tJ, tScalarFI->dnNdxn( 1 ).get_column( iNode ) ) );
            }

            // threshold tHdcAux
            tHdcAux = std::max( tHdcAux, mEpsilon );

            // compute hdc with Eqn (14) but omit factor 2 as cancels out in Eqn 12
            const real tHdcHalf = 1.0 / tHdcAux;

            // evaluate Eqn (12) without accounting for |Z|
            const real tKdc = 1.0 / tYref * std::pow( tGradNorm / tYref, tBeta - 2.0 ) * std::pow( tHdcHalf, tBeta );

            // compute stabilization parameter value
            mPPVal = { { tKdc } };
        }

        //------------------------------------------------------------------------------

        void
        SP_YZBeta_Advection::eval_dSPdMasterDOF(
                const moris::Cell< MSI::Dof_Type >& aDofTypes )
        {
            // get the dof type index
            const uint tDofIndex = mMasterGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the dof type FI
            Field_Interpolator* tFIDer =
                    mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set size for dSPdMasterDof, dTau1dDof, dTau3dDof
            mdPPdMasterDof( tDofIndex ).set_size( 1, tFIDer->get_number_of_space_time_coefficients() );

            // get the FI for scalar field
            Field_Interpolator* tScalarFI =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofScalarField );

            // get the diffusion CM
            const std::shared_ptr< Constitutive_Model >& tCMDiffusion =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::DIFFUSION ) );

            // get the beta property
            const std::shared_ptr< Property >& tPropBeta =
                    mMasterProp( static_cast< uint >( Property_Type::BETA_CONSTANT ) );

            const real tBeta = tPropBeta->val()( 0 );

            // get the reference state property
            const std::shared_ptr< Property >& tPropYref =
                    mMasterProp( static_cast< uint >( Property_Type::REFERENCE_STATE ) );

            const real tYref = tPropYref->val()( 0 );

            // compute norm of spatial gradient of scalar field; apply threshold
            const real tGradNorm = std::max( norm( tCMDiffusion->gradEnergy() ), mEpsilon );

            // get unit vector in direction of spatial gradient of scalar field
            Matrix< DDRMat > tJ = tCMDiffusion->gradEnergy() / tGradNorm;

            // compute hdc ( Eqn 14 in Bazilev (2007) )
            const uint tNumNodes = tScalarFI->dnNdxn( 1 ).n_cols();

            // auxiliary variables to compute denominator in Eqn. 14
            real tHdcAux = 0.0;

            for ( uint iNode = 0; iNode < tNumNodes; iNode++ )
            {
                tHdcAux += std::abs( dot( tJ, tScalarFI->dnNdxn( 1 ).get_column( iNode ) ) );
            }

            // threshold tHdcAux
            tHdcAux = std::max( tHdcAux, mEpsilon );

            // compute hdc with Eqn (14) but omit factor 2 as cancels out in Eqn 12
            const real tHdcHalf = 1.0 / tHdcAux;

            // if dof type is scalar field
            if ( tCMDiffusion->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative of hdc wrt scalar field dof
                Matrix< DDRMat > tdGradNormdu( 1, tFIDer->get_number_of_space_time_coefficients() );
                Matrix< DDRMat > tdHdcAuxdu( 1, tFIDer->get_number_of_space_time_coefficients(), 0.0 );
                Matrix< DDRMat > tdHdcHalfdu( 1, tFIDer->get_number_of_space_time_coefficients() );

                // compute derivative of the scalar field gradient norm (compute only derivative if not thresholded)
                if ( tGradNorm > mEpsilon )
                {
                    tdGradNormdu = trans( tCMDiffusion->gradEnergy() ) * tCMDiffusion->dGradEnergydDOF( aDofTypes ) / tGradNorm;
                }
                else
                {
                    tdGradNormdu.fill( 0.0 );
                }

                // compute derivative of unit vector in direction of spatial gradient of scalar field
                Matrix< DDRMat > tdJdu = tCMDiffusion->dGradEnergydDOF( aDofTypes ) / tGradNorm
                                       - tCMDiffusion->gradEnergy() * tdGradNormdu / tGradNorm / tGradNorm;

                // compute derivative of the abs term (compute only derivative if not thresholded)
                if ( tHdcAux > mEpsilon )
                {
                    uint tNumNodes = tScalarFI->dnNdxn( 1 ).n_cols();

                    for ( uint iNode = 0; iNode < tNumNodes; iNode++ )
                    {
                        const real tHdcAuxLoc = dot( tJ, tScalarFI->dnNdxn( 1 ).get_column( iNode ) );

                        const real tSign = tHdcAuxLoc < 0.0 ? -1.0 : 1.0;

                        tdHdcAuxdu += tSign * trans( tScalarFI->dnNdxn( 1 ).get_column( iNode ) ) * tdJdu;
                    }
                }

                // compute the derivative of hdc with Eqn (14) but omit factor 2
                tdHdcHalfdu = -1.0 / tHdcAux / tHdcAux * tdHdcAuxdu;

                // compute the derivative Eqn (12) without accounting for |Z|
                mdPPdMasterDof( tDofIndex ) =                                                       //
                        1.0 / tYref *                                                               //
                        ( ( tBeta - 2.0 ) / tYref * std::pow( tGradNorm / tYref, tBeta - 3.0 ) *    //
                                        std::pow( tHdcHalf, tBeta ) * tdGradNormdu                  //
                                + tBeta * std::pow( tGradNorm / tYref, tBeta - 2.0 ) *              //
                                          std::pow( tHdcHalf, tBeta - 1.0 ) * tdHdcHalfdu );
            }
            else
            {
                mdPPdMasterDof( tDofIndex ).fill( 0.0 );
            }

            // if beta depends on dof type
            MORIS_ASSERT( !tPropBeta->check_dof_dependency( aDofTypes ),
                    "SP_YZBeta_Advection::eval_dSPdMasterDOF - Dof depenced of beta not implemented yet.\n" );

            // if Yref depends on dof type
            MORIS_ASSERT( !tPropYref->check_dof_dependency( aDofTypes ),
                    "SP_YZBeta_Advection::eval_dSPdMasterDOF - Dof depenced of Yref not implemented yet.\n" );

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mdPPdMasterDof( tDofIndex ) ),
                    "SP_YZBeta_Advection::eval_dSPdMasterDOF - mdPPdMasterDof contains NAN or INF, exiting for tDofIndex = %d !\n",
                    tDofIndex );
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

