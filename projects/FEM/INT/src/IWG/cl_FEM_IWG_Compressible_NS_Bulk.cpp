/*
 * cl_FEM_IWG_Compressible_NS_Bulk.cpp
 *
 *  Created on: Feb 10, 2021
 *      Author: wunsch
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

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        IWG_Compressible_NS_Bulk::IWG_Compressible_NS_Bulk()
        {
            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "DynamicViscosity" ]     = static_cast< uint >( IWG_Property_Type::DYNAMIC_VISCOSITY );
            mPropertyMap[ "ThermalConductivity" ]  = static_cast< uint >( IWG_Property_Type::THERMAL_CONDUCTIVITY );
            mPropertyMap[ "BodyForce" ]            = static_cast< uint >( IWG_Property_Type::BODY_FORCE );
            mPropertyMap[ "BodyHeatLoad" ]         = static_cast< uint >( IWG_Property_Type::BODY_HEAT_LOAD );

            // set size for the material model pointer cell
            mMasterMM.resize( static_cast< uint >( IWG_Material_Type::MAX_ENUM ), nullptr );

            // populate the material map
            mMaterialMap[ "FluidMM" ] = static_cast< uint >( IWG_Material_Type::FLUID_MM );

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "FluidCM" ] = static_cast< uint >( IWG_Constitutive_Type::FLUID_CM );

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the stabilization parameter map
            mStabilizationMap[ "GenericSP" ] = static_cast< uint >( IWG_Stabilization_Type::GENERIC );
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Bulk::reset_spec_eval_flags()
        {
            // reset eval flags
            mSpaceDimEval = true;
            mNumBasesEval = true;

            mVarSetEval = true;
            mTestFuncSetEval = true;

            mFluxAMatEval = true;
            mFluxADofMatEval = true;
            mFluxKMatEval = true;
            mKijiEval = true;
            mFluxKDofMatEval = true;

            mLYEval = true;
            mLWEval = true;
            mLDofYEval = true;
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Bulk::compute_residual( real aWStar )
        {
            // check master field interpolators
#ifdef DEBUG
            this->check_field_interpolators();
#endif
            // check residual dof types
            MORIS_ASSERT( check_residual_dof_types( mResidualDofType ), 
                    "IWG_Compressible_NS_Bulk::compute_jacobian() - Only pressure or density primitive variables supported for now." );

            // assemble flux matrices
            this->assemble_variable_set();
            this->eval_A_matrices();

            // get indeces for residual dof types, indices for assembly (FIXME: assembly only for primitive vars)
            uint tMasterDof1Index      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterDof2Index      = mSet->get_dof_index_for_type( mResidualDofType( 1 ), mtk::Master_Slave::MASTER );
            uint tMasterDof3Index      = mSet->get_dof_index_for_type( mResidualDofType( 2 ), mtk::Master_Slave::MASTER );
            uint tMasterRes1StartIndex = mSet->get_res_dof_assembly_map()( tMasterDof1Index )( 0, 0 );
            uint tMasterRes1StopIndex  = mSet->get_res_dof_assembly_map()( tMasterDof1Index )( 0, 1 );
            uint tMasterRes2StartIndex = mSet->get_res_dof_assembly_map()( tMasterDof2Index )( 0, 0 );
            uint tMasterRes2StopIndex  = mSet->get_res_dof_assembly_map()( tMasterDof2Index )( 0, 1 );
            uint tMasterRes3StartIndex = mSet->get_res_dof_assembly_map()( tMasterDof3Index )( 0, 0 );
            uint tMasterRes3StopIndex  = mSet->get_res_dof_assembly_map()( tMasterDof3Index )( 0, 1 );

            // get the FIs associated with each residual dof type
            Field_Interpolator * tFIFirstDofType =  mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );
            Field_Interpolator * tFIVelocity     =  mMasterFIManager->get_field_interpolators_for_type( mDofVelocity );
            Field_Interpolator * tFIThirdDofType =  mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 2 ) );

            // get number of space dimensions
            uint tNumSpaceDims = tFIVelocity->get_number_of_fields();

            // get the material and constitutive models
            std::shared_ptr< Material_Model > tMM = mMasterMM( static_cast< uint >( IWG_Material_Type::FLUID_MM ) );
            std::shared_ptr< Constitutive_Model > tCM = mMasterCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_CM ) );

            // get the mass body force property - not used for now
            // std::shared_ptr< Property > tPropBodyForce = mMasterProp( static_cast< uint >( IWG_Property_Type::BODY_FORCE ) );
            // std::shared_ptr< Property > tPropBodyHeatLoad = mMasterProp( static_cast< uint >( IWG_Property_Type::BODY_HEAT_LOAD ) );

            // get subview for complete residual
            auto tRes = mSet->get_residual()( 0 )( { tMasterRes1StartIndex, tMasterRes3StopIndex }, { 0, 0 } );

            // A0 matrix contribution
            tRes += aWStar * trans( this->W() ) * this->A( 0 ) * this->dYdt();

            // compute the second residual (velocity)
            mSet->get_residual()( 0 )( { tMasterRes2StartIndex, tMasterRes2StopIndex }, { 0, 0 } ) += aWStar * ( 
                    trans( tCM->testStrain() ) * this->MultipMat() * tCM->flux( CM_Function_Type::MECHANICAL ) ); 

            // compute the third residual (temperature)
            mSet->get_residual()( 0 )( { tMasterRes3StartIndex, tMasterRes3StopIndex }, { 0, 0 } ) += aWStar * ( 
                    trans( tFIThirdDofType->dnNdxn( 1 ) ) * ( 
                        tCM->flux( CM_Function_Type::WORK ) - tCM->flux( CM_Function_Type::THERMAL ) ) ); 

            // loop over A-Matrices
            for ( uint iA = 1; iA < tNumSpaceDims + 1; iA++ )
            {
                // compute residual
                tRes += aWStar * trans( this->W() ) * this->A( iA ) * this->dYdx( iA - 1 );
            }

            // get the Stabilization Parameter
            const std::shared_ptr< Stabilization_Parameter > & tSP = mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::GENERIC ) );

            // add contribution of stabilization term if stabilization parameter has been set
            if ( tSP != nullptr )
            {
                // get the strong form of the residual
                Matrix< DDRMat > tStrongRes;
                this->compute_residual_strong_form( tStrongRes );

                // compute the first residual (pressure or density)
                mSet->get_residual()( 0 )( { tMasterRes1StartIndex, tMasterRes1StopIndex }, { 0, 0 } ) +=  aWStar * tFIFirstDofType->N_trans() *
                        tSP->val()( 0 ) * tStrongRes( { 0, 0 } );

                // compute the second residual (velocity)
                mSet->get_residual()( 0 )( { tMasterRes2StartIndex, tMasterRes2StopIndex }, { 0, 0 } ) +=  aWStar * tFIVelocity->N_trans() *
                        tSP->val()( 0 ) * tStrongRes( { 1, tNumSpaceDims } );

                // compute the third residual (temperature)
                mSet->get_residual()( 0 )( { tMasterRes3StartIndex, tMasterRes3StopIndex }, { 0, 0 } ) +=  aWStar * tFIThirdDofType->N_trans() *
                        tSP->val()( 0 ) * tStrongRes( { tNumSpaceDims + 1, tNumSpaceDims + 1 } );
            }

// // debug for 2D
// if ( tNumSpaceDims == 2 ) {

// // debug - construct KijYj as a matrix
// Matrix< DDRMat > tKijYj( tNumSpaceDims + 2, 2, 0.0 );
// Matrix< DDRMat > tTau = tCM->flux( CM_Function_Type::MECHANICAL );
// Matrix< DDRMat > tQ = tCM->flux( CM_Function_Type::THERMAL );
// Matrix< DDRMat > tU = tFIVelocity->val();
// tKijYj( 1, 0 ) = tTau( 0 );
// tKijYj( 1, 1 ) = tTau( 2 );
// tKijYj( 2, 0 ) = tTau( 2 );
// tKijYj( 2, 1 ) = tTau( 1 );
// Matrix< DDRMat > tEFlux = tKijYj({1,2},{0,1}) * tU - tQ;
// tKijYj( 3, 0 ) = tEFlux( 0 );
// tKijYj( 3, 1 ) = tEFlux( 1 );

// // debug - construct KijYji as vector
// Matrix< DDRMat > tKijYji( tNumSpaceDims + 2, 1, 0.0 );
// tKijYji( { 1, 2 }, { 0, 0 } ) = 
//         tCM->divflux( CM_Function_Type::MECHANICAL ).matrix_data();
// tKijYji( { 3, 3 }, { 0, 0 } ) = 
//         tCM->divflux( CM_Function_Type::WORK ) - tCM->divflux( CM_Function_Type::THERMAL );

// Matrix< DDRMat > tStrongRes;
// this->compute_residual_strong_form( tStrongRes );

// // debug - write matrices to .hdf5 file
// std::string tMorisRoot = moris::get_base_moris_dir();
// std::string tHdf5FilePath = tMorisRoot + "/tmp/flux_matrices.hdf5";
// std::cout << "Outputting flux matrices to: " << tHdf5FilePath << " ... " << std::flush;
// hid_t tFileID = create_hdf5_file( tHdf5FilePath );
// herr_t tStatus = 0;
// save_matrix_to_hdf5_file( tFileID, "Residual", mSet->get_residual()( 0 ), tStatus );
// save_matrix_to_hdf5_file( tFileID, "ResStrongForm", tStrongRes, tStatus );
// save_matrix_to_hdf5_file( tFileID, "Y", mY, tStatus );
// save_matrix_to_hdf5_file( tFileID, "dYdt", mdYdt, tStatus );
// save_matrix_to_hdf5_file( tFileID, "dYdx", mdYdx, tStatus );
// save_matrix_to_hdf5_file( tFileID, "A0", mA(0), tStatus );
// save_matrix_to_hdf5_file( tFileID, "A1", mA(1), tStatus );
// save_matrix_to_hdf5_file( tFileID, "A2", mA(2), tStatus );
// save_matrix_to_hdf5_file( tFileID, "KijYj", tKijYj, tStatus );
// save_matrix_to_hdf5_file( tFileID, "KijYji", tKijYji, tStatus );
// save_matrix_to_hdf5_file( tFileID, "Tau", tCM->flux( CM_Function_Type::MECHANICAL ), tStatus );
// save_matrix_to_hdf5_file( tFileID, "q", tCM->flux( CM_Function_Type::THERMAL ), tStatus );
// save_matrix_to_hdf5_file( tFileID, "Fw", tCM->flux( CM_Function_Type::WORK ), tStatus );
// save_matrix_to_hdf5_file( tFileID, "divTau", tCM->divflux( CM_Function_Type::MECHANICAL ), tStatus );
// save_matrix_to_hdf5_file( tFileID, "divq", tCM->divflux( CM_Function_Type::THERMAL ), tStatus );
// save_matrix_to_hdf5_file( tFileID, "divFw", tCM->divflux( CM_Function_Type::WORK ), tStatus );
// close_hdf5_file( tFileID );
// std::cout << "Done \n" << std::flush;
// }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Compressible_NS_Bulk::compute_residual - Residual contains NAN or INF, exiting!");                                 
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Bulk::compute_jacobian( real aWStar )
        {
            // check master field interpolators
#ifdef DEBUG
            this->check_field_interpolators();
#endif
            // check residual dof types
            MORIS_ASSERT( check_residual_dof_types( mResidualDofType ), 
                    "IWG_Compressible_NS_Bulk::compute_jacobian() - Only pressure or density primitive variables supported for now." );

            // assemble flux matrices
            this->assemble_variable_set();
            this->assemble_test_function_set();
            this->eval_A_matrices();
            this->eval_A_DOF_matrices();

            // get indeces for residual dof types, indices for assembly (FIXME: assembly only for primitive vars)
            uint tMasterDof1Index      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterDof2Index      = mSet->get_dof_index_for_type( mResidualDofType( 1 ), mtk::Master_Slave::MASTER );
            uint tMasterDof3Index      = mSet->get_dof_index_for_type( mResidualDofType( 2 ), mtk::Master_Slave::MASTER );
            uint tMasterRes1StartIndex = mSet->get_res_dof_assembly_map()( tMasterDof1Index )( 0, 0 );
            uint tMasterRes1StopIndex  = mSet->get_res_dof_assembly_map()( tMasterDof1Index )( 0, 1 );
            uint tMasterRes2StartIndex = mSet->get_res_dof_assembly_map()( tMasterDof2Index )( 0, 0 );
            uint tMasterRes2StopIndex  = mSet->get_res_dof_assembly_map()( tMasterDof2Index )( 0, 1 );
            uint tMasterRes3StartIndex = mSet->get_res_dof_assembly_map()( tMasterDof3Index )( 0, 0 );
            uint tMasterRes3StopIndex  = mSet->get_res_dof_assembly_map()( tMasterDof3Index )( 0, 1 );

            // get the FIs associated with each residual dof type
            Field_Interpolator * tFIFirstDofType =  mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );
            Field_Interpolator * tFIVelocity =  mMasterFIManager->get_field_interpolators_for_type( mDofVelocity );
            Field_Interpolator * tFIThirdDofType =  mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 2 ) );

            // get number of space dimensions
            uint tNumSpaceDims = tFIVelocity->get_number_of_fields();   

            // get number of bases for the elements used
            uint tNumBases = tFIVelocity->get_number_of_space_time_bases();                     

            // get the material and constitutive models
            std::shared_ptr< Material_Model > tMM = mMasterMM( static_cast< uint >( IWG_Material_Type::FLUID_MM ) );
            std::shared_ptr< Constitutive_Model > tCM = mMasterCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_CM ) );

            // get the mass body force property - not used for now
            // std::shared_ptr< Property > tPropBodyForce = mMasterProp( static_cast< uint >( IWG_Property_Type::BODY_FORCE ) );
            // std::shared_ptr< Property > tPropBodyHeatLoad = mMasterProp( static_cast< uint >( IWG_Property_Type::BODY_HEAT_LOAD ) );

            // check DoF dependencies
            MORIS_ASSERT( check_dof_dependencies( mSet, mResidualDofType, mRequestedMasterGlobalDofTypes ), 
                    "IWG_Compressible_NS_Bulk::compute_jacobian - Set of DoF dependencies not suppported." );

            // get the indeces for assembly
            sint tDofFirstDepIndex     = mSet->get_dof_index_for_type( mRequestedMasterGlobalDofTypes( 0 )( 0 ), mtk::Master_Slave::MASTER );
            sint tDofThirdDepIndex     = mSet->get_dof_index_for_type( mRequestedMasterGlobalDofTypes( 2 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterDep1StartIndex = mSet->get_jac_dof_assembly_map()( tMasterDof1Index )( tDofFirstDepIndex, 0 );
            uint tMasterDep3StopIndex  = mSet->get_jac_dof_assembly_map()( tMasterDof3Index )( tDofThirdDepIndex, 1 );                

            // add contribution to first residual dof type
            mSet->get_jacobian()(
                    { tMasterRes1StartIndex, tMasterRes1StopIndex },
                    { tMasterDep1StartIndex, tMasterDep3StopIndex } ) += aWStar * ( tFIFirstDofType->N_trans() * (
                            trans( mdYdt ) * mADOF( 0 )( 0 ) +  mA( 0 )( { 0, 0 }, { 0, tNumSpaceDims + 1 } ) * mdWdt ) ); 

            // add contribution to velocity residual dof type

            // Build the product of Y' * A_DOF for velocity DoF type
            Matrix< DDRMat > tADOFY( tNumSpaceDims, ( tNumSpaceDims + 2 ) * tNumBases );
            for (uint iDim = 0; iDim < tNumSpaceDims; iDim++)
            {
                tADOFY( { iDim, iDim }, { 0, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = trans( mdYdt ) * mADOF( 0 )( iDim + 1 );
            }

            // add contribution
            mSet->get_jacobian()(
                    { tMasterRes2StartIndex, tMasterRes2StopIndex },
                    { tMasterDep1StartIndex, tMasterDep3StopIndex } ) += aWStar * ( tFIVelocity->N_trans() * (
                            tADOFY +  mA( 0 )( { 1, tNumSpaceDims }, { 0, tNumSpaceDims + 1 } ) * mdWdt ) );                        

            // add contribution to temperature residual dof type
            mSet->get_jacobian()(
                    { tMasterRes3StartIndex, tMasterRes3StopIndex },
                    { tMasterDep1StartIndex, tMasterDep3StopIndex } ) += aWStar * ( tFIThirdDofType->N_trans() * (
                            trans( mdYdt ) * mADOF( 0 )( tNumSpaceDims + 1 ) + 
                            mA( 0 )( { tNumSpaceDims + 1, tNumSpaceDims + 1 }, { 0, tNumSpaceDims + 1 } ) * mdWdt ) );
 
            // loop over A-Matrices
            for ( uint iA = 1; iA < tNumSpaceDims + 1; iA++ )
            {
                // add contribution to first residual dof type
                mSet->get_jacobian()(
                        { tMasterRes1StartIndex, tMasterRes1StopIndex },
                        { tMasterDep1StartIndex, tMasterDep3StopIndex } ) += aWStar * ( tFIFirstDofType->N_trans() * (
                                trans( this->dYdx( iA - 1 ) ) * mADOF( iA )( 0 ) + 
                                mA( iA )( { 0, 0 }, { 0, tNumSpaceDims + 1 } ) * mdWdx( iA - 1 ) ) ); 

                // add contribution to velocity residual dof type

                // Build the product of Y' * A_DOF for velocity DoF type
                Matrix< DDRMat > tADOFY( tNumSpaceDims, ( tNumSpaceDims + 2 ) * tNumBases );
                for (uint iDim = 0; iDim < tNumSpaceDims; iDim++)
                {
                    tADOFY( { iDim, iDim }, { 0, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                            trans( this->dYdx( iA - 1 ) ) * mADOF( iA )( iDim + 1 );
                }

                // // add contribution
                mSet->get_jacobian()(
                        { tMasterRes2StartIndex, tMasterRes2StopIndex },
                        { tMasterDep1StartIndex, tMasterDep3StopIndex } ) += aWStar * ( tFIVelocity->N_trans() * (
                                tADOFY + mA( iA )( { 1, tNumSpaceDims }, { 0, tNumSpaceDims + 1 } ) * mdWdx( iA - 1 ) ) );                        

                // add contribution to temperature residual dof type
                mSet->get_jacobian()(
                        { tMasterRes3StartIndex, tMasterRes3StopIndex },
                        { tMasterDep1StartIndex, tMasterDep3StopIndex } ) += aWStar * ( tFIThirdDofType->N_trans() * (
                                trans( this->dYdx( iA - 1 ) ) * mADOF( iA )( tNumSpaceDims + 1 ) + 
                                mA( iA )( { tNumSpaceDims + 1, tNumSpaceDims + 1 }, { 0, tNumSpaceDims + 1 } ) * mdWdx( iA - 1 ) ) );
            }

            // loop over DoF dependencies for K*Y,j term
            for (uint iDof = 0; iDof < mRequestedMasterGlobalDofTypes.size(); iDof++)
            {   
                // get the treated dof type
                Cell< MSI::Dof_Type > & tDofType = mRequestedMasterGlobalDofTypes( iDof );

                // get index
                uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( iDof ), mtk::Master_Slave::MASTER );
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tMasterDepStartIndex = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 0 );
                uint tMasterDepStopIndex  = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 1 );

                // add contribution
                mSet->get_jacobian()(
                        { tMasterRes2StartIndex, tMasterRes2StopIndex },
                        { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * ( trans( tCM->testStrain() ) * 
                                this->MultipMat() * tCM->dFluxdDOF( tDofType, CM_Function_Type::MECHANICAL ) );                        

                // add contribution to temperature residual dof type
                mSet->get_jacobian()(
                        { tMasterRes3StartIndex, tMasterRes3StopIndex },
                        { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * ( trans( tFIThirdDofType->dnNdxn( 1 ) ) * (
                                tCM->dFluxdDOF( tDofType, CM_Function_Type::WORK ) - tCM->dFluxdDOF( tDofType, CM_Function_Type::THERMAL ) ) );
            }

            // get the Stabilization Parameter
            const std::shared_ptr< Stabilization_Parameter > & tSP = mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::GENERIC ) );

            // add contribution of stabilization term if stabilization parameter has been set
            if ( tSP != nullptr )
            {
                // get the strong form of the jacobian
                Matrix< DDRMat > tStrongJac;
                this->compute_jacobian_strong_form( tStrongJac );

                // compute the first residual (pressure or density)
                mSet->get_jacobian()( 
                        { tMasterRes1StartIndex, tMasterRes1StopIndex }, 
                        { tMasterDep1StartIndex, tMasterDep3StopIndex } ) +=  aWStar * tFIFirstDofType->N_trans() * tSP->val()( 0 ) * 
                                tStrongJac( { 0, 0 }, { tMasterDep1StartIndex, tMasterDep3StopIndex } );

                // compute the second residual (velocity)
                mSet->get_jacobian()( 
                        { tMasterRes2StartIndex, tMasterRes2StopIndex }, 
                        { tMasterDep1StartIndex, tMasterDep3StopIndex } ) +=  aWStar * tFIVelocity->N_trans() * tSP->val()( 0 ) * 
                                tStrongJac( { 1, tNumSpaceDims }, { tMasterDep1StartIndex, tMasterDep3StopIndex } );

                // compute the third residual (temperature)
                mSet->get_jacobian()( 
                        { tMasterRes3StartIndex, tMasterRes3StopIndex }, 
                        { tMasterDep1StartIndex, tMasterDep3StopIndex } ) +=  aWStar * tFIThirdDofType->N_trans() * tSP->val()( 0 ) * 
                                tStrongJac( { tNumSpaceDims + 1, tNumSpaceDims + 1 }, { tMasterDep1StartIndex, tMasterDep3StopIndex } );
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ) ,
                    "IWG_Compressible_NS_Bulk::compute_jacobian - Jacobian contains NAN or INF, exiting!");                     
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Bulk::compute_jacobian_and_residual( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Compressible_NS_Bulk::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Bulk::compute_dRdp( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Compressible_NS_Bulk::compute_dRdp - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Bulk::compute_residual_strong_form( Matrix< DDRMat > & aRM )
        {
            // assemble flux matrices
            this->assemble_variable_set();
            this->eval_A_matrices();

            // get the CM
            std::shared_ptr< Constitutive_Model > tCM = mMasterCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_CM ) );

            // get the FIs associated with each residual dof type
            Field_Interpolator * tFIVelocity     =  mMasterFIManager->get_field_interpolators_for_type( mDofVelocity );

            // get number of space dimensions
            uint tNumSpaceDims = tFIVelocity->get_number_of_fields();    

            // =================
            // assemble Residual 
            // =================

            // contribution from A0 matrix
            aRM = mA( 0 ) * mdYdt;

            // loop over Ai matrices and add contribution
            for ( uint iA = 1; iA < tNumSpaceDims + 1; iA++ )
            {
                aRM += this->A( iA ) * this->dYdx( iA - 1 );
            }

            // contribution from (Kij*Y,j),i
            // first residual (pressure or density) -- no contribution
            // compute contribution to the second residual (velocity)
            aRM( { 1, tNumSpaceDims }, { 0, 0 } ) -= 
                    tCM->divflux( CM_Function_Type::MECHANICAL ).matrix_data();

            // compute contribution to the third residual (temperature)
            aRM( { tNumSpaceDims + 1, tNumSpaceDims + 1 }, { 0, 0 } ) -= 
                    tCM->divflux( CM_Function_Type::WORK ) - tCM->divflux( CM_Function_Type::THERMAL );

        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Bulk::compute_jacobian_strong_form( Matrix< DDRMat > & aJM )
        {
            // assemble flux matrices
            this->assemble_variable_set();
            this->assemble_test_function_set();
            this->eval_A_matrices();
            this->eval_A_DOF_matrices();

            // get indeces for residual dof types, indices for assembly (FIXME: assembly only for primitive vars)
            uint tMasterDof1Index      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterDof3Index      = mSet->get_dof_index_for_type( mResidualDofType( 2 ), mtk::Master_Slave::MASTER );

            // get the FIs associated with each residual dof type
            Field_Interpolator * tFIVelocity =  mMasterFIManager->get_field_interpolators_for_type( mDofVelocity );

            // get number of space dimensions
            uint tNumSpaceDims = tFIVelocity->get_number_of_fields();                      

            // get the material and constitutive models
            std::shared_ptr< Constitutive_Model > tCM = mMasterCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_CM ) );

            // check DoF dependencies
            MORIS_ASSERT( check_dof_dependencies( mSet, mResidualDofType, mRequestedMasterGlobalDofTypes ), 
                    "IWG_Compressible_NS_Bulk::compute_jacobian - Set of DoF dependencies not suppported. See error message above." );

            // get the indeces for assembly
            sint tDofFirstDepIndex     = mSet->get_dof_index_for_type( mRequestedMasterGlobalDofTypes( 0 )( 0 ), mtk::Master_Slave::MASTER );
            sint tDofThirdDepIndex     = mSet->get_dof_index_for_type( mRequestedMasterGlobalDofTypes( 2 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterDep1StartIndex = mSet->get_jac_dof_assembly_map()( tMasterDof1Index )( tDofFirstDepIndex, 0 );
            uint tMasterDep3StopIndex  = mSet->get_jac_dof_assembly_map()( tMasterDof3Index )( tDofThirdDepIndex, 1 );                

            // ==================
            // assmemble Jacobian
            // ==================
            
            // add contribution of A0 matrix
            aJM = mA( 0 ) * mdWdt;
            
            // loop over rows for which the DoF derivatives are stored
            for ( uint iR = 0; iR < tNumSpaceDims + 2; iR++ )
            {
                aJM( { iR, iR }, { tMasterDep1StartIndex, tMasterDep3StopIndex } ) += 
                        trans( mdYdt ) * mADOF( 0 )( iR );
            }

            // add contribution of Ai matrices
            for ( uint iA = 1; iA < tNumSpaceDims + 1; iA++ )
            {
                aJM += mA( iA ) * mdWdx( iA - 1 );

                // loop over rows for which the DoF derivatives are stored
                for ( uint iR = 0; iR < tNumSpaceDims + 2; iR++ )
                {
                    aJM( { iR, iR }, { tMasterDep1StartIndex, tMasterDep3StopIndex } ) += 
                            trans( this->dYdx( iA - 1 ) ) * mADOF( iA )( iR );
                }
            }

            // loop over DoF dependencies for K*Y,j term
            for (uint iDof = 0; iDof < mRequestedMasterGlobalDofTypes.size(); iDof++)
            {   
                // get the treated dof type
                Cell< MSI::Dof_Type > & tDofType = mRequestedMasterGlobalDofTypes( iDof );

                // get index
                uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( iDof ), mtk::Master_Slave::MASTER );
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tMasterDepStartIndex = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 0 );
                uint tMasterDepStopIndex  = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 1 );

                // add contribution
                aJM( { 1, tNumSpaceDims }, { tMasterDepStartIndex, tMasterDepStopIndex } ) -= 
                        tCM->ddivfluxdu( tDofType, CM_Function_Type::MECHANICAL ).matrix_data();                        

                // add contribution to temperature residual dof type
                aJM( { tNumSpaceDims + 1, tNumSpaceDims + 1 }, { tMasterDepStartIndex, tMasterDepStopIndex } ) -= 
                        tCM->ddivfluxdu( tDofType, CM_Function_Type::WORK ) - tCM->ddivfluxdu( tDofType, CM_Function_Type::THERMAL ); 
            }
        }

        //------------------------------------------------------------------------------

        uint IWG_Compressible_NS_Bulk::num_space_dims()
        {
            // check if number of spatial dimensions is known
            if ( !mSpaceDimEval )
            {
                return mNumSpaceDims;
            }

            // set eval flag
            mSpaceDimEval = false;
            
            // get CM
            std::shared_ptr< Constitutive_Model > tCM = 
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_CM ) );

            // get number of spatial dimensions from CM
            mNumSpaceDims = tCM->get_num_space_dims();

            // return
            return mNumSpaceDims;
        }

        //------------------------------------------------------------------------------

        uint IWG_Compressible_NS_Bulk::num_bases()
        {
            // check if number of spatial dimensions is known
            if ( !mNumBasesEval )
            {
                return mNumBasesPerField;
            }

            // set eval flag
            mNumBasesEval = false;
            
            // get first FI
            Field_Interpolator * tFI =  mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get number of spatial dimensions from CM
            mNumBasesPerField = tFI->get_number_of_space_time_bases();

            // return
            return mNumBasesPerField;
        }

        //------------------------------------------------------------------------------

        // FIXME provided directly by the field interpolator?
        void IWG_Compressible_NS_Bulk::compute_dnNdtn(
                Matrix< DDRMat > & adnNdtn )
        {
            // get the residual dof type FI (here velocity)
            Field_Interpolator * tVelocityFI =
                    mMasterFIManager->get_field_interpolators_for_type( mDofVelocity );                 

            // init size for dnNdtn
            uint tNumRowt = tVelocityFI->get_number_of_fields();
            uint tNumColt = tVelocityFI->dnNdtn( 1 ).n_cols();
            adnNdtn.set_size( tNumRowt, tNumRowt * tNumColt , 0.0 );

            // loop over the fields
            for( uint iField = 0; iField < tNumRowt; iField++ )
            {
                // fill the matrix for each dimension
                adnNdtn( { iField, iField }, { iField * tNumColt, ( iField + 1 ) * tNumColt - 1 } ) =
                        tVelocityFI->dnNdtn( 1 ).matrix_data();
            }
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & IWG_Compressible_NS_Bulk::MultipMat()
        {
            //build multiplication matrix
            //for 2D
            if( mMasterFIManager->get_field_interpolators_for_type( mDofVelocity )->get_number_of_fields() == 2 )
            {
                return mMultipMat2D;
            }
            // for 3D
            else
            {
                return mMultipMat3D;
            }
        }       

        //------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
