/*
 * cl_FEM_IWG_Compressible_NS_Bulk.cpp
 *
 *  Created on: Feb 10, 2021
 *      Author: wunsch
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG_Compressible_NS_Bulk.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"

// debug
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

            // set size for the constitutive model pointer cell
            mMasterMM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mMaterialMap[ "FluidMM" ] = static_cast< uint >( IWG_Material_Type::FLUID_MM );

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "FluidCM" ] = static_cast< uint >( IWG_Constitutive_Type::FLUID_CM );
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Bulk::reset_spec_eval_flags()
        {

// debug
//std::cout << "Resetting Eval Flags in Comp. Flow IWG. \n" << std::flush;

            // reset eval flags
            mFluxAMatEval = true;
            mFluxADofMatEval = true;
            mFluxKMatEval = true;
            mFluxKDofMatEval = true;
            mVarVecEval = true;
            mVarDofVecEval = true;
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Bulk::compute_residual( real aWStar )
        {
            // check master field interpolators
#ifdef DEBUG
            this->check_field_interpolators();
#endif
            // check residual dof types
            MORIS_ASSERT( this->check_residual_dof_types(), 
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
            Field_Interpolator * tFIVelocity =  mMasterFIManager->get_field_interpolators_for_type( mDofVelocity );
            Field_Interpolator * tFIThirdDofType =  mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 2 ) );

            // get number of space dimensions
            uint tNumSpaceDims = tFIVelocity->get_number_of_fields();

            // get the material and constitutive models
            std::shared_ptr< Material_Model > tMM = mMasterMM( static_cast< uint >( IWG_Material_Type::FLUID_MM ) );
            std::shared_ptr< Constitutive_Model > tCM = mMasterCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_CM ) );

            // get the mass body force property - not used for now
            // std::shared_ptr< Property > tPropBodyForce = mMasterProp( static_cast< uint >( IWG_Property_Type::BODY_FORCE ) );
            // std::shared_ptr< Property > tPropBodyHeatLoad = mMasterProp( static_cast< uint >( IWG_Property_Type::BODY_HEAT_LOAD ) );

// debug
//std::cout << "IWG_Compressible_NS_Bulk::compute_residual - Norm of Residual = " << norm( mSet->get_residual()( 0 ) ) << " \n" << std::flush;

            // compute the first residual (pressure or density)
            mSet->get_residual()( 0 )( { tMasterRes1StartIndex, tMasterRes1StopIndex }, { 0, 0 } ) += aWStar * ( tFIFirstDofType->N_trans() * (
                    mA( 0 )( { 0, 0 }, { 0, tNumSpaceDims + 1 } ) * mdYdt ) );

            // compute the second residual (velocity)
            mSet->get_residual()( 0 )( { tMasterRes2StartIndex, tMasterRes2StopIndex }, { 0, 0 } ) += aWStar * ( tFIVelocity->N_trans() * (
                    mA( 0 )( { 1, tNumSpaceDims }, { 0, tNumSpaceDims + 1 } ) * mdYdt ) ); 

            // compute the third residual (temperature)
            mSet->get_residual()( 0 )( { tMasterRes3StartIndex, tMasterRes3StopIndex }, { 0, 0 } ) += aWStar * ( tFIThirdDofType->N_trans() * (
                    mA( 0 )( { tNumSpaceDims + 1, tNumSpaceDims + 1 }, { 0, tNumSpaceDims + 1 } ) * mdYdt ) ); 

            // loop over A-Matrices
            for ( uint iA = 1; iA < tNumSpaceDims + 1; iA++ )
            {
                // compute the first residual (pressure or density)
                mSet->get_residual()( 0 )( { tMasterRes1StartIndex, tMasterRes1StopIndex }, { 0, 0 } ) += aWStar * ( tFIFirstDofType->N_trans() * (
                        mA( iA )( { 0, 0 }, { 0, tNumSpaceDims + 1 } ) * mdYdx( { 0, tNumSpaceDims + 1 }, { iA - 1, iA - 1 } ) ) );

                // compute the second residual (velocity)
                mSet->get_residual()( 0 )( { tMasterRes2StartIndex, tMasterRes2StopIndex }, { 0, 0 } ) += aWStar * ( tFIVelocity->N_trans() * (
                        mA( iA )( { 1, tNumSpaceDims }, { 0, tNumSpaceDims + 1 } ) * mdYdx( { 0, tNumSpaceDims + 1 }, { iA - 1, iA - 1 } ) ) ); 

                // compute the third residual (temperature)
                mSet->get_residual()( 0 )( { tMasterRes3StartIndex, tMasterRes3StopIndex }, { 0, 0 } ) += aWStar *  ( tFIThirdDofType->N_trans() * (
                        mA( iA )( { tNumSpaceDims + 1, tNumSpaceDims + 1 }, { 0, tNumSpaceDims + 1 } ) * mdYdx( { 0, tNumSpaceDims + 1 }, { iA - 1, iA - 1 } ) ) ); 
            }

// debug - finite difference flux jacobians
// mSet->get_residual()( 0 )( { tMasterRes1StartIndex, tMasterRes1StartIndex + tNumSpaceDims + 1 }, { 0, 0 } ) += aWStar * (  
//         trans( mA( 3 )( { 4, 4 }, { 0, tNumSpaceDims + 1 } ) ) );

// mSet->get_residual()( 0 )( { tMasterRes2StartIndex, tMasterRes2StartIndex + tNumSpaceDims + 1 }, { 0, 0 } ) += aWStar * (
//         trans( mA( 0 )( { 1, 1 }, { 0, tNumSpaceDims + 1 } ) ) );

// mSet->get_residual()( 0 )( { tMasterRes3StartIndex, tMasterRes3StartIndex + tNumSpaceDims + 1 }, { 0, 0 } ) += aWStar * (
//         trans( mA( 0 )( { tNumSpaceDims + 1, tNumSpaceDims + 1 }, { 0, tNumSpaceDims + 1 } ) ) );

// print( trans( mSet->get_residual()( 0 ) ), "Residual after assembly" );   
// MORIS_ASSERT( false, "Stop here by intention." ); 

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Compressible_NS_Bulk::compute_residual - Residual contains NAN or INF, exiting!");

// debug - write to hdf5 file
// std::cout << "aWStar = " << aWStar << " \n" << std::flush;
// std::string tMorisRoot = moris::get_base_moris_dir();
// std::string tHdf5FilePath = tMorisRoot + "/tmp/residual_vectors.hdf5";
// std::cout << "Outputting HDF5 to: " << tHdf5FilePath << " \n" << std::flush;
// hid_t tFileID = create_hdf5_file( tHdf5FilePath );
// herr_t tStatus = 0;
// save_matrix_to_hdf5_file( tFileID, "N_P", tFIFirstDofType->N(), tStatus );
// save_matrix_to_hdf5_file( tFileID, "N_V", tFIVelocity->N(), tStatus );
// save_matrix_to_hdf5_file( tFileID, "N_TEMP", tFIThirdDofType->N(), tStatus );
// save_matrix_to_hdf5_file( tFileID, "mdYdt", mdYdt, tStatus );
// save_matrix_to_hdf5_file( tFileID, "Residual", mSet->get_residual()( 0 ), tStatus );
// close_hdf5_file( tFileID );                                   
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Bulk::compute_jacobian( real aWStar )
        {
            // check master field interpolators
#ifdef DEBUG
            this->check_field_interpolators();
#endif
            // check residual dof types
            MORIS_ASSERT( this->check_residual_dof_types(), 
                "IWG_Compressible_NS_Bulk::compute_jacobian() - Only pressure or density primitive variables supported for now." );

            // assemble flux matrices
            this->assemble_variable_set();
            this->assemble_variable_DOF_set();
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
            MORIS_ASSERT( this->check_dof_dependencies(), "IWG_Compressible_NS_Bulk::compute_jacobian - Set of DoF dependencies not suppported." );

            // get the indeces for assembly
            sint tDofFirstDepIndex     = mSet->get_dof_index_for_type( mRequestedMasterGlobalDofTypes( 0 )( 0 ), mtk::Master_Slave::MASTER );
            sint tDofThirdDepIndex     = mSet->get_dof_index_for_type( mRequestedMasterGlobalDofTypes( 2 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterDep1StartIndex = mSet->get_jac_dof_assembly_map()( tMasterDof1Index )( tDofFirstDepIndex, 0 );
            uint tMasterDep3StopIndex  = mSet->get_jac_dof_assembly_map()( tMasterDof3Index )( tDofThirdDepIndex, 1 );                

// debug
//std::cout << "IWG_Compressible_NS_Bulk::compute_jacobian - Norm of Jacobian = " << norm( mSet->get_jacobian() ) << " \n" << std::flush;

            // add contribution to first residual dof type
            mSet->get_jacobian()(
                    { tMasterRes1StartIndex, tMasterRes1StopIndex },
                    { tMasterDep1StartIndex, tMasterDep3StopIndex } ) += aWStar * ( tFIFirstDofType->N_trans() * (
                            trans( mdYdt ) * mADOF( 0 )( 0 ) +  mA( 0 )( { 0, 0 }, { 0, tNumSpaceDims + 1 } ) * mdYdtDOF ) ); 

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
                            tADOFY +  mA( 0 )( { 1, tNumSpaceDims }, { 0, tNumSpaceDims + 1 } ) * mdYdtDOF ) );                        

            // add contribution to temperature residual dof type
            mSet->get_jacobian()(
                    { tMasterRes3StartIndex, tMasterRes3StopIndex },
                    { tMasterDep1StartIndex, tMasterDep3StopIndex } ) += aWStar * ( tFIThirdDofType->N_trans() * (
                            trans( mdYdt ) * mADOF( 0 )( tNumSpaceDims + 1 ) + 
                            mA( 0 )( { tNumSpaceDims + 1, tNumSpaceDims + 1 }, { 0, tNumSpaceDims + 1 } ) * mdYdtDOF ) );
 
            // loop over A-Matrices
            for ( uint iA = 1; iA < tNumSpaceDims + 1; iA++ )
            {
                // add contribution to first residual dof type
                mSet->get_jacobian()(
                        { tMasterRes1StartIndex, tMasterRes1StopIndex },
                        { tMasterDep1StartIndex, tMasterDep3StopIndex } ) += aWStar * ( tFIFirstDofType->N_trans() * (
                                trans( mdYdx( { 0, tNumSpaceDims + 1 }, { iA - 1, iA - 1 } ) ) * mADOF( iA )( 0 ) + 
                                mA( iA )( { 0, 0 }, { 0, tNumSpaceDims + 1 } ) * mdYdxDOF( iA - 1 ) ) ); 

                // add contribution to velocity residual dof type

                // Build the product of Y' * A_DOF for velocity DoF type
                Matrix< DDRMat > tADOFY( tNumSpaceDims, ( tNumSpaceDims + 2 ) * tNumBases );
                for (uint iDim = 0; iDim < tNumSpaceDims; iDim++)
                {
                    tADOFY( { iDim, iDim }, { 0, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                            trans( mdYdx( { 0, tNumSpaceDims + 1 }, { iA - 1, iA - 1 } ) ) * mADOF( iA )( iDim + 1 );
                }

                // // add contribution
                mSet->get_jacobian()(
                        { tMasterRes2StartIndex, tMasterRes2StopIndex },
                        { tMasterDep1StartIndex, tMasterDep3StopIndex } ) += aWStar * ( tFIVelocity->N_trans() * (
                                tADOFY + mA( iA )( { 1, tNumSpaceDims }, { 0, tNumSpaceDims + 1 } ) * mdYdxDOF( iA - 1 ) ) );                        

                // add contribution to temperature residual dof type
                mSet->get_jacobian()(
                        { tMasterRes3StartIndex, tMasterRes3StopIndex },
                        { tMasterDep1StartIndex, tMasterDep3StopIndex } ) += aWStar * ( tFIThirdDofType->N_trans() * (
                                trans( mdYdx( { 0, tNumSpaceDims + 1 }, { iA - 1, iA - 1 } ) ) * mADOF( iA )( tNumSpaceDims + 1 ) + 
                                mA( iA )( { tNumSpaceDims + 1, tNumSpaceDims + 1 }, { 0, tNumSpaceDims + 1 } ) * mdYdxDOF( iA - 1 ) ) );
            }

// debug - finite difference flux jacobians
// mSet->get_jacobian()(
//         { tMasterRes1StartIndex, tMasterRes1StartIndex + tNumSpaceDims + 1 },
//         { tMasterDep1StartIndex, tMasterDep3StopIndex } ) += aWStar * ( mADOF( 3 )( 4 ) );

// mSet->get_jacobian()(
//         { tMasterRes2StartIndex, tMasterRes2StartIndex + tNumSpaceDims + 1 },
//         { tMasterDep1StartIndex, tMasterDep3StopIndex } ) += aWStar * 0.0 * ( mADOF( 0 )( 1 ) );

// mSet->get_jacobian()(
//         { tMasterRes3StartIndex, tMasterRes3StartIndex + tNumSpaceDims + 1 },
//         { tMasterDep1StartIndex, tMasterDep3StopIndex } ) += aWStar * ( mADOF( 0 )( tNumSpaceDims + 1 ) );   
                                                              
// print( mSet->get_jacobian()( { tMasterRes1StartIndex, tMasterRes1StartIndex + tNumSpaceDims + 1 }, { tMasterDep1StartIndex, tMasterDep3StopIndex } ), 
//         "Jacobian after assembly" );   
//MORIS_ASSERT( false, "Stop here by intention." ); 

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ) ,
                    "IWG_Compressible_NS_Bulk::compute_jacobian - Jacobian contains NAN or INF, exiting!");

// debug - write to hdf5 file
// std::cout << "aWStar = " << aWStar << " \n" << std::flush;
// std::string tMorisRoot = moris::get_base_moris_dir();
// std::string tHdf5FilePath = tMorisRoot + "/tmp/jacobian_matrices.hdf5";
// std::cout << "Outputting HDF5 to: " << tHdf5FilePath << " \n" << std::flush;
// hid_t tFileID = create_hdf5_file( tHdf5FilePath );
// herr_t tStatus = 0;
// save_matrix_to_hdf5_file( tFileID, "mdYdtDOF", mdYdtDOF, tStatus );
// save_matrix_to_hdf5_file( tFileID, "Jacobian", mSet->get_jacobian(), tStatus );
// close_hdf5_file( tFileID );                      
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

        void IWG_Compressible_NS_Bulk::compute_residual_strong_form(
                Matrix< DDRMat > & aRM,
                real             & aRC )
        {
            MORIS_ERROR( false, "IWG_Compressible_NS_Bulk::compute_residual_strong_form - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Bulk::compute_jacobian_strong_form(
                moris::Cell< MSI::Dof_Type >   aDofTypes,
                Matrix< DDRMat >             & aJM,
                Matrix< DDRMat >             & aJC )
        {
            MORIS_ERROR( false, "IWG_Compressible_NS_Bulk::compute_jacobian_strong_form - Not implemented." );
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

        bool IWG_Compressible_NS_Bulk::check_residual_dof_types()
        {
            // initialize
            bool tCheck = true;

            // FIXME: only density and pressure primitive variables supported for now

            // check that there are exactly 3 residual DoF types
            tCheck = tCheck && ( mResidualDofType.size() == 3 );
            MORIS_ASSERT( mResidualDofType.size() == 3,
                    "IWG_Compressible_NS_Bulk::check_residual_dof_types() - List of Residual DoF types must be of length 3 for pressure primitive variables." );

            // check that the right DoF types are present
            tCheck = tCheck && ( mResidualDofType( 0 ) == MSI::Dof_Type::RHO || mResidualDofType( 0 ) == MSI::Dof_Type::P );
            MORIS_ASSERT( tCheck,
                    "IWG_Compressible_NS_Bulk::check_residual_dof_types() - First DoF type must be density or pressure." );
            tCheck = tCheck && ( mResidualDofType( 1 ) == MSI::Dof_Type::VX );
            tCheck = tCheck && ( mResidualDofType( 2 ) == MSI::Dof_Type::TEMP );
            MORIS_ASSERT( tCheck,
                    "IWG_Compressible_NS_Bulk::check_residual_dof_types() - Second and third DoF types must be velocity and temperature." );

            // set variables Types
            if ( mVariableSet == fem::Variable_Set::UNDEFINED )
            {
                if ( mResidualDofType( 0 ) == MSI::Dof_Type::RHO )
                {
                    mVariableSet = fem::Variable_Set::DENSITY_PRIMITIVE;
                }
                else if ( mResidualDofType( 0 ) == MSI::Dof_Type::P )
                {
                    mVariableSet = fem::Variable_Set::PRESSURE_PRIMITIVE;
                }
                else
                {
                    MORIS_ERROR( false,
                        "IWG_Compressible_NS_Bulk::check_residual_dof_types() - Something went wrong. Only Density or Pressure supported as first DoF types." );
                }
            }
            else if ( mResidualDofType( 0 ) == MSI::Dof_Type::RHO || mResidualDofType( 0 ) == MSI::Dof_Type::P )
            {
                // do nothing, variable set already defined
            }
            else
            {
                MORIS_ERROR( false,
                    "IWG_Compressible_NS_Bulk::check_residual_dof_types() - Something went wrong; variable set not supported yet." );
            }

            // return whether check was successful or not
            return tCheck;
        }

        //------------------------------------------------------------------------------

        bool IWG_Compressible_NS_Bulk::check_dof_dependencies()
        {  
            // compute and check the number of dof dependencies
            uint tNumDofDependencies = mRequestedMasterGlobalDofTypes.size(); 
            MORIS_ERROR( tNumDofDependencies == 3, "IWG_Compressible_NS_Bulk::check_dof_dependencies() - More than three DoF dependencies." );       

            // check dof dependencies and ordering
            bool tCheck = ( mRequestedMasterGlobalDofTypes( 0 )( 0 ) == MSI::Dof_Type::RHO ) || ( mRequestedMasterGlobalDofTypes( 0 )( 0 ) == MSI::Dof_Type::P );
            tCheck = tCheck && ( mRequestedMasterGlobalDofTypes( 1 )( 0 ) == MSI::Dof_Type::VX );
            tCheck = tCheck && ( mRequestedMasterGlobalDofTypes( 2 )( 0 ) == MSI::Dof_Type::TEMP );

            MORIS_ERROR( tCheck, 
                    "IWG_Compressible_NS_Bulk::check_dof_dependencies() - Only Pressure and Density Primitive variables supported for now."
                    "DoF ordering must be (rho,u,T) or (p,u,T)." );
        
            // check that the assembly map is a connected block            
            uint tMasterDof1Index      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterDof2Index      = mSet->get_dof_index_for_type( mResidualDofType( 1 ), mtk::Master_Slave::MASTER );
            uint tMasterDof3Index      = mSet->get_dof_index_for_type( mResidualDofType( 2 ), mtk::Master_Slave::MASTER );            
            sint tDofDep1Index         = mSet->get_dof_index_for_type( mRequestedMasterGlobalDofTypes( 0 )( 0 ), mtk::Master_Slave::MASTER );
            sint tDofDep2Index         = mSet->get_dof_index_for_type( mRequestedMasterGlobalDofTypes( 1 )( 0 ), mtk::Master_Slave::MASTER );
            sint tDofDep3Index         = mSet->get_dof_index_for_type( mRequestedMasterGlobalDofTypes( 2 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterDep1StopIndex  = mSet->get_jac_dof_assembly_map()( tMasterDof1Index )( tDofDep1Index, 1 );
            uint tMasterDep2StartIndex = mSet->get_jac_dof_assembly_map()( tMasterDof2Index )( tDofDep2Index, 0 );
            uint tMasterDep2StopIndex  = mSet->get_jac_dof_assembly_map()( tMasterDof2Index )( tDofDep2Index, 1 );
            uint tMasterDep3StartIndex = mSet->get_jac_dof_assembly_map()( tMasterDof3Index )( tDofDep3Index, 0 );
            MORIS_ERROR( ( tMasterDep1StopIndex + 1 == tMasterDep2StartIndex ) && ( tMasterDep2StopIndex + 1 == tMasterDep3StartIndex ),
                    "IWG_Compressible_NS_Bulk::check_dof_dependencies() - Assembly map is not connected." );

            return tCheck;

            std::cout << "Check passed. \n" << std::flush;
        }

        //------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
