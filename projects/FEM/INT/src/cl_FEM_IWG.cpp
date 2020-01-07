/*
 * cl_FEM_IWG.cpp
 *
 *  Created on: Nov 12, 2019
 *      Author: sonne
 */

#include "cl_FEM_IWG.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

namespace moris
{
    namespace fem
    {

        void IWG::get_non_unique_dof_types( moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // set the size of the dof type list for the set
            uint tCounter = 0;

            // get number of direct master dof dependencies
            for ( uint iDOF = 0; iDOF < mMasterDofTypes.size(); iDOF++ )
            {
                tCounter += mMasterDofTypes( iDOF ).size();
            }

            // get number of direct slave dof dependencies
            for ( uint iDOF = 0; iDOF < mSlaveDofTypes.size(); iDOF++ )
            {
                tCounter += mSlaveDofTypes( iDOF ).size();
            }

            // loop over the master properties
            for ( std::shared_ptr< Property > tProperty : mMasterProp )
            {
                if ( tProperty != nullptr )
                {
                    // get property non unique dof type list
                    moris::Cell< MSI::Dof_Type > tActiveDofType;
                    tProperty->get_non_unique_dof_types( tActiveDofType );

                    //update counter
                    tCounter += tActiveDofType.size();
                }
            }

            // loop over slave properties
            for ( std::shared_ptr< Property > tProperty : mSlaveProp )
            {
                if ( tProperty != nullptr )
                {
                    // get property nn unique dof type list
                    moris::Cell< MSI::Dof_Type > tActiveDofType;
                    tProperty->get_non_unique_dof_types( tActiveDofType );

                    // update counter
                    tCounter += tActiveDofType.size();
                }
            }

            // loop over master constitutive models
            for ( std::shared_ptr< Constitutive_Model > tCM : mMasterCM )
            {
                if ( tCM != nullptr )
                {
                    // get CM non unique dof type list
                    moris::Cell< MSI::Dof_Type > tActiveDofType;
                    tCM->get_non_unique_dof_types( tActiveDofType );

                    // update counter
                    tCounter += tActiveDofType.size();
                }
            }

            // loop over slave constitutive models
            for ( std::shared_ptr< Constitutive_Model > tCM : mSlaveCM )
            {
                if( tCM != nullptr )
                {
                    // get CM non unique dof type list
                    moris::Cell< MSI::Dof_Type >  tActiveDofType;
                    tCM->get_non_unique_dof_types( tActiveDofType );

                    // update counter
                    tCounter += tActiveDofType.size();
                }
            }

            // loop over master stabilization parameters
            for ( std::shared_ptr< Stabilization_Parameter > tSP : mStabilizationParam )
            {
                if ( tSP != nullptr )
                {
                    // get SP non unique dof type list
                    moris::Cell< MSI::Dof_Type >  tActiveDofType;
                    tSP->get_non_unique_global_dof_type_list( tActiveDofType );

                    // update counter
                    tCounter += tActiveDofType.size();
                }
            }

            // reserve memory for dof type list
            aDofTypes.reserve( tCounter );

            // loop over master dof direct dependencies
            for ( uint iDOF = 0; iDOF < mMasterDofTypes.size(); iDOF++ )
            {
                // populate the dof list
                aDofTypes.append( mMasterDofTypes( iDOF ) );
            }

            // loop over slave dof direct dependencies
            for ( uint iDOF = 0; iDOF < mSlaveDofTypes.size(); iDOF++ )
            {
                //populate the dof list
                aDofTypes.append( mSlaveDofTypes( iDOF )  );
            }

            // loop over master properties
            for ( std::shared_ptr< Property > tProperty : mMasterProp )
            {
                if ( tProperty != nullptr )
                {
                    // get property non unique dof type list
                    moris::Cell< MSI::Dof_Type > tActiveDofType;
                    tProperty->get_non_unique_dof_types( tActiveDofType );

                    // populate the dof list
                    aDofTypes.append( tActiveDofType );
                }
            }

            // loop over slave properties
            for ( std::shared_ptr< Property > tProperty : mSlaveProp )
            {
                if ( tProperty != nullptr )
                {
                    // get property non unique dof type list
                    moris::Cell< MSI::Dof_Type > tActiveDofType;
                    tProperty->get_non_unique_dof_types( tActiveDofType );

                    // populate the dof list
                    aDofTypes.append( tActiveDofType );
                }
            }

            // loop over the master constitutive models
            for ( std::shared_ptr< Constitutive_Model > tCM : mMasterCM )
            {
                if ( tCM != nullptr )
                {
                    // get CM non unique dof type list
                    moris::Cell< MSI::Dof_Type > tActiveDofType;
                    tCM->get_non_unique_dof_types( tActiveDofType );

                    // populate the dof list
                    aDofTypes.append( tActiveDofType );
                }
            }

            // loop over the slave constitutive models
            for ( std::shared_ptr< Constitutive_Model > tCM : mSlaveCM )
            {
                if( tCM != nullptr )
                {
                    // get CM non unique dof type list
                    moris::Cell< MSI::Dof_Type > tActiveDofType;
                    tCM->get_non_unique_dof_types( tActiveDofType );

                    // populate the dof list
                    aDofTypes.append( tActiveDofType );
                }
            }

            // loop over the stabilization parameters
            for ( std::shared_ptr< Stabilization_Parameter > tSP : mStabilizationParam )
            {
                if ( tSP != nullptr )
                {
                    // get SP non unique master dof type list
                    moris::Cell< MSI::Dof_Type >  tActiveDofType;
                    tSP->get_non_unique_global_dof_type_list( tActiveDofType );

                    // populate the dof list
                    aDofTypes.append( tActiveDofType );
                }
            }
        }

//------------------------------------------------------------------------------
        void IWG::build_global_dof_type_list()
        {
            // MASTER-------------------------------------------------------
            // set size for the global dof type list
            uint tNumDofTypes = mSet->get_num_dof_types();

            // set size for the global dof type list
            mMasterGlobalDofTypes.reserve( tNumDofTypes );

            // set a size for the checkList ( used to avoid repeating a dof type)
            Matrix< DDSMat > tCheckList( tNumDofTypes, 1, -1 );

            // get dof type from direct dependencies
            for ( uint iDOF = 0; iDOF < mMasterDofTypes.size(); iDOF++ )
            {
                // get set index for dof type
                sint tDofTypeIndex = mSet->get_dof_index_for_type_1( mMasterDofTypes( iDOF )( 0 ) );  //FIXME'

                // put the dof type in the checklist
                tCheckList( tDofTypeIndex ) = 1;

                // put the dof type in the global type list
                mMasterGlobalDofTypes.push_back( mMasterDofTypes( iDOF ) );
            }

            // get dof type from master properties
            for ( std::shared_ptr< Property > tProperty : mMasterProp )
            {
                if ( tProperty != nullptr )
                {
                    // get dof types for property
                    moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType
                    = tProperty->get_dof_type_list();

                    // loop on property dof type
                    for ( uint iDOF = 0; iDOF < tActiveDofType.size(); iDOF++ )
                    {
                        // get set index for dof type
                        sint tDofTypeindex = mSet->get_dof_index_for_type_1( tActiveDofType( iDOF )( 0 ) );

                        // if dof enum not in the list
                        if ( tCheckList( tDofTypeindex) != 1 )
                        {
                            // put the dof type in the checklist
                            tCheckList( tDofTypeindex ) = 1;

                            // put the dof type in the global type list
                            mMasterGlobalDofTypes.push_back( tActiveDofType( iDOF ) );
                        }
                    }
                }
            }

            // get dof type from master constitutive models
            for ( std::shared_ptr< Constitutive_Model > tCM : mMasterCM )
            {
                if ( tCM != nullptr )
                {
                    // get dof types for constitutive model
                    moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType
                    = tCM->get_global_dof_type_list();

                    // loop on property dof type
                    for ( uint iDOF = 0; iDOF < tActiveDofType.size(); iDOF++ )
                    {
                        // get set index for dof type
                        sint tDofTypeindex = mSet->get_dof_index_for_type_1( tActiveDofType( iDOF )( 0 ) );

                        // if dof enum not in the list
                        if ( tCheckList( tDofTypeindex) != 1 )
                        {
                            // put the dof type in the checklist
                            tCheckList( tDofTypeindex ) = 1;

                            // put the dof type in the global type list
                            mMasterGlobalDofTypes.push_back( tActiveDofType( iDOF ) );
                        }
                    }
                }
            }

            // get dof type from master stabilization parameters
            for ( std::shared_ptr< Stabilization_Parameter > tSP : mStabilizationParam )
            {
                if ( tSP != nullptr )
                {
                    // get dof types for constitutive model
                    moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType
                    = tSP->get_global_dof_type_list( mtk::Master_Slave::MASTER );

                    // loop on property dof type
                    for ( uint iDOF = 0; iDOF < tActiveDofType.size(); iDOF++ )
                    {
                        // get set index for dof type
                        sint tDofTypeindex = mSet->get_dof_index_for_type_1( tActiveDofType( iDOF )( 0 ) );

                        // if dof enum not in the list
                        if ( tCheckList( tDofTypeindex) != 1 )
                        {
                            // put the dof type in the checklist
                            tCheckList( tDofTypeindex ) = 1;

                            // put the dof type in the global type list
                            mMasterGlobalDofTypes.push_back( tActiveDofType( iDOF ) );
                        }
                    }
                }
            }

            // reduce size of dof list to fit unique list
            mMasterGlobalDofTypes.shrink_to_fit();

            // SLAVE--------------------------------------------------------

            // set size for the global dof type list
            mSlaveGlobalDofTypes.reserve( tNumDofTypes );

            // set a size for the checkList ( used to avoid repeating a dof type)
            tCheckList.fill( -1 );

            // get dof type from slave direct dependencies
            for ( uint iDOF = 0; iDOF < mSlaveDofTypes.size(); iDOF++ )
            {
                // get set index for dof type
                sint tDofTypeindex = mSet->get_dof_index_for_type_1( mSlaveDofTypes( iDOF )( 0 ) );

                // put the dof type in the checklist
                tCheckList( tDofTypeindex ) = 1;

                // put the dof type in the global type list
                mSlaveGlobalDofTypes.push_back( mSlaveDofTypes( iDOF ) );
            }

            // get dof type from master properties
            for ( std::shared_ptr< Property > tProperty : mSlaveProp )
            {
                if ( tProperty != nullptr )
                {
                    // get dof types for property
                    moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType
                    = tProperty->get_dof_type_list();

                    // loop on property dof type
                    for ( uint iDOF = 0; iDOF < tActiveDofType.size(); iDOF++ )
                    {
                        // get set index for dof type
                        sint tDofTypeindex = mSet->get_dof_index_for_type_1( tActiveDofType( iDOF )( 0 ) );

                        // if dof enum not in the list
                        if ( tCheckList( tDofTypeindex) != 1 )
                        {
                            // put the dof type in the checklist
                            tCheckList( tDofTypeindex ) = 1;

                            // put the dof type in the global type list
                            mSlaveGlobalDofTypes.push_back( tActiveDofType( iDOF ) );
                        }
                    }
                }
            }

            // get dof type from slave constitutive models
            for ( std::shared_ptr< Constitutive_Model > tCM : mSlaveCM )
            {
                if ( tCM != nullptr )
                {
                    // get dof types for constitutive model
                    moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType
                    = tCM->get_global_dof_type_list();

                    // loop on property dof type
                    for ( uint iDOF = 0; iDOF < tActiveDofType.size(); iDOF++ )
                    {
                        // get set index for dof type
                        sint tDofTypeIndex = mSet->get_dof_index_for_type_1( tActiveDofType( iDOF )( 0 ) );

                        // if dof enum not in the list
                        if ( tCheckList( tDofTypeIndex ) != 1 )
                        {
                            // put the dof type in the checklist
                            tCheckList( tDofTypeIndex ) = 1;

                            // put the dof type in the global type list
                            mSlaveGlobalDofTypes.push_back( tActiveDofType( iDOF ) );
                        }
                    }
                }
            }

            // get dof type from stabilization parameters
            for ( std::shared_ptr< Stabilization_Parameter > tSP : mStabilizationParam )
            {
                if ( tSP != nullptr )
                {
                    // get dof types for constitutive model
                    moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType
                    = tSP->get_global_dof_type_list( mtk::Master_Slave::SLAVE );

                    // loop on property dof type
                    for ( uint iDOF = 0; iDOF < tActiveDofType.size(); iDOF++ )
                    {
                        // get set index for dof type
                        sint tDofTypeindex = mSet->get_dof_index_for_type_1( tActiveDofType( iDOF )( 0 ) );

                        // if dof enum not in the list
                        if ( tCheckList( tDofTypeindex) != 1 )
                        {
                            // put the dof type in the checklist
                            tCheckList( tDofTypeindex ) = 1;

                            // put the dof type in the global type list
                            mSlaveGlobalDofTypes.push_back( tActiveDofType( iDOF ) );
                        }
                    }
                }
            }

            // reduce size of dof list to fit unique list
             mSlaveGlobalDofTypes.shrink_to_fit();
        }

//------------------------------------------------------------------------------

void IWG::build_requested_dof_type_list( const bool aItResidual )
{
    mRequestedMasterGlobalDofTypes.clear();
    mRequestedSlaveGlobalDofTypes .clear();

    if ( aItResidual )
    {
        moris::Cell< moris::Cell< enum MSI::Dof_Type > > tRequestedDofTypes =  mSet->get_secundary_dof_types();

        mRequestedMasterGlobalDofTypes.reserve( tRequestedDofTypes.size() );
        mRequestedSlaveGlobalDofTypes .reserve( tRequestedDofTypes.size() );

        for( auto tDofTypes : tRequestedDofTypes )
        {
            for ( uint Ik = 0; Ik < mMasterGlobalDofTypes.size(); Ik++ )
            {
                if( mMasterGlobalDofTypes( Ik )( 0 ) == tDofTypes( 0 ) )
                {
                    mRequestedMasterGlobalDofTypes.push_back( mMasterGlobalDofTypes( Ik ) );

                    break;
                }
            }

            for ( uint Ik = 0; Ik < mSlaveGlobalDofTypes.size(); Ik++ )
            {
                if( mSlaveGlobalDofTypes( Ik )( 0 ) == tDofTypes( 0 ) )
                {
                    mRequestedSlaveGlobalDofTypes.push_back( mSlaveGlobalDofTypes( Ik ) );

                    break;
                }
            }
        }
    }
    else
    {
        Cell < enum MSI::Dof_Type > tRequestedDofTypes = mSet->get_requested_dof_types();

        mRequestedMasterGlobalDofTypes.reserve( tRequestedDofTypes.size() );
        mRequestedSlaveGlobalDofTypes .reserve( tRequestedDofTypes.size() );

        for( auto tDofTypes : tRequestedDofTypes )
        {
            for ( uint Ik = 0; Ik < mMasterGlobalDofTypes.size(); Ik++ )
            {
                if( mMasterGlobalDofTypes( Ik )( 0 ) == tDofTypes )
                {
                    mRequestedMasterGlobalDofTypes.push_back( mMasterGlobalDofTypes( Ik ) );
                    break;
                }
            }

            for ( uint Ik = 0; Ik < mSlaveGlobalDofTypes.size(); Ik++ )
            {
                if( mSlaveGlobalDofTypes( Ik )( 0 ) == tDofTypes )
                {
                    mRequestedSlaveGlobalDofTypes.push_back( mSlaveGlobalDofTypes( Ik ) );

                    break;
                }
            }

            if( mResidualDofType( 0 ) == tDofTypes )
            {
                mResidualDofTypeRequested = true;
            }
        }
    }

    mRequestedMasterGlobalDofTypes.shrink_to_fit();
    mRequestedSlaveGlobalDofTypes.shrink_to_fit();
}

//------------------------------------------------------------------------------

        void IWG::set_dof_field_interpolators( mtk::Master_Slave aIsMaster )
        {
            // set field interpolators for the SP
            for( std::shared_ptr< Stabilization_Parameter > tSP : this->get_stabilization_parameters() )
            {
                if ( tSP != nullptr )
                {
                    // get the list of dof types for the SP
                    moris::Cell< moris::Cell< MSI::Dof_Type > > tSPDofTypes = tSP->get_global_dof_type_list( aIsMaster );

                    // get the number of dof type for the SP
                    uint tNumDofTypes = tSPDofTypes.size();

                    // set the size of the field interpolators list for the SP
                    moris::Cell< Field_Interpolator* > tSPFIs( tNumDofTypes, nullptr );

                    // loop over the dof types
                    for( uint iDof = 0; iDof < tNumDofTypes; iDof++ )
                    {
                        // grab the field interpolator for the dof type
                        tSPFIs( iDof ) = mFieldInterpolatorManager->get_field_interpolators_for_type( tSPDofTypes( iDof )( 0 ), aIsMaster );
                    }

                    // set the field interpolators for the SP
                    tSP->set_dof_field_interpolators( tSPFIs, aIsMaster );

                    // set the field interpolator manager for the SP
                    tSP->set_field_interpolator_manager( mFieldInterpolatorManager );

                    // set th efem set pointer for the SP
                    tSP->set_set_pointer( mSet );
                }
            }

            // set field interpolators for constitutive models
            for( std::shared_ptr< Constitutive_Model > tCM : this->get_constitutive_models( aIsMaster ) )
            {
                if ( tCM != nullptr )
                {
                    // get the list of dof types for the CM
                    moris::Cell< moris::Cell< MSI::Dof_Type > > tCMDofTypes = tCM->get_global_dof_type_list();

                    // get the number of dof type for the CM
                    uint tNumDofTypes = tCMDofTypes.size();

                    // set the size of the field interpolators list for the CM
                    moris::Cell< Field_Interpolator* > tCMFIs( tNumDofTypes, nullptr );

                    // loop over the dof types
                    for( uint iDof = 0; iDof < tNumDofTypes; iDof++ )
                    {
                        // fill the field interpolators list for the CM
                        tCMFIs( iDof ) = mFieldInterpolatorManager->get_field_interpolators_for_type( tCMDofTypes( iDof )( 0 ), aIsMaster );
                    }

                    // set the field interpolators for the CM
                    tCM->set_dof_field_interpolators( tCMFIs );

                    // set the field interpolator manager for the CM
                    tCM->set_field_interpolator_manager( mFieldInterpolatorManager );

                    // set the fem set pointe for the CM
                    tCM->set_set_pointer( mSet );
                }
            }

            // set field interpolators for properties
            for( std::shared_ptr< Property > tProp : this->get_properties( aIsMaster ) )
            {
                if ( tProp != nullptr )
                {
                    // get the list of dof types for the property
                    moris::Cell< moris::Cell< MSI::Dof_Type > > tPropDofTypes = tProp->get_dof_type_list();

                    // get the number of dof type for the property
                    uint tNumDofTypes = tPropDofTypes.size();

                    // set the size of the field interpolators list for the property
                    moris::Cell< Field_Interpolator* > tPropFIs( tNumDofTypes, nullptr );

                    // loop over the dof types
                    for( uint iDof = 0; iDof < tNumDofTypes; iDof++ )
                    {
                        tPropFIs( iDof ) = mFieldInterpolatorManager->get_field_interpolators_for_type( tPropDofTypes( iDof )( 0 ), aIsMaster );
                    }

                    // set the field interpolators for the property
                    tProp->set_dof_field_interpolators( tPropFIs );

                    // set the field interpolator manager for the property
                    tProp->set_field_interpolator_manager( mFieldInterpolatorManager );

                    // set the fem set pointer for the property
                    tProp->set_set_pointer( mSet );
                }
            }
        }

//------------------------------------------------------------------------------

void IWG::check_dof_field_interpolators( mtk::Master_Slave aIsMaster )
{
    if ( aIsMaster == mtk::Master_Slave::MASTER)
    {
        // loop over the field interpolator pointers
        for( uint iFI = 0; iFI < mRequestedMasterGlobalDofTypes.size(); iFI++ )
        {
            // check that the field interpolator was set
            MORIS_ASSERT( mFieldInterpolatorManager->get_field_interpolators_for_type( mRequestedMasterGlobalDofTypes( iFI )( 0 ), aIsMaster ) != nullptr,
                    "IWG::check_dof_field_interpolators - FI missing. " );
        }
    }
    else
    {
        // loop over the field interpolator pointers
        for( uint iFI = 0; iFI < mRequestedSlaveGlobalDofTypes.size(); iFI++ )
        {
            // check that the field interpolator was set
            MORIS_ASSERT( mFieldInterpolatorManager->get_field_interpolators_for_type( mRequestedSlaveGlobalDofTypes( iFI )( 0 ), aIsMaster ) != nullptr,
                    "IWG::check_dof_field_interpolators - FI missing. " );
        }
    }
}

//------------------------------------------------------------------------------

void IWG::compute_jacobian_FD( real                                             aWStar,
                               real                                             aPerturbation,
                               moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobiansFD )
{
    // get master number of dof types
    uint tNumDofType = mRequestedMasterGlobalDofTypes.size();

    aJacobiansFD.resize( 1 );
    aJacobiansFD( 0 ).resize( tNumDofType );

    // loop over the IWG dof types
    for( uint iFI = 0; iFI < tNumDofType; iFI++ )
    {
        uint tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
        uint tDepIndex = mSet->get_dof_index_for_type( mRequestedMasterGlobalDofTypes( iFI )( 0 ), mtk::Master_Slave::MASTER );

        uint tNumRows = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 ) - mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 ) + 1;
        uint tNumCols = mSet->get_jac_dof_assembly_map()( tDofIndex )( tDepIndex, 1 ) - mSet->get_jac_dof_assembly_map()( tDofIndex )( tDepIndex, 0 ) + 1;

        aJacobiansFD( 0 )( iFI ).set_size( tNumRows, tNumCols, 0.0 );

        Field_Interpolator * tFI = mFieldInterpolatorManager->get_field_interpolators_for_type( mRequestedMasterGlobalDofTypes( iFI )( 0 ), mtk::Master_Slave::MASTER );

        // get number of master FI bases and fields
        uint tDerNumBases  = tFI->get_number_of_space_time_bases();
        uint tDerNumFields = tFI->get_number_of_fields();

        // coefficients for dof type wrt which derivative is computed
        Matrix< DDRMat > tCoeff = tFI->get_coeff();

        // init dof counter
        uint tDofCounter = 0;

        // loop over the coefficient column
        for( uint iCoeffCol = 0; iCoeffCol< tDerNumFields; iCoeffCol++ )
        {
            // loop over the coefficient row
            for( uint iCoeffRow = 0; iCoeffRow< tDerNumBases; iCoeffRow++  )
            {
                // perturbation of the coefficent
                Matrix< DDRMat > tCoeffPert = tCoeff;
                tCoeffPert( iCoeffRow, iCoeffCol ) = tCoeffPert( iCoeffRow, iCoeffCol ) + aPerturbation * tCoeffPert( iCoeffRow, iCoeffCol );

                // setting the perturbed coefficients
                tFI->set_coeff( tCoeffPert );

                // reset properties, CM and SP for IWG
                this->reset_eval_flags();

                // evaluate the residual
                mSet->get_residual().fill( 0.0 );
                this->compute_residual( aWStar );

                Matrix< DDRMat > tResidual_Plus
                =  mSet->get_residual()( { mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 ) }, { 0, 0 } );

                // perturbation of the coefficent
                tCoeffPert = tCoeff;
                tCoeffPert( iCoeffRow, iCoeffCol ) = tCoeffPert( iCoeffRow, iCoeffCol ) - aPerturbation * tCoeffPert( iCoeffRow, iCoeffCol );

                // setting the perturbed coefficients
                tFI->set_coeff( tCoeffPert );

                // reset properties, CM and SP for IWG
                this->reset_eval_flags();

                // evaluate the residual
                mSet->get_residual().fill( 0.0 );
                this->compute_residual( aWStar );

                Matrix< DDRMat > tResidual_Minus
                =  mSet->get_residual()( { mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 ) }, { 0, 0 } );

                // evaluate Jacobian
                aJacobiansFD( 0 )( iFI ).get_column( tDofCounter )
                                                   = ( tResidual_Plus - tResidual_Minus ) / ( 2.0 * aPerturbation * tCoeff( iCoeffRow, iCoeffCol ) );

                // update dof counter
                tDofCounter++;
            }
        }
        // reset the coefficients values
        tFI->set_coeff( tCoeff );
    }
}

//------------------------------------------------------------------------------

 void IWG::compute_jacobian_FD_double( real                                             aWStar,
                                       real                                             aPerturbation,
                                       moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobiansFD )
{
    // get master and slave number of dof types
    uint tMasterNumDofType = mRequestedMasterGlobalDofTypes.size();
    uint tSlaveNumDofType  = mRequestedSlaveGlobalDofTypes.size();

    aJacobiansFD.resize( 2 );
    aJacobiansFD( 0 ).resize( tMasterNumDofType + tSlaveNumDofType );
    aJacobiansFD( 1 ).resize( tMasterNumDofType + tSlaveNumDofType );

    // loop over the master dof types
    for( uint iFI = 0; iFI < tMasterNumDofType; iFI++ )
    {
        uint tDofIndexMaster = mSet->get_dof_index_for_type( mRequestedMasterGlobalDofTypes( iFI )( 0 ), mtk::Master_Slave::MASTER );
        uint tDofIndexSlave  = mSet->get_dof_index_for_type( mRequestedMasterGlobalDofTypes( iFI )( 0 ), mtk::Master_Slave::SLAVE  );

        uint tNumRowsMaster = mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 1 ) - mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 0 ) + 1;
        uint tNumColsMaster = mSet->get_jac_dof_assembly_map()( tDofIndexMaster )( tDofIndexMaster, 1 ) - mSet->get_jac_dof_assembly_map()( tDofIndexMaster )( tDofIndexMaster, 0 ) + 1;
        uint tNumRowsSlave  = mSet->get_res_dof_assembly_map()( tDofIndexSlave  )( 0 , 1 ) - mSet->get_res_dof_assembly_map()( tDofIndexSlave  )( 0 , 0 ) + 1;
        uint tNumColsSlave  = mSet->get_jac_dof_assembly_map()( tDofIndexSlave  )( tDofIndexSlave , 1 ) - mSet->get_jac_dof_assembly_map()( tDofIndexSlave  )( tDofIndexSlave , 0 ) + 1;

        aJacobiansFD( 0 )( iFI ).set_size( tNumRowsMaster, tNumColsMaster, 0.0 );
        aJacobiansFD( 1 )( iFI ).set_size( tNumRowsSlave, tNumColsSlave, 0.0 );

        Field_Interpolator * tFI = mFieldInterpolatorManager->get_field_interpolators_for_type( mRequestedMasterGlobalDofTypes( iFI )( 0 ), mtk::Master_Slave::MASTER );

        // get number of master FI bases and fields
        uint tDerNumBases  = tFI->get_number_of_space_time_bases();
        uint tDerNumFields = tFI->get_number_of_fields();

        // coefficients for dof type wrt which derivative is computed
        Matrix< DDRMat > tCoeff = tFI->get_coeff();

        // init dof counter
        uint tDofCounter = 0;

        // loop over the coefficients column
        for( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
        {
            // loop over the coefficients row
            for( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
            {
                // perturbation of the coefficent
                Matrix< DDRMat > tCoeffPert = tCoeff;
                tCoeffPert( iCoeffRow, iCoeffCol ) = tCoeffPert( iCoeffRow, iCoeffCol ) + aPerturbation * tCoeffPert( iCoeffRow, iCoeffCol );

                // setting the perturbed coefficients
                tFI->set_coeff( tCoeffPert );

                // reset properties, CM and SP for IWG
                this->reset_eval_flags();

                // evaluate the residual
                mSet->get_residual().fill( 0.0 );
                this->compute_residual( aWStar );

                Matrix< DDRMat > tResidual_Plus_Master
                =  mSet->get_residual()( { mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 1 ) }, { 0, 0 } );
                Matrix< DDRMat > tResidual_Plus_Slave
                =  mSet->get_residual()( { mSet->get_res_dof_assembly_map()( tDofIndexSlave  )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexSlave  )( 0, 1 ) }, { 0, 0 } );

                // perturbation of the coefficent
                tCoeffPert = tCoeff;
                tCoeffPert( iCoeffRow, iCoeffCol ) = tCoeffPert( iCoeffRow, iCoeffCol ) - aPerturbation * tCoeffPert( iCoeffRow, iCoeffCol );

                // setting the perturbed coefficients
                tFI->set_coeff( tCoeffPert );

                // reset properties, CM and SP for IWG
                this->reset_eval_flags();

                // evaluate the residual
                mSet->get_residual().fill( 0.0 );
                this->compute_residual( aWStar );

                Matrix< DDRMat > tResidual_Minus_Master
                      =  mSet->get_residual()( { mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 1 ) }, { 0, 0 } );
                Matrix< DDRMat > tResidual_Minus_Slave
                =  mSet->get_residual()( { mSet->get_res_dof_assembly_map()( tDofIndexSlave  )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexSlave  )( 0, 1 ) }, { 0, 0 } );

                // evaluate Jacobian
                aJacobiansFD( 0 )( iFI ).get_column( tDofCounter ) = ( tResidual_Plus_Master - tResidual_Minus_Master )/ ( 2.0 * aPerturbation * tCoeff( iCoeffRow, iCoeffCol ) );
                aJacobiansFD( 1 )( iFI ).get_column( tDofCounter ) = ( tResidual_Plus_Slave  - tResidual_Minus_Slave  )/ ( 2.0 * aPerturbation * tCoeff( iCoeffRow, iCoeffCol ) );

                // update dof counter
                tDofCounter++;
            }
        }
        // reset the coefficients values
        tFI->set_coeff( tCoeff );
    }

    // loop over the slave dof types
    for( uint iFI = 0; iFI < tSlaveNumDofType; iFI++ )
    {
        uint tDofIndexMaster = mSet->get_dof_index_for_type( mRequestedMasterGlobalDofTypes( iFI )( 0 ), mtk::Master_Slave::MASTER );
        uint tDofIndexSlave  = mSet->get_dof_index_for_type( mRequestedMasterGlobalDofTypes( iFI )( 0 ), mtk::Master_Slave::SLAVE  );

        uint tNumRowsMaster = mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 1 ) - mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 0 ) + 1;
        uint tNumColsMaster = mSet->get_jac_dof_assembly_map()( tDofIndexMaster )( tDofIndexMaster, 1 ) - mSet->get_jac_dof_assembly_map()( tDofIndexMaster )( tDofIndexMaster, 0 ) + 1;
        uint tNumRowsSlave  = mSet->get_res_dof_assembly_map()( tDofIndexSlave  )( 0 , 1 ) - mSet->get_res_dof_assembly_map()( tDofIndexSlave  )( 0 , 0 ) + 1;
        uint tNumColsSlave  = mSet->get_jac_dof_assembly_map()( tDofIndexSlave  )( tDofIndexSlave , 1 ) - mSet->get_jac_dof_assembly_map()( tDofIndexSlave  )( tDofIndexSlave , 0 ) + 1;

        aJacobiansFD( 0 )( tMasterNumDofType + iFI ).set_size( tNumRowsMaster, tNumColsMaster, 0.0 );
        aJacobiansFD( 1 )( tMasterNumDofType + iFI ).set_size( tNumRowsSlave, tNumColsSlave, 0.0 );

        Field_Interpolator * tFI = mFieldInterpolatorManager->get_field_interpolators_for_type( mRequestedMasterGlobalDofTypes( iFI )( 0 ), mtk::Master_Slave::SLAVE );

        // get number of master FI bases and fields
        uint tDerNumBases  = tFI->get_number_of_space_time_bases();
        uint tDerNumFields = tFI->get_number_of_fields();

        // coefficients for dof type wrt which derivative is computed
        Matrix< DDRMat > tCoeff = tFI->get_coeff();

        // init dof counter
        uint tDofCounter = 0;

        // loop over the coefficients columns
        for( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
        {
            // loop over the coefficients rows
            for( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
            {
                // perturbation of the coefficent
                Matrix< DDRMat > tCoeffPert = tCoeff;
                tCoeffPert( iCoeffRow, iCoeffCol ) = tCoeffPert( iCoeffRow, iCoeffCol ) + aPerturbation * tCoeffPert( iCoeffRow, iCoeffCol );

                // setting the perturbed coefficients
                tFI->set_coeff( tCoeffPert );

                // reset properties, CM and SP for IWG
                this->reset_eval_flags();

                // evaluate the residual
                mSet->get_residual().fill( 0.0 );
                this->compute_residual( aWStar );

                Matrix< DDRMat > tResidual_Plus_Master
                =  mSet->get_residual()( { mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 1 ) }, { 0, 0 } );
                Matrix< DDRMat > tResidual_Plus_Slave
                =  mSet->get_residual()( { mSet->get_res_dof_assembly_map()( tDofIndexSlave  )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexSlave  )( 0, 1 ) }, { 0, 0 } );

                // perturbation of the coefficent
                tCoeffPert = tCoeff;
                tCoeffPert( iCoeffRow, iCoeffCol ) = tCoeffPert( iCoeffRow, iCoeffCol ) - aPerturbation * tCoeffPert( iCoeffRow, iCoeffCol );

                // setting the perturbed coefficients
                tFI->set_coeff( tCoeffPert );

                // reset properties, CM and SP for IWG
                this->reset_eval_flags();

                // evaluate the residual
                mSet->get_residual().fill( 0.0 );
                this->compute_residual( aWStar );

                Matrix< DDRMat > tResidual_Minus_Master
                =  mSet->get_residual()( { mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 1 ) }, { 0, 0 } );
                Matrix< DDRMat > tResidual_Minus_Slave
                =  mSet->get_residual()( { mSet->get_res_dof_assembly_map()( tDofIndexSlave  )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexSlave  )( 0, 1 ) }, { 0, 0 } );

                // evaluate Jacobian
                aJacobiansFD( 0 )( tMasterNumDofType + iFI ).get_column( tDofCounter ) = ( tResidual_Plus_Master - tResidual_Minus_Master )/ ( 2.0 * aPerturbation * tCoeff( iCoeffRow, iCoeffCol ) );
                aJacobiansFD( 1 )( tMasterNumDofType + iFI ).get_column( tDofCounter ) = ( tResidual_Plus_Slave  - tResidual_Minus_Slave  )/ ( 2.0 * aPerturbation * tCoeff( iCoeffRow, iCoeffCol ) );

                // update dof counter
                tDofCounter++;
            }
        }
        // reset the coefficients values
        tFI->set_coeff( tCoeff );
    }
}

//------------------------------------------------------------------------------

        bool IWG::check_jacobian( real                                             aPerturbation,
                                  real                                             aEpsilon,
                                  real                                             aWStar,
                                  moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                  moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobiansFDs )
        {
            // compute jacobian with IWG
            this->compute_jacobian( aWStar );

            // compute jacobian by FD
            this->compute_jacobian_FD( aWStar, aPerturbation, aJacobiansFDs );

            // set size for comparison
            aJacobians.resize( aJacobiansFDs.size() );

            //define a boolean for check
            bool tCheckJacobian = true;

            // check each components
            for ( uint iJac = 0; iJac < aJacobiansFDs.size(); iJac++ )
            {
                // set size for comparison
                aJacobians( iJac ).resize( aJacobiansFDs( iJac ).size() );

                for( uint jJac = 0; jJac < aJacobiansFDs( iJac ).size(); jJac++ )
                {
                    uint tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
                    uint tIndexDep = mSet->get_dof_index_for_type( mRequestedMasterGlobalDofTypes( jJac )( 0 ), mtk::Master_Slave::MASTER );

                    aJacobians( iJac )( jJac )
                    = mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 ) },
                                            { mSet->get_jac_dof_assembly_map()( tDofIndex )( tIndexDep, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndex )( tIndexDep, 1 ) } );

                    for( uint iiJac = 0; iiJac < aJacobiansFDs( iJac )( jJac ).n_rows(); iiJac++ )
                    {
                        for( uint jjJac = 0; jjJac < aJacobiansFDs( iJac )( jJac ).n_cols(); jjJac++ )
                        {
                            tCheckJacobian = tCheckJacobian && ( aJacobians( iJac )( jJac )( iiJac, jjJac ) - aJacobiansFDs( iJac )( jJac )( iiJac, jjJac ) < aEpsilon );
                        }
                    }
                }
            }

            // return bool
            return tCheckJacobian;
        }

//------------------------------------------------------------------------------
        bool IWG::check_jacobian_double( real                                             aPerturbation,
                                         real                                             aEpsilon,
                                         real                                             aWStar,
                                         moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                         moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobiansFDs )
        {
            // compute jacobian with IWG
            this->compute_jacobian( aWStar );

            // compute jacobian by FD
            this->compute_jacobian_FD_double( aWStar, aPerturbation, aJacobiansFDs );

            // set jacobian size for comparison
            aJacobians.resize( 2 );

            //define a boolean for check
            bool tCheckJacobian = true;

            // check each components
            for ( uint iJac = 0; iJac < aJacobiansFDs.size(); iJac++ )
            {
                //set size for comparison
                aJacobians( iJac ).resize( aJacobiansFDs( iJac ).size() );

                uint tDofIndex;

                if( iJac == 0 )
                {
                    tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
                }
                else if( iJac == 1 )
                {
                    tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::SLAVE );
                }
                else { MORIS_ERROR( false, "only case 0 and 1 implemented" ); }

                // init dof counter
                uint tCounter = 0;

                // loop over master dof type
                for( uint jJac = 0; jJac < mRequestedMasterGlobalDofTypes.size(); jJac++ )
                {
                    // get set index for master dof type
                    uint tIndexDep = mSet->get_dof_index_for_type( mRequestedMasterGlobalDofTypes( jJac )( 0 ), mtk::Master_Slave::MASTER );

                    // fill jacobian matrix for comparison
                    aJacobians( iJac )( jJac )
                    = mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 ) },
                                            { mSet->get_jac_dof_assembly_map()( tDofIndex )( tIndexDep, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndex )( tIndexDep, 1 ) } );

                    // loop over the rows of jacobian
                    for( uint iiJac = 0; iiJac < aJacobiansFDs( iJac )( tCounter ).n_rows(); iiJac++ )
                    {
                        // loop over the columns of jacobian
                        for( uint jjJac = 0; jjJac < aJacobiansFDs( iJac )( tCounter ).n_cols(); jjJac++ )
                        {
                            // check component
                            tCheckJacobian = tCheckJacobian && ( aJacobians( iJac )( tCounter )( iiJac, jjJac ) - aJacobiansFDs( iJac )( tCounter )( iiJac, jjJac ) < aEpsilon );
                        }
                    }
                    // update dof counter
                    tCounter++;
                }

                // loop over slave dof types
                for( uint jJac = 0; jJac < mRequestedSlaveGlobalDofTypes.size(); jJac++ )
                {
                    // get set index for the slave dof type
                    uint tIndexDep = mSet->get_dof_index_for_type( mRequestedSlaveGlobalDofTypes( jJac )( 0 ), mtk::Master_Slave::SLAVE );

                    // fill jacobian matrix for comparison
                    aJacobians( iJac )( tCounter )
                    = mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 ) },
                                            { mSet->get_jac_dof_assembly_map()( tDofIndex )( tIndexDep, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndex )( tIndexDep, 1 ) } );

                    // loop over the rows of jacobian
                    for( uint iiJac = 0; iiJac < aJacobiansFDs( iJac )( tCounter ).n_rows(); iiJac++ )
                    {
                        // loop over the columns of jacobian
                        for( uint jjJac = 0; jjJac < aJacobiansFDs( iJac )( tCounter ).n_cols(); jjJac++ )
                        {
                            // check component
                            tCheckJacobian = tCheckJacobian && ( aJacobians( iJac )( tCounter )( iiJac, jjJac ) - aJacobiansFDs( iJac )( tCounter )( iiJac, jjJac ) < aEpsilon );
                        }
                    }
                    // update dof counter
                    tCounter++;
                }
            }

            // return bool
            return tCheckJacobian;
        }

//------------------------------------------------------------------------------

}   // end fem namespace
}   // end moris namespace


