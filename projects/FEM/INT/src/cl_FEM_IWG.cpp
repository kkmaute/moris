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

    for ( uint iDOF = 0; iDOF < mMasterDofTypes.size(); iDOF++ )
    {
        tCounter += mMasterDofTypes( iDOF ).size();
    }
    for ( uint iDOF = 0; iDOF < mSlaveDofTypes.size(); iDOF++ )
    {
        tCounter += mSlaveDofTypes( iDOF ).size();
    }

    for ( std::shared_ptr< Property > tProperty : mMasterProp )
    {
        moris::Cell< MSI::Dof_Type > tActiveDofType;
        tProperty->get_non_unique_dof_types( tActiveDofType );

        tCounter += tActiveDofType.size();
    }

    for ( std::shared_ptr< Property > tProperty : mSlaveProp )
    {
        moris::Cell< MSI::Dof_Type > tActiveDofType;
        tProperty->get_non_unique_dof_types( tActiveDofType );

        tCounter += tActiveDofType.size();
    }

    // get dof type from constitutive models
    for ( std::shared_ptr< Constitutive_Model > tCM : mMasterCM )
    {
        // get dof types for constitutive model
        moris::Cell< MSI::Dof_Type > tActiveDofType;
        tCM->get_non_unique_dof_types( tActiveDofType );

        tCounter += tActiveDofType.size();
    }

    // get dof type from constitutive models
    for ( std::shared_ptr< Constitutive_Model > tCM : mSlaveCM )
    {
        // get dof types for constitutive model
        moris::Cell< MSI::Dof_Type >  tActiveDofType;
        tCM->get_non_unique_dof_types( tActiveDofType );

        tCounter += tActiveDofType.size();
    }

    // get dof type from stabilization parameters
    for ( std::shared_ptr< Stabilization_Parameter > tSP : mStabilizationParam )
    {
        // get dof types for constitutive model
        moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType = tSP->get_global_dof_type_list( mtk::Master_Slave::MASTER );

        tCounter += tActiveDofType.size();
    }

    // get dof type from stabilization parameters
    for ( std::shared_ptr< Stabilization_Parameter > tSP : mStabilizationParam )
    {
        // get dof types for constitutive model
        moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType = tSP->get_global_dof_type_list( mtk::Master_Slave::SLAVE );

        tCounter += tActiveDofType.size();
    }

    aDofTypes.reserve( tCounter );

    for ( uint iDOF = 0; iDOF < mMasterDofTypes.size(); iDOF++ )
    {
        aDofTypes.append( mMasterDofTypes( iDOF ) );
    }
    for ( uint iDOF = 0; iDOF < mSlaveDofTypes.size(); iDOF++ )
    {
        aDofTypes.append( mSlaveDofTypes( iDOF )  );
    }

    for ( std::shared_ptr< Property > tProperty : mMasterProp )
    {
        moris::Cell< MSI::Dof_Type > tActiveDofType;
        tProperty->get_non_unique_dof_types( tActiveDofType );

        aDofTypes.append( tActiveDofType );
    }

    for ( std::shared_ptr< Property > tProperty : mSlaveProp )
    {
        moris::Cell< MSI::Dof_Type > tActiveDofType;
        tProperty->get_non_unique_dof_types( tActiveDofType );

        aDofTypes.append( tActiveDofType );
    }

    // get dof type from constitutive models
    for ( std::shared_ptr< Constitutive_Model > tCM : mMasterCM )
    {
        // get dof types for constitutive model
        moris::Cell< MSI::Dof_Type > tActiveDofType;
        tCM->get_non_unique_dof_types( tActiveDofType );

        aDofTypes.append( tActiveDofType );
    }

    // get dof type from constitutive models
    for ( std::shared_ptr< Constitutive_Model > tCM : mSlaveCM )
    {
        // get dof types for constitutive model
        moris::Cell< MSI::Dof_Type > tActiveDofType;
        tCM->get_non_unique_dof_types( tActiveDofType );

        aDofTypes.append( tActiveDofType );
    }

    // get dof type from stabilization parameters
    for ( std::shared_ptr< Stabilization_Parameter > tSP : mStabilizationParam )
    {
        // get dof types for constitutive model
        moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType = tSP->get_global_dof_type_list( mtk::Master_Slave::MASTER );

        // loop on property dof type
        for ( uint iDOF = 0; iDOF < tActiveDofType.size(); iDOF++ )
        {
            aDofTypes.append( tActiveDofType( iDOF ) );
        }
    }

    // get dof type from stabilization parameters
    for ( std::shared_ptr< Stabilization_Parameter > tSP : mStabilizationParam )
    {
        // get dof types for constitutive model
        moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType = tSP->get_global_dof_type_list( mtk::Master_Slave::SLAVE );

        // loop on property dof type
        for ( uint iDOF = 0; iDOF < tActiveDofType.size(); iDOF++ )
        {
            aDofTypes.append( tActiveDofType( iDOF ) );
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

        // set a size for the checkList (used to avoid repeating a dof type)
        Matrix< DDSMat > tCheckList( tNumDofTypes, 1, -1 );

        // get dof type from penalty parameter
        for ( uint iDOF = 0; iDOF < mMasterDofTypes.size(); iDOF++ )
        {
            sint tDofTypeindex = mSet->get_dof_index_for_type_1( mMasterDofTypes( iDOF )( 0 ) );  //FIXME'

            // put the dof type in the checklist
            tCheckList( tDofTypeindex ) = 1;

            // put the dof type in the global type list
            mMasterGlobalDofTypes.push_back( mMasterDofTypes( iDOF ) );
        }

        // get dof type from properties
        for ( std::shared_ptr< Property > tProperty : mMasterProp )
        {
            // get dof types for property
            moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType = tProperty->get_dof_type_list();

            // loop on property dof type
            for ( uint iDOF = 0; iDOF < tActiveDofType.size(); iDOF++ )
            {
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

        // get dof type from constitutive models
        for ( std::shared_ptr< Constitutive_Model > tCM : mMasterCM )
        {
            // get dof types for constitutive model
            moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType = tCM->get_global_dof_type_list();

            // loop on property dof type
            for ( uint iDOF = 0; iDOF < tActiveDofType.size(); iDOF++ )
            {
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

        // get dof type from stabilization parameters
        for ( std::shared_ptr< Stabilization_Parameter > tSP : mStabilizationParam )
        {
            // get dof types for constitutive model
            moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType = tSP->get_global_dof_type_list( mtk::Master_Slave::MASTER );

            // loop on property dof type
            for ( uint iDOF = 0; iDOF < tActiveDofType.size(); iDOF++ )
            {
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
        mMasterGlobalDofTypes.shrink_to_fit();

        // SLAVE--------------------------------------------------------

        // set size for the global dof type list
        mSlaveGlobalDofTypes.reserve( tNumDofTypes );

        // set a size for the checkList (used to avoid repeating a dof type)
        tCheckList.fill( -1 );

        // get dof type from penalty parameter
        for ( uint iDOF = 0; iDOF < mSlaveDofTypes.size(); iDOF++ )
        {
            sint tDofTypeindex = mSet->get_dof_index_for_type_1( mSlaveDofTypes( iDOF )( 0 ) );

            // put the dof type in the checklist
            tCheckList( tDofTypeindex ) = 1;

            // put the dof type in the global type list
            mSlaveGlobalDofTypes.push_back( mSlaveDofTypes( iDOF ) );
        }

        // get dof type from properties
        for ( std::shared_ptr< Property > tProperty : mSlaveProp )
        {
            // get dof types for property
            moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType = tProperty->get_dof_type_list();

            // loop on property dof type
            for ( uint iDOF = 0; iDOF < tActiveDofType.size(); iDOF++ )
            {
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

        // get dof type from constitutive models
        for ( std::shared_ptr< Constitutive_Model > tCM : mSlaveCM )
        {
            // get dof types for constitutive model
            moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType = tCM->get_global_dof_type_list();

            // loop on property dof type
            for ( uint iDOF = 0; iDOF < tActiveDofType.size(); iDOF++ )
            {
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

        // get dof type from stabilization parameters
        for ( std::shared_ptr< Stabilization_Parameter > tSP : mStabilizationParam )
        {
            // get dof types for constitutive model
            moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType = tSP->get_global_dof_type_list( mtk::Master_Slave::SLAVE );

            // loop on property dof type
            for ( uint iDOF = 0; iDOF < tActiveDofType.size(); iDOF++ )
            {
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
        mSlaveGlobalDofTypes.shrink_to_fit();

    }
//------------------------------------------------------------------------------

void IWG::build_requested_dof_type_list()
{
    Cell < enum MSI::Dof_Type > tRequestedDofTypes = mSet->get_requested_dof_types();

    mRequestedMasterGlobalDofTypes.clear();
    mRequestedSlaveGlobalDofTypes .clear();

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
    }

    mRequestedMasterGlobalDofTypes.shrink_to_fit();
    mRequestedSlaveGlobalDofTypes.shrink_to_fit();
}

//------------------------------------------------------------------------------

void IWG::set_residual_double( moris::Cell< Matrix< DDRMat > > & aResidual )
{
    // set the size of the residual cell
    aResidual.resize( 2 );

    Field_Interpolator * tFIMaster = mFieldInterpolatorManager->get_field_interpolators_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
    Field_Interpolator * tFISlave  = mFieldInterpolatorManager->get_field_interpolators_for_type( mResidualDofType( 0 ), mtk::Master_Slave::SLAVE );

    // set size of each residual matrix
    aResidual( 0 ).set_size( tFIMaster->get_number_of_space_time_coefficients(), 1, 0.0 );
    aResidual( 1 ).set_size( tFISlave ->get_number_of_space_time_coefficients(), 1, 0.0 );
}

//------------------------------------------------------------------------------

void IWG::set_jacobian( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians )
{

    // get number of dof types for the IWG
    uint tNumDofType = mSet->get_unique_dof_type_list().size();

    // set the size of the jacobian cell
    aJacobians.resize( 1 );
    aJacobians( 0 ).resize( tNumDofType );

    // get residual dof type number of dofs
    Field_Interpolator * tFI = mFieldInterpolatorManager->get_field_interpolators_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
    uint tNumResDof = tFI->get_number_of_space_time_coefficients();

    moris::Cell< enum MSI::Dof_Type > & tUniqueDofTypes = mSet->get_unique_dof_type_list();

    // loop over the master dof dependencies
    for( uint iDOF = 0; iDOF < tNumDofType; iDOF++ )
    {
        Field_Interpolator * tFI_2 = mFieldInterpolatorManager->get_field_interpolators_for_type( tUniqueDofTypes( iDOF ), mtk::Master_Slave::MASTER );

        if ( tFI_2 == nullptr )
        {
            continue;
        }

        uint tNumResDof_2 = tFI_2->get_number_of_space_time_coefficients();

        // set size for each jacobian matrix
        aJacobians( 0 )( iDOF ).set_size( tNumResDof, tNumResDof_2, 0.0 );
    }
}

//------------------------------------------------------------------------------

void IWG::set_jacobian_double( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians )
{
    uint tNumDofType = mSet->get_unique_dof_type_list().size();
    uint tNumDofTypeTotal = tNumDofType * 2 ;

    Field_Interpolator * tFIMaster = mFieldInterpolatorManager->get_field_interpolators_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
    Field_Interpolator * tFISlave  = mFieldInterpolatorManager->get_field_interpolators_for_type( mResidualDofType( 0 ), mtk::Master_Slave::SLAVE );

    uint tNumResDofMaster = tFIMaster->get_number_of_space_time_coefficients();
    uint tNumResDofSlave  = tFISlave->get_number_of_space_time_coefficients();

    // set the size of the jacobian cell
    aJacobians.resize( 2 );
    aJacobians( 0 ).resize( tNumDofTypeTotal );
    aJacobians( 1 ).resize( tNumDofTypeTotal );

    moris::Cell< enum MSI::Dof_Type > & tUniqueDofTypes = mSet->get_unique_dof_type_list();

    // loop over the master dof dependencies
    for( uint iMasterDOF = 0; iMasterDOF < tNumDofType; iMasterDOF++ )
    {
        Field_Interpolator * tFI_2 = mFieldInterpolatorManager->get_field_interpolators_for_type( tUniqueDofTypes( iMasterDOF ), mtk::Master_Slave::MASTER );

        if ( tFI_2 == nullptr )
        {
            continue;
        }

        uint tNumResDof = tFI_2->get_number_of_space_time_coefficients();

        // set size for each residual matrix
        aJacobians( 0 )( iMasterDOF ).set_size( tNumResDofMaster, tNumResDof, 0.0 );
        aJacobians( 1 )( iMasterDOF ).set_size( tNumResDofSlave,  tNumResDof, 0.0 );
    }

    // loop over the slave dof dependencies
    for( uint iSlaveDOF = 0; iSlaveDOF < tNumDofType; iSlaveDOF++ )
    {
        Field_Interpolator * tFI_2 = mFieldInterpolatorManager->get_field_interpolators_for_type( tUniqueDofTypes( iSlaveDOF ), mtk::Master_Slave::SLAVE );

        if ( tFI_2 == nullptr )
        {
            continue;
        }

        uint tNumResDof = tFI_2->get_number_of_space_time_coefficients();

        // set size for each residual matrix
        aJacobians( 0 )( tNumDofType + iSlaveDOF ).set_size( tNumResDofMaster, tNumResDof, 0.0 );
        aJacobians( 1 )( tNumDofType + iSlaveDOF ).set_size( tNumResDofSlave,  tNumResDof, 0.0 );
    }
}

//------------------------------------------------------------------------------

void IWG::set_dof_field_interpolators( mtk::Master_Slave aIsMaster )
{
    // set field interpolators for the SP
    for( std::shared_ptr< Stabilization_Parameter > tSP : this->get_stabilization_parameters() )
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
            tSPFIs( iDof ) = mFieldInterpolatorManager->get_field_interpolators_for_type( tSPDofTypes( iDof )( 0 ), aIsMaster );
        }

        // set the field interpolators for the SP
        tSP->set_dof_field_interpolators( tSPFIs );

        tSP->set_field_interpolator_manager( mFieldInterpolatorManager );

        tSP->set_set_pointer( mSet );
    }

    // set field interpolators for constitutive models
    for( std::shared_ptr< Constitutive_Model > tCM : this->get_constitutive_models( aIsMaster ) )
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

        tCM->set_field_interpolator_manager( mFieldInterpolatorManager );

        tCM->set_set_pointer( mSet );
    }

    // set field interpolators for properties
    for( std::shared_ptr< Property > tProp : this->get_properties( aIsMaster ) )
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

        tProp->set_field_interpolator_manager( mFieldInterpolatorManager );

        tProp->set_set_pointer( mSet );
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
        uint tDofIndex = mSet->get_dof_index_for_type( mRequestedMasterGlobalDofTypes( iFI )( 0 ), mtk::Master_Slave::MASTER );

        uint tNumRows = mSet->get_dof_assembly_map()( tDofIndex )( tDofIndex, 1 ) - mSet->get_dof_assembly_map()( tDofIndex )( tDofIndex, 0 ) + 1;
        uint tNumCols = mSet->get_dof_assembly_map()( tDofIndex )( tDofIndex, 3 ) - mSet->get_dof_assembly_map()( tDofIndex )( tDofIndex, 2 ) + 1;

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
                =  mSet->get_residual()( { mSet->get_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_dof_assembly_map()( tDofIndex )( 0, 1 ) }, { 0, 0 } );

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
                =  mSet->get_residual()( { mSet->get_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_dof_assembly_map()( tDofIndex )( 0, 1 ) }, { 0, 0 } );

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

    aJacobians.resize( 1 );

    //define a boolean for check
    bool tCheckJacobian = true;

    // check each components
    for ( uint iJac = 0; iJac < aJacobiansFDs.size(); iJac++ )
    {
        aJacobians( iJac ).resize( aJacobiansFDs( iJac ).size() );

        for( uint jJac = 0; jJac < aJacobiansFDs( iJac ).size(); jJac++ )
        {
            uint tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tIndexDep = mSet->get_dof_index_for_type( mRequestedMasterGlobalDofTypes( jJac )( 0 ), mtk::Master_Slave::MASTER );

             aJacobians( iJac )( jJac ) = mSet->get_jacobian()( { mSet->get_dof_assembly_map()( tDofIndex )( tIndexDep, 0 ), mSet->get_dof_assembly_map()( tDofIndex )( tIndexDep, 1 ) },
                                                                { mSet->get_dof_assembly_map()( tDofIndex )( tIndexDep, 2 ), mSet->get_dof_assembly_map()( tDofIndex )( tIndexDep, 3 ) } );

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


}   // end fem namespace
}   // end moris namespace


