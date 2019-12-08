/*
 * cl_FEM_Set.cpp
 *
 *  Created on: Apr 11, 2019
 *      Author: schmidt/noel
 */
#include <iostream>

#include "cl_MSI_Model_Solver_Interface.hpp" //FEM/MSI/src
#include "cl_MSI_Solver_Interface.hpp" //FEM/MSI/src
#include "cl_FEM_Set.hpp"                    //FEM/INT/src
#include "cl_FEM_Set_User_Info.hpp"                    //FEM/INT/src
#include "cl_FEM_Element_Factory.hpp"        //FEM/INT/src
#include "cl_FEM_Integrator.hpp"             //FEM/INT/src
#include "cl_FEM_Field_Interpolator_Manager.hpp"             //FEM/INT/src

#include "cl_MTK_Set.hpp"             //FEM/INT/src
#include "fn_equal_to.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
    Set::Set( moris::mtk::Set           * aMeshSet,
              fem::Set_User_Info        & aSetInfo,
              moris::Cell< Node_Base* > & aIPNodes )
    : mMeshSet( aMeshSet ),
      mNodes( aIPNodes ),
      mIWGs( aSetInfo.get_IWGs() ),
      mElementType( aSetInfo.get_set_type() )
    {
        for(  std::shared_ptr< IWG > tIWG : mIWGs )
        {
            tIWG->set_set_pointer( this );
        }
        // get mesh clusters on set
        mMeshClusterList = mMeshSet->get_clusters_on_set();

        // get number of mesh clusters on set
        uint tNumMeshClusters = mMeshClusterList.size();

        // set size for the equation objects list
        mEquationObjList.resize( tNumMeshClusters, nullptr);

        // create a fem cluster factory
        fem::Element_Factory tClusterFactory;

        // loop over mesh clusters on set
        for( luint iCluster = 0; iCluster < tNumMeshClusters; iCluster++ )
        {
            // create a fem cluster
            mEquationObjList( iCluster ) = tClusterFactory.create_cluster( mElementType,
                                                                           mMeshClusterList( iCluster ),
                                                                           mNodes,
                                                                           this );
        }

        // get spatial dimension
        mSpaceDim = mMeshSet->get_spatial_dim();

        // bool true is master IG cell are trivial
        mIsTrivialMaster = mMeshSet->is_trivial( mtk::Master_Slave::MASTER );

        // get interpolation geometry type
        mIPGeometryType = mMeshSet->get_interpolation_cell_geometry_type();

        // get integration geometry type
        mIGGeometryType = mMeshSet->get_integration_cell_geometry_type();

        // get space interpolation order for IP cells
        // FIXME if different for different fields
        mIPSpaceInterpolationOrder = mMeshSet->get_interpolation_cell_interpolation_order();

        // get space interpolation order for IG cells
        // FIXME if different for different fields
        mIGSpaceInterpolationOrder = mMeshSet->get_integration_cell_interpolation_order();

        // create geometry interpolation rule for IP cells
        Interpolation_Rule tIPGeometryInterpolationRule( mIPGeometryType,
                                                         Interpolation_Type::LAGRANGE,
                                                         mIPSpaceInterpolationOrder,
                                                         Interpolation_Type::LAGRANGE,
                                                         mtk::Interpolation_Order::LINEAR ); // FIXME not linear?

        // create geometry interpolation rule for IG cells
        Interpolation_Rule tIGGeometryInterpolationRule( mIGGeometryType,
                                                         Interpolation_Type::LAGRANGE,
                                                         mIGSpaceInterpolationOrder,
                                                         Interpolation_Type::LAGRANGE,
                                                         mtk::Interpolation_Order::LINEAR ); // FIXME not linear?

        // switch on set type
        switch ( mElementType )
        {
            // if block-set
            case ( fem::Element_Type::BULK ):
            {
                // create a geometry interpolator for IP cells
                mMasterIPGeometryInterpolator = new Geometry_Interpolator( tIPGeometryInterpolationRule, false );

                // create a geometry interpolator for IG cells
                mMasterIGGeometryInterpolator = new Geometry_Interpolator( tIGGeometryInterpolationRule, false );

                break;
            }

            // if side-set
            case( fem::Element_Type::SIDESET ):
            {
                // create a geometry interpolator for IP cells
                mMasterIPGeometryInterpolator = new Geometry_Interpolator( tIPGeometryInterpolationRule, true );

                // create a geometry interpolator for IG cells
                mMasterIGGeometryInterpolator = new Geometry_Interpolator( tIGGeometryInterpolationRule, true );

                break;
            }

            // if double side-set
            case( fem::Element_Type::DOUBLE_SIDESET ):
            {
                // bool true is slave IG cell are trivial
                mIsTrivialSlave = mMeshSet->is_trivial( mtk::Master_Slave::SLAVE );

                // create a geometry interpolator for master and slave IP cells
                mMasterIPGeometryInterpolator = new Geometry_Interpolator( tIPGeometryInterpolationRule, true );
                mSlaveIPGeometryInterpolator  = new Geometry_Interpolator( tIPGeometryInterpolationRule, true );

                // create a geometry interpolator for master and slave IG cells
                mMasterIGGeometryInterpolator = new Geometry_Interpolator( tIGGeometryInterpolationRule, true );
                mSlaveIGGeometryInterpolator  = new Geometry_Interpolator( tIGGeometryInterpolationRule, true );

                break;
            }

            // if none of the above
            default:
            {
                MORIS_ERROR(false, "Set::Set - unknown element type");
                break;
            }
        }

        // create a unique dof type list for solver
        this->create_unique_dof_type_list();

        this->create_dof_type_map_unique();

        // create a dof type list
        this->create_dof_type_list();

        // create a dof type map
        this->create_dof_type_map();

        // create an interpolation rule
        Integration_Rule tIntegrationRule = Integration_Rule( mIGGeometryType,
                                                              Integration_Type::GAUSS,
                                                              this->get_auto_integration_order( mIGGeometryType ),
                                                              Integration_Type::GAUSS,
                                                              Integration_Order::BAR_1 ); // fixme time order

        // create an integrator
        Integrator tIntegrator( tIntegrationRule );

        // get integration points
        tIntegrator.get_points( mIntegPoints );

        // get integration weights
        tIntegrator.get_weights( mIntegWeights );
    }

//------------------------------------------------------------------------------
    Set::~Set()
    {
        // delete the equation object pointers
        for( MSI::Equation_Object* tEquationObj : mEquationObjList )
        {
            delete tEquationObj;
        }
        mEquationObjList.clear();

        // delete the master interpolation geometry interpolator pointer
        if ( mMasterIPGeometryInterpolator != nullptr )
        {
            delete mMasterIPGeometryInterpolator;
        }

        // delete the slave interpolation geometry interpolator pointer
        if ( mSlaveIPGeometryInterpolator != nullptr )
        {
            delete mSlaveIPGeometryInterpolator;
        }

        // delete the master integration geometry interpolator pointer
        if ( mMasterIGGeometryInterpolator != nullptr )
        {
            delete mMasterIGGeometryInterpolator;
        }

        // delete the slave integration geometry interpolator pointer
        if ( mSlaveIGGeometryInterpolator != nullptr )
        {
            delete mSlaveIGGeometryInterpolator;
        }

        // delete the field interpolator pointers
        this->delete_pointers();

    }

//------------------------------------------------------------------------------
    void Set::delete_pointers()
    {
        //FIXME introduce Field_Interolator_Manager->delete_pointers();
    }

//------------------------------------------------------------------------------
    void Set::initialize_set( const bool aIsResidual )
    {
        if ( !mIsEmptySet )    //FIXME this flag is a hack. find better solution
        {
            this->create_residual_dof_assembly_map();

            this->create_dof_assembly_map( aIsResidual );

            this->create_requested_IWG_list();

            this->build_requested_IWG_dof_type_list( aIsResidual );

            for(  std::shared_ptr< IWG > tIWG : mRequestedIWGs )
            {
                tIWG->set_set_pointer( this );
            }
        }
    }

//------------------------------------------------------------------------------
    void Set::free_memory()
    {
        for(  std::shared_ptr< IWG > tIWG : mIWGs )
        {
            tIWG->free_memory( );
        }
    }

//------------------------------------------------------------------------------
    moris::Cell < enum MSI::Dof_Type > Set::get_requested_dof_types()
    {
        return mModelSolverInterface->get_solver_interface()
                                    ->get_requested_dof_types();
    }

//------------------------------------------------------------------------------
    moris::Cell< moris::Cell< enum MSI::Dof_Type > > Set::get_secundary_dof_types()
    {
        return mModelSolverInterface->get_solver_interface()
                                    ->get_secundary_dof_types();
    }

//------------------------------------------------------------------------------
    void Set::finalize( MSI::Model_Solver_Interface * aModelSolverInterface )
    {
        // delete the field interpolator pointers
        this->delete_pointers();

        // create the field interpolators
        this->create_field_interpolators( aModelSolverInterface );

        // set field interpolators for the IWGs
        this->set_IWG_field_interpolators();
    }

//------------------------------------------------------------------------------
    void Set::create_unique_dof_type_list()
    {
        // init dof type counter
        uint tCounter = 0;

        // loop over the IWGs
        for ( std::shared_ptr< IWG > tIWG : mIWGs )
        {
            // get an IWG non unique dof type list
            moris::Cell< MSI::Dof_Type >  tActiveDofType;
            tIWG->get_non_unique_dof_types( tActiveDofType );

            // update dof type counter
            tCounter += tActiveDofType.size();
        }

        // set max size for the dof type list
        mEqnObjDofTypeList.reserve( tCounter );

        // loop over the IWGs
        for ( std::shared_ptr< IWG > tIWG : mIWGs )
        {
            // get non unique dof type list
            moris::Cell< MSI::Dof_Type > tActiveDofType;
            tIWG->get_non_unique_dof_types( tActiveDofType );

            // populate the corresponding EqnObj dof type list
            mEqnObjDofTypeList.append( tActiveDofType );
        }

        // make the EqnObj dof type list unique
        std::sort( ( mEqnObjDofTypeList.data() ).data(),
                   ( mEqnObjDofTypeList.data() ).data() + mEqnObjDofTypeList.size());
        auto last = std::unique( ( mEqnObjDofTypeList.data() ).data(),
                                 ( mEqnObjDofTypeList.data() ).data() + mEqnObjDofTypeList.size() );
        auto pos  = std::distance( ( mEqnObjDofTypeList.data() ).data(), last );
        mEqnObjDofTypeList.resize( pos );
    }

//------------------------------------------------------------------------------
    void Set::create_dof_type_list()
    {
        uint tNumDofTypes = this->get_num_dof_types();

        // set size for the global dof type list
        mMasterDofTypes.reserve( tNumDofTypes );
        mSlaveDofTypes .reserve( tNumDofTypes );

        // create a list to check if dof type is already in the list
        Matrix< DDSMat > tMasterCheckList( tNumDofTypes, 1, -1 );
        Matrix< DDSMat > tSlaveCheckList ( tNumDofTypes, 1, -1 );

        // loop over the IWGs
        for ( std::shared_ptr< IWG > tIWG : mIWGs )
        {
            // get dof types for property
            moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypeMaster = tIWG->get_global_dof_type_list();

            // loop over the IWG active dof type
            for ( uint iDOF = 0; iDOF < tDofTypeMaster.size(); iDOF++ )
            {
                sint tDofTypeindex = this->get_dof_index_for_type_1( tDofTypeMaster( iDOF )( 0 ) );

                // if dof enum not in the list
                if ( tMasterCheckList( tDofTypeindex ) != 1 )
                {
                    // put the dof type in the checklist
                    tMasterCheckList( tDofTypeindex ) = 1;

                    // put the dof type in the global type list
                    mMasterDofTypes.push_back( tDofTypeMaster( iDOF ) );
                }
            }

            moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypeSlave = tIWG->get_global_dof_type_list( mtk::Master_Slave::SLAVE );

            // loop over the IWG active dof type
            for ( uint iDOF = 0; iDOF < tDofTypeSlave.size(); iDOF++ )
            {
                sint tDofTypeindex = this->get_dof_index_for_type_1( tDofTypeMaster( iDOF )( 0 ) );

                // if dof enum not in the list
                if ( tSlaveCheckList( tDofTypeindex ) != 1 )
                {
                    // put the dof type in the checklist
                    tSlaveCheckList( tDofTypeindex ) = 1;

                    // put the dof type in the global type list
                    mSlaveDofTypes.push_back( tDofTypeMaster( iDOF ) );
                }
            }
        }
        mMasterDofTypes.shrink_to_fit();
        mSlaveDofTypes .shrink_to_fit();
    }

//------------------------------------------------------------------------------
    void Set::create_dof_type_map()
    {
        // get number of master dof types
        uint tMasterNumDofs = this->get_number_of_field_interpolators();

        // get maximal dof type enum
        sint tMaxEnum = -1;

        // loop over the IWGs
        for ( uint iDOF = 0; iDOF < tMasterNumDofs; iDOF++ )
        {
            for ( uint Ik = 0; Ik < mMasterDofTypes( iDOF ).size(); Ik++ )
            {
                // get the highest dof type enum
                tMaxEnum = std::max( tMaxEnum, static_cast< int >( mMasterDofTypes( iDOF )( Ik ) ) );
            }
        }

        // get number of slave dof types
        uint tSlaveNumDofs =  this->get_number_of_field_interpolators( mtk::Master_Slave::SLAVE );

        // loop over the IWGs
        for ( uint iDOF = 0; iDOF < tSlaveNumDofs; iDOF++ )
        {
            for ( uint Ik = 0; Ik < mSlaveDofTypes( iDOF ).size(); Ik++ )
            {
                // get the highest dof type enum
                tMaxEnum = std::max( tMaxEnum, static_cast< int >( mSlaveDofTypes( iDOF )( Ik ) ) );
            }
        }
        // +1 since start at 0
        tMaxEnum++;

        MORIS_ASSERT( tMaxEnum != -1, "Set::create_dof_type_map(), no information to build dof type map" );

        // set size of dof type map    // FIXME replace with map
        mMasterDofTypeMap.set_size( tMaxEnum, 1, -1 );

        // loop over dof types
        for ( uint iDOF = 0; iDOF < tMasterNumDofs; iDOF++ )
        {
            mMasterDofTypeMap( static_cast< int >( mMasterDofTypes( iDOF )( 0 ) ), 0 ) = iDOF;
        }

        // set size of dof type map
        mSlaveDofTypeMap.set_size( tMaxEnum, 1, -1 );

        // loop over dof types
        for ( uint iDOF = 0; iDOF < tSlaveNumDofs; iDOF++ )
        {
            mSlaveDofTypeMap( static_cast< int >( mSlaveDofTypes( iDOF )( 0 ) ), 0 ) = iDOF;
        }
    }

//-----------------------------------------------------------------------------
    void Set::create_field_interpolators( MSI::Model_Solver_Interface * aModelSolverInterface )
    {
        mFieldInterpolatorManager = new Field_Interpolator_Manager( mMasterDofTypes,
                                                                    mSlaveDofTypes,
                                                                    this,
                                                                    aModelSolverInterface );

        mFieldInterpolatorManager->create_field_interpolators( aModelSolverInterface );
    }

//------------------------------------------------------------------------------
    void Set::set_IWG_field_interpolators()
    {
        // loop over the IWGs
        for ( std::shared_ptr< IWG > tIWG : mIWGs )
        {
            tIWG->set_field_interpolator_manager( mFieldInterpolatorManager );

            // set IWG field interpolators
            tIWG->set_dof_field_interpolators( mtk::Master_Slave::MASTER );

            // set IWG field interpolators
            tIWG->set_dof_field_interpolators( mtk::Master_Slave::SLAVE );
        }
    }

//------------------------------------------------------------------------------
    void Set::set_IWG_geometry_interpolators()
    {
        // loop over the IWGs
        for ( std::shared_ptr< IWG > tIWG : mIWGs )
        {
            //MASTER------------------------------------------------------------------------
            // set IWG geometry interpolators
            tIWG->set_geometry_interpolator( mMasterIPGeometryInterpolator );

            //SLAVE------------------------------------------------------------------------
            // set IWG field interpolators
            if( mElementType == fem::Element_Type::DOUBLE_SIDESET )
            {
                tIWG->set_geometry_interpolator( mSlaveIPGeometryInterpolator, mtk::Master_Slave::SLAVE );
            }
        }
    }

//------------------------------------------------------------------------------
    void Set::create_residual_dof_assembly_map()
    {
        // get list of requested dof types
        Cell < enum MSI::Dof_Type >  tRequestedDofTypes =  this->get_model_solver_interface()
                                                                               ->get_solver_interface()
                                                                               ->get_requested_dof_types();

        sint tMaxDofIndex = -1;

        for( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
        {
            sint tDofIndex = this->get_dof_index_for_type( tRequestedDofTypes( Ik ), mtk::Master_Slave::MASTER );

            if( tDofIndex != -1 )
            {
                tMaxDofIndex = std::max( tMaxDofIndex, tDofIndex );
            }

            tDofIndex = this->get_dof_index_for_type( tRequestedDofTypes( Ik ), mtk::Master_Slave::SLAVE );

            if( tDofIndex != -1 )
            {
                tMaxDofIndex = std::max( tMaxDofIndex, tDofIndex );
            }
        }
        tMaxDofIndex++;

        mResDofAssemblyMap.resize( tMaxDofIndex );

        for( uint Ik = 0; Ik < mResDofAssemblyMap.size(); Ik++ )
        {
            mResDofAssemblyMap( Ik ).set_size( 1, 2, -1 );
        }

        uint tCounter = 0;

        for( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
        {
            sint tDofIndex = this->get_dof_index_for_type( tRequestedDofTypes( Ik ), mtk::Master_Slave::MASTER );

            if( tDofIndex != -1 )
            {
                uint tNumCoeff = mFieldInterpolatorManager->get_field_interpolators_for_type( tRequestedDofTypes( Ik ), mtk::Master_Slave::MASTER )
                                                                          ->get_number_of_space_time_coefficients();

                mResDofAssemblyMap( tDofIndex )( 0, 0 ) = tCounter;
                mResDofAssemblyMap( tDofIndex )( 0, 1 ) = tCounter + tNumCoeff - 1;

                tCounter += tNumCoeff;
            }
        }

        for( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
        {
            sint tDofIndex = this->get_dof_index_for_type( tRequestedDofTypes( Ik ), mtk::Master_Slave::SLAVE  );

            if( tDofIndex != -1 )
            {
                uint tNumCoeff = mFieldInterpolatorManager->get_field_interpolators_for_type( tRequestedDofTypes( Ik ), mtk::Master_Slave::SLAVE )
                                                          ->get_number_of_space_time_coefficients();

                mResDofAssemblyMap( tDofIndex )( 0, 0 ) = tCounter;
                mResDofAssemblyMap( tDofIndex )( 0, 1 ) = tCounter + tNumCoeff - 1;

                tCounter += tNumCoeff;
            }
        }
    }

//------------------------------------------------------------------------------
    void Set::create_dof_assembly_map( const bool aIsResidual )
    {
        if( aIsResidual )
        {
            this->create_staggered_jacobian_dof_assembly_map();
        }
        else
        {
            this->create_jacobian_dof_assembly_map();
        }
    }

//------------------------------------------------------------------------------
    void Set::create_jacobian_dof_assembly_map()
    {
        // get list of requested dof types
        Cell < enum MSI::Dof_Type >  tRequestedDofTypes =  this->get_model_solver_interface()
                                                                               ->get_solver_interface()
                                                                               ->get_requested_dof_types();

        sint tMaxDofIndex = -1;

        for( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
        {
            sint tDofIndex = this->get_dof_index_for_type( tRequestedDofTypes( Ik ), mtk::Master_Slave::MASTER );

            if( tDofIndex != -1 )
            {
                tMaxDofIndex = std::max( tMaxDofIndex, tDofIndex );
            }

            tDofIndex = this->get_dof_index_for_type( tRequestedDofTypes( Ik ), mtk::Master_Slave::SLAVE );

            if( tDofIndex != -1 )
            {
                tMaxDofIndex = std::max( tMaxDofIndex, tDofIndex );
            }
        }

        tMaxDofIndex++;

        mJacDofAssemblyMap.resize( tMaxDofIndex );

        for( uint Ik = 0; Ik < mResDofAssemblyMap.size(); Ik++ )
        {
            mJacDofAssemblyMap( Ik ).set_size( tMaxDofIndex, 2, -1 );
        }

        for( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
        {
            sint tDofIndex = this->get_dof_index_for_type( tRequestedDofTypes( Ik ), mtk::Master_Slave::MASTER );

            if( tDofIndex != -1 )
            {
                uint tCounter_2 = 0;

                for( uint Ii = 0; Ii < tRequestedDofTypes.size(); Ii++ )
                {
                    sint tDofIndex_2 = this->get_dof_index_for_type( tRequestedDofTypes( Ii ), mtk::Master_Slave::MASTER );

                    if( tDofIndex_2 != -1 )
                    {
                        uint tNumCoeff_2 = mFieldInterpolatorManager->get_field_interpolators_for_type( tRequestedDofTypes( Ii ), mtk::Master_Slave::MASTER )
                                                                            ->get_number_of_space_time_coefficients();

                        mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 0 ) = tCounter_2;
                        mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 1 ) = tCounter_2 + tNumCoeff_2 - 1;

                        tCounter_2 += tNumCoeff_2;
                    }
                }

                for( uint Ii = 0; Ii < tRequestedDofTypes.size(); Ii++ )
                {
                    sint tDofIndex_2 = this->get_dof_index_for_type( tRequestedDofTypes( Ii ), mtk::Master_Slave::SLAVE );

                    if( tDofIndex_2 != -1 )
                    {
                        uint tNumCoeff_2 = mFieldInterpolatorManager->get_field_interpolators_for_type( tRequestedDofTypes( Ii ), mtk::Master_Slave::SLAVE )
                                                                            ->get_number_of_space_time_coefficients();

                        mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 0 ) = tCounter_2;
                        mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 1 ) = tCounter_2 + tNumCoeff_2 - 1;

                        tCounter_2 += tNumCoeff_2;
                    }
                }
            }
        }

        for( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
        {
            sint tDofIndex = this->get_dof_index_for_type( tRequestedDofTypes( Ik ), mtk::Master_Slave::SLAVE  );

            if( tDofIndex != -1 )
            {
                uint tCounter_2 = 0;

                for( uint Ii = 0; Ii < tRequestedDofTypes.size(); Ii++ )
                {
                    sint tDofIndex_2 = this->get_dof_index_for_type( tRequestedDofTypes( Ii ), mtk::Master_Slave::MASTER  );

                    if( tDofIndex_2 != -1 )
                    {
                        uint tNumCoeff_2 = mFieldInterpolatorManager->get_field_interpolators_for_type( tRequestedDofTypes( Ii ), mtk::Master_Slave::MASTER  )
                                                                    ->get_number_of_space_time_coefficients();

                        mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 0 ) = tCounter_2;
                        mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 1 ) = tCounter_2 + tNumCoeff_2 - 1;

                        tCounter_2 += tNumCoeff_2;
                    }
                }

                for( uint Ii = 0; Ii < tRequestedDofTypes.size(); Ii++ )
                {
                    sint tDofIndex_2 = this->get_dof_index_for_type( tRequestedDofTypes( Ii ), mtk::Master_Slave::SLAVE  );

                    if( tDofIndex_2 != -1 )
                    {
                        uint tNumCoeff_2 = mFieldInterpolatorManager->get_field_interpolators_for_type( tRequestedDofTypes( Ii ), mtk::Master_Slave::SLAVE  )
                                                                    ->get_number_of_space_time_coefficients();

                        mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 0 ) = tCounter_2;
                        mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 1 ) = tCounter_2 + tNumCoeff_2 - 1;

                        tCounter_2 += tNumCoeff_2;
                    }
                }
            }
        }
    }

//------------------------------------------------------------------------------
    void Set::create_staggered_jacobian_dof_assembly_map()
    {
        // get list of requested dof types
        Cell < enum MSI::Dof_Type >  tRequestedDofTypes =  this->get_model_solver_interface()
                                                                               ->get_solver_interface()
                                                                               ->get_requested_dof_types();

        Cell< Cell < enum MSI::Dof_Type > >  tSecundaryDofTypes =  this->get_model_solver_interface()
                                                                               ->get_solver_interface()
                                                                               ->get_secundary_dof_types();

        sint tMaxDofIndex = -1;

        for( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
        {
            sint tDofIndex = this->get_dof_index_for_type( tRequestedDofTypes( Ik ), mtk::Master_Slave::MASTER );

            if( tDofIndex != -1 )
            {
                tMaxDofIndex = std::max( tMaxDofIndex, tDofIndex );
            }

            tDofIndex = this->get_dof_index_for_type( tRequestedDofTypes( Ik ), mtk::Master_Slave::SLAVE );

            if( tDofIndex != -1 )
            {
                tMaxDofIndex = std::max( tMaxDofIndex, tDofIndex );
            }
        }
        tMaxDofIndex++;

        sint tMaxDofIndexSec = -1;

        for( uint Ik = 0; Ik < tSecundaryDofTypes.size(); Ik++ )
        {
            sint tDofIndex = this->get_dof_index_for_type( tSecundaryDofTypes( Ik )( 0 ), mtk::Master_Slave::MASTER );

            if( tDofIndex != -1 )
            {
                tMaxDofIndexSec = std::max( tMaxDofIndexSec, tDofIndex );
            }

            tDofIndex = this->get_dof_index_for_type( tSecundaryDofTypes( Ik )( 0 ), mtk::Master_Slave::SLAVE );

            if( tDofIndex != -1 )
            {
                tMaxDofIndexSec = std::max( tMaxDofIndexSec, tDofIndex );
            }
        }
        tMaxDofIndexSec++;

        mJacDofAssemblyMap.resize( tMaxDofIndex );

        for( uint Ik = 0; Ik < mResDofAssemblyMap.size(); Ik++ )
        {
            mJacDofAssemblyMap( Ik ).set_size( tMaxDofIndexSec, 2, -1 );
        }

        for( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
        {
            sint tDofIndex = this->get_dof_index_for_type( tRequestedDofTypes( Ik ), mtk::Master_Slave::MASTER );

            if( tDofIndex != -1 )
            {
                uint tCounter_2 = 0;

                for( uint Ii = 0; Ii < tSecundaryDofTypes.size(); Ii++ )
                {
                    sint tDofIndex_2 = this->get_dof_index_for_type( tSecundaryDofTypes( Ii )( 0 ), mtk::Master_Slave::MASTER );

                    if( tDofIndex_2 != -1 )
                    {
                        uint tNumCoeff_2 = mFieldInterpolatorManager->get_field_interpolators_for_type( tSecundaryDofTypes( Ii )( 0 ), mtk::Master_Slave::MASTER )
                                                                            ->get_number_of_space_time_coefficients();

                        mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 0 ) = tCounter_2;
                        mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 1 ) = tCounter_2 + tNumCoeff_2 - 1;

                        tCounter_2 += tNumCoeff_2;
                    }
                }

                for( uint Ii = 0; Ii < tSecundaryDofTypes.size(); Ii++ )
                {
                    sint tDofIndex_2 = this->get_dof_index_for_type(  tSecundaryDofTypes( Ii )( 0 ), mtk::Master_Slave::SLAVE );

                    if( tDofIndex_2 != -1 )
                    {
                        uint tNumCoeff_2 = mFieldInterpolatorManager->get_field_interpolators_for_type(  tSecundaryDofTypes( Ii )( 0 ), mtk::Master_Slave::SLAVE )
                                                                            ->get_number_of_space_time_coefficients();

                        mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 0 ) = tCounter_2;
                        mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 1 ) = tCounter_2 + tNumCoeff_2 - 1;

                        tCounter_2 += tNumCoeff_2;
                    }
                }
            }
        }

        for( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
        {
            sint tDofIndex = this->get_dof_index_for_type( tRequestedDofTypes( Ik ), mtk::Master_Slave::SLAVE  );

            if( tDofIndex != -1 )
            {
                uint tCounter_2 = 0;

                for( uint Ii = 0; Ii < tSecundaryDofTypes.size(); Ii++ )
                {
                    sint tDofIndex_2 = this->get_dof_index_for_type(  tSecundaryDofTypes( Ii )( 0 ), mtk::Master_Slave::MASTER  );

                    if( tDofIndex_2 != -1 )
                    {
                        uint tNumCoeff_2 = mFieldInterpolatorManager->get_field_interpolators_for_type(  tSecundaryDofTypes( Ii )( 0 ), mtk::Master_Slave::MASTER  )
                                                                    ->get_number_of_space_time_coefficients();

                        mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 0 ) = tCounter_2;
                        mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 1 ) = tCounter_2 + tNumCoeff_2 - 1;

                        tCounter_2 += tNumCoeff_2;
                    }
                }

                for( uint Ii = 0; Ii < tSecundaryDofTypes.size(); Ii++ )
                {
                    sint tDofIndex_2 = this->get_dof_index_for_type(  tSecundaryDofTypes( Ii )( 0 ), mtk::Master_Slave::SLAVE  );

                    if( tDofIndex_2 != -1 )
                    {
                        uint tNumCoeff_2 = mFieldInterpolatorManager->get_field_interpolators_for_type(  tSecundaryDofTypes( Ii )( 0 ), mtk::Master_Slave::SLAVE  )
                                                                    ->get_number_of_space_time_coefficients();

                        mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 0 ) = tCounter_2;
                        mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 1 ) = tCounter_2 + tNumCoeff_2 - 1;

                        tCounter_2 += tNumCoeff_2;
                    }
                }
            }
        }
    }

//------------------------------------------------------------------------------
    void Set::create_requested_IWG_list()
    {
        // get List of requested dof types
        Cell < enum MSI::Dof_Type >  tRequestedDofTypes =  this->get_model_solver_interface()
                                                               ->get_solver_interface()
                                                               ->get_requested_dof_types();

        mRequestedIWGs.clear();

        mRequestedIWGs.reserve( tRequestedDofTypes.size() );

        for( auto tDofType : tRequestedDofTypes )
        {
            for( uint Ik = 0; Ik < mIWGs.size(); Ik++ )
            {
                if( mIWGs( Ik )->get_residual_dof_type()( 0 ) == tDofType)
                {
                    mRequestedIWGs.push_back( mIWGs( Ik ) );

                    break;
                }
            }
        }

        mRequestedIWGs.shrink_to_fit();
    }

//------------------------------------------------------------------------------
    void Set::build_requested_IWG_dof_type_list( const bool aItResidual )
    {
        for( auto tIWG : mRequestedIWGs )
        {
            tIWG->build_requested_dof_type_list( aItResidual );
        }
    }

//------------------------------------------------------------------------------
    void Set::initialize_mJacobian()
    {
        if ( !mJacobianExist )
        {
            uint tNumCoeff = 0;

            moris::Cell < enum MSI::Dof_Type >tRequestedDofTypes = this->get_requested_dof_types();

            for( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
            {
                sint tDofIndex = this->get_dof_index_for_type( tRequestedDofTypes( Ik ), mtk::Master_Slave::MASTER );

                if( tDofIndex != -1 )
                {
                    tNumCoeff += mFieldInterpolatorManager->get_field_interpolators_for_type( tRequestedDofTypes( Ik ), mtk::Master_Slave::MASTER )
                                                                      ->get_number_of_space_time_coefficients();
                }

                tDofIndex = this->get_dof_index_for_type( tRequestedDofTypes( Ik ), mtk::Master_Slave::SLAVE  );

                if( tDofIndex != -1 )
                {
                    tNumCoeff += mFieldInterpolatorManager->get_field_interpolators_for_type( tRequestedDofTypes( Ik ), mtk::Master_Slave::SLAVE  )
                                                                                            ->get_number_of_space_time_coefficients();
                }
            }

            mJacobian.set_size( tNumCoeff, tNumCoeff, 0.0 );

            mJacobianExist = true;
        }
        else
        {
//            MORIS_ASSERT( mJacobian.numel() > 0, "Set::initialize_mJacobian() - Jacobian not properly initialized.");

            mJacobian.fill( 0.0 );
        }
    }

//------------------------------------------------------------------------------
    void Set::initialize_mResidual()
    {
        if ( !mResidualExist )
        {
            uint tNumCoeff = 0;

            moris::Cell < enum MSI::Dof_Type >tRequestedDofTypes = this->get_requested_dof_types();

            for( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
            {
                sint tDofIndex = this->get_dof_index_for_type( tRequestedDofTypes( Ik ), mtk::Master_Slave::MASTER );

                if( tDofIndex != -1 )
                {
                    tNumCoeff += mFieldInterpolatorManager->get_field_interpolators_for_type( tRequestedDofTypes( Ik ), mtk::Master_Slave::MASTER )
                                                                              ->get_number_of_space_time_coefficients();
                }

                tDofIndex = this->get_dof_index_for_type( tRequestedDofTypes( Ik ), mtk::Master_Slave::SLAVE  );

                if( tDofIndex != -1 )
                {
                    tNumCoeff += mFieldInterpolatorManager->get_field_interpolators_for_type( tRequestedDofTypes( Ik ), mtk::Master_Slave::SLAVE  )
                                                                                                    ->get_number_of_space_time_coefficients();
                }
            }

            mResidual.set_size( tNumCoeff, 1, 0.0 );

            mResidualExist = true;
        }
        else
        {
            //            MORIS_ASSERT( mResidual.numel() > 0, "Set::initialize_mResidual() - Residual not properly initialized.");

            mResidual.fill( 0.0 );
        }
    }

//------------------------------------------------------------------------------
    mtk::Interpolation_Order Set::get_auto_interpolation_order( const moris::uint        aNumVertices,
                                                                const mtk::Geometry_Type aGeometryType )
    {
        switch( aGeometryType )
        {
            case( mtk::Geometry_Type::LINE ) :
                switch( aNumVertices )
                {
                   case( 1 ) :
                       return mtk::Interpolation_Order::UNDEFINED;
                       break;

                   case( 2 ) :
                       return mtk::Interpolation_Order::LINEAR;
                       break;

                   case( 3 ) :
                       return mtk::Interpolation_Order::QUADRATIC;
                       break;

                   default :
                       MORIS_ERROR( false, " Element::get_auto_interpolation_order - not defined for LINE and number of vertices. ");
                       return mtk::Interpolation_Order::UNDEFINED;
                       break;
                }

            case( mtk::Geometry_Type::QUAD ) :
                switch( aNumVertices )
                {
                    case( 4 ) :
                        return mtk::Interpolation_Order::LINEAR;
                        break;

                    case( 8 ) :
                        return mtk::Interpolation_Order::SERENDIPITY;
                        break;

                    case( 9 ) :
                        return mtk::Interpolation_Order::QUADRATIC;
                        break;

                    case( 16 ) :
                        return mtk::Interpolation_Order::CUBIC;
                        break;

                    default :
                        MORIS_ERROR( false, " Element::get_auto_interpolation_order - not defined for QUAD and number of vertices. ");
                        return mtk::Interpolation_Order::UNDEFINED;
                        break;
                }

            case( mtk::Geometry_Type::HEX ) :
                switch( aNumVertices )
                {
                    case( 8 ) :
                        return mtk::Interpolation_Order::LINEAR;
                        break;

                    case( 20 ) :
                        return mtk::Interpolation_Order::SERENDIPITY;
                        break;

                    case( 27 ) :
                        return mtk::Interpolation_Order::QUADRATIC;
                        break;

                    case( 64 ) :
                        return mtk::Interpolation_Order::CUBIC;
                        break;

                    default :
                        MORIS_ERROR( false, " Element::get_auto_interpolation_order - not defined for HEX and number of vertices. ");
                        return mtk::Interpolation_Order::UNDEFINED;
                        break;
                }

                case( mtk::Geometry_Type::TET ) :
                switch( aNumVertices )
                {
                    case( 4 ) :
                        return mtk::Interpolation_Order::LINEAR;
                        break;

                    case( 10 ) :
                        return mtk::Interpolation_Order::QUADRATIC;
                        break;

                    case( 20 ) :
                        return mtk::Interpolation_Order::CUBIC;
                        break;

                    default :
                        MORIS_ERROR( false, " Element::get_auto_interpolation_order - not defined for TET and number of vertices. ");
                        return mtk::Interpolation_Order::UNDEFINED;
                        break;
                }

            default :
                MORIS_ERROR( false, " Element::get_auto_interpolation_order - not defined for this geometry type. ");
                return mtk::Interpolation_Order::UNDEFINED;
                break;
        }
    }

//------------------------------------------------------------------------------
    fem::Interpolation_Type Set::get_auto_time_interpolation_type( const moris::uint aNumVertices )
    {
        switch( aNumVertices )
        {
          case( 1 ) :
              return Interpolation_Type::CONSTANT;
              break;

          case( 2 ) :
          case( 3 ) :
          case( 4 ) :
              return Interpolation_Type::LAGRANGE;
              break;

          default :
              MORIS_ERROR( false, " Element::get_auto_time_interpolation_type - not defined this number of time vertices. ");
              return Interpolation_Type::UNDEFINED;
              break;
        }
    }

//------------------------------------------------------------------------------
    fem::Integration_Order Set::get_auto_integration_order( const mtk::Geometry_Type aGeometryType )
    {
        switch( aGeometryType )
        {
            case( mtk::Geometry_Type::LINE ) :
                return fem::Integration_Order::BAR_3;
                break;

            case( mtk::Geometry_Type::QUAD ) :
                 return fem::Integration_Order::QUAD_2x2;         //FIXME
                 break;

            case( mtk::Geometry_Type::HEX ) :
                return fem::Integration_Order::HEX_3x3x3;
                break;

            case( mtk::Geometry_Type::TRI ) :
                return fem::Integration_Order::TRI_6;
                break;

            case( mtk::Geometry_Type::TET ) :
                return fem::Integration_Order::TET_5;
                break;

            default :
                MORIS_ERROR( false, " Element::get_auto_integration_order - not defined for this geometry type. ");
                return Integration_Order::UNDEFINED;
                break;
        }
    }

//------------------------------------------------------------------------------
    moris::sint Set::get_dof_index_for_type_1( enum MSI::Dof_Type aDofType )
    {
        return mDofTypeMap( static_cast< int >( aDofType ), 0 );
    }

//------------------------------------------------------------------------------
    moris::uint Set::get_num_dof_types()
    {
        return this->get_unique_dof_type_list().size();
    }

//------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
