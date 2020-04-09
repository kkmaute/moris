/*
 * cl_FEM_Set.cpp
 *
 *  Created on: Apr 11, 2019
 *      Author: schmidt/noel
 */
#include <iostream>

#include "cl_MSI_Model_Solver_Interface.hpp"     //FEM/MSI/src
#include "cl_MSI_Solver_Interface.hpp"           //FEM/MSI/src
#include "cl_FEM_Model.hpp"                      //FEM/INT/src
#include "cl_FEM_Set.hpp"                        //FEM/INT/src
#include "cl_FEM_Set_User_Info.hpp"              //FEM/INT/src
#include "cl_FEM_Element_Factory.hpp"            //FEM/INT/src
#include "cl_FEM_Integrator.hpp"                 //FEM/INT/src
#include "cl_FEM_Field_Interpolator_Manager.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Element.hpp"      //FEM/INT/src
#include "cl_FEM_Cluster.hpp"                    //FEM/INT/src
#include "cl_MTK_Set.hpp"                        //FEM/INT/src
#include "fn_equal_to.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
    Set::Set(       fem::FEM_Model            * aFemModel,
                    moris::mtk::Set           * aMeshSet,
              const fem::Set_User_Info        & aSetInfo,
              const moris::Cell< Node_Base* > & aIPNodes ) : mFemModel( aFemModel ),
                                                             mMeshSet( aMeshSet ),
                                                             mNodes( aIPNodes ),
                                                             mIWGs( aSetInfo.get_IWGs() ),
                                                             mIQIs( aSetInfo.get_IQIs() )
    {
        // get the set type (BULK, SIDESET, DOUBLE_SIDESET)
        this->determine_set_type();

        // loop over the IWGs on the set
        for(  std::shared_ptr< IWG > tIWG : mIWGs )
        {
            // set the fem set pointer to the IWG
            tIWG->set_set_pointer( this );
        }

        // loop over the IQIs on the set
        for(  std::shared_ptr< IQI > tIQI : mIQIs )
        {
            // set the fem set pointer to the IQI
            tIQI->set_set_pointer( this );
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
            // init list of pointers to IP mesh cell
            moris::Cell< const mtk::Cell * > tInterpolationCell;

            // switch on set type
            switch ( mElementType )
            {
                // if bulk or sideset
                case( fem::Element_Type::BULK ):
                case( fem::Element_Type::SIDESET ):
                {
                    tInterpolationCell.resize( 1, &mMeshClusterList( iCluster )->get_interpolation_cell() );
                    break;
                }
                // if double sideset
                case( fem::Element_Type::DOUBLE_SIDESET ):
                {
                    tInterpolationCell.resize( 2 );
                    tInterpolationCell( 0 ) = &mMeshClusterList( iCluster )->get_interpolation_cell( mtk::Master_Slave::MASTER );
                    tInterpolationCell( 1 ) = &mMeshClusterList( iCluster )->get_interpolation_cell( mtk::Master_Slave::SLAVE );
                    break;
                }
                // if none of the above
                default:
                {
                    MORIS_ERROR(false, "Set::Set - unknown element type");
                    break;
                }
            }

            // create an interpolation element
            mEquationObjList( iCluster ) = new fem::Interpolation_Element( mElementType, tInterpolationCell, mNodes, this );

            // create a fem cluster
            std::shared_ptr< fem::Cluster > tCluster = std::make_shared< fem::Cluster >( mElementType,
                                                                                         mMeshClusterList( iCluster ),
                                                                                         this,
                                                                                         mEquationObjList( iCluster ) );
            // set the cluster to the interpolation element
            reinterpret_cast< fem::Interpolation_Element* >( mEquationObjList( iCluster ) )->set_cluster( tCluster, 0 );
        }

        // geometry and interpolation info----------------------------------------------
        //------------------------------------------------------------------------------
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

        // dof and dv dependencies info-------------------------------------------------
        //------------------------------------------------------------------------------
        // create a unique dof and dv type lists for solver
        this->create_unique_dof_and_dv_type_lists();

        // create a unique dof and dv type maps
        this->create_unique_dof_and_dv_type_maps();

        // create a dof and dv type lists
        this->create_dof_and_dv_type_lists();

        // create a dof and dv type maps
        this->create_dof_and_dv_type_maps();

        // integration info-------------------------------------------------------------
        //------------------------------------------------------------------------------
        // create an interpolation rule
        Integration_Rule tIntegrationRule = Integration_Rule( mIGGeometryType,
                                                              Integration_Type::GAUSS,
                                                              this->get_auto_integration_order( mIGGeometryType, mIPSpaceInterpolationOrder ),
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

        // delete the field interpolator pointers
        this->delete_pointers();
    }

//------------------------------------------------------------------------------
    void Set::initialize_set( const bool aIsResidual,
                              const bool aIsForward )
    {
        if ( !mIsEmptySet )    //FIXME this flag is a hack. find better solution
        {
            mIsResidual = aIsResidual;

            mIsForward  = aIsForward;

            this->create_residual_dof_assembly_map();

            this->create_dof_assembly_map( aIsResidual );

            this->create_dv_assembly_map();

            this->create_requested_IWG_list();

            this->create_requested_IQI_list();

            this->create_requested_IQI_type_map();

            this->build_requested_IWG_dof_type_list( aIsResidual );

            // set fem set pointer to IWGs
             for(  std::shared_ptr< IWG > tIWG : mRequestedIWGs )
             {
                 tIWG->set_set_pointer( this );
             }

             // set fem set pointer to IQIs
             for(  std::shared_ptr< IQI > tIQI : mRequestedIQIs )
             {
                 tIQI->set_set_pointer( this );
             }

//            if( !aIsForward )
//            {
//                this->create_requested_dv_assembly_map();
//            }
        }
    }

//------------------------------------------------------------------------------
    void Set::finalize( MSI::Model_Solver_Interface * aModelSolverInterface )
    {
        if ( !mIsEmptySet )    //FIXME this flag is a hack. find better solution
        {
            // delete the field interpolator pointers
            this->delete_pointers();

            // create the field interpolators
            this->create_field_interpolator_managers( aModelSolverInterface );

            // set field interpolator managers for the IWGs
            this->set_IWG_field_interpolator_managers();

            // set field interpolator managers for the IQIs
            this->set_IQI_field_interpolator_managers();
        }
    }

//------------------------------------------------------------------------------
    void Set::free_memory()
    {
        for(  std::shared_ptr< IWG > tIWG : mIWGs )
        {
            tIWG->free_memory();
        }
    }

//------------------------------------------------------------------------------
    void Set::delete_pointers()
    {
        if( mMasterFIManager != nullptr )
        {
            delete mMasterFIManager;
            mMasterFIManager = nullptr;
        }
        if( mSlaveFIManager != nullptr )
        {
            delete mSlaveFIManager;
            mSlaveFIManager = nullptr;
        }
    }

//------------------------------------------------------------------------------
    void Set::create_unique_dof_and_dv_type_lists()
    {
        // init dof and dv type counter
        uint tDofCounter = 0;
        uint tDvCounter  = 0;

        // loop over the IWGs
        for ( std::shared_ptr< IWG > tIWG : mIWGs )
        {
            // get an IWG non unique dof and dv types
            moris::Cell< MSI::Dof_Type >  tActiveDofType;
            moris::Cell< GEN_DV >         tActiveDvType;
            tIWG->get_non_unique_dof_and_dv_types( tActiveDofType, tActiveDvType );

            // update dof and dv type counters
            tDofCounter += tActiveDofType.size();
            tDvCounter  += tActiveDvType.size();
        }

        // loop over the IQIs
        for ( std::shared_ptr< IQI > tIQI : mIQIs )
        {
            // get an IWG non unique dof and dv types
            moris::Cell< MSI::Dof_Type >  tActiveDofType;
            moris::Cell< GEN_DV >   tActiveDvType;
            tIQI->get_non_unique_dof_and_dv_types( tActiveDofType, tActiveDvType );

            // update dof and dv type counter
            tDofCounter += tActiveDofType.size();
            tDvCounter  += tActiveDvType.size();
        }

        // set max size for the unique dof and dv type lists
        mUniqueDofTypeList.reserve( tDofCounter );
        mUniqueDvTypeList.reserve( tDvCounter );

        // loop over the IWGs
        for ( std::shared_ptr< IWG > tIWG : mIWGs )
        {
            // get non unique dof and dv types
            moris::Cell< MSI::Dof_Type > tActiveDofType;
            moris::Cell< GEN_DV >        tActiveDvType;
            tIWG->get_non_unique_dof_and_dv_types( tActiveDofType, tActiveDvType );

            // populate the corresponding unique dof and dv type lists
            mUniqueDofTypeList.append( tActiveDofType );
            mUniqueDvTypeList.append( tActiveDvType );
        }

        // loop over the IQIs
        for ( std::shared_ptr< IQI > tIQI : mIQIs )
        {
            // get non unique dof and dv types
            moris::Cell< MSI::Dof_Type > tActiveDofType;
            moris::Cell< GEN_DV >        tActiveDvType;
            tIQI->get_non_unique_dof_and_dv_types( tActiveDofType, tActiveDvType );

            // populate the corresponding unique dof and dv type lists
            mUniqueDofTypeList.append( tActiveDofType );
            mUniqueDvTypeList.append( tActiveDvType );
        }

        {
            // make the dof type list unique
            std::sort( ( mUniqueDofTypeList.data() ).data(),
                       ( mUniqueDofTypeList.data() ).data() + mUniqueDofTypeList.size());
            auto last = std::unique( ( mUniqueDofTypeList.data() ).data(),
                                     ( mUniqueDofTypeList.data() ).data() + mUniqueDofTypeList.size() );
            auto pos  = std::distance( ( mUniqueDofTypeList.data() ).data(), last );
            mUniqueDofTypeList.resize( pos );
        }

        {
            // make the dv type list unique
            std::sort( ( mUniqueDvTypeList.data() ).data(),
                       ( mUniqueDvTypeList.data() ).data() + mUniqueDvTypeList.size());
            auto last = std::unique( ( mUniqueDvTypeList.data() ).data(),
                                     ( mUniqueDvTypeList.data() ).data() + mUniqueDvTypeList.size() );
            auto pos  = std::distance( ( mUniqueDvTypeList.data() ).data(), last );
            mUniqueDvTypeList.resize( pos );
        }
    }

//------------------------------------------------------------------------------
    void Set::create_dof_and_dv_type_lists()
    {
        // get number of dof and dv types
        uint tNumDofTypes = this->get_num_unique_dof_types();
        uint tNumDvTypes  = this->get_num_unique_dv_types();

        // set size for the global dof type list
        mMasterDofTypes.reserve( tNumDofTypes );
        mSlaveDofTypes .reserve( tNumDofTypes );
        mMasterDvTypes.reserve( tNumDvTypes );
        mSlaveDvTypes .reserve( tNumDvTypes );

        // create a list to check if dof type is already in the list
        Matrix< DDSMat > tMasterCheckList( tNumDofTypes, 1, -1 );
        Matrix< DDSMat > tSlaveCheckList ( tNumDofTypes, 1, -1 );
        Matrix< DDSMat > tMasterDvCheckList( tNumDvTypes, 1, -1 );
        Matrix< DDSMat > tSlaveDvCheckList ( tNumDvTypes, 1, -1 );

        // loop over the IWGs
        for ( std::shared_ptr< IWG > tIWG : mIWGs )
        {
            // get master dof and dv types for the IWG
            moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypeMaster
            = tIWG->get_global_dof_type_list();
            moris::Cell< moris::Cell< GEN_DV > >  tDvTypeMaster
            = tIWG->get_global_dv_type_list();

            // loop over the IWG active master dof type
            for ( uint iDOF = 0; iDOF < tDofTypeMaster.size(); iDOF++ )
            {
                // get set index for the treated master dof type
                sint tDofTypeindex = this->get_index_from_unique_dof_type_map( tDofTypeMaster( iDOF )( 0 ) );

                // if dof enum not in the list
                if ( tMasterCheckList( tDofTypeindex ) != 1 )
                {
                    // put the dof type in the checklist
                    tMasterCheckList( tDofTypeindex ) = 1;

                    // put the dof type in the global type list
                    mMasterDofTypes.push_back( tDofTypeMaster( iDOF ) );
                }
            }

            // loop over the IWG active master dv type
            for ( uint iDv = 0; iDv < tDvTypeMaster.size(); iDv++ )
            {
                // get set index for the treated master dof type
                sint tDvTypeindex = this->get_index_from_unique_dv_type_map( tDvTypeMaster( iDv )( 0 ) );

                // if dv enum not in the list
                if ( tMasterDvCheckList( tDvTypeindex ) != 1 )
                {
                    // put the dof type in the checklist
                    tMasterDvCheckList( tDvTypeindex ) = 1;

                    // put the dof type in the global type list
                    mMasterDvTypes.push_back( tDvTypeMaster( iDv ) );
                }
            }

            // get slave dof and dv types for the IWG
            moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypeSlave
            = tIWG->get_global_dof_type_list( mtk::Master_Slave::SLAVE );
            moris::Cell< moris::Cell< GEN_DV > >  tDvTypeSlave
            = tIWG->get_global_dv_type_list( mtk::Master_Slave::SLAVE );

            // loop over the IWG active slave dof type
            for ( uint iDOF = 0; iDOF < tDofTypeSlave.size(); iDOF++ )
            {
                // get set index for the treated slave dof type
                sint tDofTypeindex = this->get_index_from_unique_dof_type_map( tDofTypeSlave( iDOF )( 0 ) );

                // if dof enum not in the list
                if ( tSlaveCheckList( tDofTypeindex ) != 1 )
                {
                    // put the dof type in the checklist
                    tSlaveCheckList( tDofTypeindex ) = 1;

                    // put the dof type in the global type list
                    mSlaveDofTypes.push_back( tDofTypeSlave( iDOF ) );
                }
            }

            // loop over the IWG active slave dv type
            for ( uint iDv = 0; iDv < tDvTypeSlave.size(); iDv++ )
            {
                // get set index for the treated slave dv type
                sint tDvTypeindex = this->get_index_from_unique_dv_type_map( tDvTypeSlave( iDv )( 0 ) );

                // if dv enum not in the list
                if ( tSlaveDvCheckList( tDvTypeindex ) != 1 )
                {
                    // put the dv type in the checklist
                    tSlaveDvCheckList( tDvTypeindex ) = 1;

                    // put the dv type in the global type list
                    mSlaveDvTypes.push_back( tDvTypeSlave( iDv ) );
                }
            }
        }

        // loop over the IQIs
        for ( std::shared_ptr< IQI > tIQI : mIQIs )
        {
            // get master dof and dv types for the IWG
            moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypeMaster
            = tIQI->get_global_dof_type_list();
            moris::Cell< moris::Cell< GEN_DV > >  tDvTypeMaster
            = tIQI->get_global_dv_type_list();

            // loop over the IQI active master dof type
            for ( uint iDOF = 0; iDOF < tDofTypeMaster.size(); iDOF++ )
            {
                // get set index for the treated master dof type
                sint tDofTypeindex = this->get_index_from_unique_dof_type_map( tDofTypeMaster( iDOF )( 0 ) );

                // if dof enum not in the list
                if ( tMasterCheckList( tDofTypeindex ) != 1 )
                {
                    // put the dof type in the checklist
                    tMasterCheckList( tDofTypeindex ) = 1;

                    // put the dof type in the global type list
                    mMasterDofTypes.push_back( tDofTypeMaster( iDOF ) );
                }
            }

            // loop over the IQI active master dv type
            for ( uint iDv = 0; iDv < tDvTypeMaster.size(); iDv++ )
            {
                // get set index for the treated master dv type
                sint tDvTypeindex = this->get_index_from_unique_dv_type_map( tDvTypeMaster( iDv )( 0 ) );

                // if dv enum not in the list
                if ( tMasterDvCheckList( tDvTypeindex ) != 1 )
                {
                    // put the dv type in the checklist
                    tMasterDvCheckList( tDvTypeindex ) = 1;

                    // put the dv type in the global type list
                    mMasterDvTypes.push_back( tDvTypeMaster( iDv ) );
                }
            }

            // get slave dof and dv types for the IWG
            moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypeSlave
            = tIQI->get_global_dof_type_list( mtk::Master_Slave::SLAVE );
            moris::Cell< moris::Cell< GEN_DV > >  tDvTypeSlave
            = tIQI->get_global_dv_type_list( mtk::Master_Slave::SLAVE );

            // loop over the IWG active slave dof type
            for ( uint iDOF = 0; iDOF < tDofTypeSlave.size(); iDOF++ )
            {
                // get set index for the treated slave dof type
                sint tDofTypeindex = this->get_index_from_unique_dof_type_map( tDofTypeSlave( iDOF )( 0 ) );

                // if dof enum not in the list
                if ( tSlaveCheckList( tDofTypeindex ) != 1 )
                {
                    // put the dof type in the checklist
                    tSlaveCheckList( tDofTypeindex ) = 1;

                    // put the dof type in the global type list
                    mSlaveDofTypes.push_back( tDofTypeSlave( iDOF ) );
                }
            }

            // loop over the IWG active slave dv type
            for ( uint iDv = 0; iDv < tDvTypeSlave.size(); iDv++ )
            {
                // get set index for the treated slave dv type
                sint tDvTypeindex = this->get_index_from_unique_dv_type_map( tDvTypeSlave( iDv )( 0 ) );

                // if dv enum not in the list
                if ( tSlaveDvCheckList( tDvTypeindex ) != 1 )
                {
                    // put the dv type in the checklist
                    tSlaveDvCheckList( tDvTypeindex ) = 1;

                    // put the dv type in the global type list
                    mSlaveDvTypes.push_back( tDvTypeSlave( iDv ) );
                }
            }
        }

        // shrink list to fit to number of unique dof and dv types
        mMasterDofTypes.shrink_to_fit();
        mSlaveDofTypes .shrink_to_fit();
        mMasterDvTypes.shrink_to_fit();
        mSlaveDvTypes .shrink_to_fit();
    }

//------------------------------------------------------------------------------
    void Set::create_dof_and_dv_type_maps()
    {
        // get number of master dof types
        uint tMasterNumDofs = this->get_dof_type_list().size();

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
        uint tSlaveNumDofs = this->get_dof_type_list( mtk::Master_Slave::SLAVE ).size();

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

        MORIS_ASSERT( tMaxEnum != -1, "Set::create_dof_and_dv_type_maps(), no information to build dof type map" );

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

        // get number of master dv types
        uint tMasterNumDvs  = this->get_dv_type_list().size();

        // get maximal dv type enum
        tMaxEnum = -1;

        // loop over the dv types
        for ( uint iDv = 0; iDv < tMasterNumDvs; iDv++ )
        {
            for ( uint Ik = 0; Ik < mMasterDvTypes( iDv ).size(); Ik++ )
            {
                // get the highest dof type enum
                tMaxEnum = std::max( tMaxEnum, static_cast< int >( mMasterDvTypes( iDv )( Ik ) ) );
            }
        }

        // get number of slave dv types
        uint tSlaveNumDvs =  this->get_dv_type_list( mtk::Master_Slave::SLAVE ).size();

        // loop over the IWGs
        for ( uint iDOF = 0; iDOF < tSlaveNumDvs; iDOF++ )
        {
            for ( uint Ik = 0; Ik < mSlaveDvTypes( iDOF ).size(); Ik++ )
            {
                // get the highest dof type enum
                tMaxEnum = std::max( tMaxEnum, static_cast< int >( mSlaveDvTypes( iDOF )( Ik ) ) );
            }
        }
        // +1 since start at 0
        tMaxEnum++;

        MORIS_ASSERT( tMaxEnum != -1, "Set::create_dv_type_map(), no information to build dv type map" );

        // set size of dv type map    // FIXME replace with map
        mMasterDvTypeMap.set_size( tMaxEnum, 1, -1 );

        // loop over dv types
        for ( uint iDv = 0; iDv < tMasterNumDvs; iDv++ )
        {
            mMasterDvTypeMap( static_cast< int >( mMasterDvTypes( iDv )( 0 ) ), 0 ) = iDv;
        }

        // set size of dv type map
        mSlaveDvTypeMap.set_size( tMaxEnum, 1, -1 );

        // loop over dv types
        for ( uint iDv = 0; iDv < tSlaveNumDvs; iDv++ )
        {
            mSlaveDvTypeMap( static_cast< int >( mSlaveDvTypes( iDv )( 0 ) ), 0 ) = iDv;
        }
    }

//-----------------------------------------------------------------------------
    void Set::create_field_interpolator_managers( MSI::Model_Solver_Interface * aModelSolverInterface )
    {
        // create the master field interpolator manager
        mMasterFIManager = new Field_Interpolator_Manager( mMasterDofTypes,
                                                           mMasterDvTypes,
                                                           this );

        // create the geometry interpolators on the master FI manager
        mMasterFIManager->create_geometry_interpolators();

        // create the field interpolators on the master FI manager
        mMasterFIManager->create_field_interpolators( aModelSolverInterface );

        // create the slave field interpolator manager
        mSlaveFIManager = new Field_Interpolator_Manager( mSlaveDofTypes,
                                                          mSlaveDvTypes,
                                                          this,
                                                          mtk::Master_Slave::SLAVE );

        // create the geometry interpolators on the slave FI manager
        mSlaveFIManager->create_geometry_interpolators();

        // create the field interpolators on the slave FI manager
        mSlaveFIManager->create_field_interpolators( aModelSolverInterface );

    }

//------------------------------------------------------------------------------
    void Set::set_IWG_field_interpolator_managers()
    {
        // loop over the IWGs
        for ( std::shared_ptr< IWG > tIWG : mIWGs )
        {
            // set the master FI manager
            tIWG->set_field_interpolator_manager( mMasterFIManager );

            // if double sideset, set slave
            if( mElementType == fem::Element_Type::DOUBLE_SIDESET )
            {
                // set IWG slave field interpolator manager
                tIWG->set_field_interpolator_manager( mSlaveFIManager, mtk::Master_Slave::SLAVE );
            }
        }
    }

//------------------------------------------------------------------------------
    void Set::set_IWG_cluster_for_stabilization_parameters( fem::Cluster * aCluster )
    {
        // loop over the IWGs
        for( auto tIWG : mIWGs)
        {
            // get the SP from the IWG
            moris::Cell< std::shared_ptr< Stabilization_Parameter > > & tSPs
            = tIWG->get_stabilization_parameters();

            // loop over the SP
            for( auto tSP : tSPs )
            {
                // check if SP is null
                if( tSP != nullptr )
                {
                    // set the fem cluster
                    tSP->set_cluster( aCluster );
                }
            }
        }
    }

//------------------------------------------------------------------------------
    void Set::set_IQI_cluster_for_stabilization_parameters( fem::Cluster * aCluster )
    {
        // loop over the IQIs
        for( auto tIQI : mIQIs)
        {
            // get the SP from the IQI
            moris::Cell< std::shared_ptr< Stabilization_Parameter > > & tSPs
            = tIQI->get_stabilization_parameters();

            // loop over the SPs
            for( auto tSP : tSPs )
            {
                // check if SP is null
                if( tSP != nullptr )
                {
                    // set the fem cluster
                    tSP->set_cluster( aCluster );
                }
            }
        }
    }


//------------------------------------------------------------------------------
    void Set::set_IQI_field_interpolator_managers()
    {
        // loop over the IQIs
        for ( std::shared_ptr< IQI > tIQI : mIQIs )
        {
            // set IQI master FI manager
            tIQI->set_field_interpolator_manager( mMasterFIManager );

            // if double sideset, set slave
            if( mElementType == fem::Element_Type::DOUBLE_SIDESET )
            {
                // set IQI slave FI manager
                tIQI->set_field_interpolator_manager( mSlaveFIManager, mtk::Master_Slave::SLAVE );
            }
        }
    }

//------------------------------------------------------------------------------
    void Set::create_residual_dof_assembly_map()
    {
        // get the list of requested dof types by the solver
        moris::Cell < enum MSI::Dof_Type > tRequestedDofTypes
        =  this->get_model_solver_interface()
               ->get_solver_interface()
               ->get_requested_dof_types();

        // init the max index for dof types
        sint tMaxDofIndex = -1;

        // loop over the requested dof types
        for( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
        {
            // get the set index for the requested master dof type
            sint tDofIndex = this->get_dof_index_for_type( tRequestedDofTypes( Ik ),
                                                           mtk::Master_Slave::MASTER );

            // if the index was set (and is different from -1)
            if( tDofIndex != -1 )
            {
                // update the max index for dof type
                tMaxDofIndex = std::max( tMaxDofIndex, tDofIndex );
            }

            // get the set index for the requested slave dof type
            tDofIndex = this->get_dof_index_for_type( tRequestedDofTypes( Ik ),
                                                      mtk::Master_Slave::SLAVE );

            // if the index was set (and is different -1)
            if( tDofIndex != -1 )
            {
                // update the max index for dof type
                tMaxDofIndex = std::max( tMaxDofIndex, tDofIndex );
            }
        }
        // add +1 to the max index for dof type (since 0 based)
        tMaxDofIndex++;

        // set size for the residual assembly map
        mResDofAssemblyMap.resize( tMaxDofIndex );

        // init the residual assembly map
        for( uint Ik = 0; Ik < mResDofAssemblyMap.size(); Ik++ )
        {
            mResDofAssemblyMap( Ik ).set_size( 1, 2, -1 );
        }

        // init dof coefficients counter
        uint tCounter = 0;

        // loop over the requested dof types
        for( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
        {
            // get the set index for the requested master dof type
            sint tDofIndex = this->get_dof_index_for_type( tRequestedDofTypes( Ik ),
                                                           mtk::Master_Slave::MASTER );

            // if the index was set (and is different from -1)
            if( tDofIndex != -1 )
            {
                // get the number of coefficients re;ated to the master dof type
                uint tNumCoeff = mMasterFIManager->get_field_interpolators_for_type( tRequestedDofTypes( Ik ) )
                                                 ->get_number_of_space_time_coefficients();

                // fill the residual assembly map with starting and ending indices for the master dof type
                mResDofAssemblyMap( tDofIndex )( 0, 0 ) = tCounter;
                mResDofAssemblyMap( tDofIndex )( 0, 1 ) = tCounter + tNumCoeff - 1;

                // update the dof coefficient counter
                tCounter += tNumCoeff;
            }
        }

        // loop over the requested dof types
        for( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
        {
            //get the set index for the requested slave dof type
            sint tDofIndex = this->get_dof_index_for_type( tRequestedDofTypes( Ik ),
                                                           mtk::Master_Slave::SLAVE );

            // if the dof type was set (its set index is different from -1)
            if( tDofIndex != -1 )
            {
                // get the number of coefficients for the slave dof type
                uint tNumCoeff = mSlaveFIManager->get_field_interpolators_for_type( tRequestedDofTypes( Ik ) )
                                                ->get_number_of_space_time_coefficients();

                // fill the residual assembly map with starting and ending indices for the slave dof type
                mResDofAssemblyMap( tDofIndex )( 0, 0 ) = tCounter;
                mResDofAssemblyMap( tDofIndex )( 0, 1 ) = tCounter + tNumCoeff - 1;

                // update the dof coefficient counter
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
        // get list of requested dof types (by the solver)
        moris::Cell < enum MSI::Dof_Type >  tRequestedDofTypes
        =  this->get_model_solver_interface()
               ->get_solver_interface()
               ->get_requested_dof_types();

        sint tMaxDofIndex = -1;

        // master
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
                        uint tNumCoeff_2 = mMasterFIManager->get_field_interpolators_for_type( tRequestedDofTypes( Ii ) )
                                                           ->get_number_of_space_time_coefficients();

                        mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 0 ) = tCounter_2;
                        mJacDofAssemblyMap( tDofIndex )( tDofIndex_2, 1 ) = tCounter_2 + tNumCoeff_2 - 1;

                        tCounter_2 += tNumCoeff_2;
                    }
                }

                // slave
                for( uint Ii = 0; Ii < tRequestedDofTypes.size(); Ii++ )
                {
                    sint tDofIndex_2 = this->get_dof_index_for_type( tRequestedDofTypes( Ii ), mtk::Master_Slave::SLAVE );

                    if( tDofIndex_2 != -1 )
                    {
                        uint tNumCoeff_2 = mSlaveFIManager->get_field_interpolators_for_type( tRequestedDofTypes( Ii ) )
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
                        uint tNumCoeff_2 = mMasterFIManager->get_field_interpolators_for_type( tRequestedDofTypes( Ii ) )
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
                        uint tNumCoeff_2 = mSlaveFIManager->get_field_interpolators_for_type( tRequestedDofTypes( Ii ) )
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
        moris::Cell < enum MSI::Dof_Type > tRequestedDofTypes
        = this->get_model_solver_interface()
              ->get_solver_interface()
              ->get_requested_dof_types();

        moris::Cell< moris::Cell < enum MSI::Dof_Type > > tSecundaryDofTypes
        = this->get_model_solver_interface()
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
                        uint tNumCoeff_2 = mMasterFIManager->get_field_interpolators_for_type( tSecundaryDofTypes( Ii )( 0 ) )
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
                        uint tNumCoeff_2 = mSlaveFIManager->get_field_interpolators_for_type(  tSecundaryDofTypes( Ii )( 0 ) )
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
                        uint tNumCoeff_2 = mMasterFIManager->get_field_interpolators_for_type(  tSecundaryDofTypes( Ii )( 0 ) )
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
                        uint tNumCoeff_2 = mSlaveFIManager->get_field_interpolators_for_type(  tSecundaryDofTypes( Ii )( 0 ) )
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
    void Set::create_dv_assembly_map()
    {
        // get the list of requested dv types by the opt solver
        // FIXME for now everything is evaluated
        //Cell < enum GEN_DV >  tRequestedDvTypes;

        // init the max index for dv types
        sint tMaxDvIndex = -1;

        // loop over the dv types
        for( uint Ik = 0; Ik < mMasterDvTypes.size(); Ik++ )
        {
            // get the set index for the requested master dof type
            sint tDvIndex = this->get_dv_index_for_type( mMasterDvTypes( Ik )( 0 ),
                                                         mtk::Master_Slave::MASTER );

            // if the index was set (and is different from -1)
            if( tDvIndex != -1 )
            {
                // update the max index for dv type
                tMaxDvIndex = std::max( tMaxDvIndex, tDvIndex );
            }

            // get the set index for the requested slave slave type
            tDvIndex = this->get_dv_index_for_type( mSlaveDvTypes( Ik )( 0 ),
                                                    mtk::Master_Slave::SLAVE );

            // if the index was set (and is different -1)
            if( tDvIndex != -1 )
            {
                // update the max index for dv type
                tMaxDvIndex = std::max( tMaxDvIndex, tDvIndex );
            }
        }
        // add +1 to the max index for dv type (since 0 based)
        tMaxDvIndex++;

        // set size for the dv assembly map
        mDvAssemblyMap.resize( tMaxDvIndex );

        // init the dv assembly map
        for( uint Ik = 0; Ik < mDvAssemblyMap.size(); Ik++ )
        {
            mDvAssemblyMap( Ik ).set_size( 1, 2, -1 );
        }

        // init dv coefficients counter
        uint tCounter = 0;

        // loop over the dv types
        for( uint Ik = 0; Ik < mMasterDvTypes.size(); Ik++ )
        {
            // get the set index for the requested master dv type
            sint tDvIndex = this->get_dv_index_for_type( mMasterDvTypes( Ik )( 0 ),
                                                         mtk::Master_Slave::MASTER );

            // if the index was set (and is different from -1)
            if( tDvIndex != -1 )
            {
                // get the number of coefficients related to the master dv type
                uint tNumCoeff = mMasterFIManager->get_field_interpolators_for_type( mMasterDvTypes( Ik )( 0 ) )
                                                 ->get_number_of_space_time_coefficients();

                // fill the dv assembly map with starting and ending indices
                // for the master dv type
                mDvAssemblyMap( tDvIndex )( 0, 0 ) = tCounter;
                mDvAssemblyMap( tDvIndex )( 0, 1 ) = tCounter + tNumCoeff - 1;

                // update the dv coefficient counter
                tCounter += tNumCoeff;
            }
        }

        // loop over the slave dv types
        for( uint Ik = 0; Ik < mSlaveDvTypes.size(); Ik++ )
        {
            //get the set index for the slave dv type
            sint tDvIndex = this->get_dv_index_for_type( mSlaveDvTypes( Ik )( 0 ),
                                                         mtk::Master_Slave::SLAVE );

            // if the dv type was set (its set index is different from -1)
            if( tDvIndex != -1 )
            {
                // get the number of coefficients for the slave dv type
                uint tNumCoeff = mSlaveFIManager->get_field_interpolators_for_type( mSlaveDvTypes( Ik )( 0 ) )
                                                ->get_number_of_space_time_coefficients();

                // fill the residual assembly map with starting and ending indices for the slave dof type
                mDvAssemblyMap( tDvIndex )( 0, 0 ) = tCounter;
                mDvAssemblyMap( tDvIndex )( 0, 1 ) = tCounter + tNumCoeff - 1;

                // update the dof coefficient counter
                tCounter += tNumCoeff;
            }
        }
    }

//------------------------------------------------------------------------------
    void Set::create_requested_IWG_list()
    {
        // get list of requested dof types from solver
        moris::Cell < enum MSI::Dof_Type > tRequestedDofTypes
        =  this->get_model_solver_interface()
               ->get_solver_interface()
               ->get_requested_dof_types();

        // clear requested IWG list
        mRequestedIWGs.clear();

        // reserve max size for requested IWG list
        mRequestedIWGs.reserve( mIWGs.size() );

        // loop over the requested dof types
        for( MSI::Dof_Type tDofType : tRequestedDofTypes )
        {
            // loop over the IWG in set IWG list
            for( uint iIWG = 0; iIWG < mIWGs.size(); iIWG++ )
            {
                // if the IWG residual dof type is requested
                if( mIWGs( iIWG )->get_residual_dof_type()( 0 ) == tDofType )
                {
                    // add the IWg to the requested IWG list
                    mRequestedIWGs.push_back( mIWGs( iIWG ) );
                }
            }
        }
        // reduce the size of requested IWG list to fit
        mRequestedIWGs.shrink_to_fit();
    }

//------------------------------------------------------------------------------
    void Set::create_requested_IQI_list()
    {
        // get list of requested IQI types from OPT
        mRequestedIQIs = mIQIs;

//        // FIXME use when MSI/GEN interface supports this
          // FIXME set requested types fomr GEN or SOL. do not use get function
//        // get list of requested IQI types from OPT through the design variable interface
//        moris::Cell < enum FEM::IQI_Type > mRequestedIQITypes;
//
//        // clear requested IQI list
//        mRequestedIQIs.clear();
//
//        // set size for requested IWG list
//        mRequestedIWGs.resize( tRequestedIQITypes.size() );
//
//        // loop over the requested IQI types
//        for( FEM::IQI_Type tIQIType : tRequestedIQITypes )
//        {
//            // loop over the IQI in set IQI list
//            for( uint iIQI = 0; iIQI < mIQIs.size(); iIQI++ )
//            {
//                // if the IQI type is requested
//                if( mIQIs( iIQI )->get_fem_IQI_type() == tIQIType )
//                {
//                    // add the IQI to the requested IQI list
//                    mRequestedIQIs.push_back( mIQIs( iIQI ) );
//                    break;
//                }
//            }
//        }
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
        // if residual not initialized before
        if ( !mJacobianExist )
        {
            // get the dof types requested by the solver
            moris::Cell< enum MSI::Dof_Type > tRequestedDofTypes
            = this->get_requested_dof_types();

            // init dof coefficient counter
            uint tNumCols = 0;

            // loop over the requested dof types
            for ( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
            {
                // get the set index for the master dof type
                sint tDofIndex = this->get_dof_index_for_type(tRequestedDofTypes( Ik ), mtk::Master_Slave::MASTER);

                // if this master dof is active
                if ( tDofIndex != -1 )
                {
                    // update number of dof coefficients
                    tNumCols += mMasterFIManager->get_field_interpolators_for_type(tRequestedDofTypes(Ik))
                                                ->get_number_of_space_time_coefficients();
                }

                // get the set index for the slave dof type
                tDofIndex = this->get_dof_index_for_type(tRequestedDofTypes(Ik), mtk::Master_Slave::SLAVE);

                // if this slave dof is active
                if ( tDofIndex != -1 )
                {
                    // update number of dof coefficients
                    tNumCols += mSlaveFIManager->get_field_interpolators_for_type(tRequestedDofTypes(Ik))
                                               ->get_number_of_space_time_coefficients();
                }
            }

            // if for residual evaluation
            if( !mIsResidual )
            {
                // set size for the jacobian matrix
                mJacobian.set_size( tNumCols, tNumCols, 0.0 );
            }
            // if for jacobian evaluation
            else
            {
                // get the secondary dof types from the solver
                moris::Cell< moris::Cell< enum MSI::Dof_Type > > tSecDofTypes
                = this->get_secundary_dof_types();

                // init dof coefficient counter for rows
                uint tNumRows = 0;

                // loop over the groups of secondary dof types
                for ( auto tSecDofTypesI : tSecDofTypes )
                {
                    // loop over the secondary dof types in group
                    for (uint Ik = 0; Ik < tSecDofTypesI.size(); Ik++)
                    {
                        // get the set index for the master dof type
                        sint tDofIndex = this->get_dof_index_for_type( tSecDofTypesI( Ik ),
                                                                       mtk::Master_Slave::MASTER );

                        // if this master dof is active
                        if ( tDofIndex != -1 )
                        {
                            // update number of dof coefficients
                            tNumRows += mMasterFIManager->get_field_interpolators_for_type( tSecDofTypesI( Ik ) )
                                                        ->get_number_of_space_time_coefficients();
                        }

                        // get the set index for the slave dof type
                        tDofIndex = this->get_dof_index_for_type( tSecDofTypesI(Ik), mtk::Master_Slave::SLAVE );

                        // if this slave dof is active
                        if ( tDofIndex != -1 )
                        {
                            // update number of dof coefficients
                            tNumRows += mSlaveFIManager->get_field_interpolators_for_type( tSecDofTypesI( Ik ) )
                                                       ->get_number_of_space_time_coefficients();
                        }
                    }
                }

                // set size for the jacobian matrix
                mJacobian.set_size( tNumCols, tNumRows, 0.0 );
            }
            // set the jacobian initialization flag to true
            mJacobianExist = true;
        }
        else
        {
//            MORIS_ASSERT( mJacobian.numel() > 0, "Set::initialize_mJacobian() - Jacobian not properly initialized.");
            // fill the jacobian matrix with zeros
            mJacobian.fill( 0.0 );
        }
    }

//------------------------------------------------------------------------------
    void Set::initialize_mResidual()
    {
        // if residual not initialized before
        if ( !mResidualExist )
        {
           // get the dof types requested by the solver
            moris::Cell < enum MSI::Dof_Type >tRequestedDofTypes
            = this->get_requested_dof_types();

            // init dof coefficient counter
            uint tNumCoeff = 0;

            // loop over the requested dof types
            for( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
            {
                // get the set index for the master dof type
                sint tDofIndex = this->get_dof_index_for_type( tRequestedDofTypes( Ik ),
                                                               mtk::Master_Slave::MASTER );

                // if this master dof is active
                if( tDofIndex != -1 )
                {
                    // update number of dof coefficients
                    tNumCoeff += mMasterFIManager->get_field_interpolators_for_type( tRequestedDofTypes( Ik ) )
                                                 ->get_number_of_space_time_coefficients();
                }

                // get the set index for the slave dof type
                tDofIndex = this->get_dof_index_for_type( tRequestedDofTypes( Ik ),
                                                          mtk::Master_Slave::SLAVE  );

                // if this slave dof is active
                if( tDofIndex != -1 )
                {
                    // update number of dof coefficients
                    tNumCoeff += mSlaveFIManager->get_field_interpolators_for_type( tRequestedDofTypes( Ik ) )
                                                ->get_number_of_space_time_coefficients();
                }
            }

            // if forward analysis (residual vector only, no dQIdu or dRdp)
            if( mIsForward )
            {
                // set size for the list of residual vector
                mResidual.resize( 1 );

                // set size for the residual vector
                mResidual( 0 ).set_size( tNumCoeff, 1, 0.0 );
            }
            // if sensitivity analysis ( dRdp + dQIdu)
            else
            {
//                //get the number of requested dv types
//                uint tNumRequestedDv  = this->get_requested_dv_types().size();

                // get the number of requested QI
                uint tNumRequestedIQI = this->get_requested_IQIs().size();

                 // set size for the list of dQIdu vectors
                mResidual.resize( tNumRequestedIQI );

                // loop over the dQIdu vectors
                for( auto & tRes : mResidual )
                {
                    // set size for the dQIdu vector
                    tRes.set_size( tNumCoeff, 1, 0.0 );
                }
            }

            // set the residual initialization flag to true
            mResidualExist = true;
        }
        // if residual initialized before
        else
        {
            // loop over the residual vectors
            for( auto & tRes : mResidual )
            {
                // fill the residual vector with zeros
                tRes.fill( 0.0 );
            }
        }
    }

//------------------------------------------------------------------------------
    void Set::initialize_mQI()
    {
        // if list of QI values not initialized before
        if ( !mQIExist )
        {
           // get the dof types requested by the solver
//            moris::Cell < moris::Cell< moris_index > > tQIAssemblyMap
//            = this->get_QI_assembly_map();
//
//            // init QI counter
//            uint tNumQI= 0;
//
//            // loop over the requested dof types
//            for( uint Ik = 0; Ik < tQIAssemblyMap.size(); Ik++ )
//            {
//                for( uint Jk = 0; Jk < tQIAssemblyMap( Ik ).size(); Jk++ )
//                {
//                    if( tQIAssemblyMap( Ik )( Jk ) != -1 )
//                    {
//                        tNumQI++;
//                    }
//                }
//            }

                   uint tNumQI= 1;

            // set size for the list of QI values
            mQI.resize( tNumQI );

            for( auto & tQI : mQI )
            {
                // set size for the QI value
                // FIXMe assumed scalar
                tQI.set_size( 1, 1, 0.0 );
            }
            // set the QI initialization flag to true
            mQIExist = true;
        }
        // if list of QI values initialized before
        else
        {
            // loop over the QI values
            for( auto & tQI : mQI )
            {
                // fill the QI value vector with zero
                tQI.fill( 0.0 );
            }
        }
    }

//------------------------------------------------------------------------------
        void Set::initialize_mdQIdp()
        {
            MORIS_ERROR( false, "Set::initialize_mdQIdp - not implemented yet." );
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
                       MORIS_ERROR( false, " Set::get_auto_interpolation_order - not defined for LINE and number of vertices. ");
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
                        MORIS_ERROR( false, " Set::get_auto_interpolation_order - not defined for QUAD and number of vertices. ");
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
                        MORIS_ERROR( false, " Set::get_auto_interpolation_order - not defined for HEX and number of vertices. ");
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
                        MORIS_ERROR( false, " Set::get_auto_interpolation_order - not defined for TET and number of vertices. ");
                        return mtk::Interpolation_Order::UNDEFINED;
                        break;
                }

            default :
                MORIS_ERROR( false, " Set::get_auto_interpolation_order - not defined for this geometry type. ");
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
    fem::Integration_Order Set::get_auto_integration_order( const mtk::Geometry_Type       aGeometryType,
                                                            const mtk::Interpolation_Order aInterpolationOrder )
    {
        switch( aGeometryType )
        {
            case( mtk::Geometry_Type::LINE ) :
            {
                switch( aInterpolationOrder )
                {
                    case( mtk::Interpolation_Order::LINEAR ):
                        return fem::Integration_Order::BAR_1;

                    case( mtk::Interpolation_Order::QUADRATIC ):
                        return fem::Integration_Order::BAR_2;

                    case( mtk::Interpolation_Order::CUBIC ):
                        return fem::Integration_Order::BAR_3;

                    default:
                        MORIS_ERROR( false, "Set::get_auto_integration_order - Unknown or unsupported interpolation order.");
                        return fem::Integration_Order::UNDEFINED;
                }
                break;
            }

            case( mtk::Geometry_Type::QUAD ) :
            {
                switch( aInterpolationOrder )
                {
                    case( mtk::Interpolation_Order::LINEAR ):
                         return fem::Integration_Order::QUAD_2x2;

                    case( mtk::Interpolation_Order::SERENDIPITY ):
                         return fem::Integration_Order::QUAD_3x3;

                    case( mtk::Interpolation_Order::QUADRATIC ):
                        return fem::Integration_Order::QUAD_3x3;

                    case( mtk::Interpolation_Order::CUBIC ):
                        return fem::Integration_Order::QUAD_4x4;

                    default:
                        MORIS_ERROR( false, "Set::get_auto_integration_order - Unknown or unsupported interpolation order.");
                        return fem::Integration_Order::UNDEFINED;
                }
                break;
            }

            case( mtk::Geometry_Type::HEX ) :
            {
                switch( aInterpolationOrder )
                {
                    case( mtk::Interpolation_Order::LINEAR ):
                         return fem::Integration_Order::HEX_2x2x2;

                    case( mtk::Interpolation_Order::SERENDIPITY ):
                        return fem::Integration_Order::HEX_3x3x3;

                    case( mtk::Interpolation_Order::QUADRATIC ):
                        return fem::Integration_Order::HEX_3x3x3;

                    case( mtk::Interpolation_Order::CUBIC ):
                        return fem::Integration_Order::HEX_4x4x4;

                    default:
                        MORIS_ERROR( false, "Set::get_auto_integration_order - Unknown or unsupported interpolation order.");
                        return fem::Integration_Order::UNDEFINED;
                }
                break;
            }

            case( mtk::Geometry_Type::TRI ) :
            {
                switch( aInterpolationOrder )
                {
                    case( mtk::Interpolation_Order::LINEAR ):
                         return fem::Integration_Order::TRI_3;

                    case( mtk::Interpolation_Order::QUADRATIC ):
                        return fem::Integration_Order::TRI_6;

                    case( mtk::Interpolation_Order::CUBIC ):
                        return fem::Integration_Order::TRI_7;

                    default:
                        MORIS_ERROR( false, "Set::get_auto_integration_order - Unknown or unsupported interpolation order.");
                        return fem::Integration_Order::UNDEFINED;
                }
                break;
            }

            case( mtk::Geometry_Type::TET ) :
            {
                switch( aInterpolationOrder )
                {
                    case( mtk::Interpolation_Order::LINEAR ):
                         return fem::Integration_Order::TET_4;

                    case( mtk::Interpolation_Order::QUADRATIC ):
                        return fem::Integration_Order::TET_11;

                    case( mtk::Interpolation_Order::CUBIC ):
                        return fem::Integration_Order::TET_15;

                    default:
                        MORIS_ERROR( false, "Set::get_auto_integration_order - Unknown or unsupported interpolation order.");
                        return fem::Integration_Order::UNDEFINED;
                }
                break;
            }

            default :
                MORIS_ERROR( false, " Set::get_auto_integration_order - Unknown or unsupported geometry type. ");
                return Integration_Order::UNDEFINED;
                break;
        }
    }

////------------------------------------------------------------------------------
//    void Set::create_requested_dv_assembly_map()
//    {
//        moris::Cell< enum GEN_DV > tRequestedDvTypes = this->get_requested_dv_types();
//
//        uint tMaxDvIndex = 0;
//
//        for( auto tDvType : tRequestedDvTypes )
//        {
//            tMaxDvIndex = std::max( tMaxDvIndex, static_cast< uint >( tDvType ) );
//        }
//
//        mQITypeMap.set_size( tMaxDvIndex+1, 1, -1 );
//
//        uint tCounter = 0;
//
//        for( auto tDvType : tRequestedDvTypes )
//        {
//            mQITypeMap( static_cast< uint >( tDvType ) ) = tCounter++;
//        }
//    }

//------------------------------------------------------------------------------

    void Set::set_visualization_set( const uint              aMeshIndex,
                                           moris::mtk::Set * aVisMeshSet,
                                     const bool              aOnlyPrimayCells)
    {
         uint tNumClustersOnSets = aVisMeshSet->get_num_clusters_on_set();

         // set vis clusters to clusters
         for( uint Ik = 0; Ik < tNumClustersOnSets; Ik++ )
         {
             // create a fem cluster
             std::shared_ptr< fem::Cluster > tCluster = std::make_shared< fem::Cluster >( mElementType,
                                                                                          aVisMeshSet->get_clusters_by_index( Ik ),
                                                                                          this,
                                                                                          mEquationObjList( Ik ));

             reinterpret_cast< fem::Interpolation_Element* >( mEquationObjList( Ik ) )->set_cluster( tCluster,
                                                                                                     aMeshIndex );
         }

         // build set element map
         uint tNumCells = aVisMeshSet->get_num_cells_on_set( aOnlyPrimayCells );

         moris::Matrix< DDSMat > tCellIndex = aVisMeshSet->get_cell_inds_on_block( aOnlyPrimayCells );

         sint tSize = std::max( ( sint )mCellAssemblyMap.size(), ( sint )aMeshIndex + 1 );

         mCellAssemblyMap.resize( tSize );
         mMtkIgCellOnSet .resize( tSize );

         mMtkIgCellOnSet( aMeshIndex )= tNumCells;

         if(tNumCells>0)
         {
             sint tMaxIndex = tCellIndex.max();
//             sint tMinIndex = tCellIndex.min();

             mCellAssemblyMap( aMeshIndex ).set_size( tMaxIndex + 1, 1, -1 );

             for( uint Ik = 0; Ik < tNumCells; Ik++ )
             {
                 mCellAssemblyMap( aMeshIndex )( tCellIndex( Ik ) ) = Ik;
             }

         }

    }

//------------------------------------------------------------------------------
    void Set::compute_quantity_of_interest( const uint              aMeshIndex,
                                            Matrix< DDRMat >      * aElementFieldValues,
                                            Matrix< DDRMat >      * aNodalFieldValues,
                                            moris::real           * aGlobalScalar,
                                            enum vis::Output_Type   aOutputType,
                                            enum vis::Field_Type    aFieldType )
    {
        mSetElementalValues = aElementFieldValues;
        mSetNodalValues     = aNodalFieldValues;
        mSetGlobalValues    = aGlobalScalar;

        mSetNodalCounter.set_size( (*mSetNodalValues).numel(), 1, 0 );

        mSetElementalValues->set_size( mMtkIgCellOnSet( aMeshIndex ), 1, 0.0 );

        for( uint Ik = 0; Ik < mEquationObjList.size(); Ik++ )
        {
            mEquationObjList( Ik )->compute_quantity_of_interest( aMeshIndex,
                                                                  aOutputType,
                                                                  aFieldType );
        }

        //FIXME I do not like this at all. someone change it
        for( uint Ik = 0; Ik < mSetNodalValues->numel(); Ik++ )
        {
            if( mSetNodalCounter(Ik) != 0)
            {
                (*mSetNodalValues)(Ik) = (*mSetNodalValues)(Ik)/mSetNodalCounter(Ik);
            }
        }
//        if( aFieldType==vis::Field_Type::NODAL || aFieldType==vis::Field_Type::NODAL_IP )
//        {
//            print( *mSetNodalValues, "mSetNodalValues" );
//        }
    }

//------------------------------------------------------------------------------
    void Set::determine_set_type()
    {
        enum moris::SetType tMtkSetType = mMeshSet->get_set_type();

        switch( tMtkSetType )
        {
            case( moris::SetType::BULK ) :
                mElementType = fem::Element_Type::BULK;
                break;

            case( moris::SetType::SIDESET ) :
                 mElementType = fem::Element_Type::SIDESET;
                 break;

            case( moris::SetType::DOUBLE_SIDED_SIDESET ) :
                mElementType = fem::Element_Type::DOUBLE_SIDESET;
                break;

            default :
                MORIS_ERROR( false, "Set::determine_set_type() - not defined for this set type. ");
                break;
        }
    }

//------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
