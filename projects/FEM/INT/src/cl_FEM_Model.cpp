
// added by christian: link to Google Perftools
#ifdef WITHGPERFTOOLS
#include <gperftools/profiler.h>
#endif

#include "cl_Stopwatch.hpp" //CHR/src

#include "cl_Map.hpp"
#include "fn_unique.hpp"
#include "fn_sum.hpp" // for check
#include "fn_iscol.hpp"
#include "fn_trans.hpp"
#include "op_equal_equal.hpp"

#include "MTK_Tools.hpp"
#include "cl_MTK_Enums.hpp"              //MTK/src
#include "cl_MTK_Mesh_Manager.hpp"       //MTK/src

#include "cl_FEM_Node_Base.hpp"          //FEM/INT/src
#include "cl_FEM_Node.hpp"               //FEM/INT/src
#include "cl_FEM_Enums.hpp"              //FEM/INT/src

#include "cl_FEM_Model.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Set_User_Info.hpp"
#include "cl_FEM_CM_Factory.hpp"
#include "cl_FEM_SP_Factory.hpp"
#include "cl_FEM_IWG_Factory.hpp"
#include "cl_FEM_IQI_Factory.hpp"

#include "cl_MSI_Equation_Object.hpp"
#include "cl_MSI_Dof_Type_Enums.hpp"

#include "fn_Exec_load_user_library.hpp"


namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
        FEM_Model::FEM_Model(       mtk::Mesh_Manager *                 aMeshManager,
                              const moris_index                       & aMeshPairIndex,
                                    moris::Cell< fem::Set_User_Info > & aSetInfo ) : mMeshManager( aMeshManager ),
                                                                                     mMeshPairIndex( aMeshPairIndex )
        {
            // get the number of sets
            uint tNumFemSets = aSetInfo.size();

            // start timer
            tic tTimer1;

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 0: initialize
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // Get pointers to interpolation and integration mesh
            mtk::Interpolation_Mesh* tInterpolationMesh = nullptr;
            mtk::Integration_Mesh*   tIntegrationMesh   = nullptr;
            mMeshManager->get_mesh_pair( mMeshPairIndex, tInterpolationMesh, tIntegrationMesh );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 1: create nodes
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // ask mesh about number of IP nodes on proc
            luint tNumOfIPNodes = tInterpolationMesh->get_num_nodes();

            // create IP node objects
            mIPNodes.resize( tNumOfIPNodes, nullptr );

            for( uint iNode = 0; iNode < tNumOfIPNodes; iNode++ )
            {
                mIPNodes( iNode ) = new fem::Node( &tInterpolationMesh->get_mtk_vertex( iNode ) );
            }

            if( par_rank() == 0)
            {
                // stop timer
                real tElapsedTime = tTimer1.toc<moris::chronos::milliseconds>().wall;

                // print output
                MORIS_LOG_INFO( "Model: created %u FEM IP nodes in %5.3f seconds.\n\n",
                                ( unsigned int ) tNumOfIPNodes,
                                ( double ) tElapsedTime / 1000 );
            }

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 2: create elements
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // start timer
            tic tTimer2;

            // create equation objects
            mFemSets.resize( tNumFemSets, nullptr );     // FIXME try to create them as equation sets

            // get the number of element to create
            mFemClusters.reserve( 100000 ); //FIXME

            //------------------------------------------------------------------------------
            // init the fem set counter
            moris::uint tFemSetCounter = 0;

            // loop over the used fem set
            for( luint iSet = 0; iSet < tNumFemSets; iSet++ )
            {
                mMeshSetToFemSetMap[ aSetInfo( iSet ).get_mesh_index() ] = tFemSetCounter;

                // create a list of clusters
                moris::mtk::Set * tMeshSet = tIntegrationMesh->get_set_by_index( aSetInfo( iSet ).get_mesh_index() );

//                std::cout<<tMeshSet->get_num_clusters_on_set()<<std::endl;
                if ( tMeshSet->get_num_clusters_on_set() !=0 )
                {
                    // create new fem set
                    mFemSets( tFemSetCounter ) = new fem::Set( this,
                                                               tMeshSet,
                                                               aSetInfo( iSet ),
                                                               mIPNodes );
                }
                else
                {
                    // FIXME why do we build empty set?
                    mFemSets( tFemSetCounter ) = new fem::Set();
                }

                // collect equation objects associated with the block-set
                mFemClusters.append( mFemSets( tFemSetCounter )->get_equation_object_list() );

                // update fem set counter
                tFemSetCounter++;
            }
            mFemClusters.shrink_to_fit();

            if( par_rank() == 0)
            {
                // stop timer
                real tElapsedTime = tTimer2.toc<moris::chronos::milliseconds>().wall;

                // print output

                MORIS_LOG_INFO( "Model: created %u FEM elements in %5.3f seconds.\n\n",
                        ( unsigned int ) mFemClusters.size(),
                        ( double ) tElapsedTime / 1000 );
            }
        }

//------------------------------------------------------------------------------
        void FEM_Model::initialize( moris::Cell< moris::Cell< ParameterList > > & aParameterList )
        {
            // save parameter list as a member data
            mParameterList = aParameterList;

            // get msi string to dof type map
            moris::map< std::string, MSI::Dof_Type > tMSIDofTypeMap
            = moris::MSI::get_msi_dof_type_map();

            // get msi string to dv type map
            moris::map< std::string, MSI::Dv_Type > tMSIDvTypeMap
            = moris::MSI::get_msi_dv_type_map();

            // create properties
            moris::Cell< std::shared_ptr< fem::Property > > tProperties;
            moris::map< std::string, uint > tPropertyMap;
            this->create_properties( tProperties, tPropertyMap, mParameterList( 0 ), tMSIDofTypeMap, tMSIDvTypeMap );

            // create constitutive models
            moris::Cell< std::shared_ptr< fem::Constitutive_Model > > tCMs;
            moris::map< std::string, uint > tCMMap;
            this->create_constitutive_models( tCMs, tCMMap, mParameterList( 1 ), tProperties, tPropertyMap, tMSIDofTypeMap, tMSIDvTypeMap );

            // create stabilization parameters
            moris::Cell< std::shared_ptr< fem::Stabilization_Parameter > > tSPs;
            moris::map< std::string, uint > tSPMap;
            this->create_stabilization_parameters( tSPs, tSPMap, mParameterList( 2 ), tProperties, tPropertyMap, tCMs, tCMMap, tMSIDofTypeMap, tMSIDvTypeMap );

            // create IWGs
            moris::Cell< std::shared_ptr< fem::IWG > > tIWGs;
            this->create_IWGs( tIWGs, mParameterList( 3 ), tProperties, tPropertyMap, tCMs, tCMMap, tSPs, tSPMap, tMSIDofTypeMap, tMSIDvTypeMap );

            // create IQIs
            moris::Cell< std::shared_ptr< fem::IQI > > tIQIs;
            this->create_IQIs( tIQIs, mParameterList( 4 ), tProperties, tPropertyMap, tCMs, tCMMap, tSPs, tSPMap, tMSIDofTypeMap, tMSIDvTypeMap );
        }

//------------------------------------------------------------------------------

        FEM_Model::~FEM_Model()
        {
            // delete fem nodes
            for( auto tIPNodes : mIPNodes )
            {
                delete tIPNodes;
            }
            mIPNodes.clear();

            // delete the fem sets
            for( auto tFemSet : mFemSets )
            {
                delete tFemSet;
            }
            mFemSets.clear();

            // delete the fem cluster
            mFemClusters.clear();
        }

//------------------------------------------------------------------------------
        void FEM_Model::create_properties
        ( moris::Cell< std::shared_ptr< fem::Property > > & aProperties,
          moris::map< std::string, uint >                 & aPropertyMap,
          moris::Cell< ParameterList >                    & aParameterList,
          moris::map< std::string, MSI::Dof_Type >        & aMSIDofTypeMap,
          moris::map< std::string, MSI::Dv_Type >         & aMSIDvTypeMap )
        {

            Library_IO tLibrary( mFilePath );

            // get the number of properties
            uint tNumProps = aParameterList.size();

            // create a list of property pointers
            aProperties.resize( tNumProps );

            // loop over the parameter lists
            for ( uint iProp = 0; iProp < tNumProps; iProp++ )
            {
                // create a property pointer
                aProperties( iProp ) = std::make_shared< fem::Property >();

                // set a name for the property
                aProperties( iProp )->set_name( aParameterList( iProp ).get< std::string >( "property_name" ) );

                // fill property map
                aPropertyMap[ aParameterList( iProp ).get< std::string >( "property_name" ) ] = iProp;

                // set dof dependencies
                moris::Cell< moris::Cell< moris::MSI::Dof_Type > > tDofTypes;
                string_to_cell_of_cell( aParameterList( iProp ).get< std::string >( "dof_dependencies" ), tDofTypes, aMSIDofTypeMap );
                aProperties( iProp )->set_dof_type_list( tDofTypes );

                // set dof dependencies
                moris::Cell< moris::Cell< moris::MSI::Dv_Type > > tDvTypes;
                string_to_cell_of_cell( aParameterList( iProp ).get< std::string >( "dv_dependencies" ), tDvTypes, aMSIDvTypeMap );
                aProperties( iProp )->set_dv_type_list( tDvTypes );

                // set function parameters
                moris::Cell< moris::Matrix< DDRMat > > tFuncParameters;
                string_to_cell_mat_2( aParameterList( iProp ).get< std::string >( "function_parameters" ), tFuncParameters );
                aProperties( iProp )->set_parameters( tFuncParameters );

                // set value function for property
                std::string tValFuncName = aParameterList( iProp ).get< std::string >( "value_function" );
                MORIS_FEM_FREE_FUNCTION tValFunction = nullptr;
                if ( tValFuncName.size() > 1 )
                {
                    tValFunction = tLibrary.load_fem_free_functions( tValFuncName );
                }
                aProperties( iProp )->set_val_function( tValFunction );

                // set dof derivative function for property
                moris::Cell< std::string > tDofDerFuncNames;
                string_to_cell( aParameterList( iProp ).get< std::string >( "dof_derivative_functions" ), tDofDerFuncNames );
                uint tNumDofDerFuncs = tDofDerFuncNames.size();
                moris::Cell< fem::PropertyFunc > tDofDerFunctions( tNumDofDerFuncs, nullptr );
                for( uint iFunc = 0; iFunc < tNumDofDerFuncs; iFunc++ )
                {
                    if( tDofDerFuncNames( iFunc ).size() > 1 )
                    {
                        MORIS_FEM_FREE_FUNCTION tValFunction = tLibrary.load_fem_free_functions( tDofDerFuncNames( iFunc ) );
                        tDofDerFunctions( iFunc ) = tValFunction;
                    }
                }
                aProperties( iProp )->set_dof_derivative_functions( tDofDerFunctions );

                // set dv derivative function for property
                moris::Cell< std::string > tDvDerFuncNames;
                string_to_cell( aParameterList( iProp ).get< std::string >( "dv_derivative_functions" ), tDvDerFuncNames );
                uint tNumDvDerFuncs = tDvDerFuncNames.size();
                moris::Cell< fem::PropertyFunc > tDvDerFunctions( tNumDvDerFuncs, nullptr );
                for( uint iFunc = 0; iFunc < tNumDvDerFuncs; iFunc++ )
                {
                    if( tDvDerFuncNames( iFunc ).size() > 1 )
                    {
                        MORIS_FEM_FREE_FUNCTION tValFunction = tLibrary.load_fem_free_functions( tDvDerFuncNames( iFunc ) );
                        tDvDerFunctions( iFunc ) = tValFunction;
                    }
                }
                aProperties( iProp )->set_dv_derivative_functions( tDvDerFunctions );
            }
        }

//------------------------------------------------------------------------------
        void FEM_Model::create_constitutive_models
        ( moris::Cell< std::shared_ptr< fem::Constitutive_Model > > & aCMs,
          moris::map< std::string, uint >                           & aCMMap,
          moris::Cell< ParameterList >                              & aParameterList,
          moris::Cell< std::shared_ptr< fem::Property > >           & aProperties,
          moris::map< std::string, uint >                           & aPropertyMap,
          moris::map< std::string, MSI::Dof_Type >                  & aMSIDofTypeMap,
          moris::map< std::string, MSI::Dv_Type >                   & aMSIDvTypeMap )
        {
            // create a constitutive model factory
            CM_Factory tCMFactory;

            // get number of constitutive models
            uint tNumCMs = aParameterList.size();

            // create a list of CMs
            aCMs.resize( tNumCMs );

            // loop over the parameter lists for CM
            for( uint iCM = 0; iCM < tNumCMs; iCM++ )
            {
                // get the constitutive type from parameter list
                fem::Constitutive_Type tCMType
                = static_cast< fem::Constitutive_Type >( aParameterList( iCM ).get< uint >( "constitutive_type" ) );

                // create a constitutive model pointer
                aCMs( iCM ) = tCMFactory.create_CM( tCMType );

                // set CM name
                aCMs( iCM )->set_name( aParameterList( iCM ).get< std::string >( "constitutive_name" ) );

                // fill property map
                aCMMap[ aParameterList( iCM ).get< std::string >( "constitutive_name" ) ] = iCM;

                // set CM dof dependencies
                moris::Cell< moris::Cell< moris::MSI::Dof_Type > > tDofTypes;
                string_to_cell_of_cell( std::get< 0 >( aParameterList( iCM ).get< std::pair< std::string, std::string > >( "dof_dependencies" ) ), tDofTypes, aMSIDofTypeMap );
                moris::Cell< std::string > tDofTypeNames;
                string_to_cell( std::get< 1 >( aParameterList( iCM ).get< std::pair< std::string, std::string > >( "dof_dependencies" ) ), tDofTypeNames );
                aCMs( iCM )->set_dof_type_list( tDofTypes, tDofTypeNames );

                // set CM dv dependencies
                moris::Cell< moris::Cell< moris::MSI::Dv_Type > > tDvTypes;
                string_to_cell_of_cell( std::get< 0 >( aParameterList( iCM ).get< std::pair< std::string, std::string > >( "dv_dependencies" ) ), tDvTypes, aMSIDvTypeMap );
                moris::Cell< std::string > tDvTypeNames;
                string_to_cell( std::get< 1 >( aParameterList( iCM ).get< std::pair< std::string, std::string > >( "dv_dependencies" ) ), tDvTypeNames );
                aCMs( iCM )->set_dv_type_list( tDvTypes, tDvTypeNames );

                // set CM properties
                moris::Cell< moris::Cell< std::string > > tPropertyNamesPair;
                string_to_cell_of_cell( aParameterList( iCM ).get< std::string >( "properties" ), tPropertyNamesPair );
                for( uint iProp = 0; iProp < tPropertyNamesPair.size(); iProp++ )
                {
                    // get property index
                    uint tPropertyIndex = aPropertyMap[ tPropertyNamesPair( iProp)( 0 ) ];

                    // set property for CM
                    aCMs( iCM )->set_property( aProperties( tPropertyIndex ), tPropertyNamesPair( iProp )( 1 ) );
                }

                // debug
                aCMs( iCM )->print_names();
            }
        }

//------------------------------------------------------------------------------
        void FEM_Model::create_stabilization_parameters
        ( moris::Cell< std::shared_ptr< fem::Stabilization_Parameter > > & aSPs,
          moris::map< std::string, uint >                                & aSPMap,
          moris::Cell< ParameterList >                                   & aParameterList,
          moris::Cell< std::shared_ptr< fem::Property > >                & aProperties,
          moris::map< std::string, uint >                                & aPropertyMap,
          moris::Cell< std::shared_ptr< fem::Constitutive_Model > >      & aCMs,
          moris::map< std::string, uint >                                & aCMMap,
          moris::map< std::string, MSI::Dof_Type >                       & aMSIDofTypeMap,
          moris::map< std::string, MSI::Dv_Type >                        & aMSIDvTypeMap )
        {
            // create a stabilization parameter factory
            SP_Factory tSPFactory;

            // get the number of stabilization parameters
            uint tNumSPs = aParameterList.size();

            // set size for the list of stabilization parameter pointer
            aSPs.resize( tNumSPs );

            // loop over the parameter list
            for ( uint iSP = 0; iSP < tNumSPs; iSP++ )
            {
                // get the stabilization type from parameter list
                fem::Stabilization_Type tSPType
                = static_cast< fem::Stabilization_Type >( aParameterList( iSP ).get< uint >( "stabilization_type" ) );

                // create a stabilization parameter pointer
                aSPs( iSP ) = tSPFactory.create_SP( tSPType );

                // set name
                aSPs( iSP )->set_name( aParameterList( iSP ).get< std::string >( "stabilization_name" ) );

                // fill stabilization map
                aSPMap[ aParameterList( iSP ).get< std::string >( "stabilization_name" ) ] = iSP;

                // set parameters
                moris::Cell< moris::Matrix< DDRMat > > tFuncParameters;
                string_to_cell_mat_2( aParameterList( iSP ).get< std::string >( "function_parameters" ), tFuncParameters );
                aSPs( iSP )->set_parameters( tFuncParameters );

                // set master dof dependencies
                moris::Cell< moris::Cell< moris::MSI::Dof_Type > > tDofTypes;
                string_to_cell_of_cell( std::get< 0 >( aParameterList( iSP ).get< std::pair< std::string, std::string > >( "master_dof_dependencies" ) ), tDofTypes, aMSIDofTypeMap );
                moris::Cell< std::string > tDofTypeNames;
                string_to_cell( std::get< 1 >( aParameterList( iSP ).get< std::pair< std::string, std::string > >( "master_dof_dependencies" ) ), tDofTypeNames );
                aSPs( iSP )->set_dof_type_list( tDofTypes, tDofTypeNames );

                // set slave dof dependencies
                string_to_cell_of_cell( std::get< 0 >( aParameterList( iSP ).get< std::pair< std::string, std::string > >( "slave_dof_dependencies" ) ), tDofTypes, aMSIDofTypeMap );
                string_to_cell( std::get< 1 >( aParameterList( iSP ).get< std::pair< std::string, std::string > >( "slave_dof_dependencies" ) ), tDofTypeNames );
                aSPs( iSP )->set_dof_type_list( tDofTypes, tDofTypeNames, mtk::Master_Slave::SLAVE );

                // set master dv dependencies
                moris::Cell< moris::Cell< moris::MSI::Dv_Type > > tDvTypes;
                string_to_cell_of_cell( std::get< 0 >( aParameterList( iSP ).get< std::pair< std::string, std::string > >( "master_dv_dependencies" ) ), tDvTypes, aMSIDvTypeMap );
                moris::Cell< std::string > tDvTypeNames;
                string_to_cell( std::get< 1 >( aParameterList( iSP ).get< std::pair< std::string, std::string > >( "master_dv_dependencies" ) ), tDvTypeNames );
                aSPs( iSP )->set_dv_type_list( tDvTypes, tDvTypeNames );

                // set slave dof dependencies
                string_to_cell_of_cell( std::get< 0 >( aParameterList( iSP ).get< std::pair< std::string, std::string > >( "slave_dv_dependencies" ) ), tDvTypes, aMSIDvTypeMap );
                string_to_cell( std::get< 1 >( aParameterList( iSP ).get< std::pair< std::string, std::string > >( "slave_dv_dependencies" ) ), tDvTypeNames );
                aSPs( iSP )->set_dv_type_list( tDvTypes, tDvTypeNames, mtk::Master_Slave::SLAVE );

                // set master properties
                moris::Cell< moris::Cell< std::string > > tMasterPropertyNamesPair;
                string_to_cell_of_cell( aParameterList( iSP ).get< std::string >( "master_properties" ), tMasterPropertyNamesPair );

                for( uint iProp = 0; iProp < tMasterPropertyNamesPair.size(); iProp++ )
                {
                    // get property index
                    uint tPropertyIndex = aPropertyMap[ tMasterPropertyNamesPair( iProp )( 0 ) ];

                    // set property for CM
                    aSPs( iSP )->set_property( aProperties( tPropertyIndex ), tMasterPropertyNamesPair( iProp )( 1 ) );
                }

                // set slave properties
                moris::Cell< moris::Cell< std::string > > tSlavePropertyNamesPair;
                string_to_cell_of_cell( aParameterList( iSP ).get< std::string >( "slave_properties" ), tSlavePropertyNamesPair );

                for( uint iProp = 0; iProp < tSlavePropertyNamesPair.size(); iProp++ )
                {
                    // get property index
                    uint tPropertyIndex = aPropertyMap[ tSlavePropertyNamesPair( iProp )( 0 ) ];

                    // set property for CM
                    aSPs( iSP )->set_property( aProperties( tPropertyIndex ), tSlavePropertyNamesPair( iProp )( 1 ), mtk::Master_Slave::SLAVE );
                }

                // set master constitutive models
                moris::Cell< moris::Cell< std::string > > tMasterCMNamesPair;
                string_to_cell_of_cell( aParameterList( iSP ).get< std::string >( "master_constitutive_models" ), tMasterCMNamesPair );

                for( uint iCM = 0; iCM < tMasterCMNamesPair.size(); iCM++ )
                {
                    // get CM index
                    uint tCMIndex = aCMMap[ tMasterCMNamesPair( iCM )( 0 ) ];

                    // set CM for SP
                    aSPs( iSP )->set_constitutive_model( aCMs( tCMIndex ), tMasterCMNamesPair( iCM )( 1 ) );
                }

                // set slave constitutive models
                 moris::Cell< moris::Cell< std::string > > tSlaveCMNamesPair;
                 string_to_cell_of_cell( aParameterList( iSP ).get< std::string >( "slave_constitutive_models" ), tSlaveCMNamesPair );

                 for( uint iCM = 0; iCM < tSlaveCMNamesPair.size(); iCM++ )
                 {
                     // get CM index
                     uint tCMIndex = aCMMap[ tSlaveCMNamesPair( iCM )( 0 ) ];

                     // set CM for SP
                     aSPs( iSP )->set_constitutive_model( aCMs( tCMIndex ), tSlaveCMNamesPair( iCM )( 1 ), mtk::Master_Slave::SLAVE );
                 }

                 // debug
                 aSPs( iSP )->print_names();
            }
        }

//------------------------------------------------------------------------------
        void FEM_Model::create_IWGs
        ( moris::Cell< std::shared_ptr< fem::IWG > >                     & aIWGs,
          moris::Cell< ParameterList >                                   & aParameterList,
          moris::Cell< std::shared_ptr< fem::Property > >                & aProperties,
          moris::map< std::string, uint >                                & aPropertyMap,
          moris::Cell< std::shared_ptr< fem::Constitutive_Model > >      & aCMs,
          moris::map< std::string, uint >                                & aCMMap,
          moris::Cell< std::shared_ptr< fem::Stabilization_Parameter > > & aSPs,
          moris::map< std::string, uint >                                & aSPMap,
          moris::map< std::string, MSI::Dof_Type >                       & aMSIDofTypeMap,
          moris::map< std::string, MSI::Dv_Type >                        & aMSIDvTypeMap )
        {
            // create an IWG factory
            IWG_Factory tIWGFactory;

            // get number of IWGs
            uint tNumIWGs = aParameterList.size();

            // create a list of IWG pointers
            aIWGs.resize( tNumIWGs );

            // loop over the parameter lists
            for( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
            {
                // get the IWG type from parameter list
                fem::IWG_Type tIWGType
                = static_cast< fem::IWG_Type >( aParameterList( iIWG ).get< uint >( "IWG_type" ) );

                // create an IWG pointer
                aIWGs( iIWG ) = tIWGFactory.create_IWG( tIWGType );

                // set name
                aIWGs( iIWG )->set_name( aParameterList( iIWG ).get< std::string >( "IWG_name" ) );

                // set residual dof type
                moris::Cell< moris::MSI::Dof_Type > tResDofTypes;
                string_to_cell( aParameterList( iIWG ).get< std::string >( "dof_residual" ),
                                tResDofTypes,
                                aMSIDofTypeMap );
                aIWGs( iIWG )->set_residual_dof_type( tResDofTypes );

                // set master dof dependencies
                moris::Cell< moris::Cell< moris::MSI::Dof_Type > > tDofTypes;
                string_to_cell_of_cell( aParameterList( iIWG ).get< std::string >( "master_dof_dependencies" ),
                                        tDofTypes,
                                        aMSIDofTypeMap );
                aIWGs( iIWG )->set_dof_type_list( tDofTypes );

                // set slave dof dependencies
                string_to_cell_of_cell( aParameterList( iIWG ).get< std::string >( "slave_dof_dependencies" ),
                                        tDofTypes,
                                        aMSIDofTypeMap );
                aIWGs( iIWG )->set_dof_type_list( tDofTypes );

                // set master dv dependencies
                moris::Cell< moris::Cell< moris::MSI::Dv_Type > > tDvTypes;
                string_to_cell_of_cell( aParameterList( iIWG ).get< std::string >( "master_dv_dependencies" ),
                                        tDvTypes,
                                        aMSIDvTypeMap );
                aIWGs( iIWG )->set_dv_type_list( tDvTypes );

                // set slave dv dependencies
                string_to_cell_of_cell( aParameterList( iIWG ).get< std::string >( "slave_dv_dependencies" ),
                                        tDvTypes,
                                        aMSIDvTypeMap );
                aIWGs( iIWG )->set_dv_type_list( tDvTypes );

                // set master properties
                moris::Cell< moris::Cell< std::string > > tMasterPropertyNamesPair;
                string_to_cell_of_cell( aParameterList( iIWG ).get< std::string >( "master_properties" ), tMasterPropertyNamesPair );

                for( uint iProp = 0; iProp < tMasterPropertyNamesPair.size(); iProp++ )
                {
                    // get property index
                    uint tPropertyIndex = aPropertyMap[ tMasterPropertyNamesPair( iProp )( 0 ) ];

                    // set property for IWG
                    aIWGs( iIWG )->set_property( aProperties( tPropertyIndex ), tMasterPropertyNamesPair( iProp )( 1 ) );
                }

                // set slave properties
                moris::Cell< moris::Cell< std::string > > tSlavePropertyNamesPair;
                string_to_cell_of_cell( aParameterList( iIWG ).get< std::string >( "slave_properties" ), tSlavePropertyNamesPair );

                for( uint iProp = 0; iProp < tSlavePropertyNamesPair.size(); iProp++ )
                {
                    // get property index
                    uint tPropertyIndex = aPropertyMap[ tMasterPropertyNamesPair( iProp )( 0 ) ];

                    // set property for IWG
                    aIWGs( iIWG )->set_property( aProperties( tPropertyIndex ), tSlavePropertyNamesPair( iProp )( 1 ), mtk::Master_Slave::SLAVE );
                }

                // set master constitutive models
                moris::Cell< moris::Cell< std::string > > tMasterCMNamesPair;
                string_to_cell_of_cell( aParameterList( iIWG ).get< std::string >( "master_constitutive_models" ), tMasterCMNamesPair );

                for( uint iCM = 0; iCM < tMasterCMNamesPair.size(); iCM++ )
                {
                    // get CM index
                    uint tCMIndex = aCMMap[ tMasterCMNamesPair( iCM )( 0 ) ];

                    // set CM for IWG
                    aIWGs( iIWG )->set_constitutive_model( aCMs( tCMIndex ), tMasterCMNamesPair( iCM )( 1 ) );
                }

                // set slave constitutive models
                moris::Cell< moris::Cell< std::string > > tSlaveCMNamesPair;
                string_to_cell_of_cell( aParameterList( iIWG ).get< std::string >( "slave_constitutive_models" ), tSlaveCMNamesPair );

                for( uint iCM = 0; iCM < tSlaveCMNamesPair.size(); iCM++ )
                {
                    // get CM index
                    uint tCMIndex = aCMMap[ tSlaveCMNamesPair( iCM )( 0 ) ];

                    // set CM for IWG
                    aIWGs( iIWG )->set_constitutive_model( aCMs( tCMIndex ), tSlaveCMNamesPair( iCM )( 1 ), mtk::Master_Slave::SLAVE );
                }

                // set stabilization parameters
                moris::Cell< moris::Cell< std::string > > tSPNamesPair;
                string_to_cell_of_cell( aParameterList( iIWG ).get< std::string >( "stabilization_parameters" ), tSPNamesPair );

                for( uint iSP = 0; iSP < tSPNamesPair.size(); iSP++ )
                {
                    // get SP index
                    uint tSPIndex = aSPMap[ tSPNamesPair( iSP )( 0 ) ];

                    // set SP for IWG
                    aIWGs( iIWG )->set_stabilization_parameter( aSPs( tSPIndex ), tSPNamesPair( iSP )( 1 ) );
                }

                // debug
                aIWGs( iIWG )->print_names();
            }
        }

//------------------------------------------------------------------------------
        void FEM_Model::create_IQIs
        ( moris::Cell< std::shared_ptr< fem::IQI > >                     & aIQIs,
          moris::Cell< ParameterList >                                   & aParameterList,
          moris::Cell< std::shared_ptr< fem::Property > >                & aProperties,
          moris::map< std::string, uint >                                & aPropertyMap,
          moris::Cell< std::shared_ptr< fem::Constitutive_Model > >      & aCMs,
          moris::map< std::string, uint >                                & aCMMap,
          moris::Cell< std::shared_ptr< fem::Stabilization_Parameter > > & aSPs,
          moris::map< std::string, uint >                                & aSPMap,
          moris::map< std::string, MSI::Dof_Type >                       & aMSIDofTypeMap,
          moris::map< std::string, MSI::Dv_Type >                        & aMSIDvTypeMap )
        {
            // create an IQI factory
            IQI_Factory tIQIFactory;

            // get number of IQIs
            uint tNumIQIs = aParameterList.size();

            // set size for list of IQI pointers
            aIQIs.resize( tNumIQIs );

            // loop over the parameter lists
            for( uint iIQI = 0; iIQI < tNumIQIs; iIQI++ )
            {
                // get the IQI type from parameter list
                fem::IQI_Type tIQIType
                = static_cast< fem::IQI_Type >( aParameterList( iIQI ).get< uint >( "IQI_type" ) );

               // create an IQI pointer
                aIQIs( iIQI ) = tIQIFactory.create_IQI( tIQIType );

                // set name
                aIQIs( iIQI )->set_name( aParameterList( iIQI ).get< std::string >( "IQI_name" ) );

                // set master dof dependencies
                moris::Cell< moris::Cell< moris::MSI::Dof_Type > > tDofTypes;
                string_to_cell_of_cell( aParameterList( iIQI ).get< std::string >( "master_dof_dependencies" ),
                                        tDofTypes,
                                        aMSIDofTypeMap );
                aIQIs( iIQI )->set_dof_type_list( tDofTypes );


                // set master dv dependencies
                moris::Cell< moris::Cell< moris::MSI::Dv_Type > > tDvTypes;
                string_to_cell_of_cell( aParameterList( iIQI ).get< std::string >( "master_dv_dependencies" ),
                                        tDvTypes,
                                        aMSIDvTypeMap );
                aIQIs( iIQI )->set_dv_type_list( tDvTypes );

                // set master properties
                moris::Cell< moris::Cell< std::string > > tMasterPropertyNamesPair;
                string_to_cell_of_cell( aParameterList( iIQI ).get< std::string >( "master_properties" ), tMasterPropertyNamesPair );

                for( uint iProp = 0; iProp < tMasterPropertyNamesPair.size(); iProp++ )
                {
                    // get property index
                    uint tPropertyIndex = aPropertyMap[ tMasterPropertyNamesPair( iProp )( 0 ) ];

                    // set property for IWG
                    aIQIs( iIQI )->set_property( aProperties( tPropertyIndex ), tMasterPropertyNamesPair( iProp )( 1 ) );
                }

                // set master constitutive models
                moris::Cell< moris::Cell< std::string > > tMasterCMNamesPair;
                string_to_cell_of_cell( aParameterList( iIQI ).get< std::string >( "master_constitutive_models" ), tMasterCMNamesPair );

                for( uint iCM = 0; iCM < tMasterCMNamesPair.size(); iCM++ )
                {
                    // get CM index
                    uint tCMIndex = aCMMap[ tMasterCMNamesPair( iCM )( 0 ) ];

                    // set CM for IWG
                    aIQIs( iIQI )->set_constitutive_model( aCMs( tCMIndex ), tMasterCMNamesPair( iCM )( 1 ) );
                }

                // set stabilization parameters
                moris::Cell< moris::Cell< std::string > > tSPNamesPair;
                string_to_cell_of_cell( aParameterList( iIQI ).get< std::string >( "stabilization_parameters" ), tSPNamesPair );

                for( uint iSP = 0; iSP < tSPNamesPair.size(); iSP++ )
                {
                    // get SP index
                    uint tSPIndex = aSPMap[ tSPNamesPair( iSP )( 0 ) ];

                    // set SP for IWG
                    aIQIs( iIQI )->set_stabilization_parameter( aSPs( tSPIndex ), tSPNamesPair( iSP )( 1 ) );
                }
                // debug
                aIQIs( iIQI )->print_names();
            }
        }

//------------------------------------------------------------------------------

    } /* namespace mdl */
} /* namespace moris */
