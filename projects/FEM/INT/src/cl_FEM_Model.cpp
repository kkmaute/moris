
// added by christian: link to Google Perftools
#ifdef WITHGPERFTOOLS
#include <gperftools/profiler.h>
#endif

//CHR/src
#include "cl_Stopwatch.hpp"
//LINALG/src
#include "cl_Map.hpp"
#include "fn_unique.hpp"
#include "fn_sum.hpp" // for check
#include "fn_iscol.hpp"
#include "fn_trans.hpp"
#include "op_equal_equal.hpp"
//MTK/src
#include "MTK_Tools.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
//FEM/INT/src
#include "cl_FEM_Node_Base.hpp"
#include "cl_FEM_Node.hpp"
#include "cl_FEM_Enums.hpp"
#include "cl_FEM_Model.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_CM_Factory.hpp"
#include "cl_FEM_SP_Factory.hpp"
#include "cl_FEM_IWG_Factory.hpp"
#include "cl_FEM_IQI_Factory.hpp"
//FEM/MSI/src
#include "cl_MSI_Equation_Object.hpp"
#include "cl_MSI_Dof_Type_Enums.hpp"
//FEM/VIS/src
#include "cl_VIS_Output_Enums.hpp"
//GEN/src
#include "cl_GEN_Pdv_Enums.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        FEM_Model::FEM_Model(
                mtk::Mesh_Manager *                 aMeshManager,
                const moris_index                 & aMeshPairIndex,
                moris::Cell< fem::Set_User_Info > & aSetInfo )
        : mMeshManager( aMeshManager ),
          mMeshPairIndex( aMeshPairIndex )
        {
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 0: unpack mesh
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // get pointers to interpolation and integration meshes
            mtk::Interpolation_Mesh* tIPMesh = nullptr;
            mtk::Integration_Mesh*   tIGMesh = nullptr;

            mMeshManager->get_mesh_pair( mMeshPairIndex, tIPMesh, tIGMesh );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 1: create nodes
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // start timer
            tic tTimer1;

            // ask IP mesh about number of IP vertices on proc
            luint tNumIPNodes = tIPMesh->get_num_nodes();

            // set size for list IP nodes
            mIPNodes.resize( tNumIPNodes, nullptr );

            // loop over IP mesh vertices
            for( uint iNode = 0; iNode < tNumIPNodes; iNode++ )
            {
                // create a new IP Node
                mIPNodes( iNode ) = new fem::Node( &tIPMesh->get_mtk_vertex( iNode ) );
            }

            if( par_rank() == 0)
            {
                // stop timer
                real tElapsedTime = tTimer1.toc<moris::chronos::milliseconds>().wall;

                // print output
                MORIS_LOG_INFO( "Model: created %u FEM IP nodes in %5.3f seconds.",
                        ( unsigned int ) tNumIPNodes,
                        ( double ) tElapsedTime / 1000 );
            }

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 2: create sets
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // start timer
            tic tTimer2;

            // get the number of sets
            uint tNumFemSets = aSetInfo.size();

            // create equation sets
            mFemSets.resize( tNumFemSets, nullptr );

            // get number of IP cells
            uint tNumIPCells = tIPMesh->get_num_elems();

            // reserve size for list of equation objects
            mFemClusters.reserve( tNumIPCells );

            // loop over the used fem set
            for( luint iSet = 0; iSet < tNumFemSets; iSet++ )
            {
                // get the mesh set name
                std::string tMeshSetName = aSetInfo( iSet).get_mesh_set_name();

                moris_index tMeshSetIndex;
                if( tMeshSetName.size() > 0 )
                {
                    // get the mesh set index from its name
                    tMeshSetIndex = tIGMesh->get_set_index_by_name( tMeshSetName );
                }
                else
                {
                    tMeshSetIndex = aSetInfo( iSet ).get_mesh_index();
                }

                // fill the mesh set index to fem set index map
                mMeshSetToFemSetMap[ std::make_tuple(
                        tMeshSetIndex,
                        aSetInfo( iSet ).get_time_continuity(),
                        aSetInfo( iSet ).get_time_boundary() ) ] = iSet;

                // get the mesh set pointer
                moris::mtk::Set * tMeshSet = tIGMesh->get_set_by_index( tMeshSetIndex );

                // if non-empty mesh set
                if ( tMeshSet->get_num_clusters_on_set() !=0 )
                {
                    // create a fem set
                    mFemSets( iSet ) = new fem::Set( this, tMeshSet, aSetInfo( iSet ), mIPNodes );
                    mFemSets( iSet )->set_equation_model( this );
                }
                // if empty mesh set
                else
                {
                    // create an empty fem set
                    mFemSets( iSet ) = new fem::Set();
                    mFemSets( iSet )->set_equation_model( this );
                }

                // collect equation objects associated with the set
                mFemClusters.append( mFemSets( iSet )->get_equation_object_list() );
            }
            // shrink list to fit size
            mFemClusters.shrink_to_fit();

            // stop timer
            real tElapsedTime = tTimer2.toc<moris::chronos::milliseconds>().wall;

            uint tNumElements       = mFemClusters.size();
            uint tGlobalNumElements = sum_all(tNumElements);

            // print output
            MORIS_LOG_INFO( "FEM_Model: created %u FEM IP elements in %5.3f seconds.",
                    tGlobalNumElements ,
                    (real) tElapsedTime / 1000.0 );
        }

        //------------------------------------------------------------------------------

        FEM_Model::FEM_Model(
                mtk::Mesh_Manager                           * aMeshManager,
                const moris_index                           & aMeshPairIndex,
                moris::Cell< moris::Cell< ParameterList > >   aParameterList,
                std::shared_ptr< Library_IO >                 aLibrary )
        : mMeshManager( aMeshManager ),
          mMeshPairIndex( aMeshPairIndex ),
          mParameterList( aParameterList )
        {
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 0: unpack fem input and mesh
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // start timer
            tic tTimer0;

            // get pointers to interpolation and integration meshes
            mtk::Interpolation_Mesh* tIPMesh = nullptr;
            mtk::Integration_Mesh*   tIGMesh = nullptr;
            mMeshManager->get_mesh_pair( mMeshPairIndex, tIPMesh, tIGMesh );

            // set the space dimension
            mSpaceDim = tIPMesh->get_spatial_dim();

            // unpack the FEM inputs
            this->initialize( aLibrary );

            if( par_rank() == 0)
            {
                // stop timer
                real tElapsedTime = tTimer0.toc<moris::chronos::milliseconds>().wall;

                // print output
                MORIS_LOG_INFO( "FEM_Model: unpack FEM input in %5.3f seconds.",
                        ( double ) tElapsedTime / 1000 );
            }

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 1: create IP nodes
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // start timer
            tic tTimer1;

            // ask mesh about number of IP nodes on proc
            luint tNumIPNodes = tIPMesh->get_num_nodes();

            // create IP node objects
            mIPNodes.resize( tNumIPNodes, nullptr );

            for( uint iNode = 0; iNode < tNumIPNodes; iNode++ )
            {
                mIPNodes( iNode ) = new fem::Node( &tIPMesh->get_mtk_vertex( iNode ) );
            }

            // stop timer
            real tElapsedTime = tTimer1.toc<moris::chronos::milliseconds>().wall;

            uint tGlobalNumNodes = sum_all(tNumIPNodes);

            // print output
            MORIS_LOG_INFO( "FEM_Model: created %u FEM IP nodes in %5.3f seconds.",
                    tGlobalNumNodes,
                    ( double ) tElapsedTime / 1000 );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 2: create fem sets
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // start timer
            tic tTimer2;

            // get number of fem sets
            uint tNumFemSets = mSetInfo.size();

            // set size for list of equation sets
            mFemSets.resize( tNumFemSets, nullptr );

            // get number of IP cells
            uint tNumIPCells = tIPMesh->get_num_elems();

            // reserve size for list of equation objects
            mFemClusters.reserve( tNumIPCells );

            // loop over the used fem set
            for( luint iSet = 0; iSet < tNumFemSets; iSet++ )
            {
                // get the mesh set name
                std::string tMeshSetName = mSetInfo( iSet ).get_mesh_set_name();

                // get the mesh set index from its name
                moris_index tMeshSetIndex = tIGMesh->get_set_index_by_name( tMeshSetName );

                // fill the mesh set index to fem set index map
                mMeshSetToFemSetMap[ std::make_tuple(
                        tMeshSetIndex,
                        mSetInfo( iSet ).get_time_continuity(),
                        mSetInfo( iSet ).get_time_boundary() ) ] = iSet;

                // get the mesh set pointer
                moris::mtk::Set * tMeshSet = tIGMesh->get_set_by_index( tMeshSetIndex );

                // if non-empty mesh set
                if ( tMeshSet->get_num_clusters_on_set() != 0 )
                {
                    // create new fem set
                    mFemSets( iSet ) = new fem::Set( this, tMeshSet, mSetInfo( iSet ), mIPNodes );

                    mFemSets( iSet )->set_equation_model( this );
                }
                // if empty mesh set
                else
                {
                    // create an empty fem set
                    mFemSets( iSet ) = new fem::Set();

                    mFemSets( iSet )->set_equation_model( this );
                }

                // collect equation objects associated with the set
                mFemClusters.append( mFemSets( iSet )->get_equation_object_list() );
            }

            // shrink to fit
            mFemClusters.shrink_to_fit();

            // stop timer
            tElapsedTime = tTimer2.toc<moris::chronos::milliseconds>().wall;

            uint tNumElements       = mFemClusters.size();
            uint tGlobalNumElements = sum_all(tNumElements);

            // print output
            MORIS_LOG_INFO( "FEM_Model: created %u FEM IP elements in %5.3f seconds.",
                    tGlobalNumElements,
                    ( double ) tElapsedTime / 1000.0 );
        }

        //------------------------------------------------------------------------------

        void FEM_Model::initialize( std::shared_ptr< Library_IO > aLibrary )
        {
            // get msi string to dof type map
            moris::map< std::string, MSI::Dof_Type > tMSIDofTypeMap =
                    moris::MSI::get_msi_dof_type_map();

            // get string to dv type map
            moris::map< std::string, PDV_Type > tMSIDvTypeMap =
                    get_pdv_type_map();

            switch( mParameterList.size() )
            {
                // without phase
                case 6:
                {
                    // create properties
                    std::map< std::string, uint > tPropertyMap;
                    this->create_properties( tPropertyMap, tMSIDofTypeMap, tMSIDvTypeMap, aLibrary );

                    // create constitutive models
                    std::map< std::string, uint > tCMMap;
                    this->create_constitutive_models( tCMMap, tPropertyMap, tMSIDofTypeMap, tMSIDvTypeMap );

                    // create stabilization parameters
                    std::map< std::string, uint > tSPMap;
                    this->create_stabilization_parameters( tSPMap, tPropertyMap, tCMMap, tMSIDofTypeMap, tMSIDvTypeMap );

                    // create IWGs
                    std::map< std::string, uint > tIWGMap;
                    this->create_IWGs( tIWGMap, tPropertyMap, tCMMap, tSPMap, tMSIDofTypeMap, tMSIDvTypeMap );

                    // create IQIs
                    std::map< std::string, uint > tIQIMap;
                    this->create_IQIs( tIQIMap, tPropertyMap, tCMMap, tSPMap, tMSIDofTypeMap, tMSIDvTypeMap );

                    // create FEM set info
                    this->create_fem_set_info();

                    break;
                }
                // with phase
                case 7:
                {
                    // create phases
                    this->create_phases();

                    // create properties
                    std::map< std::string, uint > tPropertyMap;
                    this->create_properties( tPropertyMap, tMSIDofTypeMap, tMSIDvTypeMap, aLibrary );

                    // create constitutive models
                    this->create_constitutive_models( tPropertyMap, tMSIDofTypeMap, tMSIDvTypeMap );

                    // create stabilization parameters
                    std::map< std::string, uint > tSPMap;
                    this->create_stabilization_parameters( tSPMap, tPropertyMap, tMSIDofTypeMap, tMSIDvTypeMap );

                    // create IWGs
                    this->create_IWGs( tPropertyMap, tSPMap, tMSIDofTypeMap );

                    // create IQIs
                    this->create_IQIs( tPropertyMap, tSPMap, tMSIDofTypeMap );

                    // create FEM set info
                    this->create_fem_set_info( true );

                    // get fem computation type parameter list
                    ParameterList tComputationParameterList = mParameterList( 5 )( 0 );

                    // get bool for printing physics model
                    bool tPrintPhysics =
                            tComputationParameterList.get< bool >( "print_physics_model" );

                    // if print FEM model
                    if( tPrintPhysics )
                    {
                        // phase info
                        std::cout<<"Phase info "<<std::endl;

                        // loop over phase info
                        for( uint iPhase = 0; iPhase < mPhaseInfo.size(); iPhase++ )
                        {
                            std::cout<<"%-------------------------------------------------"<<std::endl;
                            mPhaseInfo( iPhase ).print_names();
                            std::cout<<"%-------------------------------------------------"<<std::endl;
                        }

                        std::cout<<" "<<std::endl;

                        // set info
                        std::cout<<"Set info "<<std::endl;

                        // loop over set info
                        for( uint iSet = 0; iSet < mSetInfo.size(); iSet++ )
                        {
                            std::cout<<"%-------------------------------------------------"<<std::endl;
                            mSetInfo( iSet ).print_names();
                            std::cout<<"%-------------------------------------------------"<<std::endl;
                        }
                    }

                    break;
                }

                default:
                    MORIS_ERROR( false, "FEM_Model::initialize - wrong size for parameter list" );
            }
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

        void FEM_Model::finalize_equation_sets(
                MSI::Model_Solver_Interface * aModelSolverInterface )
        {
            // loop over the fem sets
            for( MSI::Equation_Set * tFemSet : mFemSets )
            {
                // finalize the fem set
                tFemSet->finalize( aModelSolverInterface );
            }
        }

        //------------------------------------------------------------------------------

        void FEM_Model::create_properties(
                std::map< std::string, uint >            & aPropertyMap,
                moris::map< std::string, MSI::Dof_Type > & aMSIDofTypeMap,
                moris::map< std::string, PDV_Type >      & aDvTypeMap,
                std::shared_ptr< Library_IO >              aLibrary )
        {
            // get the property parameter list
            moris::Cell< ParameterList > tPropParameterList = mParameterList( 0 );

            // get the number of properties
            uint tNumProps = tPropParameterList.size();

            // create a list of property pointers
            mProperties.resize( tNumProps, nullptr );

            // loop over the parameter lists
            for ( uint iProp = 0; iProp < tNumProps; iProp++ )
            {
                // get property parameter list
                ParameterList tPropParameter = tPropParameterList( iProp );

                // get property name from parameter list
                std::string tPropertyName = tPropParameter.get< std::string >( "property_name" );

                // create a property pointer
                mProperties( iProp ) = std::make_shared< fem::Property >();

                // set a name for the property
                mProperties( iProp )->set_name( tPropertyName );

                // fill property map
                aPropertyMap[ tPropertyName ] = iProp;

                // set dof dependencies
                moris::Cell< moris::Cell< moris::MSI::Dof_Type > > tDofTypes;
                string_to_cell_of_cell(
                        tPropParameter.get< std::string >( "dof_dependencies" ),
                        tDofTypes,
                        aMSIDofTypeMap );
                mProperties( iProp )->set_dof_type_list( tDofTypes );

                // set dof dependencies
                moris::Cell< moris::Cell< PDV_Type > > tDvTypes;
                string_to_cell_of_cell(
                        tPropParameter.get< std::string >( "dv_dependencies" ),
                        tDvTypes,
                        aDvTypeMap );
                mProperties( iProp )->set_dv_type_list( tDvTypes );

                // set function parameters
                moris::Cell< moris::Matrix< DDRMat > > tFuncParameters;
                string_to_cell_mat_2(
                        tPropParameter.get< std::string >( "function_parameters" ),
                        tFuncParameters );
                mProperties( iProp )->set_parameters( tFuncParameters );

                // set value function for property
                std::string tValFuncName = tPropParameter.get< std::string >( "value_function" );
                MORIS_FEM_FREE_FUNCTION tValFunction = nullptr;
                if ( tValFuncName.size() > 1 )
                {
                    tValFunction = aLibrary->load_fem_free_functions( tValFuncName );
                    mProperties( iProp )->set_val_function( tValFunction );
                }

                // set dof derivative function for property
                moris::Cell< std::string > tDofDerFuncNames;
                string_to_cell(
                        tPropParameter.get< std::string >( "dof_derivative_functions" ),
                        tDofDerFuncNames );
                uint tNumDofDerFuncs = tDofDerFuncNames.size();
                moris::Cell< fem::PropertyFunc > tDofDerFunctions( tNumDofDerFuncs, nullptr );
                for( uint iFunc = 0; iFunc < tNumDofDerFuncs; iFunc++ )
                {
                    if( tDofDerFuncNames( iFunc ).size() > 1 )
                    {
                        MORIS_FEM_FREE_FUNCTION tValFunction =
                                aLibrary->load_fem_free_functions( tDofDerFuncNames( iFunc ) );
                        tDofDerFunctions( iFunc ) = tValFunction;
                    }
                }
                mProperties( iProp )->set_dof_derivative_functions( tDofDerFunctions );

                // set dv derivative function for property
                moris::Cell< std::string > tDvDerFuncNames;
                string_to_cell(
                        tPropParameter.get< std::string >( "dv_derivative_functions" ),
                        tDvDerFuncNames );
                uint tNumDvDerFuncs = tDvDerFuncNames.size();
                moris::Cell< fem::PropertyFunc > tDvDerFunctions( tNumDvDerFuncs, nullptr );
                for( uint iFunc = 0; iFunc < tNumDvDerFuncs; iFunc++ )
                {
                    if( tDvDerFuncNames( iFunc ).size() > 1 )
                    {
                        MORIS_FEM_FREE_FUNCTION tValFunction =
                                aLibrary->load_fem_free_functions( tDvDerFuncNames( iFunc ) );
                        tDvDerFunctions( iFunc ) = tValFunction;
                    }
                }
                mProperties( iProp )->set_dv_derivative_functions( tDvDerFunctions );
            }
        }

        //------------------------------------------------------------------------------

        void FEM_Model::create_constitutive_models(
                std::map< std::string, uint >            & aPropertyMap,
                moris::map< std::string, MSI::Dof_Type > & aMSIDofTypeMap,
                moris::map< std::string, PDV_Type >      & aDvTypeMap )
        {
            // create a constitutive model factory
            CM_Factory tCMFactory;

            // get the CM parameter list
            moris::Cell< ParameterList > tCMParameterList = mParameterList( 1 );

            // get number of constitutive models
            uint tNumCMs = tCMParameterList.size();

            // loop over the parameter lists for CM
            for( uint iCM = 0; iCM < tNumCMs; iCM++ )
            {
                // get the treated CM parameter list
                ParameterList tCMParameter = tCMParameterList( iCM );

                // get the constitutive type from parameter list
                fem::Constitutive_Type tCMType =
                        static_cast< fem::Constitutive_Type >( tCMParameter.get< uint >( "constitutive_type" ) );

                // get the constitutive model name from parameter list
                std::string tCMName =
                        tCMParameter.get< std::string >( "constitutive_name" );

                // get the phase from parameter list
                std::string tPhaseName =
                        tCMParameter.get< std::string >( "phase_name" );

                // get the model type
                fem::Model_Type tCMModelType =
                        static_cast< fem::Model_Type >( tCMParameter.get< uint >( "model_type" ) );

                // create a constitutive model pointer
                std::shared_ptr< fem::Constitutive_Model > tCM =
                        tCMFactory.create_CM( tCMType );

                // set CM name
                tCM->set_name( tCMName );

                // set CM space dimension
                tCM->set_space_dim( mSpaceDim );

                // set CM model type
                if( tCMModelType != fem::Model_Type::UNDEFINED )
                {
                    tCM->set_model_type( tCMModelType );
                }

                // set CM dof dependencies
                moris::Cell< moris::Cell< moris::MSI::Dof_Type > > tDofTypes;
                string_to_cell_of_cell(
                        std::get< 0 >( tCMParameter.get< std::pair< std::string, std::string > >( "dof_dependencies" ) ),
                        tDofTypes,
                        aMSIDofTypeMap );
                moris::Cell< std::string > tDofTypeNames;
                string_to_cell(
                        std::get< 1 >( tCMParameter.get< std::pair< std::string, std::string > >( "dof_dependencies" ) ),
                        tDofTypeNames );
                tCM->set_dof_type_list( tDofTypes, tDofTypeNames );

                // set CM dv dependencies
                moris::Cell< moris::Cell< PDV_Type > > tDvTypes;
                string_to_cell_of_cell(
                        std::get< 0 >( tCMParameter.get< std::pair< std::string, std::string > >( "dv_dependencies" ) ),
                        tDvTypes,
                        aDvTypeMap );
                moris::Cell< std::string > tDvTypeNames;
                string_to_cell(
                        std::get< 1 >( tCMParameter.get< std::pair< std::string, std::string > >( "dv_dependencies" ) ),
                        tDvTypeNames );
                tCM->set_dv_type_list( tDvTypes, tDvTypeNames );

                // set CM properties
                moris::Cell< moris::Cell< std::string > > tPropertyNamesPair;
                string_to_cell_of_cell(
                        tCMParameter.get< std::string >( "properties" ),
                        tPropertyNamesPair );
                for( uint iProp = 0; iProp < tPropertyNamesPair.size(); iProp++ )
                {
                    // get the property name
                    std::string tPropertyName = tPropertyNamesPair( iProp )( 0 );

                    // check if property in the map
                    MORIS_ERROR( aPropertyMap.find( tPropertyName ) != aPropertyMap.end(),
                            "FEM_Model::create_CMs - Unknown aPropertyString : %s \n",
                            tPropertyName.c_str() );

                    // get property index
                    uint tPropertyIndex = aPropertyMap[ tPropertyName ];

                    // set property for CM
                    tCM->set_property(
                            mProperties( tPropertyIndex ),
                            tPropertyNamesPair( iProp )( 1 ) );
                }

                // set local properties
                tCM->set_local_properties();

                // check the phase exist
                MORIS_ERROR( mPhaseMap.find( tPhaseName ) != mPhaseMap.end(),
                        "FEM_Model::create_constitutive_models - Unknown tPhaseName : %s \n",
                        tPhaseName.c_str() );

                // set CM to corresponding phase
                mPhaseInfo( mPhaseMap[ tPhaseName ] ).set_CM( tCM );
            }
        }

        //------------------------------------------------------------------------------

        void FEM_Model::create_stabilization_parameters(
                std::map< std::string, uint >            & aSPMap,
                std::map< std::string, uint >            & aPropertyMap,
                moris::map< std::string, MSI::Dof_Type > & aMSIDofTypeMap,
                moris::map< std::string, PDV_Type >      & aDvTypeMap )
        {
            // create a stabilization parameter factory
            SP_Factory tSPFactory;

            // get the SP parameter list
            moris::Cell< ParameterList > tSPParameterList = mParameterList( 2 );

            // get the number of stabilization parameters
            uint tNumSPs = tSPParameterList.size();

            // set size for the list of stabilization parameter pointer
            mSPs.resize( tNumSPs, nullptr );

            // loop over the parameter list
            for ( uint iSP = 0; iSP < tNumSPs; iSP++ )
            {
                // get the stabilization parameters
                ParameterList tSPParameter = tSPParameterList( iSP );

                // get the stabilization parameter name
                std::string tSPName = tSPParameter.get< std::string >( "stabilization_name" );

                // get the stabilization type from parameter list
                fem::Stabilization_Type tSPType =
                        static_cast< fem::Stabilization_Type >( tSPParameter.get< uint >( "stabilization_type" ) );

                // create a stabilization parameter pointer
                mSPs( iSP ) = tSPFactory.create_SP( tSPType );

                // set name
                mSPs( iSP )->set_name( tSPName );

                // set SP space dimension
                mSPs( iSP )->set_space_dim( mSpaceDim );

                // fill stabilization map
                aSPMap[ tSPName ] = iSP;

                // set parameters
                moris::Cell< moris::Matrix< DDRMat > > tFuncParameters;
                string_to_cell_mat_2(
                        tSPParameter.get< std::string >( "function_parameters" ),
                        tFuncParameters );
                mSPs( iSP )->set_parameters( tFuncParameters );

                // init string for master or slave
                std::string tIsMasterString = "master";
                mtk::Master_Slave tIsMaster = mtk::Master_Slave::MASTER;

                // loop on master and slave
                for( uint iMaster = 0; iMaster <= mSPs( iSP )->get_has_slave(); iMaster++ )
                {
                    // if slave
                    if( iMaster )
                    {
                        // reset string for slave
                        tIsMasterString = "slave";
                        tIsMaster = mtk::Master_Slave::SLAVE;
                    }

                    // get the treated IWG phase
                    std::string tPhaseName =
                            tSPParameter.get< std::string >( tIsMasterString + "_phase_name" );

                    // check for unknown phase
                    MORIS_ERROR( mPhaseMap.find( tPhaseName ) != mPhaseMap.end(),
                            "FEM_Model::create_stabilization_parameters - Unknown phase name: %s \n",
                            tPhaseName.c_str() );

                    // set dof dependencies
                    moris::Cell< moris::Cell< moris::MSI::Dof_Type > > tDofTypes;
                    string_to_cell_of_cell(
                            std::get< 0 >( tSPParameter.get< std::pair< std::string, std::string > >( tIsMasterString + "_dof_dependencies" ) ),
                            tDofTypes,
                            aMSIDofTypeMap );
                    moris::Cell< std::string > tDofTypeNames;
                    string_to_cell( std::get< 1 >(
                            tSPParameter.get< std::pair< std::string, std::string > >( tIsMasterString + "_dof_dependencies" ) ),
                            tDofTypeNames );
                    mSPs( iSP )->set_dof_type_list( tDofTypes, tDofTypeNames, tIsMaster );

                    // set dv dependencies
                    moris::Cell< moris::Cell< PDV_Type > > tDvTypes;
                    string_to_cell_of_cell(
                            std::get< 0 >( tSPParameter.get< std::pair< std::string, std::string > >( tIsMasterString + "_dv_dependencies" ) ),
                            tDvTypes,
                            aDvTypeMap );
                    moris::Cell< std::string > tDvTypeNames;
                    string_to_cell(
                            std::get< 1 >( tSPParameter.get< std::pair< std::string, std::string > >( tIsMasterString + "_dv_dependencies" ) ),
                            tDvTypeNames );
                    mSPs( iSP )->set_dv_type_list( tDvTypes, tDvTypeNames, tIsMaster );

                    // set master properties
                    moris::Cell< moris::Cell< std::string > > tPropertyNamesPair;
                    string_to_cell_of_cell(
                            tSPParameter.get< std::string >( tIsMasterString + "_properties" ),
                            tPropertyNamesPair );

                    for( uint iProp = 0; iProp < tPropertyNamesPair.size(); iProp++ )
                    {
                        // get the property name
                        std::string tPropertyName = tPropertyNamesPair( iProp )( 0 );

                        // check for unknown property
                        MORIS_ERROR( aPropertyMap.find( tPropertyName ) != aPropertyMap.end(),
                                "FEM_Model::create_stabilization_parameters - Unknown %s aPropertyString : %s \n",
                                tIsMasterString.c_str(),
                                tPropertyName.c_str() );

                        // get property index
                        uint tPropertyIndex = aPropertyMap[ tPropertyName ];

                        // set property for CM
                        mSPs( iSP )->set_property(
                                mProperties( tPropertyIndex ),
                                tPropertyNamesPair( iProp )( 1 ),
                                tIsMaster );
                    }

                    // set constitutive models
                    moris::Cell< moris::Cell< std::string > > tCMNamesPair;
                    string_to_cell_of_cell(
                            tSPParameter.get< std::string >( tIsMasterString + "_constitutive_models" ),
                            tCMNamesPair );

                    // loop over CM names
                    for( uint iCM = 0; iCM < tCMNamesPair.size(); iCM++ )
                    {
                        // get the CM name
                        std::string tCMName = tCMNamesPair( iCM )( 0 );

                        // get CM from phase
                        std::shared_ptr< fem::Constitutive_Model > tCM =
                                mPhaseInfo( mPhaseMap[ tPhaseName ] ).get_CM_by_name( tCMName );

                        // set CM for SP
                        mSPs( iSP )->set_constitutive_model(
                                tCM,
                                tCMNamesPair( iCM )( 1 ) );
                    }
                }
            }
        }

        //------------------------------------------------------------------------------

        void FEM_Model::create_IWGs(
                std::map< std::string, uint >            & aPropertyMap,
                std::map< std::string, uint >            & aSPMap,
                moris::map< std::string, MSI::Dof_Type > & aMSIDofTypeMap )
        {
            // create an IWG factory
            IWG_Factory tIWGFactory;

            // get the IWG parameter list
            moris::Cell< ParameterList > tIWGParameterList = mParameterList( 3 );

            // get number of IWGs
            uint tNumIWGs = tIWGParameterList.size();

            // create a list of IWG pointers
            mIWGs.resize( tNumIWGs, nullptr );

            // loop over the parameter lists
            for( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
            {
                // get the treated IWG parameter list
                ParameterList tIWGParameter = tIWGParameterList( iIWG );

                // get the treated IWG name
                std::string tIWGName = tIWGParameter.get< std::string >( "IWG_name" );

                // get the treated IWG type
                fem::IWG_Type tIWGType =
                        static_cast< fem::IWG_Type >( tIWGParameter.get< uint >( "IWG_type" ) );

                // get the ghost order from parameter list
                uint tGhostOrder = tIWGParameter.get< uint >( "ghost_order" );

                // get the treated IWG residual dof type
                moris::Cell< moris::MSI::Dof_Type > tResDofTypes;
                string_to_cell(
                        tIWGParameter.get< std::string >( "dof_residual" ),
                        tResDofTypes,
                        aMSIDofTypeMap );

                // get the treated IWG bulk type
                fem::Element_Type tIWGBulkType =
                        static_cast< fem::Element_Type >( tIWGParameter.get< uint >( "IWG_bulk_type" ) );

                //
                bool tMasterSlave = ( tIWGBulkType == fem::Element_Type::DOUBLE_SIDESET );

                // create an IWG pointer
                mIWGs( iIWG ) = tIWGFactory.create_IWG( tIWGType );

                // set name
                mIWGs( iIWG )->set_name( tIWGName );

                // set interpolation order
                mIWGs( iIWG )->set_interpolation_order( tGhostOrder );

                // set residual dof type
                mIWGs( iIWG )->set_residual_dof_type( tResDofTypes );

                // set bulk type
                mIWGs( iIWG )->set_bulk_type( tIWGBulkType );

                // init string for master or slave
                std::string tIsMasterString = "master";
                mtk::Master_Slave tIsMaster = mtk::Master_Slave::MASTER;

                // loop on master and slave
                for( uint iMaster = 0; iMaster <= tMasterSlave; iMaster++ )
                {
                    // if slave
                    if( iMaster )
                    {
                        // reset string for slave
                        tIsMasterString = "slave";
                        tIsMaster = mtk::Master_Slave::SLAVE;
                    }

                    // get the treated IWG phase
                    std::string tPhaseName =
                            tIWGParameter.get< std::string >( tIsMasterString + "_phase_name" );

                    // check for unknown phase
                    MORIS_ERROR( mPhaseMap.find( tPhaseName ) != mPhaseMap.end(),
                            "FEM_Model::create_IWGs - Unknown phase name: %s \n",
                            tPhaseName.c_str() );

                    // set phase name
                    mIWGs( iIWG )->set_phase_name( tPhaseName, tIsMaster );

                    // get phase index
                    uint tPhaseIndex = mPhaseMap[ tPhaseName ];

                    // get dof type list from phase
                    mPhaseInfo( tPhaseIndex ).add_dof_type_to_list( tResDofTypes );

                    // set properties
                    moris::Cell< moris::Cell< std::string > > tPropertyNamesPair;
                    string_to_cell_of_cell(
                            tIWGParameter.get< std::string >( tIsMasterString + "_properties" ),
                            tPropertyNamesPair );

                    for( uint iProp = 0; iProp < tPropertyNamesPair.size(); iProp++ )
                    {
                        // get property name
                        std::string tPropertyName = tPropertyNamesPair( iProp )( 0 );

                        // check for unknown property
                        MORIS_ERROR( aPropertyMap.find( tPropertyName ) != aPropertyMap.end() ,
                                "FEM_Model::create_IWGs - Unknown %s aPropertyString: %s \n",
                                tIsMasterString.c_str(),
                                tPropertyName.c_str() );

                        // get property index
                        uint tPropertyIndex = aPropertyMap[ tPropertyName ];

                        // set property for IWG
                        mIWGs( iIWG )->set_property(
                                mProperties( tPropertyIndex ),
                                tPropertyNamesPair( iProp )( 1 ),
                                tIsMaster );
                    }

                    // set constitutive models
                    moris::Cell< moris::Cell< std::string > > tCMNamesPair;
                    string_to_cell_of_cell(
                            tIWGParameter.get< std::string >( tIsMasterString + "_constitutive_models" ),
                            tCMNamesPair );

                    // loop over constitutive models
                    for( uint iCM = 0; iCM < tCMNamesPair.size(); iCM++ )
                    {
                        // get the CM name
                        std::string tCMName = tCMNamesPair( iCM )( 0 );

                        // get CM from phase
                        std::shared_ptr< fem::Constitutive_Model > tCM =
                                mPhaseInfo( tPhaseIndex ).get_CM_by_name( tCMName );

                        // set CM for IWG
                        mIWGs( iIWG )->set_constitutive_model(
                                tCM,
                                tCMNamesPair( iCM )( 1 ),
                                tIsMaster );
                    }
                }

                // set stabilization parameters
                moris::Cell< moris::Cell< std::string > > tSPNamesPair;
                string_to_cell_of_cell(
                        tIWGParameter.get< std::string >( "stabilization_parameters" ),
                        tSPNamesPair );

                // loop over SP names
                for( uint iSP = 0; iSP < tSPNamesPair.size(); iSP++ )
                {
                    // get the SP name
                    std::string tSPName = tSPNamesPair( iSP )( 0 );

                    // check for unknown SP
                    MORIS_ERROR( aSPMap.find( tSPName ) != aSPMap.end(),
                            "FEM_Model::create_IWGs - Unknown aSPString: %s \n",
                            tSPName.c_str() );

                    // get SP index
                    uint tSPIndex = aSPMap[ tSPNamesPair( iSP )( 0 ) ];

                    // set SP for IWG
                    mIWGs( iIWG )->set_stabilization_parameter(
                            mSPs( tSPIndex ),
                            tSPNamesPair( iSP )( 1 ) );
                }
            }

            // loop over the parameter lists to set dof dependencies
            for( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
            {
                // get the treated IWG parameter list
                ParameterList tIWGParameter = tIWGParameterList( iIWG );

                // get the IWG bulk type
                fem::Element_Type tIWGBulkType = mIWGs( iIWG )->get_bulk_type();

                // get the IWG master phase name
                std::string tPhaseName =
                        mIWGs( iIWG )->get_phase_name( mtk::Master_Slave::MASTER );

                // check for unknown phase
                MORIS_ERROR( mPhaseMap.find( tPhaseName ) != mPhaseMap.end(),
                        "FEM_Model::create_IWGs - Unknown phase name: %s \n",
                        tPhaseName.c_str() );

                // get the phase index
                uint tPhaseIndex = mPhaseMap[ tPhaseName ];

                // get dof type list from phase
                moris::Cell< moris::Cell< MSI::Dof_Type > > tMasterSlaveDofTypes =
                        mPhaseInfo( tPhaseIndex ).get_dof_type_list();

                // get dof type list from phase
                moris::Cell< moris::Cell< PDV_Type > > tMasterSlavePdvTypes =
                        mPhaseInfo( tPhaseIndex ).get_dv_type_list();

                // get dof type check list
                Matrix< DDSMat > tMasterDofCheck =
                        mPhaseInfo( tPhaseIndex ).get_dof_type_check_list();

                // get pdv type check list
                Matrix< DDSMat > tMasterPdvCheck =
                        mPhaseInfo( tPhaseIndex ).get_pdv_type_check_list();

                if( tIWGBulkType == fem::Element_Type::DOUBLE_SIDESET )
                {
                    // get the IWG slave phase name
                    std::string tSlavePhaseName =
                            mIWGs( iIWG )->get_phase_name( mtk::Master_Slave::SLAVE );

                    // check for unknown phase
                    MORIS_ERROR( mPhaseMap.find( tSlavePhaseName ) != mPhaseMap.end() ,
                            "FEM_Model::create_IWGs - Unknown phase name: %s \n",
                            tSlavePhaseName.c_str() );

                    // get CM index
                    uint tSlavePhaseIndex = mPhaseMap[ tSlavePhaseName ];

                    // get dof type list from phase
                    moris::Cell< moris::Cell< MSI::Dof_Type > > tSlaveDofTypes =
                            mPhaseInfo( tSlavePhaseIndex ).get_dof_type_list();

                    // get pdv type list from phase
                    moris::Cell< moris::Cell< PDV_Type > > tSlavePdvTypes =
                            mPhaseInfo( tSlavePhaseIndex ).get_dv_type_list();

                    for ( uint iSlaveDof = 0; iSlaveDof < tSlaveDofTypes.size(); iSlaveDof++ )
                    {
                        // get the dof type index in enum list
                        uint tDofIndex = static_cast< uint >( tSlaveDofTypes( iSlaveDof )( 0 ) );

                        // if dof not yet in dof type list
                        if( tMasterDofCheck( tDofIndex ) == -1 )
                        {
                            // add dof to dof type list
                            tMasterSlaveDofTypes.push_back( tSlaveDofTypes( iSlaveDof ) );
                        }
                    }

                    for ( uint iSlavePdv = 0; iSlavePdv < tSlavePdvTypes.size(); iSlavePdv++ )
                    {
                        // get the pdv type index in enum list
                        uint tPdvIndex = static_cast< uint >( tSlavePdvTypes( iSlavePdv )( 0 ) );

                        // if pdv not yet in dof type list
                        if( tMasterPdvCheck( tPdvIndex ) == -1 )
                        {
                            // add pdv to pdv type list
                            tMasterSlavePdvTypes.push_back( tSlavePdvTypes( iSlavePdv ) );
                        }
                    }

                    // set slave dof dependencies
                    mIWGs( iIWG )->set_dof_type_list( tMasterSlaveDofTypes, mtk::Master_Slave::SLAVE );

                    // set slave dv dependencies
                    mIWGs( iIWG )->set_dv_type_list(  tMasterSlavePdvTypes, mtk::Master_Slave::SLAVE );
                }

                // set master dof dependencies
                mIWGs( iIWG )->set_dof_type_list( tMasterSlaveDofTypes, mtk::Master_Slave::MASTER );

                // set master dv dependencies
                mIWGs( iIWG )->set_dv_type_list( tMasterSlavePdvTypes, mtk::Master_Slave::MASTER );
            }
        }

        //------------------------------------------------------------------------------

        void FEM_Model::create_IQIs(
                std::map< std::string, uint >            & aPropertyMap,
                std::map< std::string, uint >            & aSPMap,
                moris::map< std::string, MSI::Dof_Type > & aMSIDofTypeMap )
        {
            // create an IQI factory
            IQI_Factory tIQIFactory;

            // get the IQI parameter list
            moris::Cell< ParameterList > tIQIParameterList = mParameterList( 4 );

            // get number of IQIs
            uint tNumIQIs = tIQIParameterList.size();

            // set size for list of IQI pointers
            mIQIs.resize( tNumIQIs, nullptr );

            // loop over the parameter lists
            for( uint iIQI = 0; iIQI < tNumIQIs; iIQI++ )
            {
                // get the treated IQI parameter list
                ParameterList tIQIParameter = tIQIParameterList( iIQI );

                // get the treated IQI name from parameter list
                std::string tIQIName =
                        tIQIParameter.get< std::string >( "IQI_name" );

                // get the IQI type from parameter list
                fem::IQI_Type tIQIType =
                        static_cast< fem::IQI_Type >( tIQIParameter.get< uint >( "IQI_type" ) );

                // get the quantity dof type from parameter list
                moris::Cell< moris::MSI::Dof_Type > tQuantityDofTypes;
                string_to_cell(
                        tIQIParameter.get< std::string >( "dof_quantity" ),
                        tQuantityDofTypes,
                        aMSIDofTypeMap );

                // get the field index from parameter list
                sint tIQIFieldIndex =
                        tIQIParameter.get< moris::sint >( "vectorial_field_index" );

                // set function parameters
                moris::Cell< moris::Matrix< DDRMat > > tFuncParameters;
                string_to_cell_mat_2(
                        tIQIParameter.get< std::string >( "function_parameters" ),
                        tFuncParameters );

                // get the treated IQI bulk type
                fem::Element_Type tIQIBulkType =
                        static_cast< fem::Element_Type >( tIQIParameter.get< uint >( "IQI_bulk_type" ) );

                // set bool to true if double sideset
                bool tMasterSlave = ( tIQIBulkType == fem::Element_Type::DOUBLE_SIDESET );

                // create an IQI pointer
                mIQIs( iIQI ) = tIQIFactory.create_IQI( tIQIType );

                // set name
                mIQIs( iIQI )->set_name( tIQIName );

                // set quantity dof type
                mIQIs( iIQI )->set_quantity_dof_type( tQuantityDofTypes );

                // set index for vectorial field
                mIQIs( iIQI )->set_output_type_index( tIQIFieldIndex );

                // set bulk type
                mIQIs( iIQI )->set_bulk_type( tIQIBulkType );

                // set constant parameters
                mIQIs( iIQI )->set_parameters( tFuncParameters );

                // init string for master or slave
                std::string tIsMasterString = "master";
                mtk::Master_Slave tIsMaster = mtk::Master_Slave::MASTER;

                // loop on master and slave
                for( uint iMaster = 0; iMaster <= tMasterSlave; iMaster++ )
                {
                    // if slave
                    if( iMaster )
                    {
                        // reset string for slave
                        tIsMasterString = "slave";
                        tIsMaster = mtk::Master_Slave::SLAVE;
                    }

                    // get the treated IWG phase
                    std::string tPhaseName = tIQIParameter.get< std::string >( tIsMasterString + "_phase_name" );

                    // check for unknown phase
                    MORIS_ERROR( mPhaseMap.find( tPhaseName ) != mPhaseMap.end() ,
                            "FEM_Model::create_IQIs - Unknown phase name: %s \n",
                            tPhaseName.c_str() );

                    // set phase name
                    mIQIs( iIQI )->set_phase_name( tPhaseName, tIsMaster );

                    // get the phase index
                    uint tPhaseIndex = mPhaseMap[ tPhaseName ];

                    // get dof type list from phase
                    moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypes =
                            mPhaseInfo( tPhaseIndex ).get_dof_type_list();

                    // get dof type list from phase
                    moris::Cell< moris::Cell< PDV_Type > > tDvTypes =
                            mPhaseInfo( tPhaseIndex ).get_dv_type_list();

                    // set master dof dependencies
                    mIQIs( iIQI )->set_dof_type_list( tDofTypes );

                    // set master dv dependencies
                    mIQIs( iIQI )->set_dv_type_list( tDvTypes );

                    // set master properties
                    moris::Cell< moris::Cell< std::string > > tPropertyNamesPair;
                    string_to_cell_of_cell(
                            tIQIParameter.get< std::string >( tIsMasterString + "_properties" ),
                            tPropertyNamesPair );

                    for( uint iProp = 0; iProp < tPropertyNamesPair.size(); iProp++ )
                    {
                        // get the property name
                        std::string tPropertyName = tPropertyNamesPair( iProp )( 0 );

                        // check for unknown property
                        MORIS_ERROR( aPropertyMap.find( tPropertyName ) != aPropertyMap.end() ,
                                "FEM_Model::create_IQIs - Unknown %s aPropertyString: %s \n",
                                tIsMasterString.c_str(),
                                tPropertyName.c_str() );

                        // get property index
                        uint tPropertyIndex = aPropertyMap[ tPropertyName ];

                        // set property for IWG
                        mIQIs( iIQI )->set_property(
                                mProperties( tPropertyIndex ),
                                tPropertyNamesPair( iProp )( 1 ),
                                tIsMaster );
                    }

                    // set master constitutive models
                    moris::Cell< moris::Cell< std::string > > tCMNamesPair;
                    string_to_cell_of_cell(
                            tIQIParameter.get< std::string >( tIsMasterString + "_constitutive_models" ),
                            tCMNamesPair );

                    for( uint iCM = 0; iCM < tCMNamesPair.size(); iCM++ )
                    {
                        // get CM name
                        std::string tCMName = tCMNamesPair( iCM )( 0 );

                        // get CM from phase
                        std::shared_ptr< fem::Constitutive_Model > tCM =
                                mPhaseInfo( tPhaseIndex ).get_CM_by_name( tCMName );

                        // set CM for IQI
                        mIQIs( iIQI )->set_constitutive_model(
                                tCM,
                                tCMNamesPair( iCM )( 1 ),
                                tIsMaster );
                    }
                }

                // set stabilization parameters
                moris::Cell< moris::Cell< std::string > > tSPNamesPair;
                string_to_cell_of_cell(
                        tIQIParameter.get< std::string >( "stabilization_parameters" ),
                        tSPNamesPair );

                for( uint iSP = 0; iSP < tSPNamesPair.size(); iSP++ )
                {
                    // get the SP name
                    std::string tSPName = tSPNamesPair( iSP )( 0 );

                    // check for unknown SP
                    MORIS_ERROR( aSPMap.find( tSPName ) != aSPMap.end(),
                            "FEM_Model::create_IQIs - Unknown aSPString: %s \n",
                            tSPName.c_str() );

                    // get SP index
                    uint tSPIndex = aSPMap[ tSPName ];

                    // set SP for IWG
                    mIQIs( iIQI )->set_stabilization_parameter(
                            mSPs( tSPIndex ),
                            tSPNamesPair( iSP )( 1 ) );
                }
            }
        }

        //------------------------------------------------------------------------------

        void FEM_Model::create_phases()
        {
            // get the phase parameter list
            moris::Cell< ParameterList > tPhaseParameterList = mParameterList( 6 );

            // get number of phases
            uint tNumPhases = tPhaseParameterList.size();

            // resize the list of phase user info
            mPhaseInfo.resize( tNumPhases );

            // loop over the parameter lists
            for( uint iPhase = 0; iPhase < tNumPhases; iPhase++ )
            {
                // get the treated phase parameter list
                ParameterList tPhaseParameter = tPhaseParameterList( iPhase );

                // get the phase name from parameter list
                std::string tPhaseName =
                        tPhaseParameter.get< std::string >( "phase_name" );

                // set phase name to phase
                mPhaseInfo( iPhase ).set_phase_name( tPhaseName );

                // get the phase index from parameter list
                moris::Matrix< moris::IndexMat > tPhaseIndices;
                string_to_mat( tPhaseParameter.get< std::string >( "phase_indices" ), tPhaseIndices );

                // set phase mesh indices to phase
                mPhaseInfo( iPhase ).set_phase_indices( tPhaseIndices );

                // fill phase map
                mPhaseMap[ tPhaseName ] = iPhase;
            }
        }

        //------------------------------------------------------------------------------

        void FEM_Model::create_fem_set_info( bool aWithPhase )
        {
            // init number of fem sets to be created
            uint tNumFEMSets = 0;

            // get the IWG and IQI parameter lists
            moris::Cell< ParameterList > tIWGParameterList = mParameterList( 3 );
            moris::Cell< ParameterList > tIQIParameterList = mParameterList( 4 );

            // get fem computation type parameter list
            ParameterList tComputationParameterList = mParameterList( 5 )( 0 );

            // get bool for analytical/finite differenec for SA
            bool tIsAnalyticalSA =
                    tComputationParameterList.get< bool >( "is_analytical_sensitivity" );

            // get enum for FD scheme
            fem::FDScheme_Type tFDSchemeForSA = static_cast< fem::FDScheme_Type >(
                    tComputationParameterList.get< uint >( "finite_difference_scheme" ) );

            // get perturbation size for FD
            real tFDPerturbation = tComputationParameterList.get< real >(
                    "finite_difference_perturbation_size" );

            // create a map of the set
            std::map< std::tuple< std::string, bool, bool >, uint > tMeshtoFemSet;

            // loop over the IWGs
            for( uint iIWG = 0; iIWG < tIWGParameterList.size(); iIWG++ )
            {
                // get the treated IWG parameter list
                ParameterList tIWGParameter = tIWGParameterList( iIWG );

                // get the IWG bulk type
                fem::Element_Type tIWGBulkType = mIWGs( iIWG )->get_bulk_type();

                // get time continuity flag
                bool tTimeContinuity = mIWGs( iIWG )->get_time_continuity();

                // get time boundary flag
                bool tTimeBoundary = mIWGs( iIWG )->get_time_boundary();

                // get bool for ghost
                bool tIsGhost = mIWGs( iIWG )->get_ghost_flag();

                // get the IWG master phase name
                std::string tMasterPhaseName =
                        mIWGs( iIWG )->get_phase_name( mtk::Master_Slave::MASTER );

                // get the IWG slave phase name
                std::string tSlavePhaseName =
                        mIWGs( iIWG )->get_phase_name( mtk::Master_Slave::SLAVE );

                // get slave phase string from IWG input
                std::string tSlavePhaseString =
                        tIWGParameter.get< std::string >( "neighbor_phases" );

                // get ordinal string from IWG input
                std::string tOrdinalString =
                        tIWGParameter.get< std::string >( "side_ordinals" );

                // get mesh set names for IWG
                moris::Cell< std::string > tMeshSetNames;
                this->get_mesh_set_names(
                        tIWGBulkType,
                        tMasterPhaseName,
                        tSlavePhaseName,
                        tSlavePhaseString,
                        tOrdinalString,
                        tIsGhost,
                        tMeshSetNames );

                // loop over the mesh set names
                for( uint iSetName = 0; iSetName < tMeshSetNames.size(); iSetName++ )
                {
                    // check if the mesh set name already in map
                    if( tMeshtoFemSet.find( std::make_tuple(
                            tMeshSetNames( iSetName ),
                            tTimeContinuity,
                            tTimeBoundary ) ) == tMeshtoFemSet.end() )
                    {
                        // add the mesh set name map
                        tMeshtoFemSet[ std::make_tuple(
                                tMeshSetNames( iSetName ),
                                tTimeContinuity,
                                tTimeBoundary ) ] = tNumFEMSets++;

                        // create a fem set info for the mesh set
                        Set_User_Info aSetUserInfo;

                        // set its mesh set name
                        aSetUserInfo.set_mesh_set_name( tMeshSetNames( iSetName ) );

                        // set its time continuity flag
                        aSetUserInfo.set_time_continuity( tTimeContinuity );

                        // set its time boundary flag
                        aSetUserInfo.set_time_boundary( tTimeBoundary );

                        // set its sensitivity analysis type flag
                        aSetUserInfo.set_is_analytical_sensitivity_analysis( tIsAnalyticalSA );

                        // set its FD scheme for sensitivity analysis
                        aSetUserInfo.set_finite_difference_scheme_for_sensitivity_analysis( tFDSchemeForSA );

                        // set its FD perturbation size for sensitivity analysis
                        aSetUserInfo.set_finite_difference_perturbation_size( tFDPerturbation );

                        // set the IWG
                        aSetUserInfo.set_IWG( mIWGs( iIWG ) );

                        // add it to the list of fem set info
                        mSetInfo.push_back( aSetUserInfo );
                    }
                    else
                    {
                        // set the IWG
                        mSetInfo( tMeshtoFemSet[ std::make_tuple(
                                tMeshSetNames( iSetName ),
                                tTimeContinuity,
                                tTimeBoundary ) ] ).set_IWG( mIWGs( iIWG ) );
                    }
                }
            }

            // loop over the IQIs
            for( uint iIQI = 0; iIQI < tIQIParameterList.size(); iIQI++ )
            {
                // get the treated IWG parameter list
                ParameterList tIQIParameter = tIQIParameterList( iIQI );

                // get the IWG bulk type
                fem::Element_Type tIQIBulkType = mIQIs( iIQI )->get_bulk_type();

                // get time continuity flag
                bool tTimeContinuity = mIQIs( iIQI )->get_time_continuity();

                // get time boundary flag
                bool tTimeBoundary = mIQIs( iIQI )->get_time_boundary();

                // get the IWG master phase name
                std::string tMasterPhaseName =
                        mIQIs( iIQI )->get_phase_name( mtk::Master_Slave::MASTER );

                // get the IWG slave phase name
                std::string tSlavePhaseName =
                        mIQIs( iIQI )->get_phase_name( mtk::Master_Slave::SLAVE );

                // get slave phase string from IQI input
                std::string tSlavePhaseString =
                        tIQIParameter.get< std::string >( "neighbor_phases" );

                // get ordinal string from IQI input
                std::string tOrdinalString =
                        tIQIParameter.get< std::string >( "side_ordinals" );

                // get mesh set names for IWG
                moris::Cell< std::string > tMeshSetNames;
                this->get_mesh_set_names(
                        tIQIBulkType,
                        tMasterPhaseName,
                        tSlavePhaseName,
                        tSlavePhaseString,
                        tOrdinalString,
                        false,
                        tMeshSetNames );

                // loop over the mesh set names
                for( uint iSetName = 0; iSetName < tMeshSetNames.size(); iSetName++ )
                {
                    // if the mesh set name not in map
                    if( tMeshtoFemSet.find( std::make_tuple(
                            tMeshSetNames( iSetName ),
                            tTimeContinuity,
                            tTimeBoundary ) ) == tMeshtoFemSet.end() )
                    {
                        // add the mesh set name map
                        tMeshtoFemSet[ std::make_tuple(
                                tMeshSetNames( iSetName ),
                                tTimeContinuity,
                                tTimeBoundary ) ] = tNumFEMSets++;

                        // create a fem set info for the mesh set
                        Set_User_Info aSetUserInfo;

                        // set its mesh set name
                        aSetUserInfo.set_mesh_set_name( tMeshSetNames( iSetName ) );

                        // set its time continuity flag
                        aSetUserInfo.set_time_continuity( tTimeContinuity );

                        // set its time boundary flag
                        aSetUserInfo.set_time_boundary( tTimeBoundary );

                        // set its sensitivity analysis type flag
                        aSetUserInfo.set_is_analytical_sensitivity_analysis( tIsAnalyticalSA );

                        // set its FD scheme for sensitivity analysis
                        aSetUserInfo.set_finite_difference_scheme_for_sensitivity_analysis( tFDSchemeForSA );

                        // set its FD perturbation size for sensitivity analysis
                        aSetUserInfo.set_finite_difference_perturbation_size( tFDPerturbation );

                        // set the IQI
                        aSetUserInfo.set_IQI( mIQIs( iIQI ) );

                        // add it to the list of fem set info
                        mSetInfo.push_back( aSetUserInfo );
                    }
                    else
                    {
                        // set the IQI
                        mSetInfo( tMeshtoFemSet[ std::make_tuple(
                                tMeshSetNames( iSetName ),
                                tTimeContinuity,
                                tTimeBoundary ) ] ).set_IQI( mIQIs( iIQI ) );
                    }
                }
            }
        }

        //------------------------------------------------------------------------------

        void FEM_Model::get_mesh_set_names(
                fem::Element_Type               aIWGBulkType,
                std::string                     aMasterPhaseName,
                std::string                     aSlavePhaseName,
                std::string                     aSlavePhaseString,
                std::string                     aOrdinalString,
                bool                            aIsGhost,
                moris::Cell< std::string >    & aMeshSetNames )
        {
            // get the master phase mesh index
            moris::Matrix< moris::IndexMat > tMasterPhaseIndices =
                    mPhaseInfo( mPhaseMap[ aMasterPhaseName ] ).get_phase_indices();

            // get the number of master phase mesh indices
            uint tNumMasterIndices = tMasterPhaseIndices.numel();

            // switch on the element type
            switch ( aIWGBulkType )
            {
                case fem::Element_Type::BULK :
                {
                    // loop over phase mesh indices
                    for( uint iMeshIndex = 0; iMeshIndex < tNumMasterIndices; iMeshIndex++ )
                    {
                        //                        // FIXME ! get mesh set names from integration mesh for index
                        //                        mMeshManager->get_integration_mesh( 0 )->
                        //                                get_block_set_names_with_color( tMasterPhaseIndices( iMeshIndex ), aMeshSetNames );

                        // add mesh set name to list
                        aMeshSetNames.push_back(
                                "HMR_dummy_c_p" +
                                std::to_string( tMasterPhaseIndices( iMeshIndex ) ) );

                        // add mesh set name to list
                        aMeshSetNames.push_back(
                                "HMR_dummy_n_p" +
                                std::to_string( tMasterPhaseIndices( iMeshIndex ) ) );
                    }
                    break;
                }
                case fem::Element_Type::SIDESET :
                {
                    // get neighbor phase names from string
                    moris::Cell< std::string > tSlavePhaseNames;
                    string_to_cell( aSlavePhaseString, tSlavePhaseNames );

                    // get number of neighbor phase
                    uint tNumSingle = tSlavePhaseNames.size();

                    // get ordinals for boundary from string
                    Matrix< DDSMat > tOrdinals;
                    string_to_mat( aOrdinalString, tOrdinals );
                    uint tNumBoundary = tOrdinals.numel();

                    // loop over master phase mesh indices
                    for( uint iMasterMeshIndex = 0; iMasterMeshIndex < tNumMasterIndices; iMasterMeshIndex++ )
                    {
                        // get single sideset
                        for( uint iSingle = 0; iSingle < tNumSingle; iSingle++ )
                        {
                            // get the neighbor phase name
                            std::string tNeighborPhaseName = tSlavePhaseNames( iSingle );

                            // get the slave phase mesh index
                            moris::Matrix< moris::IndexMat > tSlavePhaseIndices =
                                    mPhaseInfo( mPhaseMap[ tNeighborPhaseName ] ).get_phase_indices();

                            // get number of neighbor phase mesh indices
                            uint tNumNeighborIndices = tSlavePhaseIndices.numel();

                            for( uint iNeighborMeshIndex = 0; iNeighborMeshIndex < tNumNeighborIndices; iNeighborMeshIndex++ )
                            {
                                // FIXME get this info from the mesh
                                // add mesh set name to list
                                aMeshSetNames.push_back(
                                        "iside_b0_" +
                                        std::to_string( tMasterPhaseIndices( iMasterMeshIndex ) ) +
                                        "_b1_" +
                                        std::to_string( tSlavePhaseIndices( iNeighborMeshIndex ) ) );
                            }
                        }

                        // get boundary sideset
                        for( uint iBoundary = 0; iBoundary < tNumBoundary; iBoundary++ )
                        {
                            // FIXME get this info from the mesh
                            // add mesh set name to list
                            aMeshSetNames.push_back(
                                    "SideSet_" +
                                    std::to_string( tOrdinals( iBoundary ) ) +
                                    "_c_p" +
                                    std::to_string( tMasterPhaseIndices( iMasterMeshIndex ) ) );

                            // FIXME get this info from the mesh
                            // add mesh set name to list
                            aMeshSetNames.push_back(
                                    "SideSet_" +
                                    std::to_string( tOrdinals( iBoundary ) ) +
                                    "_n_p" +
                                    std::to_string( tMasterPhaseIndices( iMasterMeshIndex ) ) );
                        }
                    }
                    break;
                }
                case fem::Element_Type::DOUBLE_SIDESET :
                {
                    // if ghost
                    if ( aIsGhost )
                    {
                        // loop over master phase mesh indices
                        for( uint iMasterMeshIndex = 0; iMasterMeshIndex < tNumMasterIndices; iMasterMeshIndex++ )
                        {
                            // FIXME get this info from the mesh
                            // add mesh set name to list
                            aMeshSetNames.push_back(
                                    "ghost_p" +
                                    std::to_string( tMasterPhaseIndices( iMasterMeshIndex ) ) );
                        }
                    }
                    // if interface
                    else
                    {
                        MORIS_ERROR( aMasterPhaseName != aSlavePhaseName,
                                "FEM_Model::get_mesh_set_names - Master and slave phases are the same, FIXME case not handled yet ");

                        // get the slave phase mesh index
                        moris::Matrix< moris::IndexMat > tSlavePhaseIndices =
                                mPhaseInfo( mPhaseMap[ aSlavePhaseName ] ).get_phase_indices();

                        // get number of slave phase mesh index
                        uint tNumSlaveIndices = tSlavePhaseIndices.numel();

                        // loop over master phase mesh indices
                        for( uint iMasterMeshIndex = 0; iMasterMeshIndex < tNumMasterIndices; iMasterMeshIndex++ )
                        {
                            // get master index
                            uint tMasterPhaseIndex = tMasterPhaseIndices( iMasterMeshIndex );

                            // loop over slave phase mesh indices
                            for( uint iSlaveMeshIndex = 0; iSlaveMeshIndex < tNumSlaveIndices; iSlaveMeshIndex++ )
                            {
                                // get slave index
                                uint tSlavePhaseIndex = tSlavePhaseIndices( iSlaveMeshIndex );

                                // if master and slave index are different
                                if( tMasterPhaseIndex != tSlavePhaseIndex )
                                {
                                    // FIXME get this info from the mesh
                                    // get interface name
                                    aMeshSetNames.push_back(
                                            "dbl_iside_p0_" +
                                            std::to_string( tMasterPhaseIndices( iMasterMeshIndex ) ) +
                                            "_p1_" +
                                            std::to_string( tSlavePhaseIndices( iSlaveMeshIndex ) ) );
                                }
                            }
                        }
                    }
                    break;
                }
                default :
                {
                    MORIS_ERROR( false, "FEM_Model::get_mesh_set_names - Unknown set type" );
                }
            }
        }

        //-------------------------------------------------------------------------------------------------

        void FEM_Model::normalize_IQIs()
        {
            for (uint tRequestedIQIIndex = 0; tRequestedIQIIndex < mRequestedIQINames.size(); tRequestedIQIIndex++)
            {
                // IQI index
                uint tIQIIndex = 0;
                while (tIQIIndex < mIQIs.size()
                        and (mIQIs(tIQIIndex)->get_name() not_eq mRequestedIQINames(tRequestedIQIIndex)))
                {
                    tIQIIndex++;
                }
                MORIS_ASSERT(tIQIIndex < mIQIs.size(),
                        ("IQI was not found with the requested name " + mRequestedIQINames(tRequestedIQIIndex)).c_str());

                // Set normalization
                std::string tNormalization = mParameterList(4)(tIQIIndex).get< std::string >("normalization");
                if (tNormalization == "none")
                {
                    // Do nothing
                }
                else if (tNormalization == "time")
                {
                    MORIS_ERROR(false, "Time normalization not implemented yet for IQIs yet, implementation should go here.");
                }
                else if (tNormalization == "design")
                {
                    mIQIs(tRequestedIQIIndex)->set_reference_value(mGlobalIQIVal(tRequestedIQIIndex)(0));
                    mGlobalIQIVal(tRequestedIQIIndex)(0) = 1.0;
                }
                else
                {
                    // Try to set reference values directly
                    try
                    {
                        mIQIs(tRequestedIQIIndex)->set_reference_value(string_to_mat<DDRMat>(tNormalization)(0));
                    }
                    catch (...)
                    {
                        // create error message
                        std::string tErrMsg =
                                "FEM_Model::normalize_IQIs() - Unknown normalization: " + tNormalization +
                                ". Must be 'none', 'time', 'design', or a reference value.";

                        // error
                        MORIS_ERROR( false , tErrMsg.c_str() );
                    }
                }
            }
        }

        //-------------------------------------------------------------------------------------------------
        // FEM INPUT - old version
        //-------------------------------------------------------------------------------------------------

        void FEM_Model::create_constitutive_models(
                std::map< std::string, uint >            & aCMMap,
                std::map< std::string, uint >            & aPropertyMap,
                moris::map< std::string, MSI::Dof_Type > & aMSIDofTypeMap,
                moris::map< std::string, PDV_Type >      & aDvTypeMap )
        {
            // create a constitutive model factory
            CM_Factory tCMFactory;

            // get the CM parameter list
            moris::Cell< ParameterList > tCMParameterList = mParameterList( 1 );

            // get number of constitutive models
            uint tNumCMs = tCMParameterList.size();

            // create a list of CMs
            mCMs.resize( tNumCMs, nullptr );

            // loop over the parameter lists for CM
            for( uint iCM = 0; iCM < tNumCMs; iCM++ )
            {
                // get the constitutive type from parameter list
                fem::Constitutive_Type tCMType =
                        static_cast< fem::Constitutive_Type >( tCMParameterList( iCM ).get< uint >( "constitutive_type" ) );

                // create a constitutive model pointer
                mCMs( iCM ) = tCMFactory.create_CM( tCMType );

                // set CM name
                mCMs( iCM )->set_name( tCMParameterList( iCM ).get< std::string >( "constitutive_name" ) );

                // fill CM map
                aCMMap[ tCMParameterList( iCM ).get< std::string >( "constitutive_name" ) ] = iCM;

                // set CM space dimension
                mCMs( iCM )->set_space_dim( mSpaceDim );

                // set CM model type
                fem::Model_Type tCMModelType =
                        static_cast< fem::Model_Type >( tCMParameterList( iCM ).get< uint >( "model_type" ) );
                if( tCMModelType != fem::Model_Type::UNDEFINED )
                {
                    mCMs( iCM )->set_model_type( tCMModelType );
                }

                // set CM dof dependencies
                moris::Cell< moris::Cell< moris::MSI::Dof_Type > > tDofTypes;
                string_to_cell_of_cell(
                        std::get< 0 >( tCMParameterList( iCM ).get< std::pair< std::string, std::string > >( "dof_dependencies" ) ),
                        tDofTypes,
                        aMSIDofTypeMap );
                moris::Cell< std::string > tDofTypeNames;
                string_to_cell(
                        std::get< 1 >( tCMParameterList( iCM ).get< std::pair< std::string, std::string > >( "dof_dependencies" ) ),
                        tDofTypeNames );
                mCMs( iCM )->set_dof_type_list( tDofTypes, tDofTypeNames );

                // set CM dv dependencies
                moris::Cell< moris::Cell< PDV_Type > > tDvTypes;
                string_to_cell_of_cell(
                        std::get< 0 >( tCMParameterList( iCM ).get< std::pair< std::string, std::string > >( "dv_dependencies" ) ),
                        tDvTypes,
                        aDvTypeMap );
                moris::Cell< std::string > tDvTypeNames;
                string_to_cell(
                        std::get< 1 >( tCMParameterList( iCM ).get< std::pair< std::string, std::string > >( "dv_dependencies" ) ),
                        tDvTypeNames );
                mCMs( iCM )->set_dv_type_list( tDvTypes, tDvTypeNames );

                // set CM properties
                moris::Cell< moris::Cell< std::string > > tPropertyNamesPair;
                string_to_cell_of_cell(
                        tCMParameterList( iCM ).get< std::string >( "properties" ),
                        tPropertyNamesPair );
                for( uint iProp = 0; iProp < tPropertyNamesPair.size(); iProp++ )
                {
                    // if property name is in the property map
                    if ( aPropertyMap.find( tPropertyNamesPair( iProp )( 0 ) ) != aPropertyMap.end() )
                    {
                        // get property index
                        uint tPropertyIndex = aPropertyMap[ tPropertyNamesPair( iProp )( 0 ) ];

                        // set property for CM
                        mCMs( iCM )->set_property(
                                mProperties( tPropertyIndex ),
                                tPropertyNamesPair( iProp )( 1 ) );
                    }
                    else
                    {
                        // error message for unknown property
                        MORIS_ERROR( false,
                                "FEM_Model::create_CMs - Unknown aPropertyString : %s \n",
                                tPropertyNamesPair( iProp )( 0 ).c_str() );
                    }
                }
                // set local properties
                mCMs( iCM )->set_local_properties();

            }
        }

        //------------------------------------------------------------------------------

        void FEM_Model::create_stabilization_parameters(
                std::map< std::string, uint >            & aSPMap,
                std::map< std::string, uint >            & aPropertyMap,
                std::map< std::string, uint >            & aCMMap,
                moris::map< std::string, MSI::Dof_Type > & aMSIDofTypeMap,
                moris::map< std::string, PDV_Type >      & aDvTypeMap )
        {
            // create a stabilization parameter factory
            SP_Factory tSPFactory;

            // get the SP parameter list
            moris::Cell< ParameterList > tSPParameterList = mParameterList( 2 );

            // get the number of stabilization parameters
            uint tNumSPs = tSPParameterList.size();

            // set size for the list of stabilization parameter pointer
            mSPs.resize( tNumSPs, nullptr );

            // loop over the parameter list
            for ( uint iSP = 0; iSP < tNumSPs; iSP++ )
            {
                // get the SP parameter
                ParameterList tSPParameter = tSPParameterList( iSP );

                // get the stabilization type from parameter list
                fem::Stabilization_Type tSPType =
                        static_cast< fem::Stabilization_Type >( tSPParameter.get< uint >( "stabilization_type" ) );

                // create a stabilization parameter pointer
                mSPs( iSP ) = tSPFactory.create_SP( tSPType );

                // set name
                mSPs( iSP )->set_name( tSPParameter.get< std::string >( "stabilization_name" ) );

                // set SP space dimension
                mSPs( iSP )->set_space_dim( mSpaceDim );

                // fill stabilization map
                aSPMap[ tSPParameter.get< std::string >( "stabilization_name" ) ] = iSP;

                // set parameters
                moris::Cell< moris::Matrix< DDRMat > > tFuncParameters;
                string_to_cell_mat_2(
                        tSPParameter.get< std::string >( "function_parameters" ),
                        tFuncParameters );
                mSPs( iSP )->set_parameters( tFuncParameters );

                // init string for master or slave
                std::string tIsMasterString = "master";
                mtk::Master_Slave tIsMaster = mtk::Master_Slave::MASTER;

                // loop on master and slave
                for( uint iMaster = 0; iMaster <= mSPs( iSP )->get_has_slave(); iMaster++ )
                {
                    // if slave
                    if( iMaster )
                    {
                        // reset string for slave
                        tIsMasterString = "slave";
                        tIsMaster = mtk::Master_Slave::SLAVE;
                    }

                    // set dof dependencies
                    moris::Cell< moris::Cell< moris::MSI::Dof_Type > > tDofTypes;
                    string_to_cell_of_cell(
                            std::get< 0 >( tSPParameter.get< std::pair< std::string, std::string > >( tIsMasterString + "_dof_dependencies" ) ),
                            tDofTypes,
                            aMSIDofTypeMap );
                    moris::Cell< std::string > tDofTypeNames;
                    string_to_cell( std::get< 1 >(
                            tSPParameter.get< std::pair< std::string, std::string > >( tIsMasterString + "_dof_dependencies" ) ),
                            tDofTypeNames );
                    mSPs( iSP )->set_dof_type_list( tDofTypes, tDofTypeNames, tIsMaster );

                    // set dv dependencies
                    moris::Cell< moris::Cell< PDV_Type > > tDvTypes;
                    string_to_cell_of_cell(
                            std::get< 0 >( tSPParameter.get< std::pair< std::string, std::string > >( tIsMasterString + "_dv_dependencies" ) ),
                            tDvTypes,
                            aDvTypeMap );
                    moris::Cell< std::string > tDvTypeNames;
                    string_to_cell(
                            std::get< 1 >( tSPParameter.get< std::pair< std::string, std::string > >( tIsMasterString + "_dv_dependencies" ) ),
                            tDvTypeNames );
                    mSPs( iSP )->set_dv_type_list( tDvTypes, tDvTypeNames, tIsMaster );

                    // set master properties
                    moris::Cell< moris::Cell< std::string > > tPropertyNamesPair;
                    string_to_cell_of_cell(
                            tSPParameter.get< std::string >( tIsMasterString + "_properties" ),
                            tPropertyNamesPair );

                    for( uint iProp = 0; iProp < tPropertyNamesPair.size(); iProp++ )
                    {
                        // get the property name
                        std::string tPropertyName = tPropertyNamesPair( iProp )( 0 );

                        // check for unknown property
                        MORIS_ERROR( aPropertyMap.find( tPropertyName ) != aPropertyMap.end() ,
                                "FEM_Model::create_stabilization_parameters - Unknown master aPropertyString : %s \n",
                                tPropertyName.c_str() );

                        // get property index
                        uint tPropertyIndex = aPropertyMap[ tPropertyName ];

                        // set property for CM
                        mSPs( iSP )->set_property(
                                mProperties( tPropertyIndex ),
                                tPropertyNamesPair( iProp )( 1 ),
                                tIsMaster );
                    }

                    // set constitutive models
                    moris::Cell< moris::Cell< std::string > > tCMNamesPair;
                    string_to_cell_of_cell(
                            tSPParameter.get< std::string >( tIsMasterString + "_constitutive_models" ),
                            tCMNamesPair );

                    // loop over the CM names
                    for( uint iCM = 0; iCM < tCMNamesPair.size(); iCM++ )
                    {
                        // get the CM name
                        std::string tCMName = tCMNamesPair( iCM )( 0 );

                        // check for unknown CM
                        MORIS_ERROR( aCMMap.find( tCMName ) != aCMMap.end() ,
                                "FEM_Model::create_stabilization_parameters - Unknown master aCMString: %s \n",
                                tCMName.c_str() );

                        // get CM index
                        uint tCMIndex = aCMMap[ tCMName ];

                        // set CM for SP
                        mSPs( iSP )->set_constitutive_model(
                                mCMs( tCMIndex ),
                                tCMNamesPair( iCM )( 1 ) );
                    }
                }
            }
        }

        //------------------------------------------------------------------------------

        void FEM_Model::create_IWGs(
                std::map< std::string, uint >            & aIWGMap,
                std::map< std::string, uint >            & aPropertyMap,
                std::map< std::string, uint >            & aCMMap,
                std::map< std::string, uint >            & aSPMap,
                moris::map< std::string, MSI::Dof_Type > & aMSIDofTypeMap,
                moris::map< std::string, PDV_Type >      & aDvTypeMap )
        {
            // create an IWG factory
            IWG_Factory tIWGFactory;

            // get the IWG parameter list
            moris::Cell< ParameterList > tIWGParameterList = mParameterList( 3 );

            // get number of IWGs
            uint tNumIWGs = tIWGParameterList.size();

            // create a list of IWG pointers
            mIWGs.resize( tNumIWGs, nullptr );

            // loop over the parameter lists
            for( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
            {
                // get the treated IWG parameter list
                ParameterList tIWGParameter = tIWGParameterList( iIWG );

                // get the IWG type from parameter list
                fem::IWG_Type tIWGType =
                        static_cast< fem::IWG_Type >( tIWGParameter.get< uint >( "IWG_type" ) );

                // create an IWG pointer
                mIWGs( iIWG ) = tIWGFactory.create_IWG( tIWGType );

                // get the IWG name
                std::string tIWGName = tIWGParameter.get< std::string >( "IWG_name" );

                // set name
                mIWGs( iIWG )->set_name( tIWGName );

                // fill IWG map
                aIWGMap[ tIWGName ] = iIWG;

                // get the ghost order from parameter list
                uint tGhostOrder = tIWGParameter.get< uint >( "ghost_order" );
                mIWGs( iIWG )->set_interpolation_order( tGhostOrder );

                // set residual dof type
                moris::Cell< moris::MSI::Dof_Type > tResDofTypes;
                string_to_cell(
                        tIWGParameter.get< std::string >( "dof_residual" ),
                        tResDofTypes,
                        aMSIDofTypeMap );
                mIWGs( iIWG )->set_residual_dof_type( tResDofTypes );

                // init string for master or slave
                std::string tIsMasterString = "master";
                mtk::Master_Slave tIsMaster = mtk::Master_Slave::MASTER;

                // loop on master and slave
                for( uint iMaster = 0; iMaster <= 1; iMaster++ )
                {
                    // if slave
                    if( iMaster )
                    {
                        // reset string for slave
                        tIsMasterString = "slave";
                        tIsMaster = mtk::Master_Slave::SLAVE;
                    }

                    // set dof dependencies
                    moris::Cell< moris::Cell< moris::MSI::Dof_Type > > tDofTypes;
                    string_to_cell_of_cell(
                            tIWGParameter.get< std::string >( tIsMasterString + "_dof_dependencies" ),
                            tDofTypes,
                            aMSIDofTypeMap );
                    mIWGs( iIWG )->set_dof_type_list( tDofTypes, tIsMaster );

                    // set dv dependencies
                    moris::Cell< moris::Cell< PDV_Type > > tDvTypes;
                    string_to_cell_of_cell(
                            tIWGParameter.get< std::string >( tIsMasterString + "_dv_dependencies" ),
                            tDvTypes,
                            aDvTypeMap );
                    mIWGs( iIWG )->set_dv_type_list( tDvTypes, tIsMaster );

                    // set properties
                    moris::Cell< moris::Cell< std::string > > tPropertyNamesPair;
                    string_to_cell_of_cell(
                            tIWGParameter.get< std::string >( tIsMasterString + "_properties" ),
                            tPropertyNamesPair );

                    for( uint iProp = 0; iProp < tPropertyNamesPair.size(); iProp++ )
                    {
                        // if property name is in the property map
                        if ( aPropertyMap.find( tPropertyNamesPair( iProp )( 0 ) ) != aPropertyMap.end() )
                        {
                            // get property index
                            uint tPropertyIndex = aPropertyMap[ tPropertyNamesPair( iProp )( 0 ) ];

                            // set property for IWG
                            mIWGs( iIWG )->set_property(
                                    mProperties( tPropertyIndex ),
                                    tPropertyNamesPair( iProp )( 1 ),
                                    tIsMaster );
                        }
                        else
                        {
                            // create error message unknown property
                            MORIS_ERROR( false ,
                                    "FEM_Model::create_IWGs - Unknown %s aPropertyString: %s \n",
                                    tIsMasterString.c_str(),
                                    tPropertyNamesPair( iProp )( 0 ).c_str() );
                        }
                    }

                    // set constitutive models
                    moris::Cell< moris::Cell< std::string > > tCMNamesPair;
                    string_to_cell_of_cell(
                            tIWGParameter.get< std::string >( tIsMasterString + "_constitutive_models" ),
                            tCMNamesPair );

                    for( uint iCM = 0; iCM < tCMNamesPair.size(); iCM++ )
                    {
                        // if CM name is in the CM map
                        if ( aCMMap.find( tCMNamesPair( iCM )( 0 ) ) != aCMMap.end() )
                        {
                            // get CM index
                            uint tCMIndex = aCMMap[ tCMNamesPair( iCM )( 0 ) ];

                            // set CM for IWG
                            mIWGs( iIWG )->set_constitutive_model(
                                    mCMs( tCMIndex ),
                                    tCMNamesPair( iCM )( 1 ),
                                    tIsMaster );
                        }
                        else
                        {
                            // error message unknown CM
                            MORIS_ERROR( false ,
                                    "FEM_Model::create_IWGs - Unknown %s aCMString: %s \n",
                                    tIsMasterString.c_str(),
                                    tCMNamesPair( iCM )( 0 ).c_str() );
                        }
                    }
                }

                // set stabilization parameters
                moris::Cell< moris::Cell< std::string > > tSPNamesPair;
                string_to_cell_of_cell(
                        tIWGParameter.get< std::string >( "stabilization_parameters" ),
                        tSPNamesPair );

                for( uint iSP = 0; iSP < tSPNamesPair.size(); iSP++ )
                {
                    // if CM name is in the CM map
                    if ( aSPMap.find( tSPNamesPair( iSP )( 0 ) ) != aSPMap.end() )
                    {
                        // get SP index
                        uint tSPIndex = aSPMap[ tSPNamesPair( iSP )( 0 ) ];

                        // set SP for IWG
                        mIWGs( iIWG )->set_stabilization_parameter(
                                mSPs( tSPIndex ),
                                tSPNamesPair( iSP )( 1 ) );
                    }
                    else
                    {
                        // error message unknown SP
                        MORIS_ERROR( false ,
                                "FEM_Model::create_IWGs - Unknown aSPString: %s \n",
                                tSPNamesPair( iSP )( 0 ).c_str() );
                    }
                }

                //                // debug
                //                mIWGs( iIWG )->print_names();
            }
        }

        //------------------------------------------------------------------------------

        void FEM_Model::create_IQIs(
                std::map< std::string, uint >            & aIQIMap,
                std::map< std::string, uint >            & aPropertyMap,
                std::map< std::string, uint >            & aCMMap,
                std::map< std::string, uint >            & aSPMap,
                moris::map< std::string, MSI::Dof_Type > & aMSIDofTypeMap,
                moris::map< std::string, PDV_Type >      & aDvTypeMap )
        {
            // create an IQI factory
            IQI_Factory tIQIFactory;

            // get the IQI parameter list
            moris::Cell< ParameterList > tIQIParameterList = mParameterList( 4 );

            // get number of IQIs
            uint tNumIQIs = tIQIParameterList.size();

            // set size for list of IQI pointers
            mIQIs.resize( tNumIQIs, nullptr );

            // loop over the parameter lists
            for( uint iIQI = 0; iIQI < tNumIQIs; iIQI++ )
            {
                // get the treated IQI parameter list
                ParameterList tIQIParameter = tIQIParameterList( iIQI );

                // get name from parameter list
                std::string tIQIName = tIQIParameter.get< std::string >( "IQI_name" );

                // get the IQI type from parameter list
                fem::IQI_Type tIQIType =
                        static_cast< fem::IQI_Type >( tIQIParameter.get< uint >( "IQI_type" ) );

                // create an IQI pointer
                mIQIs( iIQI ) = tIQIFactory.create_IQI( tIQIType );

                // set name
                mIQIs( iIQI )->set_name( tIQIName );

                // fill IQI map
                aIQIMap[ tIQIName ] = iIQI;

                // get the treated IQI quantity dof type
                moris::Cell< moris::MSI::Dof_Type > tQuantityDofTypes;
                string_to_cell(
                        tIQIParameter.get< std::string >( "dof_quantity" ),
                        tQuantityDofTypes,
                        aMSIDofTypeMap );
                mIQIs( iIQI )->set_quantity_dof_type( tQuantityDofTypes );

                // set index for vectorial field
                mIQIs( iIQI )->set_output_type_index(
                        tIQIParameter.get< moris::sint >( "vectorial_field_index" ) );

                // set function parameters
                moris::Cell< moris::Matrix< DDRMat > > tFuncParameters;
                string_to_cell_mat_2(
                        tIQIParameter.get< std::string >( "function_parameters" ),
                        tFuncParameters );
                mIQIs( iIQI )->set_parameters( tFuncParameters );

                // init string for master or slave
                std::string tIsMasterString = "master";
                mtk::Master_Slave tIsMaster = mtk::Master_Slave::MASTER;

                // loop on master and slave
                for( uint iMaster = 0; iMaster <= 1; iMaster++ )
                {
                    // if slave
                    if( iMaster )
                    {
                        // reset string for slave
                        tIsMasterString = "slave";
                        tIsMaster = mtk::Master_Slave::SLAVE;
                    }

                    // set dof dependencies
                    moris::Cell< moris::Cell< moris::MSI::Dof_Type > > tDofTypes;
                    string_to_cell_of_cell(
                            tIQIParameter.get< std::string >( tIsMasterString + "_dof_dependencies" ),
                            tDofTypes,
                            aMSIDofTypeMap );
                    mIQIs( iIQI )->set_dof_type_list( tDofTypes, tIsMaster );

                    // set dv dependencies
                    moris::Cell< moris::Cell< PDV_Type > > tDvTypes;
                    string_to_cell_of_cell(
                            tIQIParameter.get< std::string >( tIsMasterString + "_dv_dependencies" ),
                            tDvTypes,
                            aDvTypeMap );
                    mIQIs( iIQI )->set_dv_type_list( tDvTypes, tIsMaster );

                    // set properties
                    moris::Cell< moris::Cell< std::string > > tPropertyNamesPair;
                    string_to_cell_of_cell(
                            tIQIParameter.get< std::string >( tIsMasterString + "_properties" ),
                            tPropertyNamesPair );

                    for( uint iProp = 0; iProp < tPropertyNamesPair.size(); iProp++ )
                    {
                        // if property name is in the property map
                        if ( aPropertyMap.find( tPropertyNamesPair( iProp )( 0 ) ) != aPropertyMap.end() )
                        {
                            // get property index
                            uint tPropertyIndex = aPropertyMap[ tPropertyNamesPair( iProp )( 0 ) ];

                            // set property for IWG
                            mIQIs( iIQI )->set_property(
                                    mProperties( tPropertyIndex ),
                                    tPropertyNamesPair( iProp )( 1 ),
                                    tIsMaster );
                        }
                        else
                        {
                            // error message unknown property
                            MORIS_ERROR( false ,
                                    "FEM_Model::create_IQIs - Unknown %s aPropertyString: %s \n",
                                    tIsMasterString.c_str(),
                                    tPropertyNamesPair( iProp )( 0 ).c_str() );
                        }
                    }

                    // set constitutive models
                    moris::Cell< moris::Cell< std::string > > tCMNamesPair;
                    string_to_cell_of_cell(
                            tIQIParameter.get< std::string >( tIsMasterString + "_constitutive_models" ),
                            tCMNamesPair );

                    for( uint iCM = 0; iCM < tCMNamesPair.size(); iCM++ )
                    {
                        // if CM name is in the CM map
                        if ( aCMMap.find( tCMNamesPair( iCM )( 0 ) ) != aCMMap.end() )
                        {
                            // get CM index
                            uint tCMIndex = aCMMap[ tCMNamesPair( iCM )( 0 ) ];

                            // set CM for IQI
                            mIQIs( iIQI )->set_constitutive_model(
                                    mCMs( tCMIndex ),
                                    tCMNamesPair( iCM )( 1 ),
                                    tIsMaster );
                        }
                        else
                        {
                            // error message unknown CM
                            MORIS_ERROR( false ,
                                    "FEM_Model::create_IQIs - Unknown %s aCMString: %s \n",
                                    tIsMasterString.c_str(),
                                    tCMNamesPair( iCM )( 0 ).c_str() );
                        }
                    }
                }

                // set stabilization parameters
                moris::Cell< moris::Cell< std::string > > tSPNamesPair;
                string_to_cell_of_cell(
                        tIQIParameter.get< std::string >( "stabilization_parameters" ),
                        tSPNamesPair );

                for( uint iSP = 0; iSP < tSPNamesPair.size(); iSP++ )
                {
                    // if SP name is in the SP map
                    if ( aSPMap.find( tSPNamesPair( iSP )( 0 ) ) != aSPMap.end() )
                    {
                        // get SP index
                        uint tSPIndex = aSPMap[ tSPNamesPair( iSP )( 0 ) ];

                        // set SP for IQI
                        mIQIs( iIQI )->set_stabilization_parameter(
                                mSPs( tSPIndex ),
                                tSPNamesPair( iSP )( 1 ) );
                    }
                    else
                    {
                        // error message unknown SP
                        MORIS_ERROR( false ,
                                "FEM_Model::create_IQIs - Unknown aSPString: %s \n",
                                tSPNamesPair( iSP )( 0 ).c_str() );
                    }
                }

                //                // debug
                //                mIQIs( iIQI )->print_names();
            }
        }

        //------------------------------------------------------------------------------

        void FEM_Model::create_fem_set_info()
        {
            // init number of fem sets to be created
            uint tNumFEMSets = 0;

            // get the IWG and IQI parameter lists
            moris::Cell< ParameterList > tIWGParameterList = mParameterList( 3 );
            moris::Cell< ParameterList > tIQIParameterList = mParameterList( 4 );

            // get fem computation type parameter list
            ParameterList tComputationParameterList = mParameterList( 5 )( 0 );

            // get bool for printing physics model
            bool tPrintPhysics =
                    tComputationParameterList.get< bool >( "print_physics_model" );

            // get bool for analytical/finite differenec for SA
            bool tIsAnalyticalSA =
                    tComputationParameterList.get< bool >( "is_analytical_sensitivity" );

            // get enum for FD scheme
            fem::FDScheme_Type tFDSchemeForSA = static_cast< fem::FDScheme_Type >(
                    tComputationParameterList.get< uint >( "finite_difference_scheme" ) );

            // get perturbation size for FD
            real tFDPerturbation = tComputationParameterList.get< real >(
                    "finite_difference_perturbation_size" );

            // create a map of the set
            std::map< std::tuple< std::string, bool, bool >, uint > tMeshtoFemSet;

            // loop over the IWGs
            for( uint iIWG = 0; iIWG < tIWGParameterList.size(); iIWG++ )
            {
                // get the mesh set names from the IWG parameter list
                moris::Cell< std::string > tMeshSetNames;
                string_to_cell( tIWGParameterList( iIWG ).get< std::string >( "mesh_set_names" ), tMeshSetNames );

                // get the time continuity flag from the IWG parameter list
                bool tTimeContinuity = tIWGParameterList( iIWG ).get< bool >( "time_continuity" );

                // get the time boundary flag from the IQI parameter list
                bool tTimeBoundary = tIWGParameterList( iIWG ).get< bool >( "time_boundary" );

                // loop over the mesh set names
                for( uint iSetName = 0; iSetName < tMeshSetNames.size(); iSetName++ )
                {
                    // check if the mesh set name already in map
                    if( tMeshtoFemSet.find( std::make_tuple(
                            tMeshSetNames( iSetName ),
                            tTimeContinuity,
                            tTimeBoundary ) ) == tMeshtoFemSet.end() )
                    {
                        // add the mesh set name map
                        tMeshtoFemSet[ std::make_tuple(
                                tMeshSetNames( iSetName ),
                                tTimeContinuity,
                                tTimeBoundary ) ] = tNumFEMSets++;

                        // create a fem set info for the mesh set
                        Set_User_Info aSetUserInfo;

                        // set its mesh set name
                        aSetUserInfo.set_mesh_set_name( tMeshSetNames( iSetName ) );

                        // set its time continuity flag
                        aSetUserInfo.set_time_continuity( tTimeContinuity );

                        // set its time boundary flag
                        aSetUserInfo.set_time_boundary( tTimeBoundary );

                        // set its sensitivity analysis type flag
                        aSetUserInfo.set_is_analytical_sensitivity_analysis( tIsAnalyticalSA );

                        // set its FD scheme for sensitivity analysis
                        aSetUserInfo.set_finite_difference_scheme_for_sensitivity_analysis( tFDSchemeForSA );

                        // set its FD perturbation size for sensitivity analysis
                        aSetUserInfo.set_finite_difference_perturbation_size( tFDPerturbation );

                        // set the IWG
                        aSetUserInfo.set_IWG( mIWGs( iIWG ) );

                        // add it to the list of fem set info
                        mSetInfo.push_back( aSetUserInfo );
                    }
                    else
                    {
                        // set the IWG
                        mSetInfo( tMeshtoFemSet[ std::make_tuple(
                                tMeshSetNames( iSetName ),
                                tTimeContinuity,
                                tTimeBoundary ) ] ).set_IWG( mIWGs( iIWG ) );
                    }
                }
            }

            // loop over the IQIs
            for( uint iIQI = 0; iIQI < tIQIParameterList.size(); iIQI++ )
            {
                // get the mesh set names from the IQI parameter list
                moris::Cell< std::string > tMeshSetNames;
                string_to_cell( tIQIParameterList( iIQI ).get< std::string >( "mesh_set_names" ), tMeshSetNames );

                // get the time continuity flag from the IQI parameter list
                bool tTimeContinuity = tIQIParameterList( iIQI ).get< bool >( "time_continuity" );

                // get the time boundary flag from the IQI parameter list
                bool tTimeBoundary = tIQIParameterList( iIQI ).get< bool >( "time_boundary" );

                // loop over the mesh set names
                for( uint iSetName = 0; iSetName < tMeshSetNames.size(); iSetName++ )
                {
                    // if the mesh set name not in map
                    if( tMeshtoFemSet.find( std::make_tuple(
                            tMeshSetNames( iSetName ),
                            tTimeContinuity,
                            tTimeBoundary ) ) == tMeshtoFemSet.end() )
                    {
                        // add the mesh set name map
                        tMeshtoFemSet[ std::make_tuple(
                                tMeshSetNames( iSetName ),
                                tTimeContinuity,
                                tTimeBoundary ) ] = tNumFEMSets++;

                        // create a fem set info for the mesh set
                        Set_User_Info aSetUserInfo;

                        // set its mesh set name
                        aSetUserInfo.set_mesh_set_name( tMeshSetNames( iSetName ) );

                        // set its time continuity flag
                        aSetUserInfo.set_time_continuity( tTimeContinuity );

                        // set its time boundary flag
                        aSetUserInfo.set_time_boundary( tTimeBoundary );

                        // set its sensitivity analysis type flag
                        aSetUserInfo.set_is_analytical_sensitivity_analysis( tIsAnalyticalSA );

                        // set its FD scheme for sensitivity analysis
                        aSetUserInfo.set_finite_difference_scheme_for_sensitivity_analysis( tFDSchemeForSA );

                        // set its FD perturbation size for sensitivity analysis
                        aSetUserInfo.set_finite_difference_perturbation_size( tFDPerturbation );

                        // set the IQI
                        aSetUserInfo.set_IQI( mIQIs( iIQI ) );

                        // add it to the list of fem set info
                        mSetInfo.push_back( aSetUserInfo );
                    }
                    else
                    {
                        // set the IQI
                        mSetInfo( tMeshtoFemSet[ std::make_tuple(
                                tMeshSetNames( iSetName ),
                                tTimeContinuity,
                                tTimeBoundary ) ] ).set_IQI( mIQIs( iIQI ) );
                    }
                }
            }

            // debug print
            if( tPrintPhysics )
            {
                for( uint iSet = 0; iSet < mSetInfo.size(); iSet++ )
                {
                    std::cout<<"%-------------------------------------------------"<<std::endl;
                    mSetInfo( iSet ).print_names();
                    std::cout<<"%-------------------------------------------------"<<std::endl;
                }
            }
        }
        //-------------------------------------------------------------------------------------------------

    } /* namespace mdl */
} /* namespace moris */
