/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Model.cpp
 *
 */


#include "cl_MTK_Integration_Rule.hpp"
#include "cl_MTK_Side_Set.hpp"
#include "cl_MTK_Integrator.hpp"
#ifdef WITHGPERFTOOLS
#include <gperftools/profiler.h>
#endif

#include <map>
#include <set>
#include <algorithm>

// LINALG/src
#include "cl_Map.hpp"
#include "cl_Matrix.hpp"
#include "fn_unique.hpp"
#include "fn_sum.hpp"    // for check
#include "fn_iscol.hpp"
#include "fn_trans.hpp"
#include "moris_typedefs.hpp"
#include "op_equal_equal.hpp"
// MTK/src
#include "MTK_Tools.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
// FEM/INT/src
#include "cl_FEM_Node_Base.hpp"
#include "cl_FEM_Node.hpp"
#include "cl_FEM_Enums.hpp"
#include "cl_FEM_Model.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_MM_Factory.hpp"
#include "cl_FEM_CM_Factory.hpp"
#include "cl_FEM_SP_Factory.hpp"
#include "cl_FEM_IWG_Factory.hpp"
#include "cl_FEM_IQI_Factory.hpp"
#include "cl_FEM_Field.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
// FEM/MSI/src
#include "cl_MSI_Equation_Object.hpp"
#include "cl_MSI_Dof_Type_Enums.hpp"
// FEM/VIS/src
#include "cl_VIS_Output_Enums.hpp"
// GEN/src
#include "GEN_Data_Types.hpp"
// Logging package
#include "cl_Logger.hpp"
#include "cl_Tracer.hpp"
#include "cl_FEM_Model_Initializer.hpp"
#include "cl_FEM_Model_Initializer_Legacy.hpp"
#include "cl_FEM_Model_Initializer_Phasebased.hpp"
#include "cl_MTK_Mesh_DataBase_IP.hpp"
#include "cl_MTK_Mesh_DataBase_IG.hpp"


#include <cl_MTK_Json_Debug_Output.hpp>

namespace moris
{
    namespace fem
    {
        FEM_Model::FEM_Model(
                std::shared_ptr< mtk::Mesh_Manager > aMeshManager,
                const moris_index                   &aMeshPairIndex,
                moris::Vector< fem::Set_User_Info > &aSetInfo )
                : mMeshManager( aMeshManager )
                , mMeshPairIndex( aMeshPairIndex )
        {
            Tracer tTracer( "FEM", "Model", "Create" );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 0: unpack mesh
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // get pointers to interpolation and integration meshes
            mtk::Interpolation_Mesh *tIPMesh = nullptr;
            mtk::Integration_Mesh   *tIGMesh = nullptr;

            mMeshManager->get_mesh_pair( mMeshPairIndex, tIPMesh, tIGMesh );

            // set the space dimension
            mSpaceDim = tIPMesh->get_spatial_dim();

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 1: create nodes
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            this->create_interpolation_nodes( tIPMesh );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 2: create sets
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            this->create_fem_sets( tIPMesh, tIGMesh, aSetInfo );
        }


        FEM_Model::FEM_Model(
                std::shared_ptr< mtk::Mesh_Manager > aMeshManager,
                const moris_index                   &aMeshPairIndex,
                moris::Vector< fem::Set_User_Info > &aSetInfo,
                MSI::Design_Variable_Interface      *aDesignVariableInterface )
                : mMeshManager( aMeshManager )
                , mMeshPairIndex( aMeshPairIndex )
        {
            Tracer tTracer( "FEM", "Model", "Create" );

            this->set_design_variable_interface( aDesignVariableInterface );

            // if no design variables have been stipulated, skip
            if ( aDesignVariableInterface == nullptr )
            {
                mFEMOnly = true;
                MORIS_LOG( "Skipping GEN, FEM Only" );
            }

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 0: unpack mesh
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // get pointers to interpolation and integration meshes
            mtk::Interpolation_Mesh *tIPMesh = nullptr;
            mtk::Integration_Mesh   *tIGMesh = nullptr;

            mMeshManager->get_mesh_pair( mMeshPairIndex, tIPMesh, tIGMesh );

            // set the space dimension
            mSpaceDim = tIPMesh->get_spatial_dim();

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 1: create interpolation nodes
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            this->create_interpolation_nodes( tIPMesh );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 2: create integration nodes
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            this->create_integration_nodes( tIGMesh );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 3: create fem sets
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            this->create_fem_sets( tIPMesh, tIGMesh, aSetInfo );
        }


        FEM_Model::FEM_Model(
                std::shared_ptr< mtk::Mesh_Manager >     aMeshManager,
                const moris_index                       &aMeshPairIndex,
                const Vector< Vector< ParameterList > > &aParameterList,
                std::shared_ptr< Library_IO >            aLibrary )
                : mMeshManager( aMeshManager )
                , mMeshPairIndex( aMeshPairIndex )
                , mParameterList( aParameterList )
        {
            Tracer tTracer( "FEM", "Model", "Create" );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 0: unpack fem input and mesh
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // get pointers to interpolation and integration meshes
            mtk::Interpolation_Mesh *tIPMesh = mMeshManager->get_interpolation_mesh( mMeshPairIndex );
            mtk::Integration_Mesh   *tIGMesh = mMeshManager->get_integration_mesh( mMeshPairIndex );


            // set the space dimension
            mSpaceDim = tIPMesh->get_spatial_dim();

            // unpack the FEM inputs
            this->initialize( aLibrary );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 1: create IP nodes
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            this->create_interpolation_nodes( tIPMesh );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 2: create fem sets
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            this->create_fem_sets( tIPMesh, tIGMesh );
        }


        FEM_Model::FEM_Model(
                std::shared_ptr< mtk::Mesh_Manager > aMeshManager,
                const moris_index                   &aMeshPairIndex,
                Vector< Vector< ParameterList > >    aParameterList,
                MSI::Design_Variable_Interface      *aDesignVariableInterface )
                : mMeshManager( aMeshManager )
                , mMeshPairIndex( aMeshPairIndex )
                , mParameterList( aParameterList )
        {
            Tracer tTracer( "FEM", "Model", "Create" );

            this->set_design_variable_interface( aDesignVariableInterface );

            // if no design variables have been stipulated, skip
            if ( aDesignVariableInterface == nullptr )
            {
                mFEMOnly = true;
                MORIS_LOG( "Skipping GEN, FEM Only" );
            }
        }


        void
        FEM_Model::create_interpolation_nodes( mtk::Interpolation_Mesh *aIPMesh )
        {
            // ask mesh about number of IP nodes on proc
            luint tNumIPNodes = aIPMesh->get_num_nodes();

            // create IP node objects
            mIPNodes.resize( tNumIPNodes, nullptr );

            for ( uint iNode = 0; iNode < tNumIPNodes; iNode++ )
            {
                mIPNodes( iNode ) = new fem::Node( &aIPMesh->get_mtk_vertex( iNode ) );
            }

            // print output
            MORIS_LOG_SPEC( "Number of IP Nodes", sum_all( tNumIPNodes ) );
        }


        void
        FEM_Model::create_integration_nodes( mtk::Integration_Mesh *aIGMesh )
        {
            // ask IG mesh about number of IG vertices on proc
            luint tNumIGNodes = aIGMesh->get_num_nodes();

            // set size for list IG nodes
            mIGNodes.resize( tNumIGNodes, nullptr );

            Vector< enum gen::PDV_Type > tGeoPdvType;
            switch ( mSpaceDim )
            {
                case 1:
                    tGeoPdvType = { gen::PDV_Type::X_COORDINATE };
                    break;
                case 2:
                    tGeoPdvType = { gen::PDV_Type::X_COORDINATE, gen::PDV_Type::Y_COORDINATE };
                    break;
                case 3:
                    tGeoPdvType = { gen::PDV_Type::X_COORDINATE, gen::PDV_Type::Y_COORDINATE, gen::PDV_Type::Z_COORDINATE };
                    break;
                default:
                    MORIS_ERROR( false,
                            "FEM_Model::create_integration_nodes - Space dimension can only be 1, 2 or 3D." );
            }

            // set the size of mIsActiveXYZ, mXYZPdvIds and mXYZLocalAssemblyIndices
            mIsActiveXYZ.set_size( tNumIGNodes, tGeoPdvType.size() );
            mXYZPdvIds.set_size( tNumIGNodes, tGeoPdvType.size() );

            // TODO: fill with default value causes nothing to crash, but it does without, why?
            mXYZLocalAssemblyIndices.set_size( tNumIGNodes, tGeoPdvType.size(), MORIS_INDEX_MAX );
            // mXYZLocalAssemblyIndices.set_size( tNumIGNodes, tGeoPdvType.size() );

            // TODO: MESHCLEANUP
            // load Mesh-GEN map
            // Vector<moris_index> tMesh_GEN_map;
            // tMesh_GEN_map.reserve(660);
            // aIGMesh->get_Mesh_GEN_map( tMesh_GEN_map );

            // loop over IG mesh vertices
            for ( uint iNode = 0; iNode < tNumIGNodes; iNode++ )
            {
                // create a new IG Node
                mIGNodes( iNode ) = new fem::Node( &aIGMesh->get_mtk_vertex( iNode ) );

                if ( mFEMOnly == false )
                {
                    // get IG node index
                    uint tVertexIndex = mIGNodes( iNode )->get_index();

                    // get the pdv values from the MSI/GEN interface
                    Matrix< IndexMat >         tVertexIndices( 1, 1, tVertexIndex );
                    Vector< Matrix< DDRMat > > tVertexCoordsFromGen( mSpaceDim );
                    Vector< Vector< bool > >   tIsActiveDv;

                    this->get_design_variable_interface()->get_ig_pdv_value(
                            tVertexIndices,
                            tGeoPdvType,
                            tVertexCoordsFromGen,
                            tIsActiveDv );

                    // set active flags for xyz
                    this->set_vertex_xyz_active_flags( tVertexIndex, tIsActiveDv );

                    // reshape the XYZ values into a cell of vectors
                    for ( uint iSpaceDim = 0; iSpaceDim < mSpaceDim; iSpaceDim++ )
                    {
                        if ( tIsActiveDv( iSpaceDim )( 0 ) )
                        {
                            // get IG node coordinates
                            Matrix< DDRMat > tVertexCoordsFromMesh;
                            mIGNodes( iNode )->get_vertex_coords( tVertexCoordsFromMesh );

                            MORIS_ERROR( equal_to(
                                                 tVertexCoordsFromGen( iSpaceDim )( 0 ),
                                                 tVertexCoordsFromMesh( iSpaceDim ),
                                                 1.0 ),
                                    "FEM_Model::create_integration_nodes - GE coordinate and MTK coordinate differ\n" );
                        }
                    }

                    // get the id associated to the pdv
                    Vector< moris::Matrix< DDSMat > > tPdvIds;
                    this->get_design_variable_interface()->get_ig_dv_ids_for_type_and_ind(
                            tVertexIndices,
                            tGeoPdvType,
                            tPdvIds );

                    // set pdv ids for xyz
                    this->set_vertex_xyz_pdv_ids( tVertexIndex, tPdvIds );
                }
            }

            // print output
            MORIS_LOG_SPEC( "Number of IG Nodes", sum_all( tNumIGNodes ) );
        }


        void
        FEM_Model::create_fem_sets(
                mtk::Interpolation_Mesh      *aIPMesh,
                mtk::Integration_Mesh        *aIGMesh,
                Vector< fem::Set_User_Info > &aSetInfo )
        {
            // get the number of sets
            uint tNumFemSets = aSetInfo.size();

            // create equation sets
            mFemSets.resize( tNumFemSets, nullptr );

            // get number of IP cells
            uint tNumIPCells = aIPMesh->get_num_elems();

            // reserve size for list of equation objects
            mFemClusters.reserve( tNumIPCells );

            // loop over the used fem set
            for ( luint iSet = 0; iSet < tNumFemSets; iSet++ )
            {
                // get the mesh set name
                std::string tMeshSetName = aSetInfo( iSet ).get_mesh_set_name();

                moris_index tMeshSetIndex;
                if ( tMeshSetName.size() > 0 )
                {
                    // get the mesh set index from its name
                    tMeshSetIndex = aIGMesh->get_set_index_by_name( tMeshSetName );
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
                moris::mtk::Set *tMeshSet = aIGMesh->get_set_by_index( tMeshSetIndex );

                // if non-empty mesh set
                if ( tMeshSet->get_num_clusters_on_set() != 0 )
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

            uint tNumElements       = mFemClusters.size();
            uint tGlobalNumElements = sum_all( tNumElements );

            // print output
            MORIS_LOG_SPEC( "IP elements", tGlobalNumElements );
        }

        void FEM_Model::update_equation_sets()
        {
            // early return if no nonconformal set is present
            if ( std::none_of( mFemSets.begin(), mFemSets.end(), []( auto const &tSet ) { return tSet->get_element_type() == fem::Element_Type::NONCONFORMAL_SIDESET; } ) )
            {
                return;
            }
            Tracer tTracer( "FEM", "Model", "Remapping" );


            // get the side set names that get used by the contact mesh editor to build nonconformal sets
            std::set< moris_index > tRequestedIGNodes;
            std::set< moris_index > tRequestedIPNodes;

            // TODO @ff: Remove! Debug only. This makes sure that every vertex displacement will be requested, not only the ones from the nonconformal sets.
//            for ( auto const &tFemSet : mFemSets )
//            {
//                if ( tFemSet->is_empty_set() || tFemSet->get_set_name().find( "ghost" ) != std::string::npos )
//                {
//                    continue;
//                }
//                auto const &tMeshSet = dynamic_cast< fem::Set *const >( tFemSet )->get_mesh_set();
//                for ( auto const &tCluster : tMeshSet->get_clusters_on_set() )
//                {
//                    for ( auto const &tCell : tCluster->get_primary_cells_in_cluster() )
//                    {
//                        for ( auto const &tVertex : tCell->get_vertex_pointers() )
//                        {
//                            tRequestedIGNodes.insert( tVertex->get_index() );
//                        }
//                    }
//                    for ( auto const &tVertex : tCluster->get_interpolation_cell().get_vertex_pointers() )
//                    {
//                        tRequestedIPNodes.insert( tVertex->get_index() );
//                    }
//                }
//            }

            std::map< moris_index, Vector< real > > tNodalDisplacements;
            for ( auto const &tSet : mFemSets )
            {
                if ( !tSet->is_empty_set() && !tRequestedIGNodes.empty() && tSet->get_set_name().find( "ghost" ) == std::string::npos )
                {
                    // skip ghost sets
                    std::map< moris::moris_index, Vector< moris::real > > tNewNodes = tSet->get_nodal_displacements( tRequestedIGNodes );
                    for ( auto const &[ tIndex, _ ] : tNewNodes )
                    {
                        tRequestedIGNodes.erase( tIndex );
                    }
                    tNodalDisplacements.merge( tNewNodes );
                }
            }
            MORIS_ASSERT( tRequestedIGNodes.size() == 0, "Not all requested nodal displacements could be found!" );


            // TODO @ff: Remove! Debug only
            //            mtk::Json_Debug_Output tDebugOutput( mMeshManager->get_integration_mesh( 0 ) );
            //            tDebugOutput.set_ig_vertex_displacements( tNodalDisplacements );
            //            uint const        tIteration = gLogger.get_iteration( "NonLinearAlgorithm", "Newton", "Solve" );
            //            std::string const tFileName  = "debug_mesh_" + std::to_string( tIteration ) + ".json";
            //            tDebugOutput.write_to_json( tFileName );


            // store the names of the mesh sets that are stored in each FEM set. This is necessary to update the newly created
            // nonconformal sets.
            Vector< std::string > tMeshSetNames;
            tMeshSetNames.reserve( mFemSets.size() );
            for ( auto const &tSet : mFemSets )
            {
                if ( tSet->is_empty_set() )
                {
                    tMeshSetNames.push_back( "" );
                }
                else
                {
                    tMeshSetNames.push_back( tSet->get_set_name() );
                }
            }

            MORIS_ASSERT( mContactMeshEditor != nullptr, "Contact mesh editor not initialized!" );
            mContactMeshEditor->update_displacements( tNodalDisplacements );
            mContactMeshEditor->update_nonconformal_side_sets();

            // loop over each fem set and check if it needs to be updated (i.e. only nonconformal sets!)
            for ( size_t iFemSet = 0; iFemSet < mFemSets.size(); ++iFemSet )
            {
                auto *const tFemSet      = dynamic_cast< fem::Set      *>( mFemSets( iFemSet ) );
                auto const &tMeshSetName = tMeshSetNames( iFemSet );
                if ( tFemSet->get_is_update_required() )
                {
                    tFemSet->set_mesh_set( mMeshManager->get_integration_mesh( 0 )->get_set_by_name( tMeshSetName ) );
                    tFemSet->update();
                }
            }
        }


        void
        FEM_Model::create_fem_sets(
                mtk::Interpolation_Mesh *aIPMesh,
                mtk::Integration_Mesh   *aIGMesh )
        {
            // get number of fem sets
            uint tNumFemSets = mSetInfo.size();

            // set size for list of equation sets
            mFemSets.resize( tNumFemSets, nullptr );

            // get number of IP cells
            uint tNumIPCells = aIPMesh->get_num_elems();

            // reserve size for list of equation objects
            mFemClusters.reserve( tNumIPCells );

            this->prepare_nonconformal_side_sets( aIGMesh );

            // loop over the used fem set
            for ( luint iSet = 0; iSet < tNumFemSets; iSet++ )
            {
                // get the mesh set name
                std::string tMeshSetName = mSetInfo( iSet ).get_mesh_set_name();

                // get the mesh set index from its name
                moris_index tMeshSetIndex = aIGMesh->get_set_index_by_name( tMeshSetName );

                // fill the mesh set index to fem set index map
                mMeshSetToFemSetMap[ std::make_tuple(
                        tMeshSetIndex,
                        mSetInfo( iSet ).get_time_continuity(),
                        mSetInfo( iSet ).get_time_boundary() ) ] = iSet;

                // get the mesh set pointer
                moris::mtk::Set *tMeshSet = aIGMesh->get_set_by_index( tMeshSetIndex );

                // get the number of clusters on the set
                uint tNumClustersOnSet = tMeshSet->get_num_clusters_on_set();

                // if clusters exist on the set, create a non-empty set
                if ( tNumClustersOnSet != 0 )
                {
                    // create new fem set
                    mFemSets( iSet ) = new fem::Set( this, tMeshSet, mSetInfo( iSet ), mIPNodes );
                }
                // if clusters don't exist on the set, create an empty set
                else
                {
                    // create an empty fem set
                    mFemSets( iSet ) = new fem::Set();
                }

                // add pointer back to this equation model to the FEM set
                mFemSets( iSet )->set_equation_model( this );

                // collect equation objects associated with the set
                Vector< MSI::Equation_Object * > const &tEquationObjectList = mFemSets( iSet )->get_equation_object_list();

                // add equation objects collected to the current set
                mFemClusters.append( tEquationObjectList );
            }

            // shrink to fit
            mFemClusters.shrink_to_fit();

            uint tNumElements       = mFemClusters.size();
            uint tGlobalNumElements = sum_all( tNumElements );

            // print output
            MORIS_LOG_SPEC( "Number of Clusters", tGlobalNumElements );
        }


        void
        FEM_Model::get_integration_xyz_active_flags(
                const Matrix< IndexMat >      &aNodeIndices,
                const Vector< gen::PDV_Type > &aRequestedPdvTypes,
                Matrix< DDSMat >              &aIsActiveDv )
        {
            // Get the number of node indices requested
            uint tNumIndices = aNodeIndices.length();

            // Get the number of dv types requested
            uint tNumTypes = aRequestedPdvTypes.size();

            // set size for is active dv types
            aIsActiveDv.set_size( tNumIndices, tNumTypes, 0 );

            // loop over the requested dv types
            for ( uint tNode = 0; tNode < tNumIndices; tNode++ )
            {
                // get node index
                uint tNodeIndex = aNodeIndices( tNode );

                // get flags from nodes
                Matrix< DDSMat > tIsActiveDvTemp;
                this->get_vertex_xyz_active_flags( mIGNodes( tNodeIndex )->get_index(), tIsActiveDvTemp, aRequestedPdvTypes );

                aIsActiveDv.get_row( tNode ) = tIsActiveDvTemp.matrix_data();
            }
        }


        void
        FEM_Model::get_integration_xyz_pdv_ids(
                const Matrix< IndexMat >      &aNodeIndices,
                const Vector< gen::PDV_Type > &aRequestedPdvTypes,
                Matrix< DDSMat >              &aXYZPdvIds )
        {
            // Get the number of node indices requested
            uint tNumIndices = aNodeIndices.length();

            // Get the number of dv types requested
            uint tNumTypes = aRequestedPdvTypes.size();

            // set size for is active dv types
            aXYZPdvIds.set_size( tNumIndices, tNumTypes, -1 );

            // loop over the requested dv types
            for ( uint tNode = 0; tNode < tNumIndices; tNode++ )
            {
                // get node index
                uint tNodeIndex = aNodeIndices( tNode );

                // get flags from nodes
                Matrix< DDSMat > tPdvIdsTemp;
                this->get_vertex_xyz_pdv_ids( mIGNodes( tNodeIndex )->get_index(), tPdvIdsTemp, aRequestedPdvTypes );

                aXYZPdvIds.get_row( tNode ) = tPdvIdsTemp.matrix_data();
            }
        }


        void
        FEM_Model::get_integration_xyz_pdv_active_flags_and_ids(
                const Matrix< IndexMat >      &aNodeIndices,
                const Vector< gen::PDV_Type > &aRequestedPdvTypes,
                Matrix< DDSMat >              &aIsActiveDv,
                Matrix< DDSMat >              &aXYZPdvIds )
        {
            // Get the number of node indices requested
            uint tNumIndices = aNodeIndices.length();

            // Get the number of dv types requested
            uint tNumTypes = aRequestedPdvTypes.size();

            // set size for is active dv types
            aIsActiveDv.set_size( tNumIndices, tNumTypes, 0 );

            // set size for is active dv types
            aXYZPdvIds.set_size( tNumIndices, tNumTypes, -1 );

            // loop over the requested dv types
            for ( uint tNode = 0; tNode < tNumIndices; tNode++ )
            {
                // get node index
                uint tNodeIndex = aNodeIndices( tNode );

                // get flags from nodes
                Matrix< DDSMat > tIsActiveDvTemp;
                this->get_vertex_xyz_active_flags( mIGNodes( tNodeIndex )->get_index(), tIsActiveDvTemp, aRequestedPdvTypes );

                aIsActiveDv.get_row( tNode ) = tIsActiveDvTemp.matrix_data();

                // get flags from nodes
                Matrix< DDSMat > tPdvIdsTemp;
                this->get_vertex_xyz_pdv_ids( mIGNodes( tNodeIndex )->get_index(), tPdvIdsTemp, aRequestedPdvTypes );
                aXYZPdvIds.get_row( tNode ) = tPdvIdsTemp.matrix_data();
            }
        }


        void
        FEM_Model::reset_integration_xyz_pdv_assembly_indices(
                const Matrix< IndexMat > &aNodeIndices )
        {
            // Get the number of node indices requested
            uint tNumIndices = aNodeIndices.length();

            // loop over the requested dv types
            for ( uint tNode = 0; tNode < tNumIndices; tNode++ )
            {
                // get node index
                uint        tNodeIndex   = aNodeIndices( tNode );
                moris_index tIgNodeIndex = mIGNodes( tNodeIndex )->get_index();

                // reset with -1
                mXYZLocalAssemblyIndices.get_row( tIgNodeIndex ).fill( -1 );
            }
        }


        void
        FEM_Model::get_integration_xyz_pdv_assembly_indices(
                const Matrix< IndexMat >      &aNodeIndices,
                const Vector< gen::PDV_Type > &aRequestedPdvTypes,
                Matrix< DDSMat >              &aXYZPdvAssemblyIndices )
        {
            // Get the number of node indices requested
            uint tNumIndices = aNodeIndices.length();

            // Get the number of dv types requested
            uint tNumTypes = aRequestedPdvTypes.size();

            // set size for is active dv types
            aXYZPdvAssemblyIndices.set_size( tNumIndices, tNumTypes, -1 );

            // loop over the requested dv types
            for ( uint tNode = 0; tNode < tNumIndices; tNode++ )
            {
                // get node index
                uint tNodeIndex = aNodeIndices( tNode );

                // get assembly index from nodes
                Matrix< DDSMat > tPdvAssemblyIndicesTemp;
                moris_index      tIgNodeIndex = mIGNodes( tNodeIndex )->get_index();
                this->get_local_xyz_pdv_assembly_indices( tIgNodeIndex, tPdvAssemblyIndicesTemp, aRequestedPdvTypes );

                aXYZPdvAssemblyIndices.get_row( tNode ) = tPdvAssemblyIndicesTemp.matrix_data();
            }
        }


        Matrix< DDSMat >
        FEM_Model::get_XYZ_local_pdv_assembly_map() const
        {
            return mXYZLocalAssemblyIndices;
        }


        void
        FEM_Model::set_integration_xyz_pdv_assembly_index(
                moris_index        aNodeIndex,
                enum gen::PDV_Type aPdvType,
                moris_index        aXYZPdvAssemblyIndex )
        {
            // get the index of the underlying node
            moris_index tIgNodeIndex = mIGNodes( aNodeIndex )->get_index();

            // sanity check the input
            MORIS_ASSERT( aXYZPdvAssemblyIndex != MORIS_INDEX_MAX,
                    "FEM_Model::set_integration_xyz_pdv_assembly_index() - "
                    "Trying to assign MORIS_INDEX_MAX in XYZ Local PDV assembly map. This shouldn't happen." );

            // assign assembly index in map
            mXYZLocalAssemblyIndices( tIgNodeIndex, (uint)aPdvType ) = aXYZPdvAssemblyIndex;
        }

        void
        FEM_Model::initialize_from_inputfile( std::shared_ptr< Library_IO > aLibrary )
        {
            Tracer tTracer( "FEM", "Model", "Initialize from Input File" );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 0: unpack fem input and mesh
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // get pointers to interpolation and integration meshes
            mtk::Interpolation_Mesh *tIPMesh = nullptr;
            mtk::Integration_Mesh   *tIGMesh = nullptr;
            mMeshManager->get_mesh_pair( mMeshPairIndex, tIPMesh, tIGMesh );

            // set the space dimension
            mSpaceDim = tIPMesh->get_spatial_dim();

            // unpack the FEM inputs
            this->initialize( aLibrary );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 1: create interpolation nodes
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            this->create_interpolation_nodes( tIPMesh );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 2: create integration nodes
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            this->create_integration_nodes( tIGMesh );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 3: create fem sets
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            this->create_fem_sets( tIPMesh, tIGMesh );
        }


        void
        FEM_Model::initialize( const std::shared_ptr< Library_IO > &aLibrary )
        {
            Tracer tTracer( "FEM", "Model", "Load and Initialize Parameters" );

            // in the INT test of the IQI strain energy, the FEM Model is initialized without a MeshManager.
            mtk::Mesh_Pair const *tMeshPair = nullptr;
            if ( mMeshManager != nullptr )
            {
                tMeshPair = &mMeshManager->get_mesh_pair( mMeshPairIndex );
            }

            /**
             * @brief The old way to input the IWGs uses the explicit mesh names without any phase definitions (e.g. iside_b0_1_b1_0, ...),
             * while the new method uses phases and phase-pairs to define the applicable sets. Old input files can be detected by the number of
             * elements in the ParameterList. If it is 8, then the legacy method was used, if it is 9, the new method was used.
             */
            std::unique_ptr< Model_Initializer > tModelInitializer;
            switch ( mParameterList.size() )
            {
                case 8:
                {
                    tModelInitializer = std::make_unique< Model_Initializer_Legacy >(
                            mParameterList,
                            aLibrary,
                            tMeshPair,
                            mSpaceDim,
                            mUseNewGhostSets,
                            mDofTypeToBsplineMeshIndex );
                    break;
                }
                case 9:
                {
                    tModelInitializer = std::make_unique< Model_Initializer_Phasebased >(
                            mParameterList,
                            aLibrary,
                            tMeshPair,
                            mSpaceDim,
                            mUseNewGhostSets,
                            mDofTypeToBsplineMeshIndex );
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "FEM_Model::initialize - wrong size for parameter list: %zu", mParameterList.size() );
                }
            }
            tModelInitializer->initialize();
            mSetInfo    = tModelInitializer->get_set_info();
            mIQIs       = tModelInitializer->get_iqis();
            mFields     = tModelInitializer->get_fields();
            mFieldTypes = tModelInitializer->get_field_types();
        }

        FEM_Model::~FEM_Model()
        {
            this->free_memory();
        }

        void
        FEM_Model::free_memory()
        {
            // delete fem nodes
            for ( auto tIPNodes : mIPNodes )
            {
                delete tIPNodes;
            }
            mIPNodes.clear();

            // delete fem nodes
            for ( auto tIGNodes : mIGNodes )
            {
                delete tIGNodes;
            }
            mIGNodes.clear();

            // delete the fem sets
            for ( auto tFemSet : mFemSets )
            {
                delete tFemSet;
            }
            mFemSets.clear();

            // delete the fem cluster
            mFemClusters.clear();

            // delete the matrices data
            mIsActiveXYZ.set_size( 0, 0 );
            mXYZPdvIds.set_size( 0, 0 );
            mXYZLocalAssemblyIndices.set_size( 0, 0 );
        }

        void
        FEM_Model::finalize_equation_sets(
                MSI::Model_Solver_Interface *aModelSolverInterface )
        {
            // loop over the fem sets
            for ( MSI::Equation_Set *tFemSet : mFemSets )
            {
                // finalize the fem set
                tFemSet->finalize( aModelSolverInterface );
            }
        }

        void
        FEM_Model::normalize_IQIs()
        {
            for ( uint tRequestedIQIIndex = 0; tRequestedIQIIndex < mRequestedIQINames.size(); tRequestedIQIIndex++ )
            {
                // IQI index
                uint tIQIIndex = 0;
                while ( tIQIIndex < mIQIs.size() and ( mIQIs( tIQIIndex )->get_name() not_eq mRequestedIQINames( tRequestedIQIIndex ) ) )
                {
                    tIQIIndex++;
                }
                MORIS_ASSERT( tIQIIndex < mIQIs.size(),
                        "IQI was not found with the requested name %s",
                        mRequestedIQINames( tRequestedIQIIndex ).c_str() );

                // Set normalization
                std::string tNormalization = mParameterList( 4 )( tIQIIndex ).get< std::string >( "normalization" );
                if ( tNormalization == "none" )
                {
                    // Do nothing
                }
                else if ( tNormalization == "time" )
                {
                    MORIS_ERROR( false, "Time normalization not implemented yet for IQIs yet, implementation should go here." );
                }
                else if ( tNormalization == "design" )
                {
                    mIQIs( tRequestedIQIIndex )->set_reference_value( mGlobalIQIVal( tRequestedIQIIndex )( 0 ) );
                    mGlobalIQIVal( tRequestedIQIIndex )( 0 ) = 1.0;
                }
                else
                {
                    // Try to set reference values directly
                    try
                    {
                        mIQIs( tRequestedIQIIndex )->set_reference_value( string_to_mat< DDRMat >( tNormalization )( 0 ) );
                    } catch ( ... )
                    {
                        // error
                        MORIS_ERROR( false,
                                "FEM_Model::normalize_IQIs() - Unknown normalization: %s. Must be 'none', 'time', 'design', or a reference value.",
                                tNormalization.c_str() );
                    }
                }
            }
        }


        const std::shared_ptr< fem::Field > &
        FEM_Model::get_field( mtk::Field_Type tFieldType )
        {
            size_t tIndex = mFieldTypes( static_cast< sint >( tFieldType ) );
            return mFields( tIndex );
        }


        Vector< std::shared_ptr< mtk::Field > >
        FEM_Model::get_fields()
        {
            Vector< std::shared_ptr< mtk::Field > > tFields( mFields.size(), nullptr );

            for ( uint Ik = 0; Ik < mFields.size(); Ik++ )
            {
                tFields( Ik ) = mFields( Ik );
            }

            return tFields;
        }

        void
        FEM_Model::populate_fields()
        {
            Tracer tTracer( "FEM", "Model", "Populate fields" );

            // check if fields exists
            if ( mFields.size() == 0 )
            {
                return;
            }

            Vector< std::shared_ptr< fem::Field > > tFieldToPopulate;
            Vector< std::string >                   tFieldIQINames;

            tFieldToPopulate.reserve( mFields.size() );
            tFieldIQINames.reserve( mFields.size() );

            for ( uint iField = 0; iField < mFields.size(); iField++ )
            {
                if ( mFields( iField )->get_populate_field_with_IQI() )
                {
                    tFieldToPopulate.push_back( mFields( iField ) );
                    tFieldIQINames.push_back( mFields( iField )->get_IQI_name() );
                }
            }

            for ( uint iSet = 0; iSet < mFemSets.size(); iSet++ )
            {
                if ( mFemSets( iSet )->get_element_type() == Element_Type::BULK )
                {
                    mFemSets( iSet )->populate_fields( tFieldToPopulate, tFieldIQINames );
                }
            }

            // output fields to file. // FIXME: better should be done somewhere else.
            for ( uint iSet = 0; iSet < mFields.size(); iSet++ )
            {
                mFields( iSet )->output_field_to_file();

                // mFields( iSet )->save_field_to_exodus( "FEM_Field.Exo" );
            }
        }


        void
        FEM_Model::create_IQI_map()
        {
            // erase the content of the map
            mIQINameToIndexMap.clear();

            uint tCounter = 0;

            // loop over all IQIs and build a name to index map
            for ( const std::shared_ptr< IQI > &tIQI : mIQIs )
            {
                std::string tIQIName = tIQI->get_name();

                mIQINameToIndexMap[ tIQIName ] = tCounter++;
            }
        }


        void
        FEM_Model::initialize_IQIs()
        {
            // create content of the map
            this->create_IQI_map();

            // number of requested IQIs for the model
            moris::uint tNumIQIs = mRequestedIQINames.size();

            // resize the Global IQI cell
            mGlobalIQIVal.resize( tNumIQIs );

            // set a counter
            uint iIQICounter = 0;

            for ( auto &tQI : mGlobalIQIVal )
            {
                // use the map to get the IQI index of requested IQI
                moris_index tIQIIndex = mIQINameToIndexMap[ mRequestedIQINames( iIQICounter ) ];

                // get the matrix dimension of the IQI
                std::pair< uint, uint > tMatrixDim = mIQIs( tIQIIndex )->get_matrix_dim();

                // set size for the QI value and initialize the value with 0
                tQI.set_size( tMatrixDim.first, tMatrixDim.second, 0.0 );

                // increase the counter by 1
                iIQICounter++;
            }
        }


        void
        FEM_Model::get_vertex_xyz_active_flags(
                moris_index                         aVertexIndex,
                Matrix< DDSMat >                   &aIsActiveDv,
                const Vector< enum gen::PDV_Type > &aPdvTypes )
        {
            // get number of requested pdv types
            uint tNumPdvTypes = aPdvTypes.size();

            // set size for active flag
            aIsActiveDv.set_size( 1, tNumPdvTypes );

            // loop over requested pdv types
            for ( uint iPdvType = 0; iPdvType < tNumPdvTypes; iPdvType++ )
            {
                // get pdv index
                uint tXYZIndex = static_cast< uint >( aPdvTypes( iPdvType ) );

                // set value
                aIsActiveDv( iPdvType ) = mIsActiveXYZ( aVertexIndex, tXYZIndex );
            }
        }


        void
        FEM_Model::set_vertex_xyz_active_flags(
                moris_index               aVertexIndex,
                Vector< Vector< bool > > &aIsActiveDv )
        {
            // get num of pdv
            uint tNumXYZPdv = aIsActiveDv.size();

            // fill mIsActiveXYZ
            for ( uint iPdvType = 0; iPdvType < tNumXYZPdv; iPdvType++ )
            {
                mIsActiveXYZ( aVertexIndex, iPdvType ) = aIsActiveDv( iPdvType )( 0 );
            }
        }


        void
        FEM_Model::set_vertex_xyz_pdv_ids(
                moris_index                 aVertexIndex,
                Vector< Matrix< DDSMat > > &aXYZPvIds )
        {
            // get num of pdv
            uint tNumXYZPdv = aXYZPvIds.size();

            // fill mXYZPdvIds
            for ( uint iPdvType = 0; iPdvType < tNumXYZPdv; iPdvType++ )
            {
                mXYZPdvIds( aVertexIndex, iPdvType ) = aXYZPvIds( iPdvType )( 0 );
            }
        }


        void
        FEM_Model::get_vertex_xyz_pdv_ids(
                moris_index                         aVertexIndex,
                Matrix< DDSMat >                   &aXYZPdvIds,
                const Vector< enum gen::PDV_Type > &aPdvTypes )
        {
            // get number of requested pdv types
            uint tNumPdvTypes = aPdvTypes.size();

            // set size for pdv ids
            aXYZPdvIds.set_size( 1, tNumPdvTypes, -1 );

            // loop over requested pdv types
            for ( uint iPdvType = 0; iPdvType < tNumPdvTypes; iPdvType++ )
            {
                // get pdv index
                uint tXYZIndex = static_cast< uint >( aPdvTypes( iPdvType ) );

                // set value
                aXYZPdvIds( iPdvType ) = mXYZPdvIds( aVertexIndex, tXYZIndex );
            }
        }


        void
        FEM_Model::get_local_xyz_pdv_assembly_indices(
                moris_index                         aVertexIndex,
                Matrix< DDSMat >                   &aXYZLocalAssemblyIndices,
                const Vector< enum gen::PDV_Type > &aPdvTypes )
        {
            // get number of requested pdv types
            uint tNumPdvTypes = aPdvTypes.size();

            // set size for pdv ids
            aXYZLocalAssemblyIndices.set_size( 1, tNumPdvTypes, -1 );

            // loop over requested pdv types
            for ( uint iPdvType = 0; iPdvType < tNumPdvTypes; iPdvType++ )
            {
                // get pdv index
                uint tXYZIndex = static_cast< uint >( aPdvTypes( iPdvType ) );

                // get the assembly index
                moris_index tAssemblyIndex = mXYZLocalAssemblyIndices( aVertexIndex, tXYZIndex );

                // sanity check the output
                MORIS_ASSERT( tAssemblyIndex != MORIS_INDEX_MAX,
                        "FEM_Model::get_local_xyz_pdv_assembly_indices() - "
                        "Local PDV assembly index is MORIS_INDEX_MAX. It has likely not been assigned." );

                // set value
                aXYZLocalAssemblyIndices( iPdvType ) = tAssemblyIndex;
            }
        }


        /**
         * @brief Set the dof type to Bspline mesh index map
         * @param aDofTypeToBsplineMeshIndex
         */
        void
        FEM_Model::set_dof_type_to_Bspline_mesh_index(
                std::unordered_map< MSI::Dof_Type, moris_index > aDofTypeToBsplineMeshIndex )
        {
            mDofTypeToBsplineMeshIndex = aDofTypeToBsplineMeshIndex;
        }


        /**
         * @brief set flag whether to use new ghost sets
         * @param aUseNewGhostSets
         */
        void
        FEM_Model::set_use_new_ghost_sets( bool aUseNewGhostSets )
        {
            mUseNewGhostSets = aUseNewGhostSets;
        }

        void FEM_Model::prepare_nonconformal_side_sets( mtk::Integration_Mesh *aIGMesh )
        {
            // early return if no nonconformal set is present
            if ( std::none_of( mSetInfo.begin(), mSetInfo.end(), []( auto &tSetInfo ) { return tSetInfo.get_mesh_set_name().find( "ncss" ) != std::string::npos; } ) )
            {
                return;
            }

            auto *tIGMesh = dynamic_cast< mtk::Integration_Mesh_DataBase_IG * >( aIGMesh );

            // TODO @ff: remove! Only for Debug.
//            std::string const tFileName = "debug_mesh_0.json";
//            mtk::Json_Debug_Output( tIGMesh ).write_to_json( tFileName );

            auto const &[ tSetNames, tCandidatePairs ] = prepare_nonconformal_candidate_pairs();

            Vector< mtk::Side_Set const * > tSideSets;
            for ( auto &tSetName : tSetNames )
            {
                tSideSets.push_back( dynamic_cast< mtk::Side_Set const * >( tIGMesh->get_set_by_name( tSetName ) ) );
            }

            mtk::Integrator tSideIntegrator = prepare_nonconformal_integrator( tIGMesh );

            auto tCMEditor = std::make_shared< mtk::Contact_Mesh_Editor >( tIGMesh, tSideIntegrator, tSideSets, tCandidatePairs );
            tCMEditor->update_nonconformal_side_sets();
            this->set_contact_mesh_editor( tCMEditor );
        }

        mtk::Integrator FEM_Model::prepare_nonconformal_integrator( mtk::Integration_Mesh const *aIGMesh )
        {
            mtk::Integration_Order tIntegrationOrder = mtk::Integration_Order::UNDEFINED;
            for ( auto const &tSetInfo : mSetInfo )
            {
                if ( tSetInfo.get_integration_order() != mtk::Integration_Order::UNDEFINED )
                {
                    tIntegrationOrder = tSetInfo.get_integration_order();
                    break;
                }
            }
            MORIS_ASSERT( tIntegrationOrder not_eq mtk::Integration_Order::UNDEFINED, "Nonconformal integration order not defined!" );
            MORIS_ASSERT( aIGMesh->get_spatial_dim() == 2, "Currently only 2D problems are supported." );
            // create a side integrator
            mtk::Integration_Rule tSideIntegRule(
                    mtk::Geometry_Type::LINE,    // TODO @ff: currently only 2d problems are supported
                    mtk::Integration_Type::GAUSS,
                    tIntegrationOrder,
                    mtk::Geometry_Type::LINE,
                    mtk::Integration_Type::GAUSS,
                    mtk::Integration_Order::BAR_1 );
            mtk::Integrator tSideIntegrator( tSideIntegRule );
            return tSideIntegrator;
        }

        std::pair< Vector< std::string >, Vector< std::pair< moris_index, moris_index > > >
        FEM_Model::prepare_nonconformal_candidate_pairs()
        {
            Vector< std::string >                           tRequiredSideSetNames;          // list of unique names that will be used to create nonconformal sets
            Vector< std::pair< moris_index, moris_index > > tNonconformalCandidatePairs;    // list of indices into the required side set names to get possible nonconformal pairs
            std::unordered_map< std::string, moris_index >  tNonconformalIndexMap;          // temporary map to store the index of the required side set names
            for ( auto &tSetInfo : mSetInfo )
            {
                if ( std::string const tMeshSetName = tSetInfo.get_mesh_set_name();
                        tMeshSetName.find( "ncss" ) != std::string::npos )
                {
                    // format 'ncss|iside_b0_0_b1_1|iside_b0_1_b1_0'
                    std::string const tName        = tMeshSetName.substr( tMeshSetName.find( '|' ) + 1 );
                    size_t const      tPos         = tName.find( '|' );
                    std::string const tLeaderSet   = tName.substr( 0, tPos );
                    std::string const tFollowerSet = tName.substr( tPos + 1 );

                    if ( tNonconformalIndexMap.find( tLeaderSet ) == tNonconformalIndexMap.end() )
                    {
                        tNonconformalIndexMap[ tLeaderSet ] = tRequiredSideSetNames.size();
                        tRequiredSideSetNames.push_back( tLeaderSet );
                    }

                    if ( tNonconformalIndexMap.find( tFollowerSet ) == tNonconformalIndexMap.end() )
                    {
                        tNonconformalIndexMap[ tFollowerSet ] = tRequiredSideSetNames.size();
                        tRequiredSideSetNames.push_back( tFollowerSet );
                    }

                    tNonconformalCandidatePairs.push_back( { tNonconformalIndexMap[ tLeaderSet ], tNonconformalIndexMap[ tFollowerSet ] } );
                }
            }
            return { tRequiredSideSetNames, tNonconformalCandidatePairs };
        }
    }    // namespace fem
} /* namespace moris */