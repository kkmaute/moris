/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Model.cpp
 *
 */

#include <map>
#ifdef WITHGPERFTOOLS
#include <gperftools/profiler.h>
#endif

// LINALG/src
#include "cl_Map.hpp"
#include "cl_Matrix.hpp"
#include "fn_unique.hpp"
#include "fn_sum.hpp"    // for check
#include "fn_iscol.hpp"
#include "fn_trans.hpp"
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
#include "cl_GEN_Pdv_Enums.hpp"
// Logging package
#include "cl_Logger.hpp"
#include "cl_Tracer.hpp"

namespace moris
{
    namespace fem
    {
        // User-defined FEM function
        typedef void ( *FEM_Function )(
                moris::Matrix< moris::DDRMat >                &aPropMatrix,
                moris::Cell< moris::Matrix< moris::DDRMat > > &aParameters,
                moris::fem::Field_Interpolator_Manager        *aFIManager );
        //------------------------------------------------------------------------------

        FEM_Model::FEM_Model(
                std::shared_ptr< mtk::Mesh_Manager > aMeshManager,
                const moris_index                   &aMeshPairIndex,
                moris::Cell< fem::Set_User_Info >   &aSetInfo )
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

        //------------------------------------------------------------------------------

        FEM_Model::FEM_Model(
                std::shared_ptr< mtk::Mesh_Manager > aMeshManager,
                const moris_index                   &aMeshPairIndex,
                moris::Cell< fem::Set_User_Info >   &aSetInfo,
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

        //------------------------------------------------------------------------------

        FEM_Model::FEM_Model(
                std::shared_ptr< mtk::Mesh_Manager >               aMeshManager,
                const moris_index                                 &aMeshPairIndex,
                const moris::Cell< moris::Cell< ParameterList > > &aParameterList,
                std::shared_ptr< Library_IO >                      aLibrary )
                : mMeshManager( aMeshManager )
                , mMeshPairIndex( aMeshPairIndex )
                , mParameterList( aParameterList )
        {
            Tracer tTracer( "FEM", "Model", "Create" );

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
            // STEP 1: create IP nodes
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            this->create_interpolation_nodes( tIPMesh );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 2: create fem sets
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            this->create_fem_sets( tIPMesh, tIGMesh );
        }

        //------------------------------------------------------------------------------

        FEM_Model::FEM_Model(
                std::shared_ptr< mtk::Mesh_Manager >        aMeshManager,
                const moris_index                          &aMeshPairIndex,
                moris::Cell< moris::Cell< ParameterList > > aParameterList,
                MSI::Design_Variable_Interface             *aDesignVariableInterface )
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

        //------------------------------------------------------------------------------

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

        //------------------------------------------------------------------------------

        void
        FEM_Model::create_integration_nodes( mtk::Integration_Mesh *aIGMesh )
        {
            // ask IG mesh about number of IG vertices on proc
            luint tNumIGNodes = aIGMesh->get_num_nodes();

            // set size for list IG nodes
            mIGNodes.resize( tNumIGNodes, nullptr );

            moris::Cell< enum PDV_Type > tGeoPdvType;
            switch ( mSpaceDim )
            {
                case 1:
                    tGeoPdvType = { PDV_Type::X_COORDINATE };
                    break;
                case 2:
                    tGeoPdvType = { PDV_Type::X_COORDINATE, PDV_Type::Y_COORDINATE };
                    break;
                case 3:
                    tGeoPdvType = { PDV_Type::X_COORDINATE, PDV_Type::Y_COORDINATE, PDV_Type::Z_COORDINATE };
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
            // moris::Cell<moris_index> tMesh_GEN_map;
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
                    Matrix< IndexMat >              tVertexIndices( 1, 1, tVertexIndex );
                    moris::Cell< Matrix< DDRMat > > tVertexCoordsFromGen( mSpaceDim );
                    moris::Cell< Matrix< DDSMat > > tIsActiveDv;

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
                    moris::Cell< moris::Matrix< DDSMat > > tPdvIds;
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

        //------------------------------------------------------------------------------

        void
        FEM_Model::create_fem_sets(
                mtk::Interpolation_Mesh           *aIPMesh,
                mtk::Integration_Mesh             *aIGMesh,
                moris::Cell< fem::Set_User_Info > &aSetInfo )
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

        //------------------------------------------------------------------------------

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
                Cell< MSI::Equation_Object * > const &tEquationObjectList = mFemSets( iSet )->get_equation_object_list();

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

        //------------------------------------------------------------------------------

        void
        FEM_Model::get_integration_xyz_active_flags(
                const Matrix< IndexMat >      &aNodeIndices,
                const moris::Cell< PDV_Type > &aRequestedPdvTypes,
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

        //------------------------------------------------------------------------------

        void
        FEM_Model::get_integration_xyz_pdv_ids(
                const Matrix< IndexMat >      &aNodeIndices,
                const moris::Cell< PDV_Type > &aRequestedPdvTypes,
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

        //------------------------------------------------------------------------------

        void
        FEM_Model::get_integration_xyz_pdv_active_flags_and_ids(
                const Matrix< IndexMat >      &aNodeIndices,
                const moris::Cell< PDV_Type > &aRequestedPdvTypes,
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

        //------------------------------------------------------------------------------

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

        //------------------------------------------------------------------------------

        void
        FEM_Model::get_integration_xyz_pdv_assembly_indices(
                const Matrix< IndexMat >      &aNodeIndices,
                const moris::Cell< PDV_Type > &aRequestedPdvTypes,
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

        //------------------------------------------------------------------------------

        Matrix< DDSMat >
        FEM_Model::get_XYZ_local_pdv_assembly_map() const
        {
            return mXYZLocalAssemblyIndices;
        }

        //------------------------------------------------------------------------------

        void
        FEM_Model::set_integration_xyz_pdv_assembly_index(
                moris_index   aNodeIndex,
                enum PDV_Type aPdvType,
                moris_index   aXYZPdvAssemblyIndex )
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

        //------------------------------------------------------------------------------

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

        //------------------------------------------------------------------------------

        void
        FEM_Model::initialize( const std::shared_ptr< Library_IO > &aLibrary )
        {
            Tracer tTracer( "FEM", "Model", "Load and Initialize Parameters" );

            // get msi string to dof type map
            moris::map< std::string, MSI::Dof_Type > tMSIDofTypeMap =
                    moris::MSI::get_msi_dof_type_map();

            // get string to dv type map
            moris::map< std::string, PDV_Type > tMSIDvTypeMap =
                    get_pdv_type_map();

            // get string to field type map
            moris::map< std::string, mtk::Field_Type > tFieldTypeMap =
                    mtk::get_field_type_map();

            using StringToIntMap = std::map< std::string, uint >;

            /**
             * @brief The old way to input the IWGs uses the explicit mesh names without any phase definitions (e.g. iside_b0_1_b1_0, ...),
             * while the new method uses phases and phase-pairs to define the applicable sets. Old input files can be detected by the number of
             * elements in the ParameterList. If it is 8, then the legacy method was used, if it is 9, the new method was used.
             */
            bool tLegacyParameterInput = false;
            switch ( mParameterList.size() )
            {
                case 8:
                {
                    tLegacyParameterInput = true;    // without phase
                    break;
                }
                case 9:
                {
                    tLegacyParameterInput = false;    // with phase
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "FEM_Model::initialize - wrong size for parameter list" );
                }
            }

            StringToIntMap tPropertyMap = this->create_properties( tMSIDofTypeMap, tMSIDvTypeMap, tFieldTypeMap, aLibrary );
            this->create_fields();
            if ( tLegacyParameterInput )
            {
                StringToIntMap tMMMap = this->create_material_models_without_phase( tPropertyMap, tMSIDofTypeMap, tMSIDvTypeMap );
                StringToIntMap tCMMap = this->create_constitutive_models_without_phase( tPropertyMap, tMMMap, tMSIDofTypeMap, tMSIDvTypeMap );
                StringToIntMap tSPMap = this->create_stabilization_parameters_without_phase( tPropertyMap, tCMMap, tMSIDofTypeMap, tMSIDvTypeMap );

                this->create_IWGs_without_phase( tPropertyMap, tMMMap, tCMMap, tSPMap, tMSIDofTypeMap, tMSIDvTypeMap, tFieldTypeMap );
                this->create_IQIs_without_phase( tPropertyMap, tCMMap, tSPMap, tMSIDofTypeMap, tMSIDvTypeMap, tFieldTypeMap );
                this->create_fem_set_info_without_phase();
            }
            else
            {
                // create phases
                StringToIntMap tSPMap = this->create_stabilization_parameters( tPropertyMap, tMSIDofTypeMap, tMSIDvTypeMap );

                this->create_phases();
                this->create_material_models( tPropertyMap, tMSIDofTypeMap, tMSIDvTypeMap );
                this->create_constitutive_models( tPropertyMap, tMSIDofTypeMap, tMSIDvTypeMap );
                this->create_IWGs( tPropertyMap, tSPMap, tMSIDofTypeMap );
                this->create_IQIs( tPropertyMap, tSPMap, tMSIDofTypeMap );
                this->create_fem_set_info();
            }
            this->print_physics_model( !tLegacyParameterInput );
        }

        void FEM_Model::print_physics_model( bool aWithPhase )
        {
            ParameterList tComputationParameterList = this->mParameterList( 5 )( 0 );
            bool          tPrintPhysics             = tComputationParameterList.get< bool >( "print_physics_model" );
            if ( tPrintPhysics && par_rank() == 0 )
            {
                if ( aWithPhase )
                {
                    std::cout << "Phase info \n";
                    for ( auto &tPhaseInfo : mPhaseInfo )
                    {
                        std::cout << "%-------------------------------------------------\n";
                        tPhaseInfo.print_names();
                        std::cout << "%-------------------------------------------------\n";
                    }
                    std::cout << " \n";
                }

                std::cout << "Set info \n";
                for ( auto &tSetInfo : mSetInfo )
                {
                    std::cout << "%-------------------------------------------------\n";
                    tSetInfo.print_names();
                    std::cout << "%-------------------------------------------------\n";
                }
            }
        }

        //------------------------------------------------------------------------------

        FEM_Model::~FEM_Model()
        {
            this->free_memory();
        }
        //------------------------------------------------------------------------------

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

        //------------------------------------------------------------------------------

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

        //------------------------------------------------------------------------------

        std::map< std::string, uint >
        FEM_Model::create_properties(
                moris::map< std::string, MSI::Dof_Type >   &aMSIDofTypeMap,
                moris::map< std::string, PDV_Type >        &aDvTypeMap,
                moris::map< std::string, mtk::Field_Type > &aFieldTypeMap,
                std::shared_ptr< Library_IO >               aLibrary )
        {
            std::map< std::string, uint > tPropertyMap;

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
                tPropertyMap[ tPropertyName ] = iProp;

                // set dof dependencies
                moris::Cell< moris::Cell< moris::MSI::Dof_Type > > tDofTypes;
                string_to_cell_of_cell(
                        tPropParameter.get< std::string >( "dof_dependencies" ),
                        tDofTypes,
                        aMSIDofTypeMap );
                mProperties( iProp )->set_dof_type_list( tDofTypes );

                // set dv dependencies
                moris::Cell< moris::Cell< PDV_Type > > tDvTypes;
                string_to_cell_of_cell(
                        tPropParameter.get< std::string >( "dv_dependencies" ),
                        tDvTypes,
                        aDvTypeMap );
                mProperties( iProp )->set_dv_type_list( tDvTypes );

                // set field dependencies
                moris::Cell< moris::Cell< mtk::Field_Type > > tFieldTypes;
                string_to_cell_of_cell(
                        tPropParameter.get< std::string >( "field_dependencies" ),
                        tFieldTypes,
                        aFieldTypeMap );
                mProperties( iProp )->set_field_type_list( tFieldTypes );

                // set function parameters
                moris::Cell< moris::Matrix< DDRMat > > tFuncParameters;
                string_to_cell_mat_2(
                        tPropParameter.get< std::string >( "function_parameters" ),
                        tFuncParameters );
                mProperties( iProp )->set_parameters( tFuncParameters );

                // set value function for property
                std::string  tValFuncName = tPropParameter.get< std::string >( "value_function" );
                FEM_Function tValFunction = nullptr;
                if ( tValFuncName.size() > 1 )
                {
                    tValFunction = aLibrary->load_function< FEM_Function >( tValFuncName );
                    mProperties( iProp )->set_val_function( tValFunction );
                }

                // set dof derivative function for property
                moris::Cell< std::string > tDofDerFuncNames;
                string_to_cell(
                        tPropParameter.get< std::string >( "dof_derivative_functions" ),
                        tDofDerFuncNames );
                uint                             tNumDofDerFuncs = tDofDerFuncNames.size();
                moris::Cell< fem::PropertyFunc > tDofDerFunctions( tNumDofDerFuncs, nullptr );
                for ( uint iFunc = 0; iFunc < tNumDofDerFuncs; iFunc++ )
                {
                    if ( tDofDerFuncNames( iFunc ).size() > 1 )
                    {
                        tDofDerFunctions( iFunc ) = aLibrary->load_function< FEM_Function >( tDofDerFuncNames( iFunc ) );
                    }
                }
                mProperties( iProp )->set_dof_derivative_functions( tDofDerFunctions );

                // set dv derivative function for property
                moris::Cell< std::string > tDvDerFuncNames;
                string_to_cell(
                        tPropParameter.get< std::string >( "dv_derivative_functions" ),
                        tDvDerFuncNames );
                uint                             tNumDvDerFuncs = tDvDerFuncNames.size();
                moris::Cell< fem::PropertyFunc > tDvDerFunctions( tNumDvDerFuncs, nullptr );
                for ( uint iFunc = 0; iFunc < tNumDvDerFuncs; iFunc++ )
                {
                    if ( tDvDerFuncNames( iFunc ).size() > 1 )
                    {
                        tDvDerFunctions( iFunc ) = aLibrary->load_function< FEM_Function >( tDvDerFuncNames( iFunc ) );
                    }
                }
                mProperties( iProp )->set_dv_derivative_functions( tDvDerFunctions );

                // set space derivative function for property
                moris::Cell< std::string > tSpaceDerFuncNames;
                string_to_cell(
                        tPropParameter.get< std::string >( "space_derivative_functions" ),
                        tSpaceDerFuncNames );
                uint tNumSpaceDerFuncs = tSpaceDerFuncNames.size();

                moris::Cell< fem::PropertyFunc > tSpaceDerFunctions( tNumSpaceDerFuncs, nullptr );
                for ( uint iFunc = 0; iFunc < tNumSpaceDerFuncs; iFunc++ )
                {
                    if ( tSpaceDerFuncNames( iFunc ).size() > 1 )
                    {
                        tSpaceDerFunctions( iFunc ) = aLibrary->load_function< FEM_Function >( tSpaceDerFuncNames( iFunc ) );
                    }
                }
                mProperties( iProp )->set_space_der_functions( tSpaceDerFunctions );
            }

            return tPropertyMap;
        }

        //------------------------------------------------------------------------------

        void
        FEM_Model::create_fields()
        {
            std::map< std::string, uint > tFieldMap;

            // get the property parameter list
            moris::Cell< ParameterList > tFieldParameterList = mParameterList( 6 );

            // get the number of properties
            sint tNumFields = tFieldParameterList.size();

            // create a list of property pointers
            mFields.resize( tNumFields, nullptr );

            // loop over the parameter lists
            for ( sint iFields = 0; iFields < tNumFields; iFields++ )
            {
                // get property parameter list
                ParameterList tFieldParameter = tFieldParameterList( iFields );

                // get property name from parameter list
                std::string tFieldName = tFieldParameter.get< std::string >( "field_name" );

                moris::map< std::string, mtk::Field_Entity_Type > tFieldEntityTypeMap =
                        mtk::get_field_entity_type_map();

                mtk::Field_Entity_Type tFieldEntityType =
                        tFieldEntityTypeMap.find( tFieldParameter.get< std::string >( "field_entity_type" ) );

                // create a property pointer
                std::shared_ptr< fem::Field > tField = std::make_shared< fem::Field >(
                        mMeshManager->get_mesh_pair( mMeshPairIndex ),
                        tFieldEntityType );

                // set a name for the property
                tField->set_label( tFieldName );

                // fill property map
                tFieldMap[ tFieldName ] = iFields;

                // set field type
                moris::map< std::string, mtk::Field_Type > tFieldTypeMap =
                        mtk::get_field_type_map();

                moris::Cell< mtk::Field_Type > tFieldTypes;
                string_to_cell(
                        tFieldParameter.get< std::string >( "field_type" ),
                        tFieldTypes,
                        tFieldTypeMap );

                // set field type
                tField->set_field_type( tFieldTypes );

                mFieldTypeMap.resize( std::max( static_cast< uint >( tFieldTypes( 0 ) ) + 1, (uint)mFieldTypeMap.size() ), -1 );
                mFieldTypeMap( static_cast< uint >( tFieldTypes( 0 ) ) ) = iFields;

                MORIS_ERROR( ( tFieldParameter.get< std::string >( "field_create_from_file" ).empty() ) or ( tFieldParameter.get< std::string >( "IQI_Name" ).empty() ),
                        "FEM_Model::create_fields(); Field must be either created based on IQI or read from file." );

                MORIS_ERROR( not( ( not tFieldParameter.get< std::string >( "field_create_from_file" ).empty() ) and ( not tFieldParameter.get< std::string >( "IQI_Name" ).empty() ) ),
                        "FEM_Model::create_fields(); Field must be either created based on IQI or read from file." );

                if ( not tFieldParameter.get< std::string >( "field_create_from_file" ).empty() )
                {
                    tField->set_field_from_file(
                            tFieldParameter.get< std::string >( "field_create_from_file" ),
                            tFieldParameter.get< sint >( "field_file_time_index" ),
                            tFieldParameter.get< sint >( "field_file_field_index" ) );
                }

                if ( not tFieldParameter.get< std::string >( "IQI_Name" ).empty() )
                {
                    tField->set_IQI_name( tFieldParameter.get< std::string >( "IQI_Name" ) );
                }

                if ( not tFieldParameter.get< std::string >( "field_output_to_file" ).empty() )
                {
                    tField->set_field_to_file( tFieldParameter.get< std::string >( "field_output_to_file" ) );
                }

                mFields( iFields ) = tField;
            }
        }

        //------------------------------------------------------------------------------

        void FEM_Model::create_material_models(
                std::map< std::string, uint >            &aPropertyMap,
                moris::map< std::string, MSI::Dof_Type > &aMSIDofTypeMap,
                moris::map< std::string, PDV_Type >      &aDvTypeMap )
        {
            // create a constitutive model factory
            MM_Factory tMMFactory;

            // get the MM parameter list
            moris::Cell< ParameterList > tMMParameterList = mParameterList( 8 );

            // get number of constitutive models
            uint tNumMMs = tMMParameterList.size();

            // loop over the parameter lists for MM
            for ( uint iMM = 0; iMM < tNumMMs; iMM++ )
            {
                // get the treated MM parameter list
                ParameterList tMMParameter = tMMParameterList( iMM );

                // get the constitutive type from parameter list
                fem::Material_Type tMMType =
                        static_cast< fem::Material_Type >( tMMParameter.get< uint >( "material_type" ) );

                // get the constitutive model name from parameter list
                std::string tMMName =
                        tMMParameter.get< std::string >( "material_name" );

                // get the phase from parameter list
                std::string tPhaseName =
                        tMMParameter.get< std::string >( "phase_name" );

                // create a constitutive model pointer
                std::shared_ptr< fem::Material_Model > tMM =
                        tMMFactory.create_MM( tMMType );

                // set MM name
                tMM->set_name( tMMName );

                // set MM space dimension
                tMM->set_space_dim( mSpaceDim );

                // set MM dof dependencies
                moris::Cell< moris::Cell< moris::MSI::Dof_Type > > tDofTypes;
                string_to_cell_of_cell(
                        std::get< 0 >( tMMParameter.get< std::pair< std::string, std::string > >( "dof_dependencies" ) ),
                        tDofTypes,
                        aMSIDofTypeMap );
                moris::Cell< std::string > tDofTypeNames;
                string_to_cell(
                        std::get< 1 >( tMMParameter.get< std::pair< std::string, std::string > >( "dof_dependencies" ) ),
                        tDofTypeNames );
                tMM->set_dof_type_list( tDofTypes, tDofTypeNames );

                // set MM properties
                moris::Cell< moris::Cell< std::string > > tPropertyNamesPair;
                string_to_cell_of_cell(
                        tMMParameter.get< std::string >( "properties" ),
                        tPropertyNamesPair );
                for ( uint iProp = 0; iProp < tPropertyNamesPair.size(); iProp++ )
                {
                    // get the property name
                    std::string tPropertyName = tPropertyNamesPair( iProp )( 0 );

                    // check if property in the map
                    MORIS_ERROR( aPropertyMap.find( tPropertyName ) != aPropertyMap.end(),
                            "FEM_Model::create_MMs - Unknown aPropertyString : %s \n",
                            tPropertyName.c_str() );

                    // get property index
                    uint tPropertyIndex = aPropertyMap[ tPropertyName ];

                    // set property for MM
                    tMM->set_property(
                            mProperties( tPropertyIndex ),
                            tPropertyNamesPair( iProp )( 1 ) );
                }

                // set local properties
                tMM->set_local_properties();

                // check the phase exist
                MORIS_ERROR( mPhaseMap.find( tPhaseName ) != mPhaseMap.end(),
                        "FEM_Model::create_material_models_without_phase - Unknown tPhaseName : %s \n",
                        tPhaseName.c_str() );

                // set MM to corresponding phase
                mPhaseInfo( mPhaseMap[ tPhaseName ] ).set_MM( tMM );
            }
        }

        //------------------------------------------------------------------------------

        void
        FEM_Model::create_constitutive_models(
                std::map< std::string, uint >            &aPropertyMap,
                moris::map< std::string, MSI::Dof_Type > &aMSIDofTypeMap,
                moris::map< std::string, PDV_Type >      &aDvTypeMap )
        {
            // create a constitutive model factory
            CM_Factory tCMFactory;

            // get the CM parameter list
            moris::Cell< ParameterList > tCMParameterList = mParameterList( 1 );

            // get number of constitutive models
            uint tNumCMs = tCMParameterList.size();

            // loop over the parameter lists for CM
            for ( uint iCM = 0; iCM < tNumCMs; iCM++ )
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

                // check for unknown phase
                MORIS_ERROR( mPhaseMap.find( tPhaseName ) != mPhaseMap.end(),
                        "FEM_Model::create_CMs - Unknown phase name: %s \n",
                        tPhaseName.c_str() );

                // get the model type
                fem::Model_Type tCMModelType =
                        static_cast< fem::Model_Type >( tCMParameter.get< uint >( "model_type" ) );

                // create a constitutive model pointer
                std::shared_ptr< fem::Constitutive_Model > tCM =
                        tCMFactory.create_CM( tCMType );

                // set CM name
                tCM->set_name( tCMName );

                // set CM model type. must come before "set_space_dim"
                // fixme: currently cannot set a plane type and tensor type at the same time from an input file
                if ( tCMModelType != fem::Model_Type::UNDEFINED )
                {
                    tCM->set_model_type( tCMModelType );
                }

                // set CM space dimension
                tCM->set_space_dim( mSpaceDim );

                // set function parameters
                moris::Cell< moris::Matrix< DDRMat > > tFuncParameters;
                string_to_cell_mat_2(
                        tCMParameter.get< std::string >( "function_parameters" ),
                        tFuncParameters );
                tCM->set_parameters( tFuncParameters );

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

                // set CM material model
                moris::Cell< moris::Cell< std::string > > tMMNamesPair;
                string_to_cell_of_cell(
                        tCMParameter.get< std::string >( "material_model" ),
                        tMMNamesPair );
                MORIS_ERROR( tMMNamesPair.size() <= 1, "FEM_Model::create_CMs() - Only one material model per CM allowed." );

                // loop over Material Model names
                for ( uint iMM = 0; iMM < tMMNamesPair.size(); iMM++ )
                {
                    // get the material name
                    std::string tMaterialName = tMMNamesPair( iMM )( 0 );

                    // get phase index
                    uint tPhaseIndex = mPhaseMap[ tPhaseName ];

                    // get MM from phase
                    std::shared_ptr< fem::Material_Model > tMM =
                            mPhaseInfo( tPhaseIndex ).get_MM_by_name( tMaterialName );

                    // set material for CM
                    tCM->set_material_model(
                            tMM,
                            tMMNamesPair( iMM )( 1 ) );
                }

                // set CM properties
                moris::Cell< moris::Cell< std::string > > tPropertyNamesPair;
                string_to_cell_of_cell(
                        tCMParameter.get< std::string >( "properties" ),
                        tPropertyNamesPair );
                for ( uint iProp = 0; iProp < tPropertyNamesPair.size(); iProp++ )
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
                        "FEM_Model::create_constitutive_models_without_phase - Unknown tPhaseName : %s \n",
                        tPhaseName.c_str() );

                // set CM to corresponding phase
                mPhaseInfo( mPhaseMap[ tPhaseName ] ).set_CM( tCM );
            }
        }

        //------------------------------------------------------------------------------

        std::map< std::string, uint >
        FEM_Model::create_stabilization_parameters( std::map< std::string, uint > &aPropertyMap, moris::map< std::string, MSI::Dof_Type > &aMSIDofTypeMap, moris::map< std::string, PDV_Type > &aDvTypeMap )
        {
            std::map< std::string, uint > tSPMap;

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
                tSPMap[ tSPName ] = iSP;

                // set parameters
                moris::Cell< moris::Matrix< DDRMat > > tFuncParameters;
                string_to_cell_mat_2(
                        tSPParameter.get< std::string >( "function_parameters" ),
                        tFuncParameters );
                mSPs( iSP )->set_parameters( tFuncParameters );

                // init string for leader or follower
                std::string          tIsLeaderString = "leader";
                mtk::Leader_Follower tIsLeader       = mtk::Leader_Follower::LEADER;

                // loop on leader and follower
                for ( uint iLeader = 0; iLeader <= mSPs( iSP )->get_has_follower(); iLeader++ )
                {
                    // if follower
                    if ( iLeader != 0u )
                    {
                        // reset string for follower
                        tIsLeaderString = "follower";
                        tIsLeader       = mtk::Leader_Follower::FOLLOWER;
                    }

                    // get the treated IWG phase
                    std::string tPhaseName =
                            tSPParameter.get< std::string >( tIsLeaderString + "_phase_name" );

                    // check for unknown phase
                    MORIS_ERROR( mPhaseMap.find( tPhaseName ) != mPhaseMap.end(),
                            "FEM_Model::create_stabilization_parameters_without_phase - Unknown phase name: %s \n",
                            tPhaseName.c_str() );

                    // set dof dependencies
                    moris::Cell< moris::Cell< moris::MSI::Dof_Type > > tDofTypes;
                    string_to_cell_of_cell(
                            std::get< 0 >( tSPParameter.get< std::pair< std::string, std::string > >( tIsLeaderString + "_dof_dependencies" ) ),
                            tDofTypes,
                            aMSIDofTypeMap );
                    moris::Cell< std::string > tDofTypeNames;
                    string_to_cell( std::get< 1 >(
                                            tSPParameter.get< std::pair< std::string, std::string > >( tIsLeaderString + "_dof_dependencies" ) ),
                            tDofTypeNames );
                    mSPs( iSP )->set_dof_type_list( tDofTypes, tDofTypeNames, tIsLeader );

                    // set dv dependencies
                    moris::Cell< moris::Cell< PDV_Type > > tDvTypes;
                    string_to_cell_of_cell(
                            std::get< 0 >( tSPParameter.get< std::pair< std::string, std::string > >( tIsLeaderString + "_dv_dependencies" ) ),
                            tDvTypes,
                            aDvTypeMap );
                    moris::Cell< std::string > tDvTypeNames;
                    string_to_cell(
                            std::get< 1 >( tSPParameter.get< std::pair< std::string, std::string > >( tIsLeaderString + "_dv_dependencies" ) ),
                            tDvTypeNames );
                    mSPs( iSP )->set_dv_type_list( tDvTypes, tDvTypeNames, tIsLeader );

                    // set leader properties
                    moris::Cell< moris::Cell< std::string > > tPropertyNamesPair;
                    string_to_cell_of_cell(
                            tSPParameter.get< std::string >( tIsLeaderString + "_properties" ),
                            tPropertyNamesPair );

                    for ( uint iProp = 0; iProp < tPropertyNamesPair.size(); iProp++ )
                    {
                        // get the property name
                        std::string tPropertyName = tPropertyNamesPair( iProp )( 0 );

                        // check for unknown property
                        MORIS_ERROR( aPropertyMap.find( tPropertyName ) != aPropertyMap.end(),
                                "FEM_Model::create_stabilization_parameters_without_phase - Unknown %s aPropertyString : %s \n",
                                tIsLeaderString.c_str(),
                                tPropertyName.c_str() );

                        // get property index
                        uint tPropertyIndex = aPropertyMap[ tPropertyName ];

                        // set property for CM
                        mSPs( iSP )->set_property(
                                mProperties( tPropertyIndex ),
                                tPropertyNamesPair( iProp )( 1 ),
                                tIsLeader );
                    }

                    // set constitutive models
                    moris::Cell< moris::Cell< std::string > > tCMNamesPair;
                    string_to_cell_of_cell(
                            tSPParameter.get< std::string >( tIsLeaderString + "_constitutive_models" ),
                            tCMNamesPair );

                    // loop over CM names
                    for ( uint iCM = 0; iCM < tCMNamesPair.size(); iCM++ )
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

                    // get the cluster measures specifications
                    moris::Cell< moris::Cell< std::string > > tClusterMeasureTypes;
                    string_to_cell_of_cell(
                            std::get< 0 >( tSPParameter.get< std::pair< std::string, std::string > >( "cluster_measures" ) ),
                            tClusterMeasureTypes );

                    // get the cluster measures names
                    moris::Cell< std::string > tClusterMeasureNames;
                    string_to_cell( std::get< 1 >( tSPParameter.get< std::pair< std::string, std::string > >( "cluster_measures" ) ),
                            tClusterMeasureNames );

                    // build a cell of tuples describing the cluster measures specifications
                    moris::Cell< std::tuple<
                            fem::Measure_Type,
                            mtk::Primary_Void,
                            mtk::Leader_Follower > >
                            tClusterMeasureTuples( tClusterMeasureNames.size() );

                    // get fem::Measure_Type, mtk::Primary_Void and mtk::Leader_Follower map
                    // to convert string to enums
                    moris::map< std::string, fem::Measure_Type >    tFemMeasureMap = fem::get_measure_type_map();
                    moris::map< std::string, mtk::Primary_Void >    tMtkPrimaryMap = mtk::get_primary_type_map();
                    moris::map< std::string, mtk::Leader_Follower > tMtkLeaderMap  = mtk::get_leader_type_map();

                    // loop over cluster measures names
                    for ( uint iCMEA = 0; iCMEA < tClusterMeasureNames.size(); iCMEA++ )
                    {
                        // check that measure type is member of map
                        MORIS_ERROR( tFemMeasureMap.key_exists( tClusterMeasureTypes( iCMEA )( 0 ) ),
                                "FEM_Model::create_stabilization_parameters_without_phase - key does not exist: %s",
                                tClusterMeasureTypes( iCMEA )( 0 ).c_str() );

                        // get fem measure type from map
                        fem::Measure_Type tFemMeasureType = tFemMeasureMap.find( tClusterMeasureTypes( iCMEA )( 0 ) );

                        // check that primary type is member of map
                        MORIS_ERROR( tMtkPrimaryMap.key_exists( tClusterMeasureTypes( iCMEA )( 1 ) ),
                                "FEM_Model::create_stabilization_parameters_without_phase - key does not exist: %s",
                                tClusterMeasureTypes( iCMEA )( 1 ).c_str() );

                        // get mtk primary type from map
                        mtk::Primary_Void tMtkPrimaryType = tMtkPrimaryMap.find( tClusterMeasureTypes( iCMEA )( 1 ) );

                        // check that leader type is member of map
                        MORIS_ERROR( tMtkLeaderMap.key_exists( tClusterMeasureTypes( iCMEA )( 2 ) ),
                                "FEM_Model::create_stabilization_parameters_without_phase - key does not exist: %s",
                                tClusterMeasureTypes( iCMEA )( 2 ).c_str() );

                        // get mtk leader type from map
                        mtk::Leader_Follower tMtkLeaderType = tMtkLeaderMap.find( tClusterMeasureTypes( iCMEA )( 2 ) );

                        // build the cluster measure specification tuple and set it in cell of tuples
                        tClusterMeasureTuples( iCMEA ) = std::make_tuple( tFemMeasureType, tMtkPrimaryType, tMtkLeaderType );
                    }

                    // set the cell of cluster measure specification tuples to the SP
                    mSPs( iSP )->set_cluster_measure_type_list(
                            tClusterMeasureTuples,
                            tClusterMeasureNames );
                }
            }
            return tSPMap;
        }

        //------------------------------------------------------------------------------

        void
        FEM_Model::create_IWGs(
                std::map< std::string, uint >            &aPropertyMap,
                std::map< std::string, uint >            &aSPMap,
                moris::map< std::string, MSI::Dof_Type > &aMSIDofTypeMap )
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
            for ( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
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
                moris::Cell< moris::Cell< moris::MSI::Dof_Type > > tResDofTypes;
                string_to_cell_of_cell(
                        tIWGParameter.get< std::string >( "dof_residual" ),
                        tResDofTypes,
                        aMSIDofTypeMap );

                // get function parameters
                moris::Cell< moris::Matrix< DDRMat > > tFuncParameters;
                string_to_cell_mat_2(
                        tIWGParameter.get< std::string >( "function_parameters" ),
                        tFuncParameters );

                // get the treated IWG bulk type
                fem::Element_Type tIWGBulkType =
                        static_cast< fem::Element_Type >( tIWGParameter.get< uint >( "IWG_bulk_type" ) );

                // set flag for leader/follower
                bool tLeaderFollower = ( tIWGBulkType == fem::Element_Type::DOUBLE_SIDESET );

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

                // set constant parameters
                mIWGs( iIWG )->set_parameters( tFuncParameters );

                // initialize string for leader or follower
                std::string          tIsLeaderString = "leader";
                mtk::Leader_Follower tIsLeader       = mtk::Leader_Follower::LEADER;

                // loop on leader and follower
                for ( uint iLeader = 0; iLeader <= tLeaderFollower; iLeader++ )
                {
                    // if follower
                    if ( iLeader )
                    {
                        // reset string for follower
                        tIsLeaderString = "follower";
                        tIsLeader       = mtk::Leader_Follower::FOLLOWER;
                    }

                    // get the treated IWG phase
                    std::string tPhaseName =
                            tIWGParameter.get< std::string >( tIsLeaderString + "_phase_name" );

                    // check for unknown phase
                    MORIS_ERROR( mPhaseMap.find( tPhaseName ) != mPhaseMap.end(),
                            "FEM_Model::create_IWGs_without_phase - Unknown phase name: %s \n",
                            tPhaseName.c_str() );

                    // set phase name
                    mIWGs( iIWG )->set_phase_name( tPhaseName, tIsLeader );

                    // get phase index
                    uint tPhaseIndex = mPhaseMap[ tPhaseName ];

                    // get dof type list from phase - ignore double-sided side sets
                    if ( !tLeaderFollower )
                    {
                        mPhaseInfo( tPhaseIndex ).add_dof_type_to_list( tResDofTypes );
                    }

                    // set properties
                    moris::Cell< moris::Cell< std::string > > tPropertyNamesPair;
                    string_to_cell_of_cell(
                            tIWGParameter.get< std::string >( tIsLeaderString + "_properties" ),
                            tPropertyNamesPair );

                    for ( uint iProp = 0; iProp < tPropertyNamesPair.size(); iProp++ )
                    {
                        // get property name
                        std::string tPropertyName = tPropertyNamesPair( iProp )( 0 );

                        // check for unknown property
                        MORIS_ERROR( aPropertyMap.find( tPropertyName ) != aPropertyMap.end(),
                                "FEM_Model::create_IWGs_without_phase - Unknown %s aPropertyString: %s \n",
                                tIsLeaderString.c_str(),
                                tPropertyName.c_str() );

                        // get property index
                        uint tPropertyIndex = aPropertyMap[ tPropertyName ];

                        // set property for IWG
                        mIWGs( iIWG )->set_property(
                                mProperties( tPropertyIndex ),
                                tPropertyNamesPair( iProp )( 1 ),
                                tIsLeader );
                    }

                    // set material model
                    moris::Cell< moris::Cell< std::string > > tMMNamesPair;
                    string_to_cell_of_cell(
                            tIWGParameter.get< std::string >( tIsLeaderString + "_material_model" ),
                            tMMNamesPair );
                    MORIS_ERROR( tMMNamesPair.size() <= 1, "FEM_Model::create_IWGs_without_phase() - Only one material model per CM allowed." );

                    // loop over material model
                    for ( uint iMM = 0; iMM < tMMNamesPair.size(); iMM++ )
                    {
                        // get the MM name
                        std::string tMMName = tMMNamesPair( iMM )( 0 );

                        // get MM from phase
                        std::shared_ptr< fem::Material_Model > tMM =
                                mPhaseInfo( tPhaseIndex ).get_MM_by_name( tMMName );

                        // set CM for IWG
                        mIWGs( iIWG )->set_material_model(
                                tMM,
                                tMMNamesPair( iMM )( 1 ),
                                tIsLeader );
                    }

                    // set constitutive models
                    moris::Cell< moris::Cell< std::string > > tCMNamesPair;
                    string_to_cell_of_cell(
                            tIWGParameter.get< std::string >( tIsLeaderString + "_constitutive_models" ),
                            tCMNamesPair );

                    // loop over constitutive models
                    for ( uint iCM = 0; iCM < tCMNamesPair.size(); iCM++ )
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
                                tIsLeader );
                    }
                }

                // set stabilization parameters
                moris::Cell< moris::Cell< std::string > > tSPNamesPair;
                string_to_cell_of_cell(
                        tIWGParameter.get< std::string >( "stabilization_parameters" ),
                        tSPNamesPair );

                // loop over SP names
                for ( uint iSP = 0; iSP < tSPNamesPair.size(); iSP++ )
                {
                    // get the SP name
                    std::string tSPName = tSPNamesPair( iSP )( 0 );

                    // check for unknown SP
                    MORIS_ERROR( aSPMap.find( tSPName ) != aSPMap.end(),
                            "FEM_Model::create_IWGs_without_phase - Unknown aSPString: %s \n",
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
            for ( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
            {
                // get the treated IWG parameter list
                ParameterList tIWGParameter = tIWGParameterList( iIWG );

                // get the IWG bulk type
                fem::Element_Type tIWGBulkType = mIWGs( iIWG )->get_bulk_type();

                // get the IWG leader phase name
                std::string tPhaseName =
                        mIWGs( iIWG )->get_phase_name( mtk::Leader_Follower::LEADER );

                // check for unknown phase
                MORIS_ERROR( mPhaseMap.find( tPhaseName ) != mPhaseMap.end(),
                        "FEM_Model::create_IWGs_without_phase - Unknown phase name: %s \n",
                        tPhaseName.c_str() );

                // get the phase index
                uint tPhaseIndex = mPhaseMap[ tPhaseName ];

                // get dof type list from leader phase
                const moris::Cell< moris::Cell< MSI::Dof_Type > > &tLeaderDofTypes =
                        mPhaseInfo( tPhaseIndex ).get_dof_type_list();

                // get dof type list from leader phase
                const moris::Cell< moris::Cell< PDV_Type > > &tLeaderPdvTypes =
                        mPhaseInfo( tPhaseIndex ).get_dv_type_list();

                // set leader dof dependencies
                mIWGs( iIWG )->set_dof_type_list( tLeaderDofTypes, mtk::Leader_Follower::LEADER );

                // set leader dv dependencies
                mIWGs( iIWG )->set_dv_type_list( tLeaderPdvTypes, mtk::Leader_Follower::LEADER );

                if ( tIWGBulkType == fem::Element_Type::DOUBLE_SIDESET )
                {
                    // get the IWG follower phase name
                    std::string tFollowerPhaseName =
                            mIWGs( iIWG )->get_phase_name( mtk::Leader_Follower::FOLLOWER );

                    // check for unknown phase
                    MORIS_ERROR( mPhaseMap.find( tFollowerPhaseName ) != mPhaseMap.end(),
                            "FEM_Model::create_IWGs_without_phase - Unknown phase name: %s \n",
                            tFollowerPhaseName.c_str() );

                    // get CM index
                    uint tFollowerPhaseIndex = mPhaseMap[ tFollowerPhaseName ];

                    // get dof type list from phase
                    const moris::Cell< moris::Cell< MSI::Dof_Type > > &tFollowerDofTypes =
                            mPhaseInfo( tFollowerPhaseIndex ).get_dof_type_list();

                    // get pdv type list from phase
                    const moris::Cell< moris::Cell< PDV_Type > > &tFollowerPdvTypes =
                            mPhaseInfo( tFollowerPhaseIndex ).get_dv_type_list();

                    // set follower dof dependencies
                    mIWGs( iIWG )->set_dof_type_list( tFollowerDofTypes, mtk::Leader_Follower::FOLLOWER );

                    // set follower dv dependencies
                    mIWGs( iIWG )->set_dv_type_list( tFollowerPdvTypes, mtk::Leader_Follower::FOLLOWER );
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        FEM_Model::create_IQIs(
                std::map< std::string, uint >            &aPropertyMap,
                std::map< std::string, uint >            &aSPMap,
                moris::map< std::string, MSI::Dof_Type > &aMSIDofTypeMap )
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
            for ( uint iIQI = 0; iIQI < tNumIQIs; iIQI++ )
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
                bool tLeaderFollower = ( tIQIBulkType == fem::Element_Type::DOUBLE_SIDESET );

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

                // init string for leader or follower
                std::string          tIsLeaderString = "leader";
                mtk::Leader_Follower tIsLeader       = mtk::Leader_Follower::LEADER;

                // loop on leader and follower
                for ( uint iLeader = 0; iLeader <= tLeaderFollower; iLeader++ )
                {
                    // if follower
                    if ( iLeader )
                    {
                        // reset string for follower
                        tIsLeaderString = "follower";
                        tIsLeader       = mtk::Leader_Follower::FOLLOWER;
                    }

                    // get the treated IWG phase
                    std::string tPhaseName = tIQIParameter.get< std::string >( tIsLeaderString + "_phase_name" );

                    // check for unknown phase
                    MORIS_ERROR( mPhaseMap.find( tPhaseName ) != mPhaseMap.end(),
                            "FEM_Model::create_IQIs_without_phase - %s: Unknown %s phase name: %s.",
                            tIQIName.c_str(),
                            tIsLeaderString.c_str(),
                            tPhaseName.c_str() );

                    // set phase name
                    mIQIs( iIQI )->set_phase_name( tPhaseName, tIsLeader );

                    // get the phase index
                    uint tPhaseIndex = mPhaseMap[ tPhaseName ];

                    // get dof type list from phase
                    const moris::Cell< moris::Cell< MSI::Dof_Type > > &tDofTypes =
                            mPhaseInfo( tPhaseIndex ).get_dof_type_list();

                    // get dof type list from phase
                    const moris::Cell< moris::Cell< PDV_Type > > &tDvTypes =
                            mPhaseInfo( tPhaseIndex ).get_dv_type_list();

                    // set leader dof dependencies
                    mIQIs( iIQI )->set_dof_type_list( tDofTypes );

                    // set leader dv dependencies
                    mIQIs( iIQI )->set_dv_type_list( tDvTypes );

                    // set leader properties
                    moris::Cell< moris::Cell< std::string > > tPropertyNamesPair;
                    string_to_cell_of_cell(
                            tIQIParameter.get< std::string >( tIsLeaderString + "_properties" ),
                            tPropertyNamesPair );

                    for ( uint iProp = 0; iProp < tPropertyNamesPair.size(); iProp++ )
                    {
                        // get the property name
                        std::string tPropertyName = tPropertyNamesPair( iProp )( 0 );

                        // check for unknown property
                        MORIS_ERROR( aPropertyMap.find( tPropertyName ) != aPropertyMap.end(),
                                "FEM_Model::create_IQIs_without_phase - Unknown %s aPropertyString: %s \n",
                                tIsLeaderString.c_str(),
                                tPropertyName.c_str() );

                        // get property index
                        uint tPropertyIndex = aPropertyMap[ tPropertyName ];

                        // set property for IWG
                        mIQIs( iIQI )->set_property(
                                mProperties( tPropertyIndex ),
                                tPropertyNamesPair( iProp )( 1 ),
                                tIsLeader );
                    }

                    // set leader constitutive models
                    moris::Cell< moris::Cell< std::string > > tCMNamesPair;
                    string_to_cell_of_cell(
                            tIQIParameter.get< std::string >( tIsLeaderString + "_constitutive_models" ),
                            tCMNamesPair );

                    for ( uint iCM = 0; iCM < tCMNamesPair.size(); iCM++ )
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
                                tIsLeader );
                    }
                }

                // set stabilization parameters
                moris::Cell< moris::Cell< std::string > > tSPNamesPair;
                string_to_cell_of_cell(
                        tIQIParameter.get< std::string >( "stabilization_parameters" ),
                        tSPNamesPair );

                for ( uint iSP = 0; iSP < tSPNamesPair.size(); iSP++ )
                {
                    // get the SP name
                    std::string tSPName = tSPNamesPair( iSP )( 0 );

                    // check for unknown SP
                    MORIS_ERROR( aSPMap.find( tSPName ) != aSPMap.end(),
                            "FEM_Model::create_IQIs_without_phase - Unknown aSPString: %s \n",
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

        void
        FEM_Model::create_phases()
        {
            // get the phase parameter list
            moris::Cell< ParameterList > tPhaseParameterList = mParameterList( 7 );

            // get number of phases
            uint tNumPhases = tPhaseParameterList.size();

            // resize the list of phase user info
            mPhaseInfo.resize( tNumPhases );

            // loop over the parameter lists
            for ( uint iPhase = 0; iPhase < tNumPhases; iPhase++ )
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

        void
        FEM_Model::create_fem_set_info()
        {
            // init number of fem sets to be created
            uint tNumFEMSets = 0;

            // get the IWG and IQI parameter lists
            moris::Cell< ParameterList > tIWGParameterList = mParameterList( 3 );
            moris::Cell< ParameterList > tIQIParameterList = mParameterList( 4 );

            // get fem computation type parameter list
            ParameterList tComputationParameterList = mParameterList( 5 )( 0 );

            // bool true for analytical forward analysis, false for finite difference
            // decide if dRdu and dQIdu are computed by A/FD
            bool tIsAnalyticalFA =
                    tComputationParameterList.get< bool >( "is_analytical_forward" );

            // get enum for FD scheme for forward analysis
            fem::FDScheme_Type tFDSchemeForFA = static_cast< fem::FDScheme_Type >(
                    tComputationParameterList.get< uint >( "finite_difference_scheme_forward" ) );

            // get perturbation size for FD for forward analysis
            real tFDPerturbationFA = tComputationParameterList.get< real >(
                    "finite_difference_perturbation_size_forward" );

            // get bool for analytical/finite difference for sensitivity analysis
            // decide if dRdp and dQIdp are computed by A/FD
            bool tIsAnalyticalSA =
                    tComputationParameterList.get< bool >( "is_analytical_sensitivity" );

            // get enum for FD scheme for sensitivity analysis
            fem::FDScheme_Type tFDSchemeForSA = static_cast< fem::FDScheme_Type >(
                    tComputationParameterList.get< uint >( "finite_difference_scheme" ) );

            // get perturbation size for FD for sensitivity analysis
            real tFDPerturbationSA = tComputationParameterList.get< real >(
                    "finite_difference_perturbation_size" );

            // get enum for perturbation strategy for finite difference
            fem::Perturbation_Type tPerturbationStrategy = static_cast< fem::Perturbation_Type >(
                    tComputationParameterList.get< uint >( "finite_difference_perturbation_strategy" ) );

            // create a map of the set
            std::map< std::tuple< std::string, bool, bool >, uint > tMeshToFemSet;

            // loop over the IWGs
            for ( uint iIWG = 0; iIWG < tIWGParameterList.size(); iIWG++ )
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

                // get the IWG leader phase name
                std::string tLeaderPhaseName =
                        mIWGs( iIWG )->get_phase_name( mtk::Leader_Follower::LEADER );

                // get the IWG follower phase name
                std::string tFollowerPhaseName =
                        mIWGs( iIWG )->get_phase_name( mtk::Leader_Follower::FOLLOWER );

                // get follower phase string from IWG input
                std::string tFollowerPhaseString =
                        tIWGParameter.get< std::string >( "neighbor_phases" );

                // get ordinal string from IWG input
                std::string tOrdinalString =
                        tIWGParameter.get< std::string >( "side_ordinals" );

                // get mesh set names for IWG
                moris::Cell< std::string > tMeshSetNames;
                this->get_mesh_set_names(
                        tIWGBulkType,
                        tLeaderPhaseName,
                        tFollowerPhaseName,
                        tFollowerPhaseString,
                        tOrdinalString,
                        tIsGhost,
                        tMeshSetNames );

                // get a representative DoF type
                MSI::Dof_Type tFirstResidualDofType = mIWGs( iIWG )->get_residual_dof_type()( 0 )( 0 );

                // loop over the mesh set names
                for ( uint iSetName = 0; iSetName < tMeshSetNames.size(); iSetName++ )
                {
                    // get the name of the set currently treated
                    std::string tMeshSetName = tMeshSetNames( iSetName );

                    // check for ghost set names and select correct B-spline mesh automatically when new ghost sets need to be used
                    this->check_and_set_ghost_set_names( tMeshSetName, tFirstResidualDofType );

                    // check if the mesh set name already in map
                    if ( tMeshToFemSet.find( std::make_tuple(
                                 tMeshSetName,
                                 tTimeContinuity,
                                 tTimeBoundary ) )
                            == tMeshToFemSet.end() )
                    {
                        // add the mesh set name map
                        tMeshToFemSet[ std::make_tuple(
                                tMeshSetName,
                                tTimeContinuity,
                                tTimeBoundary ) ] = tNumFEMSets++;

                        // create a fem set info for the mesh set
                        Set_User_Info aSetUserInfo;

                        // set its mesh set name
                        aSetUserInfo.set_mesh_set_name( tMeshSetName );

                        // set its time continuity flag
                        aSetUserInfo.set_time_continuity( tTimeContinuity );

                        // set its time boundary flag
                        aSetUserInfo.set_time_boundary( tTimeBoundary );

                        // set its forward analysis type flag
                        aSetUserInfo.set_is_analytical_forward_analysis( tIsAnalyticalFA );

                        // set its FD scheme for forward analysis
                        aSetUserInfo.set_finite_difference_scheme_for_forward_analysis( tFDSchemeForFA );

                        // set its FD perturbation size for forward analysis
                        aSetUserInfo.set_finite_difference_perturbation_size_for_forward_analysis( tFDPerturbationFA );

                        // set its sensitivity analysis type flag
                        aSetUserInfo.set_is_analytical_sensitivity_analysis( tIsAnalyticalSA );

                        // set its FD scheme for sensitivity analysis
                        aSetUserInfo.set_finite_difference_scheme_for_sensitivity_analysis( tFDSchemeForSA );

                        // set its FD perturbation size for sensitivity analysis
                        aSetUserInfo.set_finite_difference_perturbation_size( tFDPerturbationSA );

                        // set its perturbation strategy for finite difference
                        aSetUserInfo.set_perturbation_strategy( tPerturbationStrategy );

                        // set the IWG
                        aSetUserInfo.add_IWG( mIWGs( iIWG ) );

                        // add it to the list of fem set info
                        mSetInfo.push_back( aSetUserInfo );
                    }
                    else
                    {
                        // set the IWG
                        mSetInfo( tMeshToFemSet[ std::make_tuple(
                                          tMeshSetName,
                                          tTimeContinuity,
                                          tTimeBoundary ) ] )
                                .add_IWG( mIWGs( iIWG ) );
                    }
                }
            }

            // loop over the IQIs
            for ( uint iIQI = 0; iIQI < tIQIParameterList.size(); iIQI++ )
            {
                // get the treated IWG parameter list
                ParameterList tIQIParameter = tIQIParameterList( iIQI );

                // get the IWG bulk type
                fem::Element_Type tIQIBulkType = mIQIs( iIQI )->get_bulk_type();

                // get time continuity flag
                bool tTimeContinuity = mIQIs( iIQI )->get_time_continuity();

                // get time boundary flag
                bool tTimeBoundary = mIQIs( iIQI )->get_time_boundary();

                // get the IWG leader phase name
                std::string tLeaderPhaseName =
                        mIQIs( iIQI )->get_phase_name( mtk::Leader_Follower::LEADER );

                // get the IWG follower phase name
                std::string tFollowerPhaseName =
                        mIQIs( iIQI )->get_phase_name( mtk::Leader_Follower::FOLLOWER );

                // get follower phase string from IQI input
                std::string tFollowerPhaseString =
                        tIQIParameter.get< std::string >( "neighbor_phases" );

                // get ordinal string from IQI input
                std::string tOrdinalString =
                        tIQIParameter.get< std::string >( "side_ordinals" );

                // get mesh set names for IWG
                moris::Cell< std::string > tMeshSetNames;
                this->get_mesh_set_names(
                        tIQIBulkType,
                        tLeaderPhaseName,
                        tFollowerPhaseName,
                        tFollowerPhaseString,
                        tOrdinalString,
                        false,
                        tMeshSetNames );

                // loop over the mesh set names
                for ( uint iSetName = 0; iSetName < tMeshSetNames.size(); iSetName++ )
                {
                    // get the name of the set currently treated
                    std::string tMeshSetName = tMeshSetNames( iSetName );

                    // if the mesh set name not in map
                    if ( tMeshToFemSet.find( std::make_tuple(
                                 tMeshSetName,
                                 tTimeContinuity,
                                 tTimeBoundary ) )
                            == tMeshToFemSet.end() )
                    {
                        // add the mesh set name map
                        tMeshToFemSet[ std::make_tuple(
                                tMeshSetName,
                                tTimeContinuity,
                                tTimeBoundary ) ] = tNumFEMSets++;

                        // create a fem set info for the mesh set
                        Set_User_Info aSetUserInfo;

                        // set its mesh set name
                        aSetUserInfo.set_mesh_set_name( tMeshSetName );

                        // set its time continuity flag
                        aSetUserInfo.set_time_continuity( tTimeContinuity );

                        // set its time boundary flag
                        aSetUserInfo.set_time_boundary( tTimeBoundary );

                        // set its forward analysis type flag
                        aSetUserInfo.set_is_analytical_forward_analysis( tIsAnalyticalFA );

                        // set its FD scheme for forward analysis
                        aSetUserInfo.set_finite_difference_scheme_for_forward_analysis( tFDSchemeForFA );

                        // set its FD perturbation size for forward analysis
                        aSetUserInfo.set_finite_difference_perturbation_size_for_forward_analysis( tFDPerturbationFA );

                        // set its sensitivity analysis type flag
                        aSetUserInfo.set_is_analytical_sensitivity_analysis( tIsAnalyticalSA );

                        // set its FD scheme for sensitivity analysis
                        aSetUserInfo.set_finite_difference_scheme_for_sensitivity_analysis( tFDSchemeForSA );

                        // set its FD perturbation size for sensitivity analysis
                        aSetUserInfo.set_finite_difference_perturbation_size( tFDPerturbationSA );

                        // set its perturbation strategy for finite difference
                        aSetUserInfo.set_perturbation_strategy( tPerturbationStrategy );

                        // set the IQI
                        aSetUserInfo.add_IQI( mIQIs( iIQI ) );

                        // add it to the list of fem set info
                        mSetInfo.push_back( aSetUserInfo );
                    }
                    else
                    {
                        // set the IQI
                        mSetInfo( tMeshToFemSet[ std::make_tuple(
                                          tMeshSetName,
                                          tTimeContinuity,
                                          tTimeBoundary ) ] )
                                .add_IQI( mIQIs( iIQI ) );
                    }
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        FEM_Model::get_mesh_set_names(
                fem::Element_Type           aIWGBulkType,
                const std::string          &aLeaderPhaseName,
                const std::string          &aFollowerPhaseName,
                const std::string          &aFollowerPhaseString,
                const std::string          &aOrdinalString,
                bool                        aIsGhost,
                moris::Cell< std::string > &aMeshSetNames )
        {
            // get the leader phase mesh index
            moris::Matrix< moris::IndexMat > tLeaderPhaseIndices =
                    mPhaseInfo( mPhaseMap[ aLeaderPhaseName ] ).get_phase_indices();

            // get the number of leader phase mesh indices
            uint tNumLeaderIndices = tLeaderPhaseIndices.numel();

            // switch on the element type
            switch ( aIWGBulkType )
            {
                case fem::Element_Type::BULK:
                {
                    // loop over phase mesh indices
                    for ( uint iMeshIndex = 0; iMeshIndex < tNumLeaderIndices; iMeshIndex++ )
                    {
                        //                        // FIXME ! get mesh set names from integration mesh for index
                        //                        mMeshManager->get_integration_mesh( 0 )->
                        //                                get_block_set_names_with_color( tLeaderPhaseIndices( iMeshIndex ), aMeshSetNames );

                        // add mesh set name to list
                        aMeshSetNames.push_back(
                                "HMR_dummy_c_p" + std::to_string( tLeaderPhaseIndices( iMeshIndex ) ) );

                        // add mesh set name to list
                        aMeshSetNames.push_back(
                                "HMR_dummy_n_p" + std::to_string( tLeaderPhaseIndices( iMeshIndex ) ) );
                    }
                    break;
                }
                case fem::Element_Type::SIDESET:
                {
                    // get neighbor phase names from string
                    moris::Cell< std::string > tFollowerPhaseNames;
                    string_to_cell( aFollowerPhaseString, tFollowerPhaseNames );

                    // get number of neighbor phase
                    uint tNumSingle = tFollowerPhaseNames.size();

                    // get ordinals for boundary from string
                    Matrix< DDSMat > tOrdinals;
                    string_to_mat( aOrdinalString, tOrdinals );
                    uint tNumBoundary = tOrdinals.numel();

                    // loop over leader phase mesh indices
                    for ( uint iLeaderMeshIndex = 0; iLeaderMeshIndex < tNumLeaderIndices; iLeaderMeshIndex++ )
                    {
                        // get single sideset
                        for ( uint iSingle = 0; iSingle < tNumSingle; iSingle++ )
                        {
                            // get the neighbor phase name
                            std::string tNeighborPhaseName = tFollowerPhaseNames( iSingle );

                            // get the follower phase mesh index
                            moris::Matrix< moris::IndexMat > tFollowerPhaseIndices =
                                    mPhaseInfo( mPhaseMap[ tNeighborPhaseName ] ).get_phase_indices();

                            // get number of neighbor phase mesh indices
                            uint tNumNeighborIndices = tFollowerPhaseIndices.numel();

                            for ( uint iNeighborMeshIndex = 0; iNeighborMeshIndex < tNumNeighborIndices; iNeighborMeshIndex++ )
                            {
                                // FIXME get this info from the mesh
                                // add mesh set name to list
                                aMeshSetNames.push_back(
                                        "iside_b0_" + std::to_string( tLeaderPhaseIndices( iLeaderMeshIndex ) ) + "_b1_" + std::to_string( tFollowerPhaseIndices( iNeighborMeshIndex ) ) );
                            }
                        }

                        // get boundary sideset
                        for ( uint iBoundary = 0; iBoundary < tNumBoundary; iBoundary++ )
                        {
                            // FIXME get this info from the mesh
                            // add mesh set name to list
                            aMeshSetNames.push_back(
                                    "SideSet_" + std::to_string( tOrdinals( iBoundary ) ) + "_c_p" + std::to_string( tLeaderPhaseIndices( iLeaderMeshIndex ) ) );

                            // FIXME get this info from the mesh
                            // add mesh set name to list
                            aMeshSetNames.push_back(
                                    "SideSet_" + std::to_string( tOrdinals( iBoundary ) ) + "_n_p" + std::to_string( tLeaderPhaseIndices( iLeaderMeshIndex ) ) );
                        }
                    }
                    break;
                }
                case fem::Element_Type::DOUBLE_SIDESET:
                {
                    // if ghost
                    if ( aIsGhost )
                    {
                        // loop over leader phase mesh indices
                        for ( uint iLeaderMeshIndex = 0; iLeaderMeshIndex < tNumLeaderIndices; iLeaderMeshIndex++ )
                        {
                            // FIXME get this info from the mesh
                            // add mesh set name to list
                            aMeshSetNames.push_back(
                                    "ghost_p" + std::to_string( tLeaderPhaseIndices( iLeaderMeshIndex ) ) );
                        }
                    }
                    // if interface
                    else
                    {
                        MORIS_ERROR( aLeaderPhaseName != aFollowerPhaseName,
                                "FEM_Model::get_mesh_set_names - Leader and follower phases are the same, FIXME case not handled yet " );

                        // get the follower phase mesh index
                        moris::Matrix< moris::IndexMat > tFollowerPhaseIndices =
                                mPhaseInfo( mPhaseMap[ aFollowerPhaseName ] ).get_phase_indices();

                        // get number of follower phase mesh index
                        uint tNumFollowerIndices = tFollowerPhaseIndices.numel();

                        // loop over leader phase mesh indices
                        for ( uint iLeaderMeshIndex = 0; iLeaderMeshIndex < tNumLeaderIndices; iLeaderMeshIndex++ )
                        {
                            // get leader index
                            uint tLeaderPhaseIndex = tLeaderPhaseIndices( iLeaderMeshIndex );

                            // loop over follower phase mesh indices
                            for ( uint iFollowerMeshIndex = 0; iFollowerMeshIndex < tNumFollowerIndices; iFollowerMeshIndex++ )
                            {
                                // get follower index
                                uint tFollowerPhaseIndex = tFollowerPhaseIndices( iFollowerMeshIndex );

                                // if leader and follower index are different
                                if ( tLeaderPhaseIndex != tFollowerPhaseIndex )
                                {
                                    // FIXME get this info from the mesh
                                    // get interface name
                                    aMeshSetNames.push_back(
                                            "dbl_iside_p0_" + std::to_string( tLeaderPhaseIndices( iLeaderMeshIndex ) ) + "_p1_" + std::to_string( tFollowerPhaseIndices( iFollowerMeshIndex ) ) );
                                }
                            }
                        }
                    }
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "FEM_Model::get_mesh_set_names - Unknown set type" );
                }
            }
        }

        //-------------------------------------------------------------------------------------------------

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

        //-------------------------------------------------------------------------------------------------
        // FEM INPUT - old version
        //-------------------------------------------------------------------------------------------------

        std::map< std::string, uint >
        FEM_Model::create_material_models_without_phase( std::map< std::string, uint > &aPropertyMap, moris::map< std::string, MSI::Dof_Type > &aMSIDofTypeMap, moris::map< std::string, PDV_Type > &aDvTypeMap )
        {
            std::map< std::string, uint > tMMMap;

            // create a material model factory
            MM_Factory tMMFactory;

            // get the MM parameter list
            moris::Cell< ParameterList > tMMParameterList = mParameterList( 7 );

            // get number of material models
            uint tNumMMs = tMMParameterList.size();

            // create a list of MMs
            mMMs.resize( tNumMMs, nullptr );

            // loop over the parameter lists for MM
            for ( uint iMM = 0; iMM < tNumMMs; iMM++ )
            {
                // get the material type from parameter list
                fem::Material_Type tMMType =
                        static_cast< fem::Material_Type >( tMMParameterList( iMM ).get< uint >( "material_type" ) );

                // create a material model pointer
                mMMs( iMM ) = tMMFactory.create_MM( tMMType );

                // set MM name
                mMMs( iMM )->set_name( tMMParameterList( iMM ).get< std::string >( "material_name" ) );

                // fill MM map
                tMMMap[ tMMParameterList( iMM ).get< std::string >( "material_name" ) ] = iMM;

                // set MM space dimension
                mMMs( iMM )->set_space_dim( mSpaceDim );

                // set MM dof dependencies
                moris::Cell< moris::Cell< moris::MSI::Dof_Type > > tDofTypes;
                string_to_cell_of_cell(
                        std::get< 0 >( tMMParameterList( iMM ).get< std::pair< std::string, std::string > >( "dof_dependencies" ) ),
                        tDofTypes,
                        aMSIDofTypeMap );
                moris::Cell< std::string > tDofTypeNames;
                string_to_cell(
                        std::get< 1 >( tMMParameterList( iMM ).get< std::pair< std::string, std::string > >( "dof_dependencies" ) ),
                        tDofTypeNames );
                mMMs( iMM )->set_dof_type_list( tDofTypes, tDofTypeNames );

                // set MM properties
                moris::Cell< moris::Cell< std::string > > tPropertyNamesPair;
                string_to_cell_of_cell(
                        tMMParameterList( iMM ).get< std::string >( "properties" ),
                        tPropertyNamesPair );
                for ( uint iProp = 0; iProp < tPropertyNamesPair.size(); iProp++ )
                {
                    // if property name is in the property map
                    if ( aPropertyMap.find( tPropertyNamesPair( iProp )( 0 ) ) != aPropertyMap.end() )
                    {
                        // get property index
                        uint tPropertyIndex = aPropertyMap[ tPropertyNamesPair( iProp )( 0 ) ];

                        // set property for MM
                        mMMs( iMM )->set_property(
                                mProperties( tPropertyIndex ),
                                tPropertyNamesPair( iProp )( 1 ) );
                    }
                    else
                    {
                        // error message for unknown property
                        MORIS_ERROR( false,
                                "FEM_Model::create_MMs - Unknown aPropertyString : %s \n",
                                tPropertyNamesPair( iProp )( 0 ).c_str() );
                    }
                }
                // set local properties
                mMMs( iMM )->set_local_properties();
            }
            return tMMMap;
        }

        //-------------------------------------------------------------------------------------------------

        std::map< std::string, uint >
        FEM_Model::create_constitutive_models_without_phase( std::map< std::string, uint > &aPropertyMap, std::map< std::string, uint > &aMMMap, moris::map< std::string, MSI::Dof_Type > &aMSIDofTypeMap, moris::map< std::string, PDV_Type > &aDvTypeMap )
        {
            std::map< std::string, uint > tCMMap;

            // create a constitutive model factory
            CM_Factory tCMFactory;

            // get the CM parameter list
            moris::Cell< ParameterList > tCMParameterList = mParameterList( 1 );

            // get number of constitutive models
            uint tNumCMs = tCMParameterList.size();

            // create a list of CMs
            mCMs.resize( tNumCMs, nullptr );

            // loop over the parameter lists for CM
            for ( uint iCM = 0; iCM < tNumCMs; iCM++ )
            {
                // get the constitutive type from parameter list
                fem::Constitutive_Type tCMType =
                        static_cast< fem::Constitutive_Type >( tCMParameterList( iCM ).get< uint >( "constitutive_type" ) );

                // create a constitutive model pointer
                mCMs( iCM ) = tCMFactory.create_CM( tCMType );

                // set CM name
                mCMs( iCM )->set_name( tCMParameterList( iCM ).get< std::string >( "constitutive_name" ) );

                // fill CM map
                tCMMap[ tCMParameterList( iCM ).get< std::string >( "constitutive_name" ) ] = iCM;

                // set CM model type
                fem::Model_Type tCMModelType =
                        static_cast< fem::Model_Type >( tCMParameterList( iCM ).get< uint >( "model_type" ) );
                if ( tCMModelType != fem::Model_Type::UNDEFINED )
                {
                    mCMs( iCM )->set_model_type( tCMModelType );
                }

                // set CM space dimension
                mCMs( iCM )->set_space_dim( mSpaceDim );

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
                for ( uint iProp = 0; iProp < tPropertyNamesPair.size(); iProp++ )
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

                // set material model
                moris::Cell< moris::Cell< std::string > > tMMNamesPair;
                string_to_cell_of_cell(
                        tCMParameterList( iCM ).get< std::string >( "material_model" ),
                        tMMNamesPair );

                for ( uint iMM = 0; iMM < tMMNamesPair.size(); iMM++ )
                {
                    // if MM name is in the CM map
                    if ( aMMMap.find( tMMNamesPair( iMM )( 0 ) ) != aMMMap.end() )
                    {
                        // get MM index
                        uint tMMIndex = aMMMap[ tMMNamesPair( iMM )( 0 ) ];

                        // set MM for IWG
                        mCMs( iCM )->set_material_model(
                                mMMs( tMMIndex ),
                                tMMNamesPair( iMM )( 1 ) );
                    }
                    else
                    {
                        // error message unknown MM
                        MORIS_ERROR( false,
                                "FEM_Model::create_CMs - Unknown aMMString: %s \n",
                                tMMNamesPair( iMM )( 0 ).c_str() );
                    }
                }
            }
            return tCMMap;
        }

        //------------------------------------------------------------------------------

        std::map< std::string, uint >
        FEM_Model::create_stabilization_parameters_without_phase( std::map< std::string, uint > &aPropertyMap, std::map< std::string, uint > &aCMMap, moris::map< std::string, MSI::Dof_Type > &aMSIDofTypeMap, moris::map< std::string, PDV_Type > &aDvTypeMap )
        {
            std::map< std::string, uint > tSPMap;

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
                tSPMap[ tSPParameter.get< std::string >( "stabilization_name" ) ] = iSP;

                // set parameters
                moris::Cell< moris::Matrix< DDRMat > > tFuncParameters;
                string_to_cell_mat_2(
                        tSPParameter.get< std::string >( "function_parameters" ),
                        tFuncParameters );
                mSPs( iSP )->set_parameters( tFuncParameters );

                // init string for leader or follower
                std::string          tIsLeaderString = "leader";
                mtk::Leader_Follower tIsLeader       = mtk::Leader_Follower::LEADER;

                // loop on leader and follower
                for ( uint iLeader = 0; iLeader <= mSPs( iSP )->get_has_follower(); iLeader++ )
                {
                    // if follower
                    if ( iLeader )
                    {
                        // reset string for follower
                        tIsLeaderString = "follower";
                        tIsLeader       = mtk::Leader_Follower::FOLLOWER;
                    }

                    // set dof dependencies
                    moris::Cell< moris::Cell< moris::MSI::Dof_Type > > tDofTypes;
                    string_to_cell_of_cell(
                            std::get< 0 >( tSPParameter.get< std::pair< std::string, std::string > >( tIsLeaderString + "_dof_dependencies" ) ),
                            tDofTypes,
                            aMSIDofTypeMap );
                    moris::Cell< std::string > tDofTypeNames;
                    string_to_cell( std::get< 1 >(
                                            tSPParameter.get< std::pair< std::string, std::string > >( tIsLeaderString + "_dof_dependencies" ) ),
                            tDofTypeNames );
                    mSPs( iSP )->set_dof_type_list( tDofTypes, tDofTypeNames, tIsLeader );

                    // set dv dependencies
                    moris::Cell< moris::Cell< PDV_Type > > tDvTypes;
                    string_to_cell_of_cell(
                            std::get< 0 >( tSPParameter.get< std::pair< std::string, std::string > >( tIsLeaderString + "_dv_dependencies" ) ),
                            tDvTypes,
                            aDvTypeMap );
                    moris::Cell< std::string > tDvTypeNames;
                    string_to_cell(
                            std::get< 1 >( tSPParameter.get< std::pair< std::string, std::string > >( tIsLeaderString + "_dv_dependencies" ) ),
                            tDvTypeNames );
                    mSPs( iSP )->set_dv_type_list( tDvTypes, tDvTypeNames, tIsLeader );

                    // set leader properties
                    moris::Cell< moris::Cell< std::string > > tPropertyNamesPair;
                    string_to_cell_of_cell(
                            tSPParameter.get< std::string >( tIsLeaderString + "_properties" ),
                            tPropertyNamesPair );

                    for ( uint iProp = 0; iProp < tPropertyNamesPair.size(); iProp++ )
                    {
                        // get the property name
                        std::string tPropertyName = tPropertyNamesPair( iProp )( 0 );

                        // check for unknown property
                        MORIS_ERROR( aPropertyMap.find( tPropertyName ) != aPropertyMap.end(),
                                "FEM_Model::create_stabilization_parameters_without_phase - Unknown leader aPropertyString : %s \n",
                                tPropertyName.c_str() );

                        // get property index
                        uint tPropertyIndex = aPropertyMap[ tPropertyName ];

                        // set property for CM
                        mSPs( iSP )->set_property(
                                mProperties( tPropertyIndex ),
                                tPropertyNamesPair( iProp )( 1 ),
                                tIsLeader );
                    }

                    // set constitutive models
                    moris::Cell< moris::Cell< std::string > > tCMNamesPair;
                    string_to_cell_of_cell(
                            tSPParameter.get< std::string >( tIsLeaderString + "_constitutive_models" ),
                            tCMNamesPair );

                    // loop over the CM names
                    for ( uint iCM = 0; iCM < tCMNamesPair.size(); iCM++ )
                    {
                        // get the CM name
                        std::string tCMName = tCMNamesPair( iCM )( 0 );

                        // check for unknown CM
                        MORIS_ERROR( aCMMap.find( tCMName ) != aCMMap.end(),
                                "FEM_Model::create_stabilization_parameters_without_phase - Unknown leader aCMString: %s \n",
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
            return tSPMap;
        }

        //------------------------------------------------------------------------------

        void
        FEM_Model::create_IWGs_without_phase(
                std::map< std::string, uint >              &aPropertyMap,
                std::map< std::string, uint >              &aMMMap,
                std::map< std::string, uint >              &aCMMap,
                std::map< std::string, uint >              &aSPMap,
                moris::map< std::string, MSI::Dof_Type >   &aMSIDofTypeMap,
                moris::map< std::string, PDV_Type >        &aDvTypeMap,
                moris::map< std::string, mtk::Field_Type > &aFieldTypeMap )
        {
            std::map< std::string, uint > tIWGMap;
            IWG_Factory                   tIWGFactory;
            moris::Cell< ParameterList >  tIWGParameterList = mParameterList( 3 );
            uint const                    tNumIWGs          = tIWGParameterList.size();
            mIWGs.resize( tNumIWGs, nullptr );    // list of IWG pointers

            for ( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
            {
                // get the treated IWG parameter list
                ParameterList const tIWGParameter = tIWGParameterList( iIWG );

                fem::IWG_Type const tIWGType =
                        static_cast< fem::IWG_Type >( tIWGParameter.get< uint >( "IWG_type" ) );

                auto tIWG = tIWGFactory.create_IWG( tIWGType );

                std::string const tIWGName = tIWGParameter.get< std::string >( "IWG_name" );
                tIWG->set_name( tIWGName );

                // fill IWG map
                tIWGMap[ tIWGName ] = iIWG;

                // get function parameters
                auto tFuncParameters = string_to_cell_mat_2< DDRMat >( tIWGParameter.get< std::string >( "function_parameters" ) );
                tIWG->set_parameters( tFuncParameters );

                // get the ghost order from parameter list
                uint const tGhostOrder = tIWGParameter.get< uint >( "ghost_order" );
                tIWG->set_interpolation_order( tGhostOrder );

                // set residual dof type
                set_IWG_residual_dof_type( aMSIDofTypeMap, tIWGParameter, tIWG );

                // loop over leader and follower and set the appropriate properties
                moris::Cell< mtk::Leader_Follower > const tLeaderFollower{ mtk::Leader_Follower::LEADER, mtk::Leader_Follower::FOLLOWER };
                for ( auto const &tLeaderFollowerType : tLeaderFollower )
                {
                    set_IWG_dof_dependencies( aMSIDofTypeMap, tIWGParameter, tIWG, tLeaderFollowerType );
                    set_IWG_dv_dependencies( aDvTypeMap, tIWGParameter, tIWG, tLeaderFollowerType );
                    set_IWG_field_types( aFieldTypeMap, tIWGParameter, tIWG, tLeaderFollowerType );
                    set_IWG_properties( aPropertyMap, tIWGParameter, tIWG, tLeaderFollowerType );
                    set_IWG_material_models( aMMMap, tIWGParameter, tIWG, tLeaderFollowerType );
                    set_IWG_constitutive_models( aCMMap, tIWGParameter, tIWG, tLeaderFollowerType );
                }

                set_IWG_stabilization_parameters( aSPMap, tIWGParameter, tIWG );
                mIWGs( iIWG ) = tIWG;
                // debug
                // mIWGs( iIWG )->print_names();
            }
        }

        void FEM_Model::set_IWG_residual_dof_type(
                map< std::string, MSI::Dof_Type > &aMSIDofTypeMap,
                ParameterList const               &aIWGParameter,
                std::shared_ptr< IWG >            &aIWG ) const
        {
            std::string const tDofResidualString = aIWGParameter.get< std::string >( "dof_residual" );
            auto              tResDofTypes       = string_to_cell_of_cell< MSI::Dof_Type >( tDofResidualString, aMSIDofTypeMap );
            aIWG->set_residual_dof_type( tResDofTypes );
        }

        void FEM_Model::set_IWG_stabilization_parameters(
                std::map< std::string, uint > &aSPMap,
                ParameterList const           &aIWGParameter,
                std::shared_ptr< IWG >        &aIWG )
        {
            auto tSPNamesPair = string_to_cell_of_cell< std::string >( aIWGParameter.get< std::string >( "stabilization_parameters" ) );

            for ( uint iSP = 0; iSP < tSPNamesPair.size(); iSP++ )
            {
                // if CM name is in the CM map
                if ( aSPMap.find( tSPNamesPair( iSP )( 0 ) ) != aSPMap.end() )
                {
                    // get SP index
                    uint const tSPIndex = aSPMap[ tSPNamesPair( iSP )( 0 ) ];

                    // set SP for IWG
                    aIWG->set_stabilization_parameter(
                            mSPs( tSPIndex ),
                            tSPNamesPair( iSP )( 1 ) );
                }
                else
                {
                    // error message unknown SP
                    MORIS_ERROR( false,
                            "FEM_Model::create_IWGs_without_phase - Unknown aSPString: %s \n",
                            tSPNamesPair( iSP )( 0 ).c_str() );
                }
            }
        }

        void FEM_Model::set_IWG_constitutive_models(
                std::map< std::string, uint > &aCMMap,
                ParameterList const           &aIWGParameter,
                std::shared_ptr< IWG >        &aIWG,
                mtk::Leader_Follower const    &aLeaderFollowerType )
        {
            // get the prefix of the property name based on the leader or follower type (either "leader" or "follower")
            std::string const tPrefix      = mtk::get_leader_follower_string( aLeaderFollowerType );
            auto              tCMNamesPair = string_to_cell_of_cell< std::string >( aIWGParameter.get< std::string >( tPrefix + "_constitutive_models" ) );

            for ( uint iCM = 0; iCM < tCMNamesPair.size(); iCM++ )
            {
                // if CM name is in the CM map
                if ( aCMMap.find( tCMNamesPair( iCM )( 0 ) ) != aCMMap.end() )
                {
                    // get CM index
                    uint const tCMIndex = aCMMap[ tCMNamesPair( iCM )( 0 ) ];

                    // set CM for IWG
                    aIWG->set_constitutive_model(
                            mCMs( tCMIndex ),
                            tCMNamesPair( iCM )( 1 ),
                            aLeaderFollowerType );
                }
                else
                {
                    // error message unknown CM
                    MORIS_ERROR( false,
                            "FEM_Model::create_IWGs_without_phase - Unknown %s aCMString: %s \n",
                            tPrefix.c_str(),
                            tCMNamesPair( iCM )( 0 ).c_str() );
                }
            }
        }

        void FEM_Model::set_IWG_material_models(
                std::map< std::string, uint > &aMMMap,
                ParameterList const           &aIWGParameter,
                std::shared_ptr< IWG >        &aIWG,
                mtk::Leader_Follower const    &aLeaderFollowerType )
        {
            // get the prefix of the property name based on the leader or follower type (either "leader" or "follower")
            std::string const tPrefix      = mtk::get_leader_follower_string( aLeaderFollowerType );
            auto              tMMNamesPair = string_to_cell_of_cell< std::string >( aIWGParameter.get< std::string >( tPrefix + "_material_model" ) );

            for ( uint iMM = 0; iMM < tMMNamesPair.size(); iMM++ )
            {
                // if MM name is in the CM map
                if ( aMMMap.find( tMMNamesPair( iMM )( 0 ) ) != aMMMap.end() )
                {
                    // get MM index
                    uint const tMMIndex = aMMMap[ tMMNamesPair( iMM )( 0 ) ];

                    // set MM for IWG
                    aIWG->set_material_model(
                            mMMs( tMMIndex ),
                            tMMNamesPair( iMM )( 1 ),
                            aLeaderFollowerType );
                }
                else
                {
                    // error message unknown MM
                    MORIS_ERROR( false,
                            "FEM_Model::create_IWGs_without_phase - Unknown %s aMMString: %s \n",
                            tPrefix.c_str(),
                            tMMNamesPair( iMM )( 0 ).c_str() );
                }
            }
        }

        void FEM_Model::set_IWG_properties(
                std::map< std::string, uint > &aPropertyMap,
                ParameterList const           &aIWGParameter,
                std::shared_ptr< IWG >        &aIWG,
                mtk::Leader_Follower const    &aLeaderFollowerType )
        {
            // get the prefix of the property name based on the leader or follower type (either "leader" or "follower")
            std::string const tPrefix = mtk::get_leader_follower_string( aLeaderFollowerType );

            auto tPropertyNamesPair = string_to_cell_of_cell< std::string >(
                    aIWGParameter.get< std::string >( tPrefix + "_properties" ) );

            for ( uint iProp = 0; iProp < tPropertyNamesPair.size(); iProp++ )
            {
                // if property name is in the property map
                if ( aPropertyMap.find( tPropertyNamesPair( iProp )( 0 ) ) != aPropertyMap.end() )
                {
                    // get property index
                    uint const tPropertyIndex = aPropertyMap[ tPropertyNamesPair( iProp )( 0 ) ];

                    // set property for IWG
                    aIWG->set_property(
                            mProperties( tPropertyIndex ),
                            tPropertyNamesPair( iProp )( 1 ),
                            aLeaderFollowerType );
                }
                else
                {
                    // create error message unknown property
                    MORIS_ERROR( false,
                            "FEM_Model::create_IWGs_without_phase - Unknown %s aPropertyString: %s \n",
                            tPrefix.c_str(),
                            tPropertyNamesPair( iProp )( 0 ).c_str() );
                }
            }
        }

        void FEM_Model::set_IWG_field_types(
                map< std::string, mtk::Field_Type > &aFieldTypeMap,
                ParameterList const                 &aIWGParameter,
                std::shared_ptr< IWG >              &aIWG,
                mtk::Leader_Follower const          &aLeaderFollowerType ) const
        {
            // get the prefix of the property name based on the leader or follower type (either "leader" or "follower")
            std::string const tPrefix = mtk::get_leader_follower_string( aLeaderFollowerType );

            auto tFieldTypes = property_to_cell_of_cell( aIWGParameter, tPrefix + "_field_types", aFieldTypeMap );

            aIWG->set_field_type_list( tFieldTypes, aLeaderFollowerType );
        }

        void FEM_Model::set_IWG_dof_dependencies(
                map< std::string, MSI::Dof_Type > &aMSIDofTypeMap,
                ParameterList const               &aIWGParameter,
                std::shared_ptr< IWG >            &aIWG,
                mtk::Leader_Follower               aLeaderFollowerType ) const
        {
            // get the prefix of the property based on the leader or follower type
            std::string const tPrefix = mtk::get_leader_follower_string( aLeaderFollowerType );

            auto tDofTypes = property_to_cell_of_cell( aIWGParameter, tPrefix + "_dof_dependencies", aMSIDofTypeMap );

            aIWG->set_dof_type_list( tDofTypes, aLeaderFollowerType );
        }

        void FEM_Model::set_IWG_dv_dependencies(
                map< std::string, PDV_Type > &aDvTypeMap,
                ParameterList const          &aIWGParameter,
                std::shared_ptr< IWG >       &aIWG,
                mtk::Leader_Follower const   &aLeaderFollowerType ) const
        {

            // get the prefix of the property based on the leader or follower type
            std::string const tPrefix  = mtk::get_leader_follower_string( aLeaderFollowerType );
            auto              tDvTypes = property_to_cell_of_cell( aIWGParameter, tPrefix + "_dv_dependencies", aDvTypeMap );
            aIWG->set_dv_type_list( tDvTypes, aLeaderFollowerType );
        }


        //------------------------------------------------------------------------------

        void
        FEM_Model::create_IQIs_without_phase( std::map< std::string, uint > &aPropertyMap, std::map< std::string, uint > &aCMMap, std::map< std::string, uint > &aSPMap, moris::map< std::string, MSI::Dof_Type > &aMSIDofTypeMap, moris::map< std::string, PDV_Type > &aDvTypeMap, moris::map< std::string, mtk::Field_Type > &aFieldTypeMap )
        {
            std::map< std::string, uint > tIQIMap;

            // create an IQI factory
            IQI_Factory tIQIFactory;

            // get the IQI parameter list
            moris::Cell< ParameterList > tIQIParameterList = mParameterList( 4 );

            // get number of IQIs
            uint tNumIQIs = tIQIParameterList.size();

            // set size for list of IQI pointers
            mIQIs.resize( tNumIQIs, nullptr );

            // loop over the parameter lists
            for ( uint iIQI = 0; iIQI < tNumIQIs; iIQI++ )
            {
                // get the treated IQI parameter list
                ParameterList tIQIParameter = tIQIParameterList( iIQI );

                // get name from parameter list
                std::string tIQIName = tIQIParameter.get< std::string >( "IQI_name" );

                // get the IQI type from parameter list
                fem::IQI_Type tIQIType =
                        static_cast< fem::IQI_Type >( tIQIParameter.get< uint >( "IQI_type" ) );

                // get the treated IQI bulk type
                fem::Element_Type tIQIBulkType =
                        static_cast< fem::Element_Type >( tIQIParameter.get< uint >( "IQI_bulk_type" ) );

                // create an IQI pointer
                mIQIs( iIQI ) = tIQIFactory.create_IQI( tIQIType );

                // set name
                mIQIs( iIQI )->set_name( tIQIName );

                mIQIs( iIQI )->set_bulk_type( tIQIBulkType );

                // fill IQI map
                tIQIMap[ tIQIName ] = iIQI;

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

                // init string for leader or follower
                std::string          tIsLeaderString = "leader";
                mtk::Leader_Follower tIsLeader       = mtk::Leader_Follower::LEADER;

                // loop on leader and follower
                for ( uint iLeader = 0; iLeader <= 1; iLeader++ )
                {
                    // if follower
                    if ( iLeader )
                    {
                        // reset string for follower
                        tIsLeaderString = "follower";
                        tIsLeader       = mtk::Leader_Follower::FOLLOWER;
                    }

                    // set dof dependencies
                    moris::Cell< moris::Cell< moris::MSI::Dof_Type > > tDofTypes;
                    string_to_cell_of_cell(
                            tIQIParameter.get< std::string >( tIsLeaderString + "_dof_dependencies" ),
                            tDofTypes,
                            aMSIDofTypeMap );
                    mIQIs( iIQI )->set_dof_type_list( tDofTypes, tIsLeader );

                    // set dv dependencies
                    moris::Cell< moris::Cell< PDV_Type > > tDvTypes;
                    string_to_cell_of_cell(
                            tIQIParameter.get< std::string >( tIsLeaderString + "_dv_dependencies" ),
                            tDvTypes,
                            aDvTypeMap );
                    mIQIs( iIQI )->set_dv_type_list( tDvTypes, tIsLeader );

                    // set field types
                    moris::Cell< moris::Cell< moris::mtk::Field_Type > > tFieldTypes;
                    string_to_cell_of_cell(
                            tIQIParameter.get< std::string >( tIsLeaderString + "_field_types" ),
                            tFieldTypes,
                            aFieldTypeMap );
                    mIQIs( iIQI )->set_field_type_list( tFieldTypes, tIsLeader );

                    // set properties
                    moris::Cell< moris::Cell< std::string > > tPropertyNamesPair;
                    string_to_cell_of_cell(
                            tIQIParameter.get< std::string >( tIsLeaderString + "_properties" ),
                            tPropertyNamesPair );

                    for ( uint iProp = 0; iProp < tPropertyNamesPair.size(); iProp++ )
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
                                    tIsLeader );
                        }
                        else
                        {
                            // error message unknown property
                            MORIS_ERROR( false,
                                    "FEM_Model::create_IQIs_without_phase - Unknown %s aPropertyString: %s \n",
                                    tIsLeaderString.c_str(),
                                    tPropertyNamesPair( iProp )( 0 ).c_str() );
                        }
                    }

                    // set constitutive models
                    moris::Cell< moris::Cell< std::string > > tCMNamesPair;
                    string_to_cell_of_cell(
                            tIQIParameter.get< std::string >( tIsLeaderString + "_constitutive_models" ),
                            tCMNamesPair );

                    for ( uint iCM = 0; iCM < tCMNamesPair.size(); iCM++ )
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
                                    tIsLeader );
                        }
                        else
                        {
                            // error message unknown CM
                            MORIS_ERROR( false,
                                    "FEM_Model::create_IQIs_without_phase - Unknown %s aCMString: %s \n",
                                    tIsLeaderString.c_str(),
                                    tCMNamesPair( iCM )( 0 ).c_str() );
                        }
                    }
                }

                // set stabilization parameters
                moris::Cell< moris::Cell< std::string > > tSPNamesPair;
                string_to_cell_of_cell(
                        tIQIParameter.get< std::string >( "stabilization_parameters" ),
                        tSPNamesPair );

                for ( uint iSP = 0; iSP < tSPNamesPair.size(); iSP++ )
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
                        MORIS_ERROR( false,
                                "FEM_Model::create_IQIs_without_phase - Unknown aSPString: %s \n",
                                tSPNamesPair( iSP )( 0 ).c_str() );
                    }
                }

                // debug - uncomment if needed to list IQIs the code actually sees
                // mIQIs( iIQI )->print_names();
            }
        }

        //------------------------------------------------------------------------------

        void
        FEM_Model::create_fem_set_info_without_phase()
        {
            // get fem computation type parameter list
            ParameterList const tComputationParameterList = mParameterList( 5 )( 0 );

            bool const tIsAnalyticalSA = tComputationParameterList.get< bool >( "is_analytical_sensitivity" );
            auto const tFDSchemeForSA  = static_cast< fem::FDScheme_Type >( tComputationParameterList.get< uint >( "finite_difference_scheme" ) );
            real const tFDPerturbation = tComputationParameterList.get< real >( "finite_difference_perturbation_size" );

            // create a map of the set
            std::map< std::tuple< std::string, bool, bool >, uint > tMeshToFemSetIndex;

            this->create_fem_set_info_from_IWGs( tIsAnalyticalSA, tFDSchemeForSA, tFDPerturbation, tMeshToFemSetIndex );
            this->create_fem_set_info_from_IQIs( tIsAnalyticalSA, tFDSchemeForSA, tFDPerturbation, tMeshToFemSetIndex );
        }

        void FEM_Model::create_fem_set_info_from_IQIs(
                bool const                                               aIsAnalyticalSA,
                FDScheme_Type const                                     &aFDSchemeForSA,
                real const                                               aFDPerturbation,
                std::map< std::tuple< std::string, bool, bool >, uint > &aMeshToFemSetIndex )

        {
            moris::Cell< ParameterList > tIQIParameterLists = this->mParameterList( 4 );
            for ( uint iIQI = 0; iIQI < tIQIParameterLists.size(); iIQI++ )
            {
                ParameterList const          &tIQIParameterList = tIQIParameterLists( iIQI );
                std::shared_ptr< IQI > const &tIQI              = this->mIQIs( iIQI );

                bool tTimeContinuity = tIQIParameterList.get< bool >( "time_continuity" );
                bool tTimeBoundary   = tIQIParameterList.get< bool >( "time_boundary" );

                // loop over the mesh set names
                auto tMeshSetNames = string_to_cell< std::string >( tIQIParameterList.get< std::string >( "mesh_set_names" ) );
                for ( auto &tMeshSetName : tMeshSetNames )
                {
                    auto tMeshTuple = std::make_tuple( tMeshSetName, tTimeContinuity, tTimeBoundary );

                    // if the mesh set name not in map
                    if ( aMeshToFemSetIndex.find( tMeshTuple ) == aMeshToFemSetIndex.end() )
                    {
                        // add the mesh set name map
                        aMeshToFemSetIndex[ tMeshTuple ] = this->mSetInfo.size();

                        // create a fem set info for the mesh set
                        Set_User_Info aSetUserInfo;
                        aSetUserInfo.set_mesh_set_name( tMeshSetName );
                        aSetUserInfo.set_time_continuity( tTimeContinuity );
                        aSetUserInfo.set_time_boundary( tTimeBoundary );
                        aSetUserInfo.set_is_analytical_sensitivity_analysis( aIsAnalyticalSA );
                        aSetUserInfo.set_finite_difference_scheme_for_sensitivity_analysis( aFDSchemeForSA );
                        aSetUserInfo.set_finite_difference_perturbation_size( aFDPerturbation );
                        aSetUserInfo.add_IQI( tIQI );

                        // add it to the list of fem set info
                        this->mSetInfo.push_back( aSetUserInfo );
                    }
                    else
                    {
                        // set the IQI
                        this->mSetInfo( aMeshToFemSetIndex[ tMeshTuple ] ).add_IQI( tIQI );
                    }
                }
            }
        }

        void FEM_Model::create_fem_set_info_from_IWGs(
                bool const                                               aIsAnalyticalSA,
                FDScheme_Type const                                     &aFDSchemeForSA,
                real const                                               aFDPerturbation,
                std::map< std::tuple< std::string, bool, bool >, uint > &aMeshToFemSetIndex )
        {
            moris::Cell< ParameterList > tIWGParameterLists = this->mParameterList( 3 );
            for ( uint iIWG = 0; iIWG < tIWGParameterLists.size(); iIWG++ )
            {
                ParameterList const          &tIWGParameterList = tIWGParameterLists( iIWG );
                std::shared_ptr< IWG > const &tIWG              = this->mIWGs( iIWG );

                // get the time continuity and time boundary flags from the IWG parameter list to uniquely identify the fem sets
                bool const tTimeContinuity = tIWGParameterList.get< bool >( "time_continuity" );
                bool const tTimeBoundary   = tIWGParameterList.get< bool >( "time_boundary" );

                // get a representative DoF type
                MSI::Dof_Type const tFirstResidualDofType = tIWG->get_residual_dof_type()( 0 )( 0 );

                // loop over the mesh set names
                auto tMeshSetNames = string_to_cell< std::string >( tIWGParameterList.get< std::string >( "mesh_set_names" ) );
                for ( auto &tMeshSetName : tMeshSetNames )
                {
                    // check for ghost set names and select correct B-spline mesh automatically when new ghost sets need to be used
                    this->check_and_set_ghost_set_names( tMeshSetName, tFirstResidualDofType );

                    // create a tuple with the mesh set name, time continuity and time boundary flags that will be used to
                    // uniquely identify the fem set
                    auto tMeshTuple = std::make_tuple( tMeshSetName, tTimeContinuity, tTimeBoundary );

                    // check if the mesh set name already in map
                    if ( aMeshToFemSetIndex.find( tMeshTuple ) == aMeshToFemSetIndex.end() )
                    {
                        // if the set did not yet exist, create a new one and add keep track of the index
                        aMeshToFemSetIndex[ tMeshTuple ] = this->mSetInfo.size();

                        // create a fem set info for the mesh set
                        Set_User_Info aSetUserInfo;
                        aSetUserInfo.set_mesh_set_name( tMeshSetName );
                        aSetUserInfo.set_time_continuity( tTimeContinuity );
                        aSetUserInfo.set_time_boundary( tTimeBoundary );
                        aSetUserInfo.set_is_analytical_sensitivity_analysis( aIsAnalyticalSA );
                        aSetUserInfo.set_finite_difference_scheme_for_sensitivity_analysis( aFDSchemeForSA );
                        aSetUserInfo.set_finite_difference_perturbation_size( aFDPerturbation );
                        aSetUserInfo.add_IWG( tIWG );
                        this->mSetInfo.push_back( aSetUserInfo );
                    }
                    else
                    {
                        // if the fem set already exister, only add the IWG
                        this->mSetInfo( aMeshToFemSetIndex[ tMeshTuple ] ).add_IWG( tIWG );
                    }
                }
            }
        }

        //-------------------------------------------------------------------------------------------------

        const std::shared_ptr< fem::Field > &
        FEM_Model::get_field( mtk::Field_Type tFieldType )
        {
            size_t tIndex = mFieldTypeMap( static_cast< sint >( tFieldType ) );
            return mFields( tIndex );
        }

        //-------------------------------------------------------------------------------------------------

        moris::Cell< std::shared_ptr< mtk::Field > >
        FEM_Model::get_fields()
        {
            moris::Cell< std::shared_ptr< mtk::Field > > tFields( mFields.size(), nullptr );

            for ( uint Ik = 0; Ik < mFields.size(); Ik++ )
            {
                tFields( Ik ) = mFields( Ik );
            }

            return tFields;
        }
        //-------------------------------------------------------------------------------------------------

        void
        FEM_Model::populate_fields()
        {
            Tracer tTracer( "FEM", "Model", "Populate fields" );

            // check if fields exists
            if ( mFields.size() == 0 )
            {
                return;
            }

            Cell< std::shared_ptr< fem::Field > > tFieldToPopulate;
            Cell< std::string >                   tFieldIQINames;

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

        //-------------------------------------------------------------------------------------------------

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

        //-------------------------------------------------------------------------------------------------

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

        //------------------------------------------------------------------------------

        void
        FEM_Model::get_vertex_xyz_active_flags(
                moris_index                         aVertexIndex,
                Matrix< DDSMat >                   &aIsActiveDv,
                const moris::Cell< enum PDV_Type > &aPdvTypes )
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

        //-------------------------------------------------------------------------------------------------

        void
        FEM_Model::set_vertex_xyz_active_flags(
                moris_index                      aVertexIndex,
                moris::Cell< Matrix< DDSMat > > &aIsActiveDv )
        {
            // get num of pdv
            uint tNumXYZPdv = aIsActiveDv.size();

            // fill mIsActiveXYZ
            for ( uint iPdvType = 0; iPdvType < tNumXYZPdv; iPdvType++ )
            {
                mIsActiveXYZ( aVertexIndex, iPdvType ) = aIsActiveDv( iPdvType )( 0 );
            }
        }

        //------------------------------------------------------------------------------

        void
        FEM_Model::set_vertex_xyz_pdv_ids(
                moris_index                      aVertexIndex,
                moris::Cell< Matrix< DDSMat > > &aXYZPvIds )
        {
            // get num of pdv
            uint tNumXYZPdv = aXYZPvIds.size();

            // fill mXYZPdvIds
            for ( uint iPdvType = 0; iPdvType < tNumXYZPdv; iPdvType++ )
            {
                mXYZPdvIds( aVertexIndex, iPdvType ) = aXYZPvIds( iPdvType )( 0 );
            }
        }

        //------------------------------------------------------------------------------

        void
        FEM_Model::get_vertex_xyz_pdv_ids(
                moris_index                         aVertexIndex,
                Matrix< DDSMat >                   &aXYZPdvIds,
                const moris::Cell< enum PDV_Type > &aPdvTypes )
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

        //------------------------------------------------------------------------------

        void
        FEM_Model::get_local_xyz_pdv_assembly_indices(
                moris_index                         aVertexIndex,
                Matrix< DDSMat >                   &aXYZLocalAssemblyIndices,
                const moris::Cell< enum PDV_Type > &aPdvTypes )
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

        //------------------------------------------------------------------------------

        /**
         * @brief Set the dof type to Bspline mesh index map
         *
         * @param aDofTypeToBsplineMeshIndex
         */
        void
        FEM_Model::set_dof_type_to_Bspline_mesh_index(
                std::unordered_map< MSI::Dof_Type, moris_index > aDofTypeToBsplineMeshIndex )
        {
            mDofTypeToBsplineMeshIndex = aDofTypeToBsplineMeshIndex;
        }

        //------------------------------------------------------------------------------

        /**
         * @brief set flag whether to use new ghost sets
         *
         * @param aUseNewGhostSets
         */
        void
        FEM_Model::set_use_new_ghost_sets( bool aUseNewGhostSets )
        {
            mUseNewGhostSets = aUseNewGhostSets;
        }

        //------------------------------------------------------------------------------

        void
        FEM_Model::check_and_set_ghost_set_names(
                std::string       &aMeshSetName,
                enum MSI::Dof_Type aDofType )
        {
            // check whether new ghost sets should be used
            if ( mUseNewGhostSets )
            {
                if ( aMeshSetName.find( "ghost_p" ) != std::string::npos )
                {
                    // find the phase index
                    size_t tPos = aMeshSetName.find( 'p' );
                    MORIS_ERROR( tPos != std::string::npos, "FEM_Model::check_and_set_ghost_set_names() - Phase index not found in ghost set name." );
                    aMeshSetName.erase( 0, tPos + 1 );

                    moris_index tBsplineMeshIndex = mDofTypeToBsplineMeshIndex.find( aDofType )->second;

                    aMeshSetName = "ghost_B" + std::to_string( tBsplineMeshIndex ) + "_p" + aMeshSetName;
                }
            }
        }

        //------------------------------------------------------------------------------

    }    // namespace fem
} /* namespace moris */