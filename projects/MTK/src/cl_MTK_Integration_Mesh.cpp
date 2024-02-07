/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_MTK_Integration_Mesh.cpp
 *
 */

#include <string>
#include <iomanip>
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Mesh_Core.hpp"
#include "cl_MTK_Cell_Cluster.hpp"
#include "cl_MTK_Side_Cluster.hpp"
#include "cl_MTK_Double_Side_Cluster.hpp"
#include "assert.hpp"
#include "cl_Matrix.hpp"
#include "cl_MTK_Block_Set.hpp"
#include "cl_MTK_Space_Interpolator.hpp"
#include "cl_MTK_Integration_Rule.hpp"
#include "fn_norm.hpp"
#include "cl_MTK_Double_Side_Set.hpp"
#include "cl_MTK_Side_Set.hpp"

#include "HDF5_Tools.hpp"

namespace moris::mtk
{
    // ----------------------------------------------------------------------------

    Integration_Mesh::~Integration_Mesh()
    {
        for ( auto tListofBlocks : mListOfBlocks )
        {
            delete tListofBlocks;
        }
        mListOfBlocks.clear();

        for ( auto tListofSideSets : mListOfSideSets )
        {
            delete tListofSideSets;
        }
        mListOfSideSets.clear();

        for ( auto tListofDoubleSideSets : mListOfDoubleSideSets )
        {
            delete tListofDoubleSideSets;
        }
        mListOfDoubleSideSets.clear();

        for ( auto p : mDoubleSideClusters )
        {
            delete p;
        }
        mDoubleSideClusters.clear();
    }

    // ----------------------------------------------------------------------------

    void
    Integration_Mesh::get_Mesh_GEN_map( Vector< moris_index > &aMesh_GEN_map )
    {
        std::cout << "in integrtion mesh class" << std::endl;

        aMesh_GEN_map = this->mMesh_GEN_map;
    }

    // ----------------------------------------------------------------------------

    moris::uint
    Integration_Mesh::get_num_sets() const
    {
        return mListOfAllSets.size();
    }

    // ----------------------------------------------------------------------------

    moris::mtk::Set *
    Integration_Mesh::get_set_by_name( std::string aSetLabel ) const
    {
        moris_index tSetIndex = mSetNameToIndexMap.find( aSetLabel );

        return mListOfAllSets( tSetIndex );
    }

    // ----------------------------------------------------------------------------

    moris::mtk::Set *
    Integration_Mesh::get_set_by_index( moris_index aIndex ) const
    {
        return mListOfAllSets( aIndex );
    }

    // ----------------------------------------------------------------------------

    moris_index
    Integration_Mesh::get_set_index_by_name( std::string aSetLabel )
    {
        if ( !mSetNameToIndexMap.key_exists( aSetLabel ) )
        {
            MORIS_ERROR( false, "Integration_Mesh::get_set_index_by_name - Set with name does not exists: %s", aSetLabel.c_str() );
        }

        return mSetNameToIndexMap.find( aSetLabel );
    }

    // ----------------------------------------------------------------------------

    Vector< moris::mtk::Block_Set * > const &
    Integration_Mesh::get_block_sets_with_color( moris_index const &aColor )
    {
        MORIS_ASSERT( aColor <= mMaxColor, "Color above maximum color value" );
        return mColorToBlockSet( aColor );
    }

    // ----------------------------------------------------------------------------

    void
    Integration_Mesh::get_block_set_names_with_color(
            moris_index const          &aColor,
            Vector< std::string > &aSetNames )
    {
        MORIS_ASSERT( aColor <= mMaxColor, "Color above maximum color value" );

        uint tNumSets = mColorToBlockSet( aColor ).size();
        aSetNames.resize( tNumSets );
        for ( uint iS = 0; iS < tNumSets; iS++ )
        {
            aSetNames( iS ) = mColorToBlockSet( aColor )( iS )->get_set_name();
        }
    }

    // ----------------------------------------------------------------------------

    Vector< moris::mtk::Side_Set * > const &
    Integration_Mesh::get_side_sets_with_color( moris_index const &aColor )
    {
        MORIS_ASSERT( aColor <= mMaxColor, "Color above maximum color value" );
        return mColorToSideSet( aColor );
    }

    // ----------------------------------------------------------------------------

    Vector< moris::mtk::Double_Side_Set * > const &
    Integration_Mesh::get_double_side_sets_with_color( moris_index const &aColor )
    {
        MORIS_ASSERT( aColor <= mMaxColor, "Color above maximum color value" );
        return mColorToDoubleSideSet( aColor );
    }

    // ----------------------------------------------------------------------------


    Vector< moris::mtk::Set * > const &
    Integration_Mesh::get_all_sets_with_color( moris_index const &aColor )
    {
        MORIS_ASSERT( aColor <= mMaxColor, "Color above maximum color value" );
        return mColorToAllSets( aColor );
    }

    // ----------------------------------------------------------------------------
    void
    Integration_Mesh::print_sets_by_colors()
    {
        for ( moris::uint i = 0; i < mColorToBlockSet.size(); i++ )
        {
            std::cout << "\n Color: " << std::setw( 8 ) << i << std::endl;
            std::cout << "    Blocks: " << std::endl;
            for ( moris::uint iS = 0; iS < mColorToBlockSet( i ).size(); iS++ )
            {
                std::cout << "            " << mColorToBlockSet( i )( iS )->get_set_name() << std::endl;
            }

            std::cout << "    Side Sets: " << std::endl;
            for ( moris::uint iS = 0; iS < mColorToSideSet( i ).size(); iS++ )
            {
                std::cout << "            " << mColorToSideSet( i )( iS )->get_set_name() << std::endl;
            }

            std::cout << "    Double Side Sets: " << std::endl;
            for ( moris::uint iS = 0; iS < mColorToDoubleSideSet( i ).size(); iS++ )
            {
                std::cout << "            " << mColorToDoubleSideSet( i )( iS )->get_set_name() << std::endl;
            }
        }
    }
    // ----------------------------------------------------------------------------

    std::string
    Integration_Mesh::get_block_set_label( moris_index aBlockSetOrdinal ) const
    {
        MORIS_ERROR( 0, "get_block_set_label has no default implementation" );
        return "ERROR";
    }

    // ----------------------------------------------------------------------------

    moris_index
    Integration_Mesh::get_block_set_index( std::string aBlockSetLabel ) const
    {
        MORIS_ERROR( 0, "get_block_set_index has no default implementation" );
        return MORIS_INDEX_MAX;
    }

    // ----------------------------------------------------------------------------

    moris::uint
    Integration_Mesh::get_num_blocks() const
    {
        return mListOfBlocks.size();
    }

    // ----------------------------------------------------------------------------

    moris::uint
    Integration_Mesh::get_num_side_set() const
    {
        return mListOfSideSets.size();
    }

    // ----------------------------------------------------------------------------

    moris::uint
    Integration_Mesh::get_num_double_side_set() const
    {
        return mListOfDoubleSideSets.size();
    }

    // ----------------------------------------------------------------------------

    moris_index
    Integration_Mesh::get_double_sided_set_index( std::string aDoubleSideSetLabel ) const
    {
        MORIS_ERROR( false, "not implemented" );
        return 0;
    }


    // ----------------------------------------------------------------------------

    void
    Integration_Mesh::add_double_side_set( mtk::Double_Side_Set *aDblSideSet )
    {
        // add double sided set to the list
        mListOfDoubleSideSets.push_back( aDblSideSet );

        // reconstruct list of all sets by adding double sided set and a corresponding index
        this->collect_all_sets();
    }

    // ----------------------------------------------------------------------------

    void
    Integration_Mesh::collect_all_sets( bool aSetShape )
    {
        // reset
        mListOfAllSets.clear();
        mSetNameToIndexMap.clear();

        // reserve enough space
        mListOfAllSets.reserve( mListOfBlocks.size() + mListOfSideSets.size() + mListOfDoubleSideSets.size() );

        uint tCounter = 0;

        // Append block sets to list of all sets
        std::copy(
                mListOfBlocks.begin(),
                mListOfBlocks.end(),
                std::back_inserter( mListOfAllSets ) );

        // Append side sets to list of all sets
        std::copy(
                mListOfSideSets.begin(),
                mListOfSideSets.end(),
                std::back_inserter( mListOfAllSets ) );

        // FIXME implement cell topology for side set

        // Append double side sets to list of all sets
        std::copy(
                mListOfDoubleSideSets.begin(),
                mListOfDoubleSideSets.end(),
                std::back_inserter( mListOfAllSets ) );


        // iterate through all sets and register their index in the map
        for ( moris_index Ik = 0; Ik < static_cast< moris_index >( mListOfAllSets.size() ); Ik++ )
        {
            mListOfAllSets( Ik )->set_set_index( Ik );

            mSetNameToIndexMap[ mListOfAllSets( Ik )->get_set_name() ] = Ik;
        }

        // add the cell topology to the set
        for ( moris_index Ik = 0; Ik < static_cast< moris_index >( mListOfBlocks.size() ); Ik++ )
        {
            std::string tSetName = mListOfAllSets( tCounter )->get_set_name();

            // set the blockset cell topology and shape
            mListOfAllSets( tCounter )->set_cell_topology( this->get_blockset_topology( tSetName ) );
            if ( aSetShape )
            {
                mListOfAllSets( tCounter )->set_IG_cell_shape( this->get_IG_blockset_shape( tSetName ) );
                mListOfAllSets( tCounter )->set_IP_cell_shape( this->get_IP_blockset_shape( tSetName ) );
            }
            tCounter++;
        }

        if ( aSetShape )
        {
            // add the cell topology to the set
            for ( uint Ik = 0; Ik < mListOfSideSets.size(); Ik++ )
            {
                std::string tSetName = mListOfAllSets( tCounter )->get_set_name();

                // set the sideset cell topology and shape
                mListOfAllSets( tCounter )->set_IG_cell_shape( CellShape::STRAIGHT );
                mListOfAllSets( tCounter++ )->set_IP_cell_shape( CellShape::STRAIGHT );
            }

            // add the cell topology to the set
            for ( uint Ik = 0; Ik < mListOfDoubleSideSets.size(); Ik++ )
            {
                std::string tSetName = mListOfAllSets( tCounter )->get_set_name();

                // set the sideset cell topology and shape
                mListOfAllSets( tCounter )->set_IG_cell_shape( CellShape::STRAIGHT );
                mListOfAllSets( tCounter++ )->set_IP_cell_shape( CellShape::STRAIGHT );
            }
        }

        // setup color to set data
        this->setup_set_to_color();
    }

    // ----------------------------------------------------------------------------
    void
    Integration_Mesh::setup_set_to_color()
    {
        // determine maximum color for the sets (for data sizing)
        mMaxColor = 0;
        for ( moris::uint i = 0; i < mListOfAllSets.size(); i++ )
        {
            Matrix< IndexMat > const &tSetColors = mListOfAllSets( i )->get_set_colors();

            // if the max set color is greater than current max
            // replace the current max
            if ( tSetColors.numel() > 0 )
            {
                moris_index tMaxSetColor = tSetColors.max();
                if ( mMaxColor < tMaxSetColor )
                {
                    mMaxColor = tMaxSetColor;
                }
            }
        }

        // clear old data
        mColorToBlockSet.clear();
        mColorToSideSet.clear();
        mColorToDoubleSideSet.clear();
        mColorToAllSets.clear();

        // size outer cell size of member data
        mColorToBlockSet.resize( mMaxColor + 1 );
        mColorToSideSet.resize( mMaxColor + 1 );
        mColorToDoubleSideSet.resize( mMaxColor + 1 );
        mColorToAllSets.resize( mMaxColor + 1 );

        // iterate through block sets
        for ( moris::uint i = 0; i < mListOfBlocks.size(); i++ )
        {
            Matrix< IndexMat > const &tSetColors = mListOfBlocks( i )->get_set_colors();

            // iterate through the colors and add to related color grouping
            for ( moris::uint j = 0; j < tSetColors.numel(); j++ )
            {
                mColorToBlockSet( tSetColors( j ) ).push_back( mListOfBlocks( i ) );
                mColorToAllSets( tSetColors( j ) ).push_back( mListOfBlocks( i ) );
            }
        }

        // iterate through side sets
        for ( moris::uint i = 0; i < mListOfSideSets.size(); i++ )
        {
            Matrix< IndexMat > const &tSetColors = mListOfSideSets( i )->get_set_colors();

            // iterate through the colors and add to related color grouping
            for ( moris::uint j = 0; j < tSetColors.numel(); j++ )
            {
                mColorToSideSet( tSetColors( j ) ).push_back( mListOfSideSets( i ) );
                mColorToAllSets( tSetColors( j ) ).push_back( mListOfSideSets( i ) );
            }
        }

        // iterate through double side sets
        for ( moris::uint i = 0; i < mListOfDoubleSideSets.size(); i++ )
        {
            Matrix< IndexMat > const &tSetColors = mListOfDoubleSideSets( i )->get_set_colors();

            // iterate through the colors and add to related color grouping
            for ( moris::uint j = 0; j < tSetColors.numel(); j++ )
            {
                mColorToDoubleSideSet( tSetColors( j ) ).push_back( mListOfDoubleSideSets( i ) );
                mColorToAllSets( tSetColors( j ) ).push_back( mListOfDoubleSideSets( i ) );
            }
        }
    }

    // ----------------------------------------------------------------------------
    void
    Integration_Mesh::add_double_sided_cluster( mtk::Double_Side_Cluster *aDblSidedCluster )
    {
        // add double sided cluster to the double sided cluster list
        mDoubleSideClusters.push_back( aDblSidedCluster );
    }

    // ----------------------------------------------------------------------------

    void
    Integration_Mesh::save_MPC_to_hdf5( std::string aFileName )
    {
        // get number of integration vertices
        uint tNumVertices = this->get_num_nodes();

        std::cout << "Num Nodes: " << tNumVertices << std::endl;

        Matrix< IdMat > tIdentifierMat( tNumVertices, 1, -1 );

        Matrix< IdMat >                 tBSToIPMap( this->get_num_basis_functions( 0 ), 1, gNoID );
        Matrix< IdMat >                 tIPToIGMap( tNumVertices, 1, gNoID );
        Vector< Matrix< IdMat > >  tBSToIPIds( tNumVertices );
        Vector< Matrix< DDRMat > > tBSToIPWeights( tNumVertices );

        Vector< Matrix< IdMat > >  tNodeDOFs( tNumVertices );
        Vector< Matrix< DDRMat > > tMPCs( tNumVertices );

        this->create_MPC_maps( tBSToIPMap, tIPToIGMap );

        this->build_hanging_node_MPC( tBSToIPIds, tBSToIPWeights, tBSToIPMap, tIPToIGMap );

        // loop over all bulk sets
        for ( moris_index Ik = 0; Ik < static_cast< moris_index >( this->get_num_blocks() ); Ik++ )
        {
            // get num clusters on set
            uint tNumClustersOnSet = this->get_set_by_index( Ik )->get_num_clusters_on_set();

            if ( tNumClustersOnSet > 0 )
            {
                // Geometry_Type tGeoType = aCell->get_geometry_type();
                // Interpolation_Order tIPOrder = aCell->get_interpolation_order();

                // creating interpolation rule
                Interpolation_Rule tInterpolationRule( mtk::Geometry_Type::HEX,
                        Interpolation_Type::LAGRANGE,
                        mtk::Interpolation_Order::LINEAR,
                        Interpolation_Type::LAGRANGE,
                        Interpolation_Order::LINEAR );

                // create a space interpolator
                Space_Interpolator tSpaceInterpolator( tInterpolationRule );

                // loop over all clusters in set
                for ( uint Ii = 0; Ii < tNumClustersOnSet; Ii++ )
                {
                    const Cluster *tCluster = this->get_set_by_index( Ik )->get_clusters_by_index( static_cast< moris_index >( Ii ) );

                    moris::Matrix< moris::DDRMat > tLocalCoords = tCluster->get_vertices_local_coordinates_wrt_interp_cell();

                    Vector< moris::mtk::Vertex const * > tIGVertices = tCluster->get_vertices_in_cluster();

                    // get number of vertices on the treated mesh cluster
                    uint tNumNodes = tLocalCoords.n_rows();

                    const moris::mtk::Cell &tInterpolationCell = tCluster->get_interpolation_cell();

                    Vector< Vertex * > tIPVertices = tInterpolationCell.get_vertex_pointers();

                    Matrix< IdMat > tIdMat( tIGVertices.size(), 1, 0 );


                    for ( uint Iv = 0; Iv < tIPVertices.size(); Iv++ )
                    {
                        if ( tIPVertices( Iv )->get_interpolation( 0 )->get_ids().numel() > 1 )
                        {
                            tIdMat( Iv ) = 1;
                        }
                    }

                    Matrix< IndexMat > tHangingNodeIdentifier( tNumVertices, 1, 0 );
                    if ( !tCluster->is_trivial() )
                    {
                        Matrix< IndexMat > tHangingNodeIndices = tCluster->get_hanging_nodes();

                        print( tHangingNodeIndices, "tHangingNodeIndices" );

                        for ( uint Iv = 0; Iv < tHangingNodeIndices.numel(); Iv++ )
                        {
                            tHangingNodeIdentifier( tHangingNodeIndices( Iv ) ) = 1;
                        }
                    }


                    // epsilon to count T-Matrix
                    real tEpsilon = 1e-12;

                    // loop over the vertices on the treated mesh cluster
                    for ( uint Iv = 0; Iv < tNumNodes; Iv++ )
                    {
                        tSpaceInterpolator.set_space_time( trans( tLocalCoords.get_row( Iv ) ) );

                        const Matrix< DDRMat > &tN = tSpaceInterpolator.NXi();

                        moris_id tIGVertexId    = tIGVertices( Iv )->get_id();
                        moris_id tIGVertexIndex = tIGVertices( Iv )->get_index();

                        uint tNCols = tN.n_cols();

                        tMPCs( tIGVertexIndex ).set_size( tNCols, 1, -1.0 );
                        tNodeDOFs( tIGVertexIndex ).set_size( tNCols, 1, gNoID );

                        std::cout << "tIGVertexId: " << tIGVertexId << " " << tIGVertexIndex << std::endl;

                        if ( tHangingNodeIdentifier( tIGVertexIndex ) == 1 )
                        {
                            // fixme add hanging - hanging node

                            uint tCount = 0;
                            // loop over all nonzero entries
                            for ( uint i = 0; i < tNCols; i++ )
                            {
                                if ( std::abs( tN( i ) ) > tEpsilon )
                                {
                                    // copy entry of T-Matrix
                                    tMPCs( tIGVertexIndex )( tCount ) = tN( i );

                                    // copy pointer of dof and convert to mtk::Vertex
                                    tNodeDOFs( tIGVertexIndex )( tCount++ ) = tIGVertices( i )->get_id();
                                }
                            }

                            //                                tMPCs.resize( tCount, 1 );
                            //                                tNodeDOFs.resize( tCount,1 );
                        }
                        else
                        {
                            if ( tIdMat( Iv ) == 1 )
                            {
                                for ( uint Ib = 0; Ib < tBSToIPWeights( tIGVertexId ).numel(); Ib++ )
                                {
                                    tMPCs( tIGVertexIndex )( Ib )     = tBSToIPWeights( tIGVertexId )( Ib );
                                    tNodeDOFs( tIGVertexIndex )( Ib ) = tBSToIPIds( tIGVertexId )( Ib );
                                }
                            }
                            else
                            {
                                tMPCs( tIGVertexIndex )( 0 )     = 1.0;
                                tNodeDOFs( tIGVertexIndex )( 0 ) = tIGVertexId;
                            }
                            //                                tMPCs( tIGVertexIndex ).resize( 1, 1 );
                            //                                tNodeDOFs( tIGVertexIndex ).resize( 1,1 );
                        }
                    }
                }
            }
        }

        print( tNodeDOFs, "tNodeDOFs" );
        print( tMPCs, "tMPCs" );

        Matrix< IdMat >  tIDs( 1, tNumVertices );
        Matrix< IdMat >  tMPCIDs( tMPCs( 0 ).numel(), tNumVertices );
        Matrix< DDRMat > tMPCWeightss( tMPCs( 0 ).numel(), tNumVertices );

        for ( uint Ik = 0; Ik < tNumVertices; Ik++ )
        {
            tIDs( Ik ) = static_cast< int >( Ik ) + 1;

            for ( uint Ii = 0; Ii < tMPCs( 0 ).numel(); Ii++ )
            {
                tMPCIDs( Ii, Ik )      = tNodeDOFs( Ik )( Ii );
                tMPCWeightss( Ii, Ik ) = tMPCs( Ik )( Ii );
            }
        }

        //            // add order to path
        //            std::string tFilePath =    aFilePath.substr(0,aFilePath.find_last_of(".")) // base path
        //                                                                                                              + "_" + std::to_string( tMesh->get_order() ) // rank of this processor
        //            +  aFilePath.substr( aFilePath.find_last_of("."), aFilePath.length() );

        // make path parallel
        aFileName = parallelize_path( aFileName );

        // Create a new file using default properties
        hid_t tFileID = H5Fcreate( aFileName.c_str(),
                H5F_ACC_TRUNC,
                H5P_DEFAULT,
                H5P_DEFAULT );

        // error handler
        herr_t tStatus;

        // save ids to file
        save_matrix_to_hdf5_file(
                tFileID,
                "NodeID",
                tIDs,
                tStatus );

        // generate label
        std::string tLabel = "MPC_IDs";

        // save ids to file
        save_matrix_to_hdf5_file(
                tFileID,
                tLabel,
                tMPCIDs,
                tStatus );

        // generate  label
        tLabel = "MPC_weights";

        // save weights to file
        save_matrix_to_hdf5_file(
                tFileID,
                tLabel,
                tMPCWeightss,
                tStatus );

        // close file
        tStatus = H5Fclose( tFileID );
    }

    // ----------------------------------------------------------------------------

    void
    Integration_Mesh::create_MPC_maps(
            Matrix< IdMat > &tBSToIPMap,
            Matrix< IdMat > &tIPToIGMap ) const
    {
        Matrix< DDSMat > tIdentifierMatBS( tBSToIPMap.numel(), 1, 0 );
        Matrix< DDSMat > tIdentifierMatIP( tIPToIGMap.numel(), 1, 0 );

        // loop over all bulk sets
        for ( uint Ik = 0; Ik < this->get_num_blocks(); Ik++ )
        {
            // get num clusters on set
            uint tNumClustersOnSet = this->get_set_by_index( static_cast< moris_index >( Ik ) )->get_num_clusters_on_set();

            if ( tNumClustersOnSet > 0 )
            {
                // loop over all clusters in set
                for ( uint Ii = 0; Ii < tNumClustersOnSet; Ii++ )
                {
                    const Cluster *tCluster = this->get_set_by_index( static_cast< moris_index >( Ik ) )
                                                      ->get_clusters_by_index( static_cast< moris_index >( Ii ) );

                    const moris::mtk::Cell &tInterpolationCell = tCluster->get_interpolation_cell();

                    Vector< Vertex * > tIPVertices = tInterpolationCell.get_vertex_pointers();

                    Vector< moris::mtk::Vertex const * > tIGVertices = tCluster->get_vertices_in_cluster();

                    std::cout << " IGVert: " << tIGVertices.size() << " IPVert: " << tIPVertices.size() << std::endl;

                    for ( uint Ij = 0; Ij < tIPVertices.size(); Ij++ )
                    {
                        MORIS_ASSERT( tIPVertices( Ij )->has_interpolation( 0 ), "Vertex does not have interpolation" );

                        uint tIP_ID = tIPVertices( Ij )->get_id();

                        // build BSpline to iP map
                        if ( tIPVertices( Ij )->get_interpolation( 0 )->get_ids().numel() == 1 )
                        {
                            moris_id tBSId = tIPVertices( Ij )->get_interpolation( 0 )->get_ids()( 0 );

                            if ( tIdentifierMatBS( tBSId ) == 0 )
                            {
                                tBSToIPMap( tBSId ) = tIP_ID;

                                tIdentifierMatBS( tBSId ) = 1;
                            }
                        }

                        // build IP to IG map
                        // assume corner vertices are first 8
                        uint tIG_ID = tIGVertices( Ij )->get_id();

                        if ( tIdentifierMatIP( tIP_ID ) == 0 )
                        {
                            tIPToIGMap( tIP_ID ) = tIG_ID;

                            tIdentifierMatIP( tIP_ID ) = 1;
                        }
                    }
                }
            }
        }
    }

    // ----------------------------------------------------------------------------

    void
    Integration_Mesh::build_hanging_node_MPC(
            Vector< Matrix< IdMat > >  &tBSToIPIds,
            Vector< Matrix< DDRMat > > &tBSToIPWeights,
            const Matrix< IdMat >           &tBSToIPMap,
            const Matrix< IdMat >           &tIPToIGMap )
    {
        // loop over all bulk sets
        for ( uint Ik = 0; Ik < this->get_num_blocks(); Ik++ )
        {
            // get num clusters on set
            uint tNumClustersOnSet = this->get_set_by_index( Ik )->get_num_clusters_on_set();

            if ( tNumClustersOnSet > 0 )
            {
                // loop over all clusters in set
                for ( uint Ii = 0; Ii < tNumClustersOnSet; Ii++ )
                {
                    const Cluster *tCluster = this->get_set_by_index( Ik )->get_clusters_by_index( Ii );

                    const moris::mtk::Cell &tInterpolationCell = tCluster->get_interpolation_cell();

                    Vector< Vertex * > tIPVertices = tInterpolationCell.get_vertex_pointers();

                    Vector< moris::mtk::Vertex const * > tIGVertices = tCluster->get_vertices_in_cluster();

                    for ( uint Ij = 0; Ij < tIPVertices.size(); Ij++ )
                    {
                        // build BSpline to iP map
                        if ( tIPVertices( Ij )->get_interpolation( 0 )->get_ids().numel() > 1 )
                        {
                            uint tIG_ID = tIGVertices( Ij )->get_id();

                            Matrix< IdMat >  tBSId      = tIPVertices( Ij )->get_interpolation( 0 )->get_ids();
                            Matrix< DDRMat > tBSWeights = *tIPVertices( Ij )->get_interpolation( 0 )->get_weights();

                            tBSToIPIds( tIG_ID ).resize( tBSId.numel(), 1 );
                            tBSToIPWeights( tIG_ID ).resize( tBSId.numel(), 1 );

                            for ( uint Ia = 0; Ia < tBSId.numel(); Ia++ )
                            {
                                tBSToIPWeights( tIG_ID )( Ia ) = tBSWeights( Ia );
                                tBSToIPIds( tIG_ID )( Ia )     = tIPToIGMap( tBSToIPMap( tBSId( Ia ) ) );
                            }
                        }
                    }
                }
            }
        }
    }

    Vector< moris::mtk::Block_Set * > const &Integration_Mesh::get_block_sets() const
    {
        return mListOfBlocks;
    }


    Vector< moris::mtk::Side_Set * > const &Integration_Mesh::get_side_sets() const
    {
        return mListOfSideSets;
    }

    Vector< moris::mtk::Double_Side_Set * > const &Integration_Mesh::get_double_side_sets() const
    {
        return mListOfDoubleSideSets;
    }

    Vector< moris::mtk::Set * > const &Integration_Mesh::get_sets() const
    {
        return mListOfAllSets;
    }

    Vector< moris::mtk::Set * > const &
    Integration_Mesh::get_list_of_sets( SetType aSetType ) const
    {
        // Determine the set type
        switch ( aSetType )
        {
            // Bulk sets
            case SetType::BULK:
            {
                MORIS_LOG_WARNING(
                        "You are using the deprecated function Integration_Mesh::get_sets( SetType::BULK )"
                        "Consider using get_block_sets() instead" );
                return (Vector< Set * > const &)this->get_block_sets();
            }

            // Side Sets
            case SetType::SIDESET:
            {
                MORIS_LOG_WARNING(
                        "You are using the deprecated function Integration_Mesh::get_sets( SetType::SIDESET )"
                        "Consider using get_side_sets() instead" );
                return (Vector< Set * > const &)this->get_side_sets();
            }

            // Double Sided Sets
            case SetType::DOUBLE_SIDED_SIDESET:
            {
                MORIS_LOG_WARNING(
                        "You are using the deprecated function Integration_Mesh::get_sets( SetType::DOUBLE_SIDED_SIDESET )"
                        "Consider using get_double_side_sets() instead" );
                return (Vector< Set * > const &)this->get_double_side_sets();
            }

            // All sets
            case SetType::UNDEFINED:
            {
                return mListOfAllSets;
            }

            // Incorrect Input
            default:
            {
                MORIS_ERROR( 0, "Incorrect aSetType" );
                auto *tDummyCell = new Vector< moris::mtk::Set * >( 0 );
                return *tDummyCell;
            }
        }
    }

    // ----------------------------------------------------------------------------

    Vector< Vector< moris::mtk::Block_Set * > > const &
    Integration_Mesh::get_color_to_block_sets() const
    {
        return mColorToBlockSet;
    }


    Vector< Vector< moris::mtk::Side_Set * > > const &
    Integration_Mesh::get_color_to_side_sets() const
    {
        return mColorToSideSet;
    }

    Vector< Vector< moris::mtk::Double_Side_Set * > > const &
    Integration_Mesh::get_color_to_double_side_sets() const
    {
        return mColorToDoubleSideSet;
    }

    Vector< Vector< moris::mtk::Set * > > const &
    Integration_Mesh::get_color_to_sets() const
    {
        return mColorToAllSets;
    }

    // ----------------------------------------------------------------------------

    map< std::string, moris_index > const &
    Integration_Mesh::get_set_name_to_index_map() const
    {
        return mSetNameToIndexMap;
    }

    // ----------------------------------------------------------------------------

    moris_index const &
    Integration_Mesh::get_max_color() const
    {
        return mMaxColor;
    }

    // ----------------------------------------------------------------------------

    void
    Integration_Mesh::delete_visualization_sets()
    {
        // lambda function to determine if the set name has a ghost in it
        auto tIsGhostVisualization = []( moris::mtk::Set *aSet ) -> bool { return std::strstr( aSet->get_set_name().c_str(), "ghost" ); };

        // use erase remove idiom to remove sets that have the name as the ghost
        // FIXME : a proper insert function needs to be implemented that we don't use data
        mListOfBlocks.data().erase( std::remove_if( mListOfBlocks.begin(), mListOfBlocks.end(), tIsGhostVisualization ), mListOfBlocks.end() );

        mListOfSideSets.data().erase( std::remove_if( mListOfSideSets.begin(), mListOfSideSets.end(), tIsGhostVisualization ), mListOfSideSets.end() );

        mListOfAllSets.data().erase( std::remove_if( mListOfAllSets.begin(), mListOfAllSets.end(), tIsGhostVisualization ), mListOfAllSets.end() );
    }
}    // namespace moris::mtk
