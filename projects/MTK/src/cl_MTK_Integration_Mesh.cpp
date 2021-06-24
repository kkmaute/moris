#include <string>
#include <iomanip>
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Mesh_Core.hpp"
#include "cl_MTK_Cell_Cluster.hpp"
#include "cl_MTK_Side_Cluster.hpp"
#include "cl_MTK_Double_Side_Cluster.hpp"
#include "assert.hpp"
#include "cl_Matrix.hpp"
#include "cl_MTK_Block.hpp"
#include "cl_MTK_Space_Interpolator.hpp"
#include "cl_MTK_Integration_Rule.hpp"
#include "fn_norm.hpp"

#include "HDF5_Tools.hpp"

namespace moris
{
    namespace mtk
    {
        // ----------------------------------------------------------------------------

        Integration_Mesh::~Integration_Mesh()
        {
            for( auto tListofBlocks : mListofBlocks )
            {
                delete tListofBlocks;
            }
            mListofBlocks.clear();

            for( auto tListofSideSets : mListofSideSets )
            {
                delete tListofSideSets;
            }
            mListofSideSets.clear();

            for( auto tListofDoubleSideSets : mListofDoubleSideSets )
            {
                delete tListofDoubleSideSets;
            }
            mListofDoubleSideSets.clear();

            for(auto p:mDoubleSideClusters)
            {
                delete p;
            }
            mDoubleSideClusters.clear();
        }

        // ----------------------------------------------------------------------------

        moris::uint Integration_Mesh::get_num_sets() const
        {
            return mListOfAllSets.size();
        }

        // ----------------------------------------------------------------------------

        moris::mtk::Set * Integration_Mesh::get_set_by_name( std::string aSetLabel ) const
        {
            moris_index tSetIndex = mSetNameToIndexMap.find( aSetLabel );

            return mListOfAllSets( tSetIndex );
        }

        // ----------------------------------------------------------------------------

        moris::mtk::Set * Integration_Mesh::get_set_by_index( moris_index aIndex ) const
        {
            return mListOfAllSets( aIndex );
        }

        // ----------------------------------------------------------------------------

        moris_index Integration_Mesh::get_set_index_by_name( std::string aSetLabel )
        {
            if ( ! mSetNameToIndexMap.key_exists( aSetLabel ) )
            {
                std::string tErrMsg = "Integration_Mesh::get_set_index_by_name - Set with name does not exists: " + aSetLabel;

                MORIS_ERROR( false, tErrMsg.c_str() );
            }

            return mSetNameToIndexMap.find( aSetLabel );
        }

        // ----------------------------------------------------------------------------
 
        moris::Cell<moris::mtk::Set*> const &
        Integration_Mesh::get_block_sets_with_color(moris_index const & aColor)
        {
            MORIS_ASSERT(aColor <= mMaxColor,"Color above maximum color value");
            return mBlockSetToColor(aColor);
        }

        // ----------------------------------------------------------------------------

        void
        Integration_Mesh::get_block_set_names_with_color(
                moris_index const        & aColor,
                moris::Cell<std::string> & aSetNames)
        {
            MORIS_ASSERT(aColor <= mMaxColor,"Color above maximum color value");

            uint tNumSets = mBlockSetToColor(aColor).size();
            aSetNames.resize( tNumSets );
            for(uint iS=0;iS<tNumSets;iS++)
            {
                aSetNames(iS) = mBlockSetToColor(aColor)(iS)->get_set_name();
            }
        }

        // ----------------------------------------------------------------------------

        moris::Cell<moris::mtk::Set*> const &
        Integration_Mesh::get_side_sets_with_color(moris_index const & aColor)
        {
            MORIS_ASSERT(aColor <= mMaxColor,"Color above maximum color value");
            return mSideSetToColor(aColor);
        }

        // ----------------------------------------------------------------------------

        moris::Cell<moris::mtk::Set*> const &
        Integration_Mesh::get_double_side_sets_with_color(moris_index const & aColor)
        {
            MORIS_ASSERT(aColor <= mMaxColor,"Color above maximum color value");
            return mDoubleSideSetToColor(aColor);
        }

        // ----------------------------------------------------------------------------


        moris::Cell<moris::mtk::Set*> const &
        Integration_Mesh::get_all_sets_with_color(moris_index const & aColor)
        {
            MORIS_ASSERT(aColor <= mMaxColor,"Color above maximum color value");
            return mAllSetToColor(aColor);
        }

        // ----------------------------------------------------------------------------
        void
        Integration_Mesh::print_sets_by_colors()
        {
            for(moris::uint i = 0; i < mBlockSetToColor.size(); i++)
            {
                std::cout<<"\n Color: "<<std::setw(8)<<i<<std::endl;
                std::cout<<"    Blocks: "<<std::endl;
                for(moris::uint iS = 0; iS <mBlockSetToColor(i).size(); iS++ )
                {
                    std::cout<<"            "<<mBlockSetToColor(i)(iS)->get_set_name()<<std::endl;
                }

                std::cout<<"    Side Sets: "<<std::endl;
                for(moris::uint iS = 0; iS <mSideSetToColor(i).size(); iS++ )
                {
                    std::cout<<"            "<<mSideSetToColor(i)(iS)->get_set_name()<<std::endl;
                }

                std::cout<<"    Double Side Sets: "<<std::endl;
                for(moris::uint iS = 0; iS <mDoubleSideSetToColor(i).size(); iS++ )
                {
                    std::cout<<"            "<<mDoubleSideSetToColor(i)(iS)->get_set_name()<<std::endl;
                }
            }
        }
        // ----------------------------------------------------------------------------

        std::string Integration_Mesh::get_block_set_label(moris_index aBlockSetOrdinal) const
        {
            MORIS_ERROR(0, "get_block_set_label has no default implementation");
            return "ERROR";
        }

        // ----------------------------------------------------------------------------

        moris_index Integration_Mesh::get_block_set_index(std::string aBlockSetLabel) const
        {
            MORIS_ERROR(0, "get_block_set_index has no default implementation");
            return MORIS_INDEX_MAX;
        }

        // ----------------------------------------------------------------------------

        moris::uint
        Integration_Mesh::get_num_blocks() const
        {
            return mListofBlocks.size();
        }

        // ----------------------------------------------------------------------------

        moris::uint
        Integration_Mesh::get_num_side_set() const
        {
            return mListofSideSets.size();
        }

        // ----------------------------------------------------------------------------

        moris::uint
        Integration_Mesh::get_num_double_side_set() const
        {
            return mListofDoubleSideSets.size();
        }

        // ----------------------------------------------------------------------------

        moris_index Integration_Mesh::get_double_sided_set_index(std::string aDoubleSideSetLabel) const
        {
            MORIS_ERROR( false, "not implemented");
            return 0;
        }



        // ----------------------------------------------------------------------------

        void
        Integration_Mesh::add_double_side_set(mtk::Set* aDblSideSet)
        {
            //add double sided set to the list
            mListofDoubleSideSets.push_back(aDblSideSet); 

            //reconstruct list of all sets by adding double sided set and a corresponding index
            this->collect_all_sets();  
        }

        // ----------------------------------------------------------------------------

        void 
        Integration_Mesh::collect_all_sets( bool aSetShape )
        {
            // reset
            mListOfAllSets.clear();
            mSetNameToIndexMap.clear();

            uint tCounter = 0;

            // Append block sets to list of all sets
            mListOfAllSets.append( mListofBlocks );

            // Append side sets to list of all sets
            mListOfAllSets.append( mListofSideSets );

            //FIXME implement cell topology for side set

            // Append double side sets to list of all sets
            mListOfAllSets.append( mListofDoubleSideSets );

            // iterate through all sets and register their index in the map
            for( uint Ik = 0; Ik < mListOfAllSets.size(); Ik++ )
            {
                mListOfAllSets( Ik )->set_set_index( Ik );

                mSetNameToIndexMap[ mListOfAllSets( Ik )->get_set_name() ] = Ik;
            }

            // add the cell topology to the set
            for( uint Ik = 0; Ik < mListofBlocks.size(); Ik++ )
            {
                std::string tSetName = mListOfAllSets( tCounter )->get_set_name();

                // set the blockset cell topology and shape
                mListOfAllSets( tCounter )->set_cell_topology( this->get_blockset_topology( tSetName ) );
                if ( aSetShape == true )
                {
                    mListOfAllSets( tCounter )->set_IG_cell_shape( this->get_IG_blockset_shape( tSetName ) );
                    mListOfAllSets( tCounter )->set_IP_cell_shape( this->get_IP_blockset_shape( tSetName ) );
                }
                tCounter++;
            }

            if ( aSetShape == true )
            {
                // add the cell topology to the set
                for( uint Ik = 0; Ik < mListofSideSets.size(); Ik++ )
                {
                    std::string tSetName = mListOfAllSets( tCounter )->get_set_name();

                    // set the sideset cell topology and shape
                    mListOfAllSets( tCounter   )->set_IG_cell_shape( CellShape::STRAIGHT );
                    mListOfAllSets( tCounter++ )->set_IP_cell_shape( CellShape::STRAIGHT );
                }

                // add the cell topology to the set
                for( uint Ik = 0; Ik < mListofDoubleSideSets.size(); Ik++ )
                {
                    std::string tSetName = mListOfAllSets( tCounter )->get_set_name();

                    // set the sideset cell topology and shape
                    mListOfAllSets( tCounter   )->set_IG_cell_shape( CellShape::STRAIGHT );
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
            for(moris::uint i = 0; i < mListOfAllSets.size() ; i++)
            {
                Matrix<IndexMat> const & tSetColors = mListOfAllSets(i)->get_set_colors();

                // if the max set color is greater than current max
                // replace the current max
                if(tSetColors.numel() > 0)
                {
                    moris_index tMaxSetColor = tSetColors.max();
                    if(mMaxColor < tMaxSetColor)
                    {
                        mMaxColor = tMaxSetColor;
                    }
                }
            }

            // clear old data
            mBlockSetToColor.clear();
            mSideSetToColor.clear();
            mDoubleSideSetToColor.clear();
            mAllSetToColor.clear();

            // size outer cell size of member data
            mBlockSetToColor.resize(mMaxColor + 1);
            mSideSetToColor.resize(mMaxColor + 1);
            mDoubleSideSetToColor.resize(mMaxColor + 1);
            mAllSetToColor.resize(mMaxColor + 1);

            // iterate through block sets 
            for(moris::uint i = 0; i < mListofBlocks.size(); i++)
            {
                Matrix<IndexMat> const & tSetColors = mListofBlocks(i)->get_set_colors();

                // iterate through the colors and add to related color grouping
                for(moris::uint j = 0; j < tSetColors.numel(); j++)
                {
                    mBlockSetToColor(tSetColors(j)).push_back(mListofBlocks(i));
                    mAllSetToColor(tSetColors(j)).push_back(mListofBlocks(i));
                }
            }

            // iterate through side sets 
            for(moris::uint i = 0; i < mListofSideSets.size(); i++)
            {
                Matrix<IndexMat> const & tSetColors = mListofSideSets(i)->get_set_colors();

                // iterate through the colors and add to related color grouping
                for(moris::uint j = 0; j < tSetColors.numel(); j++)
                {
                    mSideSetToColor(tSetColors(j)).push_back(mListofSideSets(i));
                    mAllSetToColor(tSetColors(j)).push_back(mListofSideSets(i));
                }
            }

            // iterate through double side sets 
            for(moris::uint i = 0; i < mListofDoubleSideSets.size(); i++)
            {
                Matrix<IndexMat> const & tSetColors = mListofDoubleSideSets(i)->get_set_colors();

                // iterate through the colors and add to related color grouping
                for(moris::uint j = 0; j < tSetColors.numel(); j++)
                {
                    mDoubleSideSetToColor(tSetColors(j)).push_back(mListofDoubleSideSets(i));
                    mAllSetToColor(tSetColors(j)).push_back(mListofDoubleSideSets(i));
                }
            }
        }

        // ----------------------------------------------------------------------------
        void
        Integration_Mesh::add_double_sided_cluster(mtk::Double_Side_Cluster* aDblSidedCluster)
        {
            // add double sided cluster to the double sided cluster list
            mDoubleSideClusters.push_back(aDblSidedCluster);
        }

        // ----------------------------------------------------------------------------

        void
        Integration_Mesh::save_MPC_to_hdf5()
        {
            // get number of integration vertices
            uint tNumVertices = this->get_num_nodes();

            std::cout<<"Num Nodes: "<<tNumVertices<<std::endl;

            Matrix< IdMat >tIdentifierMat( tNumVertices, 1, -1);

            Matrix< IdMat >  tBSToIPMap( this->get_num_basis_functions( 0 ), 1, gNoID );
            Matrix< IdMat >  tIPToIGMap( tNumVertices, 1, gNoID );
            moris::Cell< Matrix< IdMat > >  tBSToIPIds( tNumVertices );
            moris::Cell< Matrix< DDRMat > >  tBSToIPWeights( tNumVertices );

            moris::Cell< Matrix< IdMat > >  tNodeDOFs( tNumVertices );
            moris::Cell< Matrix< DDRMat > >  tMPCs( tNumVertices );

            this->create_MPC_maps( tBSToIPMap, tIPToIGMap );

            this->build_hanging_node_MPC( tBSToIPIds, tBSToIPWeights, tBSToIPMap, tIPToIGMap );

            // loop over all bulk sets
            for( uint Ik = 0; Ik < this->get_num_blocks(); Ik ++ )
            {
                // get num clusters on set
                uint tNumClustersOnSet = this->get_set_by_index( Ik )->get_num_clusters_on_set();

                if( tNumClustersOnSet > 0 )
                {
                    //Geometry_Type tGeoType = aCell->get_geometry_type();
                    //Interpolation_Order tIPOrder = aCell->get_interpolation_order();

                    // creating interpolation rule
                    Interpolation_Rule tInterpolationRule(mtk::Geometry_Type::HEX,
                            Interpolation_Type::LAGRANGE,
                            mtk::Interpolation_Order::LINEAR,
                            Interpolation_Type::LAGRANGE,
                            Interpolation_Order::LINEAR);

                    // create a space interpolator
                    Space_Interpolator tSpaceInterpolator(tInterpolationRule);

                    // loop over all clusters in set
                    for( uint Ii = 0; Ii < tNumClustersOnSet; Ii++ )
                    {
                        const Cluster * tCluster = this->get_set_by_index( Ik )->get_clusters_by_index( Ii );

                        moris::Matrix<moris::DDRMat> tLocalCoords = tCluster->get_vertices_local_coordinates_wrt_interp_cell();

                        moris::Cell< moris::mtk::Vertex const * > tIGVertices = tCluster->get_vertices_in_cluster();

                        // get number of vertices on the treated mesh cluster
                        uint tNumNodes = tLocalCoords.n_rows();

                        const moris::mtk::Cell & tInterpolationCell = tCluster->get_interpolation_cell();

                        moris::Cell< Vertex *> tIPVertices = tInterpolationCell.get_vertex_pointers();

                        Matrix< IdMat > tIdentifierMat( tIGVertices.size(), 1, 0 );



                        for( uint Iv = 0; Iv < tIPVertices.size(); Iv++ )
                        {
                            if( tIPVertices( Iv )->get_interpolation( 0 )->get_ids().numel() > 1 )
                            {
                                tIdentifierMat( Iv ) = 1;
                            }
                        }

                        Matrix< IndexMat > tHangingNodeIdentifier(tNumVertices, 1, 0);
                        if( !tCluster->is_trivial() )
                        {
                            Matrix< IndexMat > tHangingNodeIndices = tCluster->get_hanging_nodes();

                            print( tHangingNodeIndices,"tHangingNodeIndices");

                            for( uint Iv = 0; Iv < tHangingNodeIndices.numel(); Iv++ )
                            {
                                tHangingNodeIdentifier( tHangingNodeIndices( Iv ) ) = 1;
                            }
                        }


                        // epsilon to count T-Matrix
                        real tEpsilon = 1e-12;

                        // loop over the vertices on the treated mesh cluster
                        for( uint Iv = 0; Iv < tNumNodes; Iv++ )
                        {
                            tSpaceInterpolator.set_space_time( trans( tLocalCoords.get_row( Iv ) ) );

                            const Matrix< DDRMat > & tN = tSpaceInterpolator.NXi();

                            moris_id tIGVertexId = tIGVertices( Iv )->get_id();
                            moris_id tIGVertexIndex = tIGVertices( Iv )->get_index();

                            uint tNCols = tN.n_cols();

                            tMPCs( tIGVertexIndex ).set_size( tNCols, 1, -1.0 );
                            tNodeDOFs( tIGVertexIndex ).set_size( tNCols, 1, gNoID );

                            std::cout<<"tIGVertexId: "<<tIGVertexId<<" "<<tIGVertexIndex<<std::endl;

                            if( tHangingNodeIdentifier( tIGVertexIndex ) == 1 )
                            {
                                // fixme add hanging - hanging node

                                uint tCount = 0;
                                // loop over all nonzero entries
                                for( uint i=0; i<tNCols; i++ )
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
                                if( tIdentifierMat( Iv ) == 1 )
                                {
                                    for( uint Ib = 0; Ib < tBSToIPWeights( tIGVertexId ).numel(); Ib++ )
                                    {
                                        tMPCs( tIGVertexIndex )(Ib) = tBSToIPWeights( tIGVertexId )(Ib);
                                        tNodeDOFs( tIGVertexIndex )(Ib) = tBSToIPIds( tIGVertexId )(Ib);
                                    }

                                }
                                else
                                {
                                    tMPCs( tIGVertexIndex )( 0 ) = 1.0;
                                    tNodeDOFs( tIGVertexIndex )( 0 ) = tIGVertexId;
                                }
//                                tMPCs( tIGVertexIndex ).resize( 1, 1 );
//                                tNodeDOFs( tIGVertexIndex ).resize( 1,1 );
                            }
                        }
                    }
                }
            }

            print(tNodeDOFs,"tNodeDOFs");
            print(tMPCs,"tMPCs");

            Matrix< IdMat >  tIDs        ( 1, tNumVertices );
            Matrix< IdMat >  tMPCIDs     ( tMPCs( 0 ).numel(), tNumVertices );
            Matrix< DDRMat > tMPCWeightss( tMPCs( 0 ).numel(), tNumVertices );

            for( uint Ik = 0; Ik < tNumVertices; Ik++ )
            {
                tIDs( Ik ) = Ik + 1;

                for( uint Ii = 0; Ii < tMPCs( 0 ).numel(); Ii++ )
                {
                    tMPCIDs( Ii, Ik )      = tNodeDOFs( Ik )( Ii );
                    tMPCWeightss( Ii, Ik ) = tMPCs( Ik )( Ii );
                }
            }

//            // add order to path
//            std::string tFilePath =    aFilePath.substr(0,aFilePath.find_last_of(".")) // base path
//                                                                                                              + "_" + std::to_string( tMesh->get_order() ) // rank of this processor
//            +  aFilePath.substr( aFilePath.find_last_of("."), aFilePath.length() );

            std::string tFilePath = "MPC.hdf5";

            // make path parallel
            tFilePath = parallelize_path( tFilePath );

            // Create a new file using default properties
            hid_t tFileID = H5Fcreate( tFilePath.c_str(),
                    H5F_ACC_TRUNC,
                    H5P_DEFAULT,
                    H5P_DEFAULT);

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
                Matrix< IdMat > & tBSToIPMap,
                Matrix< IdMat > & tIPToIGMap)
        {
            Matrix< DDSMat > tIdentifierMatBS( tBSToIPMap.numel(),1 ,0 );
            Matrix< DDSMat > tIdentifierMatIP( tIPToIGMap.numel(),1 ,0 );

            // loop over all bulk sets
            for( uint Ik = 0; Ik < this->get_num_blocks(); Ik ++ )
            {
                // get num clusters on set
                uint tNumClustersOnSet = this->get_set_by_index( Ik )->get_num_clusters_on_set();

                if( tNumClustersOnSet > 0 )
                {
                    // loop over all clusters in set
                    for( uint Ii = 0; Ii < tNumClustersOnSet; Ii++ )
                    {
                        const Cluster * tCluster = this->get_set_by_index( Ik )->get_clusters_by_index( Ii );

                        const moris::mtk::Cell & tInterpolationCell = tCluster->get_interpolation_cell();

                        moris::Cell< Vertex *> tIPVertices = tInterpolationCell.get_vertex_pointers();

                        moris::Cell< moris::mtk::Vertex const * > tIGVertices = tCluster->get_vertices_in_cluster();

                        std::cout<<" IGVert: "<<tIGVertices.size()<<" IPVert: "<<tIPVertices.size()<<std::endl;

                        for( uint Ij = 0; Ij < tIPVertices.size(); Ij++)
                        {
                            MORIS_ASSERT( tIPVertices( Ij )->has_interpolation( 0 ), "Vertex does not have interpolation");

                            uint tIP_ID = tIPVertices( Ij )->get_id();

                            // build BSpline to iP map
                            if( tIPVertices( Ij )->get_interpolation( 0 )->get_ids().numel() == 1 )
                            {
                                moris_id tBSId = tIPVertices( Ij )->get_interpolation( 0 )->get_ids()( 0 );

                                if( tIdentifierMatBS( tBSId ) == 0)
                                {
                                    tBSToIPMap( tBSId ) = tIP_ID;

                                    tIdentifierMatBS( tBSId ) = 1;
                                }
                            }

                            // build IP to IG map
                            // assume corner vertices are first 8
                            uint tIG_ID = tIGVertices( Ij )->get_id();

                            if( tIdentifierMatIP( tIP_ID ) == 0 )
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
                moris::Cell< Matrix< IdMat > > & tBSToIPIds,
                moris::Cell< Matrix< DDRMat > > & tBSToIPWeights,
                const Matrix< IdMat > & tBSToIPMap,
                const Matrix< IdMat > & tIPToIGMap)
        {
            // loop over all bulk sets
            for( uint Ik = 0; Ik < this->get_num_blocks(); Ik ++ )
            {
                // get num clusters on set
                uint tNumClustersOnSet = this->get_set_by_index( Ik )->get_num_clusters_on_set();

                if( tNumClustersOnSet > 0 )
                {
                    // loop over all clusters in set
                    for( uint Ii = 0; Ii < tNumClustersOnSet; Ii++ )
                    {
                        const Cluster * tCluster = this->get_set_by_index( Ik )->get_clusters_by_index( Ii );

                        const moris::mtk::Cell & tInterpolationCell = tCluster->get_interpolation_cell();

                        moris::Cell< Vertex *> tIPVertices = tInterpolationCell.get_vertex_pointers();

                        moris::Cell< moris::mtk::Vertex const * > tIGVertices = tCluster->get_vertices_in_cluster();

                        for( uint Ij = 0; Ij < tIPVertices.size(); Ij++)
                        {
                            // build BSpline to iP map
                            if( tIPVertices( Ij )->get_interpolation( 0 )->get_ids().numel() > 1 )
                            {
                                uint tIG_ID = tIGVertices( Ij )->get_id();

                                Matrix< IdMat > tBSId = tIPVertices( Ij )->get_interpolation( 0 )->get_ids();
                                Matrix< DDRMat > tBSWeights= *tIPVertices( Ij )->get_interpolation( 0 )->get_weights();

                                tBSToIPIds( tIG_ID ).resize( tBSId.numel(), 1 );
                                tBSToIPWeights( tIG_ID ).resize( tBSId.numel(), 1 );

                                for( uint Ia = 0; Ia < tBSId.numel(); Ia++ )
                                {
                                    tBSToIPWeights( tIG_ID )( Ia ) = tBSWeights( Ia );
                                    tBSToIPIds    ( tIG_ID )( Ia ) = tIPToIGMap( tBSToIPMap( tBSId( Ia ) ) );
                                }
                            }
                        }
                    }
                }
            }
        }


    }
}

