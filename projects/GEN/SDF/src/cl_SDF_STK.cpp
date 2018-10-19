#include "cl_SDF_STK.hpp"

#include "cl_Stopwatch.hpp"
#include "cl_Mesh_Enums.hpp"
#include "stk_impl/cl_MTK_Mesh_STK.hpp"
#include "SDF_Tools.hpp"

namespace moris
{
    namespace sdf
    {
//-------------------------------------------------------------------------------

        void
        STK::create_mesh_data(
                moris::Cell< Matrix<DDRMat> > & aFields,
                moris::Cell< std::string >    & aFieldLabels,
                const double                    aTimeStep )
        {
            // start timer
            tic tTimer;

            // get pointer to mtk mesh
            mtk::Mesh * tMesh = mMesh.get_mtk_mesh();

            // get number of nodes
            moris_index tNumberOfNodes = tMesh->get_num_nodes();

            // get number of elements
            moris_index tNumberOfElements = tMesh->get_num_elems();


            // get number of nodes per element
            uint tNumberOfNodesPerElement
                = tMesh->get_mtk_cell( 0 ).get_number_of_vertices();

            // initialize topology field
            mElementTopology.set_size( tNumberOfElements, tNumberOfNodesPerElement );

            // initialize local to global map
            mElementLocalToGlobal.set_size( tNumberOfElements, 1 );

            // initialize coordinates
            mNodeCoords.set_size( tNumberOfNodes, mNumberOfDimensions );

            // initialize local to global map for nodes
            mNodeLocalToGlobal.set_size( tNumberOfNodes, 1 );

            // initialize node ownership
            mNodeOwner.set_size( tNumberOfNodes, 1 );

            // reset field info
            mFieldsInfo.FieldsName.clear();
            mFieldsInfo.FieldsRank.clear();

            // get number of fields
            uint tNumberOfFields = aFields.size();

            for( uint k=0; k<tNumberOfFields; ++k )
            {
                mFieldsInfo.FieldsName.push_back( aFieldLabels( k ) );
                mFieldsInfo.FieldsRank.push_back( EntityRank::NODE );
            }

            uint tElementIDs = tNumberOfFields;
            uint tNodeIDs = tNumberOfFields + 1;
            Matrix< DDRMat > tEmpty;

            aFields.push_back( tEmpty );
            aFields.push_back( tEmpty );
            mFieldsInfo.FieldsName.push_back( "Cell_ID" );
            mFieldsInfo.FieldsName.push_back( "Vertex_ID" );

            mFieldsInfo.FieldsRank.push_back( EntityRank::ELEMENT );
            mFieldsInfo.FieldsRank.push_back( EntityRank::NODE );

            aFields( tElementIDs ).set_size( tNumberOfElements, 1 );
            aFields( tNodeIDs ).set_size( tNumberOfNodes, 1 );

            // loop over all elements
            for( moris_index e=0; e<tNumberOfElements; ++e )
            {

                // get node IDs
                auto tNodeIndices = tMesh->get_nodes_connected_to_element_loc_inds( e );

                // copy node IDs to topology matrix
                for( uint k=0; k<tNumberOfNodesPerElement; ++k )
                {
                    mElementTopology( e, k ) = tMesh->get_glb_entity_id_from_entity_loc_index(
                            tNodeIndices( k ), EntityRank::NODE );
                }

                // save element index in map
                mElementLocalToGlobal( e ) = tMesh->get_glb_entity_id_from_entity_loc_index(
                        e, EntityRank::ELEMENT );

                // converte entry to real
                aFields( tElementIDs )( e ) = mElementLocalToGlobal( e );
            }

            // loop over all nodes
            for( moris_index k=0; k<tNumberOfNodes; ++k )
            {

                Matrix< DDRMat > tNodeCoords = mMesh.get_node_coordinate( k );

                // copy coords to output matrix
                for( uint i=0; i<3; ++i )
                {
                    mNodeCoords( k, i ) = tNodeCoords( i );
                }

                // copy node index into map
                mNodeLocalToGlobal( k ) =
                        tMesh->get_glb_entity_id_from_entity_loc_index(
                        k, EntityRank::NODE );

                mNodeOwner( k ) = tMesh->get_entity_owner( k, EntityRank::NODE );

                // converte entry to real
                aFields( tNodeIDs )( k ) = mNodeLocalToGlobal( k );
            }

            // link mesh data object

            mMeshData.SpatialDim              = & mNumberOfDimensions;
            mMeshData.ElemConn(0)             = & mElementTopology;
            mMeshData.NodeCoords              = & mNodeCoords;
            mMeshData.EntProcOwner            = & mNodeOwner;
            mMeshData.LocaltoGlobalElemMap(0) = & mElementLocalToGlobal;
            mMeshData.LocaltoGlobalNodeMap    = & mNodeLocalToGlobal;
            mMeshData.FieldsInfo              = & mFieldsInfo;
            mFieldsInfo.FieldsData            = & aFields;

            // set timestep of mesh data object
            mMeshData.TimeStamp = aTimeStep;

            /*if( par_rank() == 1 )
            {
                std::cout << "Number of Fields" << tNumberOfFields << std::endl;
                for( uint i=0; i<tNumberOfFields; ++i )
                {
                    std::cout << mFieldsInfo.FieldsName( i ) << " " << aFields( i ).length() << std::endl;
                }

                print( mElementTopology, "topo" );
                print( mElementLocalToGlobal, "Elements" );
                print( mNodeLocalToGlobal, "Nodes" );
                print( mNodeCoords, "Coords" );
                print( mNodeOwner, "owner" );


            } */
            if ( mMesh.is_verbose() )
            {
                // stop timer
                real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;

                // print output
                std::fprintf( stdout,"%s Created MTK output object.\n               Mesh has %lu elements and %lu nodes.\n               Creation took %5.3f seconds.\n\n",
                        proc_string().c_str(),
                        ( long unsigned int ) tNumberOfElements,
                        ( long unsigned int ) tNumberOfNodes,
                        ( double ) tElapsedTime / 1000);
            }
        }

//-------------------------------------------------------------------------------

        void
        STK::save_to_file( const std::string & aFilePath )
        {
            tic tTimer;

            // create database object
            moris::mtk::Mesh_STK tMesh( mMeshData );

            // copy file path, since tMesh does not like const input
            std::string tFilePath = aFilePath;

            // save file
            tMesh.create_output_mesh( tFilePath );

            if ( mMesh.is_verbose() )
            {
                // stop timer
                real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;

                // print output
                std::fprintf( stdout,"%s Wrote MTK mesh to file.\n               Writing took %5.3f seconds.\n\n",
                        proc_string().c_str(),
                        ( double ) tElapsedTime / 1000 );
            }
        }

//-------------------------------------------------------------------------------
    }
}
