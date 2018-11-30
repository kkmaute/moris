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

            // - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 1: get number of nodes, elements etc
            // - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // get number of nodes
            uint tNumberOfNodes = tMesh->get_num_nodes();

            // get number of elements
            uint tNumberOfElements = tMesh->get_num_elems();


            // get number of nodes per element
            uint tNumberOfNodesPerElement
                = tMesh->get_mtk_cell( 0 ).get_number_of_vertices();

            // - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 2: Allocate Matrices
            // - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // initialize topology field
            mElementTopology.set_size( tNumberOfElements, tNumberOfNodesPerElement );

            // initialize local to global map
            mElementLocalToGlobal.set_size( tNumberOfElements, 1 );

            // initialize coordinates
            mNodeCoords.set_size( tNumberOfNodes, mNumberOfDimensions );

            // initialize local to global map for nodes
            mNodeLocalToGlobal.set_size( tNumberOfNodes, 1 );

            mNodeOwner.set_size( tNumberOfNodes, 1 );

            // get number of fields
            uint tNumberOfFields = aFields.size();

            // add node ids to fields
            Matrix< DDRMat > tEmpty;
            aFields.push_back( tEmpty );
            Matrix< DDRMat > & tNodeIDs = aFields( tNumberOfFields );
            tNodeIDs.set_size( tNumberOfNodes, 1 );

            // add cell ids to fields
            aFields.push_back( tEmpty );
            Matrix< DDRMat > & tCellIDs = aFields( tNumberOfFields+1 );
            tCellIDs.set_size( tNumberOfElements, 1 );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 3: Populate Matrices
            // - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // loop over all elements
            for( uint e=0; e<tNumberOfElements; ++e )
            {
                mtk::Cell & tCell = tMesh->get_mtk_cell( e );

                //Matrix< IndexMat > tNodeInds = tCell.get_vertex_inds();
                Matrix< IndexMat > tNodeIDs = tMesh->get_entity_connected_to_entity_glob_ids(
                        tCell.get_id(), EntityRank::ELEMENT, EntityRank::NODE );

                // cast copy node IDs to topology matrix
                for( uint k=0; k<tNumberOfNodesPerElement; ++k )
                {
                    mElementTopology( e, k ) = tNodeIDs( k );
                }

                // save element index in map
                mElementLocalToGlobal( e ) = tCell.get_id();

                // also save element in export array
                tCellIDs( e ) = mElementLocalToGlobal( e );
            }

            // loop over all nodes
            for( uint k=0; k<tNumberOfNodes; ++k )
            {
                mtk::Vertex & tVertex = tMesh->get_mtk_vertex( k );

                Matrix< DDRMat > tNodeCoords =  tVertex.get_coords();

                // copy coords to output matrix
                for( uint i=0; i<3; ++i )
                {
                    mNodeCoords( k, i ) = tNodeCoords( i );
                }

                // copy node Owner
                mNodeOwner( k ) = tVertex.get_owner();

                // copy node index into map
                mNodeLocalToGlobal( k ) = tVertex->get_id();

                // save vertex id
                tNodeIDs( k ) =  mNodeLocalToGlobal( k );
            }

            // - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 4: Link Fields
            // - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // Initialize scalar field data
            mFields = moris::Cell< mtk::Scalar_Field_Info<DDRMat> >( tNumberOfFields+2 );

            for( uint f=0; f<tNumberOfFields; ++f )
            {
                mFields( f ).set_field_name( aFieldLabels( f ) );
                mFields( f ).set_field_entity_rank( EntityRank::NODE );

                mFields( f ).add_field_data(
                        & mNodeLocalToGlobal,
                        & aFields( f ) );
            }

            uint tNodeFieldIndex = tNumberOfFields;

            // Add information about node id field
            mFields(tNodeFieldIndex).set_field_name( "Vertex_ID" );
            mFields(tNodeFieldIndex).set_field_entity_rank( EntityRank::NODE );
            mFields(tNodeFieldIndex).add_field_data( &mNodeLocalToGlobal, &aFields( tNodeFieldIndex ));
            mFieldsInfo.mRealScalarFields.push_back( &mFields(tNodeFieldIndex) );

            uint tElementFieldIndex = tNumberOfFields+1;

            // Add information about cell id field
            mFields(tElementFieldIndex).set_field_name( "Cell_ID" );
            mFields(tElementFieldIndex).set_field_entity_rank( EntityRank::ELEMENT );
            mFields(tElementFieldIndex).add_field_data( &mElementLocalToGlobal, &aFields( tElementFieldIndex ) );
            mFieldsInfo.mRealScalarFields.push_back( &mFields(tElementFieldIndex) );

            // link mesh data object
            mMeshData.SpatialDim                = & mNumberOfDimensions;
            mMeshData.ElemConn( 0 )             = & mElementTopology;
            mMeshData.NodeCoords                = & mNodeCoords;
            mMeshData.NodeProcOwner             = & mNodeOwner;
            mMeshData.LocaltoGlobalElemMap( 0 ) = & mElementLocalToGlobal;
            mMeshData.LocaltoGlobalNodeMap      = & mNodeLocalToGlobal;
            mMeshData.FieldsInfo                = & mFieldsInfo;
            mMeshData.SetsInfo                  = nullptr;

            // set timestep of mesh data object
            mMeshData.TimeStamp = aTimeStep;
            mMeshData.AutoAuraOptionInSTK = false;

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
                std::fprintf( stdout,"%s Wrote MTK mesh to file %s.\n               Writing took %5.3f seconds.\n\n",
                        proc_string().c_str(),
                        aFilePath.c_str(),
                        ( double ) tElapsedTime / 1000 );
            }
        }

//-------------------------------------------------------------------------------
    }
}
