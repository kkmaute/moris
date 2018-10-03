//
// Created by messe on 3/7/18.
//

#include "cl_Ge_SDF_Mesh_Data.hpp"
#include "cl_Mesh_Enums.hpp" // MTK/src
// -----------------------------------------------------------------------------

void
ge::SDF_Mesh_Data::update()
{
    // start the timer
    moris::tic tTimer;
    // nodes on proc in global index
    mNodesOnProc = mBackgroundMesh.get_entities_owned_and_shared_by_current_proc(moris::EntityRank::NODE);

    // elements on proc in global index
    mElementsOnProc = mBackgroundMesh.get_entities_owned_and_shared_by_current_proc(moris::EntityRank::ELEMENT);

    // get the number of nodes
    mNumberOfNodes = mNodesOnProc.length();

    // get number of elements
    mNumberOfElements = mElementsOnProc.length();

    // convert global node ids to local node ids
    /*mLocalNodesOnProc.set_size(mNumberOfNodes, 1);
    for ( moris::uint k=0; k<mNumberOfNodes; ++k )
    {
        mLocalNodesOnProc( k ) = mBackgroundMesh.get_local_id_on_proc_from_global_id(
                EntityRank::NODE,
                mNodesOnProc( k ) );

        //if (moris::par_rank() == 1)
        //{
        //    std:give_coordinate_from_lagrange_basis  << " " <<  mNodesOnProc( k ) << " " <<  mLocalNodesOnProc( k ) << std::endl;
        //}
    }

    // convert global element ids to local element ids
    mLocalElementsOnProc.set_size(mNumberOfElements, 1);
    for ( moris::uint k=0; k<mNumberOfElements; ++k )
    {
        mLocalElementsOnProc( k ) = mBackgroundMesh.get_local_id_on_proc_from_global_id(
                EntityRank::ELEMENT,
                mElementsOnProc( k ) );
    } */
    //moris::tic tTimer0;
    // sparse map
    //moris::Mat< moris:: uint > tNodeMap( mNodesOnProc.max()+1, 1, UINT_MAX );
    moris::map< uint, uint > tNodeMap;

    for ( moris::uint k=0; k<mNumberOfNodes; ++k )
    {
        //tNodeMap( mNodesOnProc( k ) ) = k;
        tNodeMap[ mNodesOnProc( k ) ] = k;
    }

    //moris::real tElapsedTime0
    //        = tTimer0.toc<moris::chronos::milliseconds>().wall;
    //std::fprintf(stdout, "Time for Map         : %5.3f [sec]\n",
    //             tElapsedTime0/1000);

    //moris::tic tTimer1;
    // create element topology
    moris::Matrix< moris::DDUMat > tNodes;

    mNumberOfNodesPerElement.set_size(mNumberOfElements, 1);

    // count nodes per element
    for (moris::uint k=0; k<mNumberOfElements; ++k)
    {
        mNumberOfNodesPerElement(k) = (mBackgroundMesh.get_nodes_connected_to_element(mElementsOnProc( k ))).length();
    }

    mElementTopology.set_size(mNumberOfNodesPerElement.max(), mNumberOfElements);

    //mElementTopology.clear();
    //mElementTopology.reserve(mNumberOfElements);

    for (moris::uint k=0; k<mNumberOfElements; ++k)
    {
        // get local node IDs connected to this element
        //moris::Matrix< moris::DDUMat > tNodesOfElement = mBackgroundMesh.get_local_nodes_on_proc_connected_to_local_element_on_proc(
        //        mLocalElementsOnProc( k ));

        moris::Matrix< moris::DDUMat > tNodesOfElement = mBackgroundMesh.get_nodes_connected_to_element(
                mElementsOnProc( k ));
        //std::cout << " Element " << k << " ID: " << mElementsOnProc( k ) << " dofs:" << tNodesOfElement.length() << std::endl;

         for ( moris::uint i=0; i<mNumberOfNodesPerElement(k); ++i )
        {
            // save element topology
            mElementTopology(i, k) = tNodeMap.find( tNodesOfElement( i ) );
        }

    }

    //std::cout << "SDF index of node 37380 " << tNodeMap.find(37380) << std::endl;

    //moris::real tElapsedTime1
    //        = tTimer1.toc<moris::chronos::milliseconds>().wall;

    //std::fprintf(stdout, "Time for Element Topology         : %5.3f [sec]\n",
    //              tElapsedTime1/1000);

    // allocate the matrix containing the node coordinates
    mNodeCoords.set_size(3, mNumberOfNodes);

    // temporary vector containing node coordinate
    moris::Matrix< moris::DDRMat > tNodeCoordinate(3,1);

    // initialize minimum and maximum coordinates with extreme values
    for (moris::uint i=0; i<3; ++i)
    {
        mMinNodeCoordinate(i) =  MORIS_GE_HUGE;
        mMaxNodeCoordinate(i) = -MORIS_GE_HUGE;
    }

    // loop over all nodes on proc
    for (moris::uint k = 0; k < mNumberOfNodes; ++k)
    {
        // get the coordinates of this node
        moris::Matrix< moris::DDRMat > tNodeCoordinate
        = mBackgroundMesh.get_selected_nodes_coords(mNodesOnProc.get_row(k));

        for (moris::uint i = 0; i < 3; ++i)
        {
            // save the coordinate in the matrix
            mNodeCoords(i, k) = tNodeCoordinate(i);

            // update minimal node coordinate
            mMinNodeCoordinate(i) = ge::min(mMinNodeCoordinate(i),
                    tNodeCoordinate( i ));

            // update maximal node coordinate
            mMaxNodeCoordinate(i) = ge::max(mMaxNodeCoordinate(i),
                    tNodeCoordinate( i ));
        }
    }

    // stop the timer
    // print elapsed time
    moris::real tElapsedTime
        = tTimer.toc<moris::chronos::milliseconds>().wall;

    if(moris::par_size() == 1)
    {
        std::fprintf(stdout, "Time for mesh update           : %5.3f [sec]\n",
                tElapsedTime/1000 );
    }
    else
    {
        std::fprintf(stdout, "Proc % i - Time for mesh update           : %5.3f [sec]\n",
                (int) moris::par_rank(), tElapsedTime/1000 );
    }
    // create output message
    std::fprintf(stdout,  "SDF Mesh: %i nodes and %i elements\n", (int) mNumberOfNodes, (int) mNumberOfElements );
}

// -----------------------------------------------------------------------------
