//
// Created by christopherson on 9/13/19.
//

#include <exodusII.h>
#include "cl_MTK_Reader_Exodus.hpp"
#include "cl_MTK_Mesh_Core.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"

// TODO
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Data_STK.hpp"
#include "cl_MTK_Mesh_Core_STK.hpp"
#include "cl_MTK_Interpolation_Mesh_STK.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Mesh_Tools.hpp"
#include "cl_MTK_Integration_Mesh_STK.hpp"
#include "cl_MTK_Sets_Info.hpp"

#include <iostream>

Reader_Exodus::Reader_Exodus()
{
    mExoid = 0;
}

Reader_Exodus::~Reader_Exodus()
= default;

void Reader_Exodus::set_error_options(bool abort, bool debug, bool verbose)
{
    ex_opts(abort * EX_ABORT | debug * EX_DEBUG | verbose * EX_VERBOSE);
}

void Reader_Exodus::open_file(std::string aExodusFileName, float aVersion)
{
    int tCPUWordSize = 8, tIOWordSize = 8; // TODO
    mExoid = ex_open(aExodusFileName.c_str(), EX_READ, &tCPUWordSize, &tIOWordSize, &aVersion);
}

void Reader_Exodus::close_file(bool aRename)
{
    ex_close(mExoid);
}

void Reader_Exodus::read_file(std::string aFileName)
{
    std::cout << "Start" << std::endl;
    Reader_Exodus::open_file(aFileName);
    std::cout << "Open" << std::endl;

    // Get init from Exodus
    char tDatabaseName[100] = ""; // FIXME
    moris::uint tNumNodes, tNumElements, tNumBlocks, tNumNodeSets, tNumSideSets;
    ex_get_init(mExoid, tDatabaseName, &mNumSpatialDimensions, &tNumNodes, &tNumElements, &tNumBlocks, &tNumNodeSets,
                &tNumSideSets);

    // Create and population mesh data input
    moris::mtk::MtkMeshData mMeshDataInput(tNumBlocks);
    mMeshDataInput.SpatialDim = &mNumSpatialDimensions;
    mMeshDataInput.MarkNoBlockForIO = false;
    mMeshDataInput.CreateAllEdgesAndFaces = true;
    mMeshDataInput.Verbose = true;

    std::cout << "Initialize" << std::endl;

    // Read the coordinates from the exodus file
    moris::Matrix<moris::DDRMat> tXCoordinates(tNumNodes, 1, 0.0);
    moris::Matrix<moris::DDRMat> tYCoordinates(tNumNodes, 1, 0.0);
    moris::Matrix<moris::DDRMat> tZCoordinates(tNumNodes, 1, 0.0);
    ex_get_coord(mExoid, tXCoordinates.data(), tYCoordinates.data(), tZCoordinates.data());

    // Combine the coordinates into one matrix for input
    mNumSpatialDimensions = 3; // FIXME
    mNodeCoordinates.resize(tNumNodes, mNumSpatialDimensions);
    mNodeCoordinates.set_column(0, tXCoordinates);
    mNodeCoordinates.set_column(1, tYCoordinates);
    mNodeCoordinates.set_column(2, tZCoordinates);
    mMeshDataInput.NodeCoords = &mNodeCoordinates;

    // Get the node map
    mNodeMap.resize(tNumNodes, 1);
    ex_get_id_map(mExoid, EX_NODE_MAP, mNodeMap.data());
    mMeshDataInput.LocaltoGlobalNodeMap = &mNodeMap;

    // Node ownership
    mNodeOwner.set_size(1, tNumNodes, moris::par_rank());
    mMeshDataInput.NodeProcOwner = &mNodeOwner;

    std::cout << "Nodes" << std::endl;

    // Prepare for get_block
    mBlockSetInfo.resize(tNumBlocks, moris::mtk::MtkBlockSetInfo());
    mBlockDescription.resize(tNumBlocks, "");
    mElemConn.resize(tNumBlocks, moris::Matrix<moris::IndexMat>(1, 1));
    mLocaltoGlobalElemMap.resize(tNumBlocks, moris::Matrix<moris::IdMat>(1, 1));
    moris::uint tNumElementsInBlock, tNumNodesPerElement, tNumEdgesPerElement, tNumFacesPerElement, tNumAttributesPerElement;
    char tReadName[100];

    // Loop through all of the blocks
    for (moris::uint tBlockIndex = 0; tBlockIndex < tNumBlocks; tBlockIndex++)
    {
        // Get the block parameters
        std::cout << "Reading Block..." << std::endl;
        ex_get_block(mExoid, EX_ELEM_BLOCK, tBlockIndex + 1, tReadName, &tNumElementsInBlock,
                     &tNumNodesPerElement, &tNumEdgesPerElement, &tNumFacesPerElement, &tNumAttributesPerElement);
        std::cout << tNumElementsInBlock << ", " << tNumNodesPerElement << ", " << tNumEdgesPerElement << ", " << tNumFacesPerElement << ", " <<std::endl;

        // Get connectivity
        moris::Matrix<moris::IndexMat> tBlockNodesVector(tNumElementsInBlock * tNumNodesPerElement, 1, 0);
        ex_get_conn(mExoid, EX_ELEM_BLOCK, tBlockIndex + 1, tBlockNodesVector.data(), nullptr, nullptr);

        // Reorganize connectivity
        mElemConn(tBlockIndex).resize(tNumElementsInBlock, tNumNodesPerElement);
        moris::uint tVectorIndex = 0;
        for (moris::uint tElementIndex = 0; tElementIndex < tNumElementsInBlock; tElementIndex++)
        {
            for (moris::uint tNodeIndex = 0; tNodeIndex < tNumNodesPerElement; tNodeIndex++)
            {
                mElemConn(tBlockIndex)(tElementIndex, tNodeIndex) = tBlockNodesVector(tVectorIndex);
                tVectorIndex++;
            }
        }
        mMeshDataInput.ElemConn(tBlockIndex) = &mElemConn(tBlockIndex);

        // Element map
        mLocaltoGlobalElemMap(tBlockIndex).set_size(tNumElementsInBlock, 1);
        ex_get_map(mExoid, mLocaltoGlobalElemMap(tBlockIndex).data());
        mMeshDataInput.LocaltoGlobalElemMap(tBlockIndex) = &mLocaltoGlobalElemMap(tBlockIndex);

        // Cell ids
        mBlockSetInfo(tBlockIndex).mCellIdsInSet = mMeshDataInput.LocaltoGlobalElemMap(tBlockIndex);

        // Block name
        mBlockDescription(tBlockIndex).assign(tReadName, 100);
        mBlockSetInfo(tBlockIndex).mBlockSetName = mBlockDescription(tBlockIndex);

        // Cell topology
        mBlockSetInfo(tBlockIndex).mBlockSetTopo = this->get_cell_topology(tNumNodesPerElement);



        moris::print(*mBlockSetInfo(0).mCellIdsInSet, "*mBlockSetInfo(0).mCellIdsInSet");
        std::cout << "mBlockSetInfo(0).mBlockSetName: " << mBlockSetInfo(0).mBlockSetName << std::endl;
        std::cout << "mBlockSetInfo(0).mBlockSetTopo: " << static_cast<int>(mBlockSetInfo(0).mBlockSetTopo) << std::endl;

        // Add sets to input data
        mMtkMeshSets.add_block_set(&mBlockSetInfo(tBlockIndex));


        std::cout << "End Block" << std::endl;
        moris::print(*mMeshDataInput.ElemConn(0), "mMeshDataInput.ElemConn(0)");
    }
    mMeshDataInput.SetsInfo = &mMtkMeshSets;
    moris::print(*mMeshDataInput.ElemConn(0), "mMeshDataInput.ElemConn(0)");

    std::cout << "Blocks Done" << std::endl;



    std::cout << "ElemConn: "               << mMeshDataInput.ElemConn(0)             << std::endl;
    std::cout << "LocaltoGlobalElemMap: "   << mMeshDataInput.LocaltoGlobalElemMap(0) << std::endl;
    std::cout << "CreateAllEdgesAndFaces: " << mMeshDataInput.CreateAllEdgesAndFaces  << std::endl;
    std::cout << "Verbose: "                << mMeshDataInput.Verbose                 << std::endl;
    std::cout << "SpatialDim: "             << mMeshDataInput.SpatialDim              << std::endl;
    std::cout << "NodeCoords: "             << mMeshDataInput.NodeCoords              << std::endl;
    std::cout << "NodeProcOWner: "          << mMeshDataInput.NodeProcOwner           << std::endl;
    std::cout << "LocaltoGlobalNodeMap: "   << mMeshDataInput.LocaltoGlobalNodeMap    << std::endl;
    std::cout << "SetsInfo: "               << mMeshDataInput.SetsInfo                << std::endl;
    std::cout << "MarkNoBlockForIO: "       << mMeshDataInput.MarkNoBlockForIO        << std::endl;

//    moris::print(*mMeshDataInput.ElemConn(0), "ElemConn");
//    moris::print(*mMeshDataInput.LocaltoGlobalElemMap(0), "LocaltoGlobalElemMap");
    moris::print(*mMeshDataInput.NodeCoords, "NodeCoords");
    moris::print(*mMeshDataInput.NodeProcOwner, "NodeProcOwner");
    moris::print(*mMeshDataInput.LocaltoGlobalNodeMap, "LocaltoGlobalNodeMap");
//    moris::print(*mMeshDataInput.SetsInfo, "SetsInfo");

    mMesh = moris::mtk::create_integration_mesh(MeshType::STK, mMeshDataInput);

    Reader_Exodus::close_file(false);
    std::cout << "Close" << std::endl;
}








// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------

CellTopology Reader_Exodus::get_cell_topology(int aNumNodesPerElement)
{
    switch (aNumNodesPerElement)
    {
        case 3:
            return CellTopology::TRI3;
        case 4:
            return CellTopology::TET4; // TODO QUAD4?
        case 6:
            return CellTopology::PRISM6;
        case 8:
            return CellTopology::HEX8;
        case 10:
            return CellTopology::TET10;
        default:
            MORIS_ERROR(0, "This element is invalid or it hasn't been implemented yet!");
            return CellTopology::INVALID;
    }
}
