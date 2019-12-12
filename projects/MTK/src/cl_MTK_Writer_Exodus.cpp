//
// Created by christopherson on 9/13/19.
//

#include <exodusII.h>
#include "cl_MTK_Writer_Exodus.hpp"
#include "cl_MTK_Mesh_Core.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Data_STK.hpp"
#include "cl_MTK_Mesh_Core_STK.hpp"

#include <iostream>

Writer_Exodus::Writer_Exodus(moris::mtk::Mesh* aMeshPointer)
    : mMesh(aMeshPointer)
{
    mExoid = 0;
    this->set_error_options(true, true, true);
}

//----------------------------------------------------------------------------------------------------------------------

Writer_Exodus::~Writer_Exodus()
= default;

void Writer_Exodus::set_error_options(bool abort, bool debug, bool verbose)
{
    ex_opts(abort * EX_ABORT | debug * EX_DEBUG | verbose * EX_VERBOSE);
}

void Writer_Exodus::write_mesh(std::string aFilePath, const std::string& aFileName)
{
    Writer_Exodus::create_file(aFilePath, aFileName);
    Writer_Exodus::write_nodes();
    Writer_Exodus::write_node_sets();
    Writer_Exodus::write_blocks();
    Writer_Exodus::write_side_sets();
}

//----------------------------------------------------------------------------------------------------------------------

void Writer_Exodus::set_nodal_fields(moris::Cell<std::string> aFieldNames)
{
    if (aFieldNames.size() > 0)
    {
        // Write the number of nodal fields
        ex_put_variable_param(mExoid, EX_NODAL, aFieldNames.size());

        // Write the nodal field names and store as a map
        for (moris::uint tFieldIndex = 0; tFieldIndex < aFieldNames.size(); tFieldIndex++)
        {
            ex_put_variable_name(mExoid, EX_NODAL, tFieldIndex + 1, aFieldNames(tFieldIndex).c_str());
            mNodalFieldNamesMap[aFieldNames(tFieldIndex)] = tFieldIndex;
        }
    }
}

//----------------------------------------------------------------------------------------------------------------------

void Writer_Exodus::set_elemental_fields(moris::Cell<std::string> aFieldNames)
{
    if (aFieldNames.size() > 0)
    {
        // Write the number of elemental fields
        ex_put_variable_param(mExoid, EX_ELEM_BLOCK, aFieldNames.size());

        // Write the elemental field names and store as a map
        for (moris::uint tFieldIndex = 0; tFieldIndex < aFieldNames.size(); tFieldIndex++)
        {
            ex_put_variable_name(mExoid, EX_ELEM_BLOCK, tFieldIndex + 1, aFieldNames(tFieldIndex).c_str());
            mElementalFieldNamesMap[aFieldNames(tFieldIndex)] = tFieldIndex;
        }
    }
}

//----------------------------------------------------------------------------------------------------------------------

void Writer_Exodus::set_global_variables(moris::Cell<std::string> aVariableNames)
{
    if (aVariableNames.size() > 0)
    {
        // Write the number of global fields
        ex_put_variable_param(mExoid, EX_GLOBAL, aVariableNames.size());

        // Write the global field names and store as a map
        for (moris::uint tVariableIndex = 0; tVariableIndex < aVariableNames.size(); tVariableIndex++)
        {
            ex_put_variable_name(mExoid, EX_GLOBAL, tVariableIndex + 1, aVariableNames(tVariableIndex).c_str());
            mGlobalVariableNamesMap[aVariableNames(tVariableIndex)] = tVariableIndex;
        }
    }
}

//----------------------------------------------------------------------------------------------------------------------

void Writer_Exodus::set_time(moris::real aTimeValue)
{
    ex_put_time(mExoid, ++mTimeStep, &aTimeValue);
}

//----------------------------------------------------------------------------------------------------------------------

void Writer_Exodus::write_nodal_field(std::string aFieldName, moris::Matrix<moris::DDRMat> aFieldValues)
{
    // Field name to index
    moris::uint tMapSize = mNodalFieldNamesMap.size();
    int tFieldIndex = mNodalFieldNamesMap[aFieldName];
    MORIS_ERROR(mNodalFieldNamesMap.size() == tMapSize, aFieldName.append(
            " is not a nodal field name on this mesh!").c_str());

    // Check number of field values = number of nodes
    MORIS_ERROR(aFieldValues.numel() == mMesh->get_num_nodes(), aFieldName.append(
            " field was attempted to be written with " + std::to_string(aFieldValues.numel()) + " values, but there are " +
            std::to_string(mMesh->get_num_nodes()) + " nodes in this mesh.").c_str());

    // Write the field
    ex_put_var(mExoid, mTimeStep, EX_NODAL, tFieldIndex + 1, 0, aFieldValues.numel(), aFieldValues.data());
}

//----------------------------------------------------------------------------------------------------------------------

void Writer_Exodus::write_elemental_field(std::string aBlockName, std::string aFieldName,
        moris::Matrix<moris::DDRMat> aFieldValues)
{
    // Block name to index
    moris::uint tMapSize = mBlockNamesMap.size();
    int tBlockIndex = mBlockNamesMap[aBlockName];
    MORIS_ERROR(mBlockNamesMap.size() == tMapSize, aBlockName.append(
            " is not a block name on this mesh!").c_str());

    // Field name to index
    tMapSize = mElementalFieldNamesMap.size();
    int tFieldIndex = mElementalFieldNamesMap[aFieldName];
    MORIS_ERROR(mElementalFieldNamesMap.size() == tMapSize, aFieldName.append(
            " is not an elemental field name on this mesh!").c_str());

    // Check number of field values = number of elements
    MORIS_ERROR(aFieldValues.numel() == mMesh->get_block_set_cells(aBlockName).size(), aFieldName.append(
            " field was attempted to be written with " + std::to_string(aFieldValues.numel()) + " values, but there are " +
            std::to_string(mMesh->get_block_set_cells(aBlockName).size()) + " elements in block " + aBlockName +
            ".").c_str());

    // Write the field
    ex_put_var(mExoid, mTimeStep, EX_ELEM_BLOCK, tFieldIndex + 1, tBlockIndex + 1, aFieldValues.numel(),
            aFieldValues.data());
}

//----------------------------------------------------------------------------------------------------------------------

void Writer_Exodus::write_global_variable(std::string aVariableName, moris::real aVariableValue)
{
    // Variable name to index
    moris::uint tMapSize = mGlobalVariableNamesMap.size();
    int tVariableIndex = mGlobalVariableNamesMap[aVariableName];
    MORIS_ERROR(mGlobalVariableNamesMap.size() == tMapSize, aVariableName.append(
            " is not a global variable name on this mesh!").c_str());

    // Write the variable
    ex_put_var(mExoid, mTimeStep, EX_GLOBAL, tVariableIndex + 1, 0, 1, &aVariableValue);
}

//----------------------------------------------------------------------------------------------------------------------

void Writer_Exodus::create_file(std::string aFilePath, const std::string& aFileName)
{
    // Add temporary and permanent file names to file paths
    if (!aFilePath.empty())
    {
        aFilePath += "/";
    }
    mTempFileName = aFilePath + "temp.exo";
    mPermFileName = aFilePath + aFileName;

    // Make file name parallel, if necessary
    if (moris::par_size() > 1)
    {
        mTempFileName += "." + std::to_string(moris::par_size()) + "." + std::to_string(moris::par_rank());
        mPermFileName += "." + std::to_string(moris::par_size()) + "." + std::to_string(moris::par_rank());
    }

    // Create the database
    int tCPUWordSize = sizeof(moris::real), tIOWordSize = 0;
    mExoid = ex_create(mTempFileName.c_str(), EX_CLOBBER, &tCPUWordSize, &tIOWordSize);

    // Number of dimensions
    int tNumDimensions = mMesh->get_spatial_dim();

    // Number of nodes
    int tNumNodes = mMesh->get_num_nodes();

    // Number of node sets
    moris::Cell<std::string> tNodeSetNames = mMesh->get_set_names(moris::EntityRank::NODE);
    int tNumNodeSets = tNodeSetNames.size();

    // Number of side sets
    moris::Cell<std::string> tSideSetNames = mMesh->get_set_names(moris::EntityRank::FACE);
    int tNumSideSets = tSideSetNames.size();

    // Number of elements
    int tNumElements = 0;
    moris::Cell<std::string> tBlockNames = mMesh->get_set_names(moris::EntityRank::ELEMENT);
    int tNumElementBlocks = tBlockNames.size();
    for (moris::uint tBlockIndex = 0; tBlockIndex < tBlockNames.size(); tBlockIndex++)
    {
        moris::Cell<const moris::mtk::Cell *> tElementsInBlock = mMesh->get_block_set_cells(tBlockNames(tBlockIndex));
        tNumElements += tElementsInBlock.size();
    }

    // Initialize database
    ex_put_init(mExoid, "MTK", tNumDimensions, tNumNodes, tNumElements, tNumElementBlocks, tNumNodeSets, tNumSideSets);
}

//----------------------------------------------------------------------------------------------------------------------

void Writer_Exodus::open_file(std::string aExodusFileName, bool aReadOnly, float aVersion)
{
    int tCPUWordSize = sizeof(moris::real), tIOWordSize = 0;
    if (aReadOnly)
    {
        mExoid = ex_open(aExodusFileName.c_str(), EX_READ, &tCPUWordSize, &tIOWordSize, &aVersion);
    }
    else
    {
        mExoid = ex_open(aExodusFileName.c_str(), EX_WRITE, &tCPUWordSize, &tIOWordSize, &aVersion);
    }
}

//----------------------------------------------------------------------------------------------------------------------

void Writer_Exodus::close_file(bool aRename)
{
    ex_close(mExoid);
    if (aRename)
    {
        std::rename(mTempFileName.c_str(), mPermFileName.c_str());
    }
}

//----------------------------------------------------------------------------------------------------------------------

void Writer_Exodus::write_nodes()
{
    // Get all vertices
    moris::Cell<moris::mtk::Vertex const*> tNodes = mMesh->get_all_vertices();

    // spatial dimension
    int tSpatialDim = mMesh->get_spatial_dim();
    bool tYDim = tSpatialDim >= 2;
    bool tZDim = tSpatialDim >= 3;

    // Set up coordinate and node map arrays based on the number of vertices
    MORIS_ERROR(tNodes.size() > 0, "Invalid Node Map size");
    moris::Matrix<moris::IdMat> tNodeMap(tNodes.size(), 1, 0);

    // Coordinate arrays
    moris::Matrix<moris::DDRMat> tXCoordinates(tNodes.size(), 1, 0.0);
    moris::Matrix<moris::DDRMat> tYCoordinates(tNodes.size(), 1, 0.0);
    moris::Matrix<moris::DDRMat> tZCoordinates(tNodes.size(), 1, 0.0);

    for (moris::uint tNodeIndex = 0; tNodeIndex < tNodes.size(); tNodeIndex++)
    {
        // Get coordinates
        moris::Matrix<moris::DDRMat> tNodeCoordinates = tNodes(tNodeIndex)->get_coords();

        // Place in coordinate arrays
        tXCoordinates(tNodeIndex, 0) = tNodeCoordinates(0);
        tYCoordinates(tNodeIndex, 0) = tNodeCoordinates(1 * tYDim) * tYDim;
        tZCoordinates(tNodeIndex, 0) = tNodeCoordinates(2 * tZDim) * tZDim;

        // Get global ids for id map
        tNodeMap(tNodeIndex, 0) = tNodes(tNodeIndex)->get_id();
    }

    // Write coordinates
    ex_put_coord(mExoid, tXCoordinates.data(), tYCoordinates.data(), tZCoordinates.data());

    // Write node id map
    ex_put_id_map(mExoid, EX_NODE_MAP, tNodeMap.data());
}

//----------------------------------------------------------------------------------------------------------------------

void Writer_Exodus::write_node_sets()
{
    // Get the number of node sets and their names
    moris::Cell<std::string> tNodeSetNames = mMesh->get_set_names(moris::EntityRank::NODE);
    moris::uint tNumNodeSets = tNodeSetNames.size();

    // Write each node set
    for (moris::uint tNodeSetIndex = 0; tNodeSetIndex < tNumNodeSets; tNodeSetIndex++)
    {
        moris::Matrix<moris::IndexMat> tNodeIndices = mMesh->get_set_entity_loc_inds(moris::EntityRank::NODE,
                tNodeSetNames(tNodeSetIndex));
        ex_put_set_param(mExoid, EX_NODE_SET, tNodeSetIndex + 1, tNodeIndices.numel(), 0);
        ex_put_set(mExoid, EX_NODE_SET, tNodeSetIndex + 1, tNodeIndices.data(), nullptr);
    }
}

//----------------------------------------------------------------------------------------------------------------------

void Writer_Exodus::write_blocks()
{
    // Start the element maps
    mMtkExodusElementIndexMap.set_size(mMesh->get_num_elems(), 1, 0);
    moris::Matrix<moris::IdMat> tElementIdMap(0, 1, 0);
    moris::uint tElementIdMapStartIndex = 0;

    // All of the block names
    moris::Cell<std::string> tBlockNames = mMesh->get_set_names(moris::EntityRank::ELEMENT);

    // Loop through each block
    for (moris::uint tBlockIndex = 0; tBlockIndex < tBlockNames.size(); tBlockIndex++)
    {
        // Get the block elements
        moris::Cell<const moris::mtk::Cell*> tElementsInBlock = mMesh->get_block_set_cells(tBlockNames(tBlockIndex));

        // Add name to map
        mBlockNamesMap[tBlockNames(tBlockIndex)] = tBlockIndex;

        if (tElementsInBlock.size() > 0) // Block has at least 1 element
        {
            // Resize element map
            tElementIdMap.resize(tElementIdMapStartIndex + tElementsInBlock.size(), 1);

            // Get the CellTopology of this block
            //CellTopology tMorisBlockTopology = mMesh->get_blockset_topology(tBlockNames(tBlockIndex));

            // Get a description of the type of elements in this block FIXME once we always have a CellTopology on the mesh
            CellTopology tMorisBlockTopology;
            if (mMesh->get_spatial_dim() == 2)
            {
                if (tElementsInBlock(0)->get_vertex_inds().numel() == 3)
                {
                    tMorisBlockTopology = CellTopology::TRI3;
                }
                else
                {
                    tMorisBlockTopology = CellTopology::QUAD4;
                }
            }
            else
            {
                tMorisBlockTopology = mMesh->get_blockset_topology(tBlockNames(tBlockIndex));
            }

            const char* tExodusBlockTopology = this->get_exodus_block_topology(tMorisBlockTopology);

            // Get the number of nodes/edges/faces/attributes per element
            int tNumNodesPerElement = this->get_nodes_per_element(tMorisBlockTopology);
            int tNumEdgesPerElement = 0;
            int tNumFacesPerElement = 0;
            int tNumAttributesPerElement = 0;

            // Make a block and name it
	        ex_put_block(mExoid, EX_ELEM_BLOCK, tBlockIndex + 1, tExodusBlockTopology, tElementsInBlock.size(),
	                tNumNodesPerElement, tNumEdgesPerElement, tNumFacesPerElement, tNumAttributesPerElement);
            ex_put_name(mExoid, EX_ELEM_BLOCK, tBlockIndex + 1, tBlockNames(tBlockIndex).c_str());

            // Construct matrix of node indices per element
            moris::Matrix<moris::IndexMat> tConnectivityArray(tNumNodesPerElement * tElementsInBlock.size(), 1, 0);

            // Loop through the elements in this block
            moris::uint tConnectivityIndex = 0;
            for (moris::uint tElementIndex = 0; tElementIndex < tElementsInBlock.size(); tElementIndex++)
            {
                // Get the vertex indices of this element
            	moris::Matrix<moris::IndexMat>  tNodeIndices = tElementsInBlock(tElementIndex)->get_vertex_inds();

                // Assign each vertex individually
                for (int tNodeNum = 0; tNodeNum < tNumNodesPerElement; tNodeNum++)
                {
                    tConnectivityArray(tConnectivityIndex, 0) = tNodeIndices(tNodeNum) + 1;
                    tConnectivityIndex++;
                }

                // Get the global id and index of this element and add to maps
                tElementIdMap(tElementIdMapStartIndex + tElementIndex) = tElementsInBlock(tElementIndex)->get_id();
                mMtkExodusElementIndexMap(tElementsInBlock(tElementIndex)->get_index()) = tElementIdMapStartIndex + tElementIndex;
            }

            // Update location in element map
            tElementIdMapStartIndex += tElementsInBlock.size();

            // Write connectivity
            ex_put_conn(mExoid, EX_ELEM_BLOCK, tBlockIndex + 1, tConnectivityArray.data(), nullptr, nullptr);
        }

        else // Block has no elements
        {
            ex_put_block(mExoid, EX_ELEM_BLOCK, tBlockIndex + 1, "N/A", 0, 0, 0, 0, 0);
            ex_put_name(mExoid, EX_ELEM_BLOCK, tBlockIndex + 1, tBlockNames(tBlockIndex).c_str());
        }
    }

    // Write the element map
    ex_put_id_map(mExoid, EX_ELEM_MAP, tElementIdMap.data());

}

//----------------------------------------------------------------------------------------------------------------------

void Writer_Exodus::write_side_sets()
{
    // Get side set names
    moris::Cell<std::string> tSideSetNames = mMesh->get_set_names(moris::EntityRank::FACE);

    // Write side sets
    for (moris::uint tSideSetIndex = 0; tSideSetIndex < tSideSetNames.size(); tSideSetIndex++)
    {
        // Get the side set element ids
        moris::Matrix<moris::IndexMat>  tSideSetElements;
        moris::Matrix<moris::IndexMat>  tSideSetOrdinals;
        mMesh->get_sideset_elems_loc_inds_and_ords(tSideSetNames(tSideSetIndex), tSideSetElements, tSideSetOrdinals);

        // Change element and ordinal to what Exodus wants
        for (moris::uint tElementNum = 0; tElementNum < tSideSetOrdinals.numel(); tElementNum++)
        {
            tSideSetElements(tElementNum) = mMtkExodusElementIndexMap(tSideSetElements(tElementNum)) + 1;
            tSideSetOrdinals(tElementNum)++;
        }

        // Write the side set
        ex_put_set_param(mExoid, EX_SIDE_SET, tSideSetIndex + 1, tSideSetElements.numel(), 0);
        ex_put_set(mExoid, EX_SIDE_SET, tSideSetIndex + 1, tSideSetElements.data(), tSideSetOrdinals.data());
        ex_put_name(mExoid, EX_SIDE_SET, tSideSetIndex + 1, tSideSetNames(tSideSetIndex).c_str());
    }
}

// ---------------------------------------------------------------------------------------------------------------------
// STATIC (for now)
// ---------------------------------------------------------------------------------------------------------------------

const char* Writer_Exodus::get_exodus_block_topology(CellTopology aCellTopology)
{
    switch (aCellTopology)
    {
        case CellTopology::TRI3:
            return "TRI3";
        case CellTopology::QUAD4:
            return "QUAD4";
        case CellTopology::TET4:
            return "TET4";
        case CellTopology::TET10:
            return "TET10";
        case CellTopology::HEX8:
            return "HEX8";
        case CellTopology::PRISM6:
            return "PRISM6";
        default:
            MORIS_ERROR(0, "This element is invalid or it hasn't been implemented yet!");
            return "";
    }
}

//----------------------------------------------------------------------------------------------------------------------

int Writer_Exodus::get_nodes_per_element(CellTopology aCellTopology)
{
    switch (aCellTopology)
    {
        case CellTopology::TRI3:
            return 3;
        case CellTopology::QUAD4:
            return 4;
        case CellTopology::TET4:
            return 4;
        case CellTopology::TET10:
            return 10;
        case CellTopology::HEX8:
            return 8;
        case CellTopology::PRISM6:
            return 6;
        default:
            MORIS_ERROR(0, "This element is invalid or it hasn't been implemented yet!");
            return 0;
    }
}
