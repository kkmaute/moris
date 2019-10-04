//
// Created by christopherson on 9/13/19.
//

#include <exodusII.h>
#include "cl_MTK_Writer_Exodus.hpp"
#include "cl_MTK_Mesh_Core.hpp"

Writer_Exodus::Writer_Exodus(moris::mtk::Mesh* aMeshPointer, std::string aFilePath, const std::string& aFileName) //TODO remove path
    : mMesh(aMeshPointer)
{
    mExoid = 0;
    if (!aFilePath.empty())
    {
        aFilePath += "/";
    }
    mTempFileName = aFilePath + "temp.exo";
    mPermFileName = aFilePath + aFileName;
}

Writer_Exodus::~Writer_Exodus()
= default;

void Writer_Exodus::set_error_options(bool abort, bool debug, bool verbose)
{
    ex_opts(abort * EX_ABORT | debug * EX_DEBUG | verbose * EX_VERBOSE);
}

void Writer_Exodus::write_mesh()
{
    Writer_Exodus::create_file();
    Writer_Exodus::write_nodes();
    Writer_Exodus::write_node_sets();
    Writer_Exodus::write_blocks();
    Writer_Exodus::close_file();
}

void Writer_Exodus::create_file() {
    int                         tCPUWordSize = sizeof(moris::real),
                                tIOWordSize = 0;
    moris::uint                 tNumDimensions = mMesh->get_spatial_dim(),
                                tNumNodes = mMesh->get_num_nodes(),
                                tNumElements = mMesh->get_num_elems(),
                                tNumElementBlocks,
                                tNumNodeSets,
                                tNumSideSets;
    moris::Cell<std::string>    tElementBlockNames,
                                tNodeSetNames,
                                tSideSetNames;

    // Element Blocks
    tElementBlockNames = mMesh->get_set_names(moris::EntityRank::ELEMENT);
    tNumElementBlocks = tElementBlockNames.size();

    // Node Sets
    tNodeSetNames = mMesh->get_set_names(moris::EntityRank::NODE);
    tNumNodeSets = tNodeSetNames.size();

    // Side Sets
    tSideSetNames = mMesh->get_set_names(moris::EntityRank::FACE);
    tNumSideSets = tSideSetNames.size();
    tNumSideSets = 0; //TODO

    // Create the database
    mExoid = ex_create(mTempFileName.c_str(), EX_CLOBBER, &tCPUWordSize, &tIOWordSize);

    // Initialize the database
    ex_put_init(mExoid, "MTK", tNumDimensions, tNumNodes, tNumElements, tNumElementBlocks, tNumNodeSets,
                         tNumSideSets);
}

void Writer_Exodus::open_file(const char* aExodusFileName, float aVersion)
{
    int tCPUWordSize = sizeof(moris::real);
    int tIOWordSize = 0;
    ex_open(aExodusFileName, EX_READ, &tCPUWordSize, &tIOWordSize, &aVersion);
}

void Writer_Exodus::close_file()
{
    ex_close(mExoid);
    std::rename(mTempFileName.c_str(), mPermFileName.c_str());
}

void Writer_Exodus::write_nodes()
{
    unsigned long                           tVertexNum,
                                            tCoordinateNum;
    moris::Cell<moris::mtk::Vertex const*>  tVertices;
    moris::Matrix<moris::DDRMat>            tVertexCoordinates;
    moris::real*                            tCoordinateArray[3];

    tVertexNum = mMesh->get_num_nodes();
    tCoordinateArray[0] = new moris::real[tVertexNum]();
    tCoordinateArray[1] = new moris::real[tVertexNum]();
    tCoordinateArray[2] = new moris::real[tVertexNum]();

    tVertices = mMesh->get_all_vertices();
    for (tVertexNum = 0; tVertexNum < tVertices.size(); tVertexNum++)
    {
        tVertexCoordinates = tVertices(tVertexNum)->get_coords();
        for (tCoordinateNum = 0; tCoordinateNum < std::max(tVertexCoordinates.length(), (size_t)3); tCoordinateNum++)
        {
            tCoordinateArray[tCoordinateNum][tVertexNum] = tVertexCoordinates(tCoordinateNum);
        }
    }
    ex_put_coord(mExoid, tCoordinateArray[0], tCoordinateArray[1], tCoordinateArray[2]);

    delete tCoordinateArray[0];
    delete tCoordinateArray[1];
    delete tCoordinateArray[2];
}

void Writer_Exodus::write_node_sets()
{
    int                             tNodeSetId,
                                    tNumNodeSets;
    int*                            tNodeList;
    moris::Cell<std::string>        tNodeSetNames;
    moris::Matrix<moris::IndexMat>  tNodeIndices;

    // Get the number of node sets and their names
    tNodeSetNames = mMesh->get_set_names(moris::EntityRank::NODE);
    tNumNodeSets = tNodeSetNames.size();

    for (tNodeSetId = 0; tNodeSetId < tNumNodeSets; tNodeSetId++)
    {
        tNodeIndices = mMesh->get_set_entity_loc_inds(moris::EntityRank::NODE, tNodeSetNames(tNodeSetId));
        tNodeList = new int[tNodeIndices.length()];
        ex_put_set(mExoid, EX_NODE_SET, tNodeSetId, tNodeIndices.data(), NULL);
        delete tNodeList;
    }
}

void Writer_Exodus::write_blocks()
{
    unsigned long                   tBlockNum,
                                    tElementNum,
                                    tConnectivityIndex,
                                    tVertexNum;
    int64_t                         tNumNodesPerEntry,
                                    tNumEdgesPerEntry,
                                    tNumFacesPerEntry,
                                    tNumAttributesPerEntry = 0; // TODO attributes?
    int*                            tConnectivityArray;
    const char*                     tBlockDescription;
    CellTopology                    tBlockTopology;
    moris::Cell<std::string>        tSetNames;
    moris::Matrix<moris::IdMat>     tVertexIds;

    // All of the block names
    tSetNames = mMesh->get_set_names(moris::EntityRank::ELEMENT);

    // Loop through each block
    for (tBlockNum = 0; tBlockNum < tSetNames.size(); tBlockNum++)
    {
        // Get the block elements
        moris::Cell<const moris::mtk::Cell*> tElements = mMesh->get_block_set_cells(tSetNames(tBlockNum));

        if (tElements.size() > 0)
        {
            // Get the CellTopology of this block
            tBlockTopology = mMesh->get_blockset_topology(tSetNames(tBlockNum));

            // Get a description of the type of elements in this block
            tBlockDescription = this->get_exodus_block_description(tBlockTopology);

            // Get the number of nodes/edges/faces/attributes per element
            tNumNodesPerEntry = this->get_nodes_per_element(tBlockTopology);
            tNumEdgesPerEntry = this->get_edges_per_element(tBlockTopology);
            tNumFacesPerEntry = this->get_faces_per_element(tBlockTopology);

            // Make a block
            ex_put_block(mExoid, EX_ELEM_BLOCK, tBlockNum, tBlockDescription, tElements.size(), tNumNodesPerEntry,
                         tNumEdgesPerEntry, tNumFacesPerEntry, tNumAttributesPerEntry);

            // Construct and write connectivity
            tConnectivityArray = new int[tElements.size() * tNumNodesPerEntry]();
            tConnectivityIndex = 0;
            for (tElementNum = 0; tElementNum < tElements.size(); tElementNum++)
            {
                tVertexIds = tElements(tElementNum)->get_vertex_ids();
                for (tVertexNum = 0; tVertexNum < tVertexIds.length(); tVertexNum++)
                {
                    tConnectivityArray[tConnectivityIndex] = tVertexIds(tVertexNum);
                    tConnectivityIndex++;
                }
            }
            ex_put_conn(mExoid, EX_ELEM_BLOCK, tBlockNum, tConnectivityArray, nullptr, nullptr);
            delete tConnectivityArray;
        }
        else // Block has no elements
        {
            ex_put_block(mExoid, EX_ELEM_BLOCK, tBlockNum, "N/A", tElements.size(), 0, 0, 0, tNumAttributesPerEntry);
        }

        // Name the block sets
        ex_put_names(mExoid, EX_ELEM_BLOCK, this->string_cell_to_char_array(tSetNames));
    }
}











// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------

const char* Writer_Exodus::get_exodus_element_type(CellTopology aCellTopology)
{
    switch (aCellTopology)
    {
        case CellTopology::TRI3:
            return "TRIANGLE";
        case CellTopology::QUAD4:
            return "QUAD";
        case CellTopology::TET4:
            return "TETRA";
        case CellTopology::TET10:
            return "TETRA";
        case CellTopology::HEX8:
            return "HEX";
        case CellTopology::PRISM6:
            return "WEDGE";
        default:
            MORIS_ERROR(0, "This element is invalid or it hasn't been implemented yet!");
            return "";
    }
}

const char* Writer_Exodus::get_exodus_block_description(CellTopology aCellTopology)
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

int64_t Writer_Exodus::get_nodes_per_element(CellTopology aCellTopology)
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

int64_t Writer_Exodus::get_edges_per_element(CellTopology aCellTopology)
{
    switch (aCellTopology)
    {
        case CellTopology::TRI3:
            return 3;
        case CellTopology::QUAD4:
            return 4;
        case CellTopology::TET4:
            return 0;
        case CellTopology::TET10:
            return 0;
        case CellTopology::HEX8:
            return 0;
        case CellTopology::PRISM6:
            return 0;
        default:
            MORIS_ERROR(0, "This element is invalid or it hasn't been implemented yet!");
            return 0;
    }
}

int64_t Writer_Exodus::get_faces_per_element(CellTopology aCellTopology)
{
    switch (aCellTopology)
    {
        case CellTopology::TRI3:
            return 0;
        case CellTopology::QUAD4:
            return 0;
        case CellTopology::TET4:
            return 4;
        case CellTopology::TET10:
            return 4;
        case CellTopology::HEX8:
            return 6;
        case CellTopology::PRISM6:
            return 5;
        default:
            MORIS_ERROR(0, "This element is invalid or it hasn't been implemented yet!");
            return 0;
    }
}

char** Writer_Exodus::string_cell_to_char_array(moris::Cell<std::string> aStringCell)
{
    unsigned long i;
    char** tCharacterArray = new char*[aStringCell.size()];

    for (i = 0; i < aStringCell.size(); i++)
    {
        tCharacterArray[i] = new char[aStringCell(i).size()];
        tCharacterArray[i] = const_cast<char *>(aStringCell(i).data());
    }

    return tCharacterArray;
}
