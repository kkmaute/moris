//
// Created by christopherson on 9/13/19.
//

#include <exodusII.h>
#include "cl_MTK_Writer_Exodus.hpp"

#include "cl_MTK_Mesh_Core.hpp"

Writer_Exodus::Writer_Exodus(moris::mtk::Mesh* aMeshPointer, std::string aFilePath, const std::string& aFileName)
    : mMesh(aMeshPointer)
{
    mExoid = 0;
    if (!aFilePath.empty())
    {
        aFilePath += "/";
    }
    mTempFileName = aFilePath + "temp.exo";
    mPermFileName = aFilePath + aFileName;
    mStartTime = 0;
}

Writer_Exodus::~Writer_Exodus()
= default;

void Writer_Exodus::set_error_options(bool abort, bool debug, bool verbose)
{
    ex_opts(abort * EX_ABORT|debug * EX_DEBUG|verbose * EX_VERBOSE);
}

void Writer_Exodus::write_mesh()
{
    Writer_Exodus::create_file();
    Writer_Exodus::write_coordinates();
    Writer_Exodus::close_file();
}

void Writer_Exodus::create_file() {
    int tCPUWordSize = sizeof(moris::real);
    int tIOWordSize = 0;

    uint tNumDimensions = mMesh->get_spatial_dim();
    uint tNumNodes = mMesh->get_num_nodes();
    uint tNumElements = mMesh->get_num_elems();
    uint tNumElementBlocks, tNumNodeSets, tNumSideSets;

    moris::Cell<std::string> tElementBlockNames;
    moris::Cell<std::string> tNodeSetNames;
    moris::Cell<std::string> tSideSetNames;

    // Element Blocks
    tElementBlockNames = mMesh->get_set_names(moris::EntityRank::ELEMENT);
    tNumElementBlocks = tElementBlockNames.size();

    // Node Sets
    tNodeSetNames = mMesh->get_set_names(moris::EntityRank::NODE);
    tNumNodeSets = tNodeSetNames.size();

    // Side Sets
    tSideSetNames = mMesh->get_set_names(moris::EntityRank::FACE);
    tNumSideSets = tSideSetNames.size();

    // Create the database
    mExoid = ex_create(mTempFileName.c_str(), EX_CLOBBER, &tCPUWordSize, &tIOWordSize);

    // Initialize the database
    ex_put_init(mExoid, "MTK", tNumDimensions, tNumNodes, tNumElements, tNumElementBlocks, tNumNodeSets,
                         tNumSideSets);
}

void Writer_Exodus::open_file(const char* aExodusFileName, moris::real aVersion)
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

//void Writer_Exodus::write_coordinates()
//{
//    unsigned long i;
//    moris::Cell<moris::mtk::Vertex const*> tvertices;
//    moris::Matrix<moris::DDRMat> tcoordinates;
//    moris::real x, y, z;
//
//    tvertices = mMesh->get_all_vertices();
//    for (i = 0; i < tvertices.size(); i++)
//    {
//        coordinates = vertices[i].get_coords();
//        x = coordinates(0);
//        y = coordinates(1);(0,0)
//    }
//    ex_put_coord(mExoid, x, y, z);
//}

