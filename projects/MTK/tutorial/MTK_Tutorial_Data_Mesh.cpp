/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * MTK_Tutorial_Data_Mesh.cpp
 *
 */

#include "cl_Communication_Manager.hpp" // COM/src
#include "cl_Communication_Tools.hpp"
#include "cl_Logger.hpp" // MRS/IOS/src

moris::Comm_Manager gMorisComm;
moris::Logger       gLogger;

/*
 * Include the necessary files from MTK
 */
#include "cl_MTK_Mesh.hpp" // MTK/src
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Mesh_Factory.hpp"
#include "cl_MTK_Scalar_Field_Info.hpp"

// Set namespaces to use
using namespace moris;
using namespace mtk;

/* @brief
 * This tutorial shows how to create an MTK mesh
 * from data using the mtk mesh data structure
 * in serial with multiple element types
 */
int
main( int    argc,
      char * argv[] )
{

    // Initialize the communication manager
    gMorisComm.initialize(&argc, &argv);

    // Severity level 0 - all outputs
    gLogger.initialize( 0 );

    // Assert this is a serial tutorial
    MORIS_ASSERT(par_size() == 1,"MTK_Tutorial_Data_Mesh.cpp is only a serial test");

    /*!
     * The goal of this tutorial is to generate the mesh shown in the following figures (fields not shown):
     * @image html ./figures/MTK_tutorial_1_cell.png "Mesh for MTK tutorial 1 with Cell Ids"
     * @image html ./figures/MTK_tutorial_1_nodes.png "Mesh for MTK tutorial 1 with Node Ids"
     *
     * Specify the mesh spatial dimension
     * \code{.cpp}
     *
     * // Spatial dimension of the mesh
     * uint aNumDim = 3;
     * \endcode
     *
     * Each cell to node connectivity and cell global to local map needs to
     * be grouped by cell type, this allow the mesh to create multiple types
     * of cell topologies. As shown here, we have the cell to node connectivity
     * of the hex8 and tet4 elements in their own matrix. Also, the element local
     * to global map has a matrix for hex8 and tet4s.
     * \code{.cpp}
     * // Number of element types (hex8,tet4)
     * uint aNumElemTypes = 2;
     *
     * // Specify the element to node connectivity of all hex 8's in the mesh
     * Matrix< IdMat >  aElemConnHex8  = {{ 1, 2, 4, 3, 5,  6,  8,  7},
     *                                    { 5, 6, 8, 7, 9, 10, 12, 11}};
     *
     * // Specify the element to node connectivity of all tet 4's in the mesh
     * Matrix< IdMat >  aElemConnTet4  = {{ 8, 12, 10, 13}};
     *
     * // Specify the local to global element map for all hex8s in the mesh
     * Matrix< IdMat >  aElemLocaltoGlobalHex8 = {{1,2}};
     *
     * // Specify the local to global element map for all tet4s in the mesh
     * Matrix< IdMat >  aElemLocaltoGlobalTet4 = {{3}};
     *
     * \endcode
     *
     *  The node coordinates matrix contains the nodal coordinate in the order that
     *  nodes appear in the node local to global map. For example, row 0 of the
     *  coordinate matrix is the node coordinate for node 1.
     * \code{.cpp}
     * // Node coordinate matrix
     * Matrix< DDRMat >  aCoords   = {{0.0, 0.0, 0.0},
     *                                {1.0, 0.0, 0.0},
     *                                {0.0, 1.0, 0.0},
     *                                {1.0, 1.0, 0.0},
     *                                {0.0, 0.0, 1.0},
     *                                {1.0, 0.0, 1.0},
     *                                {0.0, 1.0, 1.0},
     *                                {1.0, 1.0, 1.0},
     *                                {0.0, 0.0, 2.0},
     *                                {1.0, 0.0, 2.0},
     *                                {0.0, 1.0, 2.0},
     *                                {1.0, 1.0, 2.0},
     *                                {2.0, 1.0, 2.0}};
     *
     * // Specify the local to global map
     * Matrix< IdMat >  aNodeLocaltoGlobal     = {{1 ,2 ,3, 4, 5, 6 , 7, 8, 9, 10, 11, 12, 13}};
     *
     * \endcode
     */

    // Spatial dimension of the mesh
    uint aNumDim = 3;

    // Number of element types (hex8,tet4)
    uint aNumElemTypes = 2;

    // Specify the element to node connectivity of all hex 8's in the mesh
    Matrix< IdMat >  aElemConnHex8  = {{ 1, 2, 4, 3, 5,  6,  8,  7},
                                       { 5, 6, 8, 7, 9, 10, 12, 11}};

    // Specify the element to node connectivity of all tet 4's in the mesh
    Matrix< IdMat >  aElemConnTet4  = {{ 8, 12, 10, 13}};

    // Node coordinate matrix
    Matrix< DDRMat >  aCoords   = {{0.0, 0.0, 0.0},
                                   {1.0, 0.0, 0.0},
                                   {0.0, 1.0, 0.0},
                                   {1.0, 1.0, 0.0},
                                   {0.0, 0.0, 1.0},
                                   {1.0, 0.0, 1.0},
                                   {0.0, 1.0, 1.0},
                                   {1.0, 1.0, 1.0},
                                   {0.0, 0.0, 2.0},
                                   {1.0, 0.0, 2.0},
                                   {0.0, 1.0, 2.0},
                                   {1.0, 1.0, 2.0},
                                   {2.0, 1.0, 2.0}};

    // Specify the local to global map
    Matrix< IdMat >  aNodeLocaltoGlobal     = {{1 ,2 ,3, 4, 5, 6 ,7, 8, 9, 10, 11, 12, 13}};

    // Specify the local to global element map for all hex8s in the mesh
    Matrix< IdMat >  aElemLocaltoGlobalHex8 = {{1,2}};

    // Specify the local to global element map for all tet4s in the mesh
    Matrix< IdMat >  aElemLocaltoGlobalTet4 = {{3}};

    /*!
     * Setup the MtkSetsInfo data structure, this data structure
     * stores Cells of the various types of sets. Currently,
     * we support and their respective data structures:
     *  - 1.) Block sets - MtkBlockSetInfo
     *  - 2.) Side sets  - MtkSideSetInfo
     *  - 3.) Node/Vertex Sets - MtkNodeSetInfo
     *  \code{.cpp}
     *    // Initialize Sets information structure
     *    moris::mtk::MtkSetsInfo tMtkMeshSets;
     *  \endcode
     *
     *
     *  First, lets set up the block sets with element Id 1
     *  in the block set named "blockset_1" and element with
     *  id 2 in the block set named "blockset_2". Each block
     *  set requires a pointer to a vector of cells ids in
     *  the block set, a block set name, and a block set
     *  cell topology.
     *  \code{.cpp}
     *   // Cells in block set ids
     *   Matrix< IdMat >tCellIdsBS1({{1}});
     *
     *   moris::mtk::MtkBlockSetInfo tBlockSet1;
     *   tBlockSet1.mCellIdsInSet = &tCellIdsBS1;
     *   tBlockSet1.mBlockSetName = "blockset_1";
     *   tBlockSet1.mBlockSetTopo = CellTopology::HEX8;
     *
     *   Matrix< IdMat >tCellIdsBS2({{2}});
     *   moris::mtk::MtkBlockSetInfo tBlockSet2;
     *   tBlockSet2.mCellIdsInSet = &tCellIdsBS2;
     *   tBlockSet2.mBlockSetName = "blockset_2";
     *   tBlockSet2.mBlockSetTopo = CellTopology::HEX8;
     *
     *  \endcode
     *
     *
     *  The block sets are then added to the Mesh sets info data structure
     *  \code{.cpp}
     *  // Add block sets to mtk mesh sets
     *  tMtkMeshSets.add_block_set(&tBlockSet1);
     *  tMtkMeshSets.add_block_set(&tBlockSet2);
     *  \endcode
     *
     *
     */

     // Initialize Sets information structure
     moris::mtk::MtkSetsInfo tMtkMeshSets;

     // Declare block sets
     //////////////////////

     // Cells in block set ids
     Matrix< IdMat >tCellIdsBS1({{1}});

     moris::mtk::MtkBlockSetInfo tBlockSet1;
     tBlockSet1.mCellIdsInSet = &tCellIdsBS1;
     tBlockSet1.mBlockSetName = "blockset_1";
     tBlockSet1.mBlockSetTopo = CellTopology::HEX8;

     Matrix< IdMat >tCellIdsBS2({{2}});
     moris::mtk::MtkBlockSetInfo tBlockSet2;
     tBlockSet2.mCellIdsInSet = &tCellIdsBS2;
     tBlockSet2.mBlockSetName = "blockset_2";
     tBlockSet2.mBlockSetTopo = CellTopology::HEX8;

     Matrix< IdMat >tCellIdsBS3({{3}});
     moris::mtk::MtkBlockSetInfo tBlockSet3;
     tBlockSet3.mCellIdsInSet = &tCellIdsBS3;
     tBlockSet3.mBlockSetName = "blockset_3";
     tBlockSet3.mBlockSetTopo = CellTopology::TET4;

     // Add block sets to mtk mesh sets
     tMtkMeshSets.add_block_set(&tBlockSet1);
     tMtkMeshSets.add_block_set(&tBlockSet2);
     tMtkMeshSets.add_block_set(&tBlockSet3);

     /*!
      * Next we set up the side sets, Each side set,
      * requires the element id and the side ordinal of
      * the respective element on the side set rather than
      * declaring side sets through face ids explicitly.
      * A side set also needs a name. For this example,
      * the side sets are seen in the following figure.
      * @image html ./figures/MTK_tutorial_1_side_set.png "Mesh for MTK tutorial 1 with side sets"
      *
      * \code{.cpp}
      *   Matrix< IdMat > tElemIdsAndSideOrdsSS1  = { {1, 3},
      *                                               {1, 4},
      *                                               {1, 5},
      *                                               {2, 1},
      *                                               {2, 2} };
      *
      *   moris::mtk::MtkSideSetInfo tSideSetStruc;
      *   tSideSetStruc.mElemIdsAndSideOrds = &tElemIdsAndSideOrdsSS1;
      *   tSideSetStruc.mSideSetName        = "Sideset_1" ;
      *
      *   // Add side side set to mesh sets
      *   tMtkMeshSets.add_side_set(&tSideSetStruc);
      *
      * \endcode
      */

     // Declare side sets
     /////////////////////
     Matrix< IdMat > tElemIdsAndSideOrdsSS1  = { {1, 3},
                                                 {1, 4},
                                                 {1, 5},
                                                 {2, 1},
                                                 {2, 2} };

     moris::mtk::MtkSideSetInfo tSideSetStruc;
     tSideSetStruc.mElemIdsAndSideOrds = &tElemIdsAndSideOrdsSS1;
     tSideSetStruc.mSideSetName        = "Sideset_1" ;

     // Add side side set to mesh sets
     tMtkMeshSets.add_side_set(&tSideSetStruc);

     /*!
      * Next we set up the node sets, each node set
      * has the node ids belonging to the set and a name.
      * The node set for this example is shown below.
      *
      * @image html ./figures/MTK_tutorial_1_node_set.png "Mesh for MTK tutorial 1 with node sets"
      * \code{.cpp}
      * Matrix< IdMat > tNodeIdsNS1  = { {1}, {3}, {5}, {6} };
      *
      * moris::mtk::MtkNodeSetInfo tNodeSet1;
      * tNodeSet1.mNodeIds        = &tNodeIdsNS1;
      * tNodeSet1.mNodeSetName    = "Nodeset_1" ;
      *
      * // Add node set to Mtk mesh sets
      * tMtkMeshSets.add_node_set(&tNodeSet1);
      * \endcode
      */
     Matrix< IdMat > tNodeIdsNS1  = { {1}, {3}, {5}, {6} };

     moris::mtk::MtkNodeSetInfo tNodeSet1;
     tNodeSet1.mNodeIds        = &tNodeIdsNS1;
     tNodeSet1.mNodeSetName    = "Nodeset_1" ;

     // Add node set to Mtk mesh sets
     tMtkMeshSets.add_node_set(&tNodeSet1);

     /*
      * Now that we have all the mesh sets setup, adding fields
      * to entities is explained. First, a real scalar node field
      * setup is demonstrated where the field is defined on the
      * entire mesh. The node field name is "node_field_1".
      * If the field part name is not specified, MTK assumes that
      * the field is declared over the entire mesh.
      *
      * @image html ./figures/MTK_tutorial_1_node_field.png "Mesh for MTK tutorial 1 colored by node_field_1"
      *
      * \code{.cpp}
      *     // Initialize Scalar_Field_Info structure with a DDRMat as the template type
      *     moris::mtk::Scalar_Field_Info<DDRMat> tNodeField;
      *
      *     // number of nodes in the mesh
      *     uint tNumNodes = aNodeLocaltoGlobal.numel();
      *
      *     // Allocate field information
      *     Matrix<DDRMat> tNodeFieldData = {{0.1, 0.2, 0.3, 0.4, 0.5, 0.6 ,0.5, 0.4, 0.3, 0.2, 0.1, 0.0, -0.1}};
      *
      *     // Set the field name
      *     tNodeField.set_field_name("node_field_1");
      *
      *     // Set the entity rank associated with this field
      *     tNodeField.set_field_entity_rank(EntityRank::NODE);
      *
      *     // Give a reference to the data and entity ids associated with the data.
      *     // In this case, the local to global node map can be used since this
      *     // is a universal node field.
      *     tNodeField.add_field_data(&aNodeLocaltoGlobal, &tNodeFieldData);
      * \endcode
      *
      */

     // Initialize Scalar_Field_Info structure with a DDRMat as the template type
     moris::mtk::Scalar_Field_Info<DDRMat> tNodeField;

     // Allocate field information (1 entry for each node)
     Matrix<DDRMat> tNodeFieldData = {{0.1, 0.2, 0.3, 0.4, 0.5, 0.6 ,0.5, 0.4, 0.3, 0.2, 0.1, 0.0, -0.1}};

     // Set the field name
     tNodeField.set_field_name("node_field_1");

     // Set the entity rank associated with this field
     tNodeField.set_field_entity_rank(EntityRank::NODE);

     // Give a reference to the data and entity ids associated with the data.
     // In this case, the local to global node map can be used since this
     // is a universal node field.
     tNodeField.add_field_data(&aNodeLocaltoGlobal, &tNodeFieldData);

     /*!
      * Once all the mesh fields have been setup, they are collected
      * into the MtkFieldsInfo data structure. The data structure currently
      * is limited to scalar fields of real and sint types. For this examples,
      * the node and element fields are real scalar fields and are added to the
      * MtkFieldsInfo mRealScalarFields.
      */
     // Initialize field information container
     moris::mtk::MtkFieldsInfo tFieldsInfo;

     // Place the node field into the field info container
     tFieldsInfo.mRealScalarFields.push_back(&tNodeField);

     /*!
      *  Now that all of the data, sets and fields have been setup,
      *  references to the information is provided to a MTKMeshData
      *  struct. This structure is passed to MTK via a factory which
      *  then loads the information into the appropriate structure.
      *  \code{.cpp}
      *       moris::mtk::MtkMeshData aMeshData(aNumElemTypes);
      *       aMeshData.SpatialDim              = &aNumDim;
      *       aMeshData.ElemConn(0)             = &aElemConnHex8;
      *       aMeshData.ElemConn(1)             = &aElemConnTet4;
      *       aMeshData.NodeCoords              = &aCoords;
      *       aMeshData.SetsInfo                = &tMtkMeshSets;
      *       aMeshData.LocaltoGlobalElemMap(0) = &aElemLocaltoGlobalHex8;
      *       aMeshData.LocaltoGlobalElemMap(1) = &aElemLocaltoGlobalTet4;
      *       aMeshData.LocaltoGlobalNodeMap    = &aNodeLocaltoGlobal;
      *       aMeshData.FieldsInfo              = &tFieldsInfo;
      *
      *       // Create mesh from data with the factory
      *       moris::mtk::Mesh* tMesh = create_interpolation_mesh( MeshType::STK, aMeshData );
      *  \endcode
      *
      */

     moris::mtk::MtkMeshData aMeshData(aNumElemTypes);
     aMeshData.SpatialDim              = &aNumDim;
     aMeshData.ElemConn(0)             = &aElemConnHex8;
     aMeshData.ElemConn(1)             = &aElemConnTet4;
     aMeshData.NodeCoords              = &aCoords;
     aMeshData.SetsInfo                = &tMtkMeshSets;
     aMeshData.LocaltoGlobalElemMap(0) = &aElemLocaltoGlobalHex8;
     aMeshData.LocaltoGlobalElemMap(1) = &aElemLocaltoGlobalTet4;
     aMeshData.LocaltoGlobalNodeMap    = &aNodeLocaltoGlobal;
     aMeshData.FieldsInfo              = &tFieldsInfo;

     // Create mesh from data with the factory
     moris::mtk::Mesh* tMesh = create_interpolation_mesh( MeshType::STK, aMeshData );

     /*!
      * To access the entities in a given set, the name of the set and
      * the rank of the entities associated with the set. Based on
      * the rank of the entity provided, the mesh can deduce the set type.
      *
      * \code{.cpp}
      *   Matrix< IndexMat >  tBlockIndices1   = tMesh->get_set_entity_loc_inds( EntityRank::ELEMENT, "blockset_1" );
      *   Matrix< IndexMat >  tBlockIndices2   = tMesh->get_set_entity_loc_inds( EntityRank::ELEMENT, "blockset_2" );
      *   Matrix< IndexMat >  tNodeSetIndices1 = tMesh->get_set_entity_loc_inds( EntityRank::NODE,    "Nodeset_1" );
      *   Matrix< IndexMat >  tSideSetIndices1 = tMesh->get_set_entity_loc_inds( EntityRank::FACE,    "Sideset_1" );
      *  \endcode
      *
      */

     Matrix< IndexMat >  tBlockIndices1   = tMesh->get_set_entity_loc_inds( EntityRank::ELEMENT, "blockset_1" );
     Matrix< IndexMat >  tBlockIndices2   = tMesh->get_set_entity_loc_inds( EntityRank::ELEMENT, "blockset_2" );
     Matrix< IndexMat >  tNodeSetIndices1 = tMesh->get_set_entity_loc_inds( EntityRank::NODE,    "Nodeset_1" );
     Matrix< IndexMat >  tSideSetIndices1 = tMesh->get_set_entity_loc_inds( EntityRank::FACE,    "Sideset_1" );

     /*!
      * The mesh can be exported to an exodus file via the create_output_mesh function call.
      * This tells the underlying mesh to write an exodus file internally
      */

     std::string tMeshOutputFile = "./MTK_Tutorial_1.e";
     tMesh->create_output_mesh(tMeshOutputFile);
     std::cout<<"Mesh outputted to file: " << tMeshOutputFile<<std::endl;

    gMorisComm.finalize();
}

