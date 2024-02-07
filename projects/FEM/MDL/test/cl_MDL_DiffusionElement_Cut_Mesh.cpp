/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MDL_DiffusionElement_Cut_Mesh.cpp
 *
 */

#include "catch.hpp"

#include "cl_MTK_Vertex.hpp"    //MTK
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh.hpp"

#include "cl_MTK_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Tools.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Scalar_Field_Info.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Mesh_Data_STK.hpp"
#include "cl_MTK_Mesh_Core_STK.hpp"
#include "cl_MTK_Interpolation_Mesh_STK.hpp"
#include "cl_MTK_Integration_Mesh_STK.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Double_Side_Cluster.hpp"
#include "cl_MTK_Double_Side_Cluster_Input.hpp"
#include "cl_MTK_Side_Cluster.hpp"
#include "cl_MTK_Side_Cluster_Input.hpp"

#include "cl_Matrix.hpp"        //LINALG
#include "linalg_typedefs.hpp"
#include "fn_equal_to.hpp" // ALG/src

#include "cl_FEM_NodeProxy.hpp"                //FEM/INT/src
#include "cl_FEM_ElementProxy.hpp"             //FEM/INT/src
#include "cl_FEM_Node_Base.hpp"                //FEM/INT/src
#include "cl_FEM_Element_Factory.hpp"          //FEM/INT/src
#include "cl_FEM_IWG_Factory.hpp"              //FEM/INT/src
#include "cl_FEM_CM_Factory.hpp"              //FEM/INT/src
#include "cl_FEM_SP_Factory.hpp"              //FEM/INT/src
#include "cl_FEM_Set_User_Info.hpp"              //FEM/INT/src

#include "cl_MDL_Model.hpp"

#include "cl_HMR.hpp"
#include "cl_HMR_Background_Mesh.hpp" //HMR/src
#include "cl_HMR_BSpline_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Element.hpp" //HMR/src
#include "cl_HMR_Factory.hpp" //HMR/src
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Parameters.hpp" //HMR/src

#include "cl_DLA_Solver_Factory.hpp"
#include "cl_DLA_Solver_Interface.hpp"

#include "cl_NLA_Nonlinear_Solver_Factory.hpp"
#include "cl_NLA_Nonlinear_Solver.hpp"
#include "cl_NLA_Nonlinear_Problem.hpp"
#include "cl_NLA_Nonlinear_Algorithm.hpp"
#include "cl_MSI_Solver_Interface.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_MSI_Model_Solver_Interface.hpp"
#include "cl_DLA_Linear_Solver_Aztec.hpp"
#include "cl_DLA_Linear_Solver.hpp"
#include "cl_SOL_Warehouse.hpp"

#include "cl_TSA_Time_Solver_Factory.hpp"
#include "cl_TSA_Monolithic_Time_Solver.hpp"
#include "cl_TSA_Time_Solver.hpp"

#include "fn_norm.hpp"

#include "cl_Tracer.hpp"

namespace moris
{
namespace mdl
{

void tConstValFunction_MDLCUT
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  Vector< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 );
}

        TEST_CASE( "Diffusion_Cut", "[moris],[mdl],[Diffusion_Cut]" )
        {

        uint p_size = moris::par_size();

        if( p_size == 1 ) // specify it is a serial test only
        {
            // setup the interpolation mesh
            std::string tInterpString = "generated:1x1x4";
            moris::mtk::Interpolation_Mesh* tInterpMesh1 = moris::mtk::create_interpolation_mesh( mtk::MeshType::STK, tInterpString );

            // setup the integration mesh
            // Define the Integration Mesh (from data from xtk)
            Matrix<DDRMat>   tNodeCoordinates ={{0, 0, 0},{1, 0, 0},{0, 1, 0},{1, 1, 0},{0, 0, 1},{1, 0, 1},{0, 1, 1},{1, 1, 1},{0, 0, 2},{1, 0, 2},{0, 1, 2},{1, 1, 2},{0, 0, 3},{1, 0, 3},{0, 1, 3},{1, 1, 3},{0, 0, 4},{1, 0, 4},{0, 1, 4},{1, 1, 4},{0.5, 0, 3.5},{1, 0.5, 3.5},{0.5, 1, 3.5},{0, 0.5, 3.5},{0.5, 0.5, 3},{0.5, 0.5, 4},{0.5, 0.5, 3.5},{0.1, 0, 3.1},{0.9, 0, 3.1},{0.1, 0.1, 3.1},{0.9, 0.1, 3.1},{1, 0, 3.1},{0, 0, 3.1},{1, 0.1, 3.1},{1, 0.9, 3.1},{0.9, 0.9, 3.1},{1, 1, 3.1},{0.9, 1, 3.1},{0.1, 1, 3.1},{0.1, 0.9, 3.1},{0, 1, 3.1},{0, 0.9, 3.1},{0, 0.1, 3.1},{0.5, 0.5, 3.1}};
            Matrix<IndexMat> tLocalToGlobalNodeMap = {{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44}};

            Matrix<IndexMat> tInterpElemsAsIntegCellIds     = {{1,2,3,4}};
            Matrix<IndexMat> tInterpElemsAsIntegCellToNodes = {{1, 2, 4, 3, 5, 6, 8, 7},{5, 6, 8, 7, 9, 10, 12, 11},{9, 10, 12, 11, 13, 14, 16, 15},{13, 14, 16, 15, 17, 18, 20, 19}};

            // Tetrathedral cells in material phase 1
            Matrix<IndexMat> tCellIdsPhase0    = {{6, 8, 10, 12, 14, 16, 17, 18, 20, 31, 32, 33, 42, 43, 44, 53, 54, 55, 62, 63, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84}};
            Matrix<IndexMat> tCellToNodePhase0 = {{29, 14, 31, 32},{13, 28, 30, 33},{35, 16, 36, 37},{14, 31, 32, 34},{16, 36, 37, 38},{39, 15, 40, 41},{13, 30, 40, 43},{15, 40, 41, 42},{30, 13, 33, 43},{13, 14, 31, 29},{28, 13, 31, 29},{13, 28, 31, 30},{16, 14, 35, 36},{14, 31, 34, 36},{14, 34, 35, 36},{15, 16, 39, 40},{16, 36, 38, 40},{16, 38, 39, 40},{15, 13, 40, 43},{15, 40, 42, 43},{13, 30, 31, 44},{14, 13, 31, 44},{25, 13, 14, 44},{14, 31, 36, 44},{16, 14, 36, 44},{25, 14, 16, 44},{16, 36, 40, 44},{15, 16, 40, 44},{25, 16, 15, 44},{30, 13, 40, 44},{13, 15, 40, 44},{13, 25, 15, 44}};

            // Tetrathedral cells in material phase 1
            Matrix<IndexMat> tCellToNodeGhost0 = {{21, 27, 31, 30},{17, 18, 21, 27},{31, 27, 34, 36},{18, 20, 22, 27},{36, 27, 38, 40},{20, 19, 23, 27},{17, 24, 19, 27},{30, 27, 31, 44},{31, 27, 36, 44},{36, 27, 40, 44},{27, 30, 40, 44},{17, 26, 18, 27},{18, 26, 20, 27},{20, 26, 19, 27},{17, 19, 26, 27},{21, 28, 30, 31},{21, 28, 31, 29},{21, 29, 31, 32},{27, 21, 31, 32},{18, 21, 27, 32},{28, 21, 30, 33},{21, 27, 30, 33},{21, 17, 27, 33},{27, 22, 34, 36},{34, 22, 35, 36},{22, 35, 36, 37},{27, 22, 36, 37},{20, 22, 27, 37},{31, 27, 32, 34},{27, 18, 32, 34},{27, 22, 18, 34},{27, 23, 38, 40},{38, 23, 39, 40},{36, 27, 37, 38},{27, 20, 37, 38},{27, 23, 20, 38},{23, 39, 40, 41},{27, 23, 40, 41},{19, 23, 27, 41},{27, 24, 42, 43},{30, 27, 40, 43},{40, 27, 42, 43},{40, 27, 41, 42},{27, 19, 41, 42},{27, 24, 19, 42},{27, 30, 33, 43},{17, 27, 33, 43},{24, 27, 17, 43}};
            Matrix<IndexMat> tCellIdsGhost0 = {{5, 7, 9, 11, 13, 15, 19, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 34, 35, 36, 37, 38, 39, 40, 41, 45, 46, 47, 48, 49, 50, 51, 52, 56, 57, 58, 59, 60, 61, 64, 65, 66, 67, 68, 69, 70, 71, 72}};

            moris::mtk::MtkSetsInfo tMtkMeshSets;
            // Define side sets on the integration mesh (i.e. fixed bc, interface and ghost)

            // interface side set (all placed in a single cluster)
            moris::mtk::MtkSideSetInfo tInterfaceSideSet;
            Matrix<IndexMat> tInterfaceElemIdandSideOrd ( {{ 6, 2}, { 8, 1}, {10, 2}, {12, 1}, {14, 1}, {16, 2}, {17, 1},
                                                           {18, 1}, {20, 2}, {32, 2}, {33, 1}, {43, 1}, {44, 1}, {54, 1},
                                                           {55, 1}, {63, 1}, {73, 1}, {76, 1}, {79, 1}, {82, 2}});

            tInterfaceSideSet.mElemIdsAndSideOrds = &tInterfaceElemIdandSideOrd;
            tInterfaceSideSet.mSideSetName        = "iside" ;
            tMtkMeshSets.add_side_set(&tInterfaceSideSet);

            // Fixed bc (trivial side cluster)
            moris::mtk::MtkSideSetInfo tFixed;
            Matrix<IndexMat> tFixedElementsAndOrds = {{1, 4}};
            tFixed.mElemIdsAndSideOrds = &tFixedElementsAndOrds;
            tFixed.mSideSetName        = "fixed" ;
            tMtkMeshSets.add_side_set(&tFixed);

            // Fixed bc (trivial side cluster)
            moris::mtk::MtkSideSetInfo tGhost;
            Matrix<IndexMat> tGhostCellAndOrds = {{3,5},{4,4}};
            tGhost.mElemIdsAndSideOrds = &tGhostCellAndOrds;
            tGhost.mSideSetName        = "ghost_facets" ;
            tMtkMeshSets.add_side_set(&tGhost);

            // add block sets (Still in the mesh but not tested here)
            // Tet Cells in Omega 0
            moris::mtk::MtkBlockSetInfo tOmega0BlockSetTet;
            tOmega0BlockSetTet.mCellIdsInSet = &tCellIdsPhase0;
            tOmega0BlockSetTet.mBlockSetName = "Omega_0_tets";
            tOmega0BlockSetTet.mBlockSetTopo = mtk::CellTopology::TET4;
            tMtkMeshSets.add_block_set(&tOmega0BlockSetTet);

            // Hex Cells in Omega 0
            Matrix<IdMat> tOmega0HexCellIds = {{1,2,3}};
            moris::mtk::MtkBlockSetInfo tOmega0BlockSetHex;
            tOmega0BlockSetHex.mCellIdsInSet = &tOmega0HexCellIds;
            tOmega0BlockSetHex.mBlockSetName = "Omega_0_hex";
            tOmega0BlockSetHex.mBlockSetTopo = mtk::CellTopology::HEX8;
            tMtkMeshSets.add_block_set(&tOmega0BlockSetHex);

            // Cells in the ghost domain of omega 1
            moris::mtk::MtkBlockSetInfo tOmega0GhostBlockSetTet;
            tOmega0GhostBlockSetTet.mCellIdsInSet = &tCellIdsGhost0;
            tOmega0GhostBlockSetTet.mBlockSetName = "Omega_0_Ghost";
            tOmega0GhostBlockSetTet.mBlockSetTopo = mtk::CellTopology::TET4;
            tMtkMeshSets.add_block_set(&tOmega0GhostBlockSetTet);

            // Integration Cells for Ghost penalization only
            Matrix<IdMat> tGhostCellIds = {{3,4}};
            moris::mtk::MtkBlockSetInfo tCellsForGhost;
            tCellsForGhost.mCellIdsInSet = &tGhostCellIds;
            tCellsForGhost.mBlockSetName = "Ghost_Cells_0";
            tCellsForGhost.mBlockSetTopo = mtk::CellTopology::HEX8;
            tMtkMeshSets.add_block_set(&tCellsForGhost);

            // Mesh data input structure
            moris::mtk::MtkMeshData tMeshDataInput(3);

            moris::uint tSpatialDim   = 3;
            Matrix<IdMat> tNodeOwner(1,tNodeCoordinates.n_rows(),moris::par_rank());
            tMeshDataInput.ElemConn(0)             = &tInterpElemsAsIntegCellToNodes;
            tMeshDataInput.ElemConn(1)             = &tCellToNodePhase0;
            tMeshDataInput.ElemConn(2)             = &tCellToNodeGhost0;
            tMeshDataInput.LocaltoGlobalElemMap(0) = (&tInterpElemsAsIntegCellIds);
            tMeshDataInput.LocaltoGlobalElemMap(1) = (&tCellIdsPhase0);
            tMeshDataInput.LocaltoGlobalElemMap(2) = (&tCellIdsGhost0);
            tMeshDataInput.CreateAllEdgesAndFaces  = true;
            tMeshDataInput.Verbose                 = false;
            tMeshDataInput.SpatialDim              = &tSpatialDim;
            tMeshDataInput.NodeCoords              = &tNodeCoordinates;
            tMeshDataInput.NodeProcOwner           = &tNodeOwner;
            tMeshDataInput.LocaltoGlobalNodeMap    = &tLocalToGlobalNodeMap;
            tMeshDataInput.SetsInfo                = &tMtkMeshSets;
            tMeshDataInput.MarkNoBlockForIO        = false;

            // ---------------------------------------
            // CELL CLUSTERING
            // ---------------------------------------

            // Get the 4th cell from interpolation mesh and add a cell cluster to it)
            moris::moris_id tInterpCellIndex = 3;
            moris::mtk::Cell* tInterpCell = &tInterpMesh1->get_mtk_cell(tInterpCellIndex);

            // setup integration mesh
            // Cells and cell topology in material phase 0
            // Tetrathedral cells in material phase 1
            Matrix<IndexMat> tCellIdsCluster1Material = {{6, 8, 10, 12, 14, 16, 17, 18, 20, 31, 32, 33, 42, 43, 44, 53, 54, 55, 62, 63, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84}};

            // Tetrathedral cells in material phase 0 void
            Matrix<IndexMat> tCellIdsCluster1Void = {{}};

            // local coordinates
            // element level parameter
            Matrix<IdMat> tVertexIDsInCluster  = {{13, 14, 16, 15, 17, 18, 20, 19, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44}};
            Matrix<DDRMat> tLocalCoordinatesWrtInterpCell = {{-1, -1, -1},{ 1, -1, -1},{ 1, 1, -1},{-1, 1, -1},{-1, -1, 1},{ 1, -1, 1},{ 1, 1, 1},{-1, 1, 1},{ 0, -1, 0},{ 1, 0, 0},{ 0, 1, 0},{ -1, 0, 0},{0, 0, -1},{0, 0, 1},{0, 0, 0},{-0.8, -1, -0.8},{0.8, -1, -0.8},{-0.8, -0.8, -0.8},{0.8, -0.8, -0.8},{1, -1, -0.8},{-1, -1, -0.8},{1, -0.8, -0.8},{1, 0.8, -0.8},{0.8, 0.8, -0.8},{1, 1, -0.8},{0.8, 1, -0.8},{-0.8, 1, -0.8},{-0.8, 0.8, -0.8},{-1, 1, -0.8},{-1, 0.8, -0.8},{-1, -0.8, -0.8},{0, 0, -0.8}};

            // setup cluster data
            moris::mtk::Cell_Cluster_Input tCellClusterInput;
            tCellClusterInput.add_cluster_data(tInterpCell,&tCellIdsCluster1Material,&tCellIdsCluster1Void,&tVertexIDsInCluster,&tLocalCoordinatesWrtInterpCell);

            // add cluster to input data
            tMeshDataInput.CellClusterInput = &tCellClusterInput;

            // ---------------------------------------
            // SIDE CLUSTERING
            // NOTE: Only add non-trivial side clusters to this data structure
            // ---------------------------------------

            moris::mtk::Side_Cluster_Input tSideClusterInput;

            // register interface side set (and in turn get the index of this side set back)
            moris_index tInterfaceOrd = tSideClusterInput.add_side_set_label(tInterfaceSideSet.mSideSetName);

            Matrix<IdMat> tInterfaceVertexIDsInCluster             = {{28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44}};
            Matrix<DDRMat> tInterfaceLocalCoordinatesWrtInterpCell = {{-0.8, -1, -0.8},{0.8, -1, -0.8},{-0.8, -0.8, -0.8},{0.8, -0.8, -0.8},{1, -1, -0.8},{-1, -1, -0.8},{1, -0.8, -0.8},{1, 0.8, -0.8},{0.8, 0.8, -0.8},{1, 1, -0.8},{0.8, 1, -0.8},{-0.8, 1, -0.8},{-0.8, 0.8, -0.8},{-1, 1, -0.8},{-1, 0.8, -0.8},{-1, -0.8, -0.8},{0, 0, -0.8}};

            //moris::print(tInterfaceElemIdandSideOrd,"tInterfaceElemIdandSideOrd");
            tSideClusterInput.add_cluster_data(false,tInterfaceOrd,tInterpCell,&tInterfaceElemIdandSideOrd,&tInterfaceVertexIDsInCluster,&tInterfaceLocalCoordinatesWrtInterpCell);

            // add cluster to input data
            tMeshDataInput.SideClusterInput = & tSideClusterInput;

            // ---------------------------------------
            // DOUBLE SIDE CLUSTERING
            // ---------------------------------------
            moris::mtk::Double_Side_Cluster_Input tDoubleSideClusterInput;
            moris_index tOrd = tDoubleSideClusterInput.add_double_side_set_label(tGhost.mSideSetName);

            moris::Matrix<moris::IdMat> tDummyVerts(0,0);
            moris::Matrix<moris::DDRMat> tDummyCoords(0,0);

            // left side cluster
            moris_id         tLeftGhostInterpCellId   = 3;
            moris_index      tLeftGhostInterpCellInd  = tInterpMesh1->get_loc_entity_ind_from_entity_glb_id(tLeftGhostInterpCellId,mtk::EntityRank::ELEMENT);
            moris::mtk::Cell* tLeftInterpCell         = &tInterpMesh1->get_mtk_cell(tLeftGhostInterpCellInd);
            Matrix<IndexMat> tLeftGhostCellIdAndOrd   = {{3,5}};
            bool             tLeftTrivial = true;

            // right side cluster
            moris::mtk::Cell* tRightInterpCell         = &tInterpMesh1->get_mtk_cell(tLeftGhostInterpCellInd);
            Matrix<IndexMat> tRightGhostCellIdAndOrd = {{4,4}};
            bool             tRightTrivial = true;

            // Vertex pairing
            Matrix<IdMat> tVertexPair = {{13,13}, {14,14}, {15,15}, {16,16}};

            tDoubleSideClusterInput.add_cluster_data(tOrd,tLeftTrivial,tLeftInterpCell,&tLeftGhostCellIdAndOrd,&tDummyVerts,&tDummyCoords,
                                                     tRightTrivial,tRightInterpCell,&tRightGhostCellIdAndOrd,&tDummyVerts,&tDummyCoords, &tVertexPair);

            tMeshDataInput.DoubleSideClusterInput = &tDoubleSideClusterInput;

            moris::mtk::Integration_Mesh* tIntegMesh1  = moris::mtk::create_integration_mesh(mtk::MeshType::STK,tMeshDataInput,tInterpMesh1);

            // place the pair in mesh manager
            std::shared_ptr< mtk::Mesh_Manager > tMeshManager = std::make_shared< mtk::Mesh_Manager >();
            tMeshManager->register_mesh_pair(tInterpMesh1,tIntegMesh1);

            //------------------------------------------------------------------------------
            // create the properties
            std::shared_ptr< fem::Property > tPropConductivity = std::make_shared< fem::Property >();
            tPropConductivity->set_parameters( { {{ 1.0 }} } );
            tPropConductivity->set_val_function( tConstValFunction_MDLCUT );

            std::shared_ptr< fem::Property > tPropDirichlet = std::make_shared< fem::Property >();
            tPropDirichlet->set_parameters( { {{ 5.0 }} } );
            tPropDirichlet->set_val_function( tConstValFunction_MDLCUT );

            std::shared_ptr< fem::Property > tPropNeumann = std::make_shared< fem::Property >();
            tPropNeumann->set_parameters( { {{ 20.0 }} } );
            tPropNeumann->set_val_function( tConstValFunction_MDLCUT );

            std::shared_ptr< fem::Property > tPropTempLoad = std::make_shared< fem::Property >();
            tPropTempLoad->set_parameters( { {{ 0.0 }} } );
            tPropTempLoad->set_val_function( tConstValFunction_MDLCUT );

            // define constitutive models
            fem::CM_Factory tCMFactory;

            std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIso =
                    tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
            tCMDiffLinIso->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tCMDiffLinIso->set_property( tPropConductivity, "Conductivity" );
            tCMDiffLinIso->set_space_dim( 3 );
            tCMDiffLinIso->set_local_properties();

            // define stabilization parameters
            fem::SP_Factory tSPFactory;
            std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche = tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
            tSPDirichletNitsche->set_parameters( { {{ 1.0 }} } );
            tSPDirichletNitsche->set_property( tPropConductivity, "Material", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::Stabilization_Parameter > tSPGhost = tSPFactory.create_SP( fem::Stabilization_Type::GHOST_DISPL );
            tSPGhost->set_parameters( {{{ 1.0 }} });
            tSPGhost->set_property( tPropConductivity, "Material", mtk::Leader_Follower::LEADER );

            // define the IWGs
            fem::IWG_Factory tIWGFactory;

            std::shared_ptr< fem::IWG > tIWGBulk = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
            tIWGBulk->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGBulk->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tIWGBulk->set_constitutive_model( tCMDiffLinIso, "Diffusion", mtk::Leader_Follower::LEADER );
            tIWGBulk->set_property( tPropTempLoad, "Load", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::IWG > tIWGDirichlet = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
            tIWGDirichlet->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGDirichlet->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tIWGDirichlet->set_stabilization_parameter( tSPDirichletNitsche, "DirichletNitsche" );
            tIWGDirichlet->set_constitutive_model( tCMDiffLinIso, "Diffusion", mtk::Leader_Follower::LEADER );
            tIWGDirichlet->set_property( tPropDirichlet, "Dirichlet", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::IWG > tIWGNeumann = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_NEUMANN );
            tIWGNeumann->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGNeumann->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tIWGNeumann->set_property( tPropNeumann, "Neumann", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::IWG > tIWGGhost = tIWGFactory.create_IWG( fem::IWG_Type::GHOST_NORMAL_FIELD );
            tIWGGhost->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGGhost->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tIWGGhost->set_dof_type_list( {{ MSI::Dof_Type::TEMP }}, mtk::Leader_Follower::FOLLOWER );
            tIWGGhost->set_stabilization_parameter( tSPGhost, "GhostSP" );

//            uint iInterpOrder = 1;
//            if ( iInterpOrder > 0 )
//            {
//                std::shared_ptr< fem::Stabilization_Parameter > tSP1 = tSPFactory.create_SP( fem::Stabilization_Type::GHOST_DISPL );
//                tSP1->set_parameters( {{{ 1.0 }}, {{ 1.0 }} });
//                tSP1->set_property( tPropConductivity, "Material", mtk::Leader_Follower::LEADER );
//                tIWGGhost->set_stabilization_parameter( tSP1, "GhostDisplOrder1" );
//            }
//            if ( iInterpOrder > 1 )
//            {
//                std::shared_ptr< fem::Stabilization_Parameter > tSP2 = tSPFactory.create_SP( fem::Stabilization_Type::GHOST_DISPL );
//                tSP2->set_parameters( {{{ 1.0 }}, {{ 2.0 }} });
//                tSP2->set_property( tPropConductivity, "Material", mtk::Leader_Follower::LEADER );
//                tIWGGhost->set_stabilization_parameter( tSP2, "GhostDisplOrder2" );
//            }
//            if ( iInterpOrder > 2 )
//            {
//                std::shared_ptr< fem::Stabilization_Parameter > tSP3 = tSPFactory.create_SP( fem::Stabilization_Type::GHOST_DISPL );
//                tSP3->set_parameters( {{{ 1.0 }}, {{ 3.0 }} });
//                tSP3->set_property( tPropConductivity, "Material", mtk::Leader_Follower::LEADER );
//                tIWGGhost->set_stabilization_parameter( tSP3, "GhostDisplOrder3" );
//            }

            // define set info
            fem::Set_User_Info tSetBulk1;
            tSetBulk1.set_mesh_set_name( "Omega_0_hex" );
            tSetBulk1.set_IWGs( { tIWGBulk } );

            fem::Set_User_Info tSetBulk2;
            tSetBulk2.set_mesh_set_name( "Omega_0_tets" );
            tSetBulk2.set_IWGs( { tIWGBulk } );

            fem::Set_User_Info tSetDirichlet;
            tSetDirichlet.set_mesh_set_name( "fixed" );
            tSetDirichlet.set_IWGs( { tIWGDirichlet } );

            fem::Set_User_Info tSetNeumann;
            tSetNeumann.set_mesh_set_name( "iside" );
            tSetNeumann.set_IWGs( { tIWGNeumann } );

            fem::Set_User_Info tSetGhost;
            tSetGhost.set_mesh_set_name( tGhost.mSideSetName );
            tSetGhost.set_IWGs( { tIWGGhost } );

            // create a cell of set info
            Vector< fem::Set_User_Info > tSetInfo( 5 );
            tSetInfo( 0 ) = tSetBulk1;
            tSetInfo( 1 ) = tSetBulk2;
            tSetInfo( 2 ) = tSetDirichlet;
            tSetInfo( 3 ) = tSetNeumann;
            tSetInfo( 4 ) = tSetGhost;

            // create model
            mdl::Model * tModel = new mdl::Model(  tMeshManager,
                                                   1,
                                                   tSetInfo );

            Vector< enum MSI::Dof_Type > tDofTypes1( 1, MSI::Dof_Type::TEMP );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 1: create linear solver and algorithm
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            dla::Solver_Factory  tSolFactory;
            std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( sol::SolverType::AZTEC_IMPL );

            tLinearSolverAlgorithm->set_param("AZ_diagnostics") = AZ_none;
            tLinearSolverAlgorithm->set_param("AZ_output") = AZ_none;

            dla::Linear_Solver tLinSolver;

            tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 2: create nonlinear solver and algorithm
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            NLA::Nonlinear_Solver_Factory tNonlinFactory;
            std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );

            tNonlinearSolverAlgorithm->set_param("NLA_max_iter")   = 10;
            tNonlinearSolverAlgorithm->set_param("NLA_hard_break") = false;
            tNonlinearSolverAlgorithm->set_param("NLA_max_lin_solver_restarts") = 2;
            tNonlinearSolverAlgorithm->set_param("NLA_rebuild_jacobian") = true;

            tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );

            NLA::Nonlinear_Solver tNonlinearSolver;

            tNonlinearSolver.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 3: create time Solver and algorithm
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            tsa::Time_Solver_Factory tTimeSolverFactory;
            std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm = tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );

            tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolver );

            tsa::Time_Solver tTimeSolver;

            tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );

            sol::SOL_Warehouse tSolverWarehouse;

            tSolverWarehouse.set_solver_interface(tModel->get_solver_interface());

            tNonlinearSolver.set_solver_warehouse( &tSolverWarehouse );
            tTimeSolver.set_solver_warehouse( &tSolverWarehouse );

            tNonlinearSolver.set_dof_type_list( tDofTypes1 );
            tTimeSolver.set_dof_type_list( tDofTypes1 );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 4: Solve and check
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            tTimeSolver.solve();

            moris::Matrix< DDRMat > tSolution11;
            tTimeSolver.get_full_solution( tSolution11 );

            //delete tInterpMesh1;
            delete tIntegMesh1;
            delete tModel;
        }
        }

    }/* namespace mdl */
}/* namespace moris */

