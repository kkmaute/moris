/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_MSI_Element_Time_Sideset.cpp
 *
 */

#include "catch.hpp"
#include "fn_equal_to.hpp"
#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Communication_Tools.hpp"

#define protected public
#define private   public
#include "cl_MSI_Equation_Object.hpp"
#include "cl_MSI_Equation_Model.hpp"
#include "cl_MSI_Node_Proxy.hpp"
#include "cl_MSI_Model_Solver_Interface.hpp"
#include "cl_MSI_Dof_Manager.hpp"
#include "cl_MSI_Pdof_Host.hpp"
#include "cl_MTK_Cluster.hpp"
#include "cl_FEM_Cluster.hpp"
#include "cl_FEM_Interpolation_Element.hpp"
#include "cl_FEM_IWG.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#undef protected
#undef private

//MTK
#include "cl_MTK_Vertex.hpp"
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
//FEM/MDL
#include "cl_MDL_Model.hpp"
//FEM/MSI
#include "cl_MSI_Solver_Interface.hpp"
#include "MSI_Test_Proxy/cl_MTK_Vertex_Proxy.hpp"
#include "MSI_Test_Proxy/cl_MTK_Cell_Proxy.hpp"
#include "MSI_Test_Proxy/cl_MTK_Cluster_Proxy.hpp"
//FEM/INT
#include "cl_FEM_CM_Factory.hpp"
#include "cl_FEM_IWG_Factory.hpp"
#include "cl_FEM_IQI_Factory.hpp"
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_FEM_Geometry_Interpolator.hpp"
#include "cl_FEM_Property.hpp"
//SOL
#include "cl_SOL_Matrix_Vector_Factory.hpp"
#include "cl_SOL_Dist_Map.hpp"
#include "cl_SOL_Dist_Vector.hpp"

namespace moris
{
namespace MSI
{

void tConstValFunction_UTTimeElem
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 );
}

TEST_CASE("Element_Time_Sideset","[INT],[Element_Time_Sideset]")
{
    if(par_size() == 1 )
    {
        // setup the integration mesh
        // Define the Integration Mesh (from data from xtk)
        Matrix<DDRMat>   tNodeCoordinatesInterpolation ={{0, 0},{1, 0},{1, 1},{0, 1}};
        Matrix<IndexMat> tLocalToGlobalNodeMapInterpolation = {{1, 2, 3, 4 }};

        Matrix<IndexMat> tInterpElemsAsIntegCellIdsInterpolation     = {{1}};
        Matrix<IndexMat> tInterpElemsAsIntegCellToNodesInterpolation = {{1, 2, 3, 4}};

        // Tetrathedral cells in material phase 1
        CellTopology     tPhase0ChildTopoInterpolation  = CellTopology::QUAD4;
        Matrix<IndexMat> tCellIdsPhase0Interpolation    = {{1}};
        Matrix<IndexMat> tCellToNodePhase0Interpolation= {{1, 2, 3,4}};

        moris::mtk::MtkSetsInfo tMtkMeshSetsInterpolation;

        // add block sets (Still in the mesh but not tested here)
        // Tet Cells in Omega 0
        moris::mtk::MtkBlockSetInfo tOmega0BlockSetQuad;
        tOmega0BlockSetQuad.mCellIdsInSet = &tCellIdsPhase0Interpolation;
        tOmega0BlockSetQuad.mBlockSetName = "Omega_tets";
        tOmega0BlockSetQuad.mBlockSetTopo = CellTopology::QUAD4;
        tMtkMeshSetsInterpolation.add_block_set(&tOmega0BlockSetQuad);

        // Mesh data input structure
        moris::mtk::MtkMeshData tMeshDataInputInterpolation(1);

        moris::uint tSpatialDim   = 2;
        moris::uint tNumElemTypes = 1;
        Matrix<IdMat> tNodeOwnerInterpolation(1,tNodeCoordinatesInterpolation.n_rows(),moris::par_rank());
        tMeshDataInputInterpolation.ElemConn(0)             = &tInterpElemsAsIntegCellToNodesInterpolation;
        tMeshDataInputInterpolation.LocaltoGlobalElemMap(0) = (&tInterpElemsAsIntegCellIdsInterpolation);
        tMeshDataInputInterpolation.CreateAllEdgesAndFaces  = true;
        tMeshDataInputInterpolation.Verbose                 = false;
        tMeshDataInputInterpolation.SpatialDim              = &tSpatialDim;
        tMeshDataInputInterpolation.NodeCoords              = &tNodeCoordinatesInterpolation;
        tMeshDataInputInterpolation.NodeProcOwner           = &tNodeOwnerInterpolation;
        tMeshDataInputInterpolation.LocaltoGlobalNodeMap    = &tLocalToGlobalNodeMapInterpolation;
        tMeshDataInputInterpolation.SetsInfo                = &tMtkMeshSetsInterpolation;
        tMeshDataInputInterpolation.MarkNoBlockForIO        = false;

        moris::mtk::Interpolation_Mesh* tInterpMesh1 = moris::mtk::create_interpolation_mesh( MeshType::STK, tMeshDataInputInterpolation );

        //--------------------------------------------------------------------------------------------------------

        // setup the integration mesh
        // Define the Integration Mesh (from data from xtk)
        Matrix<DDRMat>   tNodeCoordinates = {{0, 0},{1, 0},{1, 1},{0, 1},{0.5, 0},{0, 0.5}};
        Matrix<IndexMat> tLocalToGlobalNodeMap = {{1, 2, 3, 4, 5, 6 }};

        // Tetrathedral cells in material phase 1
        CellTopology     tPhase0ChildTopo  = CellTopology::TRI3;
        Matrix<IndexMat> tCellIdsPhase0    = {{1,2,3,4}};
        Matrix<IndexMat> tCellToNodePhase0 = {{1,5,6},{2,3,5},{5,3,6},{6,3,4}};

//        CellTopology     tPhase0ChildTopo  = CellTopology::TRI3;
//        Matrix<IndexMat> tCellIdsPhase1    = {{2, 3, 4}};
//        Matrix<IndexMat> tCellToNodePhase1 = {{1,2,4},{2,5,4},{5,2,3}};

        moris::mtk::MtkSetsInfo tMtkMeshSets;

        // add block sets (Still in the mesh but not tested here)
        // Tet Cells in Omega 0
        moris::mtk::MtkBlockSetInfo tOmega0BlockSetTet;
        tOmega0BlockSetTet.mCellIdsInSet = &tCellIdsPhase0;
        tOmega0BlockSetTet.mBlockSetName = "Omega_0_tets";
        tOmega0BlockSetTet.mBlockSetTopo = CellTopology::TRI3;
        tMtkMeshSets.add_block_set(&tOmega0BlockSetTet);

        // Mesh data input structure
        moris::mtk::MtkMeshData tMeshDataInput(1);

        Matrix<IdMat> tNodeOwner(1,tNodeCoordinates.n_rows(),moris::par_rank());
        tMeshDataInput.ElemConn(0)             = &tCellToNodePhase0;
        tMeshDataInput.LocaltoGlobalElemMap(0) = (&tCellIdsPhase0);
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

        // Get the  cell from interpolation mesh and add a cell cluster to it)
        moris::moris_id tInterpCellIndex = 0;
        moris::mtk::Cell* tInterpCell = &tInterpMesh1->get_mtk_cell(tInterpCellIndex);

        // setup integration mesh
        // Cells and cell topology in material phase 0
        // Triangle cells in material phase 1
        Matrix<IndexMat> tCellIdsCluster1Material = {{1}};

        // Triangle cells in material phase 0 void
        Matrix<IndexMat> tCellIdsCluster1Void = {{2,3,4}};

        // local coordinates
        // element level parameter
        Matrix<IdMat> tVertexIDsInCluster  = {{ 1,2,3,4,5,6 }};
        Matrix<DDRMat> tLocalCoordinatesWrtInterpCell = {{-1, -1},{ 1, -1},{ 1, 1},{-1, 1},{0, -1},{-1, 0}};

        // setup cluster data
        moris::mtk::Cell_Cluster_Input tCellClusterInput;
        tCellClusterInput.add_cluster_data(tInterpCell,&tCellIdsCluster1Material,&tCellIdsCluster1Void,&tVertexIDsInCluster,&tLocalCoordinatesWrtInterpCell);

        // add cluster to input data
        tMeshDataInput.CellClusterInput = &tCellClusterInput;

        moris::mtk::Integration_Mesh* tIntegMesh1  = moris::mtk::create_integration_mesh(MeshType::STK,tMeshDataInput,tInterpMesh1);

        // place the pair in mesh manager
        std::shared_ptr< mtk::Mesh_Manager > tMeshManager = std::make_shared< mtk::Mesh_Manager >();
        tMeshManager->register_mesh_pair(tInterpMesh1,tIntegMesh1);

        // FEM inputs
        //------------------------------------------------------------------------------
        // create the properties
        std::shared_ptr< fem::Property > tPropWeightCurrent = std::make_shared< fem::Property > ();
        tPropWeightCurrent->set_parameters( { {{ 1.0 }} } );
        tPropWeightCurrent->set_val_function( tConstValFunction_UTTimeElem );

        std::shared_ptr< fem::Property > tPropWeightPrevious = std::make_shared< fem::Property > ();
        tPropWeightPrevious->set_parameters( { {{ 1.0 }} } );
        tPropWeightPrevious->set_val_function( tConstValFunction_UTTimeElem );

        // define an IWG
        fem::IWG_Factory tIWGFactory;

        std::shared_ptr< fem::IWG > tIWG = tIWGFactory.create_IWG( fem::IWG_Type::TIME_CONTINUITY_DOF );
        tIWG->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
        tIWG->set_dof_type_list( {{ MSI::Dof_Type::TEMP }}, mtk::Master_Slave::MASTER );
        tIWG->set_property( tPropWeightCurrent, "WeightCurrent", mtk::Master_Slave::MASTER );
        tIWG->set_property( tPropWeightPrevious, "WeightPrevious", mtk::Master_Slave::MASTER );

        // define set info
        moris::Cell< fem::Set_User_Info > tSetInfo( 1 );
        tSetInfo( 0 ).set_mesh_set_name( "Omega_0_tets" );
        tSetInfo( 0 ).set_time_continuity( true );
        tSetInfo( 0 ).set_IWGs( { tIWG } );

        // FEM model
        //------------------------------------------------------------------------------
        // create model
        mdl::Model * tModel = new mdl::Model( tMeshManager,
                                               0,
                                               tSetInfo );

        MSI::MSI_Solver_Interface * tSolverInterface = tModel->get_solver_interface();

        sol::Matrix_Vector_Factory tMatFactory( sol::MapType::Epetra );
        sol::Dist_Map*  mVectorMap = tMatFactory.create_map( {{0},{1},{2},{3},{4},{5},{6},{7}}, {{}} );
        sol::Dist_Vector * mVector = tMatFactory.create_vector( nullptr, mVectorMap, 1 );
        sol::Dist_Vector * mPreviousVector = tMatFactory.create_vector( nullptr, mVectorMap, 1 );

        mVector->sum_into_global_values(         {{0},{1},{2},{3},{4},{5},{6},{7}},
                                                 {{1},{2},{3},{4},{5},{6},{7},{8}});
        mPreviousVector->sum_into_global_values( {{0},{1},{2},{3},{4},{5},{6},{7}},
                                                 {{1},{2},{3},{4},{5},{6},{7},{8}});

        // FEM set
        //------------------------------------------------------------------------------
        // get the equation set from the model
        moris::Cell< MSI::Equation_Set * > tSets
        = tModel->get_fem_model()->get_equation_sets();

        // get the equation set from the model
        std::shared_ptr< MSI::Equation_Model > tEquationModel
        = tModel->get_fem_model();

        // get a working set
        MSI::Equation_Set* tWorkSet = tSets( 0 );

        // set the IWG and IQI to the set
        reinterpret_cast< fem::Set * >( tWorkSet )->mRequestedIWGs = { tIWG };

        // set size and fill the set residual assembly map
        reinterpret_cast< fem::Set * >( tWorkSet )->mResDofAssemblyMap.resize( 1 );
        reinterpret_cast< fem::Set * >( tWorkSet )->mResDofAssemblyMap( 0 ) = { { 0, 7 } };

        // set size and fill the set jacobian assembly map
        reinterpret_cast< fem::Set * >( tWorkSet )->mJacDofAssemblyMap.resize( 1 );
        reinterpret_cast< fem::Set * >( tWorkSet )->mJacDofAssemblyMap( 0 ) = { { 0, 7 } };

        // MSI Equation object
        //------------------------------------------------------------------------------
        // get list of equation objects from set
        Cell< MSI::Equation_Object * > tEqObjs = tWorkSet->get_equation_object_list();

        // get a working equation object
        MSI::Equation_Object * tWorkEqObj = tEqObjs( 0 );

        // set the solution vector
        tEquationModel->set_solution_vector( mVector );
        tEquationModel->set_previous_solution_vector( mPreviousVector );

        // set the time
        Matrix< DDRMat > tTime = { { 1.0 }, { 2.0 } };
        Matrix< DDRMat > tPreviousTime = { { 0.0 }, { 1.0 } };
        tEquationModel->set_time( tTime );
        tEquationModel->set_previous_time( tPreviousTime );

        // Init IWG and IQI for forward analysis
        //------------------------------------------------------------------------------
        // set the IWG fem set
        tIWG->set_set_pointer( reinterpret_cast< fem::Set* >( tWorkSet ) );

        // build global dof type list
        tIWG->get_global_dof_type_list();
        tIWG->get_global_dv_type_list();

        // populate the requested master dof type
        tIWG->mRequestedMasterGlobalDofTypes = {{ MSI::Dof_Type::TEMP }};

        // init residual
        tWorkSet->mResidual.resize( 1 );
        tWorkSet->mResidual( 0 ).set_size( 8, 1, 0.0 );
        tWorkSet->mResidualExist = true;

        // init jacobian
        tWorkSet->mJacobian.set_size( 8, 8, 0.0 );
        tWorkSet->mJacobianExist = true;

        // compute residual
        tWorkEqObj->compute_residual();
        print( tWorkSet->get_residual()( 0 ), "R" );

        // compute jacobian
        tWorkEqObj->compute_jacobian();
        print( tWorkSet->get_jacobian(), "dRdu" );

        delete mVectorMap;
    }

}/* END_TEST_CASE */

TEST_CASE("Element_Time_Sideset_2","[INT],[Element_Time_Sideset_2]")
{
    if(par_size() == 1 )
    {
        // setup the integration mesh
        // Define the Integration Mesh (from data from xtk)
        Matrix<DDRMat>   tNodeCoordinatesInterpolation ={{0, 0},{1, 0},{1, 1},{0, 1}};
        Matrix<IndexMat> tLocalToGlobalNodeMapInterpolation = {{1, 2, 3, 4 }};

        Matrix<IndexMat> tInterpElemsAsIntegCellIdsInterpolation     = {{1}};
        Matrix<IndexMat> tInterpElemsAsIntegCellToNodesInterpolation = {{1, 2, 3, 4}};

        // Tetrathedral cells in material phase 1
        CellTopology     tPhase0ChildTopoInterpolation  = CellTopology::QUAD4;
        Matrix<IndexMat> tCellIdsPhase0Interpolation    = {{1}};
        Matrix<IndexMat> tCellToNodePhase0Interpolation = {{1, 2, 3,4}};

        moris::mtk::MtkSetsInfo tMtkMeshSetsInterpolation;

        // add block sets (Still in the mesh but not tested here)
        // Tet Cells in Omega 0
        moris::mtk::MtkBlockSetInfo tOmega0BlockSetQuad;
        tOmega0BlockSetQuad.mCellIdsInSet = &tCellIdsPhase0Interpolation;
        tOmega0BlockSetQuad.mBlockSetName = "Omega_tets";
        tOmega0BlockSetQuad.mBlockSetTopo = CellTopology::QUAD4;
        tMtkMeshSetsInterpolation.add_block_set(&tOmega0BlockSetQuad);

        // Mesh data input structure
        moris::mtk::MtkMeshData tMeshDataInputInterpolation(1);

        moris::uint tSpatialDim   = 2;
        moris::uint tNumElemTypes = 1;
        Matrix<IdMat> tNodeOwnerInterpolation(1,tNodeCoordinatesInterpolation.n_rows(),moris::par_rank());
        tMeshDataInputInterpolation.ElemConn(0)             = &tInterpElemsAsIntegCellToNodesInterpolation;
        tMeshDataInputInterpolation.LocaltoGlobalElemMap(0) = (&tInterpElemsAsIntegCellIdsInterpolation);
        tMeshDataInputInterpolation.CreateAllEdgesAndFaces  = true;
        tMeshDataInputInterpolation.Verbose                 = false;
        tMeshDataInputInterpolation.SpatialDim              = &tSpatialDim;
        tMeshDataInputInterpolation.NodeCoords              = &tNodeCoordinatesInterpolation;
        tMeshDataInputInterpolation.NodeProcOwner           = &tNodeOwnerInterpolation;
        tMeshDataInputInterpolation.LocaltoGlobalNodeMap    = &tLocalToGlobalNodeMapInterpolation;
        tMeshDataInputInterpolation.SetsInfo                = &tMtkMeshSetsInterpolation;
        tMeshDataInputInterpolation.MarkNoBlockForIO        = false;

        moris::mtk::Interpolation_Mesh* tInterpMesh1 = moris::mtk::create_interpolation_mesh( MeshType::STK, tMeshDataInputInterpolation );

        //--------------------------------------------------------------------------------------------------------

        // setup the integration mesh
        // Define the Integration Mesh (from data from xtk)
        Matrix<DDRMat>   tNodeCoordinates = {{0, 0},{1, 0},{1, 1},{0, 1}};
        Matrix<IndexMat> tLocalToGlobalNodeMap = {{1, 2, 3, 4 }};

        // Tetrathedral cells in material phase 1
        CellTopology     tPhase0ChildTopo  = CellTopology::QUAD4;
        Matrix<IndexMat> tCellIdsPhase0    = {{1}};
        Matrix<IndexMat> tCellToNodePhase0 = {{1,2,3,4}};

        moris::mtk::MtkSetsInfo tMtkMeshSets;

        // add block sets (Still in the mesh but not tested here)
        // Tet Cells in Omega 0
        moris::mtk::MtkBlockSetInfo tOmega0BlockSetTet;
        tOmega0BlockSetTet.mCellIdsInSet = &tCellIdsPhase0;
        tOmega0BlockSetTet.mBlockSetName = "Omega_0_tets";
        tOmega0BlockSetTet.mBlockSetTopo = CellTopology::QUAD4;
        tMtkMeshSets.add_block_set(&tOmega0BlockSetTet);

        // Mesh data input structure
        moris::mtk::MtkMeshData tMeshDataInput(1);

        Matrix<IdMat> tNodeOwner(1,tNodeCoordinates.n_rows(),moris::par_rank());
        tMeshDataInput.ElemConn(0)             = &tCellToNodePhase0;
        tMeshDataInput.LocaltoGlobalElemMap(0) = (&tCellIdsPhase0);
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

        // Get the  cell from interpolation mesh and add a cell cluster to it)
        moris::moris_id tInterpCellIndex = 0;
        moris::mtk::Cell* tInterpCell = &tInterpMesh1->get_mtk_cell(tInterpCellIndex);

        // setup integration mesh
        // Cells and cell topology in material phase 0
        // Triangle cells in material phase 1
        Matrix<IndexMat> tCellIdsCluster1Material = {{1}};

        // Triangle cells in material phase 0 void
        Matrix<IndexMat> tCellIdsCluster1Void = {{}};

        // local coordinates
        // element level parameter
        Matrix<IdMat> tVertexIDsInCluster  = {{ 1,2,3,4 }};
        Matrix<DDRMat> tLocalCoordinatesWrtInterpCell = {{-1, -1},{ 1, -1},{ 1, 1},{-1, 1}};

        // setup cluster data
        moris::mtk::Cell_Cluster_Input tCellClusterInput;
        tCellClusterInput.add_cluster_data(tInterpCell,&tCellIdsCluster1Material,&tCellIdsCluster1Void,&tVertexIDsInCluster,&tLocalCoordinatesWrtInterpCell);

        // add cluster to input data
        tMeshDataInput.CellClusterInput = &tCellClusterInput;

        moris::mtk::Integration_Mesh* tIntegMesh1  = moris::mtk::create_integration_mesh(MeshType::STK,tMeshDataInput,tInterpMesh1);

        // place the pair in mesh manager
        mtk::Mesh_Manager tMeshManager;
        tMeshManager.register_mesh_pair(tInterpMesh1,tIntegMesh1);

        // FEM inputs
        //------------------------------------------------------------------------------
        // create the properties
        std::shared_ptr< fem::Property > tPropWeightCurrent = std::make_shared< fem::Property > ();
        tPropWeightCurrent->set_parameters( { {{ 1.0 }} } );
        tPropWeightCurrent->set_val_function( tConstValFunction_UTTimeElem );

        std::shared_ptr< fem::Property > tPropWeightPrevious = std::make_shared< fem::Property > ();
        tPropWeightPrevious->set_parameters( { {{ 1.0 }} } );
        tPropWeightPrevious->set_val_function( tConstValFunction_UTTimeElem );

        std::shared_ptr< fem::Property > tPropDensity = std::make_shared< fem::Property > ();
        tPropDensity->set_parameters( { {{ 1.0 }} } );
        tPropDensity->set_val_function( tConstValFunction_UTTimeElem );

        std::shared_ptr< fem::Property > tPropHeatCapacity = std::make_shared< fem::Property > ();
        tPropHeatCapacity->set_parameters( { {{ 1.0 }} } );
        tPropHeatCapacity->set_val_function( tConstValFunction_UTTimeElem );

        std::shared_ptr< fem::Property > tPropConductivity = std::make_shared< fem::Property > ();
        tPropConductivity->set_parameters( { {{ 1.0 }} } );
        tPropConductivity->set_val_function( tConstValFunction_UTTimeElem );

        std::shared_ptr< fem::Property > tPropInitCond = std::make_shared< fem::Property > ();
        tPropInitCond->set_parameters( { {{ 1.0 }} } );
        tPropInitCond->set_val_function( tConstValFunction_UTTimeElem );

        // define constitutive models
        fem::CM_Factory tCMFactory;

        std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIso =
                tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
        tCMDiffLinIso->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tCMDiffLinIso->set_property( tPropConductivity, "Conductivity" );
        tCMDiffLinIso->set_space_dim( 2 );
        tCMMasterDiffLinIso->set_local_properties();

        // define an IWG
        fem::IWG_Factory tIWGFactory;

        std::shared_ptr< fem::IWG > tIWGTime = tIWGFactory.create_IWG( fem::IWG_Type::TIME_CONTINUITY_DOF );
        tIWGTime->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
        tIWGTime->set_dof_type_list( { { MSI::Dof_Type::TEMP } }, mtk::Master_Slave::MASTER );
        tIWGTime->set_property( tPropWeightCurrent, "WeightCurrent", mtk::Master_Slave::MASTER );
        tIWGTime->set_property( tPropWeightPrevious, "WeightPrevious", mtk::Master_Slave::MASTER );
        tIWGTime->set_property( tPropInitCond, "InitialCondition", mtk::Master_Slave::MASTER );

        std::shared_ptr< fem::IWG > tIWGDiffusion = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
        tIWGDiffusion->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
        tIWGDiffusion->set_dof_type_list( { { MSI::Dof_Type::TEMP } }, mtk::Master_Slave::MASTER );
        tIWGDiffusion->set_constitutive_model( tCMDiffLinIso, "Diffusion", mtk::Master_Slave::MASTER );
        tIWGDiffusion->set_property( tPropDensity, "Density", mtk::Master_Slave::MASTER );
        tIWGDiffusion->set_property( tPropHeatCapacity, "HeatCapacity", mtk::Master_Slave::MASTER );

        // define set info
        moris::Cell< fem::Set_User_Info > tSetInfo( 2 );

        tSetInfo( 0 ).set_mesh_set_name( "Omega_0_tets" );
        tSetInfo( 0 ).set_time_continuity( true );
        tSetInfo( 0 ).set_IWGs( { tIWGTime } );

        tSetInfo( 1 ).set_mesh_set_name( "Omega_0_tets" );
        tSetInfo( 1 ).set_IWGs( { tIWGDiffusion } );

        // FEM model
        //------------------------------------------------------------------------------
        // create model
        mdl::Model * tModel = new mdl::Model( &tMeshManager,
                                               0,
                                               tSetInfo );

        MSI::MSI_Solver_Interface * tSolverInterface = tModel->get_solver_interface();

        sol::Matrix_Vector_Factory tMatFactory( sol::MapType::Epetra );
        sol::Dist_Map*  mVectorMap = tMatFactory.create_map( {{0},{1},{2},{3},{4},{5},{6},{7}}, {{}} );
        sol:: * mVector = tMatFactory.create_vector( nullptr, mVectorMap, 1 );
        sol::Dist_Vector * mPreviousVector = tMatFactory.create_vector( nullptr, mVectorMap, 1 );

        mVector->sum_into_global_values(         {{0},{1},{2},{3},{4},{5},{6},{7}},
                                                 {{1},{2},{3},{4},{5},{6},{7},{8}});
        mPreviousVector->sum_into_global_values( {{0},{1},{2},{3},{4},{5},{6},{7}},
                                                 {{1},{2},{3},{4},{5},{6},{7},{8}});

        // FEM set
        //------------------------------------------------------------------------------
        // get the equation set from the model
        moris::Cell< MSI::Equation_Set * > tSets
        = tModel->get_fem_model()->get_equation_sets();

        // get the equation set from the model
        std::shared_ptr< MSI::Equation_Model > tEquationModel
        = tModel->get_fem_model();

        // get a working set
        MSI::Equation_Set* tWorkSet0 = tSets( 0 );
        MSI::Equation_Set* tWorkSet1 = tSets( 1 );

        // set the IWG and IQI to the set
        reinterpret_cast< fem::Set * >( tWorkSet0 )->mRequestedIWGs = { tIWGTime };

        // set size and fill the set residual assembly map
        reinterpret_cast< fem::Set * >( tWorkSet0 )->mResDofAssemblyMap.resize( 1 );
        reinterpret_cast< fem::Set * >( tWorkSet0 )->mResDofAssemblyMap( 0 ) = { { 0, 7 } };

        // set size and fill the set jacobian assembly map
        reinterpret_cast< fem::Set * >( tWorkSet0 )->mJacDofAssemblyMap.resize( 1 );
        reinterpret_cast< fem::Set * >( tWorkSet0 )->mJacDofAssemblyMap( 0 ) = { { 0, 7 } };

        // set the IWG and IQI to the set
        reinterpret_cast< fem::Set * >( tWorkSet1 )->mRequestedIWGs = { tIWGDiffusion };

        // set size and fill the set residual assembly map
        reinterpret_cast< fem::Set * >( tWorkSet1 )->mResDofAssemblyMap.resize( 1 );
        reinterpret_cast< fem::Set * >( tWorkSet1 )->mResDofAssemblyMap( 0 ) = { { 0, 7 } };

        // set size and fill the set jacobian assembly map
        reinterpret_cast< fem::Set * >( tWorkSet1 )->mJacDofAssemblyMap.resize( 1 );
        reinterpret_cast< fem::Set * >( tWorkSet1 )->mJacDofAssemblyMap( 0 ) = { { 0, 7 } };

        // MSI Equation object
        //------------------------------------------------------------------------------
        // get list of equation objects from set
        Cell< MSI::Equation_Object * > tEqObjs0 = tWorkSet0->get_equation_object_list();
        Cell< MSI::Equation_Object * > tEqObjs1 = tWorkSet1->get_equation_object_list();

        // get a working equation object
        MSI::Equation_Object * tWorkEqObj0 = tEqObjs0( 0 );
        MSI::Equation_Object * tWorkEqObj1 = tEqObjs1( 0 );

        // set the solution vector
        tEquationModel->set_solution_vector( mVector );
        tEquationModel->set_previous_solution_vector( mPreviousVector );

        // set the time
        Matrix< DDRMat > tTime = { { 0.0 }, { 1.0 } };
        Matrix< DDRMat > tPreviousTime = { { 0.0 }, { 0.0 } };
        tEquationModel->set_time( tTime );
        tEquationModel->set_previous_time( tPreviousTime );

        // Init IWG and IQI for forward analysis
        //------------------------------------------------------------------------------
        // set the IWG fem set
        tIWGTime->set_set_pointer( reinterpret_cast< fem::Set* >( tWorkSet0 ) );

        // build global dof type list
        tIWGTime->get_global_dof_type_list();
        tIWGTime->get_global_dv_type_list();

        // populate the requested master dof type
        tIWGTime->mRequestedMasterGlobalDofTypes = {{ MSI::Dof_Type::TEMP }};

        // init residual
        tWorkSet0->mResidual.resize( 1 );
        tWorkSet0->mResidual( 0 ).set_size( 8, 1, 0.0 );
        tWorkSet0->mResidualExist = true;

        // init jacobian
        tWorkSet0->mJacobian.set_size( 8, 8, 0.0 );
        tWorkSet0->mJacobianExist = true;

        // compute residual
        tWorkEqObj0->compute_residual();
        print( tWorkSet0->get_residual()( 0 ), "RTime" );

        // compute jacobian
        tWorkEqObj0->compute_jacobian();
        print( tWorkSet0->get_jacobian(), "dRTimedu" );

        // Init IWG and IQI for forward analysis
        //------------------------------------------------------------------------------
        // set the IWG fem set
        tIWGDiffusion->set_set_pointer( reinterpret_cast< fem::Set* >( tWorkSet1 ) );

        // build global dof type list
        tIWGDiffusion->get_global_dof_type_list();
        tIWGDiffusion->get_global_dv_type_list();

        // populate the requested master dof type
        tIWGDiffusion->mRequestedMasterGlobalDofTypes = {{ MSI::Dof_Type::TEMP }};

        // init residual
        tWorkSet1->mResidual.resize( 1 );
        tWorkSet1->mResidual( 0 ).set_size( 8, 1, 0.0 );
        tWorkSet1->mResidualExist = true;

        // init jacobian
        tWorkSet1->mJacobian.set_size( 8, 8, 0.0 );
        tWorkSet1->mJacobianExist = true;

        // compute residual
        tWorkEqObj1->compute_residual();
        print( tWorkSet1->get_residual()( 0 ), "RDiff" );

        // compute jacobian
        tWorkEqObj1->compute_jacobian();
        print( tWorkSet1->get_jacobian(), "dRDiffdu" );

        delete mVectorMap;
    }

}/* END_TEST_CASE */

}/* end_namespace_msi */
}/* end_namespace_moris */

