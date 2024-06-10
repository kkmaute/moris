/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_Equation_Object_Pdv.cpp
 *
 */

#include "catch.hpp"
#include "fn_equal_to.hpp"
#include "moris_typedefs.hpp"
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
#include "MSI_Test_Proxy/cl_MSI_Design_Variable_Interface_Proxy.hpp"
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

void tConstValFunction_FDTest
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  Vector< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 );
}

void tFIValDvFunction_FDTest
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  Vector< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 ) * aFIManager->get_field_interpolators_for_type( gen::PDV_Type::DENSITY )->val();
}

void tFIDerDvFunction_FDTest
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  Vector< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 ) * aFIManager->get_field_interpolators_for_type( gen::PDV_Type::DENSITY )->N();
}

TEST_CASE("Eqn_Obj_pdv","[MSI],[Eqn_Obj_pdv]")
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
        Matrix<IndexMat> tCellIdsPhase0Interpolation    = {{1}};
        Matrix<IndexMat> tCellToNodePhase0Interpolation= {{1, 2, 3,4}};

        moris::mtk::MtkSetsInfo tMtkMeshSetsInterpolation;

        // add block sets (Still in the mesh but not tested here)
        // Tet Cells in Omega 0
        moris::mtk::MtkBlockSetInfo tOmega0BlockSetQuad;
        tOmega0BlockSetQuad.mCellIdsInSet = &tCellIdsPhase0Interpolation;
        tOmega0BlockSetQuad.mBlockSetName = "Omega_tets";
        tOmega0BlockSetQuad.mBlockSetTopo = mtk::CellTopology::QUAD4;
        tMtkMeshSetsInterpolation.add_block_set(&tOmega0BlockSetQuad);

        // Mesh data input structure
        moris::mtk::MtkMeshData tMeshDataInputInterpolation(1);

        moris::uint tSpatialDim   = 2;
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

        moris::mtk::Interpolation_Mesh* tInterpMesh1 = moris::mtk::create_interpolation_mesh( mtk::MeshType::STK, tMeshDataInputInterpolation );

        //--------------------------------------------------------------------------------------------------------

        // setup the integration mesh
        // Define the Integration Mesh (from data from xtk)
        Matrix<DDRMat>   tNodeCoordinates = {{0, 0},{1, 0},{1, 1},{0, 1},{0.5, 0},{0, 0.5}};
        Matrix<IndexMat> tLocalToGlobalNodeMap = {{1, 2, 3, 4, 5, 6 }};

        // Tetrathedral cells in material phase 1
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
        tOmega0BlockSetTet.mBlockSetTopo = mtk::CellTopology::TRI3;
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
        // Tetrathedral cells in material phase 1
        Matrix<IndexMat> tCellIdsCluster1Material = {{1}};

        // Tetrathedral cells in material phase 0 void
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

        moris::mtk::Integration_Mesh* tIntegMesh1  = moris::mtk::create_integration_mesh(mtk::MeshType::STK,tMeshDataInput,tInterpMesh1);

        // place the pair in mesh manager
        std::shared_ptr< mtk::Mesh_Manager > tMeshManager = std::make_shared< mtk::Mesh_Manager >();
        tMeshManager->register_mesh_pair(tInterpMesh1,tIntegMesh1);

        Design_Variable_Interface * tDesignVariableInterface = new Design_Variable_Interface_Proxy();

        // FEM inputs
        //------------------------------------------------------------------------------
        // create the properties
        std::shared_ptr< fem::Property > tPropLeaderConductivity =
                std::make_shared< fem::Property > ();
        tPropLeaderConductivity->set_parameters( { {{ 1.0 }} } );
        tPropLeaderConductivity->set_dv_type_list( {{ gen::PDV_Type::DENSITY }} );
        tPropLeaderConductivity->set_val_function( tFIValDvFunction_FDTest );
        tPropLeaderConductivity->set_dv_derivative_functions( { tFIDerDvFunction_FDTest } );

        std::shared_ptr< fem::Property > tPropLeaderTempLoad =
                std::make_shared< fem::Property > ();
        tPropLeaderTempLoad->set_parameters( { {{ 1.0 }} } );
        tPropLeaderTempLoad->set_val_function( tConstValFunction_FDTest );

        // define constitutive models
        fem::CM_Factory tCMFactory;

        std::shared_ptr< fem::Constitutive_Model > tCMLeaderDiffLinIso =
                tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
        tCMLeaderDiffLinIso->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tCMLeaderDiffLinIso->set_property( tPropLeaderConductivity, "Conductivity" );
        tCMLeaderDiffLinIso->set_space_dim( 2 );
        tCMLeaderDiffLinIso->set_local_properties();

        // define an IWG
        fem::IWG_Factory tIWGFactory;

        std::shared_ptr< fem::IWG > tIWG =
                tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
        tIWG->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
        tIWG->set_dof_type_list( {{ MSI::Dof_Type::TEMP }}, mtk::Leader_Follower::LEADER );
        tIWG->set_constitutive_model( tCMLeaderDiffLinIso, "Diffusion", mtk::Leader_Follower::LEADER );
        tIWG->set_property( tPropLeaderTempLoad, "Load", mtk::Leader_Follower::LEADER );

        // define an IQI
        fem::IQI_Factory tIQIFactory;
        std::shared_ptr< fem::IQI > tIQI =
                tIQIFactory.create_IQI( fem::IQI_Type::STRAIN_ENERGY );
        tIQI->set_dof_type_list( {{ MSI::Dof_Type::TEMP }}, mtk::Leader_Follower::LEADER );
        tIQI->set_constitutive_model( tCMLeaderDiffLinIso, "Elast", mtk::Leader_Follower::LEADER );
        tIQI->set_name("IQI_1");

        // define set info
        fem::Set_User_Info tSetBulk1;
        tSetBulk1.set_mesh_set_name( "Omega_0_tets" );
        tSetBulk1.set_IWGs( { tIWG } );
        tSetBulk1.set_IQIs( { tIQI } );
        tSetBulk1.set_is_analytical_sensitivity_analysis( false );
        tSetBulk1.set_finite_difference_scheme_for_sensitivity_analysis(
                fem::FDScheme_Type::POINT_1_FORWARD );
        tSetBulk1.set_finite_difference_perturbation_size( 1e-6 );

        // create a cell of set info
        Vector< fem::Set_User_Info > tSetInfo( 1 );
        tSetInfo( 0 ) = tSetBulk1;

        // FEM model
        //------------------------------------------------------------------------------
        // create model
        mdl::Model * tModel = new mdl::Model(
                tMeshManager,
                0,
                tSetInfo,
                tDesignVariableInterface );

        // get the solver interface
        MSI::MSI_Solver_Interface * tSolverInterface = tModel->get_solver_interface();

        // set requested dof types for solver interface
        Vector< MSI::Dof_Type > tRequestedDofTypes = { MSI::Dof_Type::TEMP };
        tSolverInterface->set_requested_dof_types( tRequestedDofTypes );

        //
        sol::Matrix_Vector_Factory tMatFactory( sol::MapType::Epetra );
        Matrix<DDSMat> tDummyMat; 
        sol::Dist_Map*  mVectorMap = tMatFactory.create_map( {{ 0},{1},{2},{3}}, tDummyMat);
        sol::Dist_Vector * mVector = tMatFactory.create_vector( nullptr, mVectorMap, 1 );

        mVector->sum_into_global_values( Matrix< DDSMat >( {{ 0},{1},{2},{3}} ), Matrix< DDRMat >( {{ 1},{2},{3},{4}} ) );

        // FEM set
        //------------------------------------------------------------------------------
        // get the equation set from the model
        Vector< MSI::Equation_Set * > tSets =
                tModel->get_fem_model()->get_equation_sets();

        // get the equation set from the model
        std::shared_ptr< MSI::Equation_Model > tEquationModel =
                tModel->get_fem_model();

        // get a working set
        MSI::Equation_Set* tWorkSet = tSets( 0 );

        tWorkSet->set_equation_model( tEquationModel.get() );

        // set the dv interface to the set
        tEquationModel->set_design_variable_interface( tDesignVariableInterface );

        // set the IWG and IQI to the set
        reinterpret_cast< fem::Set * >( tWorkSet )->mRequestedIWGs = { tIWG };
        reinterpret_cast< fem::Set * >( tWorkSet )->mRequestedIQIs = { tIQI };

        // set size and fill the set residual assembly map
        reinterpret_cast< fem::Set * >( tWorkSet )->mResDofAssemblyMap.resize( 1 );
        reinterpret_cast< fem::Set * >( tWorkSet )->mResDofAssemblyMap( 0 ) = { { 0, 3 } };

        // set size and fill the set jacobian assembly map
        reinterpret_cast< fem::Set * >( tWorkSet )->mJacDofAssemblyMap.resize( 1 );
        reinterpret_cast< fem::Set * >( tWorkSet )->mJacDofAssemblyMap( 0 ) = { { 0, 3 } };

        // set size and fill the set dv assembly map
        reinterpret_cast< fem::Set * >( tWorkSet )->mPdvMatAssemblyMap.resize( 1 );
        reinterpret_cast< fem::Set * >( tWorkSet )->mPdvMatAssemblyMap( 0 ) = { { 0, 3 } };

        // set size for the set EqnObjDofTypeList
        reinterpret_cast< fem::Set * >( tWorkSet )->mUniqueDofTypeList.resize( 4, MSI::Dof_Type::END_ENUM );

        // set size and populate the set dof type map
        reinterpret_cast< fem::Set * >( tWorkSet )->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
        reinterpret_cast< fem::Set * >( tWorkSet )->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 0;

        // set size and populate the set leader dof type map
        reinterpret_cast< fem::Set * >( tWorkSet )->mLeaderDofTypeMap.set_size( static_cast< int >(MSI::Dof_Type::END_ENUM) + 1, 1, -1 );
        reinterpret_cast< fem::Set * >( tWorkSet )->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 0;

        // set size and populate the set leader dof type map
        reinterpret_cast< fem::Set * >( tWorkSet )->mLeaderDvTypeMap.set_size( static_cast< int >( gen::PDV_Type::UNDEFINED ) + 1, 1, -1 );
        reinterpret_cast< fem::Set * >( tWorkSet )->mLeaderDvTypeMap( static_cast< int >( gen::PDV_Type::DENSITY ) ) = 0;

        // MSI Equation object
        //------------------------------------------------------------------------------
        // get list of equation objects from set
        Vector< MSI::Equation_Object * > tEqObjs = tWorkSet->get_equation_object_list();

        // get a working equation object
        MSI::Equation_Object * tWorkEqObj = tEqObjs( 0 );

        // set the solution vector
        tEquationModel->set_solution_vector( mVector );

        // set the time
        Matrix< DDRMat > tTime = { { 0.0 }, { 1.0 } };
        tEquationModel->set_time( tTime );

        tDesignVariableInterface->set_equation_model( tModel->get_fem_model() );
        tDesignVariableInterface->set_requested_IQIs( {"IQI_1"} );
        tWorkSet->create_requested_IQI_type_map();

        // Init IWG and IQI for forward analysis
        //------------------------------------------------------------------------------
        // set the IWG/IQI fem set
        tIWG->set_set_pointer( reinterpret_cast< fem::Set* >( tWorkSet ) );
        tIQI->set_set_pointer( reinterpret_cast< fem::Set* >( tWorkSet ) );

        // build global dof type list
        tIWG->get_global_dof_type_list();
        tIWG->get_global_dv_type_list();

        // populate the requested leader dof type
        tIWG->mRequestedLeaderGlobalDofTypes = {{ MSI::Dof_Type::TEMP }};

        // compute residual
        tWorkSet->mResidual.resize( 1 );
        tWorkSet->mResidual( 0 ).set_size( 4, 1, 0.0 );
        tWorkSet->mResidualExist = true;
        tWorkEqObj->compute_residual();
        //print( tWorkSet->get_residual()( 0 ), "R" );

        // compute jacobian
        tWorkSet->mJacobian.set_size( 4, 4, 0.0 );
        tWorkSet->mJacobianExist = true;
        tWorkEqObj->compute_jacobian();
        //print( tWorkSet->get_jacobian(), "dRdu" );

        // compute dRdp
        tWorkEqObj->compute_dRdp();
        //print( tWorkSet->get_drdp()( 0 ), "dRdpMat" );
        //print( tWorkSet->get_drdp()( 1 ), "dRdpGeo" );

        // compute QIs
        tWorkSet->mQI.resize( 1 );
        tWorkSet->mQI( 0 ).set_size( 1, 1, 0.0 );
        tWorkSet->mQIExist = true;
        tWorkEqObj->compute_QI();
        //print( tWorkSet->get_QI()( 0 ), "QI" );

        // compute dQIdu
        reinterpret_cast< fem::Set * >( tWorkSet )->mResidual.resize( 1 );
        reinterpret_cast< fem::Set * >( tWorkSet )->mResidual( 0 ).set_size( 4, 1, 0.0 );
        tWorkEqObj->compute_dQIdu();
        //print( tWorkSet->get_residual()( 0 ), "dQIdu" );

//        // compute dQIdp
//        tWorkEqObj->compute_dQIdp_explicit();
//        print( tWorkSet->get_dqidp()( 0 )( 0 ), "dQIdpMat" );
//        print( tWorkSet->get_dqidp()( 1 )( 0 ), "dQIdpGeo" );

//        tEquationModel->compute_explicit_dQIdp();
        delete mVectorMap;

    }

}/* END_TEST_CASE */

}/* end_namespace_msi */
}/* end_namespace_moris */

