/*
 * cl_Equation_Object_Pdv.cpp
 *
 *  Created on: Jul 14, 2018
 *      Author: schmidt
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
#include "cl_MSI_Node_Proxy.hpp"
#include "cl_MSI_Model_Solver_Interface.hpp"
#include "cl_MSI_Dof_Manager.hpp"
#include "cl_MSI_Pdof_Host.hpp"
#include "cl_MTK_Cluster.hpp"
#include "cl_FEM_Cluster.hpp"
#include "cl_FEM_Interpolation_Element.hpp"
#include "cl_FEM_IWG.hpp"         //FEM/INT/src
#include "cl_FEM_Set.hpp"         //FEM/INT/src
#include "cl_FEM_Field_Interpolator_Manager.hpp"                   //FEM//INT//src

#undef protected
#undef private

#include "cl_MTK_Vertex.hpp"    //MTK
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh.hpp"

#include "cl_Mesh_Factory.hpp"
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

#include "cl_MDL_Model.hpp"

#include "MSI_Test_Proxy/cl_MSI_Design_Variable_Interface_Proxy.hpp"
#include "MSI_Test_Proxy/cl_MTK_Vertex_Proxy.hpp"
#include "MSI_Test_Proxy/cl_MTK_Cell_Proxy.hpp"
#include "MSI_Test_Proxy/cl_MTK_Cluster_Proxy.hpp"
#include "cl_MSI_Solver_Interface.hpp"

#include "cl_FEM_CM_Factory.hpp"
#include "cl_FEM_IWG_Factory.hpp"
#include "cl_FEM_Field_Interpolator.hpp"                   //FEM//INT//src
#include "cl_FEM_Geometry_Interpolator.hpp"                   //FEM//INT//src
#include "cl_FEM_Property.hpp"                   //FEM//INT//src

#include "cl_Matrix_Vector_Factory.hpp"
#include "cl_Map_Class.hpp"
#include "cl_SOL_Dist_Vector.hpp"


namespace moris
{
    namespace MSI
    {
    moris::Matrix< moris::DDRMat > tConstValFunction_FDTest
    ( moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
      moris::fem::Field_Interpolator_Manager *         aFIManager )
    {
        return aParameters( 0 );
    }

    moris::Matrix< moris::DDRMat > tFIValDvFunction_FDTest
    ( moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
      moris::fem::Field_Interpolator_Manager *         aFIManager )
    {
        return aParameters( 0 ) * aFIManager->get_field_interpolators_for_type( moris::MSI::Dv_Type::DENSITY )->val();
    }

    moris::Matrix< moris::DDRMat > tFIDerDvFunction_FDTest
   ( moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
     moris::fem::Field_Interpolator_Manager *         aFIManager )
    {
        return aParameters( 0 ) * aFIManager->get_field_interpolators_for_type( moris::MSI::Dv_Type::DENSITY )->N();
    }

    TEST_CASE("Eqn_Obj_pdv","[MSI],[Eqn_Obj_pdv]")
    {
        if(par_size() == 1 )
        {
    	// setup the interpolation mesh
//        std::string tInterpString = "generated:1x1";
//        moris::mtk::Interpolation_Mesh* tInterpMesh1 = moris::mtk::create_interpolation_mesh( MeshType::STK, tInterpString );

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
        Matrix<DDRMat>   tNodeCoordinates ={{0, 0},{1, 0},{1, 1},{0, 1},{0.5, 0},{0, 0.5}};
        Matrix<IndexMat> tLocalToGlobalNodeMap = {{1, 2, 3, 4, 5, 6 }};

        // Tetrathedral cells in material phase 1
        CellTopology     tPhase0ChildTopo  = CellTopology::TRI3;
        Matrix<IndexMat> tCellIdsPhase0    = {{1,2,3,4}};
        Matrix<IndexMat> tCellToNodePhase0 = {{1, 5, 6},{2,3,5},{5,3,6},{6,3,4}};

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

        moris::mtk::Integration_Mesh* tIntegMesh1  = moris::mtk::create_integration_mesh(MeshType::STK,tMeshDataInput,tInterpMesh1);

        // place the pair in mesh manager
        mtk::Mesh_Manager tMeshManager;
        tMeshManager.register_mesh_pair(tInterpMesh1,tIntegMesh1);

        Design_Variable_Interface * tDesignVariableInterface = new Design_Variable_Interface_Proxy();

        /*
        MSI::Equation_Set * tSet = new fem::Set();

        tSet->set_Dv_interface( tDesignVariableInterface );

        // Create generic equation objects
        MSI::Equation_Object * EquObj = new fem::Interpolation_Element();

        EquObj->mEquationSet = tSet;
        reinterpret_cast< fem::Interpolation_Element *> ( EquObj )->mSet = reinterpret_cast< fem::Set *> (tSet);

        tSet->mMasterDofTypes = { { MSI::Dof_Type::TEMP } };


        std::shared_ptr< fem::Cluster > tFemCluster = std::make_shared< fem::Cluster >( fem::Element_Type::BULK,
                                                                                        tCluster,
                                                                                        reinterpret_cast< fem::Set *> (tSet),
                                                                                        EquObj );

        reinterpret_cast< fem::Interpolation_Element * >( EquObj )->mFemCluster = {tFemCluster};
        */


        // create the properties
        std::shared_ptr< fem::Property > tPropMasterConductivity = std::make_shared< fem::Property > ();
        tPropMasterConductivity->set_parameters( { {{ 1.0 }} } );
        tPropMasterConductivity->set_dv_type_list( {{ MSI::Dv_Type::DENSITY }} );
        tPropMasterConductivity->set_val_function( tFIValDvFunction_FDTest );
        tPropMasterConductivity->set_dv_derivative_functions( { tFIDerDvFunction_FDTest } );

        std::shared_ptr< fem::Property > tPropMasterTempLoad = std::make_shared< fem::Property > ();
        tPropMasterTempLoad->set_parameters( { {{ 1.0 }} } );
        tPropMasterTempLoad->set_val_function( tConstValFunction_FDTest );

        // define constitutive models
        fem::CM_Factory tCMFactory;

        std::shared_ptr< fem::Constitutive_Model > tCMMasterDiffLinIso = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
        tCMMasterDiffLinIso->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tCMMasterDiffLinIso->set_property( tPropMasterConductivity, "Conductivity" );
        tCMMasterDiffLinIso->set_space_dim( 2 );

        // define the IWGs
        fem::IWG_Factory tIWGFactory;

        std::shared_ptr< fem::IWG > tIWG = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
        tIWG->set_residual_dof_type( { MSI::Dof_Type::TEMP } );
        tIWG->set_dof_type_list( {{ MSI::Dof_Type::TEMP }}, mtk::Master_Slave::MASTER );
        tIWG->set_constitutive_model( tCMMasterDiffLinIso, "DiffLinIso", mtk::Master_Slave::MASTER );
        tIWG->set_property( tPropMasterTempLoad, "Load", mtk::Master_Slave::MASTER );

        // define set info
        fem::Set_User_Info tSetBulk1;
        tSetBulk1.set_mesh_index( tIntegMesh1->get_set_index_by_name("Omega_0_tets") );
        tSetBulk1.set_IWGs( { tIWG } );

        // create a cell of set info
        moris::Cell< fem::Set_User_Info > tSetInfo( 1 );
        tSetInfo( 0 ) = tSetBulk1;

        // create model
        mdl::Model * tModel = new mdl::Model( &tMeshManager,
                                               0,
                                               tSetInfo,
                                               tDesignVariableInterface);

        MSI::MSI_Solver_Interface * tSolverInterface = tModel->get_solver_interface();

        Matrix_Vector_Factory tMatFactory( MapType::Epetra );
        Map_Class * mVectorMap = tMatFactory.create_map( {{ 0},{1},{2},{3}} );
        Dist_Vector * mVector = tMatFactory.create_vector( nullptr, mVectorMap, VectorType::FREE );

        mVector->sum_into_global_values( 4, {{ 0},{1},{2},{3}}, {{ 1},{2},{3},{4}});

//        tSolverInterface->set_solution_vector( mVector );

        moris::Cell< MSI::Equation_Set * > tSets =  tModel->get_equation_sets( );

//        tSets( 0 )->set_Dv_interface( tDesignVariableInterface );

        reinterpret_cast< fem::Set * >(tSets( 0 ))->mRequestedIWGs = { tIWG };

        Cell< MSI::Equation_Object * > tEquationObject = tSets( 0 )->get_equation_object_list();

//        tEquationObject( 0 )->mPdofValues = {{ 0},{0},{0},{0}};

        tEquationObject( 0 )->mSolVec = mVector;

        tEquationObject( 0 )->set_time({{ 0},{1}});

        tIWG->set_set_pointer( reinterpret_cast< fem::Set* >( tSets( 0 ) ) );

//            // set size for the set EqnObjDofTypeList
//            tIWG->mSet->mEqnObjDofTypeList.resize( 4, MSI::Dof_Type::END_ENUM );
//
//            // set size and populate the set dof type map
//            tIWG->mSet->mDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
//            tIWG->mSet->mDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 0;
//
//            // set size and populate the set master dof type map
//            tIWG->mSet->mMasterDofTypeMap.set_size( static_cast< int >(MSI::Dof_Type::END_ENUM) + 1, 1, -1 );
//            tIWG->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 0;

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( 1 );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, 3 } };

            // set size and fill the set jacobian assembly map
            tIWG->mSet->mJacDofAssemblyMap.resize( 1 );
            tIWG->mSet->mJacDofAssemblyMap( 0 ) = { { 0, 3 } };

            // set size and fill the set jacobian assembly map
            tIWG->mSet->mDvAssemblyMap.resize( 1 );
            tIWG->mSet->mDvAssemblyMap( 0 ) = { { 0, 3 } };

            // set size and init the set residual and jacobian
            tIWG->mSet->mResidual.resize( 1 );
            tIWG->mSet->mResidual( 0 ).set_size( 4, 1, 0.0 );
            tIWG->mSet->mJacobian.set_size( 4, 4, 0.0 );

            // set requested residual dof type flag to true
            tIWG->mResidualDofTypeRequested = true;

            // build global dof type list
            tIWG->get_global_dof_type_list();
            tIWG->get_global_dv_type_list();

            // populate the requested master dof type
            tIWG->mRequestedMasterGlobalDofTypes = {{ MSI::Dof_Type::TEMP }};


        tEquationObject( 0 )->compute_dRdp();




//        // Create the pdof hosts of this equation object
//        moris::Cell < Pdof_Host * > tPdofHostList;
//        tPdofHostList.resize( 3, nullptr );
//        moris::uint tNumMaxPdofTypes = 1;
//
//        Matrix< DDSMat > tDofTypeIndexMap(4, 1, -1);
//        tDofTypeIndexMap(3, 0) = 0;
//
//        Matrix< DDUMat > tTimePerDofType(4, 1, 1);
//
//        Equation_Set tEqnBlock;
//        tEqnBlock.mEqnObjDofTypeList.resize( 1, MSI::Dof_Type::TEMP );
//        EquObj.mEquationSet = &tEqnBlock;
//
//        EquObj.create_my_pdof_hosts( tNumMaxPdofTypes, tDofTypeIndexMap, tTimePerDofType, tPdofHostList );
//
//        // Check if right pdof host was created in given pdof host list
//        CHECK( equal_to( tPdofHostList( 0 )->mNodeID, 0 ) );
//        REQUIRE( tPdofHostList( 1 ) == NULL );
//        CHECK( equal_to( tPdofHostList( 2 )->mNodeID, 2 ) );
//
//        // Check equation objects internal pdof host list
//        CHECK( equal_to( EquObj.mMyPdofHosts( 0 ).size(), 2 ) );
//        CHECK( equal_to( EquObj.mMyPdofHosts( 0 )( 0 )->mNodeID, 0 ) );
//        CHECK( equal_to( EquObj.mMyPdofHosts( 0 )( 1 )->mNodeID, 2 ) );
//        delete Node1;
//        delete Node2;
//        delete tPdofHostList(0);
//        delete tPdofHostList(1);
//        delete tPdofHostList(2);
    }
    }

    }
}


