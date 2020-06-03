/*
 * cl_Fiber_Problem.cpp
 *
 *  Created on: Oct 22, 2019
 *      Author: sonne
 */
#include "catch.hpp"

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
#include "cl_Geom_Field.hpp"
#include "typedefs.hpp"

#include "cl_MTK_Mesh_Manager.hpp"

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

#include "cl_HMR_Mesh_Interpolation.hpp"
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
#include "cl_MSI_Solver_Interface.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_MSI_Model_Solver_Interface.hpp"
#include "cl_DLA_Linear_Solver_Aztec.hpp"
#include "cl_DLA_Linear_Solver.hpp"

#include "cl_TSA_Time_Solver_Factory.hpp"
#include "cl_TSA_Monolithic_Time_Solver.hpp"
#include "cl_TSA_Time_Solver.hpp"

#include "fn_norm.hpp"

#include "cl_GEN_Circle.hpp"
#include "cl_GEN_Enums.hpp"
#include "cl_GEN_Geom_Field_HMR.hpp"
#include "cl_GEN_Geometry.hpp"
#include "cl_GEN_Geom_Data.hpp"
#include "cl_GEN_Property.hpp"


#include "cl_PRM_HMR_Parameters.hpp"


using namespace moris;
namespace ge
{

    static const uint sNumberOfRefinements = 4;
                 real gLsbwabs;

void tConstValFunction
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 );
}

void tMValFunction
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = {{ aParameters( 0 )( 0 ),                      0.0,                     0.0 },
        {                   0.0,    aParameters( 0 )( 1 ),                     0.0 },
        {                   0.0,                      0.0,    aParameters( 0 )( 2 )}};
}

moris::real LvlSetCircle_2D_outsideDomain(const moris::Matrix< moris::DDRMat > & aPoint )
{
    Matrix< DDRMat > tCenter{ { -10, -10 } };
    return    norm( aPoint - tCenter ) - 0.001;
}

int user_defined_refinement(       hmr::Element             * aElement,
        const Cell< Matrix< DDRMat > > & aElementLocalValues,
        moris::ParameterList       & aParameters )
{
    int aDoRefine = -1;

    // abs variable field threshold
    real lsth = 0.0;

    // abs variable field bandwidth (absolute)
    real lsbwabs = gLsbwabs;

    // maximum refinement level
    uint maxlevel = sNumberOfRefinements;

    // min refinement level
    uint minlevel = 0;

    // max refinement level along interface
    uint maxifcref = sNumberOfRefinements;

    // max refinement level in volume
    uint maxvolref = 0;

    // current refinement level of element
    uint curlevel = aElement->get_level();

    // refinement strategy
    if ( aElementLocalValues( 0 ).min() <= lsth + lsbwabs )
    {
        // for volume refinement
        if ( aElementLocalValues( 0 ).min() <= lsth - lsbwabs )
        {
            if( curlevel < maxvolref && curlevel < maxlevel )
            {
                aDoRefine = 1; // refine
            }
            else if ( curlevel ==  maxvolref || curlevel == minlevel )
            {
                aDoRefine = 0; // keep
            }
            else
            {
                aDoRefine = -1; // coarsen
            }
        }
        // for interface refinement
        else
        {
            if( curlevel < maxifcref && curlevel < maxlevel )
            {
                aDoRefine = 1; // refine
            }
            else if ( curlevel ==  maxifcref || curlevel == minlevel )
            {
                aDoRefine = 0; // keep
            }
            else
            {
                aDoRefine = -1; // coarsen
            }
        }
    }
    else
    {
        if( curlevel <  minlevel )
        {
            aDoRefine = 1; // refine
        }
        else if ( curlevel == minlevel )
        {
            aDoRefine = 0; // keep
        }
    }

    return aDoRefine;
}

//------------------------------------------------------------------------------

TEST_CASE("fiber_problem_test", "[GE],[fiber_test]")
{
    uint tLagrangeMeshIndex = 0;
    //  HMR Parameters setup
    moris::ParameterList tParameters = prm::create_hmr_parameter_list();

    uint tInitialMesh = 1;
    switch(tInitialMesh)
    {
    case(1) :
        {
            tParameters.set( "number_of_elements_per_dimension", std::string("80, 40, 10") );
            break;
        }
    default :
        {
            tParameters.set( "number_of_elements_per_dimension", std::string("40, 20, 5") );
        }
    }

    tParameters.set( "domain_dimensions",                std::string("40, 20, 5") );
    tParameters.set( "domain_offset",                    std::string("-0, -0, -0") );

    tParameters.set( "domain_sidesets", std::string("1, 2, 3, 4, 5, 6") );

    tParameters.set( "truncate_bsplines", 1 );
    tParameters.set( "lagrange_orders", std::string("1") );
    tParameters.set( "lagrange_pattern", std::string("0") );
    tParameters.set( "bspline_orders", std::string("1") );
    tParameters.set( "bspline_pattern", std::string("0") );

    tParameters.set( "lagrange_output_meshes", std::string("0") );
    tParameters.set( "lagrange_input_meshes", std::string("0") );

    tParameters.set( "lagrange_to_bspline", std::string("0") );

    tParameters.set( "use_multigrid", 0 );

    tParameters.set( "refinement_buffer", 1 );
    tParameters.set( "staircase_buffer", 1 );

    tParameters.insert( "initial_refinement", 0 );

    gLsbwabs=0.3;

    //  HMR Initialization
    moris::hmr::HMR tHMR( tParameters );

    auto tDatabase = tHMR.get_database(); // std::shared_ptr< Database >

    tHMR.perform_initial_refinement( 0 );
    //------------------------------------------------------------------------------

    tDatabase->update_bspline_meshes();
    tDatabase->update_lagrange_meshes();

    uint tNumberOfFibers      = 1;
    uint tNumberOfRefinements = 1;
    std::cout<<"-------------------------------------------"<<std::endl;
    std::cout<<"number of fibers being analyzed:      "<<tNumberOfFibers<<std::endl;
    std::cout<<"number of refinements performed:      "<<tNumberOfRefinements<<std::endl;
    std::cout<<"bandwidth in user-defined refinement: "<<gLsbwabs<<std::endl;
    std::cout<<"-------------------------------------------"<<std::endl;
    Cell< Matrix< DDRMat > > tFieldData( 1 );

    std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

    moris::tic tTimer0;
    for( uint k=0; k<tNumberOfRefinements; ++k )
    {
        moris::ge::GEN_CylinderWithEndCaps    tFibers( tNumberOfFibers );
        moris::Cell<moris::ge::GEN_Geometry*> tGeometryVector = { &tFibers };

        moris::ge::GEN_Phase_Table         tPhaseTable( tGeometryVector.size(),  Phase_Table_Structure::EXP_BASE_2 );
        moris::ge::Geometry_Engine     tGENGeometryEngine( tGeometryVector,tPhaseTable,3 );

        moris_index tMeshIndex = tGENGeometryEngine.register_mesh( tMesh );

        tFieldData( 0 ) = tGENGeometryEngine.get_cylinder_vals( tMeshIndex, &tFibers, tNumberOfFibers );

        tHMR.user_defined_flagging( user_defined_refinement, tFieldData, tParameters, 0 );

        tHMR.perform_refinement_based_on_working_pattern( 0, false );

        //gLsbwabs *= 0.5;
    }

    tHMR.finalize();

    moris::ge::GEN_CylinderWithEndCaps    tFibers( tNumberOfFibers );
    moris::Cell<moris::ge::GEN_Geometry*> tGeometryVector_temp = {&tFibers};

    moris::ge::GEN_Phase_Table         tPhaseTable_temp( tGeometryVector_temp.size(),  Phase_Table_Structure::EXP_BASE_2 );
    moris::ge::Geometry_Engine     tGENGeometryEngine_temp( tGeometryVector_temp, tPhaseTable_temp, 3 );

    moris_index tMeshIndex = tGENGeometryEngine_temp.register_mesh( tMesh );

    tFieldData( 0 ) = tGENGeometryEngine_temp.get_cylinder_vals( tMeshIndex, &tFibers, tNumberOfFibers );

    hmr::Interpolation_Mesh_HMR *      tInterpMesh      = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );
    moris::hmr::Integration_Mesh_HMR * tIntegrationMesh = tHMR.create_integration_mesh( 1, 0, tInterpMesh );

    mtk::Mesh_Manager tMesh1;

    real tElapsedTime0 = tTimer0.toc<moris::chronos::milliseconds>().wall;
    tElapsedTime0 /= 1000;
    std::cout<<"==============================================="<< std::endl;
    std::cout<<"Total time for evaluation: "<< tElapsedTime0 << std::endl << std::endl;
    std::cout<<"==============================================="<< std::endl;
    //------------------------------------------------------------------------------

    moris::ge::GEN_Geom_Data tGeomData( tFieldData( 0 ) );  // fiber LS values as field data
    //==============================
    real tRadius  = 10.0001;
    real tXcenter = 0.0;//-1000.0;
    real tYcenter = 0.0;//-1000.0;
    moris::ge::Circle tHoleGeom( tRadius, tXcenter, tYcenter );

    moris::Cell<moris::ge::GEN_Geometry*> tGeometryVector = { &tGeomData, &tHoleGeom };

    size_t tModelDimension = 3;
    moris::ge::GEN_Phase_Table      tPhaseTable( tGeometryVector.size(),  Phase_Table_Structure::EXP_BASE_2 );
    moris::ge::Geometry_Engine  tGENGeometryEngine( tGeometryVector, tPhaseTable, tModelDimension );

    //------------------------------------------------------------------------------
    xtk::Model                      tXTKModel( tModelDimension, tInterpMesh, &tGENGeometryEngine );
    tXTKModel.mVerbose = false;

    Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4};
    tXTKModel.decompose(tDecompositionMethods);
    //=============================== temporary ============================================
    // output solution and meshes
    bool tOutputXTKmesh = true;
    if (tOutputXTKmesh)
    {
        xtk::Output_Options tOutputOptions1;
        tOutputOptions1.mAddNodeSets = false;
        tOutputOptions1.mAddSideSets = true;
        tOutputOptions1.mAddClusters = false;

        std::string tIntegSolFieldNameUX = "UX";
        std::string tIntegSolFieldNameUY = "UY";
        std::string tIntegSolFieldNameUZ = "UZ";
        tOutputOptions1.mRealNodeExternalFieldNames = { tIntegSolFieldNameUX, tIntegSolFieldNameUY, tIntegSolFieldNameUZ };

        moris::mtk::Integration_Mesh* tIntegMesh11 = tXTKModel.get_output_mesh(tOutputOptions1);

        std::string tMeshOutputFile1 = "./fibersGeomCheckDebugging.e";
        tIntegMesh11->create_output_mesh(tMeshOutputFile1);
        delete tIntegMesh11;
    }
    //============================= end temporary ==========================================
    bool tFullProblem = false;
    if( tFullProblem )
    {
        tXTKModel.perform_basis_enrichment(EntityRank::BSPLINE,0);

        xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
        xtk::Enriched_Integration_Mesh   & tEnrIntegMesh = tXTKModel.get_enriched_integ_mesh();

        // place the pair in mesh manager
        mtk::Mesh_Manager tMeshManager;
        tMeshManager.register_mesh_pair( &tEnrInterpMesh, &tEnrIntegMesh);
        //------------------------------------------------------------------------------
        // create the material properties
        real tEModPlate  =   1.22e6;   	//    1.22*10^6 Pa
        real tEModFibers = 103.023e6; 	// 1030.23*10^6 Pa

        std::cout<<"----------------------------------"<<std::endl;
        std::cout<<"E plate:   "<<tEModPlate<<std::endl;
        std::cout<<"E fibers:  "<<tEModFibers<<std::endl;
        std::cout<<"----------------------------------"<<std::endl;

        std::shared_ptr< fem::Property > tPropEModPlate = std::make_shared< fem::Property >();
        tPropEModPlate->set_parameters( { {{ tEModPlate }} } );   // 1.22*10^6 Pa
        //    tPropEModPlate->set_parameters( { {{ 100 }} } );
        tPropEModPlate->set_val_function( tConstValFunction );

        std::shared_ptr< fem::Property > tPropNuPlate = std::make_shared< fem::Property >();
        tPropNuPlate->set_parameters( { {{ 0.2 }} } );   // Poinson's ratio = 0.4
        tPropNuPlate->set_val_function( tConstValFunction );

        std::shared_ptr< fem::Property > tPropEModFibers = std::make_shared< fem::Property >();
        tPropEModFibers->set_parameters( { {{ tEModFibers }} } );  // 1030.23*10^6 Pa
        //    tPropEModFibers->set_parameters( { {{ 100 }} } );
        tPropEModFibers->set_val_function( tConstValFunction );

        std::shared_ptr< fem::Property > tPropNuFibers = std::make_shared< fem::Property >();
        tPropNuFibers->set_parameters( { {{ 0.2 }} } );          // Poinson's ratio = 0.4
        tPropNuFibers->set_val_function( tConstValFunction );

        //======================================================
        // create the boundary conditions and loading properties
        //======================================================
        // symmetry boundary conditions
        std::shared_ptr< fem::Property > tPropDirichletUX_ss4 = std::make_shared< fem::Property >();    // fix displacement to be zero in x
        tPropDirichletUX_ss4->set_parameters( { {{ 0.0 }, {1.0}, {1.0}} } );
        //    tPropDirichletUX_ss4->set_parameters( { {{ 0.0 }, {0.0}, {0.0}} } );
        tPropDirichletUX_ss4->set_val_function( tConstValFunction );

        std::shared_ptr< fem::Property > tPropDirichletUX_ss4_select = std::make_shared< fem::Property >();    // fix displacement to be zero in x
        tPropDirichletUX_ss4_select->set_parameters( { {{ 1.0 }, {0.0}, {0.0}} } );
        //    tPropDirichletUX_ss4_select->set_parameters( { {{ 1.0 }, {1.0}, {1.0}} } );
        tPropDirichletUX_ss4_select->set_val_function( tMValFunction );

        //------------------------------------------------------------------------------
        std::shared_ptr< fem::Property > tPropDirichletUY_ss1 = std::make_shared< fem::Property >();    // fix displacement to be zero in y
        tPropDirichletUY_ss1->set_parameters( { {{ 1.0 }, {0.0}, {1.0}} } );
        tPropDirichletUY_ss1->set_val_function( tConstValFunction );

        std::shared_ptr< fem::Property > tPropDirichletUY_ss1_select = std::make_shared< fem::Property >();    // fix displacement to be zero in y
        tPropDirichletUY_ss1_select->set_parameters( { {{ 0.0 }, {1.0}, {0.0}} } );
        tPropDirichletUY_ss1_select->set_val_function( tMValFunction );

        //------------------------------------------------------------------------------
        std::shared_ptr< fem::Property > tPropDirichletUZ_ss5 = std::make_shared< fem::Property >();    // fix displacement to be zero in z
        tPropDirichletUZ_ss5->set_parameters( { {{ 1.0 }, {1.0}, {0.0}} } );
        tPropDirichletUZ_ss5->set_val_function( tConstValFunction );

        std::shared_ptr< fem::Property > tPropDirichletUZ_ss5_select = std::make_shared< fem::Property >();    // fix displacement to be zero in y
        tPropDirichletUZ_ss5_select->set_parameters( { {{ 0.0 }, {0.0}, {1.0}} } );
        tPropDirichletUZ_ss5_select->set_val_function( tMValFunction );

        //==============================
        // loading on free end
        std::shared_ptr< fem::Property > tPropNeumann = std::make_shared< fem::Property >();
        tPropNeumann->set_parameters( {{{ 1000000 } , { 0.0 }, {0.0}}} );  // 1.0 MPa traction applied in the x-direction
        tPropNeumann->set_val_function( tConstValFunction );

        //==============================
        // define constitutive models
        fem::CM_Factory tCMFactory;

        std::shared_ptr< fem::Constitutive_Model > tCMPlate = tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
        tCMPlate->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY , MSI::Dof_Type::UZ }} );
        tCMPlate->set_property( tPropEModPlate, "YoungsModulus" );
        tCMPlate->set_property( tPropNuPlate, "PoissonRatio" );
        tCMPlate->set_space_dim( 3 );

        //------------------------------------------------------------------------------
        std::shared_ptr< fem::Constitutive_Model > tCMFibers = tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
        tCMFibers->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY , MSI::Dof_Type::UZ }} );
        tCMFibers->set_property( tPropEModFibers, "YoungsModulus" );
        tCMFibers->set_property( tPropNuFibers, "PoissonRatio" );
        tCMFibers->set_space_dim( 3 );

        //==============================
        // define stabilization parameters
        fem::SP_Factory tSPFactory;

        //------------------------------------------------------------------------------
        std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitscheBCs = tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
        tSPDirichletNitscheBCs->set_parameters( { {{ 100.0 }} } );
        tSPDirichletNitscheBCs->set_property( tPropEModPlate, "Material", mtk::Master_Slave::MASTER );

        std::shared_ptr< fem::Stabilization_Parameter > tSPNitscheInterfaceFibersToPlate = tSPFactory.create_SP( fem::Stabilization_Type::NITSCHE_INTERFACE );
        tSPNitscheInterfaceFibersToPlate->set_parameters( { {{ 100.0 }} } );
        tSPNitscheInterfaceFibersToPlate->set_property( tPropEModFibers, "Material", mtk::Master_Slave::MASTER );
        tSPNitscheInterfaceFibersToPlate->set_property( tPropEModPlate, "Material", mtk::Master_Slave::SLAVE );

        std::shared_ptr< fem::Stabilization_Parameter > tSPMasterWeightInterfaceFibersToPlate = tSPFactory.create_SP( fem::Stabilization_Type::MASTER_WEIGHT_INTERFACE );
        tSPMasterWeightInterfaceFibersToPlate->set_property( tPropEModFibers, "Material", mtk::Master_Slave::MASTER );
        tSPMasterWeightInterfaceFibersToPlate->set_property( tPropEModPlate, "Material", mtk::Master_Slave::SLAVE );

        std::shared_ptr< fem::Stabilization_Parameter > tSPSlaveWeightInterfaceFibersToPlate = tSPFactory.create_SP( fem::Stabilization_Type::SLAVE_WEIGHT_INTERFACE );
        tSPSlaveWeightInterfaceFibersToPlate->set_property( tPropEModFibers, "Material", mtk::Master_Slave::MASTER );
        tSPSlaveWeightInterfaceFibersToPlate->set_property( tPropEModPlate, "Material", mtk::Master_Slave::SLAVE );


        //==============================
        // define the IWGs
        fem::IWG_Factory tIWGFactory;
        //------------------------------------------------------------------------------
        std::shared_ptr< fem::IWG > tIWGPlate = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_BULK );
        tIWGPlate->set_residual_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ } );
        tIWGPlate->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }} );
        tIWGPlate->set_constitutive_model( tCMPlate, "ElastLinIso", mtk::Master_Slave::MASTER );

        //------------------------------------------------------------------------------
        std::shared_ptr< fem::IWG > tIWGDirichletFixedUx = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE );
        tIWGDirichletFixedUx->set_residual_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ } );
        tIWGDirichletFixedUx->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }} );
        tIWGDirichletFixedUx->set_stabilization_parameter( tSPDirichletNitscheBCs, "DirichletNitsche" );
        tIWGDirichletFixedUx->set_constitutive_model( tCMPlate, "ElastLinIso", mtk::Master_Slave::MASTER );
        tIWGDirichletFixedUx->set_property( tPropDirichletUX_ss4, "Dirichlet", mtk::Master_Slave::MASTER );
        tIWGDirichletFixedUx->set_property( tPropDirichletUX_ss4_select, "Select", mtk::Master_Slave::MASTER );

        //------------------------------------------------------------------------------
        std::shared_ptr< fem::IWG > tIWGNeumann = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_NEUMANN );
        tIWGNeumann->set_residual_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ } );
        tIWGNeumann->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }} );
        tIWGNeumann->set_property( tPropNeumann, "Neumann", mtk::Master_Slave::MASTER );

        //------------------------------------------------------------------------------
        std::shared_ptr< fem::IWG > tIWGFibers = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_BULK );
        tIWGFibers->set_residual_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ } );
        tIWGFibers->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }} );
        tIWGFibers->set_constitutive_model( tCMFibers, "ElastLinIso", mtk::Master_Slave::MASTER );

        //------------------------------------------------------------------------------
        std::shared_ptr< fem::IWG > tIWGDirichletFixedUy = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE );
        tIWGDirichletFixedUy->set_residual_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ } );
        tIWGDirichletFixedUy->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }} );
        tIWGDirichletFixedUy->set_stabilization_parameter( tSPDirichletNitscheBCs, "DirichletNitsche" );
        tIWGDirichletFixedUy->set_constitutive_model( tCMPlate, "ElastLinIso", mtk::Master_Slave::MASTER );
        tIWGDirichletFixedUy->set_property( tPropDirichletUY_ss1, "Dirichlet", mtk::Master_Slave::MASTER );
        tIWGDirichletFixedUy->set_property( tPropDirichletUY_ss1_select, "Select", mtk::Master_Slave::MASTER );

        //------------------------------------------------------------------------------
        std::shared_ptr< fem::IWG > tIWGDirichletFixedUz = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE );
        tIWGDirichletFixedUz->set_residual_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ } );
        tIWGDirichletFixedUz->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }} );
        tIWGDirichletFixedUz->set_stabilization_parameter( tSPDirichletNitscheBCs, "DirichletNitsche" );
        tIWGDirichletFixedUz->set_constitutive_model( tCMPlate, "ElastLinIso", mtk::Master_Slave::MASTER );
        tIWGDirichletFixedUz->set_property( tPropDirichletUZ_ss5, "Dirichlet", mtk::Master_Slave::MASTER );
        tIWGDirichletFixedUz->set_property( tPropDirichletUZ_ss5_select, "Select", mtk::Master_Slave::MASTER );

        //------------------------------------------------------------------------------
        std::shared_ptr< fem::IWG > tIWGFiberInterfacePlateBulk = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_INTERFACE );
        tIWGFiberInterfacePlateBulk->set_residual_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ } );
        tIWGFiberInterfacePlateBulk->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }}, mtk::Master_Slave::MASTER );
        tIWGFiberInterfacePlateBulk->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }}, mtk::Master_Slave::SLAVE );
        tIWGFiberInterfacePlateBulk->set_stabilization_parameter( tSPNitscheInterfaceFibersToPlate, "NitscheInterface" );
        tIWGFiberInterfacePlateBulk->set_stabilization_parameter( tSPMasterWeightInterfaceFibersToPlate, "MasterWeightInterface" );
        tIWGFiberInterfacePlateBulk->set_stabilization_parameter( tSPSlaveWeightInterfaceFibersToPlate, "SlaveWeightInterface" );
        tIWGFiberInterfacePlateBulk->set_constitutive_model( tCMFibers, "ElastLinIso", mtk::Master_Slave::MASTER );
        tIWGFiberInterfacePlateBulk->set_constitutive_model( tCMPlate , "ElastLinIso", mtk::Master_Slave::SLAVE  );

        //==============================
        // create a list of active block-sets
        //    std::string tInterfaceSideSetName    = tEnrIntegMesh.get_interface_side_set_name( 0, 0, 0 );
        //    std::string tDblInterfaceSideSetName = tEnrIntegMesh.get_dbl_interface_side_set_name(0,1);
        //==============================
        // define set info
        // fibers - phase 1
        // hole   - phase 2
        // plate  - phase 3

        //==============================
        // bulk for plate
        fem::Set_User_Info tBulkPlate00;
        tBulkPlate00.set_mesh_index( tEnrIntegMesh.get_set_index_by_name("HMR_dummy_n_p3") );
        tBulkPlate00.set_IWGs( { tIWGPlate } );

        fem::Set_User_Info tBulkPlate01;
        tBulkPlate01.set_mesh_index( tEnrIntegMesh.get_set_index_by_name("HMR_dummy_c_p3") );
        tBulkPlate01.set_IWGs( { tIWGPlate } );

        //==============================
        // bulk for fibers
        fem::Set_User_Info tBulkFibers00;
        tBulkFibers00.set_mesh_index( tEnrIntegMesh.get_set_index_by_name("HMR_dummy_n_p1") );
        tBulkFibers00.set_IWGs( { tIWGFibers } );

        fem::Set_User_Info tBulkFibers01;
        tBulkFibers01.set_mesh_index( tEnrIntegMesh.get_set_index_by_name("HMR_dummy_c_p1") );
        tBulkFibers01.set_IWGs( { tIWGFibers } );

        //==============================
        // symmetry boundary conditions on side-set 4 ( fix Ux = 0 )
        fem::Set_User_Info tSetDirichletFixed00;
        tSetDirichletFixed00.set_mesh_index( tEnrIntegMesh.get_set_index_by_name("SideSet_4_n_p1") );
        tSetDirichletFixed00.set_IWGs( { tIWGDirichletFixedUx } );

        fem::Set_User_Info tSetDirichletFixed01;
        tSetDirichletFixed01.set_mesh_index( tEnrIntegMesh.get_set_index_by_name("SideSet_4_c_p1") );
        tSetDirichletFixed01.set_IWGs( { tIWGDirichletFixedUx } );

        fem::Set_User_Info tSetDirichletFixed02;
        tSetDirichletFixed02.set_mesh_index( tEnrIntegMesh.get_set_index_by_name("SideSet_4_n_p3") );
        tSetDirichletFixed02.set_IWGs( { tIWGDirichletFixedUx } );

        fem::Set_User_Info tSetDirichletFixed03;
        tSetDirichletFixed03.set_mesh_index( tEnrIntegMesh.get_set_index_by_name("SideSet_4_c_p3") );
        tSetDirichletFixed03.set_IWGs( { tIWGDirichletFixedUx } );

        //==============================
        // symmetry boundary conditions on side-set 1 ( fix Uy = 0 )
        fem::Set_User_Info tSetDirichletFixed04;
        tSetDirichletFixed04.set_mesh_index( tEnrIntegMesh.get_set_index_by_name("SideSet_1_n_p1") );
        tSetDirichletFixed04.set_IWGs( { tIWGDirichletFixedUy } );

        fem::Set_User_Info tSetDirichletFixed05;
        tSetDirichletFixed05.set_mesh_index( tEnrIntegMesh.get_set_index_by_name("SideSet_1_c_p1") );
        tSetDirichletFixed05.set_IWGs( { tIWGDirichletFixedUy } );

        fem::Set_User_Info tSetDirichletFixed06;
        tSetDirichletFixed06.set_mesh_index( tEnrIntegMesh.get_set_index_by_name("SideSet_1_n_p3") );
        tSetDirichletFixed06.set_IWGs( { tIWGDirichletFixedUy } );

        fem::Set_User_Info tSetDirichletFixed07;
        tSetDirichletFixed07.set_mesh_index( tEnrIntegMesh.get_set_index_by_name("SideSet_1_c_p3") );
        tSetDirichletFixed07.set_IWGs( { tIWGDirichletFixedUy } );

        //==============================
        // symmetry boundary conditions on side-set 5 ( fix Uz = 0 )
        fem::Set_User_Info tSetDirichletFixed08;
        tSetDirichletFixed08.set_mesh_index( tEnrIntegMesh.get_set_index_by_name("SideSet_5_n_p1") );
        tSetDirichletFixed08.set_IWGs( { tIWGDirichletFixedUz } );

        fem::Set_User_Info tSetDirichletFixed09;
        tSetDirichletFixed09.set_mesh_index( tEnrIntegMesh.get_set_index_by_name("SideSet_5_c_p1") );
        tSetDirichletFixed09.set_IWGs( { tIWGDirichletFixedUz } );

        fem::Set_User_Info tSetDirichletFixed10;
        tSetDirichletFixed10.set_mesh_index( tEnrIntegMesh.get_set_index_by_name("SideSet_5_n_p3") );
        tSetDirichletFixed10.set_IWGs( { tIWGDirichletFixedUz } );

        fem::Set_User_Info tSetDirichletFixed11;
        tSetDirichletFixed11.set_mesh_index( tEnrIntegMesh.get_set_index_by_name("SideSet_5_c_p3") );
        tSetDirichletFixed11.set_IWGs( { tIWGDirichletFixedUz } );

        //==============================
        // Neumann load on side-set 2
        fem::Set_User_Info tSetNeumann00;
        tSetNeumann00.set_mesh_index( tEnrIntegMesh.get_set_index_by_name("SideSet_2_n_p1") );
        tSetNeumann00.set_IWGs( { tIWGNeumann } );

        fem::Set_User_Info tSetNeumann01;
        tSetNeumann01.set_mesh_index( tEnrIntegMesh.get_set_index_by_name("SideSet_2_c_p1") );
        tSetNeumann01.set_IWGs( { tIWGNeumann } );

        fem::Set_User_Info tSetNeumann02;
        tSetNeumann02.set_mesh_index( tEnrIntegMesh.get_set_index_by_name("SideSet_2_n_p3") );
        tSetNeumann02.set_IWGs( { tIWGNeumann } );

        fem::Set_User_Info tSetNeumann03;
        tSetNeumann03.set_mesh_index( tEnrIntegMesh.get_set_index_by_name("SideSet_2_c_p3") );
        tSetNeumann03.set_IWGs( { tIWGNeumann } );

        //==============================
        // interface(s)
        fem::Set_User_Info tInterfaceFibersToPlate;
        tInterfaceFibersToPlate.set_mesh_index( tEnrIntegMesh.get_set_index_by_name("dbl_iside_p0_1_p1_3") );
        tInterfaceFibersToPlate.set_IWGs( { tIWGFiberInterfacePlateBulk } );

        fem::Set_User_Info tInterfaceFibersToBoundary00;
        tInterfaceFibersToBoundary00.set_mesh_index( tEnrIntegMesh.get_set_index_by_name("SideSet_5_n_p1") );
        tInterfaceFibersToBoundary00.set_IWGs( { tIWGDirichletFixedUz } );

        fem::Set_User_Info tInterfaceFibersToBoundary01;
        tInterfaceFibersToBoundary01.set_mesh_index( tEnrIntegMesh.get_set_index_by_name("SideSet_5_c_p1") );
        tInterfaceFibersToBoundary01.set_IWGs( { tIWGDirichletFixedUz } );

        //==============================
        // create a cell of set info
        moris::Cell< fem::Set_User_Info > tSetInfo( 23 );
        tSetInfo( 0 )  = tBulkPlate00;
        tSetInfo( 1 )  = tBulkPlate01;

        tSetInfo( 2 )  = tBulkFibers00;
        tSetInfo( 3 )  = tBulkFibers01;

        tSetInfo( 4 )  = tSetDirichletFixed00;
        tSetInfo( 5 )  = tSetDirichletFixed01;
        tSetInfo( 6 )  = tSetDirichletFixed02;
        tSetInfo( 7 )  = tSetDirichletFixed03;
        tSetInfo( 8 )  = tSetDirichletFixed04;
        tSetInfo( 9 )  = tSetDirichletFixed05;
        tSetInfo( 10 ) = tSetDirichletFixed06;
        tSetInfo( 11 ) = tSetDirichletFixed07;
        tSetInfo( 12 ) = tSetDirichletFixed08;
        tSetInfo( 13 ) = tSetDirichletFixed09;
        tSetInfo( 14 ) = tSetDirichletFixed10;
        tSetInfo( 15 ) = tSetDirichletFixed11;

        tSetInfo( 16 ) = tSetNeumann00;
        tSetInfo( 17 ) = tSetNeumann01;
        tSetInfo( 18 ) = tSetNeumann02;
        tSetInfo( 19 ) = tSetNeumann03;

        tSetInfo( 20 ) = tInterfaceFibersToPlate;

        tSetInfo( 21 ) = tInterfaceFibersToBoundary00;
        tSetInfo( 22 ) = tInterfaceFibersToBoundary01;

        //------------------------------------------------------------------------------
        // create model
        mdl::Model * tModel = new mdl::Model( &tMeshManager,
                                               0,
                                               tSetInfo,
                                               0,
                                               false );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 1: create linear solver and algorithm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        moris::Cell< enum MSI::Dof_Type > tDofTypesU( 3 );
        tDofTypesU( 0 ) = MSI::Dof_Type::UX;
        tDofTypesU( 1 ) = MSI::Dof_Type::UY;
        tDofTypesU( 2 ) = MSI::Dof_Type::UZ;

        dla::Solver_Factory  tSolFactory;

        std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm;
        bool tDirectSolve = true;
        if(tDirectSolve)
        {
            tLinearSolverAlgorithm = tSolFactory.create_solver( sol::SolverType::AMESOS_IMPL );
        }
        else
        {
            tLinearSolverAlgorithm = tSolFactory.create_solver( sol::SolverType::AZTEC_IMPL );

            tLinearSolverAlgorithm->set_param("rel_residual")   = 6e-02;
            tLinearSolverAlgorithm->set_param("AZ_diagnostics") = AZ_none;
            tLinearSolverAlgorithm->set_param("AZ_output") = AZ_all;
            tLinearSolverAlgorithm->set_param("AZ_max_iter") = 1000;
            tLinearSolverAlgorithm->set_param("AZ_solver") = AZ_gmres;
            tLinearSolverAlgorithm->set_param("AZ_kspace") = 500;

            if(par_size()==1)
            {
                tLinearSolverAlgorithm->set_param("AZ_orthog") = AZ_modified;   // only to be used in serial
            }

            uint tPreConditioner = 0;
            switch (tPreConditioner)
            {
            case 0:
            {
                tLinearSolverAlgorithm->set_param("AZ_subdomain_solve") = AZ_ilu;
                tLinearSolverAlgorithm->set_param("AZ_graph_fill") = 5;
                break;
            }
            case 1:
            {
                tLinearSolverAlgorithm->set_param("AZ_subdomain_solve") = AZ_ilut;
                tLinearSolverAlgorithm->set_param("AZ_ilut_fill") = 10.0;
                tLinearSolverAlgorithm->set_param("AZ_drop") = 1e-12;
                tLinearSolverAlgorithm->set_param("AZ_athresh") = 0.0;
                tLinearSolverAlgorithm->set_param("AZ_rthresh") = 1.0;
                break;
            }
            default:
            {
                tLinearSolverAlgorithm->set_param("Use_ML_Prec")        = true;  // precondition the system
                tLinearSolverAlgorithm->set_param("PDE equations" )     = 1;
                tLinearSolverAlgorithm->set_param("aggregation: type")  = "Uncoupled";
            }
            }
        }

        dla::Linear_Solver tLinSolver;
        tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 2: create nonlinear solver and algorithm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        NLA::Nonlinear_Solver_Factory tNonlinFactory;
        std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );
        //        std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithmMonolythicU = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );

        tNonlinearSolverAlgorithm->set_param("NLA_max_iter")          = 10;
        tNonlinearSolverAlgorithm->set_param("NLA_rel_res_norm_drop") = 1e-3;
        //        tNonlinearSolverAlgorithmMonolythic->set_param("NLA_hard_break") = false;
        //        tNonlinearSolverAlgorithmMonolythic->set_param("NLA_max_lin_solver_restarts") = 2;
        //        tNonlinearSolverAlgorithmMonolythic->set_param("NLA_rebuild_jacobian") = true;

        tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );
        //        tNonlinearSolverAlgorithmMonolythicU->set_linear_solver( &tLinSolver );

        NLA::Nonlinear_Solver tNonlinearSolverMain;
        tNonlinearSolverMain.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );
        tNonlinearSolverMain.set_dof_type_list( tDofTypesU );

        // Create solver database
        sol::SOL_Warehouse tSolverWarehouse( tModel->get_solver_interface() );
        tNonlinearSolverMain.set_solver_warehouse( &tSolverWarehouse );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 3: create time Solver and algorithm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        tsa::Time_Solver_Factory tTimeSolverFactory;
        std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm = tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );

        tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolverMain );

        tsa::Time_Solver tTimeSolver;
        tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );
        tTimeSolver.set_solver_warehouse( &tSolverWarehouse );

        tTimeSolver.set_dof_type_list( tDofTypesU );

        //------------------------------------------------------------------------------
        tTimeSolver.solve();
        //====================================================
        //    Matrix< DDRMat > tSolVec;
        //    tTimeSolver.get_full_solution(tSolVec);
        //    print( tSolVec, " Solution Vector ");
        //====================================================

        // output solution and meshes
        xtk::Output_Options tOutputOptions;
        tOutputOptions.mAddNodeSets = false;
        tOutputOptions.mAddSideSets = true;
        tOutputOptions.mAddClusters = false;

        // add solution field to integration mesh
        std::string tIntegSolFieldNameUX = "UX";
        std::string tIntegSolFieldNameUY = "UY";
        std::string tIntegSolFieldNameUZ = "UZ";
        tOutputOptions.mRealNodeExternalFieldNames = {tIntegSolFieldNameUX, tIntegSolFieldNameUY, tIntegSolFieldNameUZ};

        moris::mtk::Integration_Mesh* tIntegMesh1 = tXTKModel.get_output_mesh(tOutputOptions);

        // Write to Integration mesh for visualization
        Matrix<DDRMat> tIntegSolUX = tModel->get_solution_for_integration_mesh_output( MSI::Dof_Type::UX );
        Matrix<DDRMat> tIntegSolUY = tModel->get_solution_for_integration_mesh_output( MSI::Dof_Type::UY );
        Matrix<DDRMat> tIntegSolUZ = tModel->get_solution_for_integration_mesh_output( MSI::Dof_Type::UZ );


        Matrix<DDRMat> tSTKIntegSolUX(tIntegMesh1->get_num_entities(EntityRank::NODE),1);
        Matrix<DDRMat> tSTKIntegSolUY(tIntegMesh1->get_num_entities(EntityRank::NODE),1);
        Matrix<DDRMat> tSTKIntegSolUZ(tIntegMesh1->get_num_entities(EntityRank::NODE),1);

        uint tNumIntegNodes = tIntegMesh1->get_num_entities(EntityRank::NODE);
        for(moris::uint i = 0; i < tNumIntegNodes; i++)
        {
            moris::moris_id tID = tIntegMesh1->get_glb_entity_id_from_entity_loc_index(i,EntityRank::NODE);
            tSTKIntegSolUX(i) = tIntegSolUX(tEnrIntegMesh.get_loc_entity_ind_from_entity_glb_id(tID,EntityRank::NODE));
            tSTKIntegSolUY(i) = tIntegSolUY(tEnrIntegMesh.get_loc_entity_ind_from_entity_glb_id(tID,EntityRank::NODE));
            tSTKIntegSolUZ(i) = tIntegSolUZ(tEnrIntegMesh.get_loc_entity_ind_from_entity_glb_id(tID,EntityRank::NODE));
        }

        // add solution field to integration mesh
        tIntegMesh1->add_mesh_field_real_scalar_data_loc_inds(tIntegSolFieldNameUX,EntityRank::NODE,tSTKIntegSolUX);
        tIntegMesh1->add_mesh_field_real_scalar_data_loc_inds(tIntegSolFieldNameUY,EntityRank::NODE,tSTKIntegSolUY);
        tIntegMesh1->add_mesh_field_real_scalar_data_loc_inds(tIntegSolFieldNameUZ,EntityRank::NODE,tSTKIntegSolUZ);

        //------------------------------------------------------------------------------
        // output solution to mesh
        std::string tMeshOutputFile = "fiber_test.e";

        std::cout<<"output file name:  "<<tMeshOutputFile<<std::endl;

        tIntegMesh1->create_output_mesh(tMeshOutputFile);
        delete tIntegMesh1;
        delete tInterpMesh;
    }// end full problem if-statement
}

TEST_CASE("mesh_writing", "[GE],[fiber_write_mesh]")
{

}

}   // ge namespace
