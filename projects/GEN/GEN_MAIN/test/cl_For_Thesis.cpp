/*
 * cl_For_Thesis.cpp
 *
 *  Created on: Jan 21, 2020
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
#include "cl_MTK_Writer_Exodus.hpp"
#include "cl_MTK_Reader_Exodus.hpp"

#include "cl_Matrix.hpp"        //LINALG
#include "linalg_typedefs.hpp"
#include "fn_equal_to.hpp" // ALG/src

#include "cl_FEM_CM_Factory.hpp"
#include "cl_FEM_Element_Factory.hpp"
#include "cl_FEM_ElementProxy.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Factory.hpp"
#include "cl_FEM_IWG_Factory.hpp"
#include "cl_FEM_NodeProxy.hpp"
#include "cl_FEM_Node_Base.hpp"
#include "cl_FEM_SP_Factory.hpp"
#include "cl_FEM_Set_User_Info.hpp"

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

#include "cl_SOL_Warehouse.hpp"

#include "cl_TSA_Time_Solver_Factory.hpp"
#include "cl_TSA_Monolithic_Time_Solver.hpp"
#include "cl_TSA_Time_Solver.hpp"

#include "fn_norm.hpp"

#include "cl_GEN_Circle.hpp"
#include "cl_GEN_Enums.hpp"
#include "cl_GEN_Geometry.hpp"
#include "cl_GEN_Geom_Data.hpp"
#include "cl_GEN_Geom_Field_HMR.hpp"
#include "cl_GEN_Multi_Geometry.hpp"
#include "cl_GEN_Plane.hpp"
#include "cl_GEN_Property.hpp"
#include "cl_GEN_Sphere.hpp"

#include "cl_VIS_Factory.hpp"
#include "cl_VIS_Output_Manager.hpp"

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

void tConstValFunctionBottom
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    Matrix< DDRMat > tOutput(3,1);
    tOutput(0) = 0.0; tOutput(1) = 0.0; tOutput(2) = 0.0;

    Matrix< DDRMat > tCoords = aFIManager->get_IP_geometry_interpolator()->valx();


    if( tCoords(0) > 0.25 && tCoords(1) > 0.25 )
    {
        aPropMatrix = aParameters( 0 );
    }
    else
    {
        aPropMatrix = tOutput;
    }
}

void tMValFunctionUnique
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    Matrix< DDRMat > tIMat(3,3, 0.0);

    Matrix< DDRMat > tCoords = aFIManager->get_IP_geometry_interpolator()->valx();

    if( tCoords(0) < 0.25 && tCoords(1) < 0.25 )
    {
        tIMat(0,0) = aParameters( 0 )( 0 );
        tIMat(1,1) = aParameters( 0 )( 1 );
        tIMat(2,2) = aParameters( 0 )( 2 );
    }

    aPropMatrix = tIMat;
}

void tMValFunction
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = {{ aParameters( 0 )( 0 ),                      0.0,                     0.0 },
                   {                   0.0,    aParameters( 0 )( 1 ),                     0.0 },
                   {                   0.0,                      0.0,    aParameters( 0 )( 2 )} };
}

void tMValFunction2D
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = {{ aParameters( 0 )( 0 ),                      0.0},
                   {                   0.0,    aParameters( 0 )( 1 )}};
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

bool tSolverOutputCriteriaThesis( moris::tsa::Time_Solver * )
{
    return true;
}





//------------------------------------------------------------------------------
TEST_CASE("experiments for thesis, geom.", "[GE],[thesis_01]")
{
    if(par_size()<=1)
    {
        size_t tModelDimension  = 2;
        uint tLagrangeMeshIndex = 0;
        //  HMR Parameters setup
        moris::ParameterList tParameters = prm::create_hmr_parameter_list();

    uint tInitialMesh = 1;
    switch(tInitialMesh)
    {
    case(0) :
        {
            tParameters.set( "number_of_elements_per_dimension", std::string("40, 40, 20") );
            tParameters.set( "domain_dimensions",                std::string("2, 2, 1") );
            tParameters.set( "domain_offset",                    std::string("-0, -0, -0") );
            tParameters.set( "domain_sidesets",            std::string("1, 2, 3, 4, 5, 6") );
            break;
        }
    case(1) :
        {
        // mesh used in Moes et. al. for the first benchmark problem
            tParameters.set( "number_of_elements_per_dimension", std::string("24, 48") );
            tParameters.set( "domain_dimensions",                std::string("7, 16") );
            tParameters.set( "domain_offset",                    std::string("0, 0") );
            tParameters.set( "domain_sidesets", std::string("1, 2, 3, 4") );
        break;
        }
    default :
        {
            tParameters.set( "number_of_elements_per_dimension", std::string("20, 20") );
            tParameters.set( "domain_dimensions",                std::string(" 2,  2") );
            tParameters.set( "domain_offset",                    std::string(" 0,  0") );
            tParameters.set( "domain_sidesets", std::string("1, 2, 3, 4") );
        }
    }

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

        //  HMR Initialization
        moris::hmr::HMR tHMR( tParameters );

        auto tDatabase = tHMR.get_database(); // std::shared_ptr< Database >

        tHMR.perform_initial_refinement( 0 );
        //------------------------------------------------------------------------------
        tDatabase->update_bspline_meshes();
        tDatabase->update_lagrange_meshes();

        uint tNumberOfRefinements = 1;
        uint tStartingPoint       = 1;

        Cell< Matrix< DDRMat > > tFieldData( 1 );

        std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );
        std::cout<<"-------------------------------------"<<std::endl;
        std::cout<<"number of refinements performed:  "<<tNumberOfRefinements<<std::endl;
        std::cout<<"-------------------------------------"<<std::endl;
        for( uint k=0; k<tNumberOfRefinements; k++ )
        {
            real tCrackX;
            real tCrackY;

            switch(tInitialMesh)
            {
            case(0) :
                {
                    tCrackX = 0.0;
                    tCrackY = 0.0;
                    break;
                }
            case(1) :   // crack tip begins in middle of mesh (Moes et. al.)
                {
                    tCrackX = 3.5001;
                    tCrackY = 8.0001;
                    break;
                }
            default :
                {
                    tCrackX = 0.1502;
                    tCrackY = 1.0001;
                    break;
                }
            }

            Matrix< DDRMat > tCenters  = {{ tCrackX,  tCrackY }};
            Matrix< DDRMat > tNormals2 = {{ 1.0, 0.0 }};
            moris::ge::Plane<2> tPlane2( tCenters, tNormals2 );

//            Matrix< DDRMat > tCenter0  = {{ 0.0, tCrackY+0.02 }};
            Matrix< DDRMat > tCenter0  = {{ 0.0, tCrackY }};
            Matrix< DDRMat > tNormals0 = {{ 0.0,  1.0 }};
            moris::ge::Plane<2> tPlane0( tCenter0, tNormals0 );

//            Matrix< DDRMat > tCenter1  = {{ 0.0, tCrackY-0.02 }};
            Matrix< DDRMat > tCenter1  = {{ 0.0, tCrackY }};
            Matrix< DDRMat > tNormals1 = {{ 0.0, -1.0 }};
            moris::ge::Plane<2> tPlane1( tCenter1, tNormals1 );

            moris::Cell< moris::ge::GEN_Geometry* > tAllPlanes = { &tPlane0, &tPlane1, &tPlane2 };
            //------------------------------------------------------------------------------
//            Matrix< DDRMat > tCenters  = {{ tCrackX,  tCrackY }};
//            Matrix< DDRMat > tNormals2 = {{ 1.0, 0.0 }};
//            moris::ge::Plane<2> tPlane2( tCenters, tNormals2 );
//
//            Matrix< DDRMat > tCenter0  = {{ 0.0, tCrackY }};
//            Matrix< DDRMat > tNormals0 = {{ 0.0,  1.0 }};
//            moris::ge::Plane<2> tPlane0( tCenter0, tNormals0 );
//
//            moris::Cell< moris::ge::GEN_Geometry* > tAllPlanes = { &tPlane0, &tPlane2 };
            //------------------------------------------------------------------------------
            moris::ge::Multi_Geometry tCrack( tAllPlanes );

//            moris::ge::Circle tCircle( 0.12501, tCrackX, tCrackY );

//            moris::Cell< moris::ge::GEN_Geometry* > tTempGeometryVector = { &tCrack, &tCircle };
            moris::Cell< moris::ge::GEN_Geometry* > tTempGeometryVector = { &tCrack };

            moris::ge::GEN_Phase_Table     tTempPhaseTable( tTempGeometryVector.size(),  Phase_Table_Structure::EXP_BASE_2 );
            moris::ge::Geometry_Engine tTempGeometryEngine( tTempGeometryVector, tTempPhaseTable, tModelDimension );

            moris_index tMeshIndex = tTempGeometryEngine.register_mesh( tMesh );

            uint tNumIPNodes = tMesh->get_num_nodes();

            Matrix<DDRMat> tFieldData ( tNumIPNodes,1 );
//            Matrix<DDRMat> tFieldData0( tNumIPNodes,1 );

            tTempGeometryEngine.initialize_geometry_objects_for_background_mesh_nodes( tNumIPNodes );
            Matrix< DDRMat > tCoords( tNumIPNodes, tModelDimension );
            for( uint i = 0; i < tNumIPNodes; i++ )
            {
                tCoords.set_row( i, tMesh->get_mtk_vertex(i).get_coords() );
            }

            tTempGeometryEngine.initialize_geometry_object_phase_values( tCoords );

            for(uint i=0; i<tNumIPNodes; i++)
            {
                tFieldData ( i ) = tTempGeometryEngine.get_entity_phase_val( i, 0 );
//                tFieldData0( i ) = tTempGeometryEngine.get_entity_phase_val( i, 1 );
            }

            tHMR.based_on_field_put_elements_on_queue( tFieldData, tLagrangeMeshIndex );
//            tHMR.based_on_field_put_elements_on_queue( tFieldData0, tLagrangeMeshIndex );

            tHMR.perform_refinement_based_on_working_pattern( 0, false );
        }

        tHMR.finalize();

        hmr::Interpolation_Mesh_HMR *      tInterpMesh      = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );

        //------------------------------------------------------------------------------
        real tCrackX;
        real tCrackY;

        switch(tInitialMesh)
        {
        case(0) :
            {
                tCrackX = 0.0;
                tCrackY = 0.0;
                break;
            }
        case(1) :   // crack tip begins in middle of mesh (Moes et. al.)
            {
                tCrackX = 3.5001;
                tCrackY = 8.0001;
                break;
            }
        default :
            {
                tCrackX = 0.1502;
                tCrackY = 1.0001;
                break;
            }
        }

        Matrix< DDRMat > tCenters  = {{ tCrackX,  tCrackY }};
        Matrix< DDRMat > tNormals2 = {{ 1.0, 0.0 }};
        moris::ge::Plane<2> tPlane2( tCenters, tNormals2 );

        Matrix< DDRMat > tCenter0  = {{ 0.0, tCrackY }};
        Matrix< DDRMat > tNormals0 = {{ 0.0,  1.0 }};
        moris::ge::Plane<2> tPlane0( tCenter0, tNormals0 );

        Matrix< DDRMat > tCenter1  = {{ 0.0, tCrackY }};
        Matrix< DDRMat > tNormals1 = {{ 0.0, -1.0 }};
        moris::ge::Plane<2> tPlane1( tCenter1, tNormals1 );

        moris::Cell< moris::ge::GEN_Geometry* > tAllPlanes = { &tPlane0, &tPlane1, &tPlane2 };
        //------------------------------------------------------------------------------
//        Matrix< DDRMat > tCenters  = {{ tCrackX,  tCrackY }};
//        Matrix< DDRMat > tNormals2 = {{ 1.0, 0.0 }};
//        moris::ge::Plane<2> tPlane2( tCenters, tNormals2 );
//
//        Matrix< DDRMat > tCenter0  = {{ 0.0, tCrackY }};
//        Matrix< DDRMat > tNormals0 = {{ 0.0,  1.0 }};
//        moris::ge::Plane<2> tPlane0( tCenter0, tNormals0 );
//
//        moris::Cell< moris::ge::GEN_Geometry* > tAllPlanes = { &tPlane0, &tPlane2 };
        //------------------------------------------------------------------------------

        moris::ge::Multi_Geometry tCrack( tAllPlanes );
        //===========================================
//        moris::ge::Circle tCircle( 0.12501, tCrackX, tCrackY );
//        moris::Cell< moris::ge::GEN_Geometry* > tGeometryVector = { &tCrack, &tCircle };
        //===========================================

        moris::Cell< moris::ge::GEN_Geometry* > tGeometryVector = { &tCrack };

        moris::ge::GEN_Phase_Table      tPhaseTable( tGeometryVector.size(),  Phase_Table_Structure::EXP_BASE_2 );
        moris::ge::Geometry_Engine  tGENGeometryEngine( tGeometryVector, tPhaseTable, tModelDimension );

        //------------------------------------------------------------------------------
        xtk::Model tXTKModel( tModelDimension, tInterpMesh, &tGENGeometryEngine );

        tXTKModel.mVerbose = false;

        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3};
        tXTKModel.decompose(tDecompositionMethods);

        // ============================ output problem geometry ===================================
        bool tOutputXTKmesh = true;
        if (tOutputXTKmesh)
        {
            tXTKModel.perform_basis_enrichment(EntityRank::BSPLINE, 0);

            xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
            xtk::Enriched_Integration_Mesh   & tEnrIntegMesh = tXTKModel.get_enriched_integ_mesh();

            // Write mesh
            moris::mtk::Writer_Exodus writer(&tEnrIntegMesh);
            writer.write_mesh("", "0_geomCheckNoFibers.exo");

            // Write the fields
            writer.set_time(0.0);

            writer.close_file();
        }
        // ============================= end geometry output ======================================
        bool tFullProblem = false;
        if(tFullProblem)
        {
            tXTKModel.perform_basis_enrichment(EntityRank::BSPLINE, 0);

            xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
            xtk::Enriched_Integration_Mesh   & tEnrIntegMesh = tXTKModel.get_enriched_integ_mesh();

            // place the pair in mesh manager
            mtk::Mesh_Manager tMeshManager;
            tMeshManager.register_mesh_pair( &tEnrInterpMesh, &tEnrIntegMesh);
            //------------------------------------------------------------------------------
            // create residual dof types
            moris::Cell< MSI::Dof_Type > tResDofTypes = { MSI::Dof_Type::UX, MSI::Dof_Type::UY };
            // create the material properties
            real tEModPlate  =  100.0;
            real tNuPlate    =    0.3;
            std::cout<<"----------------------------------"<<std::endl;
            std::cout<<"E plate:   "<<tEModPlate<<std::endl;
            std::cout<<"----------------------------------"<<std::endl;
            std::cout<<"nu plate:  "<<tNuPlate<<std::endl;
            std::cout<<"----------------------------------"<<std::endl;
            //------------------------------------------------------------------------------
            std::shared_ptr< fem::Property > tPropEModPlate = std::make_shared< fem::Property >();
            tPropEModPlate->set_parameters( { {{ tEModPlate }} } );
            tPropEModPlate->set_val_function( tConstValFunction );

            std::shared_ptr< fem::Property > tPropNuPlate = std::make_shared< fem::Property >();
            tPropNuPlate->set_parameters( { {{ tNuPlate }} } );
            tPropNuPlate->set_val_function( tConstValFunction );
            //------------------------------------------------------------------------------
            // loading on top
            std::shared_ptr< fem::Property > tPropNeumannTop = std::make_shared< fem::Property >();
            tPropNeumannTop->set_parameters( {{{ 0.0 }, { 1000.0 }}} );
            tPropNeumannTop->set_val_function( tConstValFunction );
            //------------------------------------------------------------------------------
            // fixed bottom
            std::shared_ptr< fem::Property > tPropDirichlet_ss1 = std::make_shared< fem::Property >();
            tPropDirichlet_ss1->set_parameters( { {{ 0.0 }, { 0.0 }} } );
            tPropDirichlet_ss1->set_val_function( tConstValFunction );

            std::shared_ptr< fem::Property > tPropDirichlet_ss1_select = std::make_shared< fem::Property >();
            tPropDirichlet_ss1_select->set_parameters( { {{ 1.0 }, { 1.0 }} } );
            tPropDirichlet_ss1_select->set_val_function( tMValFunction2D );
            //------------------------------------------------------------------------------
            // define constitutive models
            fem::CM_Factory tCMFactory;

            std::shared_ptr< fem::Constitutive_Model > tCMPlate = tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
            tCMPlate->set_dof_type_list( { tResDofTypes } );
            tCMPlate->set_property( tPropEModPlate, "YoungsModulus" );
            tCMPlate->set_property( tPropNuPlate, "PoissonRatio" );
            tCMPlate->set_space_dim( tModelDimension );
            //------------------------------------------------------------------------------
            // define stabilization parameters
            fem::SP_Factory tSPFactory;

            // stabilization parameter for fixed bottom
            std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitscheBCs = tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
            tSPDirichletNitscheBCs->set_parameters( { {{ 10.0 }} } );
            tSPDirichletNitscheBCs->set_property( tPropEModPlate, "Material", mtk::Master_Slave::MASTER );
            //------------------------------------------------------------------------------
            // create the IQIs
            fem::IQI_Factory tIQIFactory;

            std::shared_ptr< fem::IQI > tIQIUX = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
            tIQIUX->set_output_type( vis::Output_Type::UX );
            tIQIUX->set_dof_type_list( { tResDofTypes }, mtk::Master_Slave::MASTER );
            tIQIUX->set_output_type_index( 0 );

            std::shared_ptr< fem::IQI > tIQIUY = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
            tIQIUY->set_output_type( vis::Output_Type::UY );
            tIQIUY->set_dof_type_list( { tResDofTypes }, mtk::Master_Slave::MASTER );
            tIQIUY->set_output_type_index( 1 );

            std::shared_ptr< fem::IQI > tIQIJInt = tIQIFactory.create_IQI( fem::IQI_Type::J_INTEGRAL );
            tIQIJInt->set_output_type( vis::Output_Type::J_INTEGRAL );
            tIQIJInt->set_constitutive_model( tCMPlate, "ElastLinIso", mtk::Master_Slave::MASTER );
            tIQIJInt->set_dof_type_list( { tResDofTypes }, mtk::Master_Slave::MASTER );

            //------------------------------------------------------------------------------
            // define the IWGs
            fem::IWG_Factory tIWGFactory;

            std::shared_ptr< fem::IWG > tIWGPlate = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_BULK );
            tIWGPlate->set_residual_dof_type( tResDofTypes );
            tIWGPlate->set_dof_type_list( { tResDofTypes } );
            tIWGPlate->set_constitutive_model( tCMPlate, "ElastLinIso", mtk::Master_Slave::MASTER );
            //------------------------------------------------------------------------------
            std::shared_ptr< fem::IWG > tIWGNeumannTop = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_NEUMANN );
            tIWGNeumannTop->set_residual_dof_type( tResDofTypes );
            tIWGNeumannTop->set_dof_type_list( { tResDofTypes } );
            tIWGNeumannTop->set_property( tPropNeumannTop, "Neumann", mtk::Master_Slave::MASTER );
            //------------------------------------------------------------------------------
            std::shared_ptr< fem::IWG > tIWGDirichletFixedBottom = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE );
            tIWGDirichletFixedBottom->set_residual_dof_type( tResDofTypes );
            tIWGDirichletFixedBottom->set_dof_type_list( { tResDofTypes } );
            tIWGDirichletFixedBottom->set_stabilization_parameter( tSPDirichletNitscheBCs, "DirichletNitsche" );
            tIWGDirichletFixedBottom->set_constitutive_model( tCMPlate, "ElastLinIso", mtk::Master_Slave::MASTER );
            tIWGDirichletFixedBottom->set_property( tPropDirichlet_ss1, "Dirichlet", mtk::Master_Slave::MASTER );
            tIWGDirichletFixedBottom->set_property( tPropDirichlet_ss1_select, "Select", mtk::Master_Slave::MASTER );
            //------------------------------------------------------------------------------
            //===========================================
            // bulk for plate
            fem::Set_User_Info tBulkPlate00;
            tBulkPlate00.set_mesh_set_name( "HMR_dummy_n_p3" );
                tBulkPlate00.set_IWGs( { tIWGPlate } );
                tBulkPlate00.set_IQIs( { tIQIUX, tIQIUY, tIQIJInt } );

            fem::Set_User_Info tBulkPlate01;
            tBulkPlate01.set_mesh_set_name( "HMR_dummy_c_p3" );
            tBulkPlate01.set_IWGs( { tIWGPlate } );
            tBulkPlate01.set_IQIs( { tIQIUX, tIQIUY, tIQIJInt } );
            //===========================================
            // Neumann load on side-set 3
            fem::Set_User_Info tSetNeumann00;
            tSetNeumann00.set_mesh_set_name( "SideSet_3_n_p3" );
            tSetNeumann00.set_IWGs( { tIWGNeumannTop } );

            fem::Set_User_Info tSetNeumann01;
            tSetNeumann01.set_mesh_set_name( "SideSet_3_c_p3" );
            tSetNeumann01.set_IWGs( { tIWGNeumannTop } );
            //===========================================
            // boundary conditions on side-set 1
            fem::Set_User_Info tSetDirichletFixed00;
            tSetDirichletFixed00.set_mesh_set_name( "SideSet_1_n_p3" );
            tSetDirichletFixed00.set_IWGs( { tIWGDirichletFixedBottom } );

            fem::Set_User_Info tSetDirichletFixed01;
            tSetDirichletFixed01.set_mesh_set_name( "SideSet_1_c_p3" );
            tSetDirichletFixed01.set_IWGs( { tIWGDirichletFixedBottom } );
            //===========================================
//             IQI for J-Integral
//            fem::Set_User_Info tJIntegral;
//            tJIntegral.set_mesh_set_name( "iside_g_1_b0_3_b1_2" );
//            tJIntegral.set_IQIs( { tIQIJInt } );
            //------------------------------------------------------------------------------
            // create a cell of set info
            moris::Cell< fem::Set_User_Info > tSetInfo( 6 );
            tSetInfo( 0 )  = tBulkPlate00;
            tSetInfo( 1 )  = tBulkPlate01;
            tSetInfo( 2 )  = tSetNeumann00;
            tSetInfo( 3 )  = tSetNeumann01;
            tSetInfo( 4 )  = tSetDirichletFixed00;
            tSetInfo( 5 )  = tSetDirichletFixed01;
//            tSetInfo( 6 )  = tJIntegral;
            //------------------------------------------------------------------------------
            // create model
            mdl::Model * tModel = new mdl::Model( &tMeshManager,
                                                   0,
                                                   tSetInfo,
                                                   0,
                                                   false );
            // --------------------------------------------------------------------------------------
//            moris_index tSideSetForJ = tEnrIntegMesh.get_set_index_by_name("iside_g_1_b0_3_b1_2");
//            moris_index tSideSetForJ = tEnrIntegMesh.get_set_index_by_name("HMR_dummy_n_p3");
//
//            map< moris_index, moris_index > tTempMap = tModel->get_mesh_set_to_fem_set_index_map();
//
//            auto tIter = tTempMap[ tSideSetForJ ];
//
//            Matrix< DDRMat > tDummy0(1,1,0.0);
//            Matrix< DDRMat > tDummy1(1,1,0.0);
//
//            real tJVal;

//            tModel->get_equation_sets()(tIter)->set_visualization_set( 0, tEnrIntegMesh.get_set_by_index( tSideSetForJ ), true );

//            tModel->get_equation_sets()(tIter)->compute_quantity_of_interest( 0,
//                                                    &tDummy0,
//                                                    &tDummy1,
//                                                    &tJVal,
//                                                    vis::Output_Type::J_INTEGRAL,
//                                                    vis::Field_Type::GLOBAL );
//
//            std::cout<<"value of J-Integral:  "<<tJVal<<std::endl;
            // --------------------------------------------------------------------------------------
            // Define outputs
            vis::Output_Manager tOutputData;

            tOutputData.set_outputs( 0,
                                     vis::VIS_Mesh_Type::STANDARD,
                                     "./",
                                     "aaaaaaaaaa_outputCheck2D.e",
                                     { "HMR_dummy_n_p3", "HMR_dummy_c_p3" },
                                     { "UX", "UY", "J_INTEGRAL" },
                                     { vis::Field_Type::NODAL, vis::Field_Type::NODAL, vis::Field_Type::GLOBAL},
                                     { vis::Output_Type::UX, vis::Output_Type::UY, vis::Output_Type::J_INTEGRAL} );

            tModel->set_output_manager( &tOutputData );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // SOLVER STEP 1: create linear solver and algorithm
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            moris::Cell< enum MSI::Dof_Type > tDofTypesU( 2 );
            tDofTypesU( 0 ) = MSI::Dof_Type::UX;
            tDofTypesU( 1 ) = MSI::Dof_Type::UY;

            dla::Solver_Factory  tSolFactory;

            std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm;

            bool tDirectSolve = true;
            if (tDirectSolve)
            {
                tLinearSolverAlgorithm = tSolFactory.create_solver( sol::SolverType::AMESOS_IMPL );
            }
            else
            {
                tLinearSolverAlgorithm = tSolFactory.create_solver( sol::SolverType::AZTEC_IMPL );

                tLinearSolverAlgorithm->set_param("rel_residual")   = 6e-02;
                tLinearSolverAlgorithm->set_param("AZ_diagnostics") = AZ_none;
                tLinearSolverAlgorithm->set_param("AZ_output") = AZ_all;        // AZ_none
                tLinearSolverAlgorithm->set_param("AZ_max_iter") = 1000;
                tLinearSolverAlgorithm->set_param("AZ_solver") = AZ_gmres;
                tLinearSolverAlgorithm->set_param("AZ_kspace") = 500;
                tLinearSolverAlgorithm->set_param("AZ_orthog") = AZ_modified;   // only to be used in serial
                //    tLinearSolverAlgorithm->set_param("AZ_solver") = AZ_gmres_condnum;

                uint tPreConditioner = 0;

                switch (tPreConditioner)
                {
                case 0:
                {
                    tLinearSolverAlgorithm->set_param("AZ_subdomain_solve") = AZ_ilu;
                    tLinearSolverAlgorithm->set_param("AZ_graph_fill") = 3;
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
            // SOLVER STEP 2: create nonlinear solver and algorithm
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            NLA::Nonlinear_Solver_Factory tNonlinFactory;
            std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );

            tNonlinearSolverAlgorithm->set_param("NLA_max_iter")          = 10;
            tNonlinearSolverAlgorithm->set_param("NLA_tot_res_norm_drop") = 1e-3;

            tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );

            NLA::Nonlinear_Solver tNonlinearSolverMain;
            tNonlinearSolverMain.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );
            tNonlinearSolverMain.set_dof_type_list( tDofTypesU );

            // Create solver database
            sol::SOL_Warehouse tSolverWarehouse( tModel->get_solver_interface() );
            tNonlinearSolverMain.set_solver_warehouse( &tSolverWarehouse );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // SOLVER STEP 3: create time Solver and algorithm
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            tsa::Time_Solver_Factory tTimeSolverFactory;
            std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm = tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );

            tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolverMain );

            tsa::Time_Solver tTimeSolver;
            tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );
            tTimeSolver.set_solver_warehouse( &tSolverWarehouse );

            tTimeSolver.set_dof_type_list( tDofTypesU );

            tTimeSolver.set_output( 0, tSolverOutputCriteriaThesis );
            //------------------------------------------------------------------------------
            tTimeSolver.solve();

            delete tInterpMesh;
        }   // end full problem logic statement

    } // end par size statement
}
//------------------------------------------------------------------------------
TEST_CASE("experiments for thesis", "[GE],[thesis_00]")
{
    /*
     * place holder for the implementation of the fibers into the problem
     */
    bool tRunProblem = false;

    if(tRunProblem)
    {
    if(par_size()<=1)
    {
        uint tLagrangeMeshIndex = 0;
        //  HMR Parameters setup
        moris::ParameterList tParameters = prm::create_hmr_parameter_list();

        uint tInitialMesh = 2;
        switch(tInitialMesh)
        {
        case(0) :
            {
            tParameters.set( "number_of_elements_per_dimension", std::string("40, 40, 20") );
            break;
            }
        case(1) :
            {
            tParameters.set( "number_of_elements_per_dimension", std::string("20, 20, 10") );
            }
        default :
        {
            tParameters.set( "number_of_elements_per_dimension", std::string("10, 10, 5") );
        }
        }

        tParameters.set( "domain_dimensions",                std::string("2, 2, 1") );
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

        uint tNumberOfFibers      = 5;
        uint tNumberOfRefinements = 1;
        std::cout<<"-------------------------------------------"<<std::endl;
        std::cout<<"number of fibers being analyzed:      "<<tNumberOfFibers<<std::endl;
        std::cout<<"number of refinements performed:      "<<tNumberOfRefinements<<std::endl;
        std::cout<<"bandwidth in user-defined refinement: "<<gLsbwabs<<std::endl;
        std::cout<<"-------------------------------------------"<<std::endl;
        Cell< Matrix< DDRMat > > tFieldData( 1 );

        std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

        moris::tic tTimer0;
        for( uint k=0; k<tNumberOfRefinements; k++ )
        {
            moris::ge::GEN_CylinderWithEndCaps    tFibers( tNumberOfFibers, true, 0.0499 );
            moris::Cell<moris::ge::GEN_Geometry*> tGeometryVector = { &tFibers };

            moris::ge::GEN_Phase_Table         tPhaseTable( tGeometryVector.size(),  Phase_Table_Structure::EXP_BASE_2 );
            moris::ge::Geometry_Engine     tGENGeometryEngine( tGeometryVector,tPhaseTable,3 );

            moris_index tMeshIndex = tGENGeometryEngine.register_mesh( tMesh );

            tFieldData( 0 ) = tGENGeometryEngine.get_cylinder_vals( tMeshIndex, &tFibers, tNumberOfFibers );

            tHMR.user_defined_flagging( user_defined_refinement, tFieldData, tParameters, 0 );

            tHMR.perform_refinement_based_on_working_pattern( 0, false );
        }

        tHMR.finalize();

        moris::ge::GEN_CylinderWithEndCaps    tFibers( tNumberOfFibers, true, 0.0499 );
        moris::Cell<moris::ge::GEN_Geometry*> tGeometryVector_temp = {&tFibers};

        moris::ge::GEN_Phase_Table         tPhaseTable_temp( tGeometryVector_temp.size(),  Phase_Table_Structure::EXP_BASE_2 );
        moris::ge::Geometry_Engine     tGENGeometryEngine_temp( tGeometryVector_temp, tPhaseTable_temp, 3 );

        moris_index tMeshIndex = tGENGeometryEngine_temp.register_mesh( tMesh );

        tFieldData( 0 ) = tGENGeometryEngine_temp.get_cylinder_vals( tMeshIndex, &tFibers, tNumberOfFibers );

        hmr::Interpolation_Mesh_HMR *      tInterpMesh      = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );

        real tElapsedTime0 = tTimer0.toc<moris::chronos::milliseconds>().wall;
        tElapsedTime0 /= 1000;
        std::cout<<"==============================================="<< std::endl;
        std::cout<<"Total time for evaluation: "<< tElapsedTime0 << std::endl << std::endl;
        std::cout<<"==============================================="<< std::endl;
        //------------------------------------------------------------------------------

        moris::ge::GEN_Geom_Data tFiberData( tFieldData( 0 ) );  // fiber LS values as field data

        //===========================================
        Matrix< DDRMat > tCenters0 = {{ 0.0502,  1.0001 }};
        Matrix< DDRMat > tNormals0 = {{ 0.1699, -0.9855 }};
        moris::ge::Plane<2> tPlane0( tCenters0, tNormals0 );

        Matrix< DDRMat > tCenters1 = {{ 0.0502, 1.0001 }};
        Matrix< DDRMat > tNormals1 = {{ 0.1699, 0.9855 }};
        moris::ge::Plane<2> tPlane1( tCenters1, tNormals1 );

        moris::Cell< moris::ge::GEN_Geometry* > tBothPlanes = { &tPlane0, &tPlane1 };

        moris::ge::Multi_Geometry tCrack( tBothPlanes );
        //===========================================

        moris::Cell< moris::ge::GEN_Geometry* > tGeometryVector = { &tFiberData, &tCrack };

        size_t tModelDimension = 3;
        moris::ge::GEN_Phase_Table      tPhaseTable( tGeometryVector.size(),  Phase_Table_Structure::EXP_BASE_2 );
        moris::ge::Geometry_Engine  tGENGeometryEngine( tGeometryVector, tPhaseTable, tModelDimension );

        //------------------------------------------------------------------------------
        xtk::Model tXTKModel( tModelDimension, tInterpMesh, &tGENGeometryEngine );

        tXTKModel.mVerbose = false;

        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4};
        tXTKModel.decompose(tDecompositionMethods);
        //=============================== temporary ============================================
        // output problem geometry
        bool tOutputXTKmesh = false;
        if (tOutputXTKmesh)
        {
            tXTKModel.perform_basis_enrichment(EntityRank::BSPLINE, 0);

            xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
            xtk::Enriched_Integration_Mesh   & tEnrIntegMesh = tXTKModel.get_enriched_integ_mesh();

            // Write mesh
            moris::mtk::Writer_Exodus writer(&tEnrIntegMesh);
            writer.write_mesh("", "0_geomCheck.exo");

            // Write the fields
            writer.set_time(1.0);

            writer.close_file();
        }
        delete tInterpMesh;
        //============================= end temporary ==========================================
    }   //end par_size() statement

    }   //end run statement
}



//------------------------------------------------------------------------------
}   // end ge namespace
