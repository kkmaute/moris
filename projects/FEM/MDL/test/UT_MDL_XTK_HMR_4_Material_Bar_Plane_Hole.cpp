/*
 * UT_MDL_XTK_HMR_Multi_Material_Bar_Plane.cpp
 *
 *  Created on: Oct 4, 2019
 *      Author: doble
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
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"

#include "cl_Matrix.hpp"        //LINALG
#include "linalg_typedefs.hpp"
#include "fn_equal_to.hpp" // ALG/src

#include "cl_FEM_NodeProxy.hpp"                //FEM/INT/src
#include "cl_FEM_ElementProxy.hpp"             //FEM/INT/src
#include "cl_FEM_Node_Base.hpp"                //FEM/INT/src
#include "cl_FEM_Element_Factory.hpp"          //FEM/INT/src
#include "cl_FEM_IWG_Factory.hpp"              //FEM/INT/src
#include "cl_FEM_CM_Factory.hpp"              //FEM/INT/src
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

#include "../projects/GEN/src/ripped/geometry/cl_GEN_Geometry.hpp"
#include "../projects/GEN/src/ripped/geometry/cl_GEN_Geom_Field.hpp"

#include "fn_norm.hpp"

moris::real
Plane4MatMDL1(const moris::Matrix< moris::DDRMat > & aPoint )
{
    moris::real mXC = 0.1;
    moris::real mYC = 0.1;
    moris::real mNx = 1.0;
    moris::real mNy = 0.0;
    return (mNx*(aPoint(0)-mXC) + mNy*(aPoint(1)-mYC));
}

moris::real
Plane4MatMDL2(const moris::Matrix< moris::DDRMat > & aPoint )
{
    moris::real mXC = 1.4;
    moris::real mYC = 0.0;
    moris::real mNx = 1.0;
    moris::real mNy = 0.0;
    return (mNx*(aPoint(0)-mXC) + mNy*(aPoint(1)-mYC));
}

moris::real
Plane4MatMDL3(const moris::Matrix< moris::DDRMat > & aPoint )
{
    moris::real mXC = -1.4;
    moris::real mYC = 0.0;
    moris::real mNx = 1.0;
    moris::real mNy = 0.0;
    return (mNx*(aPoint(0)-mXC) + mNy*(aPoint(1)-mYC));
}



moris::real
Circle4MatMDL(const moris::Matrix< moris::DDRMat > & aPoint )
{
    moris::real mXCenter = 0.01;
    moris::real mYCenter = 0.01;
    moris::real mRadius = 0.47334;


    return  (aPoint(0) - mXCenter) * (aPoint(0) - mXCenter)
                    + (aPoint(1) - mYCenter) * (aPoint(1) - mYCenter)
                    - (mRadius * mRadius);
}


Matrix< DDRMat > tConstValFunction2MatMDL( moris::Cell< Matrix< DDRMat > >         & aParameters,
                                           moris::Cell< fem::Field_Interpolator* > & aDofFI,
                                           moris::Cell< fem::Field_Interpolator* > & aDvFI,
                                           fem::Geometry_Interpolator              * aGeometryInterpolator )
{
    return aParameters( 0 );
}


void
run_hmr_for_multi_mat_model_2d(hmr::HMR  &                    aHMR,
                               Cell<std::shared_ptr< moris::hmr::Field >> & aFields)
{
    moris_index tLagrangeMeshIndex = 0;

    std::shared_ptr< moris::hmr::Mesh > tMesh = aHMR.create_mesh( tLagrangeMeshIndex );

    aFields.resize(2);

    // create field
    aFields(0) = tMesh->create_field( "Geom", tLagrangeMeshIndex );
    aFields(1) = tMesh->create_field( "Geom", tLagrangeMeshIndex );

    aFields(0)->evaluate_scalar_function( Circle4MatMDL );
    aFields(1)->evaluate_scalar_function( Plane4MatMDL1 );

    for( uint k=0; k<0; ++k )
    {
        aHMR.flag_surface_elements_on_working_pattern( aFields(0) );
        aHMR.flag_surface_elements_on_working_pattern( aFields(1) );
        aHMR.perform_refinement_based_on_working_pattern( 0 );
        aFields(0)->evaluate_scalar_function( Circle4MatMDL );
        aFields(1)->evaluate_scalar_function( Plane4MatMDL1 );
    }

    aHMR.finalize();

    std::shared_ptr< hmr::Interpolation_Mesh_HMR > tInterpMesh = aHMR.create_interpolation_mesh( tLagrangeMeshIndex  );

}


moris::real
MultiMat3dPlane( const moris::Matrix< moris::DDRMat > & aPoint )
{

    real mXn = 1.0;
    real mYn = 0.0;
    real mZn = 0.0;
    real mXc = 0.1;
    real mYc = 0.1;
    real mZc = 0.1;

    return mXn*(aPoint(0)-mXc) + mYn*(aPoint(1)-mYc) + mZn*(aPoint(2)-mZc);
}


moris::real
MultiMat3dCyl(const moris::Matrix< moris::DDRMat > & aPoint )
{
    moris::Matrix<moris::DDRMat> aCenter = {{0.01},{0.01},{0.0}};
    moris::Matrix<moris::DDRMat> aAxis   = {{0.0},{1.0},{0.0}};
    moris::real aRad = 0.47334;
    moris::real aLength = 10;

    MORIS_ASSERT(aCenter.numel() == 3,"Centers need to have length 3");
    MORIS_ASSERT(aAxis.numel() == 3, "axis need to have length 3");

    Cell<moris::real> relativePosition = {(aPoint(0) - aCenter(0)),(aPoint(1) - aCenter(1)),(aPoint(2) - aCenter(2))};
    moris::real lsFromLeft = (relativePosition(0)*(-aAxis(0)) + relativePosition(1)*(-aAxis(1))+ relativePosition(2)*(-aAxis(2))) - aLength/2.0;
    moris::real lsFromRight = (relativePosition(0)*(aAxis(0)) + relativePosition(1)*(aAxis(1))+ relativePosition(2)*(aAxis(2))) - aLength/2.0;

    moris::real axialCrd = (relativePosition(0)*(aAxis(0)) + relativePosition(1)*(aAxis(1))+ relativePosition(2)*(aAxis(2)));
    Cell<moris::real> radDir = {(relativePosition(0) - aAxis(0)*axialCrd), (relativePosition(1) - aAxis(1)*axialCrd),(relativePosition(2) - aAxis(2)*axialCrd)};
    moris::real radDist = std::pow(radDir(0)*radDir(0)+radDir(1)*radDir(1)+radDir(2)*radDir(2), 0.5);
    moris::real lsFromRad = radDist - aRad;

    return std::max(std::max(lsFromLeft, lsFromRight), lsFromRad);
}

void
run_hmr_for_multi_mat_model_3d(hmr::HMR  &                    aHMR,
                               Cell<std::shared_ptr< moris::hmr::Field >> & aFields)
{
    moris_index tLagrangeMeshIndex = 0;

    std::shared_ptr< moris::hmr::Mesh > tMesh = aHMR.create_mesh( tLagrangeMeshIndex );

    aFields.resize(2);

    // create field
    aFields(0) = tMesh->create_field( "Geom", tLagrangeMeshIndex );

    aFields(0)->evaluate_scalar_function( MultiMat3dCyl );

    for( uint k=0; k<1; ++k )
    {
        aHMR.flag_surface_elements_on_working_pattern( aFields(0) );
        aHMR.perform_refinement_based_on_working_pattern( 0 );
        aFields(0)->evaluate_scalar_function( MultiMat3dCyl );
    }

    aFields(1) = tMesh->create_field( "Geom", tLagrangeMeshIndex );

    aFields(1)->evaluate_scalar_function( MultiMat3dPlane );
    for( uint k=0; k<1; ++k )
    {
        aHMR.flag_surface_elements_on_working_pattern( aFields(1) );
        aHMR.perform_refinement_based_on_working_pattern( 0 );
        aFields(1)->evaluate_scalar_function( MultiMat3dPlane );
    }

    aFields(0)->evaluate_scalar_function( MultiMat3dCyl );
    aFields(1)->evaluate_scalar_function( MultiMat3dPlane );
    aHMR.finalize();


    std::shared_ptr< hmr::Interpolation_Mesh_HMR > tInterpMesh = aHMR.create_interpolation_mesh( tLagrangeMeshIndex  );

}

TEST_CASE("XTK HMR 4 Material Bar Intersected By Plane and Hole","[XTK_HMR_PLANE_BAR_HOLE_2D]")
{


    if(par_size() == 1)
    {
        std::string tFieldName = "Geometry";

         moris::uint tLagrangeOrder = 1;
         moris::uint tBsplineOrder = 1;
         moris::uint tLagrangeMeshIndex = 0;
         moris::uint tBSplineMeshIndex = 0;

         moris::hmr::Parameters tParameters;
         tParameters.set_number_of_elements_per_dimension( { {22}, {8}} );
         tParameters.set_domain_dimensions({ {6}, {2} });
         tParameters.set_domain_offset({ {-3.0}, {-1.0} });
         tParameters.set_bspline_truncation( true );
         tParameters.set_output_meshes( { {0} } );
         tParameters.set_lagrange_orders  ( { {tLagrangeOrder} });
         tParameters.set_lagrange_patterns({ {0} });
         tParameters.set_bspline_orders   ( { {tBsplineOrder} } );
         tParameters.set_bspline_patterns ( { {0} } );
         tParameters.set_side_sets({{1},{2},{3},{4} });
         tParameters.set_max_refinement_level( 2 );
         tParameters.set_union_pattern( 2 );
         tParameters.set_working_pattern( 3 );
         tParameters.set_refinement_buffer( 2 );
         tParameters.set_staircase_buffer( 2 );
         tParameters.set_lagrange_to_bspline_mesh( {{ {0} }});


         hmr::HMR tHMR(tParameters);
         Cell<std::shared_ptr< moris::hmr::Field >> tHMRFields;
         run_hmr_for_multi_mat_model_2d(tHMR, tHMRFields);


         std::shared_ptr< hmr::Interpolation_Mesh_HMR > tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex  );

         moris::ge::GEN_Geom_Field tCircleFieldAsGeom(tHMRFields(0));
         moris::ge::GEN_Geom_Field tPlaneFieldAsGeom2(tHMRFields(1));
         moris::Cell<moris::ge::GEN_Geometry*> tGeometryVector = {&tCircleFieldAsGeom,&tPlaneFieldAsGeom2};

         size_t tModelDimension = 2;
         moris::ge::GEN_Phase_Table     tPhaseTable (tGeometryVector.size(),  Phase_Table_Structure::EXP_BASE_2);
         moris::ge::GEN_Geometry_Engine tGeometryEngine(tGeometryVector,tPhaseTable,tModelDimension);
         xtk::Model           tXTKModel(tModelDimension,tInterpMesh.get(),tGeometryEngine);
         tXTKModel.mVerbose = true;

        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3};
        tXTKModel.decompose(tDecompositionMethods);

        tXTKModel.perform_basis_enrichment(EntityRank::BSPLINE_1,0);

        xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
        xtk::Enriched_Integration_Mesh   & tEnrIntegMesh = tXTKModel.get_enriched_integ_mesh();

        // place the pair in mesh manager
        mtk::Mesh_Manager tMeshManager;
        tMeshManager.register_mesh_pair(&tEnrInterpMesh, &tEnrIntegMesh);
        
        //------------------------------------------------------------------------------
        // create the properties
        std::shared_ptr< fem::Property > tPropConductivity1 = std::make_shared< fem::Property >();
        tPropConductivity1->set_parameters( { {{ 1.0 }} } );
        tPropConductivity1->set_val_function( tConstValFunction2MatMDL );

        std::shared_ptr< fem::Property > tPropConductivity2 = std::make_shared< fem::Property >();
        tPropConductivity2->set_parameters( { {{ 0.1 }} } );
        tPropConductivity2->set_val_function( tConstValFunction2MatMDL );
        
        std::shared_ptr< fem::Property > tPropDirichlet = std::make_shared< fem::Property >();
        tPropDirichlet->set_parameters( { {{ 5.0 }} } );
        tPropDirichlet->set_val_function( tConstValFunction2MatMDL );

        std::shared_ptr< fem::Property > tPropNeumann = std::make_shared< fem::Property >();
        tPropNeumann->set_parameters( { {{ 20.0 }} } );
        tPropNeumann->set_val_function( tConstValFunction2MatMDL );

        std::shared_ptr< fem::Property > tPropTempLoad1 = std::make_shared< fem::Property >();
        tPropTempLoad1->set_parameters( { {{ 0.0 }} } );
        tPropTempLoad1->set_val_function( tConstValFunction2MatMDL );

        std::shared_ptr< fem::Property > tPropTempLoad2 = std::make_shared< fem::Property >();
        tPropTempLoad2->set_parameters( { {{ 50.0 }} } );
        tPropTempLoad2->set_val_function( tConstValFunction2MatMDL );

        // define constitutive models
        fem::CM_Factory tCMFactory;

        std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIso1 = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
        tCMDiffLinIso1->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tCMDiffLinIso1->set_properties( { tPropConductivity1 } );
        tCMDiffLinIso1->set_space_dim( 2 );

        std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIso2 = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
        tCMDiffLinIso2->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tCMDiffLinIso2->set_properties( { tPropConductivity1 } );
        tCMDiffLinIso2->set_space_dim( 2 );
        
        std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIso3 = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
        tCMDiffLinIso3->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tCMDiffLinIso3->set_properties( { tPropConductivity2 } );
        tCMDiffLinIso3->set_space_dim( 2 );
        
        std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIso4 = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
        tCMDiffLinIso4->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tCMDiffLinIso4->set_properties( { tPropConductivity2 } );
        tCMDiffLinIso4->set_space_dim( 2 );
        
        std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIso5 = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
        tCMDiffLinIso5->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tCMDiffLinIso5->set_properties( { tPropConductivity1 } );
        tCMDiffLinIso5->set_space_dim( 2 );
        
        std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIso6 = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
        tCMDiffLinIso6->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tCMDiffLinIso6->set_properties( { tPropConductivity1 } );
        tCMDiffLinIso6->set_space_dim( 2 );

        // define the IWGs
        fem::IWG_Factory tIWGFactory;

        std::shared_ptr< fem::IWG > tIWGBulk1 = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
        tIWGBulk1->set_residual_dof_type( { MSI::Dof_Type::TEMP } );
        tIWGBulk1->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tIWGBulk1->set_constitutive_models( { tCMDiffLinIso1 }, mtk::Master_Slave::MASTER );
        tIWGBulk1->set_properties( { tPropTempLoad1 }, mtk::Master_Slave::MASTER );

        std::shared_ptr< fem::IWG > tIWGBulk2 = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
        tIWGBulk2->set_residual_dof_type( { MSI::Dof_Type::TEMP } );
        tIWGBulk2->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tIWGBulk2->set_constitutive_models( { tCMDiffLinIso2 }, mtk::Master_Slave::MASTER );
        tIWGBulk2->set_properties( { tPropTempLoad2 }, mtk::Master_Slave::MASTER );

        std::shared_ptr< fem::IWG > tIWGDirichlet = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_DIRICHLET );
        tIWGDirichlet->set_residual_dof_type( { MSI::Dof_Type::TEMP } );
        tIWGDirichlet->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tIWGDirichlet->set_constitutive_models( { tCMDiffLinIso1 }, mtk::Master_Slave::MASTER );
        tIWGDirichlet->set_properties( { tPropDirichlet }, mtk::Master_Slave::MASTER );

        std::shared_ptr< fem::IWG > tIWGNeumann = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_NEUMANN );
        tIWGNeumann->set_residual_dof_type( { MSI::Dof_Type::TEMP } );
        tIWGNeumann->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tIWGNeumann->set_properties( { tPropNeumann }, mtk::Master_Slave::MASTER );

        std::shared_ptr< fem::IWG > tIWGInterface1 = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_INTERFACE );
        tIWGInterface1->set_residual_dof_type( { MSI::Dof_Type::TEMP } );
        tIWGInterface1->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tIWGInterface1->set_dof_type_list( {{ MSI::Dof_Type::TEMP }},mtk::Master_Slave::SLAVE );
        tIWGInterface1->set_constitutive_models( { tCMDiffLinIso3 }, mtk::Master_Slave::MASTER );
        tIWGInterface1->set_constitutive_models( { tCMDiffLinIso4 }, mtk::Master_Slave::SLAVE );

        std::shared_ptr< fem::IWG > tIWGInterface2 = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_INTERFACE );
        tIWGInterface2->set_residual_dof_type( { MSI::Dof_Type::TEMP } );
        tIWGInterface2->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tIWGInterface2->set_dof_type_list( {{ MSI::Dof_Type::TEMP }},mtk::Master_Slave::SLAVE );
        tIWGInterface2->set_constitutive_models( { tCMDiffLinIso3 }, mtk::Master_Slave::MASTER );
        tIWGInterface2->set_constitutive_models( { tCMDiffLinIso1 }, mtk::Master_Slave::SLAVE );

        std::shared_ptr< fem::IWG > tIWGInterface3 = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_INTERFACE );
        tIWGInterface3->set_residual_dof_type( { MSI::Dof_Type::TEMP } );
        tIWGInterface3->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tIWGInterface3->set_dof_type_list( {{ MSI::Dof_Type::TEMP }},mtk::Master_Slave::SLAVE );
        tIWGInterface3->set_constitutive_models( { tCMDiffLinIso1 }, mtk::Master_Slave::MASTER );
        tIWGInterface3->set_constitutive_models( { tCMDiffLinIso2 }, mtk::Master_Slave::SLAVE );
        
        // create a list of active block-sets
        std::string tDblInterfaceSideSetName01 = tEnrIntegMesh.get_dbl_interface_side_set_name(0,1);
        std::string tDblInterfaceSideSetName02 = tEnrIntegMesh.get_dbl_interface_side_set_name(0,2);
        std::string tDblInterfaceSideSetName13 = tEnrIntegMesh.get_dbl_interface_side_set_name(1,3);
        std::string tDblInterfaceSideSetName23 = tEnrIntegMesh.get_dbl_interface_side_set_name(2,3);

        std::cout<<"tDblInterfaceSideSetName01 = "<<tDblInterfaceSideSetName01<<" | Index = "<<tEnrIntegMesh.get_double_sided_set_index(tDblInterfaceSideSetName01)<<std::endl;
        std::cout<<"tDblInterfaceSideSetName02 = "<<tDblInterfaceSideSetName02<<" | Index = "<<tEnrIntegMesh.get_double_sided_set_index(tDblInterfaceSideSetName02)<<std::endl;

        // define set info
        fem::Set_User_Info tSetBulk1;
        tSetBulk1.set_mesh_index( tEnrIntegMesh.get_block_set_index("HMR_dummy_c_p0") );
        tSetBulk1.set_set_type( fem::Element_Type::BULK );
        tSetBulk1.set_IWGs( { tIWGBulk2 } );

        fem::Set_User_Info tSetBulk2;
        tSetBulk2.set_mesh_index( tEnrIntegMesh.get_block_set_index("HMR_dummy_n_p0") );
        tSetBulk2.set_set_type( fem::Element_Type::BULK );
        tSetBulk2.set_IWGs( { tIWGBulk2 } );

        fem::Set_User_Info tSetBulk3;
        tSetBulk3.set_mesh_index( tEnrIntegMesh.get_block_set_index("HMR_dummy_c_p1") );
        tSetBulk3.set_set_type( fem::Element_Type::BULK );
        tSetBulk3.set_IWGs( { tIWGBulk2 } );

        fem::Set_User_Info tSetBulk4;
        tSetBulk4.set_mesh_index( tEnrIntegMesh.get_block_set_index("HMR_dummy_n_p1") );
        tSetBulk4.set_set_type( fem::Element_Type::BULK );
        tSetBulk4.set_IWGs( { tIWGBulk2 } );

        fem::Set_User_Info tSetBulk5;
        tSetBulk5.set_mesh_index( tEnrIntegMesh.get_block_set_index("HMR_dummy_c_p2") );
        tSetBulk5.set_set_type( fem::Element_Type::BULK );
        tSetBulk5.set_IWGs( { tIWGBulk1 } );

        fem::Set_User_Info tSetBulk6;
        tSetBulk6.set_mesh_index( tEnrIntegMesh.get_block_set_index("HMR_dummy_n_p2") );
        tSetBulk6.set_set_type( fem::Element_Type::BULK );
        tSetBulk6.set_IWGs( { tIWGBulk1 } );

        fem::Set_User_Info tSetBulk7;
        tSetBulk7.set_mesh_index( tEnrIntegMesh.get_block_set_index("HMR_dummy_c_p3") );
        tSetBulk7.set_set_type( fem::Element_Type::BULK );
        tSetBulk7.set_IWGs( { tIWGBulk1 } );

        fem::Set_User_Info tSetBulk8;
        tSetBulk8.set_mesh_index( tEnrIntegMesh.get_block_set_index("HMR_dummy_n_p3") );
        tSetBulk8.set_set_type( fem::Element_Type::BULK );
        tSetBulk8.set_IWGs( { tIWGBulk1 } );

        fem::Set_User_Info tSetDirichlet;
        tSetDirichlet.set_mesh_index( tEnrIntegMesh.get_side_set_index("SideSet_4_n_p2") );
        tSetDirichlet.set_set_type( fem::Element_Type::SIDESET );
        tSetDirichlet.set_IWGs( { tIWGDirichlet } );

        fem::Set_User_Info tSetNeumann;
        tSetNeumann.set_mesh_index( tEnrIntegMesh.get_side_set_index("SideSet_2_n_p3") );
        tSetNeumann.set_set_type( fem::Element_Type::SIDESET );
        tSetNeumann.set_IWGs( { tIWGNeumann } );

        fem::Set_User_Info tSetInterface1;
        tSetInterface1.set_mesh_index( tEnrIntegMesh.get_double_sided_set_index(tDblInterfaceSideSetName01) );
        tSetInterface1.set_set_type( fem::Element_Type::DOUBLE_SIDESET );
        tSetInterface1.set_IWGs( { tIWGInterface1 } );

        fem::Set_User_Info tSetInterface2;
        tSetInterface2.set_mesh_index( tEnrIntegMesh.get_double_sided_set_index(tDblInterfaceSideSetName02) );
        tSetInterface2.set_set_type( fem::Element_Type::DOUBLE_SIDESET );
        tSetInterface2.set_IWGs( { tIWGInterface2 } );

        fem::Set_User_Info tSetInterface3;
        tSetInterface3.set_mesh_index( tEnrIntegMesh.get_double_sided_set_index(tDblInterfaceSideSetName13) );
        tSetInterface3.set_set_type( fem::Element_Type::DOUBLE_SIDESET );
        tSetInterface3.set_IWGs( { tIWGInterface2 } );

        fem::Set_User_Info tSetInterface4;
        tSetInterface4.set_mesh_index( tEnrIntegMesh.get_double_sided_set_index(tDblInterfaceSideSetName23) );
        tSetInterface4.set_set_type( fem::Element_Type::DOUBLE_SIDESET );
        tSetInterface4.set_IWGs( { tIWGInterface3 } );

        // create a cell of set info
        moris::Cell< fem::Set_User_Info > tSetInfo( 14 );
        tSetInfo( 0 ) = tSetBulk1;
        tSetInfo( 1 ) = tSetBulk2;
        tSetInfo( 2 ) = tSetBulk3;
        tSetInfo( 3 ) = tSetBulk4;
        tSetInfo( 4 ) = tSetBulk5;
        tSetInfo( 5 ) = tSetBulk6;
        tSetInfo( 6 ) = tSetBulk7;
        tSetInfo( 7 ) = tSetBulk8;
        tSetInfo( 8 ) = tSetDirichlet;
        tSetInfo( 9 ) = tSetNeumann;
        tSetInfo( 10 ) = tSetInterface1;
        tSetInfo( 11 ) = tSetInterface2;
        tSetInfo( 12 ) = tSetInterface3;
        tSetInfo( 13 ) = tSetInterface4;

        // create model
        mdl::Model * tModel = new mdl::Model( &tMeshManager,
                                               0,
                                               tSetInfo,
                                               0,
                                               false );

        moris::Cell< enum MSI::Dof_Type > tDofTypes1( 1, MSI::Dof_Type::TEMP );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 1: create linear solver and algorithm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        dla::Solver_Factory  tSolFactory;
        std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( SolverType::AZTEC_IMPL );

        tLinearSolverAlgorithm->set_param("AZ_diagnostics") = AZ_none;
        tLinearSolverAlgorithm->set_param("AZ_output") = AZ_none;
        tLinearSolverAlgorithm->set_param("AZ_orthog") = 1;
        tLinearSolverAlgorithm->set_param("AZ_solver") = AZ_gmres_condnum;
        tLinearSolverAlgorithm->set_param("AZ_precond") = AZ_dom_decomp;
        tLinearSolverAlgorithm->set_param("AZ_ilut_fill") = 10.0;
        tLinearSolverAlgorithm->set_param("AZ_max_iter") = 100;
        tLinearSolverAlgorithm->set_param("rel_residual") = 1e-6;


        dla::Linear_Solver tLinSolver;

        tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 2: create nonlinear solver and algorithm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        NLA::Nonlinear_Solver_Factory tNonlinFactory;
        std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );

        tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );

        NLA::Nonlinear_Solver tNonlinearSolver;
        tNonlinearSolverAlgorithm->set_param("NLA_max_iter")   = 2;
        tNonlinearSolverAlgorithm->set_param("NLA_hard_break") = false;
        tNonlinearSolverAlgorithm->set_param("NLA_max_lin_solver_restarts") = 2;
        tNonlinearSolverAlgorithm->set_param("NLA_rebuild_jacobian") = true;
        tNonlinearSolverAlgorithm->set_param("NLA_rel_residual") = 1e-6;

        tNonlinearSolver.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 3: create time Solver and algorithm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        tsa::Time_Solver_Factory tTimeSolverFactory;
        std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm = tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );

        tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolver );

        tsa::Time_Solver tTimeSolver;

        tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );

        NLA::SOL_Warehouse tSolverWarehouse;

        tSolverWarehouse.set_solver_interface(tModel->get_solver_interface());

        tNonlinearSolver.set_solver_warehouse( &tSolverWarehouse );
        tTimeSolver.set_solver_warehouse( &tSolverWarehouse );

        tNonlinearSolver.set_dof_type_list( tDofTypes1 );
        tTimeSolver.set_dof_type_list( tDofTypes1 );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 4: Solve and check
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        tTimeSolver.solve();

//
//        // verify solution
//        moris::Matrix<DDRMat> tGoldSolution =
//        {{ 5.000000e+00},
//         { 2.500000e+01},
//         { 4.500000e+01},
//         { 6.500000e+01},
//         { 5.000000e+00},
//         { 2.500000e+01},
//         { 4.500000e+01},
//         { 6.500000e+01}};
//
        Matrix<DDRMat> tFullSol;
        tTimeSolver.get_full_solution(tFullSol);


        // Declare the fields related to enrichment strategy in output options
        // output solution and meshes
        xtk::Output_Options tOutputOptions;
        tOutputOptions.mAddNodeSets = false;
        tOutputOptions.mAddSideSets = true;
        tOutputOptions.mAddClusters = false;

        // add solution field to integration mesh
        std::string tIntegSolFieldName = "solution";
        tOutputOptions.mRealNodeExternalFieldNames = {tIntegSolFieldName};

        moris::mtk::Integration_Mesh* tIntegMesh1 = tXTKModel.get_output_mesh(tOutputOptions);

        // Write to Integration mesh for visualization
        Matrix<DDRMat> tIntegSol = tModel->get_solution_for_integration_mesh_output( MSI::Dof_Type::TEMP );


        Matrix<DDRMat> tSTKIntegSol(tIntegMesh1->get_num_entities(EntityRank::NODE),1);

        for(moris::uint i = 0; i < tIntegMesh1->get_num_entities(EntityRank::NODE); i++)
        {
            moris::moris_id tID = tIntegMesh1->get_glb_entity_id_from_entity_loc_index(i,EntityRank::NODE);
            moris::moris_index tMyIndex = tEnrIntegMesh.get_loc_entity_ind_from_entity_glb_id(tID,EntityRank::NODE);

            tSTKIntegSol(i) = tIntegSol(tMyIndex);
        }

        // crate field in integration mesh
        moris::moris_index tFieldIndex = tEnrIntegMesh.create_field("Solution",EntityRank::NODE);
        tEnrIntegMesh.add_field_data(tFieldIndex,EntityRank::NODE,tSTKIntegSol);

        // add solution field to integration mesh
        tIntegMesh1->add_mesh_field_real_scalar_data_loc_inds(tIntegSolFieldName,EntityRank::NODE,tSTKIntegSol);


        std::string tMeshOutputFile = "./mdl_exo/xtk_hmr_bar_plane_hole_bl_2_mat_l" + std::to_string(tLagrangeOrder) + "_b"+std::to_string(tBsplineOrder)+".e";
        tIntegMesh1->create_output_mesh(tMeshOutputFile);

        delete tModel;
        delete tIntegMesh1;
    }
}


TEST_CASE("XTK HMR 4 Material Bar Intersected By Plane and Hole 3D","[XTK_HMR_PLANE_BAR_HOLE_3D]")
{


    if(par_size() == 1)
    {
        std::string tFieldName = "Geometry";

        moris::uint tLagrangeOrder = 1;
        moris::uint tBsplineOrder = 1;
         moris::uint tLagrangeMeshIndex = 0;
         moris::uint tBSplineMeshIndex = 0;

         moris::hmr::Parameters tParameters;
         tParameters.set_number_of_elements_per_dimension( { {22}, {8},{8}} );
         tParameters.set_domain_dimensions({ {6}, {2}, {2} });
         tParameters.set_domain_offset({ {-3.0}, {-1.0},{-1} });
         tParameters.set_bspline_truncation( true );
         tParameters.set_output_meshes( { {0} } );
         tParameters.set_lagrange_orders  ( { {tLagrangeOrder} });
         tParameters.set_lagrange_patterns({ {0} });
         tParameters.set_bspline_orders   ( { {tBsplineOrder} } );
         tParameters.set_bspline_patterns ( { {0} } );
         tParameters.set_side_sets({{1},{2},{3},{4},{5},{6}});
         tParameters.set_max_refinement_level( 2 );
         tParameters.set_union_pattern( 2 );
         tParameters.set_working_pattern( 3 );
         tParameters.set_refinement_buffer( 2 );
         tParameters.set_staircase_buffer( 2 );
         tParameters.set_lagrange_to_bspline_mesh( {{ {0} }});


         hmr::HMR tHMR(tParameters);
         Cell<std::shared_ptr< moris::hmr::Field >> tHMRFields;
         run_hmr_for_multi_mat_model_3d(tHMR, tHMRFields);


         std::shared_ptr< hmr::Interpolation_Mesh_HMR > tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex  );

         moris::ge::GEN_Geom_Field tCircleFieldAsGeom(tHMRFields(0));
         moris::ge::GEN_Geom_Field tPlaneFieldAsGeom2(tHMRFields(1));
         moris::Cell<moris::ge::GEN_Geometry*> tGeometryVector = {&tCircleFieldAsGeom,&tPlaneFieldAsGeom2};

         size_t tModelDimension = 3;
         moris::ge::GEN_Phase_Table     tPhaseTable (tGeometryVector.size(),  Phase_Table_Structure::EXP_BASE_2);
         moris::ge::GEN_Geometry_Engine tGeometryEngine(tGeometryVector,tPhaseTable,tModelDimension);
         xtk::Model           tXTKModel(tModelDimension,tInterpMesh.get(),tGeometryEngine);
         tXTKModel.mVerbose = true;

        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4};
        tXTKModel.decompose(tDecompositionMethods);

        tXTKModel.perform_basis_enrichment(EntityRank::BSPLINE_1,0);

        xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
        xtk::Enriched_Integration_Mesh   & tEnrIntegMesh = tXTKModel.get_enriched_integ_mesh();

        tEnrIntegMesh.print_side_sets();

        // place the pair in mesh manager
        mtk::Mesh_Manager tMeshManager;
        tMeshManager.register_mesh_pair(&tEnrInterpMesh, &tEnrIntegMesh);

        //------------------------------------------------------------------------------
        // create the properties
        std::shared_ptr< fem::Property > tPropConductivity1 = std::make_shared< fem::Property >();
        tPropConductivity1->set_parameters( { {{ 1.0 }} } );
        tPropConductivity1->set_val_function( tConstValFunction2MatMDL );

        std::shared_ptr< fem::Property > tPropConductivity2 = std::make_shared< fem::Property >();
        tPropConductivity2->set_parameters( { {{ 0.1 }} } );
        tPropConductivity2->set_val_function( tConstValFunction2MatMDL );

        std::shared_ptr< fem::Property > tPropDirichlet = std::make_shared< fem::Property >();
        tPropDirichlet->set_parameters( { {{ 5.0 }} } );
        tPropDirichlet->set_val_function( tConstValFunction2MatMDL );

        std::shared_ptr< fem::Property > tPropNeumann = std::make_shared< fem::Property >();
        tPropNeumann->set_parameters( { {{ 20.0 }} } );
        tPropNeumann->set_val_function( tConstValFunction2MatMDL );

        std::shared_ptr< fem::Property > tPropTempLoad1 = std::make_shared< fem::Property >();
        tPropTempLoad1->set_parameters( { {{ 0.0 }} } );
        tPropTempLoad1->set_val_function( tConstValFunction2MatMDL );

        std::shared_ptr< fem::Property > tPropTempLoad2 = std::make_shared< fem::Property >();
        tPropTempLoad2->set_parameters( { {{ 50.0 }} } );
        tPropTempLoad2->set_val_function( tConstValFunction2MatMDL );

        // define constitutive models
        fem::CM_Factory tCMFactory;

        std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIso1 = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
        tCMDiffLinIso1->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tCMDiffLinIso1->set_properties( { tPropConductivity1 } );
        tCMDiffLinIso1->set_space_dim( 3 );

        std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIso2 = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
        tCMDiffLinIso2->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tCMDiffLinIso2->set_properties( { tPropConductivity2 } );
        tCMDiffLinIso2->set_space_dim( 3 );

        // define the IWGs
        fem::IWG_Factory tIWGFactory;

        std::shared_ptr< fem::IWG > tIWGBulk1 = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
        tIWGBulk1->set_residual_dof_type( { MSI::Dof_Type::TEMP } );
        tIWGBulk1->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tIWGBulk1->set_constitutive_models( { tCMDiffLinIso1 }, mtk::Master_Slave::MASTER );
        tIWGBulk1->set_properties( { tPropTempLoad1 }, mtk::Master_Slave::MASTER );

        std::shared_ptr< fem::IWG > tIWGBulk2 = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
        tIWGBulk2->set_residual_dof_type( { MSI::Dof_Type::TEMP } );
        tIWGBulk2->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tIWGBulk2->set_constitutive_models( { tCMDiffLinIso2 }, mtk::Master_Slave::MASTER );
        tIWGBulk2->set_properties( { tPropTempLoad2 }, mtk::Master_Slave::MASTER );

        std::shared_ptr< fem::IWG > tIWGDirichlet = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_DIRICHLET );
        tIWGDirichlet->set_residual_dof_type( { MSI::Dof_Type::TEMP } );
        tIWGDirichlet->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tIWGDirichlet->set_constitutive_models( { tCMDiffLinIso2 }, mtk::Master_Slave::MASTER );
        tIWGDirichlet->set_properties( { tPropDirichlet }, mtk::Master_Slave::MASTER );

        std::shared_ptr< fem::IWG > tIWGNeumann = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_NEUMANN );
        tIWGNeumann->set_residual_dof_type( { MSI::Dof_Type::TEMP } );
        tIWGNeumann->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tIWGNeumann->set_properties( { tPropNeumann }, mtk::Master_Slave::MASTER );

        std::shared_ptr< fem::IWG > tIWGInterface1 = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_INTERFACE );
        tIWGInterface1->set_residual_dof_type( { MSI::Dof_Type::TEMP } );
        tIWGInterface1->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tIWGInterface1->set_dof_type_list( {{ MSI::Dof_Type::TEMP }},mtk::Master_Slave::SLAVE );
        tIWGInterface1->set_constitutive_models( { tCMDiffLinIso2 }, mtk::Master_Slave::MASTER );
        tIWGInterface1->set_constitutive_models( { tCMDiffLinIso2 }, mtk::Master_Slave::SLAVE );

        std::shared_ptr< fem::IWG > tIWGInterface2 = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_INTERFACE );
        tIWGInterface2->set_residual_dof_type( { MSI::Dof_Type::TEMP } );
        tIWGInterface2->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tIWGInterface2->set_dof_type_list( {{ MSI::Dof_Type::TEMP }},mtk::Master_Slave::SLAVE );
        tIWGInterface2->set_constitutive_models( { tCMDiffLinIso2 }, mtk::Master_Slave::MASTER );
        tIWGInterface2->set_constitutive_models( { tCMDiffLinIso1 }, mtk::Master_Slave::SLAVE );

        std::shared_ptr< fem::IWG > tIWGInterface3 = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_INTERFACE );
        tIWGInterface3->set_residual_dof_type( { MSI::Dof_Type::TEMP } );
        tIWGInterface3->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tIWGInterface3->set_dof_type_list( {{ MSI::Dof_Type::TEMP }},mtk::Master_Slave::SLAVE );
        tIWGInterface3->set_constitutive_models( { tCMDiffLinIso1 }, mtk::Master_Slave::MASTER );
        tIWGInterface3->set_constitutive_models( { tCMDiffLinIso1 }, mtk::Master_Slave::SLAVE );


        // create a list of active block-sets
        std::string tDblInterfaceSideSetName01 = tEnrIntegMesh.get_dbl_interface_side_set_name(0,1);
        std::string tDblInterfaceSideSetName02 = tEnrIntegMesh.get_dbl_interface_side_set_name(0,2);
        std::string tDblInterfaceSideSetName13 = tEnrIntegMesh.get_dbl_interface_side_set_name(1,3);
        std::string tDblInterfaceSideSetName23 = tEnrIntegMesh.get_dbl_interface_side_set_name(2,3);

        std::cout<<"tDblInterfaceSideSetName01 = "<<tDblInterfaceSideSetName01<<" | Index = "<<tEnrIntegMesh.get_double_sided_set_index(tDblInterfaceSideSetName01)<<std::endl;
        std::cout<<"tDblInterfaceSideSetName02 = "<<tDblInterfaceSideSetName02<<" | Index = "<<tEnrIntegMesh.get_double_sided_set_index(tDblInterfaceSideSetName02)<<std::endl;

        // define set info
        fem::Set_User_Info tSetBulk1;
        tSetBulk1.set_mesh_index( tEnrIntegMesh.get_block_set_index("HMR_dummy_c_p0") );
        tSetBulk1.set_set_type( fem::Element_Type::BULK );
        tSetBulk1.set_IWGs( { tIWGBulk2 } );

        fem::Set_User_Info tSetBulk2;
        tSetBulk2.set_mesh_index( tEnrIntegMesh.get_block_set_index("HMR_dummy_n_p0") );
        tSetBulk2.set_set_type( fem::Element_Type::BULK );
        tSetBulk2.set_IWGs( { tIWGBulk2 } );

        fem::Set_User_Info tSetBulk3;
        tSetBulk3.set_mesh_index( tEnrIntegMesh.get_block_set_index("HMR_dummy_c_p1") );
        tSetBulk3.set_set_type( fem::Element_Type::BULK );
        tSetBulk3.set_IWGs( { tIWGBulk2 } );

        fem::Set_User_Info tSetBulk4;
        tSetBulk4.set_mesh_index( tEnrIntegMesh.get_block_set_index("HMR_dummy_n_p1") );
        tSetBulk4.set_set_type( fem::Element_Type::BULK );
        tSetBulk4.set_IWGs( { tIWGBulk2 } );

        fem::Set_User_Info tSetBulk5;
        tSetBulk5.set_mesh_index( tEnrIntegMesh.get_block_set_index("HMR_dummy_c_p2") );
        tSetBulk5.set_set_type( fem::Element_Type::BULK );
        tSetBulk5.set_IWGs( { tIWGBulk1 } );

        fem::Set_User_Info tSetBulk6;
        tSetBulk6.set_mesh_index( tEnrIntegMesh.get_block_set_index("HMR_dummy_n_p2") );
        tSetBulk6.set_set_type( fem::Element_Type::BULK );
        tSetBulk6.set_IWGs( { tIWGBulk1 } );

        fem::Set_User_Info tSetBulk7;
        tSetBulk7.set_mesh_index( tEnrIntegMesh.get_block_set_index("HMR_dummy_c_p3") );
        tSetBulk7.set_set_type( fem::Element_Type::BULK );
        tSetBulk7.set_IWGs( { tIWGBulk1 } );

        fem::Set_User_Info tSetBulk8;
        tSetBulk8.set_mesh_index( tEnrIntegMesh.get_block_set_index("HMR_dummy_n_p3") );
        tSetBulk8.set_set_type( fem::Element_Type::BULK );
        tSetBulk8.set_IWGs( { tIWGBulk1 } );

        fem::Set_User_Info tSetDirichlet;
        tSetDirichlet.set_mesh_index( tEnrIntegMesh.get_side_set_index("SideSet_4_n_p2") );
        tSetDirichlet.set_set_type( fem::Element_Type::SIDESET );
        tSetDirichlet.set_IWGs( { tIWGDirichlet } );

        fem::Set_User_Info tSetNeumann;
        tSetNeumann.set_mesh_index( tEnrIntegMesh.get_side_set_index("SideSet_2_n_p3") );
        tSetNeumann.set_set_type( fem::Element_Type::SIDESET );
        tSetNeumann.set_IWGs( { tIWGNeumann } );

        fem::Set_User_Info tSetInterface1;
        tSetInterface1.set_mesh_index( tEnrIntegMesh.get_double_sided_set_index(tDblInterfaceSideSetName01) );
        tSetInterface1.set_set_type( fem::Element_Type::DOUBLE_SIDESET );
        tSetInterface1.set_IWGs( { tIWGInterface1 } );

        fem::Set_User_Info tSetInterface2;
        tSetInterface2.set_mesh_index( tEnrIntegMesh.get_double_sided_set_index(tDblInterfaceSideSetName02) );
        tSetInterface2.set_set_type( fem::Element_Type::DOUBLE_SIDESET );
        tSetInterface2.set_IWGs( { tIWGInterface2 } );

        fem::Set_User_Info tSetInterface3;
        tSetInterface3.set_mesh_index( tEnrIntegMesh.get_double_sided_set_index(tDblInterfaceSideSetName13) );
        tSetInterface3.set_set_type( fem::Element_Type::DOUBLE_SIDESET );
        tSetInterface3.set_IWGs( { tIWGInterface2 } );

        fem::Set_User_Info tSetInterface4;
        tSetInterface4.set_mesh_index( tEnrIntegMesh.get_double_sided_set_index(tDblInterfaceSideSetName23) );
        tSetInterface4.set_set_type( fem::Element_Type::DOUBLE_SIDESET );
        tSetInterface4.set_IWGs( { tIWGInterface3 } );

        // create a cell of set info
        moris::Cell< fem::Set_User_Info > tSetInfo( 14 );
        tSetInfo( 0 ) = tSetBulk1;
        tSetInfo( 1 ) = tSetBulk2;
        tSetInfo( 2 ) = tSetBulk3;
        tSetInfo( 3 ) = tSetBulk4;
        tSetInfo( 4 ) = tSetBulk5;
        tSetInfo( 5 ) = tSetBulk6;
        tSetInfo( 6 ) = tSetBulk7;
        tSetInfo( 7 ) = tSetBulk8;
        tSetInfo( 8 ) = tSetDirichlet;
        tSetInfo( 9 ) = tSetNeumann;
        tSetInfo( 10 ) = tSetInterface1;
        tSetInfo( 11 ) = tSetInterface2;
        tSetInfo( 12 ) = tSetInterface3;
        tSetInfo( 13 ) = tSetInterface4;

        // create model
        mdl::Model * tModel = new mdl::Model( &tMeshManager,
                                               1,
                                               tSetInfo,
                                               0, false );

        moris::Cell< enum MSI::Dof_Type > tDofTypes1( 1, MSI::Dof_Type::TEMP );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 1: create linear solver and algorithm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


        dla::Solver_Factory  tSolFactory;
        std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( SolverType::AZTEC_IMPL );

        tLinearSolverAlgorithm->set_param("AZ_diagnostics") = AZ_none;
        tLinearSolverAlgorithm->set_param("AZ_output") = AZ_none;

        dla::Linear_Solver tLinSolver;

        tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 2: create nonlinear solver and algorithm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        NLA::Nonlinear_Solver_Factory tNonlinFactory;
        std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );

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

        NLA::SOL_Warehouse tSolverWarehouse;

        tSolverWarehouse.set_solver_interface(tModel->get_solver_interface());

        tNonlinearSolver.set_solver_warehouse( &tSolverWarehouse );
        tTimeSolver.set_solver_warehouse( &tSolverWarehouse );

        tNonlinearSolver.set_dof_type_list( tDofTypes1 );
        tTimeSolver.set_dof_type_list( tDofTypes1 );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 4: Solve and check
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        tTimeSolver.solve();

//
//        // verify solution
//        moris::Matrix<DDRMat> tGoldSolution =
//        {{ 5.000000e+00},
//         { 2.500000e+01},
//         { 4.500000e+01},
//         { 6.500000e+01},
//         { 5.000000e+00},
//         { 2.500000e+01},
//         { 4.500000e+01},
//         { 6.500000e+01}};
//
        Matrix<DDRMat> tFullSol;
        tTimeSolver.get_full_solution(tFullSol);

        print(tFullSol,"tFullSol");


        // Declare the fields related to enrichment strategy in output options
        // output solution and meshes
        xtk::Output_Options tOutputOptions;
        tOutputOptions.mAddNodeSets = false;
        tOutputOptions.mAddSideSets = true;
        tOutputOptions.mAddClusters = false;

        // add solution field to integration mesh
        std::string tIntegSolFieldName = "solution";
        tOutputOptions.mRealNodeExternalFieldNames = {tIntegSolFieldName};

        moris::mtk::Integration_Mesh* tIntegMesh1 = tXTKModel.get_output_mesh(tOutputOptions);

        // Write to Integration mesh for visualization
        Matrix<DDRMat> tIntegSol = tModel->get_solution_for_integration_mesh_output( MSI::Dof_Type::TEMP );


        Matrix<DDRMat> tSTKIntegSol(tIntegMesh1->get_num_entities(EntityRank::NODE),1);

        for(moris::uint i = 0; i < tIntegMesh1->get_num_entities(EntityRank::NODE); i++)
        {
            moris::moris_id tID = tIntegMesh1->get_glb_entity_id_from_entity_loc_index(i,EntityRank::NODE);
            moris::moris_index tMyIndex = tEnrIntegMesh.get_loc_entity_ind_from_entity_glb_id(tID,EntityRank::NODE);

            tSTKIntegSol(i) = tIntegSol(tMyIndex);
        }

        // crate field in integration mesh
        moris::moris_index tFieldIndex = tEnrIntegMesh.create_field("Solution",EntityRank::NODE);
        tEnrIntegMesh.add_field_data(tFieldIndex,EntityRank::NODE,tSTKIntegSol);

        // add solution field to integration mesh
        tIntegMesh1->add_mesh_field_real_scalar_data_loc_inds(tIntegSolFieldName,EntityRank::NODE,tSTKIntegSol);


        std::string tMeshOutputFile = "./mdl_exo/xtk_hmr_bar_plane_hole_3d_l" + std::to_string(tLagrangeOrder) + "_b"+std::to_string(tBsplineOrder)+".e";
        tIntegMesh1->create_output_mesh(tMeshOutputFile);

        delete tModel;
        delete tIntegMesh1;
    }
}

