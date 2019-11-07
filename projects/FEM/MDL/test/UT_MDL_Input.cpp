/*
 * UT_MDL_Input.cpp
 *
 *  Created on: Oct 23, 2019
 *      Author: noel
 */

#include "catch.hpp"
#include "cl_Star.hpp"
#include "cl_Circle.hpp"
#include "cl_Plane.hpp"

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

#include "cl_FEM_IWG_Factory.hpp"             //FEM/INT/src
#include "cl_FEM_CM_Factory.hpp"              //FEM/INT/src
#include "cl_FEM_Set_User_Info.hpp"           //FEM/INT/src

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

#include "../projects/GEN/src/ripped/geometry/cl_GEN_Geom_Field.hpp"

#include "fn_norm.hpp"


namespace moris
{

// define free function for properties
 Matrix< DDRMat > tPropValConstFunc( moris::Cell< Matrix< DDRMat > >         & aParameters,
                                     moris::Cell< fem::Field_Interpolator* > & aDofFieldInterpolator,
                                     moris::Cell< fem::Field_Interpolator* > & aDvFieldInterpolator,
                                     fem::Geometry_Interpolator              * aGeometryInterpolator )
 {
     return aParameters( 0 );
 }
 Matrix< DDRMat > tPropValFunc( moris::Cell< Matrix< DDRMat > >         & aParameters,
                                moris::Cell< fem::Field_Interpolator* > & aDofFieldInterpolator,
                                moris::Cell< fem::Field_Interpolator* > & aDvFieldInterpolator,
                                fem::Geometry_Interpolator              * aGeometryInterpolator )
 {
     return aParameters( 0 ) + aParameters( 1 ) * aDofFieldInterpolator( 0 )->val();
 }
 Matrix< DDRMat > tPropDerFunc( moris::Cell< Matrix< DDRMat > >         & aParameters,
                                moris::Cell< fem::Field_Interpolator* > & aDofFieldInterpolator,
                                moris::Cell< fem::Field_Interpolator* > & aDvFieldInterpolator,
                                fem::Geometry_Interpolator              * aGeometryInterpolator )
 {
     return aParameters( 1 ) * aDofFieldInterpolator( 0 )->N();
 }

 moris::real LevelSetPlaneFunction( const moris::Matrix< moris::DDRMat > & aPoint )
 {
     real mXn = 0;
     real mYn = 0;
     real mZn = 1.0;
     real mXc = 1.0;
     real mYc = 1.0;
     real mZc = 1.4;
     return mXn*(aPoint(0)-mXc) + mYn*(aPoint(1)-mYc) + mZn*(aPoint(2)-mZc);
 }

TEST_CASE("MDL Input","[MDL_Input]")
{
    if(par_size() == 1)
    {
        std::string tFieldName = "Circle";


        moris::uint tLagrangeMeshIndex = 0;
        moris::uint tBSplineMeshIndex = 0;

        moris::hmr::Parameters tParameters;

        tParameters.set_number_of_elements_per_dimension( { {1}, {1}, {4} } );
        tParameters.set_domain_dimensions({ {1}, {1}, {2} });
        tParameters.set_domain_offset({ {0.0}, {0.0}, {0.0} });
        tParameters.set_bspline_truncation( true );
        tParameters.set_side_sets({ {5}, {6} });

        tParameters.set_output_meshes( { {0} } );

        tParameters.set_lagrange_orders  ( { {1} });
        tParameters.set_lagrange_patterns({ {0} });

        tParameters.set_bspline_orders   ( { {1} } );
        tParameters.set_bspline_patterns ( { {0} } );

        tParameters.set_union_pattern( 2 );
        tParameters.set_working_pattern( 3 );

        tParameters.set_refinement_buffer( 2 );
        tParameters.set_staircase_buffer( 2);

        Cell< Matrix< DDUMat > > tLagrangeToBSplineMesh( 1 );
        tLagrangeToBSplineMesh( 0 ) = { {0} };

        tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

        hmr::HMR tHMR( tParameters );

        std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

        // create field
        std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( tFieldName, tLagrangeMeshIndex );

        tField->evaluate_scalar_function( LevelSetPlaneFunction );

        for( uint k=0; k<0; ++k )
        {
            tHMR.flag_surface_elements_on_working_pattern( tField );
            tHMR.perform_refinement_based_on_working_pattern( 0 );

            tField->evaluate_scalar_function( LevelSetPlaneFunction );
        }

        tHMR.finalize();

        tHMR.save_to_exodus( 0, "./mdl_exo/mdl_input.e" );

        std::shared_ptr< hmr::Interpolation_Mesh_HMR > tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex  );

        moris::ge::GEN_Geom_Field tFieldAsGeom(tField);

        moris::Cell<ge::GEN_Geometry*> tGeometryVector = {&tFieldAsGeom};

        // Tell the geometry engine about the discrete field mesh and how to interpret phases
        ge::GEN_Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
        ge::GEN_Geometry_Engine tGeometryEngine(tGeometryVector,tPhaseTable);

        // Tell the XTK model that it should decompose with a C_HIERARCHY_TET4, on the same mesh that the level set field is defined on.
        size_t tModelDimension = 3;
        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,Subdivision_Method::C_HIERARCHY_TET4};
        xtk::Model tXTKModel(tModelDimension,tInterpMesh.get(), tGeometryEngine);
        tXTKModel.mSameMesh = true;
        tXTKModel.mVerbose = true;

        // Do the cutting
        tXTKModel.decompose(tDecompositionMethods);

        xtk::Output_Options tOutputOptions;
        tOutputOptions.mAddNodeSets = false;
        tOutputOptions.mAddSideSets = true;
        tOutputOptions.mAddClusters = true;

        // add solution field to integration mesh
        std::string tIntegSolFieldName = "solution";
        tOutputOptions.mRealNodeExternalFieldNames = {tIntegSolFieldName};

        moris::mtk::Integration_Mesh* tIntegMesh1 = tXTKModel.get_output_mesh(tOutputOptions);

        // place the pair in mesh manager
        mtk::Mesh_Manager tMeshManager;
        tMeshManager.register_mesh_pair(tInterpMesh.get(), tIntegMesh1);

        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------
        // create the properties
        std::shared_ptr< fem::Property > tPropConductivity = std::make_shared< fem::Property >();
        tPropConductivity->set_parameters( { {{ 1.0 }}, {{ 1.0 }} } );
        tPropConductivity->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tPropConductivity->set_val_function( tPropValFunc );
        tPropConductivity->set_dof_derivative_functions( { tPropDerFunc } );

        std::shared_ptr< fem::Property > tPropDirichlet = std::make_shared< fem::Property >();
        tPropDirichlet->set_parameters( { {{ 5.0 }} } );
        tPropDirichlet->set_val_function( tPropValConstFunc );

        std::shared_ptr< fem::Property > tPropNeumann = std::make_shared< fem::Property >();
        tPropNeumann->set_parameters( { {{ 20.0 }} } );
        tPropNeumann->set_val_function( tPropValConstFunc );

        // define constitutive models
        fem::CM_Factory tCMFactory;

        std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIso = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
        tCMDiffLinIso->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} ); // FIXME through the factory?
        tCMDiffLinIso->set_properties( { tPropConductivity } );

        // define the IWGs
        fem::IWG_Factory tIWGFactory;

        std::shared_ptr< fem::IWG > tIWGBulk = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
        tIWGBulk->set_residual_dof_type( { MSI::Dof_Type::TEMP } );                          // FIXME through the factory?
        tIWGBulk->set_dof_type_list( {{ MSI::Dof_Type::TEMP }}, mtk::Master_Slave::MASTER ); // FIXME through the factory?
        tIWGBulk->set_constitutive_models( { tCMDiffLinIso }, mtk::Master_Slave::MASTER );

        std::shared_ptr< fem::IWG > tIWGDirichlet = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_DIRICHLET );
        tIWGDirichlet->set_residual_dof_type( { MSI::Dof_Type::TEMP } );                          // FIXME through the factory?
        tIWGDirichlet->set_dof_type_list( {{ MSI::Dof_Type::TEMP }}, mtk::Master_Slave::MASTER ); // FIXME through the factory?
        tIWGDirichlet->set_properties( { tPropDirichlet }, mtk::Master_Slave::MASTER );

        std::shared_ptr< fem::IWG > tIWGNeumann = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_NEUMANN );
        tIWGNeumann->set_residual_dof_type( { MSI::Dof_Type::TEMP } );                          // FIXME through the factory?
        tIWGNeumann->set_dof_type_list( {{ MSI::Dof_Type::TEMP }}, mtk::Master_Slave::MASTER ); // FIXME through the factory?
        tIWGNeumann->set_properties( { tPropNeumann }, mtk::Master_Slave::MASTER );

        // define set info
        fem::Set_User_Info tSetBulk1;
        tSetBulk1.set_mesh_index( tIntegMesh1->get_block_set_index("child_0") ); // FIXME set index within the mesh
        tSetBulk1.set_set_type( fem::Element_Type::BULK );
        tSetBulk1.set_IWGs( { tIWGBulk } );

        fem::Set_User_Info tSetBulk2;
        tSetBulk2.set_mesh_index(  tIntegMesh1->get_block_set_index("parent_0") ); // FIXME set index within the mesh
        tSetBulk2.set_set_type( fem::Element_Type::BULK );
        tSetBulk2.set_IWGs( { tIWGBulk } );

        fem::Set_User_Info tSetNeumann;
        tSetNeumann.set_mesh_index( tIntegMesh1->get_side_set_index("iside_g_0_p0_0_p1_1") ); // FIXME set index within the mesh
        tSetNeumann.set_set_type( fem::Element_Type::SIDESET );
        tSetNeumann.set_IWGs( { tIWGNeumann } );

        fem::Set_User_Info tSetDirichlet;
        tSetDirichlet.set_mesh_index( tIntegMesh1->get_side_set_index("SideSet_1") ); // FIXME set index within the mesh
        tSetDirichlet.set_set_type( fem::Element_Type::SIDESET );
        tSetDirichlet.set_IWGs( { tIWGDirichlet } );

        // create a cell of set info
        moris::Cell< fem::Set_User_Info > tSetInfo( 4 );
        tSetInfo( 0 ) = tSetBulk1;
        tSetInfo( 1 ) = tSetBulk2;
        tSetInfo( 2 ) = tSetNeumann;
        tSetInfo( 3 ) = tSetDirichlet;

        // create model
        mdl::Model * tModel = new mdl::Model( &tMeshManager,
                                               tBSplineMeshIndex,
                                               tSetInfo );
        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------

        // clean up
        delete tModel;
    }

}/* END_TEST_CASE */

}/* END_MORIS_NAMESPACE */




