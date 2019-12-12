
/*
 * UT_VIS_Visualization_Mesh.cpp
 *
 *  Created on: Dez 02, 2019
 *      Author: schmidt
 */

#include "catch.hpp"
#include "typedefs.hpp"
#include "cl_Map.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_equal_to.hpp"



#define protected public
#define private   public
#include "cl_MDL_Model.hpp"
#include "cl_VIS_Factory.hpp"
#include "cl_VIS_Visualization_Mesh.hpp"
#include "cl_VIS_Output_Manager.hpp"
#undef protected
#undef private

#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Background_Mesh.hpp" //HMR/src
#include "cl_HMR_BSpline_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Element.hpp" //HMR/src
#include "cl_HMR_Factory.hpp" //HMR/src
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Parameters.hpp" //HMR/src

#include "cl_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Tools.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Scalar_Field_Info.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
#include "../projects/GEN/src/geometry/cl_GEN_Geom_Field.hpp"
#include "../projects/GEN/src/geometry/cl_GEN_Geometry.hpp"
#include "../projects/GEN/src/geometry/cl_GEN_Geometry.hpp"

#include "cl_VIS_Factory.hpp"

#include "cl_MTK_Writer_Exodus.hpp"
#include "cl_MTK_Reader_Exodus.hpp"

#include "cl_FEM_NodeProxy.hpp"                //FEM/INT/src
#include "cl_FEM_ElementProxy.hpp"             //FEM/INT/src
#include "cl_FEM_Node_Base.hpp"                //FEM/INT/src
#include "cl_FEM_Element_Factory.hpp"          //FEM/INT/src
#include "cl_FEM_IWG_Factory.hpp"              //FEM/INT/src
#include "cl_FEM_SP_Factory.hpp"              //FEM/INT/src
#include "cl_FEM_CM_Factory.hpp"              //FEM/INT/src
#include "cl_FEM_Set_User_Info.hpp"              //FEM/INT/src
#include "cl_FEM_Set.hpp"              //FEM/INT/src

moris::real PlaneVisTest(const moris::Matrix< moris::DDRMat > & aPoint )
{
    moris::real mXC = 0.11;
    moris::real mYC = 0.11;
    moris::real mNx = 1.0;
    moris::real mNy = 0.0;
    return (mNx*(aPoint(0)-mXC) + mNy*(aPoint(1)-mYC));
}

Matrix< DDRMat > tConstValFunction_MDLDIFF( moris::Cell< Matrix< DDRMat > >         & aParameters,
                                            moris::Cell< fem::Field_Interpolator* > & aDofFI,
                                            moris::Cell< fem::Field_Interpolator* > & aDvFI,
                                            fem::Geometry_Interpolator              * aGeometryInterpolator )
{
    return aParameters( 0 );
}

namespace moris
{
    namespace vis
    {
TEST_CASE(" Output Data","[VIS],[Output_Data]")
    {
        if(par_size() == 1)
            {
                std::string tFieldName = "Geometry";

                moris::uint tLagrangeMeshIndex = 0;
                moris::uint tBSplineMeshIndex = 0;

                moris::hmr::Parameters tParameters;

                tParameters.set_number_of_elements_per_dimension( { {4}, {2} } );
                tParameters.set_domain_dimensions({ {2}, {1} });
                tParameters.set_domain_offset({ {-1.0}, {-0.0} });
                tParameters.set_bspline_truncation( true );

                tParameters.set_output_meshes( { {0} } );

                tParameters.set_lagrange_orders  ( { {1} });
                tParameters.set_lagrange_patterns({ {0} });

                tParameters.set_bspline_orders   ( { {1} } );
                tParameters.set_bspline_patterns ( { {0} } );

                tParameters.set_side_sets({{1},{2},{3},{4} });

                tParameters.set_refinement_buffer( 1 );
                tParameters.set_staircase_buffer( 1 );

                Cell< Matrix< DDUMat > > tLagrangeToBSplineMesh( 1 );
                tLagrangeToBSplineMesh( 0 ) = { {0} };

                tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

                hmr::HMR tHMR( tParameters );

                std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

                // create field
                std::shared_ptr< moris::hmr::Field > tPlaneField  = tMesh->create_field( tFieldName, tLagrangeMeshIndex );

                tPlaneField->evaluate_scalar_function( PlaneVisTest);

                for( uint k=0; k<1; ++k )
                {
                    tHMR.flag_surface_elements_on_working_pattern( tPlaneField );
                    tHMR.perform_refinement_based_on_working_pattern( 0 );
                    tPlaneField->evaluate_scalar_function( PlaneVisTest);
                }

                tHMR.finalize();

//                tHMR.save_to_exodus( 0, "./xtk_exo/xtk_hmr_2d_enr_ip2.e" );

                std::shared_ptr< hmr::Interpolation_Mesh_HMR > tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex  );

                moris::ge::GEN_Geom_Field tPlaneFieldAsGeom( tPlaneField );

                moris::Cell< moris::ge::GEN_Geometry* > tGeometryVector = {&tPlaneFieldAsGeom};

                size_t tModelDimension = 2;
                moris::ge::GEN_Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
                moris::ge::GEN_Geometry_Engine tGeometryEngine(tGeometryVector,tPhaseTable,tModelDimension);
                xtk::Model tXTKModel(tModelDimension,tInterpMesh.get(),tGeometryEngine);
                tXTKModel.mVerbose = false;

                //Specify decomposition Method and Cut Mesh ---------------------------------------
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
                std::shared_ptr< fem::Property > tPropConductivity = std::make_shared< fem::Property >();
                tPropConductivity->set_parameters( { {{ 1.0 }} } );
                tPropConductivity->set_val_function( tConstValFunction_MDLDIFF );

                std::shared_ptr< fem::Property > tPropDirichlet = std::make_shared< fem::Property >();
                tPropDirichlet->set_parameters( { {{ 5.0 }} } );
                tPropDirichlet->set_val_function( tConstValFunction_MDLDIFF );

                std::shared_ptr< fem::Property > tPropNeumann = std::make_shared< fem::Property >();
                tPropNeumann->set_parameters( { {{ 20.0 }} } );
                tPropNeumann->set_val_function( tConstValFunction_MDLDIFF );

                std::shared_ptr< fem::Property > tPropTempLoad = std::make_shared< fem::Property >();
                tPropTempLoad->set_parameters( { {{ 0.0 }} } );
                tPropTempLoad->set_val_function( tConstValFunction_MDLDIFF );

                // define constitutive models
                fem::CM_Factory tCMFactory;

                std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIso = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
                tCMDiffLinIso->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} ); // FIXME through the factory?
                tCMDiffLinIso->set_property( tPropConductivity, "Conductivity" );
                tCMDiffLinIso->set_space_dim( 3 );

                // define stabilization parameters
                fem::SP_Factory tSPFactory;

                std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche = tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
                tSPDirichletNitsche->set_parameters( { {{ 1.0 }} } );
                tSPDirichletNitsche->set_property( tPropConductivity, "Material", mtk::Master_Slave::MASTER );

                // define the IWGs
                fem::IWG_Factory tIWGFactory;

                std::shared_ptr< fem::IWG > tIWGBulk = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
                tIWGBulk->set_residual_dof_type( { MSI::Dof_Type::TEMP } );
                tIWGBulk->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
                tIWGBulk->set_constitutive_model( tCMDiffLinIso, "DiffLinIso", mtk::Master_Slave::MASTER );
                tIWGBulk->set_property( tPropTempLoad, "Load", mtk::Master_Slave::MASTER );

                std::shared_ptr< fem::IWG > tIWGDirichlet = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_DIRICHLET );
                tIWGDirichlet->set_residual_dof_type( { MSI::Dof_Type::TEMP } );
                tIWGDirichlet->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
                tIWGDirichlet->set_stabilization_parameter( tSPDirichletNitsche, "DirichletNitsche" );
                tIWGDirichlet->set_constitutive_model( tCMDiffLinIso, "DiffLinIso", mtk::Master_Slave::MASTER );
                tIWGDirichlet->set_property( tPropDirichlet, "Dirichlet", mtk::Master_Slave::MASTER );

                std::shared_ptr< fem::IWG > tIWGNeumann = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_NEUMANN );
                tIWGNeumann->set_residual_dof_type( { MSI::Dof_Type::TEMP } );
                tIWGNeumann->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
                tIWGNeumann->set_property( tPropNeumann, "Neumann", mtk::Master_Slave::MASTER );

                // define set info
                fem::Set_User_Info tSetBulk1;
                tSetBulk1.set_mesh_index( tEnrIntegMesh.get_block_set_index("HMR_dummy_c_p0") );
                tSetBulk1.set_set_type( fem::Element_Type::BULK );
                tSetBulk1.set_IWGs( { tIWGBulk } );

                fem::Set_User_Info tSetBulk2;
                tSetBulk2.set_mesh_index( tEnrIntegMesh.get_block_set_index("HMR_dummy_n_p0") );
                tSetBulk2.set_set_type( fem::Element_Type::BULK );
                tSetBulk2.set_IWGs( { tIWGBulk } );

                fem::Set_User_Info tSetBulk3;
                tSetBulk3.set_mesh_index( tEnrIntegMesh.get_block_set_index("HMR_dummy_c_p1") );
                tSetBulk3.set_set_type( fem::Element_Type::BULK );
                tSetBulk3.set_IWGs( { tIWGBulk } );

                fem::Set_User_Info tSetBulk4;
                tSetBulk4.set_mesh_index( tEnrIntegMesh.get_block_set_index("HMR_dummy_n_p1") );
                tSetBulk4.set_set_type( fem::Element_Type::BULK );
                tSetBulk4.set_IWGs( { tIWGBulk } );

                fem::Set_User_Info tSetDirichlet;
                tSetDirichlet.set_mesh_index( tEnrIntegMesh.get_side_set_index("SideSet_2_n_p1") );
                tSetDirichlet.set_set_type( fem::Element_Type::SIDESET );
                tSetDirichlet.set_IWGs( { tIWGDirichlet } );

                fem::Set_User_Info tSetNeumann;
                tSetNeumann.set_mesh_index( tEnrIntegMesh.get_side_set_index("SideSet_4_n_p0") );
                tSetNeumann.set_set_type( fem::Element_Type::SIDESET );
                tSetNeumann.set_IWGs( { tIWGNeumann } );

                // create a cell of set info
                moris::Cell< fem::Set_User_Info > tSetInfo( 6 );
                tSetInfo( 0 ) = tSetBulk1;
                tSetInfo( 1 ) = tSetBulk2;
                tSetInfo( 2 ) = tSetBulk3;
                tSetInfo( 3 ) = tSetBulk4;
                tSetInfo( 4 ) = tSetDirichlet;
                tSetInfo( 5 ) = tSetNeumann;

                // create model
                mdl::Model * tModel = new mdl::Model( &tMeshManager,
                        1,
                        tSetInfo );

                moris::Cell< MSI::Equation_Set * > tEquationSets(tModel->get_equation_sets().size());
                for(moris::uint iSet = 0; iSet < tModel->get_equation_sets().size(); iSet++ )
                {
                    tEquationSets( iSet ) = tModel->get_equation_sets()( iSet );
                }

                Output_Manager tOutputData;

                tOutputData.set_outputs( 0,
                                            VIS_Mesh_Type::OVERLAPPING_INTERFACE,
                                            "Output_Vis_Mesh.exo",
                                            { "HMR_dummy_c_p0", "HMR_dummy_n_p0", "HMR_dummy_c_p1", "HMR_dummy_n_p1"},
                                            { 0, 1, 2, 3 },
                                            { "pressure" },
                                            { Field_Type::ELEMENTAL },
                                            { Output_Type::UX } );

                tOutputData.create_visualization_mesh( 0,
                                                       &tMeshManager,
                                                        0 );

                tOutputData.set_visualization_sets( 0,
                                             tModel );

                tOutputData.write_mesh( 0,
                                        0.0);

                tOutputData.end_writing( 0 );
//
//                tOutputData.write_mesh();

//                tOutputData.write_field();



//----------------------------------------------------------------------------------------------

////                mtk::Mesh* tMesh11 = new Visualization_Mesh( &tMeshManager, 0 );
//
//                vis::Factory tVisFactory(&tMeshManager, 0);
//
//                mtk::Mesh * tVisMesh = tVisFactory.create_visualization_mesh();
//
//                Writer_Exodus writer(tVisMesh);
//                writer.write_mesh("/data/schmidt/codes/moris/build/", "Vis_Mesh_2.exo");
//
//                moris::Cell<const moris::mtk::Cell*> tElementsInBlock = tVisMesh->get_block_set_cells("HMR_dummy_c_p0");
//
//                uint tNumElements = tElementsInBlock.size();
//                moris::Matrix<moris::DDRMat> tetField(tNumElements, 1, 4);
//                moris::Cell<std::string> tElementalFieldNames(1);
//                tElementalFieldNames(0) = "pressure";
//
//                for(uint Ik = 0; Ik<tNumElements;Ik++)
//                {
//                    tetField( Ik ) = Ik;
//                }
//
//                writer.set_elemental_fields(tElementalFieldNames);
//                writer.set_time(0.0);
//                writer.write_elemental_field(0, "pressure", tetField);
//
//                writer.close_file();
            }
    }
    }
}


