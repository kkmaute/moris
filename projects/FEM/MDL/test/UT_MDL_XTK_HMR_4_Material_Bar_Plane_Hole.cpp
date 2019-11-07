/*
 * UT_MDL_XTK_HMR_Multi_Material_Bar_Plane.cpp
 *
 *  Created on: Oct 4, 2019
 *      Author: doble
 */

#include "../../../GEN/src/new/geometry/cl_GEN_Geom_Field.hpp"
#include "../../../GEN/src/new/geometry/cl_GEN_Geometry.hpp"
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
#include "cl_FEM_Property_User_Defined_Info.hpp"              //FEM/INT/src
#include "cl_FEM_IWG_User_Defined_Info.hpp"              //FEM/INT/src
#include "cl_FEM_Constitutive_User_Defined_Info.hpp"      //FEM/INT/src

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

Cell< fem::IWG_User_Defined_Info >
create_iso_diff_bulk_iwg()
{
    Cell< fem::IWG_User_Defined_Info > tBulkIWG(1);

    tBulkIWG(0) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_BULK,
                                       { MSI::Dof_Type::TEMP },
                                       {{ MSI::Dof_Type::TEMP }},
                                       {fem::Property_Type::TEMP_LOAD },
                                       { fem::Constitutive_Type::DIFF_LIN_ISO } );

    return tBulkIWG;
}

Cell< fem::IWG_User_Defined_Info >
create_iso_diff_dirichlet_iwg()
{

    Cell< fem::IWG_User_Defined_Info > tIWG(1);
    tIWG(0) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_DIRICHLET,
                                           { MSI::Dof_Type::TEMP },
                                           {{ MSI::Dof_Type::TEMP }},
                                           { fem::Property_Type::TEMP_DIRICHLET },
                                           { fem::Constitutive_Type::DIFF_LIN_ISO } );
    return tIWG;
}

Cell< fem::IWG_User_Defined_Info >
create_iso_diff_neumann_iwg()
{

    Cell< fem::IWG_User_Defined_Info > tIWG(1);
    tIWG(0) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_NEUMANN,
                                          { MSI::Dof_Type::TEMP },
                                          {{ MSI::Dof_Type::TEMP }},
                                          { fem::Property_Type::TEMP_NEUMANN },
                                          moris::Cell< fem::Constitutive_Type >( 0 ) );
    return tIWG;
}

Cell< fem::IWG_User_Defined_Info >
create_iso_diff_interface_iwg()
{

    Cell< fem::IWG_User_Defined_Info > tIWG(1);
    tIWG(0) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_INTERFACE,
                                          { MSI::Dof_Type::TEMP },
                                          {{ MSI::Dof_Type::TEMP }},
                                          Cell< fem::Property_Type >( 0 ),
                                          {fem::Constitutive_Type::DIFF_LIN_ISO },
                                          {{ MSI::Dof_Type::TEMP }},
                                          Cell< fem::Property_Type >( 0 ),
                                          {fem::Constitutive_Type::DIFF_LIN_ISO } );
    return tIWG;
}

Cell< Cell< fem::Property_User_Defined_Info > >
create_bulk_properties(fem::Property_User_Defined_Info & aConductivityProp,
                       fem::Property_User_Defined_Info & aHeatLoadProp)
{
    Cell< Cell< fem::Property_User_Defined_Info > > tProps(1);
    tProps.resize( 1 );
    tProps( 0 ).resize( 2 );
    tProps( 0 )( 0 ) = aConductivityProp;
    tProps( 0 )( 1 ) = aHeatLoadProp;

    return tProps;
}

Cell< Cell< fem::Property_User_Defined_Info > >
create_dirichlet_properties(fem::Property_User_Defined_Info & aConductivity,
                            fem::Property_User_Defined_Info & aDirchletTemp)
{
    Cell< Cell< fem::Property_User_Defined_Info > > tProps(1);
    tProps.resize( 1 );
    tProps( 0 ).resize( 2 );
    tProps( 0 )( 0 ) = aConductivity;
    tProps( 0 )( 1 ) = aDirchletTemp;

    return tProps;

}

Cell< Cell< fem::Property_User_Defined_Info > >
create_neumann_properties(fem::Property_User_Defined_Info & aNeumannFlux)
{
    Cell< Cell< fem::Property_User_Defined_Info > > tProps(1);
    tProps.resize( 1 );
    tProps( 0 ).resize( 1 );
    tProps( 0 )( 0 ) = aNeumannFlux;

    return tProps;
}

Cell< Cell< fem::Property_User_Defined_Info > >
create_interface_properties(fem::Property_User_Defined_Info & aMasterCond,
                            fem::Property_User_Defined_Info & aSlaveCond)
{
    Cell< Cell< fem::Property_User_Defined_Info > > tProps(2);
    tProps( 0 ).resize( 1 );
    tProps( 0 )( 0 ) = aMasterCond;
    tProps( 1 ).resize( 1 );
    tProps( 1 )( 0 ) = aSlaveCond;

    return tProps;
}

fem::Constitutive_User_Defined_Info
create_diff_lin_constitutive_info()
{
    return fem::Constitutive_User_Defined_Info( fem::Constitutive_Type::DIFF_LIN_ISO, {{ MSI::Dof_Type::TEMP }}, { fem::Property_Type::CONDUCTIVITY } );
}


Cell< Cell< fem::Constitutive_User_Defined_Info > >
create_bulk_diff_lin_constitutive( fem::Constitutive_User_Defined_Info & aDiffLinConst )
{
    Cell< Cell< fem::Constitutive_User_Defined_Info > > tConstitutiveUserDefInfo(1);
    tConstitutiveUserDefInfo(0).resize(1);
    tConstitutiveUserDefInfo(0)(0) = aDiffLinConst;

    return tConstitutiveUserDefInfo;
}

Cell< Cell< fem::Constitutive_User_Defined_Info > >
create_dbc_diff_lin_constitutive( fem::Constitutive_User_Defined_Info & aDiffLinConst )
{
    Cell< Cell< fem::Constitutive_User_Defined_Info > > tConstitutiveUserDefInfo(1);
    tConstitutiveUserDefInfo(0).resize(1);
    tConstitutiveUserDefInfo(0)(0) = aDiffLinConst;

    return tConstitutiveUserDefInfo;
}

Cell< Cell< fem::Constitutive_User_Defined_Info > >
create_interface_diff_lin_constitutive( fem::Constitutive_User_Defined_Info & aMasterDiffLinConst,
                                        fem::Constitutive_User_Defined_Info & aSlaveDiffLinConst)
{
    Cell< Cell< fem::Constitutive_User_Defined_Info > > tConstitutiveUserDefInfo(2);
    tConstitutiveUserDefInfo(0).resize(1);
    tConstitutiveUserDefInfo(0)(0) = aMasterDiffLinConst;

    tConstitutiveUserDefInfo(1).resize(1);
    tConstitutiveUserDefInfo(1)(0) = aSlaveDiffLinConst;

    return tConstitutiveUserDefInfo;
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

        // create IWG user defined info
        Cell< fem::IWG_User_Defined_Info > tBulkIWG = create_iso_diff_bulk_iwg();
        Cell< fem::IWG_User_Defined_Info > tDBCIWG  = create_iso_diff_dirichlet_iwg();
        Cell< fem::IWG_User_Defined_Info > tNBCIWG  = create_iso_diff_neumann_iwg();
        Cell< fem::IWG_User_Defined_Info > tIntIWG  = create_iso_diff_interface_iwg();

        Cell< Cell< fem::IWG_User_Defined_Info > > tIWGUserDefinedInfo( 14 );

        tIWGUserDefinedInfo( 0 )  = tBulkIWG;
        tIWGUserDefinedInfo( 1 )  = tBulkIWG;
        tIWGUserDefinedInfo( 2 )  = tBulkIWG;
        tIWGUserDefinedInfo( 3 )  = tBulkIWG;
        tIWGUserDefinedInfo( 4 )  = tBulkIWG;
        tIWGUserDefinedInfo( 5 )  = tBulkIWG;
        tIWGUserDefinedInfo( 6 )  = tBulkIWG;
        tIWGUserDefinedInfo( 7 )  = tBulkIWG;
        tIWGUserDefinedInfo( 8 )  = tDBCIWG;
        tIWGUserDefinedInfo( 9 )  = tNBCIWG;
        tIWGUserDefinedInfo( 10 ) = tIntIWG;
        tIWGUserDefinedInfo( 11 ) = tIntIWG;
        tIWGUserDefinedInfo( 12 ) = tIntIWG;
        tIWGUserDefinedInfo( 13 ) = tIntIWG;

        // create the property user defined infos
        fem::Property_User_Defined_Info tConductivity( fem::Property_Type::CONDUCTIVITY,
                                                       Cell< Cell< MSI::Dof_Type > >( 0 ),
                                                       {{{ 1.0 }}},
                                                       tConstValFunction2MatMDL,
                                                       Cell< fem::PropertyFunc >( 0 ) );

        fem::Property_User_Defined_Info tConductivity2( fem::Property_Type::CONDUCTIVITY,
                                                       Cell< Cell< MSI::Dof_Type > >( 0 ),
                                                       {{{ 1.0 }}},
                                                       tConstValFunction2MatMDL,
                                                       Cell< fem::PropertyFunc >( 0 ) );
        fem::Property_User_Defined_Info tTempDirichlet( fem::Property_Type::TEMP_DIRICHLET,
                                                        Cell< Cell< MSI::Dof_Type > >( 0 ),
                                                        {{{ 5.0 }}},
                                                        tConstValFunction2MatMDL,
                                                        Cell< fem::PropertyFunc >( 0 ) );
        fem::Property_User_Defined_Info tNeumannFlux( fem::Property_Type::TEMP_NEUMANN,
                                                      Cell< Cell< MSI::Dof_Type > >( 0 ),
                                                      {{{ 20.0 }}},
                                                      tConstValFunction2MatMDL,
                                                      Cell< fem::PropertyFunc >( 0 ) );

        fem::Property_User_Defined_Info tTempLoad1( fem::Property_Type::TEMP_LOAD,
                                                              Cell< Cell< MSI::Dof_Type > >( 0 ),
                                                              {{{ 0.0 }}},
                                                              tConstValFunction2MatMDL,
                                                              Cell< fem::PropertyFunc >( 0 ) );

        fem::Property_User_Defined_Info tTempLoad2( fem::Property_Type::TEMP_LOAD,
                                                              Cell< Cell< MSI::Dof_Type > >( 0 ),
                                                              {{{ 50.0 }}},
                                                              tConstValFunction2MatMDL,
                                                              Cell< fem::PropertyFunc >( 0 ) );

        // create property user defined info
        Cell< Cell< Cell< fem::Property_User_Defined_Info > > > tPropertyUserDefinedInfo( 14 );
        tPropertyUserDefinedInfo(0)  = create_bulk_properties(tConductivity2,tTempLoad2);
        tPropertyUserDefinedInfo(1)  = create_bulk_properties(tConductivity2,tTempLoad2);
        tPropertyUserDefinedInfo(2)  = create_bulk_properties(tConductivity2,tTempLoad2);
        tPropertyUserDefinedInfo(3)  = create_bulk_properties(tConductivity2,tTempLoad2);
        tPropertyUserDefinedInfo(4)  = create_bulk_properties(tConductivity,tTempLoad1);
        tPropertyUserDefinedInfo(5)  = create_bulk_properties(tConductivity,tTempLoad1);
        tPropertyUserDefinedInfo(6)  = create_bulk_properties(tConductivity,tTempLoad1);
        tPropertyUserDefinedInfo(7)  = create_bulk_properties(tConductivity,tTempLoad1);
        tPropertyUserDefinedInfo(8)  = create_dirichlet_properties(tConductivity,tTempDirichlet);
        tPropertyUserDefinedInfo(9)  = create_neumann_properties(tNeumannFlux);
        tPropertyUserDefinedInfo(10) = create_interface_properties(tConductivity2,tConductivity2);
        tPropertyUserDefinedInfo(11) = create_interface_properties(tConductivity2,tConductivity);
        tPropertyUserDefinedInfo(12) = create_interface_properties(tConductivity2,tConductivity);
        tPropertyUserDefinedInfo(13) = create_interface_properties(tConductivity,tConductivity);

        // create constitutive user defined info
        fem::Constitutive_User_Defined_Info tDiffLinIso = create_diff_lin_constitutive_info();
        // create constitutive user defined info
        Cell< Cell< Cell< fem::Constitutive_User_Defined_Info > > > tConstitutiveUserDefinedInfo( 14 );
        tConstitutiveUserDefinedInfo(0) = create_bulk_diff_lin_constitutive(tDiffLinIso);
        tConstitutiveUserDefinedInfo(1) = create_bulk_diff_lin_constitutive(tDiffLinIso);
        tConstitutiveUserDefinedInfo(2) = create_bulk_diff_lin_constitutive(tDiffLinIso);
        tConstitutiveUserDefinedInfo(3) = create_bulk_diff_lin_constitutive(tDiffLinIso);
        tConstitutiveUserDefinedInfo(4) = create_bulk_diff_lin_constitutive(tDiffLinIso);
        tConstitutiveUserDefinedInfo(5) = create_bulk_diff_lin_constitutive(tDiffLinIso);
        tConstitutiveUserDefinedInfo(6) = create_bulk_diff_lin_constitutive(tDiffLinIso);
        tConstitutiveUserDefinedInfo(7) = create_bulk_diff_lin_constitutive(tDiffLinIso);
        tConstitutiveUserDefinedInfo(8) = create_dbc_diff_lin_constitutive(tDiffLinIso);
        tConstitutiveUserDefinedInfo( 9 ).resize( 1 ); // neumann
        tConstitutiveUserDefinedInfo( 10 ) = create_interface_diff_lin_constitutive(tDiffLinIso,tDiffLinIso);
        tConstitutiveUserDefinedInfo( 11 ) = create_interface_diff_lin_constitutive(tDiffLinIso,tDiffLinIso);
        tConstitutiveUserDefinedInfo( 12 ) = create_interface_diff_lin_constitutive(tDiffLinIso,tDiffLinIso);
        tConstitutiveUserDefinedInfo( 13 ) = create_interface_diff_lin_constitutive(tDiffLinIso,tDiffLinIso);

        // create a list of active block-sets
        std::string tDblInterfaceSideSetName01 = tEnrIntegMesh.get_dbl_interface_side_set_name(0,1);
        std::string tDblInterfaceSideSetName02 = tEnrIntegMesh.get_dbl_interface_side_set_name(0,2);
        std::string tDblInterfaceSideSetName13 = tEnrIntegMesh.get_dbl_interface_side_set_name(1,3);
        std::string tDblInterfaceSideSetName23 = tEnrIntegMesh.get_dbl_interface_side_set_name(2,3);

        std::cout<<"tDblInterfaceSideSetName01 = "<<tDblInterfaceSideSetName01<<" | Index = "<<tEnrIntegMesh.get_double_sided_set_index(tDblInterfaceSideSetName01)<<std::endl;
        std::cout<<"tDblInterfaceSideSetName02 = "<<tDblInterfaceSideSetName02<<" | Index = "<<tEnrIntegMesh.get_double_sided_set_index(tDblInterfaceSideSetName02)<<std::endl;


        moris::Cell< moris_index >  tSetList = {  tEnrIntegMesh.get_block_set_index("HMR_dummy_c_p0"),
                                                  tEnrIntegMesh.get_block_set_index("HMR_dummy_n_p0"),
                                                  tEnrIntegMesh.get_block_set_index("HMR_dummy_c_p1"),
                                                  tEnrIntegMesh.get_block_set_index("HMR_dummy_n_p1"),
                                                  tEnrIntegMesh.get_block_set_index("HMR_dummy_c_p2"),
                                                  tEnrIntegMesh.get_block_set_index("HMR_dummy_n_p2"),
                                                  tEnrIntegMesh.get_block_set_index("HMR_dummy_c_p3"),
                                                  tEnrIntegMesh.get_block_set_index("HMR_dummy_n_p3"),
                                                  tEnrIntegMesh.get_side_set_index("SideSet_4_n_p2"),
                                                  tEnrIntegMesh.get_side_set_index("SideSet_2_n_p3"),
                                                  tEnrIntegMesh.get_double_sided_set_index(tDblInterfaceSideSetName01),
                                                  tEnrIntegMesh.get_double_sided_set_index(tDblInterfaceSideSetName02),
                                                  tEnrIntegMesh.get_double_sided_set_index(tDblInterfaceSideSetName13),
                                                  tEnrIntegMesh.get_double_sided_set_index(tDblInterfaceSideSetName23)};

        moris::Cell< fem::Element_Type > tSetTypeList = { fem::Element_Type::BULK,
                                                          fem::Element_Type::BULK,
                                                          fem::Element_Type::BULK,
                                                          fem::Element_Type::BULK,
                                                          fem::Element_Type::BULK,
                                                          fem::Element_Type::BULK,
                                                          fem::Element_Type::BULK,
                                                          fem::Element_Type::BULK,
                                                          fem::Element_Type::SIDESET,
                                                          fem::Element_Type::SIDESET,
                                                          fem::Element_Type::DOUBLE_SIDESET,
                                                          fem::Element_Type::DOUBLE_SIDESET,
                                                          fem::Element_Type::DOUBLE_SIDESET,
                                                          fem::Element_Type::DOUBLE_SIDESET,
                                                          };


        // create model
        mdl::Model * tModel = new mdl::Model( &tMeshManager, tBSplineMeshIndex, tSetList, tSetTypeList, tIWGUserDefinedInfo, tPropertyUserDefinedInfo, tConstitutiveUserDefinedInfo, 0, false);

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
        tLinearSolverAlgorithm->set_param("AZ_max_iter") = 10;
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
        tNonlinearSolverAlgorithm->set_param("NLA_max_iter")   = 3;
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

        // create IWG user defined info
        Cell< fem::IWG_User_Defined_Info > tBulkIWG = create_iso_diff_bulk_iwg();
        Cell< fem::IWG_User_Defined_Info > tDBCIWG  = create_iso_diff_dirichlet_iwg();
        Cell< fem::IWG_User_Defined_Info > tNBCIWG  = create_iso_diff_neumann_iwg();
        Cell< fem::IWG_User_Defined_Info > tIntIWG  = create_iso_diff_interface_iwg();

        Cell< Cell< fem::IWG_User_Defined_Info > > tIWGUserDefinedInfo( 14 );

        tIWGUserDefinedInfo( 0 )  = tBulkIWG;
        tIWGUserDefinedInfo( 1 )  = tBulkIWG;
        tIWGUserDefinedInfo( 2 )  = tBulkIWG;
        tIWGUserDefinedInfo( 3 )  = tBulkIWG;
        tIWGUserDefinedInfo( 4 )  = tBulkIWG;
        tIWGUserDefinedInfo( 5 )  = tBulkIWG;
        tIWGUserDefinedInfo( 6 )  = tBulkIWG;
        tIWGUserDefinedInfo( 7 )  = tBulkIWG;
        tIWGUserDefinedInfo( 8 )  = tDBCIWG;
        tIWGUserDefinedInfo( 9 )  = tNBCIWG;
        tIWGUserDefinedInfo( 10 ) = tIntIWG;
        tIWGUserDefinedInfo( 11 ) = tIntIWG;
        tIWGUserDefinedInfo( 12 ) = tIntIWG;
        tIWGUserDefinedInfo( 13 ) = tIntIWG;

        // create the property user defined infos
        fem::Property_User_Defined_Info tConductivity( fem::Property_Type::CONDUCTIVITY,
                                                       Cell< Cell< MSI::Dof_Type > >( 0 ),
                                                       {{{ 1.0 }}},
                                                       tConstValFunction2MatMDL,
                                                       Cell< fem::PropertyFunc >( 0 ) );

        fem::Property_User_Defined_Info tConductivity2( fem::Property_Type::CONDUCTIVITY,
                                                       Cell< Cell< MSI::Dof_Type > >( 0 ),
                                                       {{{ 0.1 }}},
                                                       tConstValFunction2MatMDL,
                                                       Cell< fem::PropertyFunc >( 0 ) );
        fem::Property_User_Defined_Info tTempDirichlet( fem::Property_Type::TEMP_DIRICHLET,
                                                        Cell< Cell< MSI::Dof_Type > >( 0 ),
                                                        {{{ 5.0 }}},
                                                        tConstValFunction2MatMDL,
                                                        Cell< fem::PropertyFunc >( 0 ) );
        fem::Property_User_Defined_Info tNeumannFlux( fem::Property_Type::TEMP_NEUMANN,
                                                      Cell< Cell< MSI::Dof_Type > >( 0 ),
                                                      {{{ 20.0 }}},
                                                      tConstValFunction2MatMDL,
                                                      Cell< fem::PropertyFunc >( 0 ) );

        fem::Property_User_Defined_Info tTempLoad1( fem::Property_Type::TEMP_LOAD,
                                                              Cell< Cell< MSI::Dof_Type > >( 0 ),
                                                              {{{ 0.0 }}},
                                                              tConstValFunction2MatMDL,
                                                              Cell< fem::PropertyFunc >( 0 ) );

        fem::Property_User_Defined_Info tTempLoad2( fem::Property_Type::TEMP_LOAD,
                                                              Cell< Cell< MSI::Dof_Type > >( 0 ),
                                                              {{{ 50.0 }}},
                                                              tConstValFunction2MatMDL,
                                                              Cell< fem::PropertyFunc >( 0 ) );

        // create property user defined info
        Cell< Cell< Cell< fem::Property_User_Defined_Info > > > tPropertyUserDefinedInfo( 14 );
        tPropertyUserDefinedInfo(0)  = create_bulk_properties(tConductivity2,tTempLoad2);
        tPropertyUserDefinedInfo(1)  = create_bulk_properties(tConductivity2,tTempLoad2);
        tPropertyUserDefinedInfo(2)  = create_bulk_properties(tConductivity2,tTempLoad2);
        tPropertyUserDefinedInfo(3)  = create_bulk_properties(tConductivity2,tTempLoad2);
        tPropertyUserDefinedInfo(4)  = create_bulk_properties(tConductivity,tTempLoad1);
        tPropertyUserDefinedInfo(5)  = create_bulk_properties(tConductivity,tTempLoad1);
        tPropertyUserDefinedInfo(6)  = create_bulk_properties(tConductivity,tTempLoad1);
        tPropertyUserDefinedInfo(7)  = create_bulk_properties(tConductivity,tTempLoad1);
        tPropertyUserDefinedInfo(8)  = create_dirichlet_properties(tConductivity,tTempDirichlet);
        tPropertyUserDefinedInfo(9)  = create_neumann_properties(tNeumannFlux);
        tPropertyUserDefinedInfo(10) = create_interface_properties(tConductivity2,tConductivity2);
        tPropertyUserDefinedInfo(11) = create_interface_properties(tConductivity2,tConductivity);
        tPropertyUserDefinedInfo(12) = create_interface_properties(tConductivity2,tConductivity);
        tPropertyUserDefinedInfo(13) = create_interface_properties(tConductivity,tConductivity);

        // create constitutive user defined info
        fem::Constitutive_User_Defined_Info tDiffLinIso = create_diff_lin_constitutive_info();
        // create constitutive user defined info
        Cell< Cell< Cell< fem::Constitutive_User_Defined_Info > > > tConstitutiveUserDefinedInfo( 14 );
        tConstitutiveUserDefinedInfo(0) = create_bulk_diff_lin_constitutive(tDiffLinIso);
        tConstitutiveUserDefinedInfo(1) = create_bulk_diff_lin_constitutive(tDiffLinIso);
        tConstitutiveUserDefinedInfo(2) = create_bulk_diff_lin_constitutive(tDiffLinIso);
        tConstitutiveUserDefinedInfo(3) = create_bulk_diff_lin_constitutive(tDiffLinIso);
        tConstitutiveUserDefinedInfo(4) = create_bulk_diff_lin_constitutive(tDiffLinIso);
        tConstitutiveUserDefinedInfo(5) = create_bulk_diff_lin_constitutive(tDiffLinIso);
        tConstitutiveUserDefinedInfo(6) = create_bulk_diff_lin_constitutive(tDiffLinIso);
        tConstitutiveUserDefinedInfo(7) = create_bulk_diff_lin_constitutive(tDiffLinIso);
        tConstitutiveUserDefinedInfo(8) = create_dbc_diff_lin_constitutive(tDiffLinIso);
        tConstitutiveUserDefinedInfo( 9 ).resize( 1 ); // neumann
        tConstitutiveUserDefinedInfo( 10 ) = create_interface_diff_lin_constitutive(tDiffLinIso,tDiffLinIso);
        tConstitutiveUserDefinedInfo( 11 ) = create_interface_diff_lin_constitutive(tDiffLinIso,tDiffLinIso);
        tConstitutiveUserDefinedInfo( 12 ) = create_interface_diff_lin_constitutive(tDiffLinIso,tDiffLinIso);
        tConstitutiveUserDefinedInfo( 13 ) = create_interface_diff_lin_constitutive(tDiffLinIso,tDiffLinIso);

        // create a list of active block-sets
        std::string tDblInterfaceSideSetName01 = tEnrIntegMesh.get_dbl_interface_side_set_name(0,1);
        std::string tDblInterfaceSideSetName02 = tEnrIntegMesh.get_dbl_interface_side_set_name(0,2);
        std::string tDblInterfaceSideSetName13 = tEnrIntegMesh.get_dbl_interface_side_set_name(1,3);
        std::string tDblInterfaceSideSetName23 = tEnrIntegMesh.get_dbl_interface_side_set_name(2,3);

        std::cout<<"tDblInterfaceSideSetName01 = "<<tDblInterfaceSideSetName01<<" | Index = "<<tEnrIntegMesh.get_double_sided_set_index(tDblInterfaceSideSetName01)<<std::endl;
        std::cout<<"tDblInterfaceSideSetName02 = "<<tDblInterfaceSideSetName02<<" | Index = "<<tEnrIntegMesh.get_double_sided_set_index(tDblInterfaceSideSetName02)<<std::endl;


        moris::Cell< moris_index >  tSetList = {  tEnrIntegMesh.get_block_set_index("HMR_dummy_c_p0"),
                                                  tEnrIntegMesh.get_block_set_index("HMR_dummy_n_p0"),
                                                  tEnrIntegMesh.get_block_set_index("HMR_dummy_c_p1"),
                                                  tEnrIntegMesh.get_block_set_index("HMR_dummy_n_p1"),
                                                  tEnrIntegMesh.get_block_set_index("HMR_dummy_c_p2"),
                                                  tEnrIntegMesh.get_block_set_index("HMR_dummy_n_p2"),
                                                  tEnrIntegMesh.get_block_set_index("HMR_dummy_c_p3"),
                                                  tEnrIntegMesh.get_block_set_index("HMR_dummy_n_p3"),
                                                  tEnrIntegMesh.get_side_set_index("SideSet_4_n_p2"),
                                                  tEnrIntegMesh.get_side_set_index("SideSet_2_n_p3"),
                                                  tEnrIntegMesh.get_double_sided_set_index(tDblInterfaceSideSetName01),
                                                  tEnrIntegMesh.get_double_sided_set_index(tDblInterfaceSideSetName02),
                                                  tEnrIntegMesh.get_double_sided_set_index(tDblInterfaceSideSetName13),
                                                  tEnrIntegMesh.get_double_sided_set_index(tDblInterfaceSideSetName23)};

        moris::Cell< fem::Element_Type > tSetTypeList = { fem::Element_Type::BULK,
                                                          fem::Element_Type::BULK,
                                                          fem::Element_Type::BULK,
                                                          fem::Element_Type::BULK,
                                                          fem::Element_Type::BULK,
                                                          fem::Element_Type::BULK,
                                                          fem::Element_Type::BULK,
                                                          fem::Element_Type::BULK,
                                                          fem::Element_Type::SIDESET,
                                                          fem::Element_Type::SIDESET,
                                                          fem::Element_Type::DOUBLE_SIDESET,
                                                          fem::Element_Type::DOUBLE_SIDESET,
                                                          fem::Element_Type::DOUBLE_SIDESET,
                                                          fem::Element_Type::DOUBLE_SIDESET,
                                                          };


        // create model
        mdl::Model * tModel = new mdl::Model( &tMeshManager, tBSplineMeshIndex, tSetList, tSetTypeList, tIWGUserDefinedInfo, tPropertyUserDefinedInfo, tConstitutiveUserDefinedInfo, 0, false);

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

