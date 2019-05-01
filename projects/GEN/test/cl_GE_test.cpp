 /*
 * cl_GE_test.cpp
 *
 *  Created on: Dec 28, 2018
 *      Author: sonne
 */
#include "../src/cl_GE_Core.hpp"
#include "catch.hpp"

// GE includes
//------------------------------------------------------------------------------
#include "cl_GE_Element.hpp"
#include "cl_GE_Factory.hpp"
#include "cl_GE_Node.hpp"

#include "cl_SDF_Generator.hpp"

//------------------------------------------------------------------------------
// HMR includes
//#include "HMR_Globals.hpp"
//#include "cl_HMR_Parameters.cpp"
//#include "cl_HMR.hpp"
//#include "cl_HMR_Field.hpp"
//#include "cl_HMR_Mesh.hpp"

//------------------------------------------------------------------------------
// MTK includes

#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Vertex.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Tools.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Scalar_Field_Info.hpp"

//------------------------------------------------------------------------------
// linalg includes

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_all_true.hpp"
#include "fn_equal_to.hpp"
#include "op_equal_equal.hpp"

//------------------------------------------------------------------------------
// FEM includes
#include "cl_FEM_Integrator.hpp" //FEM/INT/src/
#include "cl_FEM_Field_Interpolator.hpp"

//------------------------------------------------------------------------------
// other includes
#include "cl_Stopwatch.hpp"
#include "cl_Profiler.hpp"

//------------------------------------------------------------------------------

namespace moris
{
	namespace ge
	{
//------------------------------------------------------------------------------

		TEST_CASE("nodal_coordinate_test","[GE],[nodal_coordinate_test]")
				{
		        /* create GE node objects, create GE element object, check nodal coordinates */
				mtk::Vertex* tVertex1 = new Node(0.0, 0.0);
				mtk::Vertex* tVertex2 = new Node(2.0, 0.0);
				mtk::Vertex* tVertex3 = new Node(2.0, 1.0);
				mtk::Vertex* tVertex4 = new Node(0.0, 1.0);
				//------------------------------------------------------------------------------

				moris::Cell< mtk::Vertex* > Name(4);
				Name(0) = tVertex1;
				Name(1) = tVertex2;
				Name(2) = tVertex3;
				Name(3) = tVertex4;

				mtk::Cell* tElement = new Element(Name);
				//------------------------------------------------------------------------------

				Matrix< DDRMat > Name2 = tElement->get_vertex_pointers()(0)->get_coords();
				CHECK( equal_to( Name2( 0,0 ), 0.0 ) );
				CHECK( equal_to( Name2( 0,1 ), 0.0 ) );

				Matrix< DDRMat > Name3 = tElement->get_vertex_pointers()(1)->get_coords();
				CHECK( equal_to( Name3( 0,0 ), 2.0 ) );
				CHECK( equal_to( Name3( 0,1 ), 0.0 ) );

				Matrix< DDRMat > Name4 = tElement->get_vertex_pointers()(2)->get_coords();
				CHECK( equal_to( Name4( 0,0 ), 2.0 ) );
				CHECK( equal_to( Name4( 0,1 ), 1.0 ) );

				Matrix< DDRMat > Name5 = tElement->get_vertex_pointers()(3)->get_coords();
				CHECK( equal_to( Name5( 0,0 ), 0.0 ) );
				CHECK( equal_to( Name5( 0,1 ), 1.0 ) );

				//------------------------------------------------------------------------------
				delete tVertex1; delete tVertex2;
				delete tVertex3; delete tVertex4;
				delete tElement;
				}

//------------------------------------------------------------------------------

		TEST_CASE("2D_circle_LS_intersection_test","[GE],[intersection_test_2D]")
				{
		        /* create 2D 4-element mesh, add circle function, check for intersection */
				mtk::Vertex* tVertex1_1 = new Node(0.0, -1.0);
				mtk::Vertex* tVertex1_2 = new Node(1.0, -1.0);
				mtk::Vertex* tVertex1_3 = new Node(1.0,  0.0);
				mtk::Vertex* tVertex1_4 = new Node(0.0,  0.0);

				mtk::Vertex* tVertex2_3 = new Node(1.0,  1.0);
				mtk::Vertex* tVertex2_4 = new Node(0.0,  1.0);

				mtk::Vertex* tVertex3_3 = new Node(1.0,  2.0);
				mtk::Vertex* tVertex3_4 = new Node(0.0,  2.0);

				mtk::Vertex* tVertex4_3 = new Node(1.0,  3.0);
				mtk::Vertex* tVertex4_4 = new Node(0.0,  3.0);
				//------------------------------------------------------------------------------

				moris::Cell< mtk::Vertex* > Name1(4); //element 1 nodes
				Name1(0) = tVertex1_1;
				Name1(1) = tVertex1_2;
				Name1(2) = tVertex1_3;
				Name1(3) = tVertex1_4;

				moris::Cell< mtk::Vertex* > Name2(4);
				Name2(0) = tVertex1_4;
				Name2(1) = tVertex1_3;
				Name2(2) = tVertex2_3;
				Name2(3) = tVertex2_4;

				moris::Cell< mtk::Vertex* > Name3(4);
				Name3(0) = tVertex2_4;
				Name3(1) = tVertex2_3;
				Name3(2) = tVertex3_3;
				Name3(3) = tVertex3_4;

				moris::Cell< mtk::Vertex* > Name4(4);
				Name4(0) = tVertex3_4;
				Name4(1) = tVertex3_4;
				Name4(2) = tVertex4_3;
				Name4(3) = tVertex4_4;
				//------------------------------------------------------------------------------

				mtk::Cell* tElement1 = new Element(Name1);
				mtk::Cell* tElement2 = new Element(Name2);
				mtk::Cell* tElement3 = new Element(Name3);
				mtk::Cell* tElement4 = new Element(Name4);

				moris::Cell< mtk::Cell* > Elems(4);
				Elems(0) = tElement1;
				Elems(1) = tElement2;
				Elems(2) = tElement3;
				Elems(3) = tElement4;
				//------------------------------------------------------------------------------
				moris::Cell< real > tInputs(3);
				tInputs(0) = 0.0;   // x location of center
				tInputs(1) = 0.0;   // y location of center
				tInputs(2) = 1.2;   //define the radius of the circle

				//------------------------------------------------------------------------------
				Ge_Factory tFactory;

				std::shared_ptr< Geometry > type0 = tFactory.set_geometry_type(type::ANALYTIC);
				type0->set_analytical_function( type::CIRCLE );

				GE_Core geometryEngine;
				geometryEngine.set_geometry( type0 );

				moris::Cell< uint > tRefFlag;

				tRefFlag = geometryEngine.flag_element_list_for_refinement( Elems, tInputs, 0 );

				CHECK( equal_to( tRefFlag( 0 ), 1 ) );
				CHECK( equal_to( tRefFlag( 1 ), 1 ) );
				CHECK( equal_to( tRefFlag( 2 ), 1 ) );
				CHECK( equal_to( tRefFlag( 3 ), 0 ) );
				//------------------------------------------------------------------------------
				delete tVertex1_1; delete tVertex1_2;
				delete tVertex1_3; delete tVertex1_4;
				delete tVertex2_3; delete tVertex2_4;
				delete tVertex3_3; delete tVertex3_4;
				delete tVertex4_3; delete tVertex4_4;
				delete tElement1;  delete tElement2;
				delete tElement3;  delete tElement4;
				}

//------------------------------------------------------------------------------

		TEST_CASE("GE_multiple_geometry_intersection_test_3D_with_mtk_mesh","[GE],[intersection_test_3D]")
				{
				/* create 3D mtk mesh, add two sphere LS functions, loop through all elements
				 * and flag for intersection
				 */
				//------------------------------------------------------------------------------
				const std::string tFileName2 = "generated:6x6x6"; // Define background mesh size

			    // Declare scalar node field 1
			    moris::mtk::Scalar_Field_Info<DDRMat> tNodeSphereField1;
			    std::string tSphere1FieldName = "sphere1";
			    tNodeSphereField1.set_field_name(tSphere1FieldName);
			    tNodeSphereField1.set_field_entity_rank(EntityRank::NODE);
			    
			    // Declare scalar node field 2
			    moris::mtk::Scalar_Field_Info<DDRMat> tNodeSphereField2;
			    std::string tSphere2FieldName = "sphere2";
			    tNodeSphereField2.set_field_name(tSphere2FieldName);
			    tNodeSphereField2.set_field_entity_rank(EntityRank::NODE);

			    // Boolean field for flagged elements
			    moris::mtk::Scalar_Field_Info<DDRMat> tElementFlagField;
			    std::string tRefineFieldName = "refinement_flags";
			    tElementFlagField.set_field_name(tRefineFieldName);
			    tElementFlagField.set_field_entity_rank(EntityRank::ELEMENT);

			    // Initialize field information container
			    moris::mtk::MtkFieldsInfo tFieldsInfo;

			    // Place the node field into the field info container
			    add_field_for_mesh_input(&tNodeSphereField1,tFieldsInfo);
			    add_field_for_mesh_input(&tNodeSphereField2,tFieldsInfo);

			    add_field_for_mesh_input(&tElementFlagField,tFieldsInfo);

			    // Declare some supplementary fields
			    mtk::MtkMeshData tMeshData;
			    tMeshData.FieldsInfo = &tFieldsInfo;

				// Create MORIS mesh using MTK database
				mtk::Mesh* tMesh3DHexs = mtk::create_mesh( MeshType::STK, tFileName2, &tMeshData );
				//------------------------------------------------------------------------------

				// Define parameters for sphere1
				moris::Cell< real > tInput1(4);
				tInput1(0) = 4.00; // x location of center
				tInput1(1) = 4.00; // y location of center
				tInput1(2) = 4.00; // z location of center
				tInput1(3) = 2.00; // radius of sphere

				// Define parameters for sphere2
				moris::Cell< real > tInput2(4);
				tInput2(0) = 2.00;
				tInput2(1) = 2.00;
				tInput2(2) = 2.00;
				tInput2(3) = 1.50;
				//------------------------------------------------------------------------------
                Ge_Factory tFactory;
                std::shared_ptr< Geometry > type0 = tFactory.set_geometry_type(type::ANALYTIC);
                type0->set_analytical_function( type::SPHERE );

                std::shared_ptr< Geometry > type1 = tFactory.set_geometry_type(type::ANALYTIC);
                type1->set_analytical_function( type::SPHERE );

				// Compute nodal sphere values for mesh
				uint tNumNodes = tMesh3DHexs->get_num_entities(EntityRank::NODE);
				Matrix< DDRMat > tNodeSphere1Vals(1,tNumNodes);
				Matrix< DDRMat > tNodeSphere2Vals(1,tNumNodes);

				// Collect nodal sphere values
				for(uint i=0; i<tNumNodes; i++)
				{
					Matrix< DDRMat > tNodeCoord = tMesh3DHexs->get_node_coordinate(i);

					tNodeSphere1Vals(i) = type0->get_field_val_at_coordinate( tNodeCoord,tInput1 );
					tNodeSphere2Vals(i) = type1->get_field_val_at_coordinate( tNodeCoord,tInput2 );
				}
				// add nodal sphere values to mesh
				tMesh3DHexs->add_mesh_field_real_scalar_data_loc_inds(tSphere1FieldName, EntityRank::NODE, tNodeSphere1Vals);
				tMesh3DHexs->add_mesh_field_real_scalar_data_loc_inds(tSphere2FieldName, EntityRank::NODE, tNodeSphere2Vals);

				//------------------------------------------------------------------------------

				std::shared_ptr<Geometry_Engine_Interface> tInterface = std::make_shared< GE_Core >();

				tInterface->set_geometry( type0 );
				tInterface->set_geometry( type1 );

				moris::Cell< uint > tFlags0; // create cell of flags for the elements: 1=flagged, 0=unflagged
				moris::Cell< uint > tFlags1;

				tFlags0 = tInterface->check_for_intersection( tInput1, tMesh3DHexs , 0 );
				tFlags1 = tInterface->check_for_intersection( tInput2, tMesh3DHexs , 1 );

				//------------------------------------------------------------------------------

				uint tNumElems = tMesh3DHexs->get_num_elems();

				Matrix< DDRMat > tElementalFlags1(1,tNumElems); // need to cast the cell flag values into real values for STK output
				Matrix< DDRMat > tElementalFlags2(1,tNumElems);

				for(uint t=0; t<tNumElems; t++)
				{
					tElementalFlags1(0,t) = (real)tFlags0(t);
					tElementalFlags2(0,t) = (real)tFlags1(t);
				}

				Matrix< DDRMat > tElementalFlags_total(1,tNumElems); // create total flag field due to both geometries
				tElementalFlags_total = tElementalFlags1 + tElementalFlags2;

				// add total flag field to mesh
				tMesh3DHexs->add_mesh_field_real_scalar_data_loc_inds(tRefineFieldName, EntityRank::ELEMENT, tElementalFlags_total);

//				std::string tOutputFile = "./ge_test5.exo";
//				tMesh3DHexs->create_output_mesh(tOutputFile);
				//------------------------------------------------------------------------------

				/* fixme need to add checks for test */

				delete tMesh3DHexs;
//				delete type0;       delete type1;
				}
//------------------------------------------------------------------------------
		TEST_CASE("2D_quad4_edge_normal_test_with_specific_mtk_mesh","[GE],[normal_test_2D_quad4]")
				{
				/*	create a 2D MORIS mesh of quad4's using mtk database and determine the edge normals */
				//------------------------------------------------------------------------------
				uint aNumElemTypes = 1;		// quad
				uint aNumDim = 2;		// specify number of spatial dimensions

				Matrix< IdMat > aElementConnQuad = {{ 1, 2, 5, 8 },
													{ 2, 3, 4, 5 },
													{ 8, 5, 6, 7 },
													{ 5, 4, 9, 6 }};		// specify element connectivity of quad for mesh

				Matrix< IdMat > aElemLocalToGlobalQuad = {{ 1, 2, 3, 4 }};		// specify the local to global element map for quads

				Matrix< DDRMat > aCoords = {{ 0.0, 0.0 },
											{ 1.0, 0.0 },
											{ 2.0, 0.0 },
											{ 2.0, 1.0 },
											{ 1.0, 1.0 },
											{ 1.0, 2.0 },
											{ 0.0, 2.0 },
											{ 0.0, 1.0 },
											{ 2.0, 2.0 }};		// Node coordinate matrix

				Matrix< IdMat > aNodeLocalToGlobal = {{ 1, 2, 3, 4, 5, 6, 7, 8, 9 }};		// specify the local to global map
				//------------------------------------------------------------------------------
				// create MORIS mesh using MTK database
				mtk::MtkMeshData aMeshData( aNumElemTypes );
				aMeshData.CreateAllEdgesAndFaces = true;
				aMeshData.SpatialDim = & aNumDim;
				aMeshData.ElemConn(0) = & aElementConnQuad;
				aMeshData.NodeCoords = & aCoords;
				aMeshData.LocaltoGlobalElemMap(0) = & aElemLocalToGlobalQuad;
				aMeshData.LocaltoGlobalNodeMap = & aNodeLocalToGlobal;

				mtk::Mesh* tMesh2D_Quad4 = create_mesh( MeshType::STK, aMeshData );
				//------------------------------------------------------------------------------

				GE_Core geometryEngine;

				Matrix< DDRMat > tNormal( 16, 2 );

				uint count = 0;
				for(uint f=1; f<5; f++)
				{
					for (uint e=0; e<4; e++)
					{
						Matrix< DDRMat > tTemp = geometryEngine.get_edge_normal_for_straight_edge_quad4( f, e, tMesh2D_Quad4 );
						tNormal( count, 0 ) = tTemp( 0, 0 ); tNormal( count, 1 ) = tTemp( 0, 1 );
						++count;
					}
				}
				//------------------------------------------------------------------------------
				Matrix< DDRMat > tCheckMatrix ( 16, 2 );
				tCheckMatrix(  0, 0 ) =  0.0;			tCheckMatrix(  0, 1 ) = -1.0;
				tCheckMatrix(  1, 0 ) =  1.0;			tCheckMatrix(  1, 1 ) =  0.0;
				tCheckMatrix(  2, 0 ) =  0.0;			tCheckMatrix(  2, 1 ) =  1.0;
				tCheckMatrix(  3, 0 ) = -1.0;			tCheckMatrix(  3, 1 ) =  0.0;
				tCheckMatrix(  4, 0 ) =  0.0;			tCheckMatrix(  4, 1 ) = -1.0;
				tCheckMatrix(  5, 0 ) =  1.0;			tCheckMatrix(  5, 1 ) =  0.0;
				tCheckMatrix(  6, 0 ) =  0.0;			tCheckMatrix(  6, 1 ) =  1.0;
				tCheckMatrix(  7, 0 ) = -1.0;			tCheckMatrix(  7, 1 ) =  0.0;
				tCheckMatrix(  8, 0 ) =  0.0;			tCheckMatrix(  8, 1 ) = -1.0;
				tCheckMatrix(  9, 0 ) =  1.0;			tCheckMatrix(  9, 1 ) =  0.0;
				tCheckMatrix( 10, 0 ) =  0.0;			tCheckMatrix( 10, 1 ) =  1.0;
				tCheckMatrix( 11, 0 ) = -1.0;			tCheckMatrix( 11, 1 ) =  0.0;
				tCheckMatrix( 12, 0 ) =  0.0;			tCheckMatrix( 12, 1 ) = -1.0;
				tCheckMatrix( 13, 0 ) =  1.0;			tCheckMatrix( 13, 1 ) =  0.0;
				tCheckMatrix( 14, 0 ) =  0.0;			tCheckMatrix( 14, 1 ) =  1.0;
				tCheckMatrix( 15, 0 ) = -1.0;			tCheckMatrix( 15, 1 ) =  0.0;

				bool tNormalMatrixMatch = all_true( tNormal == tCheckMatrix );
				CHECK( tNormalMatrixMatch );
				//------------------------------------------------------------------------------
				delete tMesh2D_Quad4;
				}
//------------------------------------------------------------------------------

				TEST_CASE("GE_calculate_phi_values_at_nodes","[GE],[L2_projection_test]")
				{
                    if(par_size()<=1)
                    {
                        // create T-Matrix to be used for test (9x9)
                        //------------------------------------------------------------------------------
                        Matrix< DDRMat > tTMat( 9, 9);
                        tTMat(0,0) = 0.2500;        tTMat(0,1) = 0.0000;        tTMat(0,2) = 0.0000;
                        tTMat(0,3) = 0.0000;        tTMat(0,4) = 0.2500;        tTMat(0,5) = 0.0000;
                        tTMat(0,6) = 0.0000;        tTMat(0,7) = 0.2500;        tTMat(0,8) = 0.2500;

                        tTMat(1,0) = 0.0000;        tTMat(1,1) = 0.2500;        tTMat(1,2) = 0.0000;
                        tTMat(1,3) = 0.0000;        tTMat(1,4) = 0.2500;        tTMat(1,5) = 0.2500;
                        tTMat(1,6) = 0.0000;        tTMat(1,7) = 0.0000;        tTMat(1,8) = 0.2500;

                        tTMat(2,0) = 0.0000;        tTMat(2,1) = 0.0000;        tTMat(2,2) = 0.2500;
                        tTMat(2,3) = 0.0000;        tTMat(2,4) = 0.0000;        tTMat(2,5) = 0.2500;
                        tTMat(2,6) = 0.2500;        tTMat(2,7) = 0.0000;        tTMat(2,8) = 0.2500;

                        tTMat(3,0) = 0.0000;        tTMat(3,1) = 0.0000;        tTMat(3,2) = 0.0000;
                        tTMat(3,3) = 0.2500;        tTMat(3,4) = 0.0000;        tTMat(3,5) = 0.0000;
                        tTMat(3,6) = 0.2500;        tTMat(3,7) = 0.2500;        tTMat(3,8) = 0.2500;

                        tTMat(4,0) = 0.0625;        tTMat(4,1) = 0.0625;        tTMat(4,2) = 0.0000;
                        tTMat(4,3) = 0.0000;        tTMat(4,4) = 0.3750;        tTMat(4,5) = 0.0625;
                        tTMat(4,6) = 0.0000;        tTMat(4,7) = 0.0625;        tTMat(4,8) = 0.3750;

                        tTMat(5,0) = 0.0000;        tTMat(5,1) = 0.0625;        tTMat(5,2) = 0.0625;
                        tTMat(5,3) = 0.0000;        tTMat(5,4) = 0.0625;        tTMat(5,5) = 0.3750;
                        tTMat(5,6) = 0.0625;        tTMat(5,7) = 0.0000;        tTMat(5,8) = 0.3750;

                        tTMat(6,0) = 0.0000;        tTMat(6,1) = 0.0000;        tTMat(6,2) = 0.0625;
                        tTMat(6,3) = 0.0625;        tTMat(6,4) = 0.0000;        tTMat(6,5) = 0.0625;
                        tTMat(6,6) = 0.3750;        tTMat(6,7) = 0.0625;        tTMat(6,8) = 0.3750;

                        tTMat(7,0) = 0.0625;        tTMat(7,1) = 0.0000;        tTMat(7,2) = 0.0000;
                        tTMat(7,3) = 0.0625;        tTMat(7,4) = 0.0625;        tTMat(7,5) = 0.0000;
                        tTMat(7,6) = 0.0625;        tTMat(7,7) = 0.3750;        tTMat(7,8) = 0.3750;

                        tTMat(8,0) = 0.0156;        tTMat(8,1) = 0.0156;        tTMat(8,2) = 0.0156;
                        tTMat(8,3) = 0.0156;        tTMat(8,4) = 0.0938;        tTMat(8,5) = 0.0938;
                        tTMat(8,6) = 0.0938;        tTMat(8,7) = 0.0938;        tTMat(8,8) = 0.5625;
                        //------------------------------------------------------------------------------

                        // create mtk mesh to be used for test
                        //------------------------------------------------------------------------------
                        uint aNumElemTypes = 1;     // quad
                        uint aNumDim = 2;       // specify number of spatial dimensions

                        Matrix< IdMat > aElementConnQuad9 = {{ 1, 2, 3, 4, 5, 6, 7, 8, 9 }};        // specify element connectivity of quad for mesh

                        Matrix< IdMat > aElemLocalToGlobalQuad9 = {{ 1 }};      // specify the local to global element map for quads

                        Matrix< DDRMat > aCoords = {{ 0.0, 0.0 },
                                                    { 5.0, 0.0 },
                                                    { 5.0, 5.0 },
                                                    { 0.0, 5.0 },
                                                    { 2.5, 0.0 },
                                                    { 5.0, 2.5 },
                                                    { 2.5, 5.0 },
                                                    { 0.0, 2.5 },
                                                    { 2.5, 2.5 }};      // Node coordinate matrix

                        Matrix< IdMat > aNodeLocalToGlobal = {{ 1, 2, 3, 4, 5, 6, 7, 8, 9 }};       // specify the local to global map

                        // create MORIS mesh using MTK database
                        mtk::MtkMeshData aMeshData( aNumElemTypes );
                        aMeshData.CreateAllEdgesAndFaces = true;
                        aMeshData.SpatialDim = & aNumDim;
                        aMeshData.ElemConn(0) = & aElementConnQuad9;
                        aMeshData.NodeCoords = & aCoords;
                        aMeshData.LocaltoGlobalElemMap(0) = & aElemLocalToGlobalQuad9;
                        aMeshData.LocaltoGlobalNodeMap = & aNodeLocalToGlobal;

                        mtk::Mesh* tMesh2D_Quad9 = create_mesh( MeshType::STK, aMeshData);
                        //------------------------------------------------------------------------------

                        // setup the geometry engine and specify LS function parameters (circle)
                        //------------------------------------------------------------------------------
                        Ge_Factory tFactory;
                        std::shared_ptr< Geometry > tGeomType = tFactory.set_geometry_type(type::ANALYTIC);
                        tGeomType->set_analytical_function( type::CIRCLE );

                        tGeomType->set_mesh_and_t_matrix(tMesh2D_Quad9, tTMat);

                        GE_Core tGeometryEngine;
                        tGeometryEngine.set_geometry( tGeomType );

                        moris::Cell< real > tCircleInputs(3);
                        tCircleInputs(0) = 0.00; // x location of center
                        tCircleInputs(1) = 0.00; // y location of center
                        tCircleInputs(2) = 2.50; // radius of circle
                        //------------------------------------------------------------------------------

                        // call function to determine nodal ADVs and use it to check phi values at nodes
                        //------------------------------------------------------------------------------
                        uint tWhichGeom = 0;
                        Matrix< DDRMat > tPhiHat = tGeometryEngine.compute_nodal_advs( tWhichGeom,
                                                                                       tCircleInputs,
                                                                                       mtk::Geometry_Type::QUAD,
                                                                                       fem::Integration_Type::GAUSS,
                                                                                       fem::Integration_Order::QUAD_3x3,
                                                                                       fem::Interpolation_Type::LAGRANGE,
                                                                                       mtk::Interpolation_Order::QUADRATIC,
                                                                                       fem::Integration_Type::GAUSS,
                                                                                       fem::Integration_Order::BAR_1,
                                                                                       fem::Interpolation_Type::LAGRANGE,
                                                                                       mtk::Interpolation_Order::LINEAR,
                                                                                       fem::Interpolation_Type::CONSTANT,
                                                                                       mtk::Interpolation_Order::CONSTANT);

                        Matrix< DDRMat > tPhi = tTMat*tPhiHat;

                        //------------------------------------------------------------------------------

                        // check values
                        //------------------------------------------------------------------------------
                        real tEpsilon     = 0.000000000001;    // error tolerance 1e-12

                        Matrix< DDRMat > tCheckMatrix(9,1);
                        tCheckMatrix(0,0) = -6.249999999999949;
                        tCheckMatrix(1,0) = 18.750000000000099;
                        tCheckMatrix(2,0) = 43.750000000000121;
                        tCheckMatrix(3,0) = 18.750000000000085;
                        tCheckMatrix(4,0) = -0.000000000000035;
                        tCheckMatrix(5,0) = 24.999999999999936;
                        tCheckMatrix(6,0) = 24.999999999999943;
                        tCheckMatrix(7,0) = -0.000000000000039;
                        tCheckMatrix(8,0) =  6.250000000000020;

                        Matrix< DDRMat > tDiffMat = tPhi - tCheckMatrix;

                        for (uint i=0; i<9; i++)
                        {
                            REQUIRE( tDiffMat(i) < tEpsilon);
                        }
                        //------------------------------------------------------------------------------

//                        std::string tOutputFile = "./phiVals.exo";
//                        tMesh2D_Quad9->create_output_mesh(tOutputFile);
                        //------------------------------------------------------------------------------
                        delete tMesh2D_Quad9;
                    }
				}

//------------------------------------------------------------------------------

//				TEST_CASE("GE_SDF_generator","[GE],[SDF_generator_test]")
//						{
////						Profiler tProf("/home/sonne/Desktop/temp_profile");
//
//						hmr::ParameterList tParameters = hmr::create_hmr_parameter_list();
//					    tParameters.set( "number_of_elements_per_dimension", "10, 10, 10" );
//					    tParameters.set( "domain_dimensions",                "5.6, 2.6, 3.4" );
//					    tParameters.set( "domain_offset",                    "-4.9, 3.25, -1.7" );
//					    tParameters.set( "verbose", 1 );
//
//
//						std::string tObjectPath = "/projects/GEN/test/objfiles/designbox.obj";
//						tObjectPath = std::getenv("MORISROOT") + tObjectPath;
//						sdf::SDF_Generator tSdfGen( tObjectPath );
//
//						moris::hmr::HMR tHMR( tParameters );
//					    auto tMesh = tHMR.create_mesh();
//					    for( uint k=0; k<3; ++k )
//					    {
//					       // matrices with surface element IDs
//					       Matrix< IndexMat > tSurfaceElements;
//
//					       tSdfGen.raycast( tMesh, tSurfaceElements );
//					       // get number of surface elements
//					       uint tNumberOfSurfaceElements = tSurfaceElements.length();
//
//					       // loop over all elements
//					       for( uint e=0; e<tNumberOfSurfaceElements; ++e )
//					       {
//					           // manually flag element
//					           tHMR.flag_element( tSurfaceElements( e ) );
//					       }
//
//					       // refine
//					       tHMR.perform_refinement( moris::hmr::RefinementMode::SIMPLE  );
//					    }
//
//					    // calculate T-Matrices etc
//					    tHMR.finalize();
//
//					    // calculate SDF
//					    auto tField = tMesh->create_field( "SDF", 1);
//
//					//------------------------------------------------------------------------------
//
//					    tic tTimer;
//
//					    tSdfGen.calculate_sdf( tMesh, tField->get_node_values() );
//
//					    real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;
//
//					    std::cout<<"time for SDF generator         : "<<tElapsedTime/1000<<" [sec]"<<std::endl;
//					    tHMR.save_to_exodus( "genTestSDF.exo" );
//
////						tProf.stop();
//						}


	} /* namespace ge */
} /* namespace moris */

