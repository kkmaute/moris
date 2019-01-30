 /*
 * cl_GE_test.cpp
 *
 *  Created on: Dec 28, 2018
 *      Author: sonne
 */

#include "catch.hpp"

// GE includes
//------------------------------------------------------------------------------
#include "cl_GE_Factory.hpp"
#include "cl_GE_Main.hpp"
#include "cl_GE_Element.hpp"
#include "cl_GE_Node.hpp"
#include "fn_equal_to.hpp"

#include "cl_SDF_Generator.hpp"

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
#include "op_plus.hpp"
#include "cl_Matrix.hpp"
#include "fn_print.hpp"
#include "fn_all_true.hpp"
#include "op_equal_equal.hpp"

//------------------------------------------------------------------------------
// other includes


//------------------------------------------------------------------------------

namespace moris
{
real
circle_function( const Matrix< DDRMat > & aPoint, Cell< real > inputs )
{
	// inputs(0) = x location of center
	// inputs(1) = y location of center
	// inputs(2) = radius
	Matrix< DDRMat > tCenterVec(1,2);
	tCenterVec(0,0) = inputs(0);
	tCenterVec(0,1) = inputs(1);
	return norm( aPoint - tCenterVec ) - inputs(2);
}

real
linear_function( const Matrix< DDRMat > & aPoint, moris::Cell< real> inputs )
{
	// inputs(0) = slope
	// inputs(1) = y_intercept
	return inputs(0)*aPoint(0,0) + inputs(1) - aPoint(0,1);
}

real
hyperboloid_function( const Matrix< DDRMat > & aPoint, Cell< real > inputs )
{
	// inputs(0) = offset in x
	// inputs(1) = curvature in x
	// inputs(2) = offset in y
	// inputs(3) = curvature in y
	// inputs(4) = offset in z
	// inputs(5) = curvature in z
	return std::pow(( aPoint(0,0) - inputs(0) )/inputs(1),2.0) + std::pow(( aPoint(0,1) - inputs(2))/inputs(3),2.0) - std::pow(( aPoint(0,2) - inputs(4))/inputs(5),2.0);
}

real
sphere_function( const Matrix< DDRMat > & aPoint, Cell< real > inputs )
{
	// inputs(0) = x location of center
	// inputs(1) = y location of center
	// inputs(2) = z location of center
	// inputs(3) = radius
	Matrix< DDRMat > tCenterVec(1,3);
	tCenterVec(0,0) = inputs(0);
	tCenterVec(0,1) = inputs(1);
	tCenterVec(0,2) = inputs(2);
	return norm( aPoint - tCenterVec ) - inputs(3);
}

real
plane_function( const Matrix< DDRMat > & aPoint, Cell< real > inputs )
{
	return (aPoint(0,0) - inputs(0)) + (aPoint(0,1) - inputs(1)) + (aPoint(0,2) - inputs(2));
}

	namespace ge
	{
//------------------------------------------------------------------------------

		TEST_CASE("GE test1","[GE],[GE_test1]")
				{
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

		TEST_CASE("GE test2","[GE],[GE_test2]")
				{
					mtk::Vertex* tVertex1 = new Node(0.0, 0.0); //create node objects
					mtk::Vertex* tVertex2 = new Node(1.0, 0.0);
					mtk::Vertex* tVertex3 = new Node(1.0, 1.0);
					mtk::Vertex* tVertex4 = new Node(0.0, 1.0);
					mtk::Vertex* tVertex5 = new Node(0.5, 0.0);
					mtk::Vertex* tVertex6 = new Node(1.0, 0.5);
					mtk::Vertex* tVertex7 = new Node(0.5, 1.0);
					mtk::Vertex* tVertex8 = new Node(0.0, 0.5);
					//------------------------------------------------------------------------------

					moris::Cell< mtk::Vertex* > Name(8);
					Name(0) = tVertex1; //collect node objects into cell
					Name(1) = tVertex2;
					Name(2) = tVertex3;
					Name(3) = tVertex4;
					Name(4) = tVertex5;
					Name(5) = tVertex6;
					Name(6) = tVertex7;
					Name(7) = tVertex8;
					//------------------------------------------------------------------------------

					mtk::Cell* tElement = new Element(Name); //create element object

					moris::Cell< mtk::Cell* > Elems(1);
					Elems(0) = tElement;
					//------------------------------------------------------------------------------

					moris::Cell< real > tInputs(3);
					tInputs(0) = 0.0;   // x location of center
					tInputs(1) = 0.0;   // y location of center
					tInputs(2) = 1.2;   // radius

					Ge_Factory tFactory;
					Geometry* type0 = tFactory.pick_flag(flagType::Analytical);
					type0->set_analytical_function( circle_function );

					GE geometryEngine;
					geometryEngine.set_geometry( type0 );

					moris::Cell< double > tThreshVals(1); // cell with list of threshold values
					tThreshVals(0) = 0.0;

					geometryEngine.set_threshold( tThreshVals );	// set threshold value for LS

					moris::Cell< uint > tRefFlag;

					tRefFlag = geometryEngine.flag_element_list_for_refinement( Elems, tInputs, 0 , 0);
					CHECK( equal_to( tRefFlag( 0 ), 1 ) );
					//------------------------------------------------------------------------------
					delete tVertex1; delete tVertex2;
					delete tVertex3; delete tVertex4;
					delete tVertex5; delete tVertex6;
					delete tVertex7; delete tVertex8;
					delete tElement;
				}
//------------------------------------------------------------------------------

		TEST_CASE("GE test3","[GE],[GE_test3]")
				{
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

				Ge_Factory tFactory;
				Geometry* type0 = tFactory.pick_flag(flagType::Analytical);
				type0->set_analytical_function( circle_function );

				GE geometryEngine;
				geometryEngine.set_geometry( type0 );

				moris::Cell< double > tThreshVals(2); // cell with list of threshold values
				tThreshVals(0) = 0.0;

				geometryEngine.set_threshold( tThreshVals );

				moris::Cell< uint > tRefFlag;

				tRefFlag = geometryEngine.flag_element_list_for_refinement( Elems, tInputs, 0 , 0);

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

		TEST_CASE("GE test4","[GE],[GE_test4]")
				{
					moris::Cell< real > tInput1(3);		// define inputs for circle function
					tInput1(0) = 0.0;   // x location of center
					tInput1(1) = 0.0;   // y location of center
					tInput1(2) = 1.5; //radius of circle 1
					//------------------------------------------------------------------------------

					moris::Cell< real > linInputs(2); //define parameters for the linear function
					linInputs(0) = 1.0;
					linInputs(1) = 0.5;
					//------------------------------------------------------------------------------

					mtk::Vertex* tVertex1_1 = new Node(0.0, -1.0);
					mtk::Vertex* tVertex1_2 = new Node(1.0, -1.0);
					mtk::Vertex* tVertex1_3 = new Node(1.0,  0.0);
					mtk::Vertex* tVertex1_4 = new Node(0.0,  0.0);

					mtk::Vertex* tVertex2_2 = new Node(2.0, -1.0);
					mtk::Vertex* tVertex2_3 = new Node(2.0,  0.0);

					mtk::Vertex* tVertex3_3 = new Node(2.0,  1.0);
					mtk::Vertex* tVertex3_4 = new Node(1.0,  1.0);

					mtk::Vertex* tVertex4_4 = new Node(0.0,  1.0);
					//------------------------------------------------------------------------------

					moris::Cell< mtk::Vertex* > nodes1(4);		// define the nodes in each element
					nodes1(0) = tVertex1_1; nodes1(1) = tVertex1_2;
					nodes1(2) = tVertex1_3; nodes1(3) = tVertex1_4;

					moris::Cell< mtk::Vertex* > nodes2(4);
					nodes2(0) = tVertex1_2; nodes2(1) = tVertex2_2;
					nodes2(2) = tVertex2_3; nodes2(3) = tVertex1_3;

					moris::Cell< mtk::Vertex* > nodes3(4);
					nodes3(0) = tVertex1_3; nodes3(1) = tVertex2_3;
					nodes3(2) = tVertex3_3; nodes3(3) = tVertex3_4;

					moris::Cell< mtk::Vertex* > nodes4(4);
					nodes4(0) = tVertex1_4; nodes4(1) = tVertex1_3;
					nodes4(2) = tVertex3_4; nodes4(3) = tVertex4_4;
					//------------------------------------------------------------------------------

					mtk::Cell* tElement1 = new Element(nodes1);
					mtk::Cell* tElement2 = new Element(nodes2);
					mtk::Cell* tElement3 = new Element(nodes3);
					mtk::Cell* tElement4 = new Element(nodes4);

					moris::Cell< mtk::Cell* > Elems(4);
					Elems(0) = tElement1; Elems(1) = tElement2;
					Elems(2) = tElement3; Elems(3) = tElement4;
					//------------------------------------------------------------------------------

					Ge_Factory tFactory;
					Geometry* type0 = tFactory.pick_flag(flagType::Analytical);
					type0->set_analytical_function(circle_function);

					Geometry* type1 = tFactory.pick_flag(flagType::Analytical);
					type1->set_analytical_function(linear_function);

					GE geometryEngine;

					geometryEngine.set_geometry( type0 );
					geometryEngine.set_geometry( type1 );

					moris::Cell< double > tThreshVals(1); // cell with list of threshold values
					tThreshVals(0) = 0.0;

					geometryEngine.set_threshold( tThreshVals );

					moris::Cell< uint > tRefFlag0;
					moris::Cell< uint > tRefFlag1;

					tRefFlag0 = geometryEngine.flag_element_list_for_refinement( Elems, tInput1, 0, 0 );
					tRefFlag1 = geometryEngine.flag_element_list_for_refinement( Elems, linInputs, 1, 0 );
					//------------------------------------------------------------------------------
					CHECK( equal_to( tRefFlag0( 0 ), 0 ) );			CHECK( equal_to( tRefFlag1( 0 ), 0 ) );
					CHECK( equal_to( tRefFlag0( 1 ), 1 ) );			CHECK( equal_to( tRefFlag1( 1 ), 0 ) );
					CHECK( equal_to( tRefFlag0( 2 ), 1 ) );			CHECK( equal_to( tRefFlag1( 2 ), 0 ) );
					CHECK( equal_to( tRefFlag0( 3 ), 0 ) );			CHECK( equal_to( tRefFlag1( 3 ), 1 ) );
					//------------------------------------------------------------------------------
					delete tVertex1_1; delete tVertex1_2;
					delete tVertex1_3; delete tVertex1_4;
					delete tVertex2_2; delete tVertex2_3;
					delete tVertex3_3; delete tVertex3_4;
					delete tVertex4_4;
					delete tElement1;  delete tElement2;
					delete tElement3;  delete tElement4;
				}
//------------------------------------------------------------------------------

		TEST_CASE("GE test5","[GE],[GE_test5]")
				{
				/* create 3D mtk mesh, add two sphere LS functions, loop through all elements and flag for intersection */
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

				// Compute nodal sphere values for mesh
				uint tNumNodes = tMesh3DHexs->get_num_entities(EntityRank::NODE);
				Matrix< DDRMat > tNodeSphere1Vals(1,tNumNodes);
				Matrix< DDRMat > tNodeSphere2Vals(1,tNumNodes);

				// Collect nodal sphere values
				for(uint i=0; i<tNumNodes; i++)
				{
					Matrix< DDRMat > tNodeCoord = tMesh3DHexs->get_node_coordinate(i);

					tNodeSphere1Vals(i) = sphere_function(tNodeCoord,tInput1);

					tNodeSphere2Vals(i) = sphere_function(tNodeCoord,tInput2);
				}
				// add nodal sphere values to mesh
				tMesh3DHexs->add_mesh_field_real_scalar_data_loc_inds(tSphere1FieldName, EntityRank::NODE, tNodeSphere1Vals);
				tMesh3DHexs->add_mesh_field_real_scalar_data_loc_inds(tSphere2FieldName, EntityRank::NODE, tNodeSphere2Vals);

				//------------------------------------------------------------------------------

				Ge_Factory tFactory;
				Geometry* type0 = tFactory.pick_flag(flagType::Analytical);
				type0->set_analytical_function(sphere_function);

				Geometry* type1 = tFactory.pick_flag(flagType::Analytical);
				type1->set_analytical_function(sphere_function);

				GE geometryEngine;

				geometryEngine.set_geometry( type0 );
				geometryEngine.set_geometry( type1 );

				geometryEngine.set_mesh( tMesh3DHexs );

				moris::Cell< uint > tFlags0; // create cell of flags for the elements: 1=flagged, 0=unflagged
				moris::Cell< uint > tFlags1;

				moris::Cell< double > tThreshVals(1); // cell with list of threshold values
				tThreshVals(0) = 0.0;

				geometryEngine.set_threshold( tThreshVals );

				tFlags0 = geometryEngine.check_for_intersection( tInput1, 0 , 0, 0 );
				tFlags1 = geometryEngine.check_for_intersection( tInput2, 1 , 0, 0 );

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

				std::string tOutputFile = "./ge_test5.exo";
				tMesh3DHexs->create_output_mesh(tOutputFile);
				//------------------------------------------------------------------------------

				/* add checks for test */

				delete tMesh3DHexs;
				}
//------------------------------------------------------------------------------
		TEST_CASE("GE test6","[GE],[GE_test6]")
				{
				/*	create a 2D MORIS mesh of quad4's using MTK database and determine the edge normals */
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

				GE geometryEngine;

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
				TEST_CASE("GE test7","[GE],[GE_test7]")
						{
					if(par_size()<=1)
					{
						/*	create a 3D Hex MORIS mesh using MTK database and determine the edge normals */
						//------------------------------------------------------------------------------
						const std::string tFileName2 = "generated:1x1x1";
					    // Initialize field information container
					    moris::mtk::MtkFieldsInfo tFieldsInfo;
					    // Declare some supplementary fields
					    mtk::MtkMeshData tMeshData;
					    tMeshData.FieldsInfo = &tFieldsInfo;
						// Create MORIS mesh using MTK database
						mtk::Mesh* tMesh3DHexs_norms = mtk::create_mesh( MeshType::STK, tFileName2, &tMeshData );
						//------------------------------------------------------------------------------
						uint tInd = 1;
						Matrix< IdMat > tTemp= tMesh3DHexs_norms->get_nodes_connected_to_element_glob_ids( tInd );
						print( tTemp, "tempMatrix: ");

						for( uint w=0; w<6; w++ )
						{
						Matrix< DDRMat > tNodeCoords = tMesh3DHexs_norms->get_node_coordinate( tTemp( 0, w ) );
						print( tNodeCoords, "node coords: ");
						}
						//------------------------------------------------------------------------------
//						std::string tOutputFile = "./ge_test7_edgeNorms.exo";
//						tMesh3DHexs_norms->create_output_mesh(tOutputFile);
						delete tMesh3DHexs_norms;
					}
						}
//------------------------------------------------------------------------------
				TEST_CASE("GE test8","[GE],[GE_test8]")
						{
					if(par_size()<=1)
					{
						// Define background mesh size
						const std::string tFileName = "generated:6x6x6";

						// create field for SDF
					    moris::mtk::Scalar_Field_Info<DDRMat> tElementFlagField;
					    std::string tRefineFieldName = "SDF";
					    tElementFlagField.set_field_name(tRefineFieldName);
					    tElementFlagField.set_field_entity_rank(EntityRank::ELEMENT);

					    // Initialize field information container
					    moris::mtk::MtkFieldsInfo tFieldsInfo;

					    // Place the node field into the field info container
					    add_field_for_mesh_input(&tElementFlagField,tFieldsInfo);

					    // Declare some supplementary fields
					    mtk::MtkMeshData tMeshData;
					    tMeshData.FieldsInfo = &tFieldsInfo;

						// Create MORIS mesh using MTK database
						mtk::Mesh* tMesh = mtk::create_mesh( MeshType::STK, tFileName, &tMeshData );
						//------------------------------------------------------------------------------

					    std::string tObjectPath = "/projects/GEN/test/hemisphere.obj";
					    // get path for STL file to load
					    tObjectPath = std::getenv("MORISROOT") + tObjectPath;

					    // create SDF generator
					    sdf::SDF_Generator tSdfGen( tObjectPath );

					    for( uint k=0; k<3; k++ )
					    {
					    //matrices with surface element IDs
					    Matrix< IndexMat > tSurfaceElements;
					    tSdfGen.raycast( tMesh, tSurfaceElements );		// perform a raycast on the mesh and determine number of surface elements
					    uint tNumberOfSurfaceElements = tSurfaceElements.length();
					    }
					    //------------------------------------------------------------------------------

					    // calculate SDF
					    uint tNumNodes = tMesh->get_num_entities(EntityRank::NODE);
					    for(uint i=0; i<tNumNodes; i++)
						{
							Matrix< DDRMat > tNodeCoord = tMesh->get_node_coordinate(i);
							tSdfGen.calculate_sdf( tMesh, tNodeCoord );
						}


					    /* This test is not completed yet! */

//						std::string tOutputFile = "./ge_test7.exo";
//						tMesh->create_output_mesh(tOutputFile);

						delete tMesh;
					}
						}
		//------------------------------------------------------------------------------

	} /* namespace ge */
} /* namespace moris */

