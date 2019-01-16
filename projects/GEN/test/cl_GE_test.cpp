 /*
 * cl_GE_test.cpp
 *
 *  Created on: Dec 28, 2018
 *      Author: sonne
 */

#include "catch.hpp"

#include "cl_GE_Factory.hpp"
#include "cl_GE_Main.hpp"

#include "cl_GE_Element.hpp"
#include "cl_GE_Node.hpp"
#include "fn_equal_to.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Vertex.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Tools.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Scalar_Field_Info.hpp"

#include "op_plus.hpp"


// linalg includes
//#include "cl_Matrix.hpp"
//#include "linalg_typedefs.hpp"

//#include "fn_print.hpp"
//#include "op_equal_equal.hpp"
//#include "fn_all_true.hpp"


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
//	return pow(aPoint(0,0),2.0) + pow(aPoint(0,0),2.0) - pow(aPoint(0,0),2.0);
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
				//    if(par_size()<=1)
				//    {
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
				//    }
				}
//------------------------------------------------------------------------------

		TEST_CASE("GE test2","[GE],[GE_test2]")
				{
				//    if(par_size()<=1)
				//    {
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
				tInputs(2) = 1.2; // radius

				Ge_Factory tFactory;
				Geometry* type0 = tFactory.pick_flag(flagType::Analytical);
				type0->set_analytical_function( circle_function );

				GE geometryEngine;
				geometryEngine.set_geometry( type0 );

				moris::Cell< double > tThreshVals(1); // cell with list of threshold values
				tThreshVals(0) = 0.0;

				geometryEngine.set_threshold( tThreshVals );

				moris::Cell< uint > tRefFlag;

				tRefFlag = geometryEngine.flag_element_list_for_refinement( Elems, tInputs, 0 , 0);

				//------------------------------------------------------------------------------
				delete tVertex1; delete tVertex2;
				delete tVertex3; delete tVertex4;
				delete tVertex5; delete tVertex6;
				delete tVertex7; delete tVertex8;
				delete tElement;

				//    }
				}
//------------------------------------------------------------------------------

		TEST_CASE("GE test3","[GE],[GE_test3]")
				{
				//    if(par_size()<=1)
				//    {
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
				tInputs(0) = 0.0; // x location of center
				tInputs(1) = 0.0;  // y location of center
				tInputs(2) = 1.2; //define the radius of the circle

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

				/*  elements 1, 2, 3 should be flagged
				 * 	and element 4 should not */
				//------------------------------------------------------------------------------
				delete tVertex1_1; delete tVertex1_2;
				delete tVertex1_3; delete tVertex1_4;
				delete tVertex2_3; delete tVertex2_4;
				delete tVertex3_3; delete tVertex3_4;
				delete tVertex4_3; delete tVertex4_4;
				delete tElement1;  delete tElement2;
				delete tElement3;  delete tElement4;
				//    }
				}
//------------------------------------------------------------------------------

		TEST_CASE("GE test4","[GE],[GE_test4]")
				{
				//    if(par_size()<=1)
				//    {
				moris::Cell< real > tInput1(3);
				tInput1(0) = 0.0;   // x location of center
				tInput1(1) = 0.0;   // y location of center
				tInput1(2) = 1.5; //radius of circle 1
	            //------------------------------------------------------------------------------

				moris::Cell< real > tInput2(3);
				tInput2(0) = 0.0;
				tInput2(1) = 0.0;
				tInput2(2) = 0.5;
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

				moris::Cell< mtk::Vertex* > nodes1(4);
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
				delete tVertex1_1; delete tVertex1_2;
				delete tVertex1_3; delete tVertex1_4;
				delete tVertex2_2; delete tVertex2_3;
				delete tVertex3_3; delete tVertex3_4;
				delete tVertex4_4;
				delete tElement1;  delete tElement2;
				delete tElement3;  delete tElement4;

				//    }
				}
//------------------------------------------------------------------------------

		TEST_CASE("GE test5","[GE],[GE_test5]")
				{
				//    if(par_size()<=1)
				//    {

				// Define background mesh size
				const std::string tFileName2 = "generated:6x6x6";

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
				tInput2(0) = 2.00; // x location of center
				tInput2(1) = 2.00; // y location of center
				tInput2(2) = 2.00; // z location of center
				tInput2(3) = 1.50; // radius of sphere
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

				delete tMesh3DHexs;
				//    }
				}
//------------------------------------------------------------------------------
		TEST_CASE("GE test6","[GE],[GE_test6]")
				{
				//    if(par_size()<=1)
				//	  {

				/* create two separate background meshes and multiple geometries
				 *
				 * place one geometry in each mesh
				 *
				 * Geometries:
				 * 1st - sphere
				 * 2nd - hyperboloid
				 *
				 * Meshs:
				 * 1st - 10x10x10
				 * 2nd - 5x5x5
				 *
				 */

				// Define background mesh size
				const std::string tFileName0 = "generated:10x10x10";
				const std::string tFileName1 = "generated:5x5x5";

				//------------------------------------------------------------------------------
				// Declare scalar node fieldS
				moris::mtk::Scalar_Field_Info< DDRMat > tNodeSphereField;
				std::string tSphereFieldName = "sphere";
				tNodeSphereField.set_field_name(tSphereFieldName);
				tNodeSphereField.set_field_entity_rank(EntityRank::NODE);

				moris::mtk::Scalar_Field_Info< DDRMat > tNodeParaboloidField;
				std::string tHyperboloidFieldName = "hyperboloid";
				tNodeParaboloidField.set_field_name(tHyperboloidFieldName);
				tNodeParaboloidField.set_field_entity_rank(EntityRank::NODE);

			    moris::mtk::Scalar_Field_Info< DDRMat > tElementFlagField;
			    std::string tRefineFieldName = "refinement_flags";
			    tElementFlagField.set_field_name(tRefineFieldName);
			    tElementFlagField.set_field_entity_rank(EntityRank::ELEMENT);

			    //------------------------------------------------------------------------------
			    // Initialize field information container
			    moris::mtk::MtkFieldsInfo tFieldsInfo0; // container for first mesh

			    moris::mtk::MtkFieldsInfo tFieldsInfo1; // container for second mesh

			    // Add fields to mesh
			    add_field_for_mesh_input( & tNodeSphereField, tFieldsInfo0 );
			    add_field_for_mesh_input( & tNodeParaboloidField, tFieldsInfo0 );
			    add_field_for_mesh_input( & tElementFlagField, tFieldsInfo0 );

			    // Declare some supplementary fields
			    mtk::MtkMeshData tMeshData0;
				tMeshData0.FieldsInfo = & tFieldsInfo0;

				mtk::MtkMeshData tMeshData1;
				tMeshData1.FieldsInfo = & tFieldsInfo1;
				// Create MORIS mesh using MTK database
				mtk::Mesh* tMesh3DHexs0 = mtk::create_mesh( MeshType::STK, tFileName0, & tMeshData0 );

				mtk::Mesh* tMesh3DHexs1 = mtk::create_mesh( MeshType::STK, tFileName1, & tMeshData1 );
				//------------------------------------------------------------------------------
				moris::Cell< real > tInput1(4); //sphere
				tInput1(0) = 1.50;
				tInput1(1) = 1.50;
				tInput1(2) = 1.50;
				tInput1(3) = 1.50;

				moris::Cell< real > tInput2(6); // hyperboloid
				tInput2(0) = 5.00;
				tInput2(1) = 1.00;
				tInput2(2) = 5.00;
				tInput2(3) = 1.00;
				tInput2(4) = 5.00;
				tInput2(5) = 1.00;

				//------------------------------------------------------------------------------
				// Compute nodal sphere values for mesh
				uint tNumNodes = tMesh3DHexs0->get_num_entities(EntityRank::NODE);
				Matrix< DDRMat > tNodeSphereVals(1,tNumNodes);
				Matrix< DDRMat > tNodeParaboloidVals(1,tNumNodes);

				// Collect nodal sphere values
				for(uint i=0; i<tNumNodes; i++)
				{
					Matrix< DDRMat > tNodeCoord = tMesh3DHexs0->get_node_coordinate(i);

					tNodeSphereVals(i) = sphere_function(tNodeCoord,tInput1);

					tNodeParaboloidVals(i) = sphere_function(tNodeCoord,tInput2);

				}
				// add nodal values to mesh
				tMesh3DHexs0->add_mesh_field_real_scalar_data_loc_inds(tSphereFieldName, EntityRank::NODE, tNodeSphereVals);

				tMesh3DHexs0->add_mesh_field_real_scalar_data_loc_inds(tHyperboloidFieldName, EntityRank::NODE, tNodeParaboloidVals);
				//------------------------------------------------------------------------------

				Ge_Factory tFactory;

				Geometry* type0 = tFactory.pick_flag(flagType::Analytical);
				type0->set_analytical_function(sphere_function);

				Geometry* type1 = tFactory.pick_flag(flagType::Analytical);
				type1->set_analytical_function(hyperboloid_function);

				GE geometryEngine;

				geometryEngine.set_geometry(type0);
				geometryEngine.set_geometry(type1);

				geometryEngine.set_mesh(tMesh3DHexs0);

				moris::Cell< double > tThreshVals(2); // cell with list of threshold values
				tThreshVals(0) = 0.0;
				tThreshVals(1) = 1.0;

				geometryEngine.set_threshold( tThreshVals );
				//------------------------------------------------------------------------------

				moris::Cell< uint > tFlags0; // create cell of flags for the elements: 1=flagged, 0=unflagged
				moris::Cell< uint > tFlags1;

				tFlags0 = geometryEngine.check_for_intersection( tInput1, 0 , 0, 0); //sphere
				tFlags1 = geometryEngine.check_for_intersection( tInput2, 1 , 0, 1); //hyperboloid
				//------------------------------------------------------------------------------

				uint tNumElems = tMesh3DHexs0->get_num_elems();

				Matrix< DDRMat > tElemFlag0(1,tNumElems);
				Matrix< DDRMat > tElemFlag1(1,tNumElems);

				for (uint i=0; i<tNumElems; i++)
				{
					tElemFlag0(0,i) = (real)tFlags0(i);
					tElemFlag1(0,i) = (real)tFlags1(i);
				}

				Matrix< DDRMat > tElemFlag_total(1,tNumElems);
				tElemFlag_total = tElemFlag0 + tElemFlag1;

				tMesh3DHexs0->add_mesh_field_real_scalar_data_loc_inds(tRefineFieldName, EntityRank::ELEMENT, tElemFlag_total);

				std::string tOutputFile = "./ge_test6.exo";
				tMesh3DHexs0->create_output_mesh(tOutputFile);
				//------------------------------------------------------------------------------

				delete tMesh3DHexs0;
				//    }
				}
//------------------------------------------------------------------------------


	} /* namespace ge */
} /* namespace moris */


