//------------------------------------------------------------------------------
// moris core includes
#include "cl_GE_Core.hpp"
#include "cl_Communication_Manager.hpp"
#include "cl_Communication_Tools.hpp"
#include "typedefs.hpp"
//------------------------------------------------------------------------------
// necessary includes for tutorial
#include "catch.hpp"

#include "cl_GE_Factory.hpp"
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

//------------------------------------------------------------------------------

using namespace moris;
using namespace ge;

//------------------------------------------------------------------------------

moris::Comm_Manager gMorisComm;

//------------------------------------------------------------------------------

/*!
 * the following function is used as this level set for the tutorial
 *
 * \code{.cpp}
 * real
 * sphere_function( const Matrix< DDRMat > & aPoint, Cell< real > inputs )
 * {
 *	Matrix< DDRMat > tCenterVec(1,3);
 *	tCenterVec(0,0) = inputs(0);
 *	tCenterVec(0,1) = inputs(1);
 *	tCenterVec(0,2) = inputs(2);
 *	return norm( aPoint - tCenterVec ) - inputs(3);
 *  }
 * \endcode
 *
 */

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

//------------------------------------------------------------------------------

int
main(
		int    argc,
		char * argv[] )
{
	gMorisComm = moris::Comm_Manager( &argc, &argv );

	//------------------------------------------------------------------------------
	/*!
	 * <b> Step 1: setup the background mesh </b>
	 */

	/*!
	 * For this tutorial, a 6x6x6 STK mesh is defined.
	 *
	 * For more information on setting up a background mesh, please see
	 * the relevant documentation. Here, the steps are defined only briefly.
	 */

	/*!
	 * define the background mesh size:
	 *
	 * \code{.cpp}
	 * const std::string tFileName2 = "generated:6x6x6";
	 * \endcode
	 *
	 */

	const std::string tFileName2 = "generated:6x6x6";

	/*!
	 * declare scalar node field 1:
	 *
	 * \code{.cpp}
	 * moris::mtk::Scalar_Field_Info<DDRMat> tNodeSphereField1;
	 * std::string tSphere1FieldName = "sphere1";
	 * tNodeSphereField1.set_field_name(tSphere1FieldName);
	 * tNodeSphereField1.set_field_entity_rank(EntityRank::NODE);
	 * \endcode
	 *
	 */

	moris::mtk::Scalar_Field_Info<DDRMat> tNodeSphereField1;
	std::string tSphere1FieldName = "sphere1";
	tNodeSphereField1.set_field_name(tSphere1FieldName);
	tNodeSphereField1.set_field_entity_rank(EntityRank::NODE);

	/*!
	 * declare scalar node field 2
	 *
	 * \code{.cpp}
	 * moris::mtk::Scalar_Field_Info<DDRMat> tNodeSphereField2;
	 * std::string tSphere2FieldName = "sphere2";
	 * tNodeSphereField2.set_field_name(tSphere2FieldName);
	 * tNodeSphereField2.set_field_entity_rank(EntityRank::NODE);
	 * \endcode
	 *
	 */

    moris::mtk::Scalar_Field_Info<DDRMat> tNodeSphereField2;
    std::string tSphere2FieldName = "sphere2";
    tNodeSphereField2.set_field_name(tSphere2FieldName);
    tNodeSphereField2.set_field_entity_rank(EntityRank::NODE);

    /*!
     * boolean field for flagged elements
     *
     * \code{.cpp}
     * moris::mtk::Scalar_Field_Info<DDRMat> tElementFlagField;
     * std::string tRefineFieldName = "refinement_flags";
     * tElementFlagField.set_field_name(tRefineFieldName);
     * tElementFlagField.set_field_entity_rank(EntityRank::ELEMENT);
     * \endcode
     *
     */

    moris::mtk::Scalar_Field_Info<DDRMat> tElementFlagField;
    std::string tRefineFieldName = "refinement_flags";
    tElementFlagField.set_field_name(tRefineFieldName);
    tElementFlagField.set_field_entity_rank(EntityRank::ELEMENT);

    /*!
     * initialize field information container
     *
     * \code{.cpp}
     * moris::mtk::MtkFieldsInfo tFieldsInfo;
     * \endcode
     *
     */

    moris::mtk::MtkFieldsInfo tFieldsInfo;

    /*!
     * place the node fields and element flag field into the field info container
     *
     * \code{.cpp}
     * add_field_for_mesh_input(&tNodeSphereField1,tFieldsInfo);
     * add_field_for_mesh_input(&tNodeSphereField2,tFieldsInfo);
     * add_field_for_mesh_input(&tElementFlagField,tFieldsInfo);
     * \endcode
     *
     */

    add_field_for_mesh_input(&tNodeSphereField1,tFieldsInfo);
    add_field_for_mesh_input(&tNodeSphereField2,tFieldsInfo);

    add_field_for_mesh_input(&tElementFlagField,tFieldsInfo);

    /*!
     * declare some supplementary fields
     *
     * \code{.cpp}
     * mtk::MtkMeshData tMeshData;
     * tMeshData.FieldsInfo = &tFieldsInfo;
     * \endcode
     *
     */

    mtk::MtkMeshData tMeshData;
    tMeshData.FieldsInfo = &tFieldsInfo;

	//------------------------------------------------------------------------------
    /*!
     * <b> Step 2: create mesh and define geometry parameters </b>
     */

    /*!
     * create MORIS mesh using MTK database
     *
     * \code{.cpp}
     * mtk::Mesh* tMesh3DHexs = mtk::create_mesh( MeshType::STK, tFileName2, &tMeshData );
     * \endcode
     *
     * this results in the following background mesh:
     *
     * @image html ./figures/tutorial_1_backgroundMesh.png "6x6x6 Background Mesh"
     *
     */

	mtk::Mesh* tMesh3DHexs = mtk::create_mesh( MeshType::STK, tFileName2, &tMeshData );

	/*!
	 * define parameters for sphere1
	 *
	 * \code{.cpp}
	 * moris::Cell< real > tInput1(4);
	 * tInput1(0) = 4.00;
	 * tInput1(1) = 4.00;
	 * tInput1(2) = 4.00;
	 * tInput1(3) = 2.00;
	 * \endcode
	 *
	 */

	moris::Cell< real > tInput1(4);
	tInput1(0) = 4.00;
	tInput1(1) = 4.00;
	tInput1(2) = 4.00;
	tInput1(3) = 2.00;

	/*!
	 * define parameters for sphere 2
	 *
	 * \code{.cpp}
	 * moris::Cell< real > tInput2(4);
	 * tInput2(0) = 2.00;
	 * tInput2(1) = 2.00;
	 * tInput2(2) = 2.00;
	 * tInput2(3) = 1.50;
	 * \endcode
	 *
	 */

	moris::Cell< real > tInput2(4);
	tInput2(0) = 2.00;
	tInput2(1) = 2.00;
	tInput2(2) = 2.00;
	tInput2(3) = 1.50;

	/*!
	 * compute and collect nodal sphere values
	 *
	 * \code{.cpp}
	 *
	 * uint tNumNodes = tMesh3DHexs->get_num_entities(EntityRank::NODE);
	 * Matrix< DDRMat > tNodeSphere1Vals(1,tNumNodes);
	 * Matrix< DDRMat > tNodeSphere2Vals(1,tNumNodes);
	 *
	 * for(uint i=0; i<tNumNodes; i++)
	 * {
	 * 		Matrix< DDRMat > tNodeCoord = tMesh3DHexs->get_node_coordinate(i);
	 * 		tNodeSphere1Vals(i) = sphere_function(tNodeCoord,tInput1);
	 * 		tNodeSphere2Vals(i) = sphere_function(tNodeCoord,tInput2);
	 * }
	 * \endcode
	 *
	 */

	uint tNumNodes = tMesh3DHexs->get_num_entities(EntityRank::NODE);
	Matrix< DDRMat > tNodeSphere1Vals(1,tNumNodes);
	Matrix< DDRMat > tNodeSphere2Vals(1,tNumNodes);

	for(uint i=0; i<tNumNodes; i++)
	{
		Matrix< DDRMat > tNodeCoord = tMesh3DHexs->get_node_coordinate(i);

		tNodeSphere1Vals(i) = sphere_function(tNodeCoord,tInput1);

		tNodeSphere2Vals(i) = sphere_function(tNodeCoord,tInput2);
	}

	/*!
	 * add nodal sphere values to mesh
	 *
	 * \code{.cpp}
	 * tMesh3DHexs->add_mesh_field_real_scalar_data_loc_inds(tSphere1FieldName, EntityRank::NODE, tNodeSphere1Vals);
	 * tMesh3DHexs->add_mesh_field_real_scalar_data_loc_inds(tSphere2FieldName, EntityRank::NODE, tNodeSphere2Vals);
	 * \endcode
	 *
	 * @image html ./figures/tutorial_1_firstGeometry.png "First Sphere"
	 *
	 * @image html ./figures/tutorial_1_bothGeometries.png "Both Spheres"
	 *
	 */

	tMesh3DHexs->add_mesh_field_real_scalar_data_loc_inds(tSphere1FieldName, EntityRank::NODE, tNodeSphere1Vals);
	tMesh3DHexs->add_mesh_field_real_scalar_data_loc_inds(tSphere2FieldName, EntityRank::NODE, tNodeSphere2Vals);

	//------------------------------------------------------------------------------
	/*!
	 * <b> Step 3: setup geometry engine </b>
	 */

	/*!
	 * create the factory object, define flagging method and level set function
	 *
	 * \code{.cpp}
	 * Ge_Factory tFactory;
	 *
	 * Geometry* type0 = tFactory.pick_flag(flagType::Analytical);
	 * type0->set_analytical_function(sphere_function);
	 *
	 * Geometry* type1 = tFactory.pick_flag(flagType::Analytical);
	 * type1->set_analytical_function(sphere_function);
	 * \endcode
	 *
	 */

	Ge_Factory tFactory;
	Geometry* type0 = tFactory.set_geometry_type(type::ANALYTIC);
	type0->set_analytical_function(sphere_function);

	Geometry* type1 = tFactory.set_geometry_type(type::ANALYTIC);
	type1->set_analytical_function(sphere_function);

	/*!
	 * create the geometry engine object, set the geometries and mesh
	 *
	 * \code{.cpp}
	 *
	 * GE geometryEngine;
	 *
	 * geometryEngine.set_geometry( type0 );
	 * geometryEngine.set_geometry( type1 );
	 *
	 * geometryEngine.set_mesh( tMesh3DHexs );
	 * \endcode
	 *
	 * note: the geometry engine allows for multiple geometries and meshs to be set,
	 * for this tutorial, both of the geometries are the same (sphere) and there is
	 * only one mesh (tMesh3DHexs)
	 *
	 */

	GE geometryEngine;

	geometryEngine.set_geometry( type0 );
	geometryEngine.set_geometry( type1 );

	geometryEngine.set_mesh( tMesh3DHexs );

	/*!
	 * now that the geometries and mesh are defined, we will set the threshold
	 * for the level set function
	 *
	 * in most cases, this threshold is set to zero
	 *
	 * however, the geometry engine can set a different threshold for each geometry
	 * if desired
	 *
	 * for this tutorial, we will use a single threshold for both geometries with a
	 * value of zero
	 *
	 * \code{.cpp}
	 *
	 * moris::Cell< double > tThreshVals(1);
	 * tThreshVals(0) = 0.0;
	 *
	 * geometryEngine.set_threshold( tThreshVals );
	 * \endcode
	 *
	 */

	moris::Cell< double > tThreshVals(1); // cell with list of threshold values
	tThreshVals(0) = 0.0;

	geometryEngine.set_threshold( tThreshVals );

	//------------------------------------------------------------------------------

	/*!
	 * create cell of flags for the elements and call the intersection function
	 *
	 * \code{.cpp}
	 *
	 * moris::Cell< uint > tFlags0;
	 * moris::Cell< uint > tFlags1;
	 *
	 * tFlags0 = geometryEngine.check_for_intersection( tInput1, 0, 0, 0 );
	 * tFlags1 = geometryEngine.check_for_intersection( tInput2, 1, 0, 0 );
	 *
	 * \endcode
	 *
	 * note:
	 * the arguments for the check_for_intersection( 1, 2, 3, 4 ) function are:
	 *
	 * 1)	level set function specific parameters
	 *
	 * 2)	which level set function from list
	 *
	 * 3)	which mesh from list
	 *
	 * 4)	which threshold from list
	 *
	 */

	moris::Cell< uint > tFlags0;
	moris::Cell< uint > tFlags1;

	tFlags0 = geometryEngine.check_for_intersection( tInput1, 0, 0, 0 );
	tFlags1 = geometryEngine.check_for_intersection( tInput2, 1, 0, 0 );

	//------------------------------------------------------------------------------

	/*!
	 * <b> Step 4: create flag field and add to mesh </b>
	 */

	/*!
	 * need to cast the cell flag values to real for STK output
	 *
	 * \code{.cpp}
	 *
	 * uint tNumElems = tMesh3DHexs->get_num_elems();
	 * Matrix< DDRMat > tElementalFlags1(1,tNumElems);
	 * Matrix< DDRMat > tElementalFlags2(1,tNumElems);
	 *
	 * for(uint t=0; t<tNumElems; t++)
	 * {
	 * 		tElementalFlags1(0,t) = (real)tFlags0(t);
	 * 		tElementalFlags2(0,t) = (real)tFlags0(t);
	 * }
	 * \endcode
	 *
	 */

	uint tNumElems = tMesh3DHexs->get_num_elems();

	Matrix< DDRMat > tElementalFlags1(1,tNumElems);
	Matrix< DDRMat > tElementalFlags2(1,tNumElems);

	for(uint t=0; t<tNumElems; t++)
	{
		tElementalFlags1(0,t) = (real)tFlags0(t);
		tElementalFlags2(0,t) = (real)tFlags1(t);
	}

	/*!
	 * create the total flag field from both geomtries and add to mesh
	 *
	 * \code{.cpp}
	 *
	 * Matrix< DDRMat > tElementalFlags_total(1,tNumElems);
	 * tElementalFlags_total = tElementalFlags1 + tElementalFlags2;
	 *
	 * tMesh3DHexs->add_mesh_field_real_scalar_data_loc_inds(tRefineFieldName, EntityRank::ELEMENT, tElementalFlags_total);
	 * \endcode
	 *
	 * resulting in the flag field as:
	 * @image html ./figures/tutorial_1_elementFlags.png "Elemental Flag Field"
	 *
	 * in the above figure, the elements which are NOT flagged for refinement are solid
	 *
	 * the elemental flag list can now be sent to the refinement module
	 *
	 */

	Matrix< DDRMat > tElementalFlags_total(1,tNumElems);
	tElementalFlags_total = tElementalFlags1 + tElementalFlags2;

	tMesh3DHexs->add_mesh_field_real_scalar_data_loc_inds(tRefineFieldName, EntityRank::ELEMENT, tElementalFlags_total);

	/*!
	 * <b> Step 5: output exodus file </b>
	 */

	/*!
	 * output the exodus file
	 *
	 * \code{.cpp}
	 * std::string tOutputFile = "./exodusOutputFile.exo";
	 * tMesh3DHexs->create_output_mesh(tOutputFile);
	 * \endcode
	 *
	 */

	std::string tOutputFile = "./exodusOutputFile.exo";
	tMesh3DHexs->create_output_mesh(tOutputFile);

	//------------------------------------------------------------------------------

	/*!
	 * <b> Step 6: clean up your mess </b>
	 */

	/*!
	 * \code{.cpp}
	 * delete tMesh3DHexs
	 * \endcode
	 *
	 */

	delete tMesh3DHexs;

	//------------------------------------------------------------------------------
	gMorisComm.finalize();

	return 0;

}
