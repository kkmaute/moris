/*
 * cl_GE_test_discrete.cpp
 *
 *  Created on: Apr 22, 2019
 *      Author: sonne
 */

#include "../../src/cl_GE_Core.hpp"
#include "catch.hpp"

//------------------------------------------------------------------------------
// GE includes
#include "cl_GE_Element.hpp"
#include "cl_GE_Factory.hpp"
#include "cl_GE_Node.hpp"

// linalg includes
#include "cl_Matrix.hpp"
#include "fn_all_true.hpp"
#include "fn_equal_to.hpp"
#include "linalg_typedefs.hpp"
#include "op_equal_equal.hpp"

// MTK includes
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Vertex.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Tools.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Scalar_Field_Info.hpp"
//------------------------------------------------------------------------------

using namespace moris;
using namespace ge;

TEST_CASE("discrete_functionalities_test","[GE],[discrete_functionalities]")
{
    if(par_size()<=1)
    {
        /*
         * 1) create single-element mesh
         * 2) discretize a circle function onto the mesh
         * 3) determine LS values and sensitivity at nodes
         * 4) determine intersection locations along edges
         *
         *         [3]
         *   (0,1)      (1,1)
         *      x--------x
         *      |        |
         *      O        |
         * [4]  |        |  [2]
         *      |        |
         *      x-----O--x
         *   (0,0)      (1,0)
         *         [1]
         *
         */
        Matrix< DDRMat > tTMatrix( 4,4, 0.0 ); //T-matrix to be used for L2 projection
        tTMatrix(0,0) = 0.25;
        tTMatrix(1,1) = 0.25;
        tTMatrix(2,2) = 0.25;
        tTMatrix(3,3) = 0.25;
        /* *************************************************************************
         * create the 2D mesh with a single element and initialize scalar node field
         * *************************************************************************
         */
        uint aNumElemTypes = 1;     // quad
        uint aNumDim = 2;           // specify number of spatial dimensions

        Matrix< IdMat > aElementConnQuad = {{ 1, 2, 3, 4 }};   // specify element connectivity of quad for mesh

        Matrix< IdMat > aElemLocalToGlobalQuad = {{ 1 }};      // specify the local to global element map for quads

        Matrix< DDRMat > aCoords = {{ 0.0, 0.0 },
                                    { 1.0, 0.0 },
                                    { 1.0, 1.0 },
                                    { 0.0, 1.0 }};             // Node coordinate matrix

        Matrix< IdMat > aNodeLocalToGlobal = {{ 1, 2, 3, 4 }}; // specify the local to global map
        //------------------------------------------------------------------------------
        // create MORIS mesh using MTK database
        mtk::MtkMeshData tMeshData( aNumElemTypes );
        tMeshData.CreateAllEdgesAndFaces = true;
        tMeshData.SpatialDim = & aNumDim;
        tMeshData.ElemConn(0) = & aElementConnQuad;
        tMeshData.NodeCoords = & aCoords;
        tMeshData.LocaltoGlobalElemMap(0) = & aElemLocalToGlobalQuad;
        tMeshData.LocaltoGlobalNodeMap = & aNodeLocalToGlobal;
        //------------------------------------------------------------------------------
        // declare scalar node field for the circle LS
        moris::mtk::Scalar_Field_Info<DDRMat> tNodeField1;
        std::string tFieldName = "circle";
        tNodeField1.set_field_name(tFieldName);
        tNodeField1.set_field_entity_rank(EntityRank::NODE);

        moris::mtk::Scalar_Field_Info<DDRMat> tNodeField2;
        std::string tFieldName1 = "projectionVals";
        tNodeField2.set_field_name(tFieldName1);
        tNodeField2.set_field_entity_rank(EntityRank::NODE);

        // initialize field information container
        moris::mtk::MtkFieldsInfo tFieldsInfo;
        // Place the node field into the field info container
        add_field_for_mesh_input(&tNodeField1,tFieldsInfo);
        add_field_for_mesh_input(&tNodeField2,tFieldsInfo);

        // declare some supplementary fields
        tMeshData.FieldsInfo = &tFieldsInfo;
        //------------------------------------------------------------------------------

        mtk::Mesh* tMesh2D_Quad4 = create_mesh( MeshType::STK, tMeshData );

        /* **************************************
         * discretize field onto mesh
         * **************************************
         */
        // input parameters for circle
        moris::Cell< real > tInputs(3);
        tInputs(0) = 0.0;   // x center
        tInputs(1) = 0.0;   // y center
        tInputs(2) = 0.6;   // radius

        // compute nodal circle values for mesh
        uint tNumNodes = tMesh2D_Quad4->get_num_entities(EntityRank::NODE);
        Matrix< DDRMat > tNodeVals(1,tNumNodes);

        // collect nodal circle values
        for(uint i=0; i<tNumNodes; i++)
        {
            tNodeVals(i) = circle_function( tMesh2D_Quad4->get_node_coordinate(i), tInputs );
        }
        // add nodal circle values to mesh
        tMesh2D_Quad4->add_mesh_field_real_scalar_data_loc_inds(tFieldName, EntityRank::NODE, tNodeVals);
        //------------------------------------------------------------------------------

        /* **************************************
         * build GE and ask questions
         * **************************************
         */
        Ge_Factory tFactory;
        std::shared_ptr< Geometry > circle = tFactory.set_geometry_type(type::DISCRETE);

        Cell<std::string> tFields(1);   // cell of field names
        tFields(0) = tFieldName;

        circle->set_member_variables(tMesh2D_Quad4, tFields);

        circle->set_mesh_and_t_matrix(tMesh2D_Quad4, tTMatrix);

        GE_Core tGeometryEngine;                    // create geometry engine and set geometry
        tGeometryEngine.set_geometry( circle );

        // determine LS values at the nodes
        Matrix< DDRMat > tLSVals(4,1,0.0);

        for(uint i=0; i<tNumNodes; i++)
        {
            tLSVals(i,0) = tGeometryEngine.get_geometry_pointer(0)->access_field_value_with_entity_index(i,EntityRank::NODE);
        }
        /*
         * *****************************************
         * check field values
         * *****************************************
         */
        CHECK( equal_to( tLSVals( 0,0 ), -0.36 ) );
        CHECK( equal_to( tLSVals( 1,0 ),  0.64 ) );
        CHECK( equal_to( tLSVals( 2,0 ),  1.64 ) );
        CHECK( equal_to( tLSVals( 3,0 ),  0.64 ) );



        Matrix< DDRMat > tADVs = tGeometryEngine.compute_nodal_advs( 0,
                                                                     tInputs,
                                                                     mtk::Geometry_Type::QUAD,
                                                                     fem::Integration_Type::GAUSS,
                                                                     fem::Integration_Order::QUAD_2x2,
                                                                     fem::Interpolation_Type::LAGRANGE,
                                                                     mtk::Interpolation_Order::LINEAR );
        Matrix< DDRMat > tProjectedVals = tTMatrix*tADVs;
        tMesh2D_Quad4->add_mesh_field_real_scalar_data_loc_inds(tFieldName1, EntityRank::NODE, tProjectedVals);

        /*
         * *****************************************
         * determine intersection location
         * *****************************************
         */
        // edge (1)




        //------------------------------------------------------------------------------
//        std::string tOutputFile = "./discrete_functionalities.exo";
//        tMesh2D_Quad4->create_output_mesh(tOutputFile);
        delete tMesh2D_Quad4;
    }
}
//------------------------------------------------------------------------------
