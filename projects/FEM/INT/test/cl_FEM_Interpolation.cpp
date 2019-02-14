#include <string>

#include <catch.hpp>
#include "fn_equal_to.hpp"

#include "typedefs.hpp" //MRS/COR/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_FEM_Interpolation_Function.hpp"
#include "cl_FEM_Interpolation_Rule.hpp"

#include <ctime>

using namespace moris;
using namespace fem;

//TEST_CASE( "Interpolation", "[moris],[fem],[Interpolation2]" )
//{
//	mtk::Vertex* tNode0 = new NodeProxy(0.0, 0.0, 0);
//	mtk::Vertex* tNode1 = new NodeProxy(2.0, 0.0, 1);
//	mtk::Vertex* tNode2 = new NodeProxy(2.0, 2.0, 2);
//	mtk::Vertex* tNode3 = new NodeProxy(0.0, 2.0, 3);
//
//	Matrix< DDRMat > tCoords = tNode0->get_coords();
//	print(tCoords, "node coordinates");
//	CHECK( equal_to( tCoords( 0,0 ), 0.0 ) );
//
//	moris::Cell < mtk::Vertex* > tNodeList(4);
//	tNodeList(0) = tNode0;
//	tNodeList(1) = tNode1;
//	tNodeList(2) = tNode2;
//	tNodeList(3) = tNode3;
//	enum mtk::Geometry_Type tGeometryType = mtk::Geometry_Type::QUAD;
//	mtk::Cell* tElement = new ElementProxy(tNodeList,tGeometryType);
//
//	//Building an interpolation rule
////	enum fem::Interpolation_Type tSpaceInterpolationType = fem::Interpolation_Type::LAGRANGE;
////	enum mtk::Interpolation_Order tSpaceInterpolationOrder = mtk::Interpolation_Order::LINEAR;
////	enum fem::Interpolation_Type tTimeInterpolationType = fem::Interpolation_Type::LAGRANGE;
////	enum mtk::Interpolation_Order tTimeInterpolationOrder = mtk::Interpolation_Order::LINEAR;
////	Interpolation_Rule tInterpolationRule = Interpolation_Rule(tGeometryType,
////															   tSpaceInterpolationType,tSpaceInterpolationOrder,
////															   tTimeInterpolationType,tTimeInterpolationOrder);
//
//	enum fem::Interpolation_Type tSpaceTimeInterpolationType   = fem::Interpolation_Type::LAGRANGE;
//	enum mtk::Interpolation_Order tSpaceTimeInterpolationOrder = mtk::Interpolation_Order::LINEAR;
//	Interpolation_Rule tInterpolationRule = Interpolation_Rule(tGeometryType,
//															   tSpaceTimeInterpolationType,tSpaceTimeInterpolationOrder);
//	std::cout<<"Has 2 rules ? "; std::cout<<tInterpolationRule.has_two_rules()<<std::endl;
//
//	Interpolation_Function_Base* tInterpolation = tInterpolationRule.create_space_time_interpolation_function();
//	uint tNDim = tInterpolation->get_number_of_dimensions();
//	std::cout<<"Number of dimensions = "; std::cout<<tNDim<<std::endl;
//	mtk::Interpolation_Order tIntOrder = tInterpolation->get_interpolation_order();
//	std::cout<<"Interpolation order = "; std::cout<<static_cast<std::underlying_type<mtk::Interpolation_Order>::type>(tIntOrder)<<std::endl;
//	fem::Interpolation_Type tIntType = tInterpolation->get_interpolation_type();
//	std::cout<<"Interpolation type = "; std::cout<<static_cast<std::underlying_type<fem::Interpolation_Type>::type>(tIntType)<<std::endl;
//
//	//Building a FEM element
//	fem::IWG* tIWG = new IWG_L2(0.0);
//
//	moris::Cell< fem::Node_Base* > tNodeList3(4);
//	fem::Node tNode02(tNode0);
//	fem::Node tNode12(tNode1);
//	fem::Node tNode22(tNode2);
//	fem::Node tNode32(tNode3);
//	tNodeList3(0) = &tNode02;
//	tNodeList3(1) = &tNode12;
//	tNodeList3(2) = &tNode22;
//	tNodeList3(3) = &tNode32;
//
//	fem::Element* tFEMElement = new Element(tElement, tIWG, tNodeList3);
//
//	//Building a geometry interpolator
//	//Geometry_Interpolator tGeomInterpolator = Geometry_Interpolator(tFEMElement, tInterpolationRule);
//	Geometry_Interpolator tGeomInterpolator = Geometry_Interpolator(tInterpolationRule);
//	uint tNBasis = tGeomInterpolator.get_number_of_basis();
//	std::cout<<"Number of basis = "; std::cout<<tNBasis<<std::endl;
//	//Setting the evaluation point (xi_1,xi_2)
//	Matrix< DDRMat > tXi(1,2); tXi(0)=1.0; tXi(1)=1.0;
//	print(tXi,"Xi");
//	//Evaluation of the shape function N(xi_1,xi_2)
//	fem::Interpolation_Matrix tN(2,1,1,4);
//    tGeomInterpolator.eval_N(tN,tXi);
//    print(tN.matrix(),"tN");
//	//Evaluation of the 1st derivatives of the shape functions wrt xi dNdXi(xi_1,xi_2)
//    fem::Interpolation_Matrix tdNdXi(2,1,2,4);
//    tGeomInterpolator.eval_dNdXi(tdNdXi,tXi);
//    print(tdNdXi.matrix(),"tdNdXi");
//	//Evaluation of the 2nd derivatives of the shape functions wrt xi d2NdXi2(xi_1,xi_2)
//	fem::Interpolation_Matrix td2NdXi2(2,1,3,4);
//	tGeomInterpolator.eval_d2NdXi2(td2NdXi2,tXi);
//	print(td2NdXi2.matrix(),"td2NdXi2");
//	//Evaluation of the Jacobian(xi_1,xi_2)
//	Matrix< DDRMat > tJt;
//	Matrix< DDRMat > tXhat(4,2);
//	tXhat(0,0) = tNode0->get_coords()(0); tXhat(0,1) = tNode0->get_coords()(1);
//	tXhat(1,0) = tNode1->get_coords()(0); tXhat(1,1) = tNode1->get_coords()(1);
//	tXhat(2,0) = tNode2->get_coords()(0); tXhat(2,1) = tNode2->get_coords()(1);
//	tXhat(3,0) = tNode3->get_coords()(0); tXhat(3,1) = tNode3->get_coords()(1);
//	print(tXhat,"tXhat");
//	tGeomInterpolator.eval_jacobian(tJt,tdNdXi,tXhat);
//	print(tJt,"tJt");
//	//Evaluation of the Jacobian and matrices for 2nd order derivatives(xi_1,xi_2)
//	Matrix< DDRMat > tJtbis;
//	Matrix< DDRMat > tKt; //Need to define their size, not the case for the Jacobian, Maybe because he fills tLt component per component.
//	Matrix< DDRMat > tLt(3,3);
//	tGeomInterpolator.eval_jacobian_and_matrices_for_second_derivatives(tJtbis,tKt,tLt,tdNdXi,td2NdXi2,tXhat);
//
//
////
////	//moris::Cell <mtk::Vertex*> tNodeList2 = tElement->get_vertex_pointers();
////	//Matrix< DDRMat > nodeList2Coord = tNodeList2(0)->get_coords();
////	//print(nodeList2Coord, "node coordinates 2");
////	//CHECK( equal_to( nodeList2Coord( 0,0 ), 0.0 ) );
////	//CHECK( equal_to( nodeList2Coord( 0,1 ), 0.0 ) );
////
////
////	//------------------------------------------------------------------------------
//}

TEST_CASE( "Interpolation1", "[moris],[fem],[Interpolation1]" )
{
    enum mtk::Geometry_Type tGeometryType = mtk::Geometry_Type::QUAD;
    enum fem::Interpolation_Type tSpaceTimeInterpolationType   = fem::Interpolation_Type::LAGRANGE;
    enum mtk::Interpolation_Order tSpaceTimeInterpolationOrder = mtk::Interpolation_Order::LINEAR;

    Interpolation_Rule tInterpolationRule = Interpolation_Rule(tGeometryType, tSpaceTimeInterpolationType,tSpaceTimeInterpolationOrder);
    uint tNDim =0;
    auto start = std::chrono::high_resolution_clock::now();

    for ( moris::uint Ik =0; Ik <= 100000; Ik++)
    {
         Interpolation_Function_Base* tInterpolation = tInterpolationRule.create_space_time_interpolation_function();

          tNDim = tInterpolation->get_number_of_dimensions();

         delete(tInterpolation);
    }
    auto finish = std::chrono::high_resolution_clock::now();
    std::cout<<"Number of dimensions = "<<tNDim<<std::endl;

    auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>(finish-start);
    std::cout << microseconds.count() << "µs\n"<<std::endl;
    //-------------------------------------------------------------------------------------------------
    start = std::chrono::high_resolution_clock::now();

    for ( moris::uint Ik =0; Ik <= 1000000; Ik++)
    {
         Interpolation_Function< Interpolation_Type::LAGRANGE, 2, 4 >* tInterpolation = new Interpolation_Function< Interpolation_Type::LAGRANGE, 2, 4 >();

         tNDim = tInterpolation->get_number_of_dimensions();

         delete( tInterpolation );
    }
    finish = std::chrono::high_resolution_clock::now();
    std::cout<<"Number of dimensions = "<<tNDim<<std::endl;

    microseconds = std::chrono::duration_cast<std::chrono::microseconds>(finish-start);
    std::cout << microseconds.count() << "µs\n"<<std::endl;;
}

