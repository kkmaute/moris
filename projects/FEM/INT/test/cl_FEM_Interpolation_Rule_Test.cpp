#include <string>

#include <catch.hpp>
#include "fn_equal_to.hpp"

#include "typedefs.hpp" //MRS/COR/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_FEM_Node_Proxy.hpp"
#include "cl_FEM_Element_Proxy.hpp"
#include "cl_MTK_Enums.hpp" //MTK/src

#include "cl_FEM_Enums.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Rule_Bis.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Base.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Matrix.hpp" //FEM/INT/src

using namespace moris;
using namespace fem;

TEST_CASE( "Interpolation_Rule", "[moris],[fem],[InterRule]" )
{

	//Building an interpolation rule
	Interpolation_Rule_Bis tInterpolationRule = Interpolation_Rule_Bis(mtk::Geometry_Type::QUAD,
																	   Interpolation_Type::LAGRANGE,
																	   mtk::Interpolation_Order::LINEAR,
																	   Interpolation_Type::LAGRANGE,
																	   mtk::Interpolation_Order::LINEAR);
	//Sample point in space Xi
	Matrix< DDRMat > tXi(1,2);  tXi(0)  = 1; tXi(1) = 1;
	Matrix< DDRMat > tTau(1,1); tTau(0) = 1;

	Interpolation_Function_Base* tSpaceIntFunction = tInterpolationRule.create_space_interpolation_function();
	auto tSpaceN = tSpaceIntFunction->create_matrix( 2, 0, 0 );
	auto tSpacedNdXi = tSpaceIntFunction->create_matrix( 2, 1, 0 );
	auto tSpaced2NdXi2 = tSpaceIntFunction->create_matrix( 2, 2, 0 );
	print(tSpaceN.n_rows(),"NRow");
	print(tSpaceN.n_cols(),"NCol");
	tSpaceIntFunction->eval_N( tSpaceN, tXi );
	print(tSpaceN.matrix(),"tSpaceN");
	tSpaceIntFunction->eval_dNdXi( tSpacedNdXi, tXi );
	print(tSpacedNdXi.matrix(),"tSpacedNdXi");
	tSpaceIntFunction->eval_d2NdXi2(tSpaced2NdXi2, tXi);
	print(tSpaced2NdXi2.matrix(),"tSpaced2NdXi2");

	Interpolation_Function_Base* tTimeIntFunction = tInterpolationRule.create_time_interpolation_function();
	auto tTimeN = tTimeIntFunction->create_matrix( 1, 0, 0 );
	auto tTimedNdXi = tTimeIntFunction->create_matrix( 1, 0, 1 );
	auto tTimed2NdXi2 = tTimeIntFunction->create_matrix( 1, 0, 2 );
	tTimeIntFunction->eval_N( tTimeN, tTau );
	print(tTimeN.matrix(),"tTimeN");
	tTimeIntFunction->eval_dNdXi( tTimedNdXi, tTau );
	print(tTimedNdXi.matrix(),"tTimedNdXi");
	tTimeIntFunction->eval_d2NdXi2(tTimed2NdXi2, tTau);
	print(tTimed2NdXi2.matrix(),"tTimed2NdXi2");

	//Interpolation_Function_Base* tSpaceTimeIntFunction  = tInterpolationRule.create_space_time_interpolation_function();
	auto tSpacetimeN = trans( tSpaceN.matrix() )*tTimeN.matrix();
	print(reshape(tSpacetimeN,1,tSpaceN.n_cols()*tTimeN.n_cols()),"tSpacetimeN");

//	Interpolation_Matrix* tSpacetimeN;
//	tSpacetimeN->operator*( tSpaceN, tTimeN );
	//print(tSpacetimeN,"tSpacetimeN");
//------------------------------------------------------------------------------
}
