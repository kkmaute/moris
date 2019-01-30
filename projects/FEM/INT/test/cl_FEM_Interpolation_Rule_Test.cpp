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
#include "cl_FEM_Integration_Rule.hpp" //FEM/INT/src
#include "cl_FEM_Interpolator.hpp" //FEM/INT/src

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
    Matrix< DDRMat > tXi(1,2);  tXi(0)  = 1; tXi(1) = 0;
    Matrix< DDRMat > tTau(1,1); tTau(0) = 1;

    Interpolation_Function_Base* tSpaceIntFunction = tInterpolationRule.create_space_interpolation_function();
    Interpolation_Matrix tSpaceN = tSpaceIntFunction->create_matrix( 2, 0, 0 );
    Interpolation_Matrix tSpacedNdXi = tSpaceIntFunction->create_matrix( 2, 1, 0 );
    Interpolation_Matrix tSpaced2NdXi2 = tSpaceIntFunction->create_matrix( 2, 2, 0 );

    print(tSpaceN.matrix(),"NRow");
    print(tSpaceN.n_cols(),"NCol");

    tSpaceIntFunction->eval_N( tSpaceN, tXi );
    print(tSpaceN.matrix(),"tSpaceN");

    tSpaceIntFunction->eval_dNdXi( tSpacedNdXi, tXi );
    print(tSpacedNdXi.matrix(),"tSpacedNdXi");

    tSpaceIntFunction->eval_d2NdXi2(tSpaced2NdXi2, tXi);
    print(tSpaced2NdXi2.matrix(),"tSpaced2NdXi2");

    Interpolation_Function_Base* tTimeIntFunction = tInterpolationRule.create_time_interpolation_function();
    Interpolation_Matrix tTimeN = tTimeIntFunction->create_matrix( 1, 0, 0 );
    Interpolation_Matrix tTimedNdXi = tTimeIntFunction->create_matrix( 1, 0, 1 );
    Interpolation_Matrix tTimed2NdXi2 = tTimeIntFunction->create_matrix( 1, 0, 2 );

    tTimeIntFunction->eval_N( tTimeN, tTau );
    print(tTimeN.matrix(),"tTimeN");

    tTimeIntFunction->eval_dNdXi( tTimedNdXi, tTau );
    print(tTimedNdXi.matrix(),"tTimedNdXi");

    tTimeIntFunction->eval_d2NdXi2(tTimed2NdXi2, tTau);
    print(tTimed2NdXi2.matrix(),"tTimed2NdXi2");

    //Interpolation_Function_Base* tSpaceTimeIntFunction  = tInterpolationRule.create_space_time_interpolation_function();
    auto tSpacetimeN = trans( tSpaceN.matrix() )*tTimeN.matrix();
    print(reshape(tSpacetimeN,1,tSpaceN.n_cols()*tTimeN.n_cols()),"tSpacetimeN");

//    Interpolation_Matrix* tSpacetimeN;
//    tSpacetimeN->operator*( tSpaceN, tTimeN );
    //print(tSpacetimeN,"tSpacetimeN");
//------------------------------------------------------------------------------
}

TEST_CASE( "Interpolation_Rule2", "[moris],[fem],[InterRule2]" )
{
    moris::Matrix< DDRMat > tNodeCoords( 4, 2 );
    tNodeCoords( 0, 0 ) = 0.0; tNodeCoords( 0, 1 ) = 0.0;
    tNodeCoords( 1, 0 ) = 1.0; tNodeCoords( 1, 1 ) = 0.0;
    tNodeCoords( 2, 0 ) = 1.0; tNodeCoords( 2, 1 ) = 1.0;
    tNodeCoords( 3, 0 ) = 0.0; tNodeCoords( 3, 1 ) = 1.0;

    moris::Matrix< DDRMat > tPoint( 1, 2 );
    tPoint( 0, 0 ) = -1.0; tPoint( 0, 1 ) = 1.0;

    moris::Matrix< DDRMat > tTime( 1, 1 );
    tTime( 0, 0 ) = 1.0;

    //Building an interpolation rule
    Interpolation_Rule tInterpolationRule ( mtk::Geometry_Type::QUAD,
                                            Interpolation_Type::LAGRANGE,
                                            mtk::Interpolation_Order::LINEAR,
                                            Interpolation_Type::LAGRANGE,
                                            mtk::Interpolation_Order::LINEAR );

    Integration_Rule tIntegrationRule ( mtk::Geometry_Type::QUAD,
                                        Integration_Type::GAUSS,
                                        Integration_Order::QUAD_2x2,
                                        Integration_Type::GAUSS,
                                        Integration_Order::BAR_2);

    Interpolator tInterpolator( tNodeCoords,
                                1,
                                tInterpolationRule,
                                tInterpolationRule,
                                tIntegrationRule );

    Interpolation_Matrix tSpaceTimeN;

    tInterpolator.eval_N( tSpaceTimeN,
                          tPoint,
                          tTime );

    CHECK( equal_to( tSpaceTimeN.matrix()( 0, 0 ), 0 ) );
    CHECK( equal_to( tSpaceTimeN.matrix()( 0, 1 ), 0 ) );
    CHECK( equal_to( tSpaceTimeN.matrix()( 0, 2 ), 0 ) );
    CHECK( equal_to( tSpaceTimeN.matrix()( 0, 3 ), 0 ) );
    CHECK( equal_to( tSpaceTimeN.matrix()( 0, 4 ), 0 ) );
    CHECK( equal_to( tSpaceTimeN.matrix()( 0, 5 ), 0 ) );
    CHECK( equal_to( tSpaceTimeN.matrix()( 0, 6 ), 0 ) );
    CHECK( equal_to( tSpaceTimeN.matrix()( 0, 7 ), 1 ) );



//------------------------------------------------------------------------------
}
