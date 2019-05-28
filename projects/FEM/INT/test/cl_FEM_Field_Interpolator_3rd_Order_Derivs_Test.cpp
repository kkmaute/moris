#include "catch.hpp"

#define protected public
#define private   public
#include "cl_FEM_Field_Interpolator.hpp" //FEM/INT/src
#undef protected
#undef private

using namespace moris;
using namespace fem;

TEST_CASE( "Field_Interpolator_Derivatives", "[moris],[fem],[FieldInterpolator_Derivatives]" )
{

    // define an epsilon environment
    double tEpsilon = 1E-12;

    SECTION( "Field Interpolator : 2D space - 3rd Derivatives" )
    {
        // space and time geometry interpolator
        //------------------------------------------------------------------------------
        // create a QUAD9 space element with domain (x,y) e [-3,-1]x[-1,3]
        Matrix< DDRMat > tXHat( 9, 2 );
        tXHat( 0, 0 ) = -1.0;   tXHat( 0, 1 ) = -1.0;
        tXHat( 1, 0 ) = -3.0;   tXHat( 1, 1 ) = -1.0;
        tXHat( 2, 0 ) = -3.0;   tXHat( 2, 1 ) =  3.0;
        tXHat( 3, 0 ) = -1.0;   tXHat( 3, 1 ) =  3.0;
        tXHat( 4, 0 ) = -2.0;   tXHat( 4, 1 ) = -1.0;
        tXHat( 5, 0 ) = -3.0;   tXHat( 5, 1 ) =  1.0;
        tXHat( 6, 0 ) = -2.0;   tXHat( 6, 1 ) =  3.0;
        tXHat( 7, 0 ) = -1.0;   tXHat( 7, 1 ) =  1.0;
        tXHat( 8, 0 ) = -2.0;   tXHat( 8, 1 ) =  1.0;

        //create a line time element
        Matrix< DDRMat > tTHat( 2, 1 );
        tTHat( 0 ) = 0.0;
        tTHat( 1 ) = 5.0;

        // create a space and time geometry interpolation rule
        Interpolation_Rule tGeomInterpRule( mtk::Geometry_Type::QUAD,
                                            Interpolation_Type::LAGRANGE,
                                            mtk::Interpolation_Order::QUADRATIC,
                                            Interpolation_Type::LAGRANGE,
                                            mtk::Interpolation_Order::LINEAR );

        // create a space and time geometry interpolator
        Geometry_Interpolator tGeomInterpolator( tGeomInterpRule, true );

        //set the coefficients xHat, tHat
        tGeomInterpolator.set_coeff( tXHat, tTHat );

        // field interpolator
        //------------------------------------------------------------------------------
        //create a space time interpolation rule
        Interpolation_Rule tInterpolationRule ( mtk::Geometry_Type::QUAD,
                                                Interpolation_Type::LAGRANGE,
                                                mtk::Interpolation_Order::QUADRATIC,
                                                Interpolation_Type::LAGRANGE,
                                                mtk::Interpolation_Order::LINEAR );

        //create a field interpolator
        uint tNumberOfFields = 1;
        Field_Interpolator tFieldInterpolator( tNumberOfFields,
                                               tInterpolationRule,
                                               tGeomInterpolator );











        // put values of field f(x,y) = x^2*y^2 - 2*x^2*y on nodes
        Matrix< DDRMat > tUHat( 9, 1 );
        tUHat( 0, 0 ) =  3.0;
        tUHat( 1, 0 ) = 27.0;
        tUHat( 2, 0 ) = 27.0;
        tUHat( 3, 0 ) =  3.0;
        tUHat( 4, 0 ) = 12.0;
        tUHat( 5, 0 ) = -9.0;
        tUHat( 6, 0 ) = 12.0;
        tUHat( 7, 0 ) = -1.0;
        tUHat( 8, 0 ) = -4.0;

        // specify test point in element: x = (-2.2, -0.2) ==> xi = (0.2, -0.6)
        Matrix< DDRMat > tXiPtest(1,2);
        tXiPtest( 0, 0 ) =   0.2;
        tXiPtest( 0, 1 ) = - 0.6;

        // set nominal values ---------------

        // first derivatives
        Matrix< DDRMat > tdPhidx_nominal(2,1);
        tdPhidx_nominal(0,0) = -  1.936;
        tdPhidx_nominal(1,0) = - 11.616;

        // second derivatives
        Matrix< DDRMat > td2Phidx2_nominal(3,1);
        td2Phidx2_nominal(0,0) =  0.88;
        td2Phidx2_nominal(1,0) =  9.68;
        td2Phidx2_nominal(2,0) = 10.56;

        // third derivatives
        Matrix< DDRMat > td3Phidx3_nominal(4,1);
        td3Phidx3_nominal(0,0) =   0.0;
        td3Phidx3_nominal(1,0) =   0.0;
        td3Phidx3_nominal(2,0) = - 4.8;
        td3Phidx3_nominal(2,0) = - 8.8;
    }


//    SECTION( "Field Interpolator : 3D space - 3rd Derivatives" )
//    {
//        // space and time geometry interpolator
//        //------------------------------------------------------------------------------
//        // create a HEX27 space element
//        Matrix< DDRMat > tXHat( 27, 3 );
//        tXHat(  0, 0 ) = ;   tXHat(  0, 1 ) = ;   tXHat(  0, 2 ) = ;
//        tXHat(  1, 0 ) = ;   tXHat(  1, 1 ) = ;   tXHat(  1, 2 ) = ;
//        tXHat(  2, 0 ) = ;   tXHat(  2, 1 ) = ;   tXHat(  2, 2 ) = ;
//        tXHat(  3, 0 ) = ;   tXHat(  3, 1 ) = ;   tXHat(  3, 2 ) = ;
//        tXHat(  4, 0 ) = ;   tXHat(  4, 1 ) = ;   tXHat(  4, 2 ) = ;
//        tXHat(  5, 0 ) = ;   tXHat(  5, 1 ) = ;   tXHat(  5, 2 ) = ;
//        tXHat(  6, 0 ) = ;   tXHat(  6, 1 ) = ;   tXHat(  6, 2 ) = ;
//        tXHat(  7, 0 ) = ;   tXHat(  7, 1 ) = ;   tXHat(  7, 2 ) = ;
//        tXHat(  8, 0 ) = ;   tXHat(  8, 1 ) = ;   tXHat(  8, 2 ) = ;
//        tXHat(  9, 0 ) = ;   tXHat(  9, 1 ) = ;   tXHat(  9, 2 ) = ;
//        tXHat( 10, 0 ) = ;   tXHat( 10, 1 ) = ;   tXHat( 10, 2 ) = ;
//        tXHat( 11, 0 ) = ;   tXHat( 11, 1 ) = ;   tXHat( 11, 2 ) = ;
//        tXHat( 12, 0 ) = ;   tXHat( 12, 1 ) = ;   tXHat( 12, 2 ) = ;
//        tXHat( 13, 0 ) = ;   tXHat( 13, 1 ) = ;   tXHat( 13, 2 ) = ;
//        tXHat( 14, 0 ) = ;   tXHat( 14, 1 ) = ;   tXHat( 14, 2 ) = ;
//        tXHat( 15, 0 ) = ;   tXHat( 15, 1 ) = ;   tXHat( 15, 2 ) = ;
//        tXHat( 16, 0 ) = ;   tXHat( 16, 1 ) = ;   tXHat( 16, 2 ) = ;
//        tXHat( 17, 0 ) = ;   tXHat( 17, 1 ) = ;   tXHat( 17, 2 ) = ;
//        tXHat( 18, 0 ) = ;   tXHat( 18, 1 ) = ;   tXHat( 18, 2 ) = ;
//        tXHat( 19, 0 ) = ;   tXHat( 19, 1 ) = ;   tXHat( 19, 2 ) = ;
//        tXHat( 20, 0 ) = ;   tXHat( 20, 1 ) = ;   tXHat( 20, 2 ) = ;
//        tXHat( 21, 0 ) = ;   tXHat( 21, 1 ) = ;   tXHat( 21, 2 ) = ;
//        tXHat( 22, 0 ) = ;   tXHat( 22, 1 ) = ;   tXHat( 22, 2 ) = ;
//        tXHat( 23, 0 ) = ;   tXHat( 23, 1 ) = ;   tXHat( 23, 2 ) = ;
//        tXHat( 24, 0 ) = ;   tXHat( 24, 1 ) = ;   tXHat( 24, 2 ) = ;
//        tXHat( 25, 0 ) = ;   tXHat( 25, 1 ) = ;   tXHat( 25, 2 ) = ;
//        tXHat( 26, 0 ) = ;   tXHat( 26, 1 ) = ;   tXHat( 26, 2 ) = ;
//
//
//
//    }


//------------------------------------------------------------------------------
}


