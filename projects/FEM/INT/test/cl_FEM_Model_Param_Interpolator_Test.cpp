#include "catch.hpp"

//#define protected public
//#define private   public
//#include "cl_FEM_Field_Interpolator.hpp" //FEM/INT/src
//#undef protected
//#undef private

#include "cl_FEM_Model_Parameter_Interpolator.hpp" //FEM/INT/src


using namespace moris;
using namespace fem;

// a function pointer
real force_1d_func( const Matrix< DDRMat > aPhysPoint )
{
    real tUHat = 0;
    // evaluate nodal coeff at a physical point
    if ( aPhysPoint( 0 ) > 1.0 )
    {
        tUHat = std::sin( aPhysPoint( 1 ) );
    }
    return tUHat;
}

real force_2d_func( const Matrix< DDRMat > aPhysPoint )
{
    real tUHat = 0;
    // evaluate nodal coeff at a physical point
    if ( aPhysPoint( 0 ) > 1.0 )
    {
        tUHat = std::sin( aPhysPoint( 2 ) );
    }
    return tUHat;
}

TEST_CASE( "Model_Parameter_Interpolator", "[moris],[fem],[MPInterpolator]" )
{

    // define an epsilon environment
    real tEpsilon = 1E-12;

    SECTION( "Field Interpolator : Space bar2 - Time bar2" )
    {
        // space time geometry interpolator
        //------------------------------------------------------------------------------
        //create a bar2 space element
        Matrix< DDRMat > tXHat( 2, 1 );
        tXHat( 0, 0 ) =  0.0;
        tXHat( 1, 0 ) =  2.0;

        //create a bar2 time element
        Matrix< DDRMat > tTHat( 3, 1 );
        tTHat( 0 ) = 0.0;
        tTHat( 1 ) = 5.0;
        tTHat( 2 ) = 2.5;

        //create a space geometry interpolation rule
        Interpolation_Rule tGeomInterpRule( mtk::Geometry_Type::LINE,
                                            Interpolation_Type::LAGRANGE,
                                            mtk::Interpolation_Order::LINEAR,
                                            Interpolation_Type::LAGRANGE,
                                            mtk::Interpolation_Order::QUADRATIC );

        //create a space and a time geometry interpolator
        Geometry_Interpolator* tGeomInterpolator = new Geometry_Interpolator( tGeomInterpRule );

        //set the coefficients xHat, tHat
        tGeomInterpolator->set_coeff( tXHat, tTHat );

        // field interpolator
        //------------------------------------------------------------------------------
        //create a space time interpolation rule
        Interpolation_Rule tInterpolationRule ( mtk::Geometry_Type::LINE,
                                                Interpolation_Type::LAGRANGE,
                                                mtk::Interpolation_Order::LINEAR,
                                                Interpolation_Type::LAGRANGE,
                                                mtk::Interpolation_Order::QUADRATIC );

        //create a field interpolator
        uint tNumberOfFields = 1;
        Model_Parameter_Interpolator tMPInterpolator( tNumberOfFields,
                                                      tInterpolationRule,
                                                      tGeomInterpolator,
                                                      force_1d_func);
        tMPInterpolator.set_coeff();

        Matrix< DDRMat > tForceHat = tMPInterpolator.get_coeff();
        //print(tForceHat,"tForceHat");


        // clean up
        delete tGeomInterpolator;
    }

    SECTION( "Field Interpolator : Space QUAD4 - Time BAR3" )
       {
           // space time geometry interpolator
           //------------------------------------------------------------------------------
           //create a bar2 space element
           Matrix< DDRMat > tXHat( 4, 2 );
           tXHat = {{ 0.0, 0.0 },
                    { 2.0, 0.0 },
                    { 2.0, 2.0 },
                    { 0.0, 2.0 }};


           //create a bar2 time element
           Matrix< DDRMat > tTHat( 3, 1 );
           tTHat( 0 ) = 0.0;
           tTHat( 1 ) = 5.0;
           tTHat( 2 ) = 2.5;


           //create a space geometry interpolation rule
           Interpolation_Rule tGeomInterpRule( mtk::Geometry_Type::QUAD,
                                               Interpolation_Type::LAGRANGE,
                                               mtk::Interpolation_Order::LINEAR,
                                               Interpolation_Type::LAGRANGE,
                                               mtk::Interpolation_Order::QUADRATIC );

           //create a space and a time geometry interpolator
           Geometry_Interpolator* tGeomInterpolator = new Geometry_Interpolator( tGeomInterpRule );

           //set the coefficients xHat, tHat
           tGeomInterpolator->set_coeff( tXHat, tTHat );

           // field interpolator
           //------------------------------------------------------------------------------
           //create a space time interpolation rule
           Interpolation_Rule tInterpolationRule ( mtk::Geometry_Type::QUAD,
                                                   Interpolation_Type::LAGRANGE,
                                                   mtk::Interpolation_Order::LINEAR,
                                                   Interpolation_Type::LAGRANGE,
                                                   mtk::Interpolation_Order::QUADRATIC );

           //create a field interpolator
           uint tNumberOfFields = 1;
           Model_Parameter_Interpolator tMPInterpolator( tNumberOfFields,
                                                         tInterpolationRule,
                                                         tGeomInterpolator,
                                                         force_2d_func);
           tMPInterpolator.set_coeff();

           Matrix< DDRMat > tForceHat = tMPInterpolator.get_coeff();
           //print(tForceHat,"tForceHat");


           // clean up
           delete tGeomInterpolator;
       }

//------------------------------------------------------------------------------
}


