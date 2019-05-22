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
real force_func( const Matrix< DDRMat > aSpacePhysPoint,
                 const Matrix< DDRMat > aTimePhysPoint )
{
    real tUHat = 0;
    // evaluate nodal coeff at a physical point
    if ( aSpacePhysPoint( 0 ) > 1.0 )
    {
        tUHat = std::sin( aTimePhysPoint( 0 ) );
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
        fem::Mp_Type tMpType = fem::Mp_Type::TEMP_DIRICHLET;
        Model_Parameter_Interpolator tMPInterpolator( tNumberOfFields,
                                                      tInterpolationRule,
                                                      tGeomInterpolator,
                                                      tMpType,
                                                      force_func);
        tMPInterpolator.set_coeff();

        Matrix< DDRMat > tForceHat = tMPInterpolator.get_coeff();

        // expected coeff values for check
        Matrix< DDRMat > tForceHatExpected = {{+0.000000000000000e+00},
                                              {+0.000000000000000e+00},
                                              {+0.000000000000000e+00},
                                              {-9.589242746631385e-01},
                                              {+0.000000000000000e+00},
                                              {+5.984721441039565e-01}};

        // boolean to check coeff values
        bool tCoeffCheck = true;

        // tol to check coeff values
        real tEpsilon = 1E-12;

        // loop over coeff values
        for( uint iCheck = 0; iCheck < tForceHat.numel(); iCheck++)
        {
            // check coeff values
            tCoeffCheck = tCoeffCheck && ( std::abs( tForceHatExpected( iCheck ) - tForceHat( iCheck ) ) < tEpsilon );
        }
        REQUIRE( tCoeffCheck );

        Matrix< DDRMat > tForceHatGradDof = tMPInterpolator.val_gradDof( MSI::Dof_Type::UX );
        print(tForceHatGradDof, "tForceHatGradDof");

        // clean up
        delete tGeomInterpolator;
    }

    SECTION( "Field Interpolator : Space QUAD4 - Time BAR4" )
       {
           // space time geometry interpolator
           //------------------------------------------------------------------------------
           //create a bar2 space element
           Matrix< DDRMat > tXHat = {{ 0.0, 0.0 },
                                     { 2.0, 0.0 },
                                     { 2.0, 2.0 },
                                     { 0.0, 2.0 }};

           //create a bar2 time element
           Matrix< DDRMat > tTHat( 4, 1 );
           tTHat( 0 ) = 0.0;
           tTHat( 1 ) = 5.0;
           tTHat( 2 ) = 5.0/3.0;
           tTHat( 3 ) = 2.0*5.0/3.0;


           //create a space geometry interpolation rule
           Interpolation_Rule tGeomInterpRule( mtk::Geometry_Type::QUAD,
                                               Interpolation_Type::LAGRANGE,
                                               mtk::Interpolation_Order::LINEAR,
                                               Interpolation_Type::LAGRANGE,
                                               mtk::Interpolation_Order::CUBIC );

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
                                                   mtk::Interpolation_Order::CUBIC );

           //create a field interpolator
           uint tNumberOfFields = 1;
           fem::Mp_Type tMpType = fem::Mp_Type::TEMP_DIRICHLET;
           Model_Parameter_Interpolator tMPInterpolator( tNumberOfFields,
                                                         tInterpolationRule,
                                                         tGeomInterpolator,
                                                         tMpType,
                                                         force_func);
           tMPInterpolator.set_coeff();

           Matrix< DDRMat > tForceHat = tMPInterpolator.get_coeff();

           // expected coeff values for check
           Matrix< DDRMat > tForceHatExpected = {{ +0.000000000000000e+00 }, { +0.000000000000000e+00 }, { +0.000000000000000e+00 }, { +0.000000000000000e+00 },
                                                 { +0.000000000000000e+00 }, { -9.589242746631385e-01 }, { -9.589242746631385e-01 }, { +0.000000000000000e+00 },
                                                 { +0.000000000000000e+00 }, { +9.954079577517649e-01 }, { +9.954079577517649e-01 }, { +0.000000000000000e+00 },
                                                 { +0.000000000000000e+00 }, { -1.905679628754854e-01 }, { -1.905679628754854e-01 }, { +0.000000000000000e+00 }};

           // boolean to check coeff values
           bool tCoeffCheck = true;

           // tol to check coeff values
           real tEpsilon = 1E-12;

           // loop over coeff values
           for( uint iCheck = 0; iCheck < tForceHat.numel(); iCheck++)
           {
               // check coeff values
               tCoeffCheck = tCoeffCheck && ( std::abs( tForceHatExpected( iCheck ) - tForceHat( iCheck ) ) < tEpsilon );
           }
           REQUIRE( tCoeffCheck );

           Matrix< DDRMat > tForceHatGradDof = tMPInterpolator.val_gradDof( MSI::Dof_Type::UX );
           print(tForceHatGradDof, "tForceHatGradDof");

           // clean up
           delete tGeomInterpolator;
       }

//------------------------------------------------------------------------------
}


