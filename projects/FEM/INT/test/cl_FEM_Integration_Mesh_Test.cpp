#include "catch.hpp"
#include "cl_FEM_Geometry_Interpolator.hpp" //FEM/INT/sr
#include "cl_FEM_Integrator.hpp" //FEM/INT/sr

using namespace moris;
using namespace fem;

TEST_CASE( "Intergration_Mesh", "[moris],[fem],[IntegMesh]" )
{
    // define an epsilon environment
    double tEpsilon = 1E-12;

    SECTION( "Interpolation mesh QUAD4 - Integration mesh QUAD4 " )
    {
        // define an interpolation mesh
        //------------------------------------------------------------------------------
        // define a QUAD4 space element, i.e. space coordinates xHat
        Matrix< DDRMat > tXHat = {{ 0.0, 0.0 },
                                  { 1.0, 0.0 },
                                  { 1.0, 1.0 },
                                  { 0.0, 1.0 } };

        // define a line time element, i.e. time coordinates tHat
        Matrix< DDRMat > tTHat = {{ 0.0 },
                                  { 1.0 },
                                  { 0.5 }};

        // create a space and time geometry interpolation rule
        Interpolation_Rule tGeomInterpRule( mtk::Geometry_Type::QUAD,
                                            Interpolation_Type::LAGRANGE,
                                            mtk::Interpolation_Order::LINEAR,
                                            Interpolation_Type::LAGRANGE,
                                            mtk::Interpolation_Order::QUADRATIC );

        // create a space and time geometry interpolator
        Geometry_Interpolator tGeomInterpolator( tGeomInterpRule, true );

        //set the coefficients xHat, tHat
        tGeomInterpolator.set_space_coeff( tXHat );
        tGeomInterpolator.set_time_coeff(  tTHat );

        // space and time geometry interpolations
        Interpolation_Function_Base * tSpaceInterpolation = tGeomInterpRule.create_space_interpolation_function();
        Interpolation_Function_Base * tTimeInterpolation  = tGeomInterpRule.create_time_interpolation_function();

        uint tNumSpaceBases    = tSpaceInterpolation->get_number_of_bases();
        uint tNumTimeBases     = tTimeInterpolation->get_number_of_bases();
        uint tNumSpaceDim      = tSpaceInterpolation->get_number_of_dimensions();
        uint tNumTimeDim       = tTimeInterpolation->get_number_of_dimensions();
        uint tNumParamSpaceDim = tSpaceInterpolation->get_number_of_param_dimensions();

        // define an integration mesh
        //------------------------------------------------------------------------------
        // integration mesh geometry type
        mtk::Geometry_Type tIntegGeometryType = mtk::Geometry_Type::QUAD;

        // define a QUAD4 integration element, i.e. space param coordinates xiHat
        Matrix< DDRMat > tXiHat = {{ -1.0, -1.0 },
                                   {  0.0, -1.0 },
                                   {  0.0,  1.0 },
                                   { -1.0,  1.0 }};

//        Matrix< DDRMat > tXiHat = {{ -1.0, -1.0 },
//                                   {  1.0, -1.0 },
//                                   {  1.0,  1.0 },
//                                   { -1.0,  1.0 }};

//        Matrix< DDRMat > tXiHat = {{ -0.5, -0.5 },
//                                   {  0.5, -0.5 },
//                                   {  0.5,  0.5 },
//                                   { -0.5,  0.5 }};

        // integration mesh interpolation rule
        Interpolation_Rule tIntegInterpRule( tIntegGeometryType,
                                             Interpolation_Type::LAGRANGE,
                                             mtk::Interpolation_Order::LINEAR,
                                             Interpolation_Type::LAGRANGE,
                                             mtk::Interpolation_Order::QUADRATIC );

        // space geometry interpolations for integration mesh
        Interpolation_Function_Base * tIntegSpaceInterpolation = tIntegInterpRule.create_space_interpolation_function();
        uint tIntegNumSpaceBases    = tIntegSpaceInterpolation->get_number_of_bases();
        uint tIntegNumSpaceDim      = tIntegSpaceInterpolation->get_number_of_dimensions();
        uint tIntegNumParamSpaceDim = tIntegSpaceInterpolation->get_number_of_param_dimensions();

        // create a integration rule
        Integration_Rule tIntegrationRule( tIntegGeometryType,
                                           Integration_Type::GAUSS,
                                           Integration_Order::QUAD_3x3,
                                           Integration_Type::GAUSS,
                                           Integration_Order::BAR_3 );

        // create a side integrator
        Integrator tIntegrator( tIntegrationRule );

        //get number of integration points, integration points and weights
        uint             tNumOfIntegPoints = tIntegrator.get_number_of_points();
        Matrix< DDRMat > tIntegPoints      = tIntegrator.get_points();
        Matrix< DDRMat > tIntegWeights     = tIntegrator.get_weights();

        // boolean for surface check
        bool tSurfaceCheck = true;

        // init the surface of the integration mesh
        real tSurface = 0;

        // loop over the integration points
        for ( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
        {
            // get the treated integration point location in integration space
            Matrix< DDRMat > tIntegPointI = tIntegPoints.get_column( iGP );

            // get the treated integration point location in the interpolation space
            //------------------------------------------------------------------------------
            // unpack the space and time param coords of the integration point
            Matrix< DDRMat > tXi = tIntegPointI( { 0, tIntegNumParamSpaceDim-1 }, { 0, 0 } );
            Matrix< DDRMat > tTau( 1, 1, tIntegPointI( tIntegNumParamSpaceDim ) );

            // evaluate space interpolation shape functions at integration point
            Matrix< DDRMat > tNIntegSpace = tIntegSpaceInterpolation->eval_N( tXi );

            // evaluate time interpolation shape functions at aParamPoint
            Matrix< DDRMat > tNTime = tTimeInterpolation->eval_N( tTau );

            // build space time interpolation functions for integration mesh
            Matrix< DDRMat > tN = reshape( trans( tNIntegSpace ) * tNTime, 1, tNumTimeBases*tIntegNumSpaceBases );

            // get the parametric coordinates of the integration mesh in interpolation mesh
            Matrix< DDRMat > tInterpParamCoords( tNumTimeDim+tNumParamSpaceDim, tNumTimeBases*tIntegNumSpaceBases );

            // get the time parametric coordinates
            Matrix< DDRMat > tTimeParamCoords  = tTimeInterpolation->get_param_coords();

            // get a vector of ones
            Matrix< DDRMat > tOnes( 1, tIntegNumSpaceBases, 1.0 );

            // loop on the time bases
            for( uint i = 0; i < tNumTimeBases; i++ )
            {
                // fill the space time parametric coordinates matrix with space coordinates
                tInterpParamCoords( { 0, tIntegNumParamSpaceDim-1 }, { i * tIntegNumSpaceBases, ( i + 1 ) * tIntegNumSpaceBases-1 })
                    = trans( tXiHat );

                // fill the space time parametric coordinates matrix with time coordinates
                tInterpParamCoords( { tIntegNumParamSpaceDim, tIntegNumParamSpaceDim }, { i * tIntegNumSpaceBases, ( i + 1 ) * tIntegNumSpaceBases-1 })
                    = tTimeParamCoords( i ) * tOnes;
             }

            // compute the parametric coordinates of the SideParamPoint in the parent reference element
            Matrix< DDRMat > tRefIntegPointI = tInterpParamCoords * trans( tN );

            // evaluate detJ
            //------------------------------------------------------------------------------
            real tDetJ1 = tGeomInterpolator.det_J( tRefIntegPointI );

            // get the space jacobian
            Matrix <DDRMat> tdNSpacedXi = tIntegSpaceInterpolation->eval_dNdXi( tXi );
            Matrix< DDRMat > tSpaceJt   = tdNSpacedXi * tXiHat ;
            // get the time Jacobian
            Matrix< DDRMat > tTauHat     = tTimeInterpolation->get_param_coords();
            Matrix< DDRMat > tdNTimedTau = tTimeInterpolation->eval_dNdXi( tTau );
            Matrix< DDRMat > tTimeJt     = tdNTimedTau * trans( tTauHat ) ;
            real tDetJ2 = det( tSpaceJt ) * det( tTimeJt );

            real tDetJ = tDetJ1 * tDetJ2;

            // add contribution to the surface
            tSurface = tSurface + tDetJ * tIntegWeights( iGP );
        }

        // check the surface value
        tSurfaceCheck = tSurfaceCheck && ( std::abs( tSurface - 0.5 ) < tEpsilon );
        std::cout<<tSurface<<std::endl;
        REQUIRE( tSurfaceCheck );

    }

    SECTION( "Interpolation mesh QUAD4 - Integration mesh TRI3 " )
       {
           // define an interpolation mesh
           //------------------------------------------------------------------------------
           // define a QUAD4 space element, i.e. space coordinates xHat
           Matrix< DDRMat > tXHat = {{ 0.0, 0.0 },
                                     { 1.0, 0.0 },
                                     { 1.0, 1.0 },
                                     { 0.0, 1.0 }};

           // define a line time element, i.e. time coordinates tHat
           Matrix< DDRMat > tTHat = {{ 0.0 },
                                     { 1.0 },
                                     { 0.5 }};

           // create a space and time geometry interpolation rule
           Interpolation_Rule tGeomInterpRule( mtk::Geometry_Type::QUAD,
                                               Interpolation_Type::LAGRANGE,
                                               mtk::Interpolation_Order::LINEAR,
                                               Interpolation_Type::LAGRANGE,
                                               mtk::Interpolation_Order::QUADRATIC );

           // create a space and time geometry interpolator
           Geometry_Interpolator tGeomInterpolator( tGeomInterpRule, true );

           //set the coefficients xHat, tHat
           tGeomInterpolator.set_space_coeff( tXHat );
           tGeomInterpolator.set_time_coeff(  tTHat );

           // space and time geometry interpolations
           Interpolation_Function_Base * tSpaceInterpolation = tGeomInterpRule.create_space_interpolation_function();
           Interpolation_Function_Base * tTimeInterpolation  = tGeomInterpRule.create_time_interpolation_function();

           uint tNumSpaceBases = tSpaceInterpolation->get_number_of_bases();
           uint tNumTimeBases  = tTimeInterpolation->get_number_of_bases();
           uint tNumSpaceDim   = tSpaceInterpolation->get_number_of_dimensions();
           uint tNumTimeDim    = tTimeInterpolation->get_number_of_dimensions();
           uint tNumParamSpaceDim = tSpaceInterpolation->get_number_of_param_dimensions();

           // define an integration mesh
           //------------------------------------------------------------------------------
           // integration mesh geometry type
           mtk::Geometry_Type tIntegGeometryType = mtk::Geometry_Type::TRI;

           // integration mesh interpolation rule
           Interpolation_Rule tIntegInterpRule( tIntegGeometryType,
                                                Interpolation_Type::LAGRANGE,
                                                mtk::Interpolation_Order::LINEAR,
                                                Interpolation_Type::LAGRANGE,
                                                mtk::Interpolation_Order::QUADRATIC );

           // define a QUAD4 integration element, i.e. space param coordinates xiHat
           Matrix< DDRMat > tXiHat = {{ -1.0, -1.0 },
                                      {  0.0, -1.0 },
                                      { -1.0,  1.0 }};

           // space geometry interpolations for integration mesh
           Interpolation_Function_Base * tIntegSpaceInterpolation = tIntegInterpRule.create_space_interpolation_function();
           uint tIntegNumSpaceBases    = tIntegSpaceInterpolation->get_number_of_bases();
           uint tIntegNumSpaceDim      = tIntegSpaceInterpolation->get_number_of_dimensions();
           uint tIntegNumParamSpaceDim = tIntegSpaceInterpolation->get_number_of_param_dimensions();

           // create a integration rule
           Integration_Rule tIntegrationRule( tIntegGeometryType,
                                              Integration_Type::GAUSS,
                                              Integration_Order::TRI_3,
                                              Integration_Type::GAUSS,
                                              Integration_Order::BAR_3 );

           // create a side integrator
           Integrator tIntegrator( tIntegrationRule );

           //get number of integration points, integration points and weights
           uint             tNumOfIntegPoints = tIntegrator.get_number_of_points();
           Matrix< DDRMat > tIntegPoints      = tIntegrator.get_points();
           Matrix< DDRMat > tIntegWeights     = tIntegrator.get_weights();

           // boolean for surface check
           bool tSurfaceCheck = true;

           // init the surface of the integration mesh
           real tSurface = 0;

           // loop over the integration points
           for ( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
           {
               // get the treated integration point location in integration space
               Matrix< DDRMat > tIntegPointI = tIntegPoints.get_column( iGP );

               // get the treated integration point location in the interpolation space
               //------------------------------------------------------------------------------
               // unpack the space and time param coords of the integration point
               Matrix< DDRMat > tXi = tIntegPointI( { 0, tIntegNumParamSpaceDim-1 }, { 0, 0 } );
               Matrix< DDRMat > tTau( 1, 1, tIntegPointI( tIntegNumParamSpaceDim ) );
               //print(tXi,"tXi");
               //print(tTau,"tTau");

               // evaluate space interpolation shape functions at integration point
               Matrix< DDRMat > tNIntegSpace = tIntegSpaceInterpolation->eval_N( tXi );

               // evaluate time interpolation shape functions at aParamPoint
               Matrix< DDRMat > tNTime = tTimeInterpolation->eval_N( tTau );

               // build space time interpolation functions for integration mesh
               Matrix< DDRMat > tN = reshape( trans( tNIntegSpace ) * tNTime, 1, tNumTimeBases*tIntegNumSpaceBases );

               // get the parametric coordinates of the integration mesh in interpolation mesh
               Matrix< DDRMat > tInterpParamCoords( tNumTimeDim+tNumParamSpaceDim, tNumTimeBases*tIntegNumSpaceBases );

               // get the time parametric coordinates
               Matrix< DDRMat > tTimeParamCoords  = tTimeInterpolation->get_param_coords();

               // get a vector of ones
               Matrix< DDRMat > tOnes( 1, tIntegNumSpaceBases, 1.0 );

               // loop on the time bases
               for( uint i = 0; i < tNumTimeBases; i++ )
               {
                   // fill the space time parametric coordinates matrix with space coordinates
                   tInterpParamCoords( { 0, tNumParamSpaceDim-1 }, { i * tIntegNumSpaceBases, ( i + 1 ) * tIntegNumSpaceBases-1 })
                       = trans( tXiHat );

                   // fill the space time parametric coordinates matrix with time coordinates
                   tInterpParamCoords( { tNumParamSpaceDim, tNumParamSpaceDim }, { i * tIntegNumSpaceBases, ( i + 1 ) * tIntegNumSpaceBases-1 })
                       = tTimeParamCoords( i ) * tOnes;
                }
                //print(tInterpParamCoords,"tInterpParamCoords");

                // compute the parametric coordinates of the SideParamPoint in the parent reference element
                Matrix< DDRMat > tRefIntegPointI = tInterpParamCoords * trans( tN );

                // evaluate detJ
                //------------------------------------------------------------------------------
                real tDetJ1 = tGeomInterpolator.det_J( tRefIntegPointI );

                // get the space jacobian
                Matrix <DDRMat> tdNSpacedXi = tIntegSpaceInterpolation->eval_dNdXi( tXi );
                Matrix< DDRMat > tSpaceJt   = tdNSpacedXi * tXiHat;
                // get the time Jacobian
                Matrix< DDRMat > tTauHat     = tTimeInterpolation->get_param_coords();
                Matrix< DDRMat > tdNTimedTau = tTimeInterpolation->eval_dNdXi( tTau );
                Matrix< DDRMat > tTimeJt     = tdNTimedTau * trans( tTauHat );

                Matrix< DDRMat > tSpaceJt2( tIntegNumParamSpaceDim, tIntegNumParamSpaceDim, 1.0 );
                tSpaceJt2({ 1, tIntegNumParamSpaceDim-1 },{ 0, tIntegNumParamSpaceDim-1 }) = trans( tSpaceJt );
                real detJSpace = det( tSpaceJt2 ) / 2.0;
                real tDetJ2 = detJSpace * det( tTimeJt );

                real tDetJ = tDetJ1 * tDetJ2;

                // add contribution to the surface
                tSurface = tSurface + tDetJ * tIntegWeights( iGP );
            }

            // check the surface value
            tSurfaceCheck = tSurfaceCheck && ( std::abs( tSurface - 0.25 ) < tEpsilon );
            std::cout<<tSurface<<std::endl;
            REQUIRE( tSurfaceCheck );

       }


    SECTION( "Interpolation mesh HEX8 - Integration mesh HEX8 " )
    {
        // define an interpolation mesh
        //------------------------------------------------------------------------------
        // define a HEX8 space element, i.e. space coordinates xHat
        Matrix< DDRMat > tXHat = {{ 0.0, 0.0, 0.0 },
                                  { 1.0, 0.0, 0.0 },
                                  { 1.0, 1.0, 0.0 },
                                  { 0.0, 1.0, 0.0 },
                                  { 0.0, 0.0, 1.0 },
                                  { 1.0, 0.0, 1.0 },
                                  { 1.0, 1.0, 1.0 },
                                  { 0.0, 1.0, 1.0 }};

        // define a line time element, i.e. time coordinates tHat
        Matrix< DDRMat > tTHat = {{ 0.0 },
                                  { 1.0 },
                                  { 0.5 }};

        // create a space and time geometry interpolation rule
        Interpolation_Rule tGeomInterpRule( mtk::Geometry_Type::HEX,
                                            Interpolation_Type::LAGRANGE,
                                            mtk::Interpolation_Order::LINEAR,
                                            Interpolation_Type::LAGRANGE,
                                            mtk::Interpolation_Order::QUADRATIC );

        // create a space and time geometry interpolator
        Geometry_Interpolator tGeomInterpolator( tGeomInterpRule, true );

        //set the coefficients xHat, tHat
        tGeomInterpolator.set_space_coeff( tXHat );
        tGeomInterpolator.set_time_coeff(  tTHat );

        // space and time geometry interpolations
        Interpolation_Function_Base * tSpaceInterpolation = tGeomInterpRule.create_space_interpolation_function();
        Interpolation_Function_Base * tTimeInterpolation  = tGeomInterpRule.create_time_interpolation_function();

        uint tNumSpaceBases    = tSpaceInterpolation->get_number_of_bases();
        uint tNumTimeBases     = tTimeInterpolation->get_number_of_bases();
        uint tNumSpaceDim      = tSpaceInterpolation->get_number_of_dimensions();
        uint tNumTimeDim       = tTimeInterpolation->get_number_of_dimensions();
        uint tNumParamSpaceDim = tSpaceInterpolation->get_number_of_param_dimensions();

        // define an integration mesh
        //------------------------------------------------------------------------------
        // integration mesh geometry type
        mtk::Geometry_Type tIntegGeometryType = mtk::Geometry_Type::HEX;

        // define a HEX8 integration element, i.e. space param coordinates xiHat
        Matrix< DDRMat > tXiHat = {{ -1.0, -1.0, -1.0 },
                                   {  0.0, -1.0, -1.0 },
                                   {  0.0,  1.0, -1.0 },
                                   { -1.0,  1.0, -1.0 },
                                   { -1.0, -1.0,  1.0 },
                                   {  0.0, -1.0,  1.0 },
                                   {  0.0,  1.0,  1.0 },
                                   { -1.0,  1.0,  1.0 }};

        // space geometry interpolations for integration mesh
        Interpolation_Function_Base * tIntegSpaceInterpolation = tGeomInterpRule.create_space_interpolation_function();
        uint tIntegNumSpaceBases    = tIntegSpaceInterpolation->get_number_of_bases();
        uint tIntegNumSpaceDim      = tIntegSpaceInterpolation->get_number_of_dimensions();
        uint tIntegNumParamSpaceDim = tIntegSpaceInterpolation->get_number_of_param_dimensions();

        // create a integration rule
        Integration_Rule tIntegrationRule( tIntegGeometryType,
                                           Integration_Type::GAUSS,
                                           Integration_Order::HEX_3x3x3,
                                           Integration_Type::GAUSS,
                                           Integration_Order::BAR_3 );

        // create a side integrator
        Integrator tIntegrator( tIntegrationRule );

        //get number of integration points, integration points and weights
        uint             tNumOfIntegPoints = tIntegrator.get_number_of_points();
        Matrix< DDRMat > tIntegPoints      = tIntegrator.get_points();
        Matrix< DDRMat > tIntegWeights     = tIntegrator.get_weights();

        // boolean for surface check
        bool tSurfaceCheck = true;

        // init the surface of the integration mesh
        real tSurface = 0;

        // loop over the integration points
        for ( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
        {
            // get the treated integration point location in integration space
            Matrix< DDRMat > tIntegPointI = tIntegPoints.get_column( iGP );

            // get the treated integration point location in the interpolation space
            //------------------------------------------------------------------------------
            // unpack the space and time param coords of the integration point
            Matrix< DDRMat > tXi = tIntegPointI( { 0, tIntegNumParamSpaceDim-1 }, { 0, 0 } );
            Matrix< DDRMat > tTau( 1, 1, tIntegPointI( tIntegNumParamSpaceDim ) );

            // evaluate space interpolation shape functions at integration point
            Matrix< DDRMat > tNIntegSpace = tIntegSpaceInterpolation->eval_N( tXi );

            // evaluate time interpolation shape functions at aParamPoint
            Matrix< DDRMat > tNTime = tTimeInterpolation->eval_N( tTau );

            // build space time interpolation functions for integration mesh
            Matrix< DDRMat > tN = reshape( trans( tNIntegSpace ) * tNTime, 1, tNumTimeBases*tIntegNumSpaceBases );

            // get the parametric coordinates of the integration mesh in interpolation mesh
            Matrix< DDRMat > tInterpParamCoords( tNumTimeDim+tNumParamSpaceDim, tNumTimeBases*tIntegNumSpaceBases );

            // get the time parametric coordinates
            Matrix< DDRMat > tTimeParamCoords  = tTimeInterpolation->get_param_coords();

            // get a vector of ones
            Matrix< DDRMat > tOnes( 1, tIntegNumSpaceBases, 1.0 );

            // loop on the time bases
            for( uint i = 0; i < tNumTimeBases; i++ )
            {
                // fill the space time parametric coordinates matrix with space coordinates
                tInterpParamCoords( { 0, tIntegNumParamSpaceDim-1 }, { i * tIntegNumSpaceBases, ( i + 1 ) * tIntegNumSpaceBases-1 })
                    = trans( tXiHat );

                // fill the space time parametric coordinates matrix with time coordinates
                tInterpParamCoords( { tIntegNumParamSpaceDim, tIntegNumParamSpaceDim }, { i * tIntegNumSpaceBases, ( i + 1 ) * tIntegNumSpaceBases-1 })
                    = tTimeParamCoords( i ) * tOnes;
             }

            // compute the parametric coordinates of the SideParamPoint in the parent reference element
            Matrix< DDRMat > tRefIntegPointI = tInterpParamCoords * trans( tN );

            // evaluate detJ
            //------------------------------------------------------------------------------
            real tDetJ1 = tGeomInterpolator.det_J( tRefIntegPointI );

            // get the space jacobian
            Matrix <DDRMat> tdNSpacedXi = tIntegSpaceInterpolation->eval_dNdXi( tXi );
            Matrix< DDRMat > tSpaceJt   = tdNSpacedXi * tXiHat ;
            // get the time Jacobian
            Matrix< DDRMat > tTauHat     = tTimeInterpolation->get_param_coords();
            Matrix< DDRMat > tdNTimedTau = tTimeInterpolation->eval_dNdXi( tTau );
            Matrix< DDRMat > tTimeJt     = tdNTimedTau * trans( tTauHat ) ;
            real tDetJ2 = det( tSpaceJt ) * det( tTimeJt );

            real tDetJ = tDetJ1 * tDetJ2;

            // add contribution to the surface
            tSurface = tSurface + tDetJ * tIntegWeights( iGP );
        }

        // check the surface value
        tSurfaceCheck = tSurfaceCheck && ( std::abs( tSurface - 0.5 ) < tEpsilon );
        std::cout<<tSurface<<std::endl;
        REQUIRE( tSurfaceCheck );

    }


    SECTION( "Interpolation mesh HEX8 - Integration mesh TET4 " )
       {
           // define an interpolation mesh
           //------------------------------------------------------------------------------
           // define a HEX8 space element, i.e. space coordinates xHat
                   Matrix< DDRMat > tXHat = {{ 0.0, 0.0, 0.0 },
                                             { 1.0, 0.0, 0.0 },
                                             { 1.0, 1.0, 0.0 },
                                             { 0.0, 1.0, 0.0 },
                                             { 0.0, 0.0, 1.0 },
                                             { 1.0, 0.0, 1.0 },
                                             { 1.0, 1.0, 1.0 },
                                             { 0.0, 1.0, 1.0 }};

           // define a line time element, i.e. time coordinates tHat
           Matrix< DDRMat > tTHat = {{ 0.0 },
                                     { 1.0 },
                                     { 0.5 }};

           // create a space and time geometry interpolation rule
           Interpolation_Rule tGeomInterpRule( mtk::Geometry_Type::HEX,
                                               Interpolation_Type::LAGRANGE,
                                               mtk::Interpolation_Order::LINEAR,
                                               Interpolation_Type::LAGRANGE,
                                               mtk::Interpolation_Order::QUADRATIC );

           // create a space and time geometry interpolator
           Geometry_Interpolator tGeomInterpolator( tGeomInterpRule, true );

           //set the coefficients xHat, tHat
           tGeomInterpolator.set_space_coeff( tXHat );
           tGeomInterpolator.set_time_coeff(  tTHat );

           // space and time geometry interpolations
           Interpolation_Function_Base * tSpaceInterpolation = tGeomInterpRule.create_space_interpolation_function();
           Interpolation_Function_Base * tTimeInterpolation  = tGeomInterpRule.create_time_interpolation_function();

           uint tNumSpaceBases = tSpaceInterpolation->get_number_of_bases();
           uint tNumTimeBases  = tTimeInterpolation->get_number_of_bases();
           uint tNumSpaceDim   = tSpaceInterpolation->get_number_of_dimensions();
           uint tNumTimeDim    = tTimeInterpolation->get_number_of_dimensions();
           uint tNumParamSpaceDim = tSpaceInterpolation->get_number_of_param_dimensions();

           // define an integration mesh
           //------------------------------------------------------------------------------
           // integration mesh geometry type
           mtk::Geometry_Type tIntegGeometryType = mtk::Geometry_Type::TET;

           // integration mesh interpolation rule
           Interpolation_Rule tIntegInterpRule( tIntegGeometryType,
                                                Interpolation_Type::LAGRANGE,
                                                mtk::Interpolation_Order::LINEAR,
                                                Interpolation_Type::LAGRANGE,
                                                mtk::Interpolation_Order::QUADRATIC );

           // define a TET4 integration element, i.e. space param coordinates xiHat
           Matrix< DDRMat > tXiHat = {{ -1.0, -1.0, -1.0 },
                                      {  0.0, -1.0, -1.0 },
                                      { -1.0,  1.0, -1.0 },
                                      { -1.0, -1.0,  1.0 }};

           // space geometry interpolations for integration mesh
           Interpolation_Function_Base * tIntegSpaceInterpolation = tIntegInterpRule.create_space_interpolation_function();
           uint tIntegNumSpaceBases    = tIntegSpaceInterpolation->get_number_of_bases();
           uint tIntegNumSpaceDim      = tIntegSpaceInterpolation->get_number_of_dimensions();
           uint tIntegNumParamSpaceDim = tIntegSpaceInterpolation->get_number_of_param_dimensions();

           // create a integration rule
           Integration_Rule tIntegrationRule( tIntegGeometryType,
                                              Integration_Type::GAUSS,
                                              Integration_Order::TET_5,
                                              Integration_Type::GAUSS,
                                              Integration_Order::BAR_3 );

           // create a side integrator
           Integrator tIntegrator( tIntegrationRule );

           //get number of integration points, integration points and weights
           uint             tNumOfIntegPoints = tIntegrator.get_number_of_points();
           Matrix< DDRMat > tIntegPoints      = tIntegrator.get_points();
           Matrix< DDRMat > tIntegWeights     = tIntegrator.get_weights();

           // boolean for surface check
           bool tSurfaceCheck = true;

           // init the surface of the integration mesh
           real tSurface = 0;

           // loop over the integration points
           for ( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
           {
               // get the treated integration point location in integration space
               Matrix< DDRMat > tIntegPointI = tIntegPoints.get_column( iGP );

               // get the treated integration point location in the interpolation space
               //------------------------------------------------------------------------------
               // unpack the space and time param coords of the integration point
               Matrix< DDRMat > tXi = tIntegPointI( { 0, tIntegNumParamSpaceDim-1 }, { 0, 0 } );
               Matrix< DDRMat > tTau( 1, 1, tIntegPointI( tIntegNumParamSpaceDim ) );
               //print(tXi,"tXi");
               //print(tTau,"tTau");

               // evaluate space interpolation shape functions at integration point
               Matrix< DDRMat > tNIntegSpace = tIntegSpaceInterpolation->eval_N( tXi );

               // evaluate time interpolation shape functions at aParamPoint
               Matrix< DDRMat > tNTime = tTimeInterpolation->eval_N( tTau );

               // build space time interpolation functions for integration mesh
               Matrix< DDRMat > tN = reshape( trans( tNIntegSpace ) * tNTime, 1, tNumTimeBases*tIntegNumSpaceBases );

               // get the parametric coordinates of the integration mesh in interpolation mesh
               Matrix< DDRMat > tInterpParamCoords( tNumTimeDim+tNumParamSpaceDim, tNumTimeBases*tIntegNumSpaceBases );

               // get the time parametric coordinates
               Matrix< DDRMat > tTimeParamCoords  = tTimeInterpolation->get_param_coords();

               // get a vector of ones
               Matrix< DDRMat > tOnes( 1, tIntegNumSpaceBases, 1.0 );

               // loop on the time bases
               for( uint i = 0; i < tNumTimeBases; i++ )
               {
                   // fill the space time parametric coordinates matrix with space coordinates
                   tInterpParamCoords( { 0, tNumParamSpaceDim-1 }, { i * tIntegNumSpaceBases, ( i + 1 ) * tIntegNumSpaceBases-1 })
                       = trans( tXiHat );

                   // fill the space time parametric coordinates matrix with time coordinates
                   tInterpParamCoords( { tNumParamSpaceDim, tNumParamSpaceDim }, { i * tIntegNumSpaceBases, ( i + 1 ) * tIntegNumSpaceBases-1 })
                       = tTimeParamCoords( i ) * tOnes;
                }
                //print(tInterpParamCoords,"tInterpParamCoords");

                // compute the parametric coordinates of the SideParamPoint in the parent reference element
                Matrix< DDRMat > tRefIntegPointI = tInterpParamCoords * trans( tN );

                // evaluate detJ
                //------------------------------------------------------------------------------
                real tDetJ1 = tGeomInterpolator.det_J( tRefIntegPointI );

                // get the space jacobian
                Matrix <DDRMat> tdNSpacedXi = tIntegSpaceInterpolation->eval_dNdXi( tXi );
                Matrix< DDRMat > tSpaceJt   = tdNSpacedXi * tXiHat;
                // get the time Jacobian
                Matrix< DDRMat > tTauHat     = tTimeInterpolation->get_param_coords();
                Matrix< DDRMat > tdNTimedTau = tTimeInterpolation->eval_dNdXi( tTau );
                Matrix< DDRMat > tTimeJt     = tdNTimedTau * trans( tTauHat );

                Matrix< DDRMat > tSpaceJt2( tIntegNumParamSpaceDim, tIntegNumParamSpaceDim, 1.0 );
                tSpaceJt2({ 1, tIntegNumParamSpaceDim-1 },{ 0, tIntegNumParamSpaceDim-1 }) = trans( tSpaceJt );
                real detJSpace = det( tSpaceJt2 ) / 6.0;
                real tDetJ2 = detJSpace * det( tTimeJt );

                real tDetJ = tDetJ1 * tDetJ2;

                // add contribution to the surface
                tSurface = tSurface + tDetJ * tIntegWeights( iGP );
            }

            // check the surface value
            tSurfaceCheck = tSurfaceCheck && ( std::abs( tSurface - 0.083333333334 ) < tEpsilon );
            std::cout<<tSurface<<std::endl;
            REQUIRE( tSurfaceCheck );

       }
}

