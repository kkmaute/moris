// Third-party header files.
#include <catch.hpp>
#include <utility>
#include <iostream>

// MORIS project header files.
#include "algorithms.hpp"
#include "linalg.hpp"
#include "chronos.hpp"
#include "ios.hpp"
// ----------------------------------------------------------------------------

TEST_CASE(
        "moris::Tensor",
        "[linalgebra],[Tensor]" )
        {

    moris::Tensor< moris::real, 1, 3, true > T13t; // this creates a 1st order, 3D tensor.
    T13t(0) = 1.; T13t(1) = 2.; T13t(2) = 4.;

    moris::Tensor< moris::real, 1, 2, true > T12t;
    T12t(0) = 1.; T12t(1) = 2.;

    // symmetric 2nd order, 3D tensor called T23t.
    #include "linalg/cl_Tensor/cl_Tensor_23sym.inc"

    moris::Tensor< moris::real, 2, 2, false > T22f; // unsymmetric 2nd order 2D tensor
    T22f(0,0) = 1.; T22f(0,1) = 2.;
    T22f(1,0) = 3.; T22f(1,1) = 4.;

    moris::Tensor< moris::real ,2, 3, false > T23f; // unsymmetric 2nd order 3D tensor
    T23f(0,0) = 1.; T23f(0,1) = 2.; T23f(0,2) = 3.;
    T23f(1,0) = 2.; T23f(1,1) = 4.; T23f(1,2) = 5.;
    T23f(2,0) = 3.; T23f(2,1) = 5.; T23f(2,2) = 6.;

    moris::Tensor< moris::real, 3, 3, true > T33t; // symmetric (default) 3rd order, 3D tensor
    T33t(0,0,0) = 1.; T33t(1,1,0) = 1.; T33t(2,2,0) = 1.; T33t(1,2,0) = 1.; T33t(0,2,0) = 1.; T33t(0,1,0) = 1.;
    T33t(0,0,1) = 1.; T33t(1,1,1) = 1.; T33t(2,2,1) = 1.; T33t(1,2,1) = 1.; T33t(0,2,1) = 1.; T33t(0,1,1) = 1.;
    T33t(0,0,2) = 1.; T33t(1,1,2) = 1.; T33t(2,2,2) = 1.; T33t(1,2,2) = 1.; T33t(0,2,2) = 1.; T33t(0,1,2) = 1.;

    moris::Tensor< moris::real, 4, 3, true > T43t; // symmetric 4th order 3D tensor
    T43t(0,0,0,0) = 2.; T43t(1,1,0,0) = 1.; T43t(2,2,0,0) = 1.; T43t(1,2,0,0) = 1.; T43t(0,2,0,0) = 1.; T43t(0,1,0,0) = 1.;
    T43t(0,0,1,1) = 1.; T43t(1,1,1,1) = 3.; T43t(2,2,1,1) = 1.; T43t(1,2,1,1) = 1.; T43t(0,2,1,1) = 1.; T43t(0,1,1,1) = 1.;
    T43t(0,0,2,2) = 1.; T43t(1,1,2,2) = 1.; T43t(2,2,2,2) = 4.; T43t(1,2,2,2) = 1.; T43t(0,2,2,2) = 1.; T43t(0,1,2,2) = 1.;
    T43t(0,0,1,2) = 1.; T43t(1,1,1,2) = 1.; T43t(2,2,1,2) = 1.; T43t(1,2,1,2) = 1.; T43t(0,2,1,2) = 1.; T43t(0,1,1,2) = 7.;
    T43t(0,0,0,2) = 1.; T43t(1,1,0,2) = 1.; T43t(2,2,0,2) = 1.; T43t(1,2,0,2) = 1.; T43t(0,2,0,2) = 1.; T43t(0,1,0,2) = 1.;
    T43t(0,0,0,1) = 5.; T43t(1,1,0,1) = 1.; T43t(2,2,0,1) = 1.; T43t(1,2,0,1) = 6.; T43t(0,2,0,1) = 1.; T43t(0,1,0,1) = 1.;

    moris::Tensor< moris::real, 4, 2, false > T42f; // unsymmetric 4th order 2D tensor
    T42f(0,0,0,0) = 1.; T42f(1,1,0,0) = 1.; T42f(0,1,0,0) = 1.; T42f(1,0,0,0) = 1.;
    T42f(0,0,1,1) = 1.; T42f(1,1,1,1) = 1.; T42f(0,1,1,1) = 1.; T42f(1,0,1,1) = 1.;
    T42f(0,0,0,1) = 1.; T42f(1,1,0,1) = 1.; T42f(0,1,0,1) = 1.; T42f(1,0,0,1) = 1.;
    T42f(0,0,1,0) = 1.; T42f(1,1,1,0) = 1.; T42f(0,1,1,0) = 1.; T42f(1,0,1,0) = 1.;

    moris::Tensor< moris::real, 4, 2, true > T42t; // symmetric 4th order 2D tensor
    T42t(0,0,0,0) = 1.; T42t(1,1,0,0) = 2.; T42t(0,1,0,0) = 3.;
    T42t(0,0,1,1) = 4.; T42t(1,1,1,1) = 5.; T42t(0,1,1,1) = 6.;
    T42t(0,0,0,1) = 7.; T42t(1,1,0,1) = 8.; T42t(0,1,0,1) = 9.;

    /*
    SECTION( "moris::Tensor failing accessors")
    {
        // Because of how these functions are defined, these "tests"
        // fail at compile time, i.e. this section does not compile.
        // These are thus, commented out and can be periodically
        // un-commented to ensure they still fail to compile.

        // 1st order tensors should not have 2+ index accessors.
        REQUIRE_THROWS( T13t(0,0) ); REQUIRE_THROWS( T13t(0,0,0) ); REQUIRE_THROWS( T13t(0,0,0,0) );

        // 2nd order tensors should not have  1 nor 3+ index accessors.
        REQUIRE_THROWS( T22f(0) );   REQUIRE_THROWS( T22f(0,0,0) ); REQUIRE_THROWS( T22f(0,0,0,0) );

        // 3rd order tensors should not have 1, 2 nor 4 index accessors.
        REQUIRE_THROWS( T33t(0) );   REQUIRE_THROWS( T33t(0,0) );   REQUIRE_THROWS( T33t(0,0,0,0) );

        // 4th order tensors should not have 3- index accessors.
        REQUIRE_THROWS( T42t(0) );   REQUIRE_THROWS( T42t(0,0) );   REQUIRE_THROWS( T42t(0,0,0) );
    }
    */

    SECTION( " moris::Tensor copy")
    {
        // copy of 2nd order tensor
        // checks that the order and dim were copied correctly and that
        // the components were stored correctly.

        moris::Tensor< moris::real, 2, 3, true > B(T23t); // this copies the tensor, checking the order
        REQUIRE( B.order() == 2 ); REQUIRE( B.dim() == 3 );
        REQUIRE( B(0,0) == 1. );   REQUIRE(B(2,1) == B(1,2) ); REQUIRE( B(2,2) == 6 );

        moris::Tensor< moris::real, 2, 2, false > uB = T22f;
        REQUIRE( uB.order() == 2 ); REQUIRE( uB.dim() == 2 );
        REQUIRE( moris::equal_to( uB(0,0), 1.0) ); REQUIRE( moris::equal_to( uB(0,1), 2.0) );
        REQUIRE( moris::equal_to( uB(1,0), 3.0) ); REQUIRE( moris::equal_to( uB(1,1), 4.0) );

        moris::Tensor< moris::real, 4, 2, true> D(T42t);
        REQUIRE( D.order() == 4 ); REQUIRE( D.dim() == 2 );
        REQUIRE( moris::equal_to( D(0,0,0,0), 1.0 ) ); REQUIRE( moris::equal_to( D(1,1,0,0), 2.0 ) );
        REQUIRE( moris::equal_to( D(1,1,0,1), 8.0 ) ); REQUIRE( moris::equal_to( D(0,1,0,1), 9.0 ) );

    } // ends section

    SECTION( "moris::Tensor", "inner products")
    {
#include "linalg/cl_Tensor/cl_Tensor_prod22.inc"
        REQUIRE( moris::equal_to( prod22, 129. ) );

        moris::real prod22us = T23f*T23f;
        REQUIRE( moris::equal_to( prod22us, 129. ) );

        moris::real prod11 = T13t*T13t;
        REQUIRE( moris::equal_to( prod11, 21.) );

    }

    SECTION( "moris::Tensor", "tensor 4-2 product")
    {
#include "linalg/cl_Tensor/cl_Tensor_prod42.inc"
        REQUIRE( moris::equal_to( prod42(0,0), 48. ) ); REQUIRE( moris::equal_to( prod42(1,1), 39. ) );
        REQUIRE( moris::equal_to( prod42(2,2), 49. ) ); REQUIRE( moris::equal_to( prod42(1,2), 51. ) );
        REQUIRE( moris::equal_to( prod42(0,2), 31. ) ); REQUIRE( moris::equal_to( prod42(1,0), 91. ) );

        moris::Tensor< moris::real, 2, 2, false > prod42us = T42f*T22f;
        REQUIRE( moris::equal_to( prod42us(0,0), 10. ) ); REQUIRE( moris::equal_to( prod42us(0,1), 10. ) );
        REQUIRE( moris::equal_to( prod42us(1,0), 10. ) ); REQUIRE( moris::equal_to( prod42us(1,1), 10. ) );
    }

    SECTION( "moris::Tensor", "tensor 3-1 product")
    {
        moris::Tensor< moris::real, 2, 3, true > prod31(T33t*T13t);
        REQUIRE( moris::equal_to( prod31(0,0), 7. ) ); REQUIRE( moris::equal_to( prod31(0,1), 7. ) );
        REQUIRE( moris::equal_to( prod31(1,1), 7. ) ); REQUIRE( moris::equal_to( prod31(0,2), 7. ) );
        REQUIRE( moris::equal_to( prod31(2,2), 7. ) ); REQUIRE( moris::equal_to( prod31(2,1), 7. ) );
    }

    SECTION( "moris::Tensor", "tensor 2-1 product")
    {
        moris::Tensor< moris::real, 1, 3, true > prod21 = T23t*T13t;
        REQUIRE( prod21.order() == 1 ); REQUIRE( prod21.dim() == 3 );
        REQUIRE( moris::equal_to( prod21(0), 17. ) );
        REQUIRE( moris::equal_to( prod21(1), 30. ) );
        REQUIRE( moris::equal_to( prod21(2), 37. ) );

        moris::Tensor< moris::real, 1, 3, true > prod21u = T23f*T13t;
        REQUIRE( prod21u.order() == 1 ); REQUIRE( prod21u.dim() == 3 );
        REQUIRE( moris::equal_to( prod21u(0), 17. ) );
        REQUIRE( moris::equal_to( prod21u(1), 30. ) );
        REQUIRE( moris::equal_to( prod21u(2), 37. ) );

        moris::Tensor< moris::real, 1, 3 , true> prod21f = prod(T23f, T13t);
        REQUIRE( prod21f.order() == 1 ); REQUIRE( prod21f.dim() == 3 );
        REQUIRE( moris::equal_to( prod21f(0), 17. ) );
        REQUIRE( moris::equal_to( prod21f(1), 30. ) );
        REQUIRE( moris::equal_to( prod21f(2), 37. ) );
    }

    SECTION( "moris::Tensor", "3D Speed Test")
    {
        // Performs the product of a 3-by-3 matrix and a 3-by-1 matrix
        // (i.e.a 2nd order tensor and a 1st order tensor). The Tensor
        // product should not take more time.

        moris::uint numOfOps = 1000;


        moris::Tensor< moris::real, 1, 3, true > prod21;

        // using moris::Tensor
        moris::tic timer2;

        for(moris::uint i = 0; i < numOfOps; i++)
        {
            moris::Mat<moris::real > tProd = T23t*T13t;

//            prod21 = T23t*T13t;
        }

        moris::real time2 = timer2.toc<moris::chronos::nanoseconds>().wall;


        // using moris::Mat
        moris::Mat< moris::real > tMat2 = { { 1.,2.,3.},{2.,4.,5.},{3.,5.,6.} };
        moris::Mat< moris::real > tMat1 = { { 1.}, {2.},{4.} };


        moris::tic timer1;

        for(moris::uint i = 0; i < numOfOps; i++)
        {
            moris::Mat<moris::real > tProd = tMat2*tMat1;
        }

        moris::real time1 = timer1.toc<moris::chronos::nanoseconds>().wall;


        // using prod() function, which wraps operator* for Tensors
        moris::tic timer3;

        for(moris::uint i = 0; i < numOfOps; i++)
        {
            moris::Mat<moris::real > tProd = prod(T23t,T13t);
        }

        moris::real time3 = timer3.toc<moris::chronos::nanoseconds>().wall;
#if defined(NDEBUG) || !defined(DEBUG)
        // timing test only makes sense with debug flags off
        // Checks that time with Tensor is less than 5% more than time with Mat
        CHECK( (time2 - time1) < 0.1*time1 );

        // operator wrapping should not cost more than a 5% overhead
        CHECK( (time3 - time2) < 0.10*time2 );
#endif
    }


    SECTION( "moris::Tensor", "2D Speed Test")
    {
        // Performs the product of a 2-by-2 matrix and a 2-by-1 matrix
        // (i.e.a 2nd order tensor and a 1st order tensor). The Tensor
        // product should not take more time.

        moris::uint numOfOps = 10000;

        // using moris::Mat
        moris::Mat< moris::real > tMat2 = { { 1.,2.},{3.,4.} };
        moris::Mat< moris::real > tMat1 = { { 1.}, {2.} };

        moris::tic timer1;

        for(moris::uint i = 0; i < numOfOps; i++)
        {
            moris::Mat<moris::real > tProd = tMat2*tMat1;
        }

        moris::real time1 = timer1.toc<moris::chronos::nanoseconds>().wall;

        //using moris::Tensor
        moris::tic timer2;

        for(moris::uint i = 0; i < numOfOps; i++)
        {
            moris::Tensor< moris::real, 1, 2, true > prod21 = T22f*T12t;
        }

        moris::real time2 = timer2.toc<moris::chronos::nanoseconds>().wall;
#if defined(NDEBUG) || !defined(DEBUG)
        // timing test only makes sense with debug flags off
        // Checks that time with Tensor is less than 5% more than time with Mat
        CHECK( (time2 - time1) < 0.10*time1 );
#endif
    }
} // ends Test case

