#include <string>
#include <catch.hpp>

#include "typedefs.hpp" //MRS/COR/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "assert.hpp"

#include "cl_MTK_Enums.hpp" //MTK/src

#include "cl_FEM_Enums.hpp"                       //FEM/INT/src
#include "cl_FEM_Field_Interpolator.hpp"          //FEM//INT//src
#include "cl_FEM_IWG_Olsson_CLS_Bulk.hpp"         //FEM//INT//src
#include "cl_FEM_IWG_Olsson_CLS_Interface.hpp"    //FEM//INT//src

#include "op_equal_equal.hpp"

using namespace moris;
using namespace fem;

TEST_CASE( "IWG_Olsson_CLS", "[moris],[fem],[IWG_OCLS]" )
{
    // space and time geometry interpolators
    //------------------------------------------------------------------------------
    //create a space geometry interpolation rule
    Interpolation_Rule tGeomInterpRule( mtk::Geometry_Type::QUAD,
                                        Interpolation_Type::LAGRANGE,
                                        mtk::Interpolation_Order::LINEAR,
                                        Interpolation_Type::LAGRANGE,
                                        mtk::Interpolation_Order::LINEAR);

    //create a space and a time geometry interpolator
    Geometry_Interpolator tGeomInterpolator( tGeomInterpRule );

    //create space coeff xHat
    Matrix< DDRMat > tXHat( 4, 2 );
    tXHat( 0, 0 ) = 0.0; tXHat( 0, 1 ) = 0.0;
    tXHat( 1, 0 ) = 3.0; tXHat( 1, 1 ) = 1.25;
    tXHat( 2, 0 ) = 4.5; tXHat( 2, 1 ) = 4.0;
    tXHat( 3, 0 ) = 1.0; tXHat( 3, 1 ) = 3.25;

    //create time coeff tHat
    Matrix< DDRMat > tTHat( 2, 1 );
    tTHat( 0 ) = 0.0;
    tTHat( 1 ) = 5.0;

    //set the coefficients xHat, tHat
    tGeomInterpolator.set_coeff( tXHat, tTHat);

    // field interpolator
    //------------------------------------------------------------------------------
    //create a space time interpolation rule
    Interpolation_Rule tInterpolationRule ( mtk::Geometry_Type::QUAD,
                                            Interpolation_Type::LAGRANGE,
                                            mtk::Interpolation_Order::LINEAR,
                                            Interpolation_Type::LAGRANGE,
                                            mtk::Interpolation_Order::LINEAR );

    //set the number of interpolated fields
    uint tNumberOfFields = 1;

    //create a field interpolator
    Field_Interpolator * tFieldInterpolator = new Field_Interpolator( tNumberOfFields,
                                                                      tInterpolationRule,
                                                                      tGeomInterpolator );

    //get the number of basis for space time
    uint tNSpaceTimeBases = tFieldInterpolator->get_number_of_space_time_bases();

    //create field coeff uHat
    Matrix< DDRMat > tUHat( tNSpaceTimeBases , tNumberOfFields );
    tUHat( 0 ) = -1.0;
    tUHat( 1 ) =  7.0;
    tUHat( 2 ) = -6.0;
    tUHat( 3 ) = 82.0;
    tUHat( 4 ) = -1.0;
    tUHat( 5 ) =  7.0;
    tUHat( 6 ) = -6.0;
    tUHat( 7 ) = 82.0;

    //set the coefficients uHat
    tFieldInterpolator->set_coeff( tUHat );

    // create evaluation point xi, tau
    Matrix< DDRMat > tXi( 2, 1 );
    tXi( 0, 0 ) =  0.35; tXi( 1, 0 ) = -0.25;
    Matrix< DDRMat > tTau( 1, 1 );
    tTau( 0, 0 ) = 0.0;

    //set the evaluation point xi, tau
    tFieldInterpolator->set_space_time( tXi, tTau );

    // define an epsilon environment
    double tEpsilon = 1E-6;

    SECTION( "IWG_Olsson_CLS_Bulk : check residual and jacobian" )
    {
        // IWG
        //------------------------------------------------------------------------------
        // set field upper and lower bounds
        double tFieldUpperbound =  1.0;
        double tFieldLowerbound = -1.0;

        // set a method parameter epsilon
        double tEpsilonParameter = 1.0;

        // create an IWG Olsson CLS Bulk
        IWG_Olsson_CLS_Bulk tIWGOCLSBulk( tFieldInterpolator,
                                          tFieldUpperbound,
                                          tFieldLowerbound,
                                          tEpsilonParameter );

        // check evaluation of the residual for IWG Olsson CLS Bulk ?
        //------------------------------------------------------------------------------
        // create a velocity field at evaluation point
        Matrix< DDRMat > tFieldNormal( 2, 1, 1.0);

        // evaluate the residual from IWG_Olsson_CLS_Bulk
        Matrix< DDRMat > tResidualOCLSBulk;
        tIWGOCLSBulk.compute_residual( tResidualOCLSBulk,
                                       tFieldNormal );

        // check evaluation of the jacobian for IWG Olsson CLS Bulk by FD
        //------------------------------------------------------------------------------
        // evaluate the jacobian from IWG_Olsson_CLS_Bulk
        Matrix< DDRMat > tJacobianOCLSBulk;
        tIWGOCLSBulk.compute_jacobian( tJacobianOCLSBulk,
                                       tFieldNormal );

        //define a boolean for check
        bool tJacobianBulk = true;

        //evaluate the jacobian by FD
        Matrix< DDRMat > tUHatPert, tResidualOCLSBulkPert, tJacobianRow;
        real tPert;
        for( uint k=0; k<tNSpaceTimeBases; k++ )
        {
            //set the perturbed values of uHat
            tUHatPert = tUHat;
            tPert = 1e-6 * tUHatPert( k );
            tUHatPert( k ) = tUHatPert( k ) + tPert;

            //set the coefficients uHatPert
            tFieldInterpolator->set_coeff( tUHatPert );

            // compute the perturbed residual
            tIWGOCLSBulk.compute_residual( tResidualOCLSBulkPert,
                                           tFieldNormal );

            // compute the jacobian by FD for the kth uHat
            tJacobianRow = ( tResidualOCLSBulkPert - tResidualOCLSBulk ) / tPert;

            // check the value of the jacobian evaluated from the IWG with FD
            for( uint i=0; i<tNSpaceTimeBases; i++ )
            {
                tJacobianBulk = tJacobianBulk && ( tJacobianRow( i ) - tJacobianOCLSBulk( i, k ) < tEpsilon );
            }
        }
        REQUIRE( tJacobianBulk );

        //set the coefficients uHat
        tFieldInterpolator->set_coeff( tUHat );
    }

    SECTION( "IWG_Olsson_CLS_Interface : check residual and jacobian" )
    {
        // check evaluation of the residual for IWG Olsson_CLS_Interface ?
        //------------------------------------------------------------------------------
        // set field upper and lower bounds
        double tFieldUpperbound =  1.0;
        double tFieldLowerbound = -1.0;

        // set a method parameter epsilon
        double tEpsilonParameter = 1.0;

        // create an IWG HamiltonJacobi Bulk
        IWG_Olsson_CLS_Interface tIWGOCLSInterface( tFieldInterpolator,
                                                    tFieldUpperbound,
                                                    tFieldLowerbound,
                                                    tEpsilonParameter );

         // create a velocity field at evaluation point
         Matrix< DDRMat > tFieldNormal( 2, 1, 1.0);

         //create an interface normal at evaluation point
         Matrix< DDRMat > tInterfaceNormal( 2, 1, 1.0);

         // evaluate the residual from IWG_Olsson_CLS_Interface
         Matrix< DDRMat > tResidualOCLSInterface;
         tIWGOCLSInterface.compute_residual( tResidualOCLSInterface,
                                             tFieldNormal,
                                             tInterfaceNormal );

         // check evaluation of the jacobian for IWG_Olsson_CLS_Interface by FD
         //------------------------------------------------------------------------------
         // evaluate the jacobian from IWG_Olsson_CLS_Interface
         Matrix< DDRMat > tJacobianOCLSInterface;
         tIWGOCLSInterface.compute_jacobian( tJacobianOCLSInterface,
                                             tFieldNormal,
                                             tInterfaceNormal );

        //define a boolean for check
        bool tCheckJacobianInterface = true;

        //evaluate the jacobian by FD
        Matrix< DDRMat > tUHatPert, tResidualOCLSInterfacePert, tJacobianRow;
        real tPert;
        for( uint k=0; k<tNSpaceTimeBases; k++ )
        {
            //set the perturbed values of uHat
            tUHatPert = tUHat;
            tPert = 1e-6 * tUHatPert( k );
            tUHatPert( k ) = tUHatPert( k ) + tPert;

            //set the coefficients uHatPert
            tFieldInterpolator->set_coeff( tUHatPert );

            // compute the perturbed residual
            tIWGOCLSInterface.compute_residual( tResidualOCLSInterfacePert,
                                                tFieldNormal,
                                                tInterfaceNormal );

            // compute the jacobian by FD for the kth uHat
            tJacobianRow = ( tResidualOCLSInterfacePert - tResidualOCLSInterface ) / tPert;

            // check the value of the jacobian evaluated from the IWG with FD
            for( uint i=0; i<tNSpaceTimeBases; i++ )
            {
                tCheckJacobianInterface = tCheckJacobianInterface && ( tJacobianRow( i ) - tJacobianOCLSInterface( i, k ) < tEpsilon );
            }
        }
        REQUIRE( tCheckJacobianInterface );

        tFieldInterpolator->set_coeff( tUHat );
    }


    SECTION( "IWG_Olsson_CLS_Bulk : check simultaneous evaluation of residual and jacobian" )
    {

        // define an epsilon environment
        double tEpsilon2 = 1E-12;

        // set field upper and lower bounds
        double tFieldUpperbound =  1.0;
        double tFieldLowerbound = -1.0;

        // set a method parameter epsilon
        double tEpsilonParameter = 1.0;

        // create an IWG Olsson CLS Bulk
        IWG_Olsson_CLS_Bulk tIWGOCLSBulk( tFieldInterpolator,
                                          tFieldUpperbound,
                                          tFieldLowerbound,
                                          tEpsilonParameter );

        // check evaluation of the residual for IWG Olsson CLS Bulk ?
        //------------------------------------------------------------------------------
        // create a velocity field at evaluation point
        Matrix< DDRMat > tFieldNormal( 2, 1, 1.0);

        // evaluate the residual from IWG_Olsson_CLS_Bulk
        Matrix< DDRMat > tResidualOCLSBulk;
        tIWGOCLSBulk.compute_residual( tResidualOCLSBulk,
                                       tFieldNormal );

        // check evaluation of the jacobian for IWG Olsson CLS Bulk by FD
        //------------------------------------------------------------------------------
        // evaluate the jacobian from IWG_Olsson_CLS_Bulk
        Matrix< DDRMat > tJacobianOCLSBulk;
        tIWGOCLSBulk.compute_jacobian( tJacobianOCLSBulk,
                                       tFieldNormal );

        // evaluate the residual and jacobian from IWG_Olsson_CLS_Bulk
        Matrix< DDRMat > tResidualOCLSBulk2, tJacobianOCLSBulk2;
        tIWGOCLSBulk.compute_jacobian_and_residual( tJacobianOCLSBulk2,
                                                    tResidualOCLSBulk2,
                                                    tFieldNormal );

        // check simulatneous evaluation of the jacobian and residual
        // for IWG Olsson CLS Bulk

        //define a boolean for check
        bool tCheckResidualBulkSim = true;
        bool tCheckJacobianBulkSim = true;

        for( uint k=0; k<tNSpaceTimeBases; k++ )
        {
            // check the value of the residual
            tCheckResidualBulkSim = tCheckResidualBulkSim && ( tResidualOCLSBulk2( k ) - tResidualOCLSBulk( k ) < tEpsilon2 );

            // check the value of the jacobian
           for( uint i=0; i<tNSpaceTimeBases; i++ )
           {
               tCheckJacobianBulkSim = tCheckJacobianBulkSim && ( tJacobianOCLSBulk2( i, k ) - tJacobianOCLSBulk( i, k ) < tEpsilon2 );
           }
        }
        REQUIRE( tCheckResidualBulkSim );
        REQUIRE( tCheckJacobianBulkSim );
    }

    SECTION( "IWG_Olsson_CLS_Interface : check simultaneous evaluation of residual and jacobian" )
    {
        // define an epsilon environment
        double tEpsilon2 = 1E-12;

        // set field upper and lower bounds
        double tFieldUpperbound =  1.0;
        double tFieldLowerbound = -1.0;

        // set a method parameter epsilon
        double tEpsilonParameter = 1.0;

        // create an IWG HamiltonJacobi Bulk
        IWG_Olsson_CLS_Interface tIWGOCLSInterface( tFieldInterpolator,
                                                    tFieldUpperbound,
                                                    tFieldLowerbound,
                                                    tEpsilonParameter );

         // create a velocity field at evaluation point
         Matrix< DDRMat > tFieldNormal( 2, 1, 1.0);

         //create an interface normal at evaluation point
         Matrix< DDRMat > tInterfaceNormal( 2, 1, 1.0);

         // evaluate the residual from IWG_Olsson_CLS_Interface
         Matrix< DDRMat > tResidualOCLSInterface;
         tIWGOCLSInterface.compute_residual( tResidualOCLSInterface,
                                             tFieldNormal,
                                             tInterfaceNormal );

         // evaluate the jacobian from IWG_Olsson_CLS_Interface
         Matrix< DDRMat > tJacobianOCLSInterface;
         tIWGOCLSInterface.compute_jacobian( tJacobianOCLSInterface,
                                             tFieldNormal,
                                             tInterfaceNormal );

        // evaluate the residual and jacobian from IWG_Olsson_CLS_Interface
        Matrix< DDRMat > tResidualOCLSInterface2, tJacobianOCLSInterface2;
        tIWGOCLSInterface.compute_jacobian_and_residual( tJacobianOCLSInterface2,
                                                         tResidualOCLSInterface2,
                                                         tFieldNormal,
                                                         tInterfaceNormal );

        // check simulatneous evaluation of the jacobian and residual
        // for IWG Olsson CLS Interface

        //define a boolean for check
        bool tCheckResidualInterfaceSim = true;
        bool tCheckJacobianInterfaceSim = true;

        for( uint k=0; k<tNSpaceTimeBases; k++ )
        {
            // check the value of the residual
            tCheckResidualInterfaceSim = tCheckResidualInterfaceSim && ( tResidualOCLSInterface2( k ) - tResidualOCLSInterface( k ) < tEpsilon2 );

            // check the value of the jacobian
           for( uint i=0; i<tNSpaceTimeBases; i++ )
           {
               tCheckJacobianInterfaceSim = tCheckJacobianInterfaceSim && ( tJacobianOCLSInterface2( i, k ) - tJacobianOCLSInterface( i, k ) < tEpsilon2 );
           }
        }
        REQUIRE( tCheckResidualInterfaceSim );
        REQUIRE( tCheckJacobianInterfaceSim );
    }
}
