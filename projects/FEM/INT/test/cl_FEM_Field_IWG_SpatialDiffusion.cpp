#include <string>
#include <catch.hpp>

#include "assert.hpp"

#include "cl_MTK_Enums.hpp" //MTK/src
#include "cl_FEM_Enums.hpp"                       //FEM//INT/src
#include "cl_FEM_Field_Interpolator.hpp"          //FEM//INT//src
#include "cl_FEM_IWG_Helmholtz_Bulk.hpp"          //FEM//INT//src
#include "cl_FEM_IWG_Isotropic_Spatial_Diffusion_Dirichlet.hpp"          //FEM//INT//src

#include "op_equal_equal.hpp"

using namespace moris;
using namespace fem;

TEST_CASE( "IWG_Helmholtz", "[moris],[fem],[IWG_SpatialDiff]" )
{

    // space and time geometry interpolators
    //------------------------------------------------------------------------------
    //create a space geometry interpolation rule
    Interpolation_Rule tGeomInterpRule( mtk::Geometry_Type::HEX,
                                        Interpolation_Type::LAGRANGE,
                                        mtk::Interpolation_Order::LINEAR,
                                        Interpolation_Type::LAGRANGE,
                                        mtk::Interpolation_Order::LINEAR );


    //create a space time geometry interpolator
    Geometry_Interpolator* tGeomInterpolator = new Geometry_Interpolator( tGeomInterpRule );

    //create space coeff xHat
    Matrix< DDRMat > tXHat = {{ 0.0, 0.0, 0.0 },
                              { 1.0, 0.0, 0.0 },
                              { 1.0, 1.0, 0.0 },
                              { 0.0, 1.0, 0.0 },
                              { 0.0, 0.0, 1.0 },
                              { 1.0, 0.0, 1.0 },
                              { 1.0, 1.0, 1.0 },
                              { 0.0, 1.0, 1.0 }};

    //create time coeff tHat
    Matrix< DDRMat > tTHat = {{ 0.0 }, { 1.0 }};

    //set the coefficients xHat, tHat
    tGeomInterpolator->set_coeff( tXHat, tTHat );

    // field interpolator
    //------------------------------------------------------------------------------
    //create a space time interpolation rule
    Interpolation_Rule tInterpRule ( mtk::Geometry_Type::HEX,
                                     Interpolation_Type::LAGRANGE,
                                     mtk::Interpolation_Order::LINEAR,
                                     Interpolation_Type::CONSTANT,
                                     mtk::Interpolation_Order::CONSTANT );

    //set the number of interpolated fields
    uint tNumOfFields = 1;

    //create a field interpolator
    Field_Interpolator* tTemp = new Field_Interpolator( tNumOfFields,
                                                        tInterpRule,
                                                        tGeomInterpolator );

    //get the number of basis for space time
    uint tNumOfBases = tTemp->get_number_of_space_time_bases();

    //create field coeff uHat
    Matrix< DDRMat > tTempHat( tNumOfBases , tNumOfFields );
    tTempHat( 0 ) = -1.0;
    tTempHat( 1 ) =  7.0;
    tTempHat( 2 ) = -6.0;
    tTempHat( 3 ) = 82.0;
    tTempHat( 4 ) = -1.0;
    tTempHat( 5 ) =  7.0;
    tTempHat( 6 ) = -6.0;
    tTempHat( 7 ) = 82.0;

    //set the coefficients uHat
    tTemp->set_coeff( tTempHat );

    // create evaluation point xi, tau
    Matrix< DDRMat > tXi = {{ 0.35}, {-0.25}, { 0.75}};
    Matrix< DDRMat > tTau( 1, 1, 0.0 );
    Matrix< DDRMat > tParamPoint = {{ 0.35}, {-0.25}, { 0.75}, { 0.0 }};

    //set the evaluation point xi, tau
    tTemp->set_space_time( tParamPoint );

    // define an epsilon environment
    double tEpsilon = 1E-6;

    SECTION( "IWG_Helmholtz_Bulk : check residual and jacobian" )
    {
        // IWG
        //------------------------------------------------------------------------------
        Cell< Field_Interpolator * > tFieldInterpolators( 1 , nullptr);
        tFieldInterpolators( 0 ) = tTemp;

        // create an IWG Helmholtz Bulk
        IWG_Isotropic_Spatial_Diffusion_Dirichlet tIWG;

        // set the nodal weak bcs
        Matrix< DDRMat > tNodalWeakBCs( tNumOfBases, 1, 1.0 );
        tIWG.set_nodal_weak_bcs( tNodalWeakBCs );

        // set the normal
        Matrix< DDRMat > tNormal( 3, 1, 1.0 );
        tIWG.set_normal( tNormal );

        // check evaluation of the residual for IWG Helmholtz Bulk ?
        //------------------------------------------------------------------------------
        // evaluate the residual
        Matrix< DDRMat > tResidual;
        tIWG.compute_residual( tResidual, tFieldInterpolators );

        // check evaluation of the jacobian  by FD
        //------------------------------------------------------------------------------
        // evaluate the jacobian
        Cell< Matrix< DDRMat > > tJacobians;
        tIWG.compute_jacobian( tJacobians,
                               tFieldInterpolators );

        //define a boolean for check
        bool tCheckJacobian = true;

        //evaluate the jacobian by FD
        Matrix< DDRMat > tTempHatPert, tResidualPert, tJacobianRow;
        real tPert;
        for( uint k = 0; k < tNumOfBases; k++ )
        {
            //set the perturbed values of uHat
        	tTempHatPert = tTempHat;
            tPert = 1e-4 * tTempHatPert( k );
            tTempHatPert( k ) = tTempHatPert( k ) + tPert;

            //set the coefficients uHatPert
            tFieldInterpolators( 0 )->set_coeff( tTempHatPert );

            // compute the perturbed residual
            tIWG.compute_residual( tResidualPert,
                                   tFieldInterpolators );

            // compute the jacobian by FD for the kth uHat
            tJacobianRow = ( tResidualPert - tResidual ) / tPert;

            // check the value of the jacobian evaluated from the IWG with FD
            for( uint i = 0; i < tNumOfBases; i++ )
            {
                tCheckJacobian = tCheckJacobian && ( std::abs( tJacobianRow( i ) - tJacobians( 0 )( i, k ) ) < tEpsilon );
            }
        }
        REQUIRE( tCheckJacobian );
    }

    delete tTemp;
    delete tGeomInterpolator;
}
