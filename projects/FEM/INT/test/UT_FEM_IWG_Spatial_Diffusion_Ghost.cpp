#include <string>
#include <catch.hpp>
#include "assert.hpp"

#include "cl_MTK_Enums.hpp" //MTK/src
#include "cl_FEM_Enums.hpp"                                 //FEM//INT/src
#include "cl_FEM_CM_Factory.hpp"                            //FEM//INT/src
#include "cl_FEM_IWG_Isotropic_Spatial_Diffusion_Ghost.hpp" //FEM//INT//src

#include "op_equal_equal.hpp"

moris::Matrix< moris::DDRMat > tFIValFunction_UTGhost( moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                                                       moris::Cell< moris::fem::Field_Interpolator* > & aDofFI,
                                                       moris::Cell< moris::fem::Field_Interpolator* > & aDvFI,
                                                       moris::fem::Geometry_Interpolator              * aGI )
{
    return aParameters( 0 ) * aDofFI( 0 )->val();
}

moris::Matrix< moris::DDRMat > tFIDerFunction_UTGhost( moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                                                       moris::Cell< moris::fem::Field_Interpolator* > & aDofFI,
                                                       moris::Cell< moris::fem::Field_Interpolator* > & aDvFI,
                                                       moris::fem::Geometry_Interpolator              * aGI )
{
    return aParameters( 0 ) * aDofFI( 0 )->N();
}

using namespace moris;
using namespace fem;

// This UT tests the isotropic spatial diffusion ghost IWG
// for QUAD, HEX geometry type
// for LINEAR, QUADRATIC and CUBIC interpolation order
TEST_CASE( "IWG_SpatialDiff_Ghost", "[moris],[fem],[IWG_SpatialDiff_Ghost]" )
{

    // define an epsilon environment
    real tEpsilon = 1E-4;

    // define a perturbation relative size
    real tPerturbation = 1E-4;

    for( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // set geometry inputs
        //------------------------------------------------------------------------------
    	// create geometry type
        mtk::Geometry_Type tGeometryType;

        // create space coeff xHat
        Matrix< DDRMat > tXHat;

        // create evaluation point xi, tau
        Matrix< DDRMat > tParamPoint;

        // create list with number of coeffs
        Matrix< DDRMat > tNumCoeffs;

        // create the normal
        Matrix< DDRMat > tNormal;

        switch( iSpaceDim )
        {
            case( 2 ):
            {
                // set geometry type
                tGeometryType = mtk::Geometry_Type::QUAD;

                // fill space coeff xHat
                tXHat = {{ 0.0, 0.0 },
                         { 1.0, 0.0 },
                        { 1.0, 1.0 },
                        { 0.0, 1.0 }};

               // fill evaluation point xi, tau
               tParamPoint = {{ 0.35}, {-0.25}, { 0.0 }};

               // number of coefficients
               tNumCoeffs = {{ 4 },{ 9 },{ 16 }};

               // set the normal
               tNormal = {{1.0},{0.0}};

               break;
            }
            case( 3 ):
            {
                // set geometry type
                tGeometryType = mtk::Geometry_Type::HEX;

                // fill space coeff xHat
                tXHat = {{ 0.0, 0.0, 0.0 },
                         { 1.0, 0.0, 0.0 },
                         { 1.0, 1.0, 0.0 },
                         { 0.0, 1.0, 0.0 },
                         { 0.0, 0.0, 1.0 },
                         { 1.0, 0.0, 1.0 },
                         { 1.0, 1.0, 1.0 },
                         { 0.0, 1.0, 1.0 }};

                // fill evaluation point xi, tau
                tParamPoint = {{ 0.35}, {-0.25}, { 0.75}, { 0.0 }};

                // number of coefficients
                tNumCoeffs = {{ 8 },{ 27 },{ 64 }};

                // set the normal
                tNormal = {{1.0},{0.0},{0.0}};

                break;
            }
            default:
            {
                MORIS_ERROR( false, " QUAD or HEX only." );
                break;
            }
        }

        // create a spatial diffusion virtual work ghost IWG
        //------------------------------------------------------------------------------
        // create an IWG Spatial Difffusion Ghost
        IWG_Isotropic_Spatial_Diffusion_Ghost tIWG;

        // set residual dof type
        tIWG.set_residual_dof_type( { MSI::Dof_Type::TEMP } );

        // set master dof type
        tIWG.set_dof_type_list( {{ MSI::Dof_Type::TEMP }});

        // set slave dof type
        tIWG.set_dof_type_list( {{ MSI::Dof_Type::TEMP }}, mtk::Master_Slave::SLAVE );

        // set master constitutive type
        tIWG.set_property_type_list( { fem::Property_Type::CONDUCTIVITY } );

        // set IWG normal
        tIWG.set_normal( tNormal );

        // space and time geometry interpolators
        //------------------------------------------------------------------------------
        // create a space geometry interpolation rule
        Interpolation_Rule tGIRule( tGeometryType,
                                    Interpolation_Type::LAGRANGE,
                                    mtk::Interpolation_Order::LINEAR,
                                    Interpolation_Type::LAGRANGE,
                                    mtk::Interpolation_Order::LINEAR );

        // create a space time geometry interpolator
        Geometry_Interpolator tGI = Geometry_Interpolator( tGIRule );

        // create time coeff tHat
        Matrix< DDRMat > tTHat = {{ 0.0 }, { 1.0 }};

        // set the coefficients xHat, tHat
        tGI.set_coeff( tXHat, tTHat );

        // set the evaluation point
        tGI.set_space_time( tParamPoint );

        for( uint iInterpOrder = 1; iInterpOrder < 4; iInterpOrder++ )
        {

            // field interpolators
            //------------------------------------------------------------------------------
            // create an interpolation order
            mtk::Interpolation_Order tInterpolationOrder;

            // create random coefficients for master FI
            arma::Mat< double > tMasterMatrix;
            arma::Mat< double > tSlaveMatrix;

            // switch on interpolation order
            switch( iInterpOrder )
            {
                case ( 1 ):
                {
                    // set interpolation type
                    tInterpolationOrder = mtk::Interpolation_Order::LINEAR;

                    // set interpolation order for IWG
                    tIWG.set_interpolation_order( 1 );

                    // create random coefficients for master FI
                    tMasterMatrix.randu( tNumCoeffs( 0 ), 1 );
                    // create random coefficients for slave FI
                    tSlaveMatrix.randu( tNumCoeffs( 0 ), 1 );
                    break;
                }
                case ( 2 ):
                {
                    // set interpolation type
                    tInterpolationOrder = mtk::Interpolation_Order::QUADRATIC;

                    // set interpolation order for IWG
                    tIWG.set_interpolation_order( 2 );

                    // create random coefficients for master FI
                    tMasterMatrix.randu( tNumCoeffs( 1 ), 1 );
                    // create random coefficients for slave FI
                    tSlaveMatrix.randu( tNumCoeffs( 1 ), 1 );
                    break;
                }
                case ( 3 ):
                {
                    // set interpolation type
                    tInterpolationOrder = mtk::Interpolation_Order::CUBIC;

                    // set interpolation order for IWG
                    tIWG.set_interpolation_order( 3 );

                    // create random coefficients for master FI
                    tMasterMatrix.randu( tNumCoeffs( 2 ), 1 );
                    // create random coefficients for slave FI
                    tSlaveMatrix.randu( tNumCoeffs( 2 ), 1 );
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "LINEAR, QUADRATIC or CUBIC only.");
                    break;
                }
            }

            //create a space time interpolation rule
            Interpolation_Rule tFIRule ( tGeometryType,
                                         Interpolation_Type::LAGRANGE,
                                         tInterpolationOrder,
                                         Interpolation_Type::CONSTANT,
                                         mtk::Interpolation_Order::CONSTANT );

            // fill random coefficients for master FI
            Matrix< DDRMat > tMasterDOFHat;
            tMasterDOFHat.matrix_data() = 10.0 * tMasterMatrix;

            // fill random coefficients for slave FI
            Matrix< DDRMat > tSlaveDOFHat;
            tSlaveDOFHat.matrix_data() = 10.0 * tSlaveMatrix;

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tMasterFIs( 1 );

            // get the number of field
            uint tMasterNumOfFields = 1;

            // create the field interpolator
            tMasterFIs( 0 ) = new Field_Interpolator( tMasterNumOfFields,
                                                      tFIRule,
                                                      &tGI,
                                                      {MSI::Dof_Type::TEMP} );

            // set the coefficients uHat
            tMasterFIs( 0 )->set_coeff( tMasterDOFHat );

            //set the evaluation point xi, tau
            tMasterFIs( 0 )->set_space_time( tParamPoint );

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tSlaveFIs( 1 );

            // get the number of fields
            uint tSlaveNumOfFields = 1;

            // create the field interpolator
            tSlaveFIs( 0 ) = new Field_Interpolator( tSlaveNumOfFields,
                                                     tFIRule,
                                                     &tGI,
                                                     {MSI::Dof_Type::TEMP} );

            // set the coefficients uHat
            tSlaveFIs( 0 )->set_coeff( tSlaveDOFHat );

            //set the evaluation point xi, tau
            tSlaveFIs( 0 )->set_space_time( tParamPoint );

            // properties
            //------------------------------------------------------------------------------
            // create property coefficients
            Cell< Matrix< DDRMat > > tPropCoeff = { {{1.0}} };

            // create a cell of properties for IWG
            Cell< Property* > tMasterProps( 1 );

            // create a property
            tMasterProps( 0 ) = new Property( fem::Property_Type::CONDUCTIVITY,
                                              {{ MSI::Dof_Type::TEMP }},
                                              tPropCoeff,
                                              tFIValFunction_UTGhost,
                                              { tFIDerFunction_UTGhost },
                                              &tGI );

            // set IWG field interpolators
            tMasterProps( 0 )->set_dof_field_interpolators( tMasterFIs );

            // set IWG properties
            Cell< Constitutive_Model* > tMasterCMs;
            tIWG.set_constitutive_models( tMasterCMs );
            Cell< Constitutive_Model* > tSlaveCMs;
            tIWG.set_constitutive_models( tSlaveCMs, mtk::Master_Slave::SLAVE );

            // set IWG properties
            tIWG.set_properties( tMasterProps );
            Cell< Property* > tSlaveProps;
            tIWG.set_properties( tSlaveProps, mtk::Master_Slave::SLAVE );

            // set IWG field interpolators
            tIWG.set_dof_field_interpolators( tMasterFIs );
            tIWG.set_dof_field_interpolators( tSlaveFIs, mtk::Master_Slave::SLAVE );

            // check evaluation of the residual
            //------------------------------------------------------------------------------
            // evaluate the residual
            Cell< Matrix< DDRMat > > tResidual;
            tIWG.compute_residual( tResidual );

            // check evaluation of the jacobian  by FD
            //------------------------------------------------------------------------------
            // init the jacobian for IWG and FD evaluation
            Cell< Cell< Matrix< DDRMat > > > tJacobians;
            Cell< Cell< Matrix< DDRMat > > > tJacobiansFD;

            // check jacobian by FD
            bool tCheckJacobian = tIWG.check_jacobian_double( tPerturbation,
                                                              tEpsilon,
                                                              tJacobians,
                                                              tJacobiansFD );

//            // print for debug
//            print( tJacobians( 0 )( 0 ),"tJacobians00");
//            print( tJacobiansFD( 0 )( 0 ),"tJacobiansFD00");
//
//            print( tJacobians( 0 )( 1 ),"tJacobians01");
//            print( tJacobiansFD( 0 )( 1 ),"tJacobiansFD01");
//
//            print( tJacobians( 1 )( 0 ),"tJacobians10");
//            print( tJacobiansFD( 1 )( 0 ),"tJacobiansFD10");
//
//            print( tJacobians( 1 )( 1 ),"tJacobians11");
//            print( tJacobiansFD( 1 )( 1 ),"tJacobiansFD11");

            // print the treated case
            std::cout<<"Case: Geometry "<<static_cast<int>(tGeometryType)<<" Order "<<static_cast<int>(tInterpolationOrder)<<std::endl;

            // require check is true
            REQUIRE( tCheckJacobian );

            // clean up
            for( Property* tProp : tMasterProps )
            {
                delete tProp;
            }
            tMasterProps.clear();

            for( Field_Interpolator* tFI : tMasterFIs )
            {
                delete tFI;
            }
            tMasterFIs.clear();

            for( Field_Interpolator* tFI : tSlaveFIs )
            {
                delete tFI;
            }
            tSlaveFIs.clear();
        }
    }

}/* END_TEST_CASE */
