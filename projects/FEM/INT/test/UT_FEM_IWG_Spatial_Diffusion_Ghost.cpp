#include <string>
#include <catch.hpp>
#include "assert.hpp"

#include "cl_MTK_Enums.hpp" //MTK/src
#include "cl_FEM_Enums.hpp"                                 //FEM//INT/src
#include "cl_FEM_IWG_Factory.hpp"                            //FEM//INT/src
#include "cl_FEM_CM_Factory.hpp"                            //FEM//INT/src
#include "cl_FEM_SP_Factory.hpp"                            //FEM//INT/src

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
TEST_CASE( "IWG_Diff_Ghost", "[moris],[fem],[IWG_Diff_Ghost]" )
{

    // define an epsilon environment
    real tEpsilon = 1E-4;

    // define a perturbation relative size
    real tPerturbation = 1E-4;

    // loop on the space dimension
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

        // loop on the interpolation order
        for( uint iInterpOrder = 1; iInterpOrder < 3; iInterpOrder++ )
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

            // create the field interpolator
            tMasterFIs( 0 ) = new Field_Interpolator( 1, tFIRule, &tGI, { MSI::Dof_Type::TEMP } );

            // set the coefficients uHat
            tMasterFIs( 0 )->set_coeff( tMasterDOFHat );

            //set the evaluation point xi, tau
            tMasterFIs( 0 )->set_space_time( tParamPoint );

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tSlaveFIs( 1 );

            // create the field interpolator
            tSlaveFIs( 0 ) = new Field_Interpolator( 1, tFIRule, &tGI, { MSI::Dof_Type::TEMP } );

            // set the coefficients uHat
            tSlaveFIs( 0 )->set_coeff( tSlaveDOFHat );

            //set the evaluation point xi, tau
            tSlaveFIs( 0 )->set_space_time( tParamPoint );

            // create the properties
            std::shared_ptr< fem::Property > tPropMasterConductivity = std::make_shared< fem::Property > ();
            tPropMasterConductivity->set_parameters( { {{ 1.0 }} } );
            tPropMasterConductivity->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tPropMasterConductivity->set_val_function( tFIValFunction_UTGhost );
            tPropMasterConductivity->set_dof_derivative_functions( { tFIDerFunction_UTGhost } );

            // define stabilization parameters
            fem::SP_Factory tSPFactory;

            // define the IWGs
            fem::IWG_Factory tIWGFactory;

            std::shared_ptr< fem::IWG > tIWG = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_GHOST );

            if ( iInterpOrder > 0 )
            {
                std::shared_ptr< fem::Stabilization_Parameter > tSP1 = tSPFactory.create_SP( fem::Stabilization_Type::GHOST_DISPL );
                tSP1->set_parameters( {{{ 1.0 }}, {{ 1.0 }} });
                tSP1->set_property( tPropMasterConductivity, "Material", mtk::Master_Slave::MASTER );
                tIWG->set_stabilization_parameter( tSP1, "GhostDisplOrder1" );
            }
            if ( iInterpOrder > 1 )
            {
                std::shared_ptr< fem::Stabilization_Parameter > tSP2 = tSPFactory.create_SP( fem::Stabilization_Type::GHOST_DISPL );
                tSP2->set_parameters( {{{ 1.0 }}, {{ 2.0 }} });
                tSP2->set_property( tPropMasterConductivity, "Material", mtk::Master_Slave::MASTER );
                tIWG->set_stabilization_parameter( tSP2, "GhostDisplOrder2" );
            }
            if ( iInterpOrder > 2 )
            {
                std::shared_ptr< fem::Stabilization_Parameter > tSP3 = tSPFactory.create_SP( fem::Stabilization_Type::GHOST_DISPL );
                tSP3->set_parameters( {{{ 1.0 }}, {{ 3.0 }} });
                tSP3->set_property( tPropMasterConductivity, "Material", mtk::Master_Slave::MASTER );
                tIWG->set_stabilization_parameter( tSP3, "GhostDisplOrder3" );
            }

            tIWG->set_residual_dof_type( { MSI::Dof_Type::TEMP } );
            tIWG->set_dof_type_list( {{ MSI::Dof_Type::TEMP }}, mtk::Master_Slave::MASTER );
            tIWG->set_dof_type_list( {{ MSI::Dof_Type::TEMP }}, mtk::Master_Slave::SLAVE );

            // set IWG normal
            tIWG->set_normal( tNormal );

            // build global property type list
            tIWG->build_global_dof_type_list();

            // set IWG field interpolators
            tIWG->set_dof_field_interpolators( tMasterFIs );
            tIWG->set_dof_field_interpolators( tSlaveFIs, mtk::Master_Slave::SLAVE );

            // set IWG geometry interpolator
            tIWG->set_geometry_interpolator( &tGI );
            tIWG->set_geometry_interpolator( &tGI, mtk::Master_Slave::SLAVE );

            // check evaluation of the residual
            //------------------------------------------------------------------------------
            // evaluate the residual
            Cell< Matrix< DDRMat > > tResidual;
            tIWG->compute_residual( tResidual );

            // check evaluation of the jacobian  by FD
            //------------------------------------------------------------------------------
            // init the jacobian for IWG and FD evaluation
            Cell< Cell< Matrix< DDRMat > > > tJacobians;
            Cell< Cell< Matrix< DDRMat > > > tJacobiansFD;

            // check jacobian by FD
            bool tCheckJacobian = tIWG->check_jacobian_double( tPerturbation,
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

//            // print the treated case
//            std::cout<<"Case: Geometry "<<iSpaceDim<<" Order "<<iInterpOrder<<std::endl;


            // require check is true
            REQUIRE( tCheckJacobian );

            // clean up
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
