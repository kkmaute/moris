#include <string>
#include <catch.hpp>
#include <memory>
#include "assert.hpp"

#define protected public
#define private   public
//FEM//INT//src
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG.hpp"
#include "cl_FEM_Set.hpp"
#undef protected
#undef private

//LINALG/src
#include "op_equal_equal.hpp"
//MTK/src
#include "cl_MTK_Enums.hpp"
//FEM//INT//src
#include "cl_FEM_Enums.hpp"
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_FEM_Property.hpp"
#include "cl_FEM_CM_Factory.hpp"
#include "cl_FEM_SP_Factory.hpp"
#include "cl_FEM_IWG_Factory.hpp"

void tConstValFunction_UTVelocityDirichlet
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 );
}

using namespace moris;
using namespace fem;

TEST_CASE( "IWG_Incompressible_NS_Dirichlet_Nitsche", "[IWG_Incompressible_NS_Dirichlet_Nitsche]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-3;

    // define a perturbation relative size
    real tPerturbation = 1E-4;

    // loop on the space dimension
    for( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // set geometry inputs
        //------------------------------------------------------------------------------
        // create geometry type
        mtk::Geometry_Type tGeometryType = mtk::Geometry_Type::UNDEFINED;

        // create space coeff xHat
        Matrix< DDRMat > tXHat;

        // create evaluation point xi, tau
        Matrix< DDRMat > tParamPoint;

        // create list with number of coeffs
        Matrix< DDRMat > tNumCoeffs;

        // create the normal
        Matrix< DDRMat > tNormal;

        // dof type list
        moris::Cell< MSI::Dof_Type > tVelDofTypes;
        moris::Cell< MSI::Dof_Type > tPDofTypes = { MSI::Dof_Type::P };

        // gravity
        Matrix< DDRMat > tVelocity( iSpaceDim, 1, 10.0 );

        // switch on space dimension
        switch( iSpaceDim )
        {
            case 2 :
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
               tNumCoeffs = {{ 8 },{ 18 },{ 32 }};


               // set the normal
               tNormal = {{ 1.0 }, { 0.0 }};

               // set velocity dof types
               tVelDofTypes = { MSI::Dof_Type::VX, MSI::Dof_Type::VY };

               break;
            }
            case 3 :
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
                tParamPoint = {{ 0.35 }, {-0.25}, { 0.75}, { 0.0 }};

                // number of coefficients
                tNumCoeffs = {{ 16 },{ 54 },{ 128 }};


                // set the normal
                tNormal = {{ 1.0 }, { 0.0 }, { 0.0 }};

                // set velocity dof types
                tVelDofTypes = { MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::VZ };

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
        for( uint iInterpOrder = 1; iInterpOrder < 4; iInterpOrder++ )
        {
            // field interpolators
            //------------------------------------------------------------------------------
            // create an interpolation order
            mtk::Interpolation_Order tInterpolationOrder;

            // create random coefficients for master FI
            arma::Mat< double > tMasterMatrixVel;
            arma::Mat< double > tMasterMatrixP;

            // get number of dof
            int tNumDofVel = 0;
            int tNumDofP = 0;

            // switch on interpolation order
            switch( iInterpOrder )
            {
                case 1 :
                {
                    // set interpolation type
                    tInterpolationOrder = mtk::Interpolation_Order::LINEAR;

                    // get number of dof
                    tNumDofVel = tNumCoeffs( 0 ) * iSpaceDim;
                    tNumDofP = tNumCoeffs( 0 );

                    // create random coefficients for master FI
                    tMasterMatrixVel.randu( tNumCoeffs( 0 ), iSpaceDim );
                    tMasterMatrixP.randu( tNumCoeffs( 0 ), 1 );
                    break;
                }
                case 2 :
                {
                    // set interpolation type
                    tInterpolationOrder = mtk::Interpolation_Order::QUADRATIC;

                    // get number of dof
                    tNumDofVel = tNumCoeffs( 1 ) * iSpaceDim;
                    tNumDofP = tNumCoeffs( 1 );

                    // create random coefficients for master FI
                    tMasterMatrixVel.randu( tNumCoeffs( 1 ), iSpaceDim );
                    tMasterMatrixP.randu( tNumCoeffs( 1 ), 1 );
                    break;
                }
                case 3 :
                {
                    // set interpolation type
                    tInterpolationOrder = mtk::Interpolation_Order::CUBIC;

                    // get number of dof
                    tNumDofVel = tNumCoeffs( 2 ) * iSpaceDim;
                    tNumDofP = tNumCoeffs( 2 );

                    // create random coefficients for master FI
                    tMasterMatrixVel.randu( tNumCoeffs( 2 ), iSpaceDim );
                    tMasterMatrixP.randu( tNumCoeffs( 2 ), 1 );
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
                                         Interpolation_Type::LAGRANGE,
                                         mtk::Interpolation_Order::LINEAR );

            // fill random coefficients for master FI
            Matrix< DDRMat > tMasterDOFHatVel;
            tMasterDOFHatVel.matrix_data() = 10.0 * tMasterMatrixVel;
            Matrix< DDRMat > tMasterDOFHatP;
            tMasterDOFHatP.matrix_data() = 10.0 * tMasterMatrixP;

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tMasterFIs( 2 );

            // create the field interpolator velocity
            tMasterFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelDofTypes );
            tMasterFIs( 0 )->set_coeff( tMasterDOFHatVel );
            tMasterFIs( 0 )->set_space_time( tParamPoint );

            // create the field interpolator pressure
            tMasterFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tPDofTypes );
            tMasterFIs( 1 )->set_coeff( tMasterDOFHatP );
            tMasterFIs( 1 )->set_space_time( tParamPoint );

            // create the properties
            std::shared_ptr< fem::Property > tPropViscosity = std::make_shared< fem::Property >();
            tPropViscosity->set_parameters( { {{ 1.0 }} } );
            tPropViscosity->set_val_function( tConstValFunction_UTVelocityDirichlet );

            std::shared_ptr< fem::Property > tPropDensity = std::make_shared< fem::Property >();
            tPropDensity->set_parameters( { {{ 1.0 }} } );
            tPropDensity->set_val_function( tConstValFunction_UTVelocityDirichlet );

            std::shared_ptr< fem::Property > tPropVelocity = std::make_shared< fem::Property >();
            tPropVelocity->set_parameters( { tVelocity } );
            tPropVelocity->set_val_function( tConstValFunction_UTVelocityDirichlet );

            // define constitutive models
            fem::CM_Factory tCMFactory;

            std::shared_ptr< fem::Constitutive_Model > tCMMasterIncFluid =
                    tCMFactory.create_CM( fem::Constitutive_Type::FLUID_INCOMPRESSIBLE );
            tCMMasterIncFluid->set_dof_type_list( { tVelDofTypes, tPDofTypes } );
            tCMMasterIncFluid->set_property( tPropViscosity, "Viscosity" );
            tCMMasterIncFluid->set_property( tPropDensity, "Density" );
            tCMMasterIncFluid->set_space_dim( iSpaceDim );

            // define stabilization parameters
            fem::SP_Factory tSPFactory;

            std::shared_ptr< fem::Stabilization_Parameter > tSPNitsche =
                    tSPFactory.create_SP( fem::Stabilization_Type::VELOCITY_DIRICHLET_NITSCHE );
            tSPNitsche->set_dof_type_list( { tVelDofTypes }, mtk::Master_Slave::MASTER );
            tSPNitsche->set_property( tPropDensity, "Density", mtk::Master_Slave::MASTER );
            tSPNitsche->set_property( tPropViscosity, "Viscosity", mtk::Master_Slave::MASTER );
            tSPNitsche->set_parameters( { {{ 1.0 }}, {{ 1.0 }} } );

            // define the IWGs
            fem::IWG_Factory tIWGFactory;

            std::shared_ptr< fem::IWG > tIWGVelocity =
                    tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_SYMMETRIC_NITSCHE );
            tIWGVelocity->set_residual_dof_type( tVelDofTypes );
            tIWGVelocity->set_dof_type_list( { tVelDofTypes, tPDofTypes}, mtk::Master_Slave::MASTER );
            tIWGVelocity->set_property( tPropVelocity, "Dirichlet" );
            tIWGVelocity->set_constitutive_model( tCMMasterIncFluid, "IncompressibleFluid" );
            tIWGVelocity->set_stabilization_parameter( tSPNitsche, "DirichletNitsche" );

            std::shared_ptr< fem::IWG > tIWGPressure =
                    tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_SYMMETRIC_NITSCHE );
            tIWGPressure->set_residual_dof_type( tPDofTypes );
            tIWGPressure->set_dof_type_list( { tVelDofTypes, tPDofTypes }, mtk::Master_Slave::MASTER );
            tIWGPressure->set_property( tPropVelocity, "Dirichlet" );
            tIWGPressure->set_constitutive_model( tCMMasterIncFluid, "IncompressibleFluid" );

            // set a fem set pointer
            MSI::Equation_Set * tSet = new fem::Set();
            tIWGVelocity->set_set_pointer( static_cast< fem::Set* >( tSet ) );
            tIWGPressure->set_set_pointer( static_cast< fem::Set* >( tSet ) );

            // set size for the set EqnObjDofTypeList
            tIWGVelocity->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

            // set size and populate the set dof type map
            tIWGVelocity->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
            tIWGVelocity->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )   = 0;
            tIWGVelocity->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )    = 1;

            // set size and populate the set master dof type map
            tIWGVelocity->mSet->mMasterDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
            tIWGVelocity->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )   = 0;
            tIWGVelocity->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )    = 1;

            // set size and fill the set residual assembly map
            tIWGVelocity->mSet->mResDofAssemblyMap.resize( 2 );
            tIWGVelocity->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofVel-1 } };
            tIWGVelocity->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofVel, tNumDofVel + tNumDofP - 1 } };

            // set size and fill the set jacobian assembly map
            tIWGVelocity->mSet->mJacDofAssemblyMap.resize( 2 );
            tIWGVelocity->mSet->mJacDofAssemblyMap( 0 ) = { { 0, tNumDofVel - 1 }, { tNumDofVel, tNumDofVel + tNumDofP - 1 } };
            tIWGVelocity->mSet->mJacDofAssemblyMap( 1 ) = { { 0, tNumDofVel - 1 }, { tNumDofVel, tNumDofVel + tNumDofP - 1 } };

            // set size and init the set residual and jacobian
            tIWGVelocity->mSet->mResidual.resize( 1 );
            tIWGVelocity->mSet->mResidual( 0 ).set_size( tNumDofVel + tNumDofP, 1, 0.0 );
            tIWGVelocity->mSet->mJacobian.set_size( tNumDofVel + tNumDofP, tNumDofVel + tNumDofP, 0.0 );

            // set the normal
            tIWGVelocity->set_normal( tNormal );
            tIWGPressure->set_normal( tNormal );

            // build global dof type list
            tIWGVelocity->get_global_dof_type_list();
            tIWGPressure->get_global_dof_type_list();

            // populate the requested master dof type
            tIWGVelocity->mRequestedMasterGlobalDofTypes = { { MSI::Dof_Type::VX }, { MSI::Dof_Type::P } };
            tIWGPressure->mRequestedMasterGlobalDofTypes = { { MSI::Dof_Type::VX }, { MSI::Dof_Type::P } };

            // create a field interpolator manager
            moris::Cell< moris::Cell< enum MSI::Dof_Type > > tDummyDof;
            moris::Cell< moris::Cell< enum PDV_Type > > tDummyDv;
            Field_Interpolator_Manager tFIManager( tDummyDof, tDummyDv, tSet );

            // populate the field interpolator manager
            tFIManager.mFI = tMasterFIs;
            tFIManager.mIPGeometryInterpolator = &tGI;
            tFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tIWGVelocity->mSet->mMasterFIManager = &tFIManager;

            // set IWG field interpolator manager
            tIWGVelocity->set_field_interpolator_manager( &tFIManager );
            tIWGPressure->set_field_interpolator_manager( &tFIManager );

            // check evaluation of the residual/jacobian for IWG Velocity
            //------------------------------------------------------------------------------
            // reset residual and jacobian
            tIWGVelocity->mSet->mResidual( 0 ).fill( 0.0 );
            tIWGVelocity->mSet->mJacobian.fill( 0.0 );

            // evaluate the residual
            tIWGVelocity->compute_residual( 1.0 );

            // init the jacobian for IWG and FD evaluation
            Matrix< DDRMat > tVelocityJacobian;
            Matrix< DDRMat > tVelocityJacobianFD;

            // check jacobian by FD
            bool tCheckVelocityJacobian = tIWGVelocity->check_jacobian(
                    tPerturbation,
                    tEpsilon,
                    1.0,
                    tVelocityJacobian,
                    tVelocityJacobianFD );
//            print( tVelocityJacobian(   { 0, tNumDofVel-1 }, { 0, tNumDofVel-1 } ), "tJacobianVV" );
//            print( tVelocityJacobianFD( { 0, tNumDofVel-1 }, { 0, tNumDofVel-1 } ), "tJacobianFDVV" );
//            print( tVelocityJacobian(   { 0, tNumDofVel-1 }, { tNumDofVel, tNumDofVel + tNumDofP - 1 }), "tJacobianVP" );
//            print( tVelocityJacobianFD( { 0, tNumDofVel-1 }, { tNumDofVel, tNumDofVel + tNumDofP - 1 }), "tJacobianFDVP" );

            // require check is true
            REQUIRE( tCheckVelocityJacobian );

            // check evaluation of the residual/jacobian for IWG Pressure
            //------------------------------------------------------------------------------
            // reset residual and jacobian
            tIWGPressure->mSet->mResidual( 0 ).fill( 0.0 );
            tIWGPressure->mSet->mJacobian.fill( 0.0 );

            // evaluate the residual
            tIWGPressure->compute_residual( 1.0 );

            // init the jacobian for IWG and FD evaluation
            Matrix< DDRMat > tPressureJacobian;
            Matrix< DDRMat > tPressureJacobianFD;

            // check jacobian by FD
            bool tCheckPressureJacobian = tIWGPressure->check_jacobian(
                    tPerturbation,
                    tEpsilon,
                    1.0,
                    tPressureJacobian,
                    tPressureJacobianFD );

//            print( tPressureJacobian(   { 0, tNumDofP-1 }, { 0, tNumDofVel-1 } ), "tJacobianPV" );
//            print( tPressureJacobianFD( { 0, tNumDofP-1 }, { 0, tNumDofVel-1 } ), "tJacobianFDPV" );
//            print( tPressureJacobian(   { 0, tNumDofP-1 }, { tNumDofVel, tNumDofVel + tNumDofP - 1 }), "tJacobianPP" );
//            print( tPressureJacobianFD( { 0, tNumDofP-1 }, { tNumDofVel, tNumDofVel + tNumDofP - 1 }), "tJacobianFDPP" );

//            std::cout<<"Case: Geometry "<<iSpaceDim<<" Order "<<iInterpOrder<<std::endl;

            // require check is true
            REQUIRE( tCheckPressureJacobian );

            // clean up
            tMasterFIs.clear();
        }
    }
}/*END_TEST_CASE*/
