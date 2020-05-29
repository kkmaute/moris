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

void tConstValFunction_NSVBULK
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 );
}

void tTEMPFIValFunction_NSVBULK
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 ) * aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::TEMP )->val();
}

void tTEMPFIDerFunction_NSVBULK
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 ) * aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::TEMP )->N();
}

void tPFIValFunction_NSVBULK
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 ) * aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::P )->val();
}

void tPFIDerFunction_NSVBULK
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 ) * aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::P )->N();
}

using namespace moris;
using namespace fem;

TEST_CASE( "IWG_Incompressible_NS_Velocity_Bulk", "[IWG_Incompressible_NS_Velocity_Bulk]" )
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

        // dof type list
        moris::Cell< MSI::Dof_Type > tVelDofTypes;
        moris::Cell< MSI::Dof_Type > tPDofTypes    = { MSI::Dof_Type::P };
        moris::Cell< MSI::Dof_Type > tTEMPDofTypes = { MSI::Dof_Type::TEMP };
        moris::Cell< MSI::Dof_Type > tVisDofTypes  = { MSI::Dof_Type::VISCOSITY };

        // gravity
        Matrix< DDRMat > tGravity( iSpaceDim, 1, 10.0 );

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
            arma::Mat< double > tMasterMatrixTEMP;
            arma::Mat< double > tMasterMatrixVis;

            // get number of dof
            int tNumDofVel  = 0;
            int tNumDofP    = 0;
            int tNumDofTEMP = 0;
            int tNumDofVis  = 0;

            // switch on interpolation order
            switch( iInterpOrder )
            {
                case 1 :
                {
                    // set interpolation type
                    tInterpolationOrder = mtk::Interpolation_Order::LINEAR;

                    // get number of dof
                    tNumDofVel  = tNumCoeffs( 0 ) * iSpaceDim;
                    tNumDofP    = tNumCoeffs( 0 );
                    tNumDofTEMP = tNumCoeffs( 0 );
                    tNumDofVis  = tNumCoeffs( 0 );

                    // create random coefficients for master FI
                    tMasterMatrixVel.randu( tNumCoeffs( 0 ), iSpaceDim );
                    tMasterMatrixP.randu( tNumCoeffs( 0 ), 1 );
                    tMasterMatrixTEMP.randu( tNumCoeffs( 0 ), 1 );
                    tMasterMatrixVis.randu( tNumCoeffs( 0 ), 1 );
                    break;
                }
                case 2 :
                {
                    // set interpolation type
                    tInterpolationOrder = mtk::Interpolation_Order::QUADRATIC;

                    // get number of dof
                    tNumDofVel = tNumCoeffs( 1 ) * iSpaceDim;
                    tNumDofP = tNumCoeffs( 1 );
                    tNumDofTEMP = tNumCoeffs( 1 );
                    tNumDofVis = tNumCoeffs( 1 );

                    // create random coefficients for master FI
                    tMasterMatrixVel.randu( tNumCoeffs( 1 ), iSpaceDim );
                    tMasterMatrixP.randu( tNumCoeffs( 1 ), 1 );
                    tMasterMatrixTEMP.randu( tNumCoeffs( 1 ), 1 );
                    tMasterMatrixVis.randu( tNumCoeffs( 1 ), 1 );
                    break;
                }
                case 3 :
                {
                    // set interpolation type
                    tInterpolationOrder = mtk::Interpolation_Order::CUBIC;

                    // get number of dof
                    tNumDofVel = tNumCoeffs( 2 ) * iSpaceDim;
                    tNumDofP = tNumCoeffs( 2 );
                    tNumDofTEMP = tNumCoeffs( 2 );
                    tNumDofVis = tNumCoeffs( 2 );

                    // create random coefficients for master FI
                    tMasterMatrixVel.randu( tNumCoeffs( 2 ), iSpaceDim );
                    tMasterMatrixP.randu( tNumCoeffs( 2 ), 1 );
                    tMasterMatrixTEMP.randu( tNumCoeffs( 2 ), 1 );
                    tMasterMatrixVis.randu( tNumCoeffs( 2 ), 1 );
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
            Matrix< DDRMat > tMasterDOFHatTEMP;
            tMasterDOFHatTEMP.matrix_data() = 10.0 * tMasterMatrixTEMP;
            Matrix< DDRMat > tMasterDOFHatVis;
            tMasterDOFHatVis.matrix_data() = 10.0 * tMasterMatrixVis;

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tMasterFIs( 4 );

            // create the field interpolator velocity
            tMasterFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelDofTypes );
            tMasterFIs( 0 )->set_coeff( tMasterDOFHatVel );
            tMasterFIs( 0 )->set_space_time( tParamPoint );

            // create the field interpolator pressure
            tMasterFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tPDofTypes );
            tMasterFIs( 1 )->set_coeff( tMasterDOFHatP );
            tMasterFIs( 1 )->set_space_time( tParamPoint );

            // create the field interpolator temperature
            tMasterFIs( 2 ) = new Field_Interpolator( 1, tFIRule, &tGI, tTEMPDofTypes );
            tMasterFIs( 2 )->set_coeff( tMasterDOFHatTEMP );
            tMasterFIs( 2 )->set_space_time( tParamPoint );

            // create the field interpolator viscosity
            tMasterFIs( 3 ) = new Field_Interpolator( 1, tFIRule, &tGI, tVisDofTypes );
            tMasterFIs( 3 )->set_coeff( tMasterDOFHatVis );
            tMasterFIs( 3 )->set_space_time( tParamPoint );

            // create the properties
            std::shared_ptr< fem::Property > tPropViscosity = std::make_shared< fem::Property >();
            tPropViscosity->set_parameters( { {{ 1.0 }} } );
            tPropViscosity->set_dof_type_list( { tPDofTypes } );
            tPropViscosity->set_val_function( tPFIValFunction_NSVBULK );
            tPropViscosity->set_dof_derivative_functions( { tPFIDerFunction_NSVBULK } );

            std::shared_ptr< fem::Property > tPropDensity = std::make_shared< fem::Property >();
            tPropDensity->set_parameters( { {{ 2.0 }} } );
            tPropDensity->set_dof_type_list( { tPDofTypes } );
            tPropDensity->set_val_function( tPFIValFunction_NSVBULK );
            tPropDensity->set_dof_derivative_functions( { tPFIDerFunction_NSVBULK } );

            std::shared_ptr< fem::Property > tPropGravity = std::make_shared< fem::Property >();
            tPropGravity->set_parameters( { tGravity } );
            tPropGravity->set_dof_type_list( { tPDofTypes } );
            tPropGravity->set_val_function( tPFIValFunction_NSVBULK );
            tPropGravity->set_dof_derivative_functions( { tPFIDerFunction_NSVBULK } );

            std::shared_ptr< fem::Property > tPropThermalExp = std::make_shared< fem::Property >();
            tPropThermalExp->set_parameters( { {{ 23.0 }} } );
            tPropThermalExp->set_dof_type_list( { tTEMPDofTypes } );
            tPropThermalExp->set_val_function( tTEMPFIValFunction_NSVBULK );
            tPropThermalExp->set_dof_derivative_functions( { tTEMPFIDerFunction_NSVBULK } );

            std::shared_ptr< fem::Property > tPropRefTemp = std::make_shared< fem::Property >();
            tPropRefTemp->set_parameters( { {{ 15.0 }} } );
            tPropRefTemp->set_dof_type_list( { tTEMPDofTypes } );
            tPropRefTemp->set_val_function( tTEMPFIValFunction_NSVBULK );
            tPropRefTemp->set_dof_derivative_functions( { tTEMPFIDerFunction_NSVBULK } );

            // define constitutive models
            fem::CM_Factory tCMFactory;

            std::shared_ptr< fem::Constitutive_Model > tCMMasterIncFluid =
                    tCMFactory.create_CM( fem::Constitutive_Type::FLUID_INCOMPRESSIBLE );
            tCMMasterIncFluid->set_dof_type_list( { tVelDofTypes, tPDofTypes } );
            tCMMasterIncFluid->set_property( tPropViscosity, "Viscosity" );
            tCMMasterIncFluid->set_property( tPropDensity, "Density" );
            tCMMasterIncFluid->set_space_dim( iSpaceDim );

            std::shared_ptr< fem::Constitutive_Model > tCMMasterTurbulence =
                              tCMFactory.create_CM( fem::Constitutive_Type::FLUID_TURBULENCE );
            tCMMasterTurbulence->set_dof_type_list( { tVelDofTypes, tVisDofTypes } );
            tCMMasterTurbulence->set_property( tPropViscosity, "Viscosity" );
            tCMMasterTurbulence->set_property( tPropDensity, "Density" );
            tCMMasterTurbulence->set_space_dim( iSpaceDim );

            // define stabilization parameters
            fem::SP_Factory tSPFactory;

            std::shared_ptr< fem::Stabilization_Parameter > tSPIncFlow =
                    tSPFactory.create_SP( fem::Stabilization_Type::INCOMPRESSIBLE_FLOW );
            tSPIncFlow->set_dof_type_list( { tVelDofTypes, tPDofTypes }, mtk::Master_Slave::MASTER );
            tSPIncFlow->set_property( tPropDensity, "Density", mtk::Master_Slave::MASTER );
            tSPIncFlow->set_property( tPropViscosity, "Viscosity", mtk::Master_Slave::MASTER );
            tSPIncFlow->set_parameters( { {{ 36.0 }} } );
            tSPIncFlow->set_space_dim( iSpaceDim );

            // define the IWGs
            fem::IWG_Factory tIWGFactory;

            std::shared_ptr< fem::IWG > tIWG = tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_BULK );
            tIWG->set_residual_dof_type( tVelDofTypes );
            tIWG->set_dof_type_list( { tVelDofTypes, tPDofTypes, tTEMPDofTypes, tVisDofTypes }, mtk::Master_Slave::MASTER );
            tIWG->set_property( tPropDensity, "Density" );
            tIWG->set_property( tPropGravity, "Gravity" );
            tIWG->set_property( tPropThermalExp, "ThermalExpansion" );
            tIWG->set_property( tPropRefTemp, "ReferenceTemp" );
            tIWG->set_constitutive_model( tCMMasterIncFluid, "IncompressibleFluid" );
            tIWG->set_constitutive_model( tCMMasterTurbulence, "TurbulenceFluid" );
            tIWG->set_stabilization_parameter( tSPIncFlow, "IncompressibleFlow" );

            // set a fem set pointer
            MSI::Equation_Set * tSet = new fem::Set();
            tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

            // set size for the set EqnObjDofTypeList
            tIWG->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

            // set size and populate the set dof type map
            tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
            tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )        = 0;
            tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )         = 1;
            tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) )      = 2;
            tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 3;

            // set size and populate the set master dof type map
            tIWG->mSet->mMasterDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
            tIWG->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )        = 0;
            tIWG->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )         = 1;
            tIWG->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) )      = 2;
            tIWG->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 3;

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( 4 );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofVel-1 } };
            tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofVel, tNumDofVel + tNumDofP - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 2 ) = { { tNumDofVel + tNumDofP, tNumDofVel + tNumDofP + tNumDofTEMP - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 3 ) = { { tNumDofVel + tNumDofP + tNumDofTEMP, tNumDofVel + tNumDofP + tNumDofTEMP + tNumDofVis - 1 } };

            // set size and fill the set jacobian assembly map
            tIWG->mSet->mJacDofAssemblyMap.resize( 4 );
            tIWG->mSet->mJacDofAssemblyMap( 0 ) = {
                    { 0, tNumDofVel - 1 },
                    { tNumDofVel, tNumDofVel + tNumDofP - 1 },
                    { tNumDofVel + tNumDofP, tNumDofVel + tNumDofP + tNumDofTEMP - 1 },
                    { tNumDofVel + tNumDofP + tNumDofTEMP, tNumDofVel + tNumDofP + tNumDofTEMP + tNumDofVis - 1 } };
            tIWG->mSet->mJacDofAssemblyMap( 1 ) = {
                    { 0, tNumDofVel - 1 },
                    { tNumDofVel, tNumDofVel + tNumDofP - 1 },
                    { tNumDofVel + tNumDofP, tNumDofVel + tNumDofP + tNumDofTEMP - 1 },
                    { tNumDofVel + tNumDofP + tNumDofTEMP, tNumDofVel + tNumDofP + tNumDofTEMP + tNumDofVis - 1 } };
            tIWG->mSet->mJacDofAssemblyMap( 2 ) = {
                    { 0, tNumDofVel - 1 },
                    { tNumDofVel, tNumDofVel + tNumDofP - 1 },
                    { tNumDofVel + tNumDofP, tNumDofVel + tNumDofP + tNumDofTEMP - 1 },
                    { tNumDofVel + tNumDofP + tNumDofTEMP, tNumDofVel + tNumDofP + tNumDofTEMP + tNumDofVis - 1 } };
            tIWG->mSet->mJacDofAssemblyMap( 3 ) = {
                    { 0, tNumDofVel - 1 },
                    { tNumDofVel, tNumDofVel + tNumDofP - 1 },
                    { tNumDofVel + tNumDofP, tNumDofVel + tNumDofP + tNumDofTEMP - 1 },
                    { tNumDofVel + tNumDofP + tNumDofTEMP, tNumDofVel + tNumDofP + tNumDofTEMP + tNumDofVis - 1 } };

            // set size and init the set residual and jacobian
            tIWG->mSet->mResidual.resize( 1 );
            tIWG->mSet->mResidual( 0 ).set_size(
                    tNumDofVel + tNumDofP + tNumDofTEMP + tNumDofVis,
                    1,
                    0.0 );
            tIWG->mSet->mJacobian.set_size(
                    tNumDofVel + tNumDofP + tNumDofTEMP + tNumDofVis,
                    tNumDofVel + tNumDofP + tNumDofTEMP + tNumDofVis,
                    0.0 );

            // build global dof type list
            tIWG->get_global_dof_type_list();

            // populate the requested master dof type
            tIWG->mRequestedMasterGlobalDofTypes = {
                    { MSI::Dof_Type::VX },
                    { MSI::Dof_Type::P },
                    { MSI::Dof_Type::TEMP },
                    { MSI::Dof_Type::VISCOSITY } };

            // create a field interpolator manager
            moris::Cell< moris::Cell< enum MSI::Dof_Type > > tDummyDof;
            moris::Cell< moris::Cell< enum PDV_Type > > tDummyDv;
            Field_Interpolator_Manager tFIManager( tDummyDof, tDummyDv, tSet );

            // populate the field interpolator manager
            tFIManager.mFI = tMasterFIs;
            tFIManager.mIPGeometryInterpolator = &tGI;
            tFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tIWG->mSet->mMasterFIManager = &tFIManager;

            // set IWG field interpolator manager
            tIWG->set_field_interpolator_manager( &tFIManager );

            // check evaluation of the residual for IWG
            //------------------------------------------------------------------------------
            // evaluate the residual
            tIWG->compute_residual( 1.0 );

            // check evaluation of the jacobian by FD
            //------------------------------------------------------------------------------
            // init the jacobian for IWG and FD evaluation
            Matrix< DDRMat > tJacobian;
            Matrix< DDRMat > tJacobianFD;

            // check jacobian by FD
            bool tCheckJacobian = tIWG->check_jacobian( tPerturbation,
                                                        tEpsilon,
                                                        1.0,
                                                        tJacobian,
                                                        tJacobianFD );

//            print( tJacobian(   { 0, tNumDofVel-1 }, { 0, tNumDofVel-1 } ), "tJacobianVV" );
//            print( tJacobianFD( { 0, tNumDofVel-1 }, { 0, tNumDofVel-1 } ), "tJacobianFDVv" );
//            print( tJacobian(   { 0, tNumDofVel-1 }, { tNumDofVel, tNumDofVel + tNumDofP - 1 }), "tJacobianVP" );
//            print( tJacobianFD( { 0, tNumDofVel-1 }, { tNumDofVel, tNumDofVel + tNumDofP - 1 }), "tJacobianFDVP" );
//            print( tJacobian(   { 0, tNumDofVel-1 }, { tNumDofVel + tNumDofP, tNumDofVel + tNumDofP + tNumDofTEMP - 1 }), "tJacobianVTEMP" );
//            print( tJacobianFD( { 0, tNumDofVel-1 }, { tNumDofVel + tNumDofP, tNumDofVel + tNumDofP + tNumDofTEMP - 1 }), "tJacobianFDVTEMP" );
//            print( tJacobian(   { 0, tNumDofVel-1 }, { tNumDofVel + tNumDofP + tNumDofTEMP, tNumDofVel + tNumDofP + tNumDofTEMP + tNumDofVis - 1 }), "tJacobianVTEMP" );
//            print( tJacobianFD( { 0, tNumDofVel-1 }, { tNumDofVel + tNumDofP + tNumDofTEMP, tNumDofVel + tNumDofP + tNumDofTEMP + tNumDofVis - 1 }), "tJacobianFDVTEMP" );

//            std::cout<<"Case: Geometry "<<iSpaceDim<<" Order "<<iInterpOrder<<std::endl;

            // require check is true
            REQUIRE( tCheckJacobian );

            // clean up
            tMasterFIs.clear();
        }
    }
}/*END_TEST_CASE*/
