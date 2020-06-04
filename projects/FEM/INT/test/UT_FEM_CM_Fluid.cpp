
#include "catch.hpp"

#define protected public
#define private   public
//FEM/INT/src
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Constitutive_Model.hpp"
#include "cl_FEM_Set.hpp"
#undef protected
#undef private

//LINALG/src
#include "fn_equal_to.hpp"
//FEM/INT/src
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_FEM_Property.hpp"
#include "cl_FEM_CM_Factory.hpp"

void tValFunc_CMFluid(
        moris::Matrix< moris::DDRMat >                 & aPropMatrix,
        moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
        moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix =
            aParameters( 0 )
            + aParameters( 0 ) * aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::TEMP )->val();
}

void tConstValFunc_CMFluid(
        moris::Matrix< moris::DDRMat >                 & aPropMatrix,
        moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
        moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 );
}

void tDerFunc_CMFluid(
        moris::Matrix< moris::DDRMat >                 & aPropMatrix,
        moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
        moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix =
            aParameters( 0 ) *
            aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::TEMP )->N();
}

bool tCheckFunc(
        moris::Matrix< moris::DDRMat > & aMatrixCheck,
        moris::Matrix< moris::DDRMat > & aMatrixRef,
        moris::real                      aEpsilon )
{
    //define a boolean for check
    bool tCheckBool = true;

    // define a real for absolute difference
    moris::real tAbsolute = 0.0;

    // define a real for relative difference
    moris::real tRelative = 0.0;

    for( uint ii = 0; ii < aMatrixCheck.n_rows(); ii++ )
    {
        for( uint jj = 0; jj < aMatrixCheck.n_cols(); jj++ )
        {
            // get absolute difference
            tAbsolute = std::abs( aMatrixCheck( ii, jj ) - aMatrixRef( ii, jj ) );

            // get relative difference
            tRelative = std::abs( ( aMatrixRef( ii, jj ) - aMatrixCheck( ii, jj ) ) / aMatrixRef( ii, jj ) );

            // update check value
            tCheckBool = tCheckBool && ( ( tAbsolute < aEpsilon ) || ( tRelative < aEpsilon ) );

            // for debug
            if( ( ( tAbsolute < aEpsilon ) || ( tRelative < aEpsilon ) ) == false )
            {
                std::cout<<"ii "<<ii<<" - jj "<<jj<<"\n"<<std::flush;
                std::cout<<"aMatrixCheck( ii, jj ) "<<aMatrixCheck( ii, jj )<<"\n"<<std::flush;
                std::cout<<"aMatrixRef( ii, jj ) "<<aMatrixRef( ii, jj )<<"\n"<<std::flush;
                std::cout<<"Absolute difference "<<tAbsolute<<"\n"<<std::flush;
                std::cout<<"Relative difference "<<tRelative<<"\n"<<std::flush;
            }
        }
    }

    // return bool
    return tCheckBool;
}

using namespace moris;
using namespace fem;

TEST_CASE( "CM_Fluid", "[CM_Fluid]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-3;

    // define a perturbation relative size
    real tPerturbation = 1E-6;

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
        Matrix< DDRMat > tNormal;

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

                // fill normal
                tNormal = {{0.87},{-0.23}};

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

                // fill normal
                tNormal = {{0.87},{-0.23},{0.12}};

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
            arma::Mat< double > tMasterMatrixVis;

            // switch on interpolation order
            switch( iInterpOrder )
            {
                case 1 :
                {
                    // set interpolation type
                    tInterpolationOrder = mtk::Interpolation_Order::LINEAR;

                    // create random coefficients for master FI
                    tMasterMatrixVel.randu( tNumCoeffs( 0 ), iSpaceDim );
                    tMasterMatrixP.randu( tNumCoeffs( 0 ), 1 );
                    tMasterMatrixVis.randu( tNumCoeffs( 0 ), 1 );
                    break;
                }
                case 2 :
                {
                    // set interpolation type
                    tInterpolationOrder = mtk::Interpolation_Order::QUADRATIC;

                    // create random coefficients for master FI
                    tMasterMatrixVel.randu( tNumCoeffs( 1 ), iSpaceDim );
                    tMasterMatrixP.randu( tNumCoeffs( 1 ), 1 );
                    tMasterMatrixVis.randu( tNumCoeffs( 1 ), 1 );
                    break;
                }
                case 3 :
                {
                    // set interpolation type
                    tInterpolationOrder = mtk::Interpolation_Order::CUBIC;

                    // create random coefficients for master FI
                    tMasterMatrixVel.randu( tNumCoeffs( 2 ), iSpaceDim );
                    tMasterMatrixP.randu( tNumCoeffs( 2 ), 1 );
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
            Matrix< DDRMat > tMasterDOFHatVis;
            tMasterDOFHatVis.matrix_data() = 10.0 * tMasterMatrixVis;

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tMasterFIs( 3 );

            // create the field interpolator velocity
            tMasterFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelDofTypes );
            tMasterFIs( 0 )->set_coeff( tMasterDOFHatVel );
            tMasterFIs( 0 )->set_space_time( tParamPoint );

            // create the field interpolator pressure
            tMasterFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tPDofTypes );
            tMasterFIs( 1 )->set_coeff( tMasterDOFHatP );
            tMasterFIs( 1 )->set_space_time( tParamPoint );

            // create the field interpolator viscosity
            tMasterFIs( 2 ) = new Field_Interpolator( 1, tFIRule, &tGI, tVisDofTypes );
            tMasterFIs( 2 )->set_coeff( tMasterDOFHatVis );
            tMasterFIs( 2 )->set_space_time( tParamPoint );

            // create the properties
            std::shared_ptr< fem::Property > tPropViscosity = std::make_shared< fem::Property >();
            tPropViscosity->set_parameters( { {{ 1.0 }} } );
            //tPropViscosity->set_dof_type_list( { tPDofTypes } );
            tPropViscosity->set_val_function( tConstValFunc_CMFluid );
            //tPropViscosity->set_dof_derivative_functions( { tPFIDerFunction_NSVBULK } );

            std::shared_ptr< fem::Property > tPropDensity = std::make_shared< fem::Property >();
            tPropDensity->set_parameters( { {{ 2.0 }} } );
            //tPropDensity->set_dof_type_list( { tPDofTypes } );
            tPropDensity->set_val_function( tConstValFunc_CMFluid );
            //tPropDensity->set_dof_derivative_functions( { tPFIDerFunction_NSVBULK } );

            // define constitutive models
            fem::CM_Factory tCMFactory;

            std::shared_ptr< fem::Constitutive_Model > tCMMasterFluid =
                    tCMFactory.create_CM( fem::Constitutive_Type::FLUID_INCOMPRESSIBLE );
            tCMMasterFluid->set_dof_type_list( { tVelDofTypes, tPDofTypes } );
            tCMMasterFluid->set_property( tPropViscosity, "Viscosity" );
            tCMMasterFluid->set_property( tPropDensity, "Density" );
            tCMMasterFluid->set_space_dim( iSpaceDim );

            std::shared_ptr< fem::Constitutive_Model > tCMMasterTurbulence =
                    tCMFactory.create_CM( fem::Constitutive_Type::FLUID_TURBULENCE );
            tCMMasterTurbulence->set_dof_type_list( { tVelDofTypes, tVisDofTypes } );
            tCMMasterTurbulence->set_property( tPropViscosity, "Viscosity" );
            tCMMasterTurbulence->set_property( tPropDensity, "Density" );
            tCMMasterTurbulence->set_space_dim( iSpaceDim );

            // set a fem set pointer
            MSI::Equation_Set * tSet = new fem::Set();
            tCMMasterFluid->set_set_pointer( static_cast< fem::Set* >( tSet ) );
            tCMMasterTurbulence->set_set_pointer( static_cast< fem::Set* >( tSet ) );

            // set size for the set EqnObjDofTypeList
            tCMMasterFluid->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

            // set size and populate the set dof type map
            tCMMasterFluid->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
            tCMMasterFluid->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )        = 0;
            tCMMasterFluid->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )         = 1;
            tCMMasterFluid->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 2;

            // set size and populate the set master dof type map
            tCMMasterFluid->mSet->mMasterDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
            tCMMasterFluid->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )        = 0;
            tCMMasterFluid->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )         = 1;
            tCMMasterFluid->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 2;

            // build global dof type list
            tCMMasterFluid->get_global_dof_type_list();
            tCMMasterTurbulence->get_global_dof_type_list();

            // create a field interpolator manager
            moris::Cell< moris::Cell< enum MSI::Dof_Type > > tDummyDof;
            moris::Cell< moris::Cell< enum PDV_Type > > tDummyDv;
            Field_Interpolator_Manager tFIManager( tDummyDof, tDummyDv, tSet );

            // populate the field interpolator manager
            tFIManager.mFI = tMasterFIs;
            tFIManager.mIPGeometryInterpolator = &tGI;
            tFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tCMMasterFluid->mSet->mMasterFIManager = &tFIManager;

            // set CM field interpolator manager
            tCMMasterFluid->set_field_interpolator_manager( &tFIManager );
            tCMMasterTurbulence->set_field_interpolator_manager( &tFIManager );

//            // print for debug
//            std::cout<<"Case: Geometry "<<iSpaceDim<<" Order "<<iInterpOrder<<std::endl;

            // populate the requested master dof type for CM
            moris::Cell< moris::Cell< MSI::Dof_Type > > tRequestedMasterGlobalDofTypes =
                    tCMMasterFluid->get_global_dof_type_list();

            // loop over requested dof type
            for( uint iRequestedDof = 0; iRequestedDof < tRequestedMasterGlobalDofTypes.size(); iRequestedDof++ )
            {
                // derivative dof type
                Cell< MSI::Dof_Type > tDofDerivative = tRequestedMasterGlobalDofTypes( iRequestedDof );

                // flux
                //------------------------------------------------------------------------------
                // evaluate dfluxdu
                Matrix< DDRMat > tdfluxdu = tCMMasterFluid->dFluxdDOF( tDofDerivative );

                // evaluate dfluxdu by FD
                Matrix< DDRMat > tdfluxduFD;
                tCMMasterFluid->eval_dFluxdDOF_FD( tDofDerivative, tdfluxduFD, tPerturbation );

                // check that analytical and FD match
                bool tCheckFluxFluid = tCheckFunc( tdfluxdu, tdfluxduFD, tEpsilon );
                REQUIRE( tCheckFluxFluid );

                // strain
                //------------------------------------------------------------------------------
                // evaluate dstraindu
                Matrix< DDRMat > tdstraindu = tCMMasterFluid->dStraindDOF( tDofDerivative );

                // evaluate dstraindu by FD
                Matrix< DDRMat > tdstrainduFD;
                tCMMasterFluid->eval_dStraindDOF_FD( tDofDerivative, tdstrainduFD, tPerturbation );

                // check that analytical and FD match
                bool tCheckStrainFluid = tCheckFunc( tdstraindu, tdstrainduFD, tEpsilon );
                REQUIRE( tCheckStrainFluid );

                // traction
                //------------------------------------------------------------------------------
                // evaluate dtractiondu
                Matrix< DDRMat > tdtractiondu = tCMMasterFluid->dTractiondDOF( tDofDerivative, tNormal );

                // evaluate dtractiondu by FD
                Matrix< DDRMat > tdtractionduFD;
                tCMMasterFluid->eval_dtractiondu_FD( tDofDerivative, tdtractionduFD, tPerturbation, tNormal );

                // check that analytical and FD match
                bool tCheckTractionFluid = tCheckFunc( tdtractiondu, tdtractionduFD, tEpsilon );
                REQUIRE( tCheckTractionFluid );

                // div flux
                //------------------------------------------------------------------------------
                // evaluate ddivfluxdu
                Matrix< DDRMat > tddivfluxdu = tCMMasterFluid->ddivfluxdu( tDofDerivative );

                // evaluate ddivfluxdu by FD
                Matrix< DDRMat > tddivfluxduFD;
                tCMMasterFluid->eval_ddivfluxdu_FD( tDofDerivative, tddivfluxduFD, tPerturbation );

                // check that analytical and FD match
                bool tCheckDivFluxFluid = tCheckFunc( tddivfluxdu, tddivfluxduFD, tEpsilon );
                REQUIRE( tCheckDivFluxFluid );

                // div strain
                //------------------------------------------------------------------------------
                // evaluate ddivstraindu
                Matrix< DDRMat > tddivstraindu = tCMMasterFluid->ddivstraindu( tDofDerivative );

                // evaluate ddivstraindu by FD
                Matrix< DDRMat > tddivstrainduFD;
                tCMMasterFluid->eval_ddivstraindu_FD( tDofDerivative, tddivstrainduFD, tPerturbation );

                // check that analytical and FD match
                bool tCheckDivStrainFluid = tCheckFunc( tddivstraindu, tddivstrainduFD, tEpsilon );
                REQUIRE( tCheckDivStrainFluid );
            }

            // populate the requested master dof type for CM
            tRequestedMasterGlobalDofTypes =
                    tCMMasterTurbulence->get_global_dof_type_list();

            // loop over requested dof type
            for( uint iRequestedDof = 0; iRequestedDof < tRequestedMasterGlobalDofTypes.size(); iRequestedDof++ )
            {
                // derivative dof type
                Cell< MSI::Dof_Type > tDofDerivative = tRequestedMasterGlobalDofTypes( iRequestedDof );

                // flux
                //------------------------------------------------------------------------------
                // evaluate dfluxdu
                Matrix< DDRMat > tdfluxdu = tCMMasterTurbulence->dFluxdDOF( tDofDerivative );

                // evaluate dfluxdu by FD
                Matrix< DDRMat > tdfluxduFD;
                tCMMasterTurbulence->eval_dFluxdDOF_FD( tDofDerivative, tdfluxduFD, tPerturbation );

                // check that analytical and FD match
                bool tCheckFluxTurbulence = tCheckFunc( tdfluxdu, tdfluxduFD, tEpsilon );
                REQUIRE( tCheckFluxTurbulence );

                // strain
                //------------------------------------------------------------------------------
                // evaluate dstraindu
                Matrix< DDRMat > tdstraindu = tCMMasterTurbulence->dStraindDOF( tDofDerivative );

                // evaluate dstraindu by FD
                Matrix< DDRMat > tdstrainduFD;
                tCMMasterTurbulence->eval_dStraindDOF_FD( tDofDerivative, tdstrainduFD, tPerturbation );

                // check that analytical and FD match
                bool tCheckStrainTurbulence = tCheckFunc( tdstraindu, tdstrainduFD, tEpsilon );
                REQUIRE( tCheckStrainTurbulence );

                // traction
                //------------------------------------------------------------------------------
                // evaluate dtractiondu
                Matrix< DDRMat > tdtractiondu = tCMMasterTurbulence->dTractiondDOF( tDofDerivative, tNormal );

                // evaluate dtractiondu by FD
                Matrix< DDRMat > tdtractionduFD;
                tCMMasterTurbulence->eval_dtractiondu_FD( tDofDerivative, tdtractionduFD, tPerturbation, tNormal );

                // check that analytical and FD match
                bool tCheckTractionTurbulence = tCheckFunc( tdtractiondu, tdtractionduFD, tEpsilon );
                REQUIRE( tCheckTractionTurbulence );

                // div flux
                //------------------------------------------------------------------------------
                // evaluate ddivfluxdu
                Matrix< DDRMat > tddivfluxdu = tCMMasterTurbulence->ddivfluxdu( tDofDerivative );

                // evaluate ddivfluxdu by FD
                Matrix< DDRMat > tddivfluxduFD;
                tCMMasterTurbulence->eval_ddivfluxdu_FD( tDofDerivative, tddivfluxduFD, tPerturbation );

                // check that analytical and FD match
                bool tCheckDivFluxTurbulence = tCheckFunc( tddivfluxdu, tddivfluxduFD, tEpsilon );
                REQUIRE( tCheckDivFluxTurbulence );

                // div strain
                //------------------------------------------------------------------------------
                // evaluate ddivstraindu
                Matrix< DDRMat > tddivstraindu = tCMMasterTurbulence->ddivstraindu( tDofDerivative );

                // evaluate ddivstraindu by FD
                Matrix< DDRMat > tddivstrainduFD;
                tCMMasterTurbulence->eval_ddivstraindu_FD( tDofDerivative, tddivstrainduFD, tPerturbation );

                // check that analytical and FD match
                bool tCheckDivStrainTurbulence = tCheckFunc( tddivstraindu, tddivstrainduFD, tEpsilon );
                REQUIRE( tCheckDivStrainTurbulence );
            }
            // clean up
            tMasterFIs.clear();
        }
    }
}/*END_TEST_CASE*/

