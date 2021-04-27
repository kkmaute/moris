#include <string>
#include <catch.hpp>
#include <memory>
#include "assert.hpp"

#define protected public
#define private   public
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Cluster.hpp"
#undef protected
#undef private
//MTK/src
#include "cl_MTK_Enums.hpp"
//LINALG/src
#include "op_equal_equal.hpp"
//FEM/INT/src
#include "cl_FEM_Enums.hpp"
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_FEM_Property.hpp"
#include "cl_FEM_Stabilization_Parameter.hpp"
#include "cl_FEM_CM_Factory.hpp"
#include "cl_FEM_IWG_Factory.hpp"
#include "cl_FEM_SP_Factory.hpp"

void tConstValFunction_UTIWGGGLSDIFFBULK
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
        moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
        moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 );
}

void tGeoValFunction_UTIWGGGLSDIFFBULK
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
        moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
        moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 ) * aFIManager->get_IP_geometry_interpolator()->valx()( 0 );
}

void tFIValFunction_UTIWGGGLSDIFFBULK
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
        moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
        moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 ) + 0.1 * aParameters( 0 ) * aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::TEMP )->val();
}

void tFIDerFunction_UTIWGGGLSDIFFBULK
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
        moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
        moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = 0.1 * aParameters( 0 ) * aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::TEMP )->N();
}

void tFIDer0Function_UTIWGGGLSDIFFBULK
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
        moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
        moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = 0.0 * aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::TEMP )->N();
}

namespace moris
{
    namespace fem
    {

        moris::Cell<bool> test_IWG_Diffusion_Phase_Change_GGLS(
                Matrix< DDRMat > aXHat,
                Matrix< DDRMat > aTHat,
                mtk::Interpolation_Rule aGIRule,
                mtk::Interpolation_Rule aFIRule,
                Matrix< DDRMat > aDOFHat,
                Matrix< DDRMat > aParamPoint,
                uint aNumDOFs,
                uint aSpatialDim = 2 )
        {
            // initialize cell of checks
            moris::Cell<bool> tChecks( 1, false );

            // define an epsilon environment
            real tEpsilonRel = 1.0E-6;

            // define a perturbation relative size
            real tPerturbation = 1.0E-6;

            // create the properties ------------------------------------------------------------------- //

            std::shared_ptr< fem::Property > tPropMasterConductivity = std::make_shared< fem::Property > ();
            tPropMasterConductivity->set_parameters( {{{ 1.1 }}, {{ 1.1 }}} );
            //    tPropMasterConductivity->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            //    tPropMasterConductivity->set_val_function( tFIValFunction_UTIWGDIFFBULK );
            //    tPropMasterConductivity->set_dof_derivative_functions( { tFIDerFunction_UTIWGGGLSDIFFBULK } );
            tPropMasterConductivity->set_val_function( tConstValFunction_UTIWGGGLSDIFFBULK );

            std::shared_ptr< fem::Property > tPropMasterDensity = std::make_shared< fem::Property > ();
            tPropMasterDensity->set_parameters( { {{ 1.2 }} } );
            tPropMasterDensity->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tPropMasterDensity->set_val_function( tFIValFunction_UTIWGGGLSDIFFBULK );
            tPropMasterDensity->set_dof_derivative_functions( { tFIDerFunction_UTIWGGGLSDIFFBULK } );
//            tPropMasterDensity->set_val_function( tConstValFunction_UTIWGGGLSDIFFBULK );
//            tPropMasterDensity->set_dof_derivative_functions( { tFIDer0Function_UTIWGGGLSDIFFBULK } );

            std::shared_ptr< fem::Property > tPropMasterHeatCapacity = std::make_shared< fem::Property > ();
            tPropMasterHeatCapacity->set_parameters( { {{ 0.3 }} } );
            tPropMasterHeatCapacity->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tPropMasterHeatCapacity->set_val_function( tFIValFunction_UTIWGGGLSDIFFBULK );
            tPropMasterHeatCapacity->set_dof_derivative_functions( { tFIDerFunction_UTIWGGGLSDIFFBULK } );
//            tPropMasterHeatCapacity->set_val_function( tConstValFunction_UTIWGGGLSDIFFBULK );
//            tPropMasterHeatCapacity->set_dof_derivative_functions( { tFIDer0Function_UTIWGGGLSDIFFBULK } );

            // latent heat
            std::shared_ptr< fem::Property > tPropMasterLatentHeat = std::make_shared< fem::Property >();
            tPropMasterLatentHeat->set_parameters( {{{ 100.0}}} );
            //tPropMasterLatentHeat->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tPropMasterLatentHeat->set_val_function( tConstValFunction_UTIWGGGLSDIFFBULK );

            // phase change temp
            std::shared_ptr< fem::Property > tPropMasterTmelt = std::make_shared< fem::Property >();
            tPropMasterTmelt->set_parameters( {{{ 5.2 }}} );
            //tPropMasterTupper->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tPropMasterTmelt->set_val_function( tConstValFunction_UTIWGGGLSDIFFBULK );

            // phase change constant
            std::shared_ptr< fem::Property > tPropMasterPCconst = std::make_shared< fem::Property >();
            tPropMasterPCconst->set_parameters( {{{ 2.2 }}} );
            //tPropMasterPCconst->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tPropMasterPCconst->set_val_function( tConstValFunction_UTIWGGGLSDIFFBULK );

            // phase state function type
            std::shared_ptr< fem::Property > tPropMasterPCfunction = std::make_shared< fem::Property >();
            tPropMasterPCfunction->set_parameters( {{{ 2 }}} );
            //tPropMasterPCfunction->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tPropMasterPCfunction->set_val_function( tConstValFunction_UTIWGGGLSDIFFBULK );

            // temperature load
            std::shared_ptr< fem::Property > tPropMasterBodyLoad = nullptr;

            // define constitutive model ---------------------------------------------------------------- //
            fem::CM_Factory tCMFactory;

            std::shared_ptr< fem::Constitutive_Model > tCMMasterDiffLinIsoPC = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO_PC );
            tCMMasterDiffLinIsoPC->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tCMMasterDiffLinIsoPC->set_property( tPropMasterConductivity, "Conductivity" );
            tCMMasterDiffLinIsoPC->set_property( tPropMasterDensity     , "Density" );
            tCMMasterDiffLinIsoPC->set_property( tPropMasterHeatCapacity, "HeatCapacity" );
            tCMMasterDiffLinIsoPC->set_property( tPropMasterLatentHeat  , "LatentHeat" );
            tCMMasterDiffLinIsoPC->set_property( tPropMasterTmelt       , "PCTemp" );
            tCMMasterDiffLinIsoPC->set_property( tPropMasterPCfunction  , "PhaseStateFunction" );
            tCMMasterDiffLinIsoPC->set_property( tPropMasterPCconst     , "PhaseChangeConst" );
            tCMMasterDiffLinIsoPC->set_space_dim( 3 );
            tCMMasterDiffLinIsoPC->set_local_properties();

            // define stabilization parameter ----------------------------------------------------------- //
            fem::SP_Factory tSPFactory;

            std::shared_ptr< fem::Stabilization_Parameter > tSPGGLSParam = tSPFactory.create_SP( fem::Stabilization_Type::GGLS_DIFFUSION );
            tSPGGLSParam->set_dof_type_list( {{ MSI::Dof_Type::TEMP }}, mtk::Master_Slave::MASTER );
            tSPGGLSParam->set_parameters( { {{ 1.0 }} });
            tSPGGLSParam->set_property( tPropMasterConductivity, "Conductivity", mtk::Master_Slave::MASTER );
            tSPGGLSParam->set_property( tPropMasterDensity, "Density", mtk::Master_Slave::MASTER );
            tSPGGLSParam->set_property( tPropMasterHeatCapacity, "HeatCapacity", mtk::Master_Slave::MASTER );

            // create a dummy fem cluster and set it to SP
            fem::Cluster * tCluster = new fem::Cluster();
            tSPGGLSParam->set_cluster( tCluster );

            // define the IWGs
            fem::IWG_Factory tIWGFactory;

            std::shared_ptr< fem::IWG > tIWG = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
            tIWG->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWG->set_dof_type_list( {{ MSI::Dof_Type::TEMP }}, mtk::Master_Slave::MASTER );
            tIWG->set_constitutive_model( tCMMasterDiffLinIsoPC, "Diffusion", mtk::Master_Slave::MASTER );
            tIWG->set_stabilization_parameter( tSPGGLSParam, "GGLSParam");
            tIWG->set_property( tPropMasterBodyLoad, "Load", mtk::Master_Slave::MASTER );

            // space and time geometry interpolators
            //------------------------------------------------------------------------------

            // create a space time geometry interpolator
            Geometry_Interpolator tGI( aGIRule );

            // set the coefficients xHat, tHat
            tGI.set_coeff( aXHat, aTHat );

            // set the evaluation point
            tGI.set_space_time( aParamPoint );


            // field interpolators
            //------------------------------------------------------------------------------

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tFIs( 1 );

            // create the field interpolator
            tFIs( 0 ) = new Field_Interpolator( 1, aFIRule, &tGI, { MSI::Dof_Type::TEMP } );

            // set the coefficients uHat
            tFIs( 0 )->set_coeff( aDOFHat );

            //set the evaluation point xi, tau
            tFIs( 0 )->set_space_time( aParamPoint );

            // set a fem set pointer
            MSI::Equation_Set * tSet = new fem::Set();
            static_cast<fem::Set*>(tSet)->set_set_type( fem::Element_Type::BULK );
            tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

            // set size for the set EqnObjDofTypeList
            tIWG->mSet->mUniqueDofTypeList.resize( 4, MSI::Dof_Type::END_ENUM );

            // set size and populate the set dof type map
            tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
            tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 0;

            // set size and populate the set master dof type map
            tIWG->mSet->mMasterDofTypeMap.set_size( static_cast< int >(MSI::Dof_Type::END_ENUM) + 1, 1, -1 );
            tIWG->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 0;

            int aInt = (aNumDOFs-1);

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( 1 );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, aInt } };

            // set size and fill the set jacobian assembly map
            tIWG->mSet->mJacDofAssemblyMap.resize( 1 );
            tIWG->mSet->mJacDofAssemblyMap( 0 ) = { { 0, aInt } };

            // set size and init the set residual and jacobian
            tIWG->mSet->mResidual.resize( 1 );
            tIWG->mSet->mResidual( 0 ).set_size( aNumDOFs, 1, 0.0 );
            tIWG->mSet->mJacobian.set_size( aNumDOFs, aNumDOFs, 0.0 );

            // build global dof type list
            tIWG->get_global_dof_type_list();

            // populate the requested master dof type
            tIWG->mRequestedMasterGlobalDofTypes = {{ MSI::Dof_Type::TEMP }};

            // create a field interpolator manager
            moris::Cell< moris::Cell< enum MSI::Dof_Type > > tDummy;
            Field_Interpolator_Manager tFIManager( tDummy, tSet );

            // populate the field interpolator manager
            tFIManager.mFI = tFIs;
            tFIManager.mIPGeometryInterpolator = &tGI;
            tFIManager.mIGGeometryInterpolator = &tGI;

            // set IWG field interpolator manager
            tIWG->set_field_interpolator_manager( &tFIManager );

            // check evaluation of the residual for IWG
            //------------------------------------------------------------------------------
            // evaluate the residual
            tIWG->compute_residual( 1.0 );

            // check evaluation of the jacobian  by FD
            //------------------------------------------------------------------------------
            // init the jacobian for IWG and FD evaluation
            Matrix< DDRMat > tJacobian;
            Matrix< DDRMat > tJacobianFD;

            // check jacobian by FD
            bool tCheckJacobian = tIWG->check_jacobian( tPerturbation,
                    tEpsilonRel,
                    1.0,
                    tJacobian,
                    tJacobianFD );


            // require check is true
            //REQUIRE( tCheckJacobian );
            tChecks(0) = tCheckJacobian;

            // debug
            // moris::Matrix<DDRMat> test1 = tJacobianFD-tJacobian;
            // real tMax = test1.max();
            // print( tJacobian,   "tJacobian" );
            // print( tJacobianFD, "tJacobianFD" );
            //print( test1, "JacobianDifference" );
            // std::cout << "Maximum difference = " << tMax << " \n" << std::flush;

            return tChecks;

        } // end test function

        // ------------------------------------------------------------------------------------- //
        // ------------------------------------------------------------------------------------- //
        TEST_CASE( "IWG_GGLS_Diffusion_Phase_Change_HEX8", "[moris],[fem],[IWG_GGLS_Diffusion_Phase_Change_HEX8]" )
        {
            //create a quad4 space element
            Matrix< DDRMat > tXHat = {
                    { 0.0, 0.0, 0.0 },
                    { 1.0, 0.0, 0.0 },
                    { 1.0, 1.0, 0.0 },
                    { 0.0, 1.0, 0.0 },
                    { 0.0, 0.0, 1.0 },
                    { 1.0, 0.0, 1.0 },
                    { 1.0, 1.0, 1.0 },
                    { 0.0, 1.0, 1.0 }};

            //create a line time element
            Matrix< DDRMat > tTHat( 2, 1 );
            tTHat( 0 ) = 1.0e-3;
            tTHat( 1 ) = 1.1e-3;

            //create a space geometry interpolation rule
            mtk::Interpolation_Rule tGeomInterpRule( mtk::Geometry_Type::HEX,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // create an interpolation rule
            mtk::Interpolation_Rule tIPRule (
                    mtk::Geometry_Type::HEX,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // set coefficients for field interpolators
            Matrix< DDRMat > tUHat0 = {
                    {4.7},{4.8},{4.9},{5.3},{5.2},{5.3},{5.4},{4.5},
                    {4.9},{5.2},{4.9},{5.3},{5.7},{5.1},{4.8},{4.8}};
            uint tNumDOFs = 16;
            Matrix< DDRMat > tParametricPoint = {{ 0.35}, {-0.25}, { 0.75}, { 0.4 }};

            // run test
            moris::Cell<bool> tChecks = test_IWG_Diffusion_Phase_Change_GGLS(
                    tXHat,
                    tTHat,
                    tGeomInterpRule,
                    tIPRule,
                    tUHat0,
                    tParametricPoint,
                    tNumDOFs,
                    3);

            // checks
            bool tCheckJacobian = tChecks(0);
            REQUIRE( tCheckJacobian );

        } // end TEST_CASE

        // ------------------------------------------------------------------------------------- //
        // ------------------------------------------------------------------------------------- //
        TEST_CASE( "IWG_GGLS_Diffusion_Phase_Change_HEX27", "[moris],[fem],[IWG_GGLS_Diffusion_Phase_Change_HEX27]" )
        {
            //create a HEX27 space element
            Matrix< DDRMat > tXHat = {
                    { 0.0, 0.0, 0.0}, { 4.0, 0.0, 0.0}, { 4.0, 1.0, 0.0}, { 0.0, 1.0, 0.0},
                    { 0.0, 0.0, 3.0}, { 4.0, 0.0, 3.0}, { 4.0, 1.0, 3.0}, { 0.0, 1.0, 3.0},
                    { 2.0, 0.0, 0.0}, { 4.0, 0.5, 0.0}, { 2.0, 1.0, 0.0}, { 0.0, 0.5, 0.0},
                    { 0.0, 0.0, 1.5}, { 4.0, 0.0, 1.5}, { 4.0, 1.0, 1.5}, { 0.0, 1.0, 1.5},
                    { 2.0, 0.0, 3.0}, { 4.0, 0.5, 3.0}, { 2.0, 1.0, 3.0}, { 0.0, 0.5, 3.0},
                    { 2.0, 0.5, 1.5}, { 2.0, 0.5, 0.0}, { 2.0, 0.5, 3.0},
                    { 2.0, 0.0, 1.5}, { 4.0, 0.5, 1.5}, { 2.0, 1.0, 1.5}, { 0.0, 0.5, 1.5}};

            //create a line time element
            Matrix< DDRMat > tTHat( 3, 1 );
            tTHat( 0 ) = 1.00e-3;
            tTHat( 2 ) = 1.05e-3;
            tTHat( 1 ) = 1.10e-3;

            //create a space geometry interpolation rule
            mtk::Interpolation_Rule tGeomInterpRule( mtk::Geometry_Type::HEX,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::QUADRATIC,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::QUADRATIC );

            // create an interpolation rule
            mtk::Interpolation_Rule tIPRule (
                    mtk::Geometry_Type::HEX,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::QUADRATIC,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::QUADRATIC );

            // set coefficients for field interpolators
            Matrix< DDRMat > tDOFHat = {
                    {5.1},{5.2},{5.3},{5.4},{5.5},{4.6},{4.7},{5.4},{4.9},{5.1},{5.2},{5.3},{5.4},{5.5},{4.6},{4.7},{4.8},{5.1},{5.3},{5.2},{5.3},{5.4},{4.5},{4.6},{4.7},{4.8},{4.9},
                    {5.1},{5.2},{5.3},{5.3},{5.5},{4.8},{4.7},{4.8},{4.9},{5.3},{5.2},{5.3},{5.4},{4.5},{5.2},{4.7},{4.8},{4.9},{5.1},{5.4},{5.3},{4.6},{5.5},{4.9},{4.7},{4.8},{4.9},
                    {5.4},{5.2},{5.3},{5.4},{5.2},{5.2},{4.7},{5.1},{5.3},{5.1},{5.1},{5.1},{5.2},{5.3},{5.2},{4.7},{4.6},{4.5},{5.6},{5.2},{5.3},{5.4},{5.3},{4.5},{5.3},{5.2},{5.4}};

            uint tNumDOFs = 27 * 3;
            Matrix< DDRMat > tParametricPoint = {{ 0.35}, {-0.25}, { 0.75}, { 0.4 }};

            // run test
            moris::Cell<bool> tChecks = test_IWG_Diffusion_Phase_Change_GGLS(
                    tXHat,
                    tTHat,
                    tGeomInterpRule,
                    tIPRule,
                    tDOFHat,
                    tParametricPoint,
                    tNumDOFs,
                    3);

            // checks
            bool tCheckJacobian = tChecks(0);
            REQUIRE( tCheckJacobian );

        } // end TEST_CASE

    } // end namespace fem
} // end namespace moris
