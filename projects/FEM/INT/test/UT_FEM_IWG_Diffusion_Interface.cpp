#include <string>
#include <catch.hpp>
#include "assert.hpp"

#define protected public
#define private   public
//FEM//INT//src
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Cluster.hpp"
#undef protected
#undef private
//MTK/src
#include "cl_MTK_Enums.hpp"
//FEM//INT//src
#include "cl_FEM_Enums.hpp"
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_FEM_Property.hpp"
#include "cl_FEM_IWG_Factory.hpp"
#include "cl_FEM_CM_Factory.hpp"
#include "cl_FEM_SP_Factory.hpp"
#include "FEM_Test_Proxy/cl_FEM_Inputs_for_Diffusion_UT.cpp"
//LINALG/src
#include "op_equal_equal.hpp"
#include "fn_norm.hpp"

using namespace moris;
using namespace fem;

void UT_IWG_Diff_Interface_Core ( enum fem::IWG_Type tIWGType )
{
    // define an epsilon environment
    real tEpsilon = 1E-6;

    // define a perturbation relative size
    real tPerturbation = 1E-6;

    // init geometry inputs
    //------------------------------------------------------------------------------
    // create geometry type
    mtk::Geometry_Type tGeometryType = mtk::Geometry_Type::UNDEFINED;

    // create space coeff xHat
    Matrix< DDRMat > tXHat;

    // create list of interpolation orders
    moris::Cell< mtk::Interpolation_Order > tInterpolationOrders = {
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Order::QUADRATIC,
            mtk::Interpolation_Order::CUBIC };

    // create list of integration orders
    moris::Cell< mtk::Integration_Order > tIntegrationOrders = {
            mtk::Integration_Order::QUAD_2x2,
            mtk::Integration_Order::HEX_2x2x2 };

    // create list with number of coeffs
    Matrix< DDRMat > tNumCoeffs = {{ 8, 18, 32 },{ 16, 54, 128 }};

    // dof type list
    moris::Cell< moris::Cell< MSI::Dof_Type > > tTempDofTypes  = { { MSI::Dof_Type::TEMP } };

    // init IWG
    //------------------------------------------------------------------------------
    // create the properties
    std::shared_ptr< fem::Property > tPropMasterConductivity = std::make_shared< fem::Property > ();
    tPropMasterConductivity->set_parameters( { {{ 1.0 }} } );
    tPropMasterConductivity->set_val_function( tConstValFunc_Diff );
    //            tPropMasterConductivity->set_dof_type_list( { tTempDofTypes } );
    //            tPropMasterConductivity->set_val_function( tTEMPFIValFunc_Diff );
    //            tPropMasterConductivity->set_dof_derivative_functions( { tTEMPFIDerFunc_Diff } );

    std::shared_ptr< fem::Property > tPropSlaveConductivity = std::make_shared< fem::Property > ();
    tPropSlaveConductivity->set_parameters( { {{ 5.0 }} } );
    tPropSlaveConductivity->set_val_function( tConstValFunc_Diff );
    //            tPropSlaveConductivity->set_dof_type_list( { tTempDofTypes } );
    //            tPropSlaveConductivity->set_val_function( tTEMPFIValFunc_Diff );
    //            tPropSlaveConductivity->set_dof_derivative_functions( { tTEMPFIDerFunc_Diff } );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMMasterDiffLinIso =
            tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
    tCMMasterDiffLinIso->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
    tCMMasterDiffLinIso->set_property( tPropMasterConductivity, "Conductivity" );
    tCMMasterDiffLinIso->set_space_dim( 3 );
    tCMMasterDiffLinIso->set_local_properties();

    std::shared_ptr< fem::Constitutive_Model > tCMSlaveDiffLinIso = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
    tCMSlaveDiffLinIso->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
    tCMSlaveDiffLinIso->set_property( tPropSlaveConductivity, "Conductivity" );
    tCMSlaveDiffLinIso->set_space_dim( 3 );
    tCMSlaveDiffLinIso->set_local_properties();

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPNitscheInterface =
            tSPFactory.create_SP( fem::Stabilization_Type::NITSCHE_INTERFACE );
    tSPNitscheInterface->set_parameters( { {{ 1.0 }} } );
    tSPNitscheInterface->set_property( tPropMasterConductivity, "Material", mtk::Master_Slave::MASTER );
    tSPNitscheInterface->set_property( tPropSlaveConductivity, "Material", mtk::Master_Slave::SLAVE );

    // create a dummy fem cluster and set it to SP
    fem::Cluster * tCluster = new fem::Cluster();
    tSPNitscheInterface->set_cluster( tCluster );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG =
            tIWGFactory.create_IWG( tIWGType );
    tIWG->set_residual_dof_type( tTempDofTypes );
    tIWG->set_dof_type_list( tTempDofTypes, mtk::Master_Slave::MASTER );
    tIWG->set_dof_type_list( tTempDofTypes, mtk::Master_Slave::SLAVE );
    tIWG->set_constitutive_model( tCMMasterDiffLinIso, "Diffusion", mtk::Master_Slave::MASTER );
    tIWG->set_constitutive_model( tCMSlaveDiffLinIso, "Diffusion", mtk::Master_Slave::SLAVE );
    tIWG->set_stabilization_parameter( tSPNitscheInterface, "NitscheInterface" );

    // init set info
    //------------------------------------------------------------------------------
    // set a fem set pointer
    MSI::Equation_Set * tSet = new fem::Set();
    static_cast<fem::Set*>(tSet)->set_set_type( fem::Element_Type::DOUBLE_SIDESET );
    tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tIWG->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 0;

    // set size and populate the set master dof type map
    tIWG->mSet->mMasterDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 0;

    // set size and populate the set slave dof type map
    tIWG->mSet->mSlaveDofTypeMap .set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mSlaveDofTypeMap ( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 0;

    // loop on the space dimension
    for( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // set geometry inputs
        //------------------------------------------------------------------------------

        // create and set normal
        Matrix< DDRMat > tNormal( iSpaceDim, 1, 0.5 );
        tNormal = tNormal / norm( tNormal );
        tIWG->set_normal( tNormal );

        // set space dimension to CM, SP
        tCMMasterDiffLinIso->set_space_dim( iSpaceDim );
        tCMSlaveDiffLinIso->set_space_dim( iSpaceDim );
        tSPNitscheInterface->set_space_dim( iSpaceDim );

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
        mtk::Interpolation_Rule tGIRule( tGeometryType,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::LINEAR,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::LINEAR );

        // create a space time geometry interpolator
        Geometry_Interpolator tGI = Geometry_Interpolator( tGIRule );

        // create time coeff tHat
        Matrix< DDRMat > tTHat = {{ 0.0 }, { 1.0 }};

        // set the coefficients xHat, tHat
        tGI.set_coeff( tXHat, tTHat );

        // loop on the interpolation order
        for( uint iInterpOrder = 1; iInterpOrder < 4; iInterpOrder++ )
        {
            // integration points
            //------------------------------------------------------------------------------
            // get an integration order
            mtk::Integration_Order tIntegrationOrder = tIntegrationOrders( iSpaceDim - 2 );

            // create an integration rule
            mtk::Integration_Rule tIntegrationRule(
                    tGeometryType,
                    mtk::Integration_Type::GAUSS,
                    tIntegrationOrder,
                    mtk::Geometry_Type::LINE,
                    mtk::Integration_Type::GAUSS,
                    mtk::Integration_Order::BAR_1 );

            // create an integrator
            mtk::Integrator tIntegrator( tIntegrationRule );

            // get integration points
            Matrix< DDRMat > tIntegPoints;
            tIntegrator.get_points( tIntegPoints );

            // field interpolators
            //------------------------------------------------------------------------------
            // create an interpolation order
            mtk::Interpolation_Order tInterpolationOrder = tInterpolationOrders( iInterpOrder - 1 );

            // number of dof for interpolation order
            uint tNumCoeff = tNumCoeffs( iSpaceDim - 2, iInterpOrder - 1 );

            // get number of dof per type
            int tNumDofTemp = tNumCoeff;

            //create a space time interpolation rule
            mtk::Interpolation_Rule tFIRule ( tGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    tInterpolationOrder,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // fill coefficients for master FI
            Matrix< DDRMat > tMasterDOFHatTemp;
            fill_that( tMasterDOFHatTemp, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tMasterFIs( tTempDofTypes.size() );

            // create the field interpolator temperature
            tMasterFIs( 0 ) = new Field_Interpolator( 1, tFIRule, &tGI, tTempDofTypes( 0 ) );
            tMasterFIs( 0 )->set_coeff( tMasterDOFHatTemp );

            // fill random coefficients for slave FI
            Matrix< DDRMat > tSlaveDOFHatTemp;
            fill_that( tSlaveDOFHatTemp, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tSlaveFIs( tTempDofTypes.size() );

            // create the field interpolator temperature
            tSlaveFIs( 0 ) = new Field_Interpolator( 1, tFIRule, &tGI, tTempDofTypes( 0 ) );
            tSlaveFIs( 0 )->set_coeff( tSlaveDOFHatTemp );

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( 2 * tTempDofTypes.size() );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofTemp - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofTemp, 2 * tNumDofTemp - 1 } };

            // set size and fill the set jacobian assembly map
            Matrix< DDSMat > tJacAssembly = { { 0, tNumDofTemp-1 }, { tNumDofTemp, 2 * tNumDofTemp - 1 } };
            tIWG->mSet->mJacDofAssemblyMap.resize( 2 * tTempDofTypes.size() );
            tIWG->mSet->mJacDofAssemblyMap( 0 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 1 ) = tJacAssembly;

            // set size and init the set residual and jacobian
            tIWG->mSet->mResidual.resize( 1 );
            tIWG->mSet->mResidual( 0 ).set_size( 2 * tNumDofTemp, 1, 0.0 );
            tIWG->mSet->mJacobian.set_size( 2 * tNumDofTemp, 2 * tNumDofTemp, 0.0 );

            // build global dof type list
            tIWG->get_global_dof_type_list();

            // populate the requested master dof type
            tIWG->mRequestedMasterGlobalDofTypes = tTempDofTypes;
            tIWG->mRequestedSlaveGlobalDofTypes  = tTempDofTypes;

            // create a field interpolator manager
            moris::Cell< moris::Cell< enum PDV_Type > > tDummyDv;
            moris::Cell< moris::Cell< enum mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager tMasterFIManager( tTempDofTypes, tDummyDv, tDummyField, tSet );
            Field_Interpolator_Manager tSlaveFIManager( tTempDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tMasterFIManager.mFI = tMasterFIs;
            tMasterFIManager.mIPGeometryInterpolator = &tGI;
            tMasterFIManager.mIGGeometryInterpolator = &tGI;
            tSlaveFIManager.mFI = tSlaveFIs;
            tSlaveFIManager.mIPGeometryInterpolator = &tGI;
            tSlaveFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tIWG->mSet->mMasterFIManager = &tMasterFIManager;
            tIWG->mSet->mSlaveFIManager = &tSlaveFIManager;

            // set IWG field interpolator manager
            tIWG->set_field_interpolator_manager( &tMasterFIManager, mtk::Master_Slave::MASTER );
            tIWG->set_field_interpolator_manager( &tSlaveFIManager, mtk::Master_Slave::SLAVE );

            // loop over inetgration points
            uint tNumGPs = tIntegPoints.n_cols();
            for( uint iGP = 0; iGP < tNumGPs; iGP ++ )
            {
                // reset IWG evaluation flags
                tIWG->reset_eval_flags();

                // create evaluation point xi, tau
                Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );

                // set integration point
                tIWG->mSet->mMasterFIManager->set_space_time( tParamPoint );
                tIWG->mSet->mSlaveFIManager->set_space_time( tParamPoint );

                // check evaluation of the residual for IWG
                //------------------------------------------------------------------------------
                // reset residual
                tIWG->mSet->mResidual( 0 ).fill( 0.0 );

                // compute residual
                tIWG->compute_residual( 1.0 );

                // check evaluation of the jacobian by FD
                //------------------------------------------------------------------------------
                // reset jacobian
                tIWG->mSet->mJacobian.fill( 0.0 );

                // init the jacobian for IWG and FD evaluation
                Matrix< DDRMat > tJacobian;
                Matrix< DDRMat > tJacobianFD;

                // check jacobian by FD
                bool tCheckJacobian = tIWG->check_jacobian(
                        tPerturbation,
                        tEpsilon,
                        1.0,
                        tJacobian,
                        tJacobianFD,
                        true );

                // check that jacobian is symmetric for symmetric version
                if ( tIWGType == fem::IWG_Type::SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE )
                {
                    real tRelError = norm( tJacobian - trans( tJacobian ) ) / norm ( tJacobian );
                    REQUIRE( tRelError < 1e-12 );
                }

                // print for debug
                if( !tCheckJacobian )
                {
                    std::cout<<"Case: Geometry "<<iSpaceDim<<" Order "<<iInterpOrder<<" iGP "<<iGP<<std::endl;
                }

                // require check is true
                REQUIRE( tCheckJacobian );
            }

            // clean up
            tMasterFIs.clear();
            tSlaveFIs.clear();
        }
    }
}

TEST_CASE( "IWG_Diff_Interface_Symmetric", "[moris],[fem],[IWG_Diff_Interface_Symmetric]" )
{
    UT_IWG_Diff_Interface_Core(  fem::IWG_Type::SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE );

}

TEST_CASE( "IWG_Diff_Interface_Unsymmetric", "[moris],[fem],[IWG_Diff_Interface_Unsymmetric]" )
{
    UT_IWG_Diff_Interface_Core( fem::IWG_Type::SPATIALDIFF_INTERFACE_UNSYMMETRIC_NITSCHE );
}
