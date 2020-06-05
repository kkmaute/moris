#include <string>
#include <catch.hpp>
#include "assert.hpp"

#define protected public
#define private   public
 //FEM//INT//src
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG.hpp"
#include "cl_FEM_Set.hpp"
#undef protected
#undef private

//MTK/src
#include "cl_MTK_Enums.hpp"
//FEM/INT/src
#include "cl_FEM_Enums.hpp"
#include "cl_FEM_IWG_Factory.hpp"
#include "cl_FEM_CM_Factory.hpp"
#include "cl_FEM_SP_Factory.hpp"
//LINALG
#include "op_equal_equal.hpp"

void tFIConstValFunction_UTViscosityInterface
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 );
}

void tFIValFunction_UTViscosityInterface
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 ) * sum( aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::VISCOSITY )->val() );
}

void tFIDerFunction_UTViscosityInterface
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    moris::Matrix< moris::DDRMat > tTemp = trans( aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::VISCOSITY )->N() );
    moris::Matrix< moris::DDRMat > tReturn( tTemp.n_rows(), 1 );
    for( uint i = 0; i < tTemp.n_rows(); i++ )
    {
        tReturn( i ) = sum( tTemp.get_row( i ) );
    }
    aPropMatrix = aParameters( 0 )( 0, 0 ) * trans( tReturn );
}

using namespace moris;
using namespace fem;

// This UT tests the velocity interface IWG for incompressible NS
// for QUAD, HEX geometry type
// for LINEAR, QUADRATIC and CUBIC interpolation order
TEST_CASE( "IWG_Spalart_Allmaras_Turbulence_Interface", "[moris],[fem],[IWG_Spalart_Allmaras_Turbulence_Interface]" )
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

        // dof type list
        Cell< MSI::Dof_Type > tVisDofTypes = { MSI::Dof_Type::VISCOSITY };

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
               tNormal = {{ 1.0 }, { 0.0 }};

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
        for( uint iInterpOrder = 1; iInterpOrder < 4; iInterpOrder++ )
        {
            // field interpolators
            //------------------------------------------------------------------------------
            // create an interpolation order
            mtk::Interpolation_Order tInterpolationOrder;

            // create random coefficients for master FI
            arma::Mat< double > tMasterMatrixVis;
            arma::Mat< double > tSlaveMatrixVis;

            // get number of dof
            int tNumDofVis  = 0;

            // switch on interpolation order
            switch( iInterpOrder )
            {
                case ( 1 ):
                {
                    // set interpolation type
                    tInterpolationOrder = mtk::Interpolation_Order::LINEAR;

                    // get number of dof
                    tNumDofVis = tNumCoeffs( 0 );

                    // create random coefficients for master FI
                    tMasterMatrixVis.randu( tNumCoeffs( 0 ), 1 );

                    // create random coefficients for slave FI
                    tSlaveMatrixVis.randu( tNumCoeffs( 0 ), 1 );
                    break;
                }
                case ( 2 ):
                {
                    // set interpolation type
                    tInterpolationOrder = mtk::Interpolation_Order::QUADRATIC;

                    // get number of dof
                    tNumDofVis = tNumCoeffs( 1 );

                    // create random coefficients for master FI
                    tMasterMatrixVis.randu( tNumCoeffs( 1 ), 1 );

                    // create random coefficients for slave FI
                    tSlaveMatrixVis.randu( tNumCoeffs( 1 ), 1 );
                    break;
                }
                case ( 3 ):
                {
                    // set interpolation type
                    tInterpolationOrder = mtk::Interpolation_Order::CUBIC;

                    // get number of dof
                    tNumDofVis = tNumCoeffs( 2 );

                    // create random coefficients for master FI
                    tMasterMatrixVis.randu( tNumCoeffs( 2 ), 1 );

                    // create random coefficients for slave FI
                    tSlaveMatrixVis.randu( tNumCoeffs( 2 ), 1 );
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
            Matrix< DDRMat > tMasterDOFHatVis;
            tMasterDOFHatVis.matrix_data() = 10.0 * tMasterMatrixVis;

            // fill random coefficients for slave FI
            Matrix< DDRMat > tSlaveDOFHatVis;
            tSlaveDOFHatVis.matrix_data() = 10.0 * tSlaveMatrixVis;

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tMasterFIs( 1 );

            // create the field interpolator
            tMasterFIs( 0 ) = new Field_Interpolator( 1, tFIRule, &tGI, tVisDofTypes );

            // set the coefficients uHat
            tMasterFIs( 0 )->set_coeff( tMasterDOFHatVis );

            //set the evaluation point xi, tau
            tMasterFIs( 0 )->set_space_time( tParamPoint );

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tSlaveFIs( 1 );

            // create the field interpolator
            tSlaveFIs( 0 ) = new Field_Interpolator( 1, tFIRule, &tGI, tVisDofTypes );

            // set the coefficients uHat
            tSlaveFIs( 0 )->set_coeff( tSlaveDOFHatVis );

            //set the evaluation point xi, tau
            tSlaveFIs( 0 )->set_space_time( tParamPoint );

            // create the master properties
            std::shared_ptr< fem::Property > tPropMasterViscosity = std::make_shared< fem::Property >();
            tPropMasterViscosity->set_parameters( { {{ 1.0 }} } );
            tPropMasterViscosity->set_val_function( tFIConstValFunction_UTViscosityInterface );

            // create the slave properties
            std::shared_ptr< fem::Property > tPropSlaveViscosity = std::make_shared< fem::Property >();
            tPropSlaveViscosity->set_parameters( { {{ 1.0 }} } );
            tPropSlaveViscosity->set_val_function( tFIConstValFunction_UTViscosityInterface );

            // define stabilization parameters
            fem::SP_Factory tSPFactory;

            std::shared_ptr< fem::Stabilization_Parameter > tSPNitsche =
                    tSPFactory.create_SP( fem::Stabilization_Type::NITSCHE_INTERFACE );
            tSPNitsche->set_parameters( { {{ 1.0 }} } );
            tSPNitsche->set_property( tPropMasterViscosity, "Material", mtk::Master_Slave::MASTER );
            tSPNitsche->set_property( tPropSlaveViscosity, "Material", mtk::Master_Slave::SLAVE );

            std::shared_ptr< fem::Stabilization_Parameter > tSPMasterWeight =
                    tSPFactory.create_SP( fem::Stabilization_Type::MASTER_WEIGHT_INTERFACE );
            tSPMasterWeight->set_property( tPropMasterViscosity, "Material", mtk::Master_Slave::MASTER );
            tSPMasterWeight->set_property( tPropSlaveViscosity, "Material", mtk::Master_Slave::SLAVE );

            std::shared_ptr< fem::Stabilization_Parameter > tSPSlaveWeight =
                    tSPFactory.create_SP( fem::Stabilization_Type::SLAVE_WEIGHT_INTERFACE );
            tSPSlaveWeight->set_property( tPropMasterViscosity, "Material", mtk::Master_Slave::MASTER );
            tSPSlaveWeight->set_property( tPropSlaveViscosity, "Material", mtk::Master_Slave::SLAVE );

            // define the IWGs
            fem::IWG_Factory tIWGFactory;

            std::shared_ptr< fem::IWG > tIWG =
                    tIWGFactory.create_IWG( fem::IWG_Type::SPALART_ALLMARAS_TURBULENCE_INTERFACE );
            tIWG->set_residual_dof_type( tVisDofTypes );
            tIWG->set_dof_type_list( { tVisDofTypes }, mtk::Master_Slave::MASTER );
            tIWG->set_dof_type_list( { tVisDofTypes }, mtk::Master_Slave::SLAVE );
            tIWG->set_property( tPropMasterViscosity, "Viscosity", mtk::Master_Slave::MASTER );
            tIWG->set_property( tPropSlaveViscosity, "Viscosity", mtk::Master_Slave::SLAVE );
            tIWG->set_stabilization_parameter( tSPNitsche, "NitscheInterface" );
            tIWG->set_stabilization_parameter( tSPMasterWeight, "MasterWeightInterface" );
            tIWG->set_stabilization_parameter( tSPSlaveWeight, "SlaveWeightInterface" );

            // create and set the fem set for the IWG
            MSI::Equation_Set * tSet = new fem::Set();

            // set pointer for IWG
            tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

            // set size for the set EqnObjDofTypeList
            tIWG->mSet->mUniqueDofTypeList.resize( 4, MSI::Dof_Type::END_ENUM );

            // set size and populate the set dof type map
            tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
            tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 0;

            // set size and populate the set master and slave dof type map
            tIWG->mSet->mMasterDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
            tIWG->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 0;

            tIWG->mSet->mSlaveDofTypeMap .set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
            tIWG->mSet->mSlaveDofTypeMap ( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 0;

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( 2 );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofVis-1 } };
            tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofVis, 2 * tNumDofVis - 1 } };

            // set size and fill the set jacobian assembly map
            Matrix< DDSMat > tJacAssembly = {
                    { 0, tNumDofVis-1 },
                    { tNumDofVis, 2 * tNumDofVis - 1 } };
            tIWG->mSet->mJacDofAssemblyMap.resize( 2 );
            tIWG->mSet->mJacDofAssemblyMap( 0 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 1 ) = tJacAssembly;

            // set IWG normal
            tIWG->set_normal( tNormal );

            // build global property type list
            tIWG->build_global_dof_and_dv_type_list();

            // populate the requested master dof type
            tIWG->mRequestedMasterGlobalDofTypes = { tVisDofTypes };
            tIWG->mRequestedSlaveGlobalDofTypes  = { tVisDofTypes };

            // create a field interpolator manager
            moris::Cell< moris::Cell< enum MSI::Dof_Type > > tDummy;
            Field_Interpolator_Manager tMasterFIManager( tDummy, tSet, mtk::Master_Slave::MASTER );
            Field_Interpolator_Manager tSlaveFIManager( tDummy, tSet, mtk::Master_Slave::SLAVE );

            // populate the field interpolator manager
            tMasterFIManager.mFI = tMasterFIs;
            tMasterFIManager.mIPGeometryInterpolator = &tGI;
            tMasterFIManager.mIGGeometryInterpolator = &tGI;
            tSlaveFIManager.mFI  = tSlaveFIs;
            tSlaveFIManager.mIPGeometryInterpolator = &tGI;
            tSlaveFIManager.mIGGeometryInterpolator = &tGI;

            // set IWG field interpolator manager
            tIWG->set_field_interpolator_manager( &tMasterFIManager );
            tIWG->set_field_interpolator_manager( &tSlaveFIManager, mtk::Master_Slave::SLAVE );

            // reset residual and jacobian
            //------------------------------------------------------------------------------
            tIWG->mSet->mResidual.resize( 1 );
            tIWG->mSet->mResidual( 0 ).set_size( 2 * tNumDofVis, 1 , 0.0 );
            tIWG->mSet->mJacobian.set_size( 2 * tNumDofVis, 2 * tNumDofVis, 0.0 );

            // check evaluation of the residual
            //------------------------------------------------------------------------------
            // evaluate the residual
            tIWG->compute_residual( 1.0 );

            // check evaluation of the jacobian  by FD
            //------------------------------------------------------------------------------
            // init the jacobian for IWG and FD evaluation
            Matrix< DDRMat > tJacobians;
            Matrix< DDRMat > tJacobiansFD;

            // check jacobian by FD
            bool tCheckJacobian = tIWG->check_jacobian_double(
                    tPerturbation,
                    tEpsilon,
                    1.0,
                    tJacobians,
                    tJacobiansFD );

//            // print for debug
//            print( tJacobians,"tJacobians");
//            print( tJacobiansFD,"tJacobiansFD");

            // require check is true
            REQUIRE( tCheckJacobian );

//            // print the treated case
//            std::cout<<"Case: Geometry "<<iSpaceDim<<" Order "<<iInterpOrder<<std::endl;

            // clean up
            tMasterFIs.clear();
            tSlaveFIs.clear();
        }
    }

}/* END_TEST_CASE */
