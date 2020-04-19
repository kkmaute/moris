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

void tFIConstValFunction_UTPressureGhost
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 );
}

void tFIValFunction_UTPressureGhost
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 ) * sum( aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::VX )->val() );
}

void tFIDerFunction_UTPressureGhost
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    moris::Matrix< moris::DDRMat > tTemp = trans( aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::VX )->N() );
    moris::Matrix< moris::DDRMat > tReturn( tTemp.n_rows(), 1 );
    for( uint i = 0; i < tTemp.n_rows(); i++ )
    {
        tReturn( i ) = sum( tTemp.get_row( i ) );
    }
    aPropMatrix = aParameters( 0 )( 0, 0 ) * trans( tReturn );
}

using namespace moris;
using namespace fem;

// This UT tests the pressure ghost IWG for incompressible NS
// for QUAD, HEX geometry type
// for LINEAR, QUADRATIC and CUBIC interpolation order
TEST_CASE( "IWG_Incompressible_NS_Pressure_Ghost", "[moris],[fem],[IWG_Incompressible_NS_Pressure_Ghost]" )
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
        Cell< MSI::Dof_Type > tVelocityDofTypes;
        Cell< MSI::Dof_Type > tPressureDofTypes;

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

               // set dof type list
               tVelocityDofTypes = { MSI::Dof_Type::VX, MSI::Dof_Type::VY };
               tPressureDofTypes = { MSI::Dof_Type::P };

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

                // set dof type list
                tVelocityDofTypes = { MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::VZ };
                tPressureDofTypes = { MSI::Dof_Type::P };
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
            arma::Mat< double > tMasterMatrix;
            arma::Mat< double > tMasterPressureMatrix;
            arma::Mat< double > tSlaveMatrix;
            arma::Mat< double > tSlavePressureMatrix;

            // get number of dof
            int tNumDof = 0;
            int tNumPDof = 0;

            // switch on interpolation order
            switch( iInterpOrder )
            {
                case ( 1 ):
                {
                    // set interpolation type
                    tInterpolationOrder = mtk::Interpolation_Order::LINEAR;

                    // get number of dof
                    tNumDof  = tNumCoeffs( 0 ) * iSpaceDim;
                    tNumPDof = tNumCoeffs( 0 );

                    // create random coefficients for master FIs
                    tMasterMatrix.randu( tNumCoeffs( 0 ), iSpaceDim );
                    tMasterPressureMatrix.randu( tNumCoeffs( 0 ), 1 );

                    // create random coefficients for slave FIs
                    tSlaveMatrix.randu( tNumCoeffs( 0 ), iSpaceDim );
                    tSlavePressureMatrix.randu( tNumCoeffs( 0 ), 1 );
                    break;
                }
                case ( 2 ):
                {
                    // set interpolation type
                    tInterpolationOrder = mtk::Interpolation_Order::QUADRATIC;

                    // get number of dof
                    tNumDof  = tNumCoeffs( 1 ) * iSpaceDim;
                    tNumPDof = tNumCoeffs( 1 );

                    // create random coefficients for master FI
                    tMasterMatrix.randu( tNumCoeffs( 1 ), iSpaceDim );
                    tMasterPressureMatrix.randu( tNumCoeffs( 1 ), 1 );

                    // create random coefficients for slave FI
                    tSlaveMatrix.randu( tNumCoeffs( 1 ), iSpaceDim );
                    tSlavePressureMatrix.randu( tNumCoeffs( 1 ), 1 );
                    break;
                }
                case ( 3 ):
                {
                    // set interpolation type
                    tInterpolationOrder = mtk::Interpolation_Order::CUBIC;

                    // get number of dof
                    tNumDof  = tNumCoeffs( 2 ) * iSpaceDim;
                    tNumPDof = tNumCoeffs( 2 );

                    // create random coefficients for master FI
                    tMasterMatrix.randu( tNumCoeffs( 2 ), iSpaceDim );
                    tMasterPressureMatrix.randu( tNumCoeffs( 2 ), 1 );

                    // create random coefficients for slave FI
                    tSlaveMatrix.randu( tNumCoeffs( 2 ), iSpaceDim );
                    tSlavePressureMatrix.randu( tNumCoeffs( 2 ), 1 );
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
            Matrix< DDRMat > tMasterPressureDOFHat;
            tMasterPressureDOFHat.matrix_data() = 10.0 * tMasterPressureMatrix;

            // fill random coefficients for slave FI
            Matrix< DDRMat > tSlaveDOFHat;
            tSlaveDOFHat.matrix_data() = 10.0 * tSlaveMatrix;
            Matrix< DDRMat > tSlavePressureDOFHat;
            tSlavePressureDOFHat.matrix_data() = 10.0 * tSlavePressureMatrix;

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tMasterFIs( 2 );

            // create the field interpolator
            tMasterFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelocityDofTypes );
            tMasterFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tPressureDofTypes );

            // set the coefficients uHat
            tMasterFIs( 0 )->set_coeff( tMasterDOFHat );
            tMasterFIs( 1 )->set_coeff( tMasterPressureDOFHat );

            //set the evaluation point xi, tau
            tMasterFIs( 0 )->set_space_time( tParamPoint );
            tMasterFIs( 1 )->set_space_time( tParamPoint );

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tSlaveFIs( 2 );

            // create the field interpolator
            tSlaveFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelocityDofTypes );
            tSlaveFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tPressureDofTypes );

            // set the coefficients uHat
            tSlaveFIs( 0 )->set_coeff( tSlaveDOFHat );
            tSlaveFIs( 1 )->set_coeff( tSlavePressureDOFHat );

            //set the evaluation point xi, tau
            tSlaveFIs( 0 )->set_space_time( tParamPoint );
            tSlaveFIs( 1 )->set_space_time( tParamPoint );

            // create the properties
            std::shared_ptr< fem::Property > tPropMasterViscosity = std::make_shared< fem::Property >();
            tPropMasterViscosity->set_parameters( { {{ 1.0 }} } );
            tPropMasterViscosity->set_val_function( tFIConstValFunction_UTPressureGhost );

            std::shared_ptr< fem::Property > tPropMasterDensity = std::make_shared< fem::Property >();
            tPropMasterDensity->set_parameters( { {{ 1.0 }} } );
            tPropMasterDensity->set_val_function( tFIConstValFunction_UTPressureGhost );

            // define stabilization parameters
            fem::SP_Factory tSPFactory;

            std::shared_ptr< fem::Stabilization_Parameter > tSPPressureGhost
            = tSPFactory.create_SP( fem::Stabilization_Type::PRESSURE_GHOST );
            tSPPressureGhost->set_dof_type_list( { tVelocityDofTypes }, mtk::Master_Slave::MASTER );
            tSPPressureGhost->set_parameters( { {{ 1.0 }}, {{ 1.0 }} });
            tSPPressureGhost->set_property( tPropMasterViscosity, "Viscosity", mtk::Master_Slave::MASTER );
            tSPPressureGhost->set_property( tPropMasterDensity, "Density", mtk::Master_Slave::MASTER );

            // define the IWGs
            fem::IWG_Factory tIWGFactory;

            std::shared_ptr< fem::IWG > tIWGPressure
            = tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_GHOST );
            tIWGPressure->set_residual_dof_type( tPressureDofTypes );
            tIWGPressure->set_dof_type_list( { tVelocityDofTypes, tPressureDofTypes }, mtk::Master_Slave::MASTER );
            tIWGPressure->set_dof_type_list( { tVelocityDofTypes, tPressureDofTypes }, mtk::Master_Slave::SLAVE );
            tIWGPressure->set_stabilization_parameter( tSPPressureGhost, "PressureGhost" );

            // create and set the fem set for the IWG
            MSI::Equation_Set * tSet = new fem::Set();

            // set pointer for IWG
            tIWGPressure->set_set_pointer( static_cast< fem::Set* >( tSet ) );

            // set size for the set EqnObjDofTypeList
            tIWGPressure->mSet->mUniqueDofTypeList.resize( 4, MSI::Dof_Type::END_ENUM );

            // set size and populate the set dof type map
            tIWGPressure->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
            tIWGPressure->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) ) = 0;
            tIWGPressure->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )  = 1;

            // set size and populate the set master and slave dof type map
            tIWGPressure->mSet->mMasterDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
            tIWGPressure->mSet->mSlaveDofTypeMap .set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
            tIWGPressure->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) ) = 0;
            tIWGPressure->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )  = 1;
            tIWGPressure->mSet->mSlaveDofTypeMap ( static_cast< int >( MSI::Dof_Type::VX ) ) = 0;
            tIWGPressure->mSet->mSlaveDofTypeMap ( static_cast< int >( MSI::Dof_Type::P ) )  = 1;

            // set size and fill the set residual assembly map
            tIWGPressure->mSet->mResDofAssemblyMap.resize( 4 );
            tIWGPressure->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDof-1 } };
            tIWGPressure->mSet->mResDofAssemblyMap( 1 ) = { { tNumDof, tNumPDof + tNumDof - 1 } };
            tIWGPressure->mSet->mResDofAssemblyMap( 2 ) = { { tNumPDof + tNumDof, ( 2 * tNumDof ) + tNumDof - 1 } };
            tIWGPressure->mSet->mResDofAssemblyMap( 3 ) = { { ( 2 * tNumDof ) + tNumPDof, ( 2 * ( tNumPDof + tNumDof ) ) - 1 } };

            // set size and fill the set jacobian assembly map
            tIWGPressure->mSet->mJacDofAssemblyMap.resize( 4 );
            tIWGPressure->mSet->mJacDofAssemblyMap( 0 ) = { { 0, tNumDof-1 },
                                              { tNumDof, tNumDof + tNumPDof - 1 },
                                              { tNumDof + tNumPDof, ( 2 * tNumDof ) + tNumPDof - 1 },
                                              { ( 2 * tNumDof ) + tNumPDof, ( 2 * ( tNumDof + tNumPDof ) ) - 1 } };
            tIWGPressure->mSet->mJacDofAssemblyMap( 1 ) = { { 0, tNumDof-1 },
                                              { tNumDof, tNumDof + tNumPDof - 1 },
                                              { tNumDof + tNumPDof, ( 2 * tNumDof ) + tNumPDof - 1 },
                                              { ( 2 * tNumDof ) + tNumPDof, ( 2 * ( tNumDof + tNumPDof ) ) - 1 } };
            tIWGPressure->mSet->mJacDofAssemblyMap( 2 ) = { { 0, tNumDof-1 },
                                              { tNumDof, tNumDof + tNumPDof - 1 },
                                              { tNumDof + tNumPDof, ( 2 * tNumDof ) + tNumPDof - 1 },
                                              { ( 2 * tNumDof ) + tNumPDof, ( 2 * ( tNumDof + tNumPDof ) ) - 1 } };
            tIWGPressure->mSet->mJacDofAssemblyMap( 3 ) = { { 0, tNumDof-1 },
                                              { tNumDof, tNumDof + tNumPDof - 1 },
                                              { tNumDof + tNumPDof, ( 2 * tNumDof ) + tNumPDof - 1 },
                                              { ( 2 * tNumDof ) + tNumPDof, ( 2 * ( tNumDof + tNumPDof ) ) - 1 } };

            // set IWG normal
            tIWGPressure->set_normal( tNormal );

            // build global property type list
            tIWGPressure->build_global_dof_and_dv_type_list();

            // populate the requested master dof type
            tIWGPressure->mRequestedMasterGlobalDofTypes = { tVelocityDofTypes, tPressureDofTypes };
            tIWGPressure->mRequestedSlaveGlobalDofTypes  = { tVelocityDofTypes, tPressureDofTypes };

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

            // set the interpolator manager to the set
            tIWGPressure->mSet->mMasterFIManager = &tMasterFIManager;
            tIWGPressure->mSet->mSlaveFIManager  = &tSlaveFIManager;

            // set IWG field interpolator manager
            tIWGPressure->set_field_interpolator_manager( &tMasterFIManager );
            tIWGPressure->set_field_interpolator_manager( &tSlaveFIManager, mtk::Master_Slave::SLAVE );

            // reset residual and jacobian
            //------------------------------------------------------------------------------
            tIWGPressure->mSet->mResidual.resize( 1 );
            tIWGPressure->mSet->mResidual( 0 ).set_size( 2 * ( tNumPDof + tNumDof ), 1 , 0.0 );
            tIWGPressure->mSet->mJacobian.set_size( 2 * ( tNumPDof + tNumDof ), 2 * ( tNumPDof + tNumDof ), 0.0 );

            // check evaluation of the residual
            //------------------------------------------------------------------------------
            // evaluate the residual
            tIWGPressure->compute_residual( 1.0 );

            // check evaluation of the jacobian  by FD
            //------------------------------------------------------------------------------
            // init the jacobian for IWG and FD evaluation
            Matrix< DDRMat > tJacobians;
            Matrix< DDRMat > tJacobiansFD;

            // check jacobian by FD
            bool tCheckJacobian = tIWGPressure->check_jacobian_double( tPerturbation,
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
