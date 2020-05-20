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

void tFIConstValFunction_UTVelocityGhost
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 );
}

void tFIValFunction_UTVelocityGhost
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 ) * sum( aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::VX )->val() );
}

void tFIDerFunction_UTVelocityGhost
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

// This UT tests the velocity ghost IWG for incompressible NS
// for QUAD, HEX geometry type
// for LINEAR, QUADRATIC and CUBIC interpolation order
TEST_CASE( "IWG_Incompressible_NS_Velocity_Ghost", "[moris],[fem],[IWG_Incompressible_NS_Velocity_Ghost]" )
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
        Cell< MSI::Dof_Type > tDofTypes;

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
               tDofTypes = { MSI::Dof_Type::VX, MSI::Dof_Type::VY };

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
                tDofTypes = { MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::VZ };

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
            arma::Mat< double > tSlaveMatrix;

            // get number of dof
            int tNumDof = 0;

            // switch on interpolation order
            switch( iInterpOrder )
            {
                case ( 1 ):
                {
                    // set interpolation type
                    tInterpolationOrder = mtk::Interpolation_Order::LINEAR;

                    // get number of dof
                    tNumDof = tNumCoeffs( 0 ) * iSpaceDim;

                    // create random coefficients for master FI
                    tMasterMatrix.randu( tNumCoeffs( 0 ), iSpaceDim );

                    // create random coefficients for slave FI
                    tSlaveMatrix.randu( tNumCoeffs( 0 ), iSpaceDim );
                    break;
                }
                case ( 2 ):
                {
                    // set interpolation type
                    tInterpolationOrder = mtk::Interpolation_Order::QUADRATIC;

                    // get number of dof
                    tNumDof = tNumCoeffs( 1 ) * iSpaceDim;

                    // create random coefficients for master FI
                    tMasterMatrix.randu( tNumCoeffs( 1 ), iSpaceDim );

                    // create random coefficients for slave FI
                    tSlaveMatrix.randu( tNumCoeffs( 1 ), iSpaceDim );
                    break;
                }
                case ( 3 ):
                {
                    // set interpolation type
                    tInterpolationOrder = mtk::Interpolation_Order::CUBIC;

                    // get number of dof
                    tNumDof = tNumCoeffs( 2 ) * iSpaceDim;

                    // create random coefficients for master FI
                    tMasterMatrix.randu( tNumCoeffs( 2 ), iSpaceDim );

                    // create random coefficients for slave FI
                    tSlaveMatrix.randu( tNumCoeffs( 2 ), iSpaceDim );
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
            tMasterFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tDofTypes );

            // set the coefficients uHat
            tMasterFIs( 0 )->set_coeff( tMasterDOFHat );

            //set the evaluation point xi, tau
            tMasterFIs( 0 )->set_space_time( tParamPoint );

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tSlaveFIs( 1 );

            // create the field interpolator
            tSlaveFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tDofTypes );

            // set the coefficients uHat
            tSlaveFIs( 0 )->set_coeff( tSlaveDOFHat );

            //set the evaluation point xi, tau
            tSlaveFIs( 0 )->set_space_time( tParamPoint );

            // create the properties
            std::shared_ptr< fem::Property > tPropMasterViscosity = std::make_shared< fem::Property >();
            tPropMasterViscosity->set_parameters( { {{ 1.0 }} } );
            tPropMasterViscosity->set_val_function( tFIConstValFunction_UTVelocityGhost );

            std::shared_ptr< fem::Property > tPropMasterDensity = std::make_shared< fem::Property >();
            tPropMasterDensity->set_parameters( { {{ 1.0 }} } );
            tPropMasterDensity->set_val_function( tFIConstValFunction_UTVelocityGhost );

            // define stabilization parameters
            fem::SP_Factory tSPFactory;

            std::shared_ptr< fem::Stabilization_Parameter > tSPViscousGhost =
                    tSPFactory.create_SP( fem::Stabilization_Type::VISCOUS_GHOST );
            tSPViscousGhost->set_parameters( {{{ 1.0 }} });
            tSPViscousGhost->set_property( tPropMasterViscosity, "Viscosity", mtk::Master_Slave::MASTER );

            std::shared_ptr< fem::Stabilization_Parameter > tSPTimeGhost =
                    tSPFactory.create_SP( fem::Stabilization_Type::TIME_VELOCITY_GHOST );
            tSPTimeGhost->set_parameters( {{{ 1.0 }}, {{ 1.0 }} });
            tSPTimeGhost->set_property( tPropMasterDensity, "Density", mtk::Master_Slave::MASTER );

            std::shared_ptr< fem::Stabilization_Parameter > tSPConvectiveGhost =
                    tSPFactory.create_SP( fem::Stabilization_Type::CONVECTIVE_GHOST );
            tSPConvectiveGhost->set_dof_type_list( { tDofTypes }, mtk::Master_Slave::MASTER );
            tSPConvectiveGhost->set_parameters( {{{ 1.0 }} });
            tSPConvectiveGhost->set_property( tPropMasterDensity, "Density", mtk::Master_Slave::MASTER );

            // define the IWGs
            fem::IWG_Factory tIWGFactory;

            std::shared_ptr< fem::IWG > tIWGViscous =
                    tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_VISCOUS_VELOCITY_GHOST );
            tIWGViscous->set_residual_dof_type( tDofTypes );
            tIWGViscous->set_dof_type_list( { tDofTypes }, mtk::Master_Slave::MASTER );
            tIWGViscous->set_dof_type_list( { tDofTypes }, mtk::Master_Slave::SLAVE );
            tIWGViscous->set_stabilization_parameter( tSPViscousGhost, "ViscousGhost" );
            tIWGViscous->set_stabilization_parameter( tSPTimeGhost, "TimeVelocityGhost" );

            std::shared_ptr< fem::IWG > tIWGConvective =
                    tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_CONVECTIVE_VELOCITY_GHOST );
            tIWGConvective->set_residual_dof_type( tDofTypes );
            tIWGConvective->set_dof_type_list( { tDofTypes }, mtk::Master_Slave::MASTER );
            tIWGConvective->set_dof_type_list( { tDofTypes }, mtk::Master_Slave::SLAVE );
            tIWGConvective->set_stabilization_parameter( tSPConvectiveGhost, "ConvectiveGhost" );

            // create and set the fem set for the IWG
            MSI::Equation_Set * tSet = new fem::Set();

            // set pointer for IWG
            tIWGViscous->set_set_pointer( static_cast< fem::Set* >( tSet ) );
            tIWGConvective->set_set_pointer( static_cast< fem::Set* >( tSet ) );

            // set size for the set EqnObjDofTypeList
            tIWGViscous->mSet->mUniqueDofTypeList.resize( 4, MSI::Dof_Type::END_ENUM );

            // set size and populate the set dof type map
            tIWGViscous->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
            tIWGViscous->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) ) = 0;

            // set size and populate the set master and slave dof type map
            tIWGViscous->mSet->mMasterDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
            tIWGViscous->mSet->mSlaveDofTypeMap .set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
            tIWGViscous->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) ) = 0;
            tIWGViscous->mSet->mSlaveDofTypeMap ( static_cast< int >( MSI::Dof_Type::VX ) ) = 0;

            // set size and fill the set residual assembly map
            tIWGViscous->mSet->mResDofAssemblyMap.resize( 2 );
            tIWGViscous->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDof-1 } };
            tIWGViscous->mSet->mResDofAssemblyMap( 1 ) = { { tNumDof, ( 2 * tNumDof )-1 } };

            // set size and fill the set jacobian assembly map
            tIWGViscous->mSet->mJacDofAssemblyMap.resize( 2 );
            tIWGViscous->mSet->mJacDofAssemblyMap( 0 ) = { { 0, tNumDof-1 },{ tNumDof, ( 2 * tNumDof )-1 } };
            tIWGViscous->mSet->mJacDofAssemblyMap( 1 ) = { { 0, tNumDof-1 },{ tNumDof, ( 2 * tNumDof )-1 } };

            // set IWG normal
            tIWGViscous->set_normal( tNormal );
            tIWGConvective->set_normal( tNormal );

            // build global property type list
            tIWGViscous->build_global_dof_and_dv_type_list();
            tIWGConvective->build_global_dof_and_dv_type_list();

            // populate the requested master dof type
            tIWGViscous->mRequestedMasterGlobalDofTypes = { tDofTypes };
            tIWGViscous->mRequestedSlaveGlobalDofTypes  = { tDofTypes };
            tIWGConvective->mRequestedMasterGlobalDofTypes = { tDofTypes };
            tIWGConvective->mRequestedSlaveGlobalDofTypes  = { tDofTypes };

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
            tIWGViscous->set_field_interpolator_manager( &tMasterFIManager );
            tIWGViscous->set_field_interpolator_manager( &tSlaveFIManager, mtk::Master_Slave::SLAVE );
            tIWGConvective->set_field_interpolator_manager( &tMasterFIManager );
            tIWGConvective->set_field_interpolator_manager( &tSlaveFIManager, mtk::Master_Slave::SLAVE );

            // reset residual and jacobian
            //------------------------------------------------------------------------------
            tIWGViscous->mSet->mResidual.resize( 1 );
            tIWGViscous->mSet->mResidual( 0 ).set_size( 2 * tNumDof, 1 , 0.0 );
            tIWGViscous->mSet->mJacobian.set_size( 2 * tNumDof, 2 * tNumDof, 0.0 );

            // check evaluation of the residual
            //------------------------------------------------------------------------------
            // evaluate the residual
            tIWGViscous->compute_residual( 1.0 );

            // check evaluation of the jacobian  by FD
            //------------------------------------------------------------------------------
            // init the jacobian for IWG and FD evaluation
            Matrix< DDRMat > tJacobians;
            Matrix< DDRMat > tJacobiansFD;

            // check jacobian by FD
            bool tCheckJacobian = tIWGViscous->check_jacobian_double( tPerturbation,
                                                                      tEpsilon,
                                                                      1.0,
                                                                      tJacobians,
                                                                      tJacobiansFD );

//            // print for debug
//            print( tJacobians,"tJacobians");
//            print( tJacobiansFD,"tJacobiansFD");

            // require check is true
            REQUIRE( tCheckJacobian );

            // reset residual and jacobian
            //------------------------------------------------------------------------------
            tIWGConvective->mSet->mResidual.resize( 1 );
            tIWGConvective->mSet->mResidual( 0 ).set_size( 2 * tNumDof, 1 , 0.0 );
            tIWGConvective->mSet->mJacobian.set_size( 2 * tNumDof, 2 * tNumDof, 0.0 );

            // check evaluation of the residual
            //------------------------------------------------------------------------------
            // evaluate the residual
            tIWGConvective->compute_residual( 1.0 );

            // init the jacobian for IWG and FD evaluation
            Matrix< DDRMat > tConvectiveJacobians;
            Matrix< DDRMat > tConvectiveJacobiansFD;

            // check jacobian by FD
            bool tCheckConvectiveJacobian
            = tIWGConvective->check_jacobian_double( tPerturbation,
                                                     tEpsilon,
                                                     1.0,
                                                     tConvectiveJacobians,
                                                     tConvectiveJacobiansFD );

//            // print for debug
//            print( tConvectiveJacobians,"tConvectiveJacobians");
//            print( tConvectiveJacobiansFD,"tConvectiveJacobiansFD");

            // require check is true
            REQUIRE( tCheckConvectiveJacobian );

//            // print the treated case
//            std::cout<<"Case: Geometry "<<iSpaceDim<<" Order "<<iInterpOrder<<std::endl;


            // clean up
            tMasterFIs.clear();
            tSlaveFIs.clear();
        }
    }

}/* END_TEST_CASE */
