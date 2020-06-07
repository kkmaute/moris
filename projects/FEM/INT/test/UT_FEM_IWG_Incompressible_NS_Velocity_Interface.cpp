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

void tFIConstValFunction_UTVelocityInterface
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 );
}

void tFIValFunction_UTVelocityInterface
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 ) * sum( aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::VX )->val() );
}

void tFIDerFunction_UTVelocityInterface
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

// This UT tests the velocity interface IWG for incompressible NS
// for QUAD, HEX geometry type
// for LINEAR, QUADRATIC and CUBIC interpolation order
TEST_CASE( "IWG_Incompressible_NS_Velocity_Interface", "[moris],[fem],[IWG_Incompressible_NS_Velocity_Interface]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-4;

    // define a perturbation relative size
    real tPerturbation = 1E-6;

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
        arma::Mat< double > tNormalRand;
        tNormalRand.randu( iSpaceDim, 1 );
        tNormal.matrix_data() = tNormalRand;

        // dof type list
        Cell< MSI::Dof_Type > tVelDofTypes;
        Cell< MSI::Dof_Type > tPDofTypes = { MSI::Dof_Type::P };
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

               // set dof type list
               tVelDofTypes = { MSI::Dof_Type::VX, MSI::Dof_Type::VY };

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

                // set dof type list
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
            arma::Mat< double > tMasterMatrixVis;
            arma::Mat< double > tSlaveMatrixVel;
            arma::Mat< double > tSlaveMatrixP;
            arma::Mat< double > tSlaveMatrixVis;

            // get number of dof
            int tNumDofVel  = 0;
            int tNumDofP    = 0;
            int tNumDofVis  = 0;

            // switch on interpolation order
            switch( iInterpOrder )
            {
                case ( 1 ):
                {
                    // set interpolation type
                    tInterpolationOrder = mtk::Interpolation_Order::LINEAR;

                    // get number of dof
                    tNumDofVel = tNumCoeffs( 0 ) * iSpaceDim;
                    tNumDofP   = tNumCoeffs( 0 );
                    tNumDofVis = tNumCoeffs( 0 );

                    // create random coefficients for master FI
                    tMasterMatrixVel.randu( tNumCoeffs( 0 ), iSpaceDim );
                    tMasterMatrixP.randu( tNumCoeffs( 0 ), 1 );
                    tMasterMatrixVis.randu( tNumCoeffs( 0 ), 1 );

                    // create random coefficients for slave FI
                    tSlaveMatrixVel.randu( tNumCoeffs( 0 ), iSpaceDim );
                    tSlaveMatrixP.randu( tNumCoeffs( 0 ), 1 );
                    tSlaveMatrixVis.randu( tNumCoeffs( 0 ), 1 );
                    break;
                }
                case ( 2 ):
                {
                    // set interpolation type
                    tInterpolationOrder = mtk::Interpolation_Order::QUADRATIC;

                    // get number of dof
                    tNumDofVel = tNumCoeffs( 1 ) * iSpaceDim;
                    tNumDofP   = tNumCoeffs( 1 );
                    tNumDofVis = tNumCoeffs( 1 );

                    // create random coefficients for master FI
                    tMasterMatrixVel.randu( tNumCoeffs( 1 ), iSpaceDim );
                    tMasterMatrixP.randu( tNumCoeffs( 1 ), 1 );
                    tMasterMatrixVis.randu( tNumCoeffs( 1 ), 1 );

                    // create random coefficients for slave FI
                    tSlaveMatrixVel.randu( tNumCoeffs( 1 ), iSpaceDim );
                    tSlaveMatrixP.randu( tNumCoeffs( 1 ), 1 );
                    tSlaveMatrixVis.randu( tNumCoeffs( 1 ), 1 );
                    break;
                }
                case ( 3 ):
                {
                    // set interpolation type
                    tInterpolationOrder = mtk::Interpolation_Order::CUBIC;

                    // get number of dof
                    tNumDofVel = tNumCoeffs( 2 ) * iSpaceDim;
                    tNumDofP   = tNumCoeffs( 2 );
                    tNumDofVis = tNumCoeffs( 2 );

                    // create random coefficients for master FI
                    tMasterMatrixVel.randu( tNumCoeffs( 2 ), iSpaceDim );
                    tMasterMatrixP.randu( tNumCoeffs( 2 ), 1 );
                    tMasterMatrixVis.randu( tNumCoeffs( 2 ), 1 );

                    // create random coefficients for slave FI
                    tSlaveMatrixVel.randu( tNumCoeffs( 2 ), iSpaceDim );
                    tSlaveMatrixP.randu( tNumCoeffs( 2 ), 1 );
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
            Matrix< DDRMat > tMasterDOFHatVel;
            tMasterDOFHatVel.matrix_data() = 10.0 * tMasterMatrixVel;
            Matrix< DDRMat > tMasterDOFHatP;
            tMasterDOFHatP.matrix_data() = 10.0 * tMasterMatrixP;
            Matrix< DDRMat > tMasterDOFHatVis;
            tMasterDOFHatVis.matrix_data() = 10.0 * tMasterMatrixVis;

            // fill random coefficients for slave FI
            Matrix< DDRMat > tSlaveDOFHatVel;
            tSlaveDOFHatVel.matrix_data() = 10.0 * tSlaveMatrixVel;
            Matrix< DDRMat > tSlaveDOFHatP;
            tSlaveDOFHatP.matrix_data() = 10.0 * tSlaveMatrixP;
            Matrix< DDRMat > tSlaveDOFHatVis;
            tSlaveDOFHatVis.matrix_data() = 10.0 * tSlaveMatrixVis;

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tMasterFIs( 3 );

            // create the field interpolator
            tMasterFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelDofTypes );
            tMasterFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tPDofTypes );
            tMasterFIs( 2 ) = new Field_Interpolator( 1, tFIRule, &tGI, tVisDofTypes );

            // set the coefficients uHat
            tMasterFIs( 0 )->set_coeff( tMasterDOFHatVel );
            tMasterFIs( 1 )->set_coeff( tMasterDOFHatP );
            tMasterFIs( 2 )->set_coeff( tMasterDOFHatVis );

            //set the evaluation point xi, tau
            tMasterFIs( 0 )->set_space_time( tParamPoint );
            tMasterFIs( 1 )->set_space_time( tParamPoint );
            tMasterFIs( 2 )->set_space_time( tParamPoint );

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tSlaveFIs( 3 );

            // create the field interpolator
            tSlaveFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelDofTypes );
            tSlaveFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tPDofTypes );
            tSlaveFIs( 2 ) = new Field_Interpolator( 1, tFIRule, &tGI, tVisDofTypes );

            // set the coefficients uHat
            tSlaveFIs( 0 )->set_coeff( tSlaveDOFHatVel );
            tSlaveFIs( 1 )->set_coeff( tSlaveDOFHatP );
            tSlaveFIs( 2 )->set_coeff( tSlaveDOFHatVis );

            //set the evaluation point xi, tau
            tSlaveFIs( 0 )->set_space_time( tParamPoint );
            tSlaveFIs( 1 )->set_space_time( tParamPoint );
            tSlaveFIs( 2 )->set_space_time( tParamPoint );

            // create the master properties
            std::shared_ptr< fem::Property > tPropMasterViscosity = std::make_shared< fem::Property >();
            tPropMasterViscosity->set_parameters( { {{ 1.0 }} } );
            tPropMasterViscosity->set_val_function( tFIConstValFunction_UTVelocityInterface );

            std::shared_ptr< fem::Property > tPropMasterDensity = std::make_shared< fem::Property >();
            tPropMasterDensity->set_parameters( { {{ 1.0 }} } );
            tPropMasterDensity->set_val_function( tFIConstValFunction_UTVelocityInterface );

            // create the slave properties
            std::shared_ptr< fem::Property > tPropSlaveViscosity = std::make_shared< fem::Property >();
            tPropSlaveViscosity->set_parameters( { {{ 1.0 }} } );
            tPropSlaveViscosity->set_val_function( tFIConstValFunction_UTVelocityInterface );

            std::shared_ptr< fem::Property > tPropSlaveDensity = std::make_shared< fem::Property >();
            tPropSlaveDensity->set_parameters( { {{ 1.0 }} } );
            tPropSlaveDensity->set_val_function( tFIConstValFunction_UTVelocityInterface );

            // define constitutive models
            fem::CM_Factory tCMFactory;

            std::shared_ptr< fem::Constitutive_Model > tCMMasterFluid =
                    tCMFactory.create_CM( fem::Constitutive_Type::FLUID_INCOMPRESSIBLE );
            tCMMasterFluid->set_dof_type_list( { tVelDofTypes, tPDofTypes } );
            tCMMasterFluid->set_property( tPropMasterViscosity, "Viscosity" );
            tCMMasterFluid->set_property( tPropMasterDensity, "Density" );
            tCMMasterFluid->set_space_dim( iSpaceDim );

            std::shared_ptr< fem::Constitutive_Model > tCMMasterTurbulence =
                    tCMFactory.create_CM( fem::Constitutive_Type::FLUID_TURBULENCE );
            tCMMasterTurbulence->set_dof_type_list( { tVelDofTypes, tVisDofTypes } );
            tCMMasterTurbulence->set_property( tPropMasterViscosity, "Viscosity" );
            tCMMasterTurbulence->set_property( tPropMasterDensity, "Density" );
            tCMMasterTurbulence->set_space_dim( iSpaceDim );

            std::shared_ptr< fem::Constitutive_Model > tCMSlaveFluid =
                    tCMFactory.create_CM( fem::Constitutive_Type::FLUID_INCOMPRESSIBLE );
            tCMSlaveFluid->set_dof_type_list( { tVelDofTypes, tPDofTypes } );
            tCMSlaveFluid->set_property( tPropSlaveViscosity, "Viscosity" );
            tCMSlaveFluid->set_property( tPropSlaveDensity, "Density" );
            tCMSlaveFluid->set_space_dim( iSpaceDim );

            std::shared_ptr< fem::Constitutive_Model > tCMSlaveTurbulence =
                    tCMFactory.create_CM( fem::Constitutive_Type::FLUID_TURBULENCE );
            tCMSlaveTurbulence->set_dof_type_list( { tVelDofTypes, tVisDofTypes } );
            tCMSlaveTurbulence->set_property( tPropSlaveViscosity, "Viscosity" );
            tCMSlaveTurbulence->set_property( tPropSlaveDensity, "Density" );
            tCMSlaveTurbulence->set_space_dim( iSpaceDim );

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
                    tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_INTERFACE );
            tIWG->set_residual_dof_type( tVelDofTypes );
            tIWG->set_dof_type_list( { tVelDofTypes, tPDofTypes, tVisDofTypes }, mtk::Master_Slave::MASTER );
            tIWG->set_dof_type_list( { tVelDofTypes, tPDofTypes, tVisDofTypes }, mtk::Master_Slave::SLAVE );
            tIWG->set_constitutive_model( tCMMasterFluid, "IncompressibleFluid", mtk::Master_Slave::MASTER );
            tIWG->set_constitutive_model( tCMSlaveFluid, "IncompressibleFluid", mtk::Master_Slave::SLAVE );
            tIWG->set_constitutive_model( tCMMasterTurbulence, "TurbulenceFluid", mtk::Master_Slave::MASTER );
            tIWG->set_constitutive_model( tCMSlaveTurbulence, "TurbulenceFluid", mtk::Master_Slave::SLAVE );
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
            tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )        = 0;
            tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )         = 1;
            tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 2;

            // set size and populate the set master and slave dof type map
            tIWG->mSet->mMasterDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
            tIWG->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )        = 0;
            tIWG->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )         = 1;
            tIWG->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 2;

            tIWG->mSet->mSlaveDofTypeMap .set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
            tIWG->mSet->mSlaveDofTypeMap ( static_cast< int >( MSI::Dof_Type::VX ) )        = 0;
            tIWG->mSet->mSlaveDofTypeMap ( static_cast< int >( MSI::Dof_Type::P ) )         = 1;
            tIWG->mSet->mSlaveDofTypeMap ( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 2;

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( 6 );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofVel-1 } };
            tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofVel, tNumDofVel + tNumDofP - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 2 ) = { { tNumDofVel + tNumDofP, tNumDofVel + tNumDofP + tNumDofVis - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 3 ) = { { tNumDofVel + tNumDofP + tNumDofVis, 2 * tNumDofVel + tNumDofP + tNumDofVis -1 } };
            tIWG->mSet->mResDofAssemblyMap( 4 ) = { { 2 * tNumDofVel + tNumDofP + tNumDofVis, 2 * tNumDofVel + 2 * tNumDofP + tNumDofVis - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 5 ) = { { 2 * tNumDofVel + 2 * tNumDofP + tNumDofVis, 2 * tNumDofVel + 2 * tNumDofP + 2 * tNumDofVis - 1 } };

            // set size and fill the set jacobian assembly map
            Matrix< DDSMat > tJacAssembly = {
                    { 0, tNumDofVel-1 },
                    { tNumDofVel, tNumDofVel + tNumDofP - 1 },
                    { tNumDofVel + tNumDofP, tNumDofVel + tNumDofP + tNumDofVis - 1 },
                    { tNumDofVel + tNumDofP + tNumDofVis, 2 * tNumDofVel + tNumDofP + tNumDofVis -1 },
                    { 2 * tNumDofVel + tNumDofP + tNumDofVis, 2 * tNumDofVel + 2 * tNumDofP + tNumDofVis - 1 },
                    { 2 * tNumDofVel + 2 * tNumDofP + tNumDofVis, 2 * tNumDofVel + 2 * tNumDofP + 2 * tNumDofVis - 1 } };
            tIWG->mSet->mJacDofAssemblyMap.resize( 6 );
            tIWG->mSet->mJacDofAssemblyMap( 0 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 1 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 2 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 3 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 4 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 5 ) = tJacAssembly;

            // set IWG normal
            tIWG->set_normal( tNormal );

            // build global property type list
            tIWG->build_global_dof_and_dv_type_list();

            // populate the requested master dof type
            tIWG->mRequestedMasterGlobalDofTypes = { tVelDofTypes, tPDofTypes, tVisDofTypes };
            tIWG->mRequestedSlaveGlobalDofTypes  = { tVelDofTypes, tPDofTypes, tVisDofTypes };

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
            tIWG->mSet->mResidual( 0 ).set_size( 2 * ( tNumDofVel + tNumDofP + tNumDofVis ), 1 , 0.0 );
            tIWG->mSet->mJacobian.set_size( 2 * ( tNumDofVel + tNumDofP + tNumDofVis ), 2 * ( tNumDofVel + tNumDofP + tNumDofVis ), 0.0 );

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
