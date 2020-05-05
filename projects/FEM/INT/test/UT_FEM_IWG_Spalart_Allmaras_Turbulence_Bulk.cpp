#include <string>
#include <catch.hpp>
#include <memory>

#include "assert.hpp"

#include "cl_MTK_Enums.hpp" //MTK/src
#include "cl_FEM_Enums.hpp"                                //FEM//INT/src

#include "op_equal_equal.hpp"

#define protected public
#define private   public
#include "cl_FEM_Field_Interpolator_Manager.hpp"                   //FEM//INT//src
#include "cl_FEM_IWG.hpp"         //FEM/INT/src
#include "cl_FEM_Set.hpp"         //FEM/INT/src
#undef protected
#undef private

#include "cl_FEM_Field_Interpolator.hpp"                   //FEM//INT//src
#include "cl_FEM_Property.hpp"                   //FEM//INT//src
#include "cl_FEM_CM_Factory.hpp"                   //FEM//INT//src
#include "cl_FEM_SP_Factory.hpp"                   //FEM//INT//src
#include "cl_FEM_IWG_Factory.hpp"                   //FEM//INT//src

void tConstValFunction_SATurbulenceBulk
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 );
}

void tTEMPFIValFunction_SATurbulenceBulk
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 ) * aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::VISCOSITY )->val();
}

void tTEMPFIDerFunction_SATurbulenceBulk
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 ) * aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::VISCOSITY )->N();
}

using namespace moris;
using namespace fem;

TEST_CASE( "IWG_Spalart_Allmaras_Turbulence_Bulk", "[IWG_Spalart_Allmaras_Turbulence_Bulk]" )
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
        moris::Cell< MSI::Dof_Type > tVisDofTypes = { MSI::Dof_Type::VISCOSITY };

        // switch on space dimension
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
               tNumCoeffs = {{ 8 },{ 18 },{ 32 }};

               // set velocity dof types
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
            arma::Mat< double > tMasterMatrixVis;

            // get number of dof
            int tNumDofVel = 0;
            int tNumDofVis = 0;

            // switch on interpolation order
            switch( iInterpOrder )
            {
                case ( 1 ):
                {
                    // set interpolation type
                    tInterpolationOrder = mtk::Interpolation_Order::LINEAR;

                    // get number of dof
                    tNumDofVel = tNumCoeffs( 0 ) * iSpaceDim;
                    tNumDofVis = tNumCoeffs( 0 );

                    // create random coefficients for master FI
                    tMasterMatrixVel.randu( tNumCoeffs( 0 ), iSpaceDim );
                    tMasterMatrixVis.randu( tNumCoeffs( 0 ), 1 );
                    break;
                }
                case ( 2 ):
                {
                    // set interpolation type
                    tInterpolationOrder = mtk::Interpolation_Order::QUADRATIC;

                    // get number of dof
                    tNumDofVel = tNumCoeffs( 1 ) * iSpaceDim;
                    tNumDofVis = tNumCoeffs( 1 );

                    // create random coefficients for master FI
                    tMasterMatrixVel.randu( tNumCoeffs( 1 ), iSpaceDim );
                    tMasterMatrixVis.randu( tNumCoeffs( 1 ), 1 );

                    break;
                }
                case ( 3 ):
                {
                    // set interpolation type
                    tInterpolationOrder = mtk::Interpolation_Order::CUBIC;

                    // get number of dof
                    tNumDofVel = tNumCoeffs( 2 ) * iSpaceDim;
                    tNumDofVis = tNumCoeffs( 2 );

                    // create random coefficients for master FI
                    tMasterMatrixVel.randu( tNumCoeffs( 2 ), iSpaceDim );
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
            Matrix< DDRMat > tMasterDOFHatVis;
            tMasterDOFHatVis.matrix_data() = 10.0 * tMasterMatrixVis;

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tMasterFIs( 2 );

            // create the field interpolator velocity
            tMasterFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelDofTypes );
            tMasterFIs( 0 )->set_coeff( tMasterDOFHatVel );
            tMasterFIs( 0 )->set_space_time( tParamPoint );

            // create the field interpolator viscosity
            tMasterFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tVisDofTypes );
            tMasterFIs( 1 )->set_coeff( tMasterDOFHatVis );
            tMasterFIs( 1 )->set_space_time( tParamPoint );

            // create the properties
            std::shared_ptr< fem::Property > tPropWallDistance = std::make_shared< fem::Property >();
            tPropWallDistance->set_parameters( { {{ 1.0 }} } );
            tPropWallDistance->set_val_function( tConstValFunction_SATurbulenceBulk );

            std::shared_ptr< fem::Property > tPropViscosity = std::make_shared< fem::Property >();
            tPropViscosity->set_parameters( { {{ 2.0 }} } );
            tPropViscosity->set_val_function( tConstValFunction_SATurbulenceBulk );

            // define the IWGs
            fem::IWG_Factory tIWGFactory;

            std::shared_ptr< fem::IWG > tIWG = tIWGFactory.create_IWG( fem::IWG_Type::SPALART_ALLMARAS_TURBULENCE_BULK );
            tIWG->set_residual_dof_type( tVisDofTypes );
            tIWG->set_dof_type_list( { tVelDofTypes, tVisDofTypes }, mtk::Master_Slave::MASTER );
            tIWG->set_property( tPropWallDistance, "WallDistance" );
            tIWG->set_property( tPropViscosity, "Viscosity" );

            // set a fem set pointer
            MSI::Equation_Set * tSet = new fem::Set();
            tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

            // set size for the set EqnObjDofTypeList
            tIWG->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

            // set size and populate the set dof type map
            tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
            tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) ) = 0;
            tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 1;

            // set size and populate the set master dof type map
            tIWG->mSet->mMasterDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
            tIWG->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) ) = 0;
            tIWG->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 1;

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( 2 );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofVel-1 } };
            tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofVel, tNumDofVel + tNumDofVis - 1 } };

            // set size and fill the set jacobian assembly map
            tIWG->mSet->mJacDofAssemblyMap.resize( 2 );
            tIWG->mSet->mJacDofAssemblyMap( 0 ) = { { 0, tNumDofVel-1 }, { tNumDofVel, tNumDofVel + tNumDofVis - 1 } };
            tIWG->mSet->mJacDofAssemblyMap( 1 ) = { { 0, tNumDofVel-1 }, { tNumDofVel, tNumDofVel + tNumDofVis - 1 } };

            // set size and init the set residual and jacobian
            tIWG->mSet->mResidual.resize( 1 );
            tIWG->mSet->mResidual( 0 ).set_size( tNumDofVel + tNumDofVis, 1, 0.0 );
            tIWG->mSet->mJacobian.set_size( tNumDofVel + tNumDofVis, tNumDofVel + tNumDofVis, 0.0 );

            // build global dof type list
            tIWG->get_global_dof_type_list();

            // populate the requested master dof type
            tIWG->mRequestedMasterGlobalDofTypes = { { MSI::Dof_Type::VX }, { MSI::Dof_Type::VISCOSITY }};

            // create a field interpolator manager
            moris::Cell< moris::Cell< enum MSI::Dof_Type > > tDummyDof;
            moris::Cell< moris::Cell< enum PDV > > tDummyDv;
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

//            print( tJacobian(   { 0, tNumDofVis-1 }, { 0, tNumDofVel-1 } ), "tJacobianVU" );
//            print( tJacobianFD( { 0, tNumDofVis-1 }, { 0, tNumDofVel-1 } ), "tJacobianFDVU" );
//            print( tJacobian(   { 0, tNumDofVis-1 }, { tNumDofVel, tNumDofVel + tNumDofVis - 1 }), "tJacobianVV" );
//            print( tJacobianFD( { 0, tNumDofVis-1 }, { tNumDofVel, tNumDofVel + tNumDofVis - 1 }), "tJacobianFDVV" );

            // require check is true
            REQUIRE( tCheckJacobian );

            // clean up
            tMasterFIs.clear();
        }
    }

}/*END_TEST_CASE*/
