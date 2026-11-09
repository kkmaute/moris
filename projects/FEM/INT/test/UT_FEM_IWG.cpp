/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_IWG.cpp
 *
 */

#include <memory>
#include "catch.hpp"
#include "fn_equal_to.hpp"
#include "fn_norm.hpp"
#include "fn_trans.hpp"
#include "cl_FEM_Property.hpp"                 //FEM/INT/src
#include "cl_FEM_Geometry_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_Field_Interpolator.hpp"       //FEM/INT/src
#include "cl_FEM_CM_Factory.hpp"               //FEM/INT/src
#include "cl_FEM_Constitutive_Model.hpp"       //FEM/INT/src

#define protected public
#define private public
#include "cl_FEM_Field_Interpolator_Manager.hpp"    //FEM//INT//src
#include "cl_FEM_IWG.hpp"                           //FEM/INT/src
#include "cl_MSI_Dof_Manager.hpp"                   //FEM/INT/src
#include "cl_MSI_Model_Solver_Interface.hpp"        //FEM/INT/src
#include "cl_FEM_Set.hpp"                           //FEM/INT/src
#undef protected
#undef private

#include "cl_FEM_IWG_Factory.hpp"    //FEM/INT/src

namespace moris::fem
{

    void tValFunction_UTIWG( moris::Matrix< moris::DDRMat >& aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >&        aParameters,
            moris::fem::Field_Interpolator_Manager*          aFIManager )
    {
        moris::Matrix< moris::DDRMat > tPropertyVal( 1, 1, 1.0 );
        aPropMatrix = tPropertyVal;
    }

    // This test case tests all the member functions of the IWG,
    // except the residual and jacobian related ones.
    TEST_CASE( "IWG_dof_dv", "[moris],[fem],[IWG_dof_dv]" )
    {
        // create the properties
        std::shared_ptr< fem::Property > tPropLeader1 = std::make_shared< fem::Property >();
        tPropLeader1->set_dof_type_list( { { MSI::Dof_Type::TEMP }, { MSI::Dof_Type::VX } } );
        tPropLeader1->set_val_function( tValFunction_UTIWG );

        std::shared_ptr< fem::Property > tPropLeader2 = std::make_shared< fem::Property >();
        tPropLeader2->set_dof_type_list( { { MSI::Dof_Type::TEMP }, { MSI::Dof_Type::LS1 } } );
        //        tPropLeader2->set_dv_type_list( {{ gen::PDV_Type::LS1 }, { gen::PDV_Type::LS2 }} );
        tPropLeader2->set_val_function( tValFunction_UTIWG );

        std::shared_ptr< fem::Property > tPropFollower1 = std::make_shared< fem::Property >();
        tPropFollower1->set_val_function( tValFunction_UTIWG );

        std::shared_ptr< fem::Property > tPropFollower2 = std::make_shared< fem::Property >();
        tPropFollower2->set_dof_type_list( { { MSI::Dof_Type::TEMP }, { MSI::Dof_Type::LS1 } } );
        //        tPropFollower2->set_dv_type_list( {{ gen::PDV_Type::LS1 }, { gen::PDV_Type::LS2 }} );
        tPropFollower2->set_val_function( tValFunction_UTIWG );

        // define constitutive models
        fem::CM_Factory tCMFactory;

        std::shared_ptr< fem::Constitutive_Model > tCMLeader1 =
                tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
        tCMLeader1->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
        tCMLeader1->set_property( tPropLeader2, "Conductivity" );
        tCMLeader1->set_space_dim( 3 );
        tCMLeader1->set_local_properties();

        std::shared_ptr< fem::Constitutive_Model > tCMFollower1 =
                tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
        tCMFollower1->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
        tCMFollower1->set_property( tPropFollower2, "Conductivity" );
        tCMFollower1->set_space_dim( 3 );
        tCMFollower1->set_local_properties();

        // create leader dof field interpolators
        uint                          tNumberOfFields = 1;
        Vector< Field_Interpolator* > tLeaderDofFI( 4, nullptr );
        tLeaderDofFI( 1 ) = new Field_Interpolator( tNumberOfFields, { MSI::Dof_Type::TEMP } );
        tLeaderDofFI( 2 ) = new Field_Interpolator( tNumberOfFields, { MSI::Dof_Type::LS1 } );
        tLeaderDofFI( 3 ) = new Field_Interpolator( tNumberOfFields, { MSI::Dof_Type::VX } );

        //        // create leader dv field interpolators
        //        Vector< Field_Interpolator* > tLeaderDvFI( 3, nullptr );
        //        tLeaderDvFI( 0 ) = new Field_Interpolator ( tNumberOfFields, { gen::PDV_Type::DENSITY } );
        //        tLeaderDvFI( 1 ) = new Field_Interpolator ( tNumberOfFields, { gen::PDV_Type::LS1 } );
        //        tLeaderDvFI( 2 ) = new Field_Interpolator ( tNumberOfFields, { gen::PDV_Type::LS2 } );

        //        // create follower dof field interpolators
        //        Vector< Field_Interpolator* > tFollowerDofFI( 4, nullptr );
        //        tFollowerDofFI( 0 ) = new Field_Interpolator ( tNumberOfFields, { MSI::Dof_Type::UX } );
        //        tFollowerDofFI( 1 ) = new Field_Interpolator ( tNumberOfFields, { MSI::Dof_Type::TEMP } );
        //        tFollowerDofFI( 2 ) = new Field_Interpolator ( tNumberOfFields, { MSI::Dof_Type::LS1 } );

        //        // create follower dv field interpolators
        //        Vector< Field_Interpolator* > tFollowerDvFI( 2, nullptr );
        //        tFollowerDvFI( 0 ) = new Field_Interpolator ( tNumberOfFields, { gen::PDV_Type::LS1 } );
        //        tFollowerDvFI( 1 ) = new Field_Interpolator ( tNumberOfFields, { gen::PDV_Type::LS2 } );

        // define the IWGs
        fem::IWG_Factory tIWGFactory;

        std::shared_ptr< moris::fem::IWG > tIWG = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
        tIWG->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
        tIWG->set_dof_type_list( { { MSI::Dof_Type::TEMP } }, mtk::Leader_Follower::LEADER );
        tIWG->set_dof_type_list( { { MSI::Dof_Type::TEMP }, { MSI::Dof_Type::UX } }, mtk::Leader_Follower::FOLLOWER );
        // tIWG->set_dv_type_list( {{ gen::PDV_Type::DENSITY }, {gen::PDV_Type::LS1 }}, mtk::Leader_Follower::LEADER );
        // tIWG->set_dv_type_list( {{gen::PDV_Type::LS1 }}, mtk::Leader_Follower::FOLLOWER );
        tIWG->set_constitutive_model( tCMLeader1, "Diffusion", mtk::Leader_Follower::LEADER );
        // tIWG->set_constitutive_model( tCMFollower1, "Diffusion", mtk::Leader_Follower::FOLLOWER );
        tIWG->set_property( tPropLeader1, "Load", mtk::Leader_Follower::LEADER );
        // tIWG->set_property( tPropFollower1, "Load", mtk::Leader_Follower::FOLLOWER );

        tIWG->mRequestedLeaderGlobalDofTypes = { { MSI::Dof_Type::TEMP }, { MSI::Dof_Type::LS1 }, { MSI::Dof_Type::VX } };
        // tIWG->mRequestedFollowerGlobalDofTypes = {{ MSI::Dof_Type::UX },{ MSI::Dof_Type::TEMP},{ MSI::Dof_Type::LS1}};
        MSI::Equation_Set* tSet = new fem::Set();

        static_cast< fem::Set* >( tSet )->set_set_type( fem::Element_Type::BULK );

        tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

        tIWG->mSet->mUniqueDofTypeList.resize( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, MSI::Dof_Type::END_ENUM );

        tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
        tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) )   = 0;
        tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 1;
        tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::LS1 ) )  = 2;
        tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )   = 3;

        tIWG->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
        tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) )   = 0;
        tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 1;
        tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::LS1 ) )  = 2;
        tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )   = 3;

        //        tIWG->mSet->mFollowerDofTypeMap.set_size( static_cast< int >(MSI::Dof_Type::UNDEFINED) + 1, 1, -1 );
        //        tIWG->mSet->mFollowerDofTypeMap( static_cast< int >(MSI::Dof_Type::UX) ) = 0;
        //        tIWG->mSet->mFollowerDofTypeMap( static_cast< int >(MSI::Dof_Type::TEMP) ) = 1;
        //        tIWG->mSet->mFollowerDofTypeMap( static_cast< int >(MSI::Dof_Type::LS1) ) = 2;
        //        tIWG->mSet->mFollowerDofTypeMap( static_cast< int >(MSI::Dof_Type::VX) ) = 3;

        tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

        Vector< Vector< enum MSI::Dof_Type > > tDummy;
        Field_Interpolator_Manager             tFIManager( tDummy, tSet );
        tFIManager.mFI = tLeaderDofFI;
        //        tFIManager.mFollowerFI = tFollowerDofFI;

        // build leader and follower global dof type list
        tIWG->get_global_dof_type_list();

        // set IWG field interpolators
        tIWG->mLeaderFIManager = &tFIManager;
        //
        //        // build leader and follower global dv type list
        //        tIWG->build_global_dv_type_list();
        //
        //        // set IWG leader and follower dv field interpolators
        //        tIWG->set_dv_field_interpolators( tLeaderDvFI );
        //        tIWG->set_dv_field_interpolators( tFollowerDvFI, mtk::Leader_Follower::FOLLOWER );

        // dof check--------------------------------------------------------------------
        // check leader global dof list size
        CHECK( equal_to( tIWG->mLeaderGlobalDofTypes.size(), 3 ) );

        // check leader global dof list content
        CHECK( equal_to( static_cast< uint >( tIWG->mLeaderGlobalDofTypes( 0 )( 0 ) ), 3 ) );
        CHECK( equal_to( static_cast< uint >( tIWG->mLeaderGlobalDofTypes( 1 )( 0 ) ), 14 ) );
        CHECK( equal_to( static_cast< uint >( tIWG->mLeaderGlobalDofTypes( 2 )( 0 ) ), 6 ) );

        //        // check follower global dof list size
        //        CHECK( equal_to( tIWG->mFollowerGlobalDofTypes.size(), 3 ));
        //
        //        // check follower global dof list content
        //        CHECK( equal_to( static_cast< uint >( tIWG->mFollowerGlobalDofTypes( 0 )( 0 ) ), 3 ) );
        //        CHECK( equal_to( static_cast< uint >( tIWG->mFollowerGlobalDofTypes( 1 )( 0 ) ), 0 ) );
        //        CHECK( equal_to( static_cast< uint >( tIWG->mFollowerGlobalDofTypes( 2 )( 0 ) ), 6 ) );

        // check dof field interpolators
        tIWG->check_field_interpolators();
        tIWG->check_field_interpolators( mtk::Leader_Follower::FOLLOWER );

        //        // dv check---------------------------------------------------------------------
        //        // check leader global dv list size
        //        CHECK( equal_to( tIWG->mLeaderGlobalDvTypes.size(), 3 ));
        //
        //        // check leader global dv list content
        //        CHECK( equal_to( static_cast< uint >( tIWG->mLeaderGlobalDvTypes( 0 )( 0 ) ), 2 ) );
        //        CHECK( equal_to( static_cast< uint >( tIWG->mLeaderGlobalDvTypes( 1 )( 0 ) ), 0 ) );
        //        CHECK( equal_to( static_cast< uint >( tIWG->mLeaderGlobalDvTypes( 2 )( 0 ) ), 1 ) );
        //
        //        // check follower global dv list size
        //        CHECK( equal_to( tIWG->mFollowerGlobalDvTypes.size(), 2 ));
        //
        //        // check leader global dv list content
        //        CHECK( equal_to( static_cast< uint >( tIWG->mFollowerGlobalDvTypes( 0 )( 0 ) ), 0 ) );
        //        CHECK( equal_to( static_cast< uint >( tIWG->mFollowerGlobalDvTypes( 1 )( 0 ) ), 1 ) );
        //
        //        // check dv field interpolators
        //        tIWG->check_dv_field_interpolators();
        //        tIWG->check_dv_field_interpolators( mtk::Leader_Follower::FOLLOWER );

        // clean up---------------------------------------------------------------------
        // delete the dof field interpolator pointers
        tLeaderDofFI.clear();

        //        // delete the dof field interpolator pointers
        //        tFollowerDofFI.clear();

        // delete the dv field interpolator pointers
        //        for( Field_Interpolator* tFI : tLeaderDvFI )
        //        {
        //            delete tFI;
        //        }
        //        tLeaderDvFI.clear();

        // delete the dv field interpolator pointers
        //        for( Field_Interpolator* tFI : tFollowerDvFI )
        //        {
        //            delete tFI;
        //        }
        //        tFollowerDvFI.clear();

        delete ( tSet );

    } /* TEST_CASE */
}    // namespace moris::fem

namespace moris::fem
{
    TEST_CASE( "IWG Spatially Lumped Mass Construction", "[fem],[IWG_Lumped_Mass]" )
    {
        // Compare the shared builder with explicit dense block row sums. Different
        // current/previous points make the cross operator nonsymmetric and exercise
        // temporal off-diagonal blocks for scalar and vector QUAD/HEX fields.
        for ( uint tDimension : { 2u, 3u } )
        {
            const auto              tGeometry = tDimension == 2 ? mtk::Geometry_Type::QUAD : mtk::Geometry_Type::HEX;
            mtk::Interpolation_Rule tGeometryRule( tGeometry, mtk::Interpolation_Type::LAGRANGE, mtk::Interpolation_Order::LINEAR, mtk::Interpolation_Type::LAGRANGE, mtk::Interpolation_Order::LINEAR );

            Geometry_Interpolator tGI( tGeometryRule );
            Matrix< DDRMat >      tXHat = tDimension == 2
                                                ? Matrix< DDRMat >{ { 0.0, 0.0 }, { 1.0, 0.0 }, { 1.0, 1.0 }, { 0.0, 1.0 } }
                                                : Matrix< DDRMat >{ { 0.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 }, { 1.0, 1.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 }, { 1.0, 0.0, 1.0 }, { 1.0, 1.0, 1.0 }, { 0.0, 1.0, 1.0 } };
            tGI.set_coeff( tXHat, { { 0.0 }, { 1.0 } } );

            for ( auto tSpaceOrder : { mtk::Interpolation_Order::LINEAR, mtk::Interpolation_Order::QUADRATIC } )
            {
                for ( auto tTimeOrder : { mtk::Interpolation_Order::CONSTANT, mtk::Interpolation_Order::LINEAR } )
                {
                    for ( uint tNumFields : { 1u, 3u } )
                    {
                        mtk::Interpolation_Rule tRule( tGeometry, mtk::Interpolation_Type::LAGRANGE, tSpaceOrder, tTimeOrder == mtk::Interpolation_Order::CONSTANT ? mtk::Interpolation_Type::CONSTANT : mtk::Interpolation_Type::LAGRANGE, tTimeOrder );

                        const Vector< MSI::Dof_Type > tDofs = tNumFields == 1
                                                                    ? Vector< MSI::Dof_Type >{ MSI::Dof_Type::TEMP }
                                                                    : Vector< MSI::Dof_Type >{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ };
                        Field_Interpolator            tCurrent( tNumFields, tRule, &tGI, tDofs );
                        Field_Interpolator            tPrevious( tNumFields, tRule, &tGI, tDofs );
                        const uint                    tNumSpace        = tCurrent.get_number_of_space_bases();
                        const uint                    tNumTime         = tCurrent.get_number_of_time_bases();
                        const uint                    tNumCoefficients = tNumFields * tNumTime * tNumSpace;
                        Matrix< DDRMat >              tCoefficients( tNumTime * tNumSpace, tNumFields, 0.7 );
                        tCurrent.set_coeff( tCoefficients );
                        tPrevious.set_coeff( tCoefficients );

                        for ( bool tInterior : { false, true } )
                        {
                            INFO( "Dimension=" << tDimension << " spatial order=" << static_cast< uint >( tSpaceOrder )
                                               << " time bases=" << tNumTime << " fields=" << tNumFields << " interior=" << tInterior );
                            Matrix< DDRMat > tCurrentPoint( tDimension + 1, 1 );
                            Matrix< DDRMat > tPreviousPoint( tDimension + 1, 1 );

                            for ( uint i = 0; i < tDimension; ++i )
                            {
                                tCurrentPoint( i )  = 0.2 - 0.1 * i;
                                tPreviousPoint( i ) = -0.3 + 0.1 * i;
                            }

                            tCurrentPoint( tDimension )  = tInterior ? 0.2 : -1.0;
                            tPreviousPoint( tDimension ) = tInterior ? -0.4 : 1.0;
                            tCurrent.set_space_time( tCurrentPoint );
                            tPrevious.set_space_time( tPreviousPoint );

                            const Matrix< DDRMat > tDenseCurrent = trans( tCurrent.N() ) * tCurrent.N();
                            const Matrix< DDRMat > tDenseCross   = trans( tCurrent.N() ) * tPrevious.N();
                            Matrix< DDRMat >       tReferenceCurrent( tNumCoefficients, tNumCoefficients, 0.0 );
                            Matrix< DDRMat >       tReferenceCross( tNumCoefficients, tNumCoefficients, 0.0 );

                            for ( uint iField = 0; iField < tNumFields; ++iField )
                            {
                                for ( uint iTime = 0; iTime < tNumTime; ++iTime )
                                {
                                    for ( uint jTime = 0; jTime < tNumTime; ++jTime )
                                    {
                                        const uint tRowStart    = ( iField * tNumTime + iTime ) * tNumSpace;
                                        const uint tColumnStart = ( iField * tNumTime + jTime ) * tNumSpace;

                                        for ( uint iSpace = 0; iSpace < tNumSpace; ++iSpace )
                                        {
                                            for ( uint jSpace = 0; jSpace < tNumSpace; ++jSpace )
                                            {
                                                tReferenceCurrent( tRowStart + iSpace, tColumnStart + iSpace ) +=
                                                        tDenseCurrent( tRowStart + iSpace, tColumnStart + jSpace );
                                                tReferenceCross( tRowStart + iSpace, tColumnStart + iSpace ) +=
                                                        tDenseCross( tRowStart + iSpace, tColumnStart + jSpace );
                                            }
                                        }
                                    }
                                }
                            }

                            IWG::Lumped_Matrix tCurrentLump;
                            IWG::Lumped_Matrix tPreviousLump;
                            IWG::compute_lumped_matrix(
                                    &tCurrent, &tPrevious, tPreviousLump );

                            Matrix< DDRMat > tCross( tNumCoefficients, tNumCoefficients, 0.0 );
                            tCross += tPreviousLump;

                            CHECK( norm( tCross - tReferenceCross ) < 1e-12 );

                            Matrix< DDRMat > tInitialValues( tNumFields, 1 );

                            for ( uint iField = 0; iField < tNumFields; ++iField )
                            {
                                tInitialValues( iField ) = -0.5 + 1.3 * iField;
                            }

                            Matrix< DDRMat > tInitialCoefficients;
                            IWG::compute_lumped_matrices(
                                    &tCurrent, &tPrevious, tCurrentLump, tPreviousLump, tInitialCoefficients, tInitialValues, true );

                            Matrix< DDRMat > tActualCurrent( tNumCoefficients, tNumCoefficients, 0.0 );
                            Matrix< DDRMat > tActualPrevious = tActualCurrent;
                            tActualCurrent += tCurrentLump;
                            tActualPrevious += tPreviousLump;

                            CHECK( norm( tActualCurrent - tReferenceCurrent ) < 1e-12 );
                            CHECK( norm( tActualPrevious - ( tNumTime == 1 ? tReferenceCurrent : tReferenceCross ) ) < 1e-12 );
                            CHECK( tCurrentLump.num_stored_entries() == tNumSpace * tNumFields * tNumTime * tNumTime );
                            CHECK( tInitialCoefficients.n_rows() == tNumTime * tNumSpace );
                            CHECK( tInitialCoefficients.n_cols() == tNumFields );

                            for ( uint iField = 0; iField < tNumFields; ++iField )
                            {
                                for ( uint iBasis = 0; iBasis < tNumSpace * tNumTime; ++iBasis )
                                {
                                    CHECK( tInitialCoefficients( iBasis, iField ) == Approx( tInitialValues( iField ) ) );
                                }
                            }

                            const Matrix< DDRMat > tInitialProduct = tPreviousLump * vectorize( tInitialCoefficients );

                            CHECK( norm( tInitialProduct - tActualPrevious * vectorize( tInitialCoefficients ) ) < 1e-12 );

                            // No initial data or previous interpolator is needed for
                            // the constant-time previous=current convention.
                            IWG::compute_lumped_matrices(
                                    &tCurrent, tNumTime == 1 ? nullptr : &tPrevious, tCurrentLump, tPreviousLump, tInitialCoefficients, {}, false );

                            CHECK( tInitialCoefficients.numel() == 0 );
                            REQUIRE_THROWS( IWG::compute_lumped_matrix(
                                    nullptr, &tPrevious, tPreviousLump ) );
                            REQUIRE_THROWS( IWG::compute_lumped_matrices(
                                    &tCurrent, &tPrevious, tCurrentLump, tPreviousLump, tInitialCoefficients, {}, true ) );
                        }
                    }
                }
            }
        }
    }

}    // namespace moris::fem

namespace moris::fem
{
    TEST_CASE( "IWG Lumped Mass Validation", "[fem],[IWG_Lumped_Mass]" )
    {
        // Reject incompatible field layouts and unsupported temporal orders before
        // evaluating shape functions or allocating the compact operator.
        mtk::Interpolation_Rule tLinearRule( mtk::Geometry_Type::QUAD, mtk::Interpolation_Type::LAGRANGE, mtk::Interpolation_Order::LINEAR, mtk::Interpolation_Type::LAGRANGE, mtk::Interpolation_Order::LINEAR );
        mtk::Interpolation_Rule tConstantRule( mtk::Geometry_Type::QUAD, mtk::Interpolation_Type::LAGRANGE, mtk::Interpolation_Order::LINEAR, mtk::Interpolation_Type::CONSTANT, mtk::Interpolation_Order::CONSTANT );
        mtk::Interpolation_Rule tQuadraticSpaceRule( mtk::Geometry_Type::QUAD, mtk::Interpolation_Type::LAGRANGE, mtk::Interpolation_Order::QUADRATIC, mtk::Interpolation_Type::LAGRANGE, mtk::Interpolation_Order::LINEAR );
        mtk::Interpolation_Rule tQuadraticTimeRule( mtk::Geometry_Type::QUAD, mtk::Interpolation_Type::LAGRANGE, mtk::Interpolation_Order::LINEAR, mtk::Interpolation_Type::LAGRANGE, mtk::Interpolation_Order::QUADRATIC );

        Geometry_Interpolator   tGI( tLinearRule );
        Field_Interpolator      tScalar( 1, tLinearRule, &tGI, { MSI::Dof_Type::TEMP } );
        Field_Interpolator      tVector( 3, tLinearRule, &tGI, { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ } );
        Field_Interpolator      tConstant( 1, tConstantRule, &tGI, { MSI::Dof_Type::TEMP } );
        Field_Interpolator      tQuadraticSpace( 1, tQuadraticSpaceRule, &tGI, { MSI::Dof_Type::TEMP } );
        Field_Interpolator      tQuadraticTime( 1, tQuadraticTimeRule, &tGI, { MSI::Dof_Type::TEMP } );
        IWG::Lumped_Matrix      tLump;

        REQUIRE_THROWS( IWG::compute_lumped_matrix(
                &tScalar, &tVector, tLump ) );
        REQUIRE_THROWS( IWG::compute_lumped_matrix(
                &tScalar, &tConstant, tLump ) );
        REQUIRE_THROWS( IWG::compute_lumped_matrix(
                &tScalar, &tQuadraticSpace, tLump ) );
        REQUIRE_THROWS( IWG::compute_lumped_matrix(
                &tQuadraticTime, &tQuadraticTime, tLump ) );
        CHECK( tLump.num_stored_entries() == 0 );
    }

}    // namespace moris::fem
