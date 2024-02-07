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
#include "cl_FEM_Property.hpp"              //FEM/INT/src
#include "cl_FEM_Geometry_Interpolator.hpp"         //FEM/INT/src
#include "cl_FEM_Field_Interpolator.hpp"         //FEM/INT/src
#include "cl_FEM_CM_Factory.hpp"         //FEM/INT/src
#include "cl_FEM_Constitutive_Model.hpp"         //FEM/INT/src

#define protected public
#define private   public
#include "cl_FEM_Field_Interpolator_Manager.hpp"                   //FEM//INT//src
#include "cl_FEM_IWG.hpp"         //FEM/INT/src
#include "cl_MSI_Dof_Manager.hpp"         //FEM/INT/src
#include "cl_MSI_Model_Solver_Interface.hpp"         //FEM/INT/src
#include "cl_FEM_Set.hpp"         //FEM/INT/src
#undef protected
#undef private

#include "cl_FEM_IWG_Factory.hpp"         //FEM/INT/src

void tValFunction_UTIWG
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Vector< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    moris::Matrix< moris::DDRMat > tPropertyVal( 1, 1, 1.0);
    aPropMatrix = tPropertyVal;
}

namespace moris
{
    namespace fem
    {

    // This test case tests all the member functions of the IWG,
    // except the residual and jacobian related ones.
    TEST_CASE( "IWG_dof_dv", "[moris],[fem],[IWG_dof_dv]" )
    {
        // create the properties
        std::shared_ptr< fem::Property > tPropLeader1 = std::make_shared< fem::Property > ();
        tPropLeader1->set_dof_type_list( {{ MSI::Dof_Type::TEMP }, { MSI::Dof_Type::VX }} );
        tPropLeader1->set_val_function( tValFunction_UTIWG );

        std::shared_ptr< fem::Property > tPropLeader2 = std::make_shared< fem::Property > ();
        tPropLeader2->set_dof_type_list( {{ MSI::Dof_Type::TEMP }, { MSI::Dof_Type::LS1 }} );
//        tPropLeader2->set_dv_type_list( {{ gen::PDV_Type::LS1 }, { gen::PDV_Type::LS2 }} );
        tPropLeader2->set_val_function( tValFunction_UTIWG );

        std::shared_ptr< fem::Property > tPropFollower1 = std::make_shared< fem::Property > ();
        tPropFollower1->set_val_function( tValFunction_UTIWG );

        std::shared_ptr< fem::Property > tPropFollower2 = std::make_shared< fem::Property > ();
        tPropFollower2->set_dof_type_list( {{ MSI::Dof_Type::TEMP }, { MSI::Dof_Type::LS1 }} );
//        tPropFollower2->set_dv_type_list( {{ gen::PDV_Type::LS1 }, { gen::PDV_Type::LS2 }} );
        tPropFollower2->set_val_function( tValFunction_UTIWG );

        // define constitutive models
        fem::CM_Factory tCMFactory;

        std::shared_ptr< fem::Constitutive_Model > tCMLeader1 =
                tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
        tCMLeader1->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tCMLeader1->set_property( tPropLeader2, "Conductivity" );
        tCMLeader1->set_space_dim( 3 );
        tCMLeader1->set_local_properties();

        std::shared_ptr< fem::Constitutive_Model > tCMFollower1 =
                tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
        tCMFollower1->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tCMFollower1->set_property( tPropFollower2, "Conductivity" );
        tCMFollower1->set_space_dim( 3 );
        tCMFollower1->set_local_properties();

        // create leader dof field interpolators
        uint tNumberOfFields = 1;
        Vector< Field_Interpolator* > tLeaderDofFI( 4, nullptr );
        tLeaderDofFI( 1 ) = new Field_Interpolator ( tNumberOfFields, { MSI::Dof_Type::TEMP } );
        tLeaderDofFI( 2 ) = new Field_Interpolator ( tNumberOfFields, { MSI::Dof_Type::LS1 } );
        tLeaderDofFI( 3 ) = new Field_Interpolator ( tNumberOfFields, { MSI::Dof_Type::VX } );

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
        tIWG->set_dof_type_list( {{ MSI::Dof_Type::TEMP }}, mtk::Leader_Follower::LEADER );
        tIWG->set_dof_type_list( {{ MSI::Dof_Type::TEMP },{ MSI::Dof_Type::UX }}, mtk::Leader_Follower::FOLLOWER );
        //tIWG->set_dv_type_list( {{ gen::PDV_Type::DENSITY }, {gen::PDV_Type::LS1 }}, mtk::Leader_Follower::LEADER );
        //tIWG->set_dv_type_list( {{gen::PDV_Type::LS1 }}, mtk::Leader_Follower::FOLLOWER );
        tIWG->set_constitutive_model( tCMLeader1, "Diffusion", mtk::Leader_Follower::LEADER );
        //tIWG->set_constitutive_model( tCMFollower1, "Diffusion", mtk::Leader_Follower::FOLLOWER );
        tIWG->set_property( tPropLeader1, "Load", mtk::Leader_Follower::LEADER );
        //tIWG->set_property( tPropFollower1, "Load", mtk::Leader_Follower::FOLLOWER );

        tIWG->mRequestedLeaderGlobalDofTypes = {{ MSI::Dof_Type::TEMP },{ MSI::Dof_Type::LS1},{ MSI::Dof_Type::VX}};
        //tIWG->mRequestedFollowerGlobalDofTypes = {{ MSI::Dof_Type::UX },{ MSI::Dof_Type::TEMP},{ MSI::Dof_Type::LS1}};
        MSI::Equation_Set * tSet = new fem::Set();

        static_cast<fem::Set*>(tSet)->set_set_type( fem::Element_Type::BULK );

        tIWG->set_set_pointer(static_cast<fem::Set*>(tSet));

        tIWG->mSet->mUniqueDofTypeList.resize( static_cast< int >(MSI::Dof_Type::END_ENUM) + 1, MSI::Dof_Type::END_ENUM );

        tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >(MSI::Dof_Type::END_ENUM) + 1, 1, -1 );
        tIWG->mSet->mUniqueDofTypeMap( static_cast< int >(MSI::Dof_Type::UX) ) = 0;
        tIWG->mSet->mUniqueDofTypeMap( static_cast< int >(MSI::Dof_Type::TEMP) ) = 1;
        tIWG->mSet->mUniqueDofTypeMap( static_cast< int >(MSI::Dof_Type::LS1) ) = 2;
        tIWG->mSet->mUniqueDofTypeMap( static_cast< int >(MSI::Dof_Type::VX) ) = 3;

        tIWG->mSet->mLeaderDofTypeMap.set_size( static_cast< int >(MSI::Dof_Type::END_ENUM) + 1, 1, -1 );
        tIWG->mSet->mLeaderDofTypeMap( static_cast< int >(MSI::Dof_Type::UX) ) = 0;
        tIWG->mSet->mLeaderDofTypeMap( static_cast< int >(MSI::Dof_Type::TEMP) ) = 1;
        tIWG->mSet->mLeaderDofTypeMap( static_cast< int >(MSI::Dof_Type::LS1) ) = 2;
        tIWG->mSet->mLeaderDofTypeMap( static_cast< int >(MSI::Dof_Type::VX) ) = 3;

//        tIWG->mSet->mFollowerDofTypeMap.set_size( static_cast< int >(MSI::Dof_Type::UNDEFINED) + 1, 1, -1 );
//        tIWG->mSet->mFollowerDofTypeMap( static_cast< int >(MSI::Dof_Type::UX) ) = 0;
//        tIWG->mSet->mFollowerDofTypeMap( static_cast< int >(MSI::Dof_Type::TEMP) ) = 1;
//        tIWG->mSet->mFollowerDofTypeMap( static_cast< int >(MSI::Dof_Type::LS1) ) = 2;
//        tIWG->mSet->mFollowerDofTypeMap( static_cast< int >(MSI::Dof_Type::VX) ) = 3;

        tIWG->set_set_pointer(static_cast<fem::Set*>(tSet));

        moris::Vector< moris::Vector< enum MSI::Dof_Type > > tDummy;
        Field_Interpolator_Manager tFIManager( tDummy, tSet );
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
        CHECK( equal_to( tIWG->mLeaderGlobalDofTypes.size(), 3 ));

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

        delete (tSet);

    }/* TEST_CASE */

    }/* namespace fem */
}/* namespace moris */

