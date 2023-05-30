/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_CM.cpp
 *
 */

#include "catch.hpp"
#include "fn_equal_to.hpp"
#include "cl_FEM_Property.hpp"              //FEM/INT/src
#include "cl_FEM_Geometry_Interpolator.hpp" //FEM/INT/src
#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_CM_Factory.hpp"            //FEM/INT/src

#define protected public
#define private   public
#include "cl_FEM_Constitutive_Model.hpp" //FEM/INT/src
#include "cl_FEM_Field_Interpolator_Manager.hpp"                   //FEM//INT//src
#include "cl_FEM_Set.hpp"                   //FEM//INT//src
#undef protected
#undef private

void tValFunctionCM
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 ) + aParameters( 1 ) * aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::TEMP )->val();
}

void tDerFunctionCM
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 1 ) * aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::TEMP )->N();
}

namespace moris
{
    namespace fem
    {
        TEST_CASE( "Constitutive_Model", "[moris],[fem],[CM]" )
        {

            // create the properties
            std::shared_ptr< fem::Property > tProp = std::make_shared< fem::Property > ();
            tProp->set_parameters( {{{ 1.0}}, {{1.0 }}} );
            tProp->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tProp->set_val_function( tValFunctionCM );
            tProp->set_dof_derivative_functions( { tDerFunctionCM } );

            // create a constitutive model
            CM_Factory tCMFactory;
            std::shared_ptr< fem::Constitutive_Model > tCM = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );

            // set space dim
            tCM->set_space_dim( 2 );

            // set dof types
            tCM->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );

            // set dv types
            tCM->set_dv_type_list( {{ PDV_Type::LS1 }} );

            // set property
            tCM->set_property( tProp, "Conductivity" );

            //create a space and a time geometry interpolator
            Geometry_Interpolator tGI;

            // create a dof field interpolator
            Cell< Field_Interpolator* > tDofFIs( 1, nullptr );
            tDofFIs( 0 ) = new Field_Interpolator ( 1, { MSI::Dof_Type::TEMP } );

            // create a dv field interpolator
            Cell< Field_Interpolator* > tDvFIs( 2, nullptr );
            tDvFIs( 0 ) = new Field_Interpolator ( 1, { PDV_Type::LS1 } );
            tDvFIs( 1 ) = new Field_Interpolator ( 1, { PDV_Type::LS2 } );

            // set a fem set pointer
            MSI::Equation_Set * tSet = new fem::Set();

            // set size and populate the set dof type map
            reinterpret_cast<fem::Set*>(tSet)->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
            reinterpret_cast<fem::Set*>(tSet)->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 0;

            // set size and populate the set leader dof type map
            tSet->mLeaderDofTypeMap.set_size( static_cast< int >(MSI::Dof_Type::END_ENUM) + 1, 1, -1 );
            tSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 0;

            // create a field interpolator manager
            moris::Cell< moris::Cell< enum MSI::Dof_Type > > tDummy;
            Field_Interpolator_Manager tFIManager( tDummy, tSet );

            // populate the field interpolator manager
            tFIManager.mFI = tDofFIs;
            tFIManager.mDvFI = tDvFIs;
            tFIManager.mIPGeometryInterpolator = &tGI;
            tFIManager.mIGGeometryInterpolator = &tGI;

            // set IWG field interpolator manager
            tCM->set_field_interpolator_manager( &tFIManager );

            // dof map
            //------------------------------------------------------------------------------
            // check global dof map size
            CHECK( equal_to( tCM->get_dof_type_map().n_cols(), 1 ));
            CHECK( equal_to( tCM->get_dof_type_map().n_rows(), 4 ));

            // check global dof map content
            CHECK( equal_to( tCM->get_dof_type_map()( 0, 0 ), -1 ));
            CHECK( equal_to( tCM->get_dof_type_map()( 1, 0 ), -1 ));
            CHECK( equal_to( tCM->get_dof_type_map()( 2, 0 ), -1 ));
            CHECK( equal_to( tCM->get_dof_type_map()( 3, 0 ), 0 ));

            // global dof list and map
            //------------------------------------------------------------------------------
            // check global dof list size
            CHECK( equal_to( tCM->get_global_dof_type_list().size(), 1 ));

            // check global dof list content
            CHECK( equal_to( static_cast< uint >( tCM->get_global_dof_type_list()( 0 )( 0 ) ), 3 ) ); //TEMP

            // check global dof map size
            CHECK( equal_to( tCM->get_global_dof_type_map().n_cols(), 1 ));
            CHECK( equal_to( tCM->get_global_dof_type_map().n_rows(), 4 ));

            // check global dof map content
            CHECK( equal_to( tCM->get_global_dof_type_map()( 0, 0 ), -1 ));
            CHECK( equal_to( tCM->get_global_dof_type_map()( 1, 0 ), -1 ));
            CHECK( equal_to( tCM->get_global_dof_type_map()( 2, 0 ), -1 ));
            CHECK( equal_to( tCM->get_global_dof_type_map()( 3, 0 ), 0 ));

            // dof dependency
            //------------------------------------------------------------------------------
            // check dependency on TEMP
            CHECK( tCM->check_dof_dependency({ MSI::Dof_Type::TEMP }) );

            // check dependency on UX
            CHECK( !tCM->check_dof_dependency({ MSI::Dof_Type::UX }) );

            // clean up
            //------------------------------------------------------------------------------

            // delete the field interpolator pointers
            tDofFIs.clear();

        }/* TEST_CASE */

    }/* namespace fem */
}/* namespace moris */

