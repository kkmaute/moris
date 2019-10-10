
#include "catch.hpp"
#include "fn_equal_to.hpp"
#include "cl_FEM_Property.hpp"              //FEM/INT/src
#include "cl_FEM_Geometry_Interpolator.hpp"         //FEM/INT/src
#include "cl_FEM_Field_Interpolator.hpp"         //FEM/INT/src
#include "cl_FEM_CM_Factory.hpp"         //FEM/INT/src
#include "cl_FEM_Constitutive_Model.hpp"         //FEM/INT/src
#include "cl_FEM_CM_Diffusion_Linear_Isotropic.hpp"         //FEM/INT/src

#define protected public
#define private   public
#include "cl_FEM_IWG.hpp"         //FEM/INT/src
#undef protected
#undef private

#include "cl_FEM_IWG_Factory.hpp"         //FEM/INT/src

moris::Matrix< moris::DDRMat > tValFunction_UTIWG( moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                                                   moris::Cell< moris::fem::Field_Interpolator* > & aDofFI,
                                                   moris::Cell< moris::fem::Field_Interpolator* > & aDvFI,
                                                   moris::fem::Geometry_Interpolator              * aGeometryInterpolator )
{
    moris::Matrix< moris::DDRMat > tPropertyVal( 1, 1, 1.0);
    return tPropertyVal;
}

namespace moris
{
    namespace fem
    {

    // This test case tests all the member functions of the IWG,
    // except the residual and jacobian related ones.
    TEST_CASE( "IWG_dof_dv", "[moris],[fem],[IWG_dof_dv]" )
    {
        // create a geometry interpolator
        fem::Geometry_Interpolator tGI;

        // create IWG master properties
        Cell< fem::Property* > tMasterIWGProp( 2, nullptr );
        tMasterIWGProp( 0 ) = new Property( fem::Property_Type::TEMP_NEUMANN,
                                         {{ MSI::Dof_Type::TEMP }, { MSI::Dof_Type::VX }},
                                         Cell< Matrix< DDRMat > >( 0 ),
                                         tValFunction_UTIWG,
                                         Cell< fem::PropertyFunc >( 0 ),
                                         & tGI );

        tMasterIWGProp( 1 ) = new Property( fem::Property_Type::CONDUCTIVITY,
                                         {{ MSI::Dof_Type::TEMP }, { MSI::Dof_Type::LS1 }},
                                         {{ MSI::Dv_Type::LS1 }, { MSI::Dv_Type::LS2 }},
                                         Cell< Matrix< DDRMat > >( 0 ),
                                         tValFunction_UTIWG,
                                         Cell< fem::PropertyFunc >( 0 ),
                                         Cell< fem::PropertyFunc >( 0 ),
                                         & tGI );
        // create CM properties
        Cell< fem::Property* > tMasterCMProp( 1, nullptr );
        tMasterCMProp( 0 ) = tMasterIWGProp( 1 );

        // create IWG slave properties
        Cell< fem::Property* > tSlaveIWGProp( 2, nullptr );
        tSlaveIWGProp( 0 ) = new Property( fem::Property_Type::TEMP_LOAD,
                                      Cell< Cell< MSI::Dof_Type > >( 0 ),
                                      Cell< Matrix< DDRMat > >( 0 ),
                                      tValFunction_UTIWG,
                                      Cell< fem::PropertyFunc >( 0 ),
                                      & tGI );

        tSlaveIWGProp( 1 ) = new Property( fem::Property_Type::CONDUCTIVITY,
                                      {{ MSI::Dof_Type::TEMP }, { MSI::Dof_Type::LS1 }},
                                      {{ MSI::Dv_Type::LS1 }, { MSI::Dv_Type::LS2 }},
                                      Cell< Matrix< DDRMat > >( 0 ),
                                      tValFunction_UTIWG,
                                      Cell< fem::PropertyFunc >( 0 ),
                                      Cell< fem::PropertyFunc >( 0 ),
                                      & tGI );
        // create CM properties
        Cell< fem::Property* > tSlaveCMProp( 1, nullptr );
        tSlaveCMProp( 0 ) = tSlaveIWGProp( 1 );

        // create master constitutive model
        Cell< fem::Constitutive_Model* > tMasterCMs( 1 );
        tMasterCMs( 0 ) = new CM_Diffusion_Linear_Isotropic();
        tMasterCMs( 0 )->set_space_dim( 3 );
        tMasterCMs( 0 )->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tMasterCMs( 0 )->set_property_type_list( { fem::Property_Type::CONDUCTIVITY } );
        tMasterCMs( 0 )->set_properties( tMasterCMProp );

        // create slave constitutive model
        Cell< fem::Constitutive_Model* > tSlaveCMs( 1 );
        tSlaveCMs( 0 ) = new CM_Diffusion_Linear_Isotropic();
        tSlaveCMs( 0 )->set_space_dim( 3 );
        tSlaveCMs( 0 )->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tSlaveCMs( 0 )->set_property_type_list( { fem::Property_Type::CONDUCTIVITY } );
        tSlaveCMs( 0 )->set_properties( tSlaveCMProp );

        // create master dof field interpolators
        uint tNumberOfFields = 1;
        Cell< Field_Interpolator* > tMasterDofFI( 3, nullptr );
        tMasterDofFI( 0 ) = new Field_Interpolator ( tNumberOfFields, { MSI::Dof_Type::TEMP } );
        tMasterDofFI( 1 ) = new Field_Interpolator ( tNumberOfFields, { MSI::Dof_Type::VX } );
        tMasterDofFI( 2 ) = new Field_Interpolator ( tNumberOfFields, { MSI::Dof_Type::LS1 } );

        // create master dv field interpolators
        Cell< Field_Interpolator* > tMasterDvFI( 3, nullptr );
        tMasterDvFI( 0 ) = new Field_Interpolator ( tNumberOfFields, { MSI::Dv_Type::DENSITY } );
        tMasterDvFI( 1 ) = new Field_Interpolator ( tNumberOfFields, { MSI::Dv_Type::LS1 } );
        tMasterDvFI( 2 ) = new Field_Interpolator ( tNumberOfFields, { MSI::Dv_Type::LS2 } );

        // create slave dof field interpolators
        Cell< Field_Interpolator* > tSlaveDofFI( 3, nullptr );
        tSlaveDofFI( 0 ) = new Field_Interpolator ( tNumberOfFields, { MSI::Dof_Type::TEMP } );
        tSlaveDofFI( 1 ) = new Field_Interpolator ( tNumberOfFields, { MSI::Dof_Type::UX } );
        tSlaveDofFI( 2 ) = new Field_Interpolator ( tNumberOfFields, { MSI::Dof_Type::LS1 } );

        // create slave dv field interpolators
        Cell< Field_Interpolator* > tSlaveDvFI( 2, nullptr );
        tSlaveDvFI( 0 ) = new Field_Interpolator ( tNumberOfFields, { MSI::Dv_Type::LS1 } );
        tSlaveDvFI( 1 ) = new Field_Interpolator ( tNumberOfFields, { MSI::Dv_Type::LS2 } );

        // create an IWG
        fem::IWG tIWG;

        // set residual dof type
        tIWG.set_residual_dof_type( { MSI::Dof_Type::TEMP } );

        // set active master and slave dof types
        tIWG.set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tIWG.set_dof_type_list( {{ MSI::Dof_Type::TEMP },{ MSI::Dof_Type::UX }}, mtk::Master_Slave::SLAVE );

        // set active master and slave dv types
        tIWG.set_dv_type_list( {{ MSI::Dv_Type::DENSITY }, {MSI::Dv_Type::LS1 }} );
        tIWG.set_dv_type_list( {{MSI::Dv_Type::LS1 }}, mtk::Master_Slave::SLAVE );

        // set active master and slave constitutive types
        tIWG.set_constitutive_type_list( { fem::Constitutive_Type::DIFF_LIN_ISO } );
        tIWG.set_constitutive_type_list( { fem::Constitutive_Type::DIFF_LIN_ISO }, mtk::Master_Slave::SLAVE );

        // set active master and slave property type
        tIWG.set_property_type_list( { fem::Property_Type::TEMP_NEUMANN } );
        tIWG.set_property_type_list( { fem::Property_Type::TEMP_LOAD }, mtk::Master_Slave::SLAVE );

        // set master and slave constitutive models
        tIWG.set_constitutive_models( tMasterCMs );
        tIWG.set_constitutive_models( tSlaveCMs, mtk::Master_Slave::SLAVE );

        // set master and slave properties
        tIWG.set_properties( tMasterIWGProp );
        tIWG.set_properties( tSlaveIWGProp, mtk::Master_Slave::SLAVE );

        // set IWG master and slave dof field interpolators
        tIWG.set_dof_field_interpolators( tMasterDofFI );
        tIWG.set_dof_field_interpolators( tSlaveDofFI, mtk::Master_Slave::SLAVE );

        // set IWG master and slave dv field interpolators
        tIWG.set_dv_field_interpolators( tMasterDvFI );
        tIWG.set_dv_field_interpolators( tSlaveDvFI, mtk::Master_Slave::SLAVE );

        // constitutive models check-----------------------------------------------------
        // check constitutive model pointers
        tIWG.check_constitutive_models();
        tIWG.check_constitutive_models( mtk::Master_Slave::SLAVE );

        // properties check--------------------------------------------------------------
        // check master global property list size
        CHECK( equal_to( tIWG.mMasterGlobalPropTypes.size(), 2 ));

        //check master global property list content
        CHECK( equal_to( static_cast< uint >( tIWG.mMasterGlobalPropTypes( 0 ) ), 2 ) );
        CHECK( equal_to( static_cast< uint >( tIWG.mMasterGlobalPropTypes( 1 ) ), 3 ) );

        // check slave global property list size
        CHECK( equal_to( tIWG.mSlaveGlobalPropTypes.size(), 2 ));

        //check slave global property list content
        CHECK( equal_to( static_cast< uint >( tIWG.mSlaveGlobalPropTypes( 0 ) ), 4 ) );
        CHECK( equal_to( static_cast< uint >( tIWG.mSlaveGlobalPropTypes( 1 ) ), 3 ) );

        // check property pointers
        tIWG.check_properties();
        tIWG.check_properties( mtk::Master_Slave::SLAVE );

        // dof check--------------------------------------------------------------------
        // check master global dof list size
        CHECK( equal_to( tIWG.mMasterGlobalDofTypes.size(), 3 ));

        // check master global dof list content
        CHECK( equal_to( static_cast< uint >( tIWG.mMasterGlobalDofTypes( 0 )( 0 ) ), 3 ) );
        CHECK( equal_to( static_cast< uint >( tIWG.mMasterGlobalDofTypes( 1 )( 0 ) ), 11 ) );
        CHECK( equal_to( static_cast< uint >( tIWG.mMasterGlobalDofTypes( 2 )( 0 ) ), 6 ) );

        // check slave global dof list size
        CHECK( equal_to( tIWG.mSlaveGlobalDofTypes.size(), 3 ));

        // check slave global dof list content
        CHECK( equal_to( static_cast< uint >( tIWG.mSlaveGlobalDofTypes( 0 )( 0 ) ), 3 ) );
        CHECK( equal_to( static_cast< uint >( tIWG.mSlaveGlobalDofTypes( 1 )( 0 ) ), 0 ) );
        CHECK( equal_to( static_cast< uint >( tIWG.mSlaveGlobalDofTypes( 2 )( 0 ) ), 6 ) );

        // check dof field interpolators
        tIWG.check_dof_field_interpolators();
        tIWG.check_dof_field_interpolators( mtk::Master_Slave::SLAVE );

        // dv check---------------------------------------------------------------------
        // check master global dv list size
        CHECK( equal_to( tIWG.mMasterGlobalDvTypes.size(), 3 ));

        // check master global dv list content
        CHECK( equal_to( static_cast< uint >( tIWG.mMasterGlobalDvTypes( 0 )( 0 ) ), 2 ) );
        CHECK( equal_to( static_cast< uint >( tIWG.mMasterGlobalDvTypes( 1 )( 0 ) ), 0 ) );
        CHECK( equal_to( static_cast< uint >( tIWG.mMasterGlobalDvTypes( 2 )( 0 ) ), 1 ) );

        // check slave global dv list size
        CHECK( equal_to( tIWG.mSlaveGlobalDvTypes.size(), 2 ));

        // check master global dv list content
        CHECK( equal_to( static_cast< uint >( tIWG.mSlaveGlobalDvTypes( 0 )( 0 ) ), 0 ) );
        CHECK( equal_to( static_cast< uint >( tIWG.mSlaveGlobalDvTypes( 1 )( 0 ) ), 1 ) );

        // check dv field interpolators
        tIWG.check_dv_field_interpolators();
        tIWG.check_dv_field_interpolators( mtk::Master_Slave::SLAVE );

        // clean up---------------------------------------------------------------------
        // delete the property pointers
        for( Property* tProp : tMasterIWGProp )
        {
            delete tProp;
        }
        tMasterIWGProp.clear();

        // delete the property pointers
        for( Property* tProp : tSlaveIWGProp )
        {
            delete tProp;
        }
        tSlaveIWGProp.clear();

        // delete the dof field interpolator pointers
        for( Field_Interpolator* tFI : tMasterDofFI )
        {
            delete tFI;
        }
        tMasterDofFI.clear();

        // delete the dof field interpolator pointers
        for( Field_Interpolator* tFI : tSlaveDofFI )
        {
            delete tFI;
        }
        tSlaveDofFI.clear();

        // delete the dv field interpolator pointers
        for( Field_Interpolator* tFI : tMasterDvFI )
        {
            delete tFI;
        }
        tMasterDvFI.clear();

        // delete the dv field interpolator pointers
        for( Field_Interpolator* tFI : tSlaveDvFI )
        {
            delete tFI;
        }
        tSlaveDvFI.clear();

        // delete the constitutive model pointers
        for( Constitutive_Model* tCM : tMasterCMs )
        {
            delete tCM;
        }
        tMasterCMs.clear();

        // delete the constitutive model pointers
        for( Constitutive_Model* tCM : tSlaveCMs )
        {
            delete tCM;
        }
        tSlaveCMs.clear();

    }/* TEST_CASE */

    }/* namespace fem */
}/* namespace moris */
