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
        // create the properties
        std::shared_ptr< fem::Property > tPropMaster1 = std::make_shared< fem::Property > ();
        tPropMaster1->set_dof_type_list( {{ MSI::Dof_Type::TEMP }, { MSI::Dof_Type::VX }} );
        tPropMaster1->set_val_function( tValFunction_UTIWG );

        std::shared_ptr< fem::Property > tPropMaster2 = std::make_shared< fem::Property > ();
        tPropMaster2->set_dof_type_list( {{ MSI::Dof_Type::TEMP }, { MSI::Dof_Type::LS1 }} );
        tPropMaster2->set_dv_type_list( {{ MSI::Dv_Type::LS1 }, { MSI::Dv_Type::LS2 }} );
        tPropMaster2->set_val_function( tValFunction_UTIWG );

        std::shared_ptr< fem::Property > tPropSlave1 = std::make_shared< fem::Property > ();
        tPropSlave1->set_val_function( tValFunction_UTIWG );

        std::shared_ptr< fem::Property > tPropSlave2 = std::make_shared< fem::Property > ();
        tPropSlave2->set_dof_type_list( {{ MSI::Dof_Type::TEMP }, { MSI::Dof_Type::LS1 }} );
        tPropSlave2->set_dv_type_list( {{ MSI::Dv_Type::LS1 }, { MSI::Dv_Type::LS2 }} );
        tPropSlave2->set_val_function( tValFunction_UTIWG );

        // define constitutive models
        fem::CM_Factory tCMFactory;

        std::shared_ptr< fem::Constitutive_Model > tCMMaster1 = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
        tCMMaster1->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tCMMaster1->set_properties( { tPropMaster2 } );
        tCMMaster1->set_space_dim( 3 );

        std::shared_ptr< fem::Constitutive_Model > tCMSlave1 = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
        tCMSlave1->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tCMSlave1->set_properties( { tPropSlave2 } );
        tCMSlave1->set_space_dim( 3 );

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

        // define the IWGs
        fem::IWG_Factory tIWGFactory;

        std::shared_ptr< moris::fem::IWG > tIWG = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
        tIWG->set_residual_dof_type( { MSI::Dof_Type::TEMP } );
        tIWG->set_dof_type_list( {{ MSI::Dof_Type::TEMP }}, mtk::Master_Slave::MASTER );
        tIWG->set_dof_type_list( {{ MSI::Dof_Type::TEMP },{ MSI::Dof_Type::UX }}, mtk::Master_Slave::SLAVE );
        tIWG->set_dv_type_list( {{ MSI::Dv_Type::DENSITY }, {MSI::Dv_Type::LS1 }}, mtk::Master_Slave::MASTER );
        tIWG->set_dv_type_list( {{MSI::Dv_Type::LS1 }}, mtk::Master_Slave::SLAVE );
        tIWG->set_constitutive_models( { tCMMaster1 }, mtk::Master_Slave::MASTER );
        tIWG->set_constitutive_models( { tCMSlave1 }, mtk::Master_Slave::SLAVE );
        tIWG->set_properties( { tPropMaster1 }, mtk::Master_Slave::MASTER );
        tIWG->set_properties( { tPropSlave1 }, mtk::Master_Slave::SLAVE );

        // build master and slave global dof type list
        tIWG->build_global_dof_type_list();

        // set IWG master and slave dof field interpolators
        tIWG->set_dof_field_interpolators( tMasterDofFI );
        tIWG->set_dof_field_interpolators( tSlaveDofFI, mtk::Master_Slave::SLAVE );

        // build master and slave global dv type list
        tIWG->build_global_dv_type_list();

        // set IWG master and slave dv field interpolators
        tIWG->set_dv_field_interpolators( tMasterDvFI );
        tIWG->set_dv_field_interpolators( tSlaveDvFI, mtk::Master_Slave::SLAVE );

        // dof check--------------------------------------------------------------------
        // check master global dof list size
        CHECK( equal_to( tIWG->mMasterGlobalDofTypes.size(), 3 ));

        // check master global dof list content
        CHECK( equal_to( static_cast< uint >( tIWG->mMasterGlobalDofTypes( 0 )( 0 ) ), 3 ) );
        CHECK( equal_to( static_cast< uint >( tIWG->mMasterGlobalDofTypes( 1 )( 0 ) ), 11 ) );
        CHECK( equal_to( static_cast< uint >( tIWG->mMasterGlobalDofTypes( 2 )( 0 ) ), 6 ) );

        // check slave global dof list size
        CHECK( equal_to( tIWG->mSlaveGlobalDofTypes.size(), 3 ));

        // check slave global dof list content
        CHECK( equal_to( static_cast< uint >( tIWG->mSlaveGlobalDofTypes( 0 )( 0 ) ), 3 ) );
        CHECK( equal_to( static_cast< uint >( tIWG->mSlaveGlobalDofTypes( 1 )( 0 ) ), 0 ) );
        CHECK( equal_to( static_cast< uint >( tIWG->mSlaveGlobalDofTypes( 2 )( 0 ) ), 6 ) );

        // check dof field interpolators
        tIWG->check_dof_field_interpolators();
        tIWG->check_dof_field_interpolators( mtk::Master_Slave::SLAVE );

        // dv check---------------------------------------------------------------------
        // check master global dv list size
        CHECK( equal_to( tIWG->mMasterGlobalDvTypes.size(), 3 ));

        // check master global dv list content
        CHECK( equal_to( static_cast< uint >( tIWG->mMasterGlobalDvTypes( 0 )( 0 ) ), 2 ) );
        CHECK( equal_to( static_cast< uint >( tIWG->mMasterGlobalDvTypes( 1 )( 0 ) ), 0 ) );
        CHECK( equal_to( static_cast< uint >( tIWG->mMasterGlobalDvTypes( 2 )( 0 ) ), 1 ) );

        // check slave global dv list size
        CHECK( equal_to( tIWG->mSlaveGlobalDvTypes.size(), 2 ));

        // check master global dv list content
        CHECK( equal_to( static_cast< uint >( tIWG->mSlaveGlobalDvTypes( 0 )( 0 ) ), 0 ) );
        CHECK( equal_to( static_cast< uint >( tIWG->mSlaveGlobalDvTypes( 1 )( 0 ) ), 1 ) );

        // check dv field interpolators
        tIWG->check_dv_field_interpolators();
        tIWG->check_dv_field_interpolators( mtk::Master_Slave::SLAVE );

        // clean up---------------------------------------------------------------------
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

    }/* TEST_CASE */

    }/* namespace fem */
}/* namespace moris */
