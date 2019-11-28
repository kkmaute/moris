
#include "catch.hpp"
#include "fn_equal_to.hpp"
#include "cl_FEM_Property.hpp"              //FEM/INT/src
#include "cl_FEM_Geometry_Interpolator.hpp" //FEM/INT/src
#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_CM_Factory.hpp"            //FEM/INT/src

#define protected public
#define private   public
#include "cl_FEM_Constitutive_Model.hpp" //FEM/INT/src
#undef protected
#undef private

moris::Matrix< moris::DDRMat > tValFunctionCM( moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                                               moris::Cell< moris::fem::Field_Interpolator* > & aDofFI,
                                               moris::Cell< moris::fem::Field_Interpolator* > & aDvFI,
                                               moris::fem::Geometry_Interpolator              * aGeometryInterpolator )
{
    return aParameters( 0 ) + aParameters( 1 ) * aDofFI( 0 )->val();
}

moris::Matrix< moris::DDRMat > tDerFunctionCM( moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                                               moris::Cell< moris::fem::Field_Interpolator* > & aDofFI,
                                               moris::Cell< moris::fem::Field_Interpolator* > & aDvFI,
                                               moris::fem::Geometry_Interpolator              * aGeometryInterpolator )
{
    return aParameters( 1 ) * aDofFI( 0 )->N();
}

namespace moris
{
    namespace fem
    {
        TEST_CASE( "Constitutive_Model", "[moris],[fem],[CM]" )
        {
            // real for check
            real tEpsilon = 1E-6;

            // create the properties
            std::shared_ptr< fem::Property > tProp = std::make_shared< fem::Property > ();
            tProp->set_parameters( {{{ 1.0}}, {{1.0 }}} );
            tProp->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tProp->set_val_function( tValFunctionCM );
            tProp->set_dof_derivative_functions( { tDerFunctionCM } );

            // create a constitutive model
            Constitutive_Model tCM;

            // set constitutive type
            tCM.set_constitutive_type( fem::Constitutive_Type::DIFF_LIN_ISO );

            // set space dim
            tCM.set_space_dim( 2 );

            // set dof types
            tCM.set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );

            // set dv types
            tCM.set_dv_type_list( {{ MSI::Dv_Type::LS1 }} );

            // set property
            tCM.set_properties( { tProp } );

            //create a space and a time geometry interpolator
            Geometry_Interpolator tGI;

            // create a dof field interpolator
            Cell< Field_Interpolator* > tDofFIs( 1, nullptr );
            tDofFIs( 0 ) = new Field_Interpolator ( 1, { MSI::Dof_Type::TEMP } );

            // create a dv field interpolator
            Cell< Field_Interpolator* > tDvFIs( 2, nullptr );
            tDvFIs( 0 ) = new Field_Interpolator ( 1, { MSI::Dv_Type::LS1 } );
            tDvFIs( 1 ) = new Field_Interpolator ( 1, { MSI::Dv_Type::LS2 } );

            // set FI for constitutive model
            tCM.set_dof_field_interpolators( tDofFIs );

            // set geo FI for constitutive model
            tCM.set_geometry_interpolator( &tGI );

            // check dof field interpolator
            //------------------------------------------------------------------------------
            tCM.check_dof_field_interpolators();

            // dof map
            //------------------------------------------------------------------------------
            // check global dof map size
            CHECK( equal_to( tCM.get_dof_type_map().n_cols(), 1 ));
            CHECK( equal_to( tCM.get_dof_type_map().n_rows(), 4 ));

            // check global dof map content
            CHECK( equal_to( tCM.get_dof_type_map()( 0, 0 ), -1 ));
            CHECK( equal_to( tCM.get_dof_type_map()( 1, 0 ), -1 ));
            CHECK( equal_to( tCM.get_dof_type_map()( 2, 0 ), -1 ));
            CHECK( equal_to( tCM.get_dof_type_map()( 3, 0 ), 0 ));

            // global dof list and map
            //------------------------------------------------------------------------------
            // check global dof list size
            CHECK( equal_to( tCM.get_global_dof_type_list().size(), 1 ));

            // check global dof list content
            CHECK( equal_to( static_cast< uint >( tCM.get_global_dof_type_list()( 0 )( 0 ) ), 3 ) ); //TEMP

            // check global dof map size
            CHECK( equal_to( tCM.get_global_dof_type_map().n_cols(), 1 ));
            CHECK( equal_to( tCM.get_global_dof_type_map().n_rows(), 4 ));

            // check global dof map content
            CHECK( equal_to( tCM.get_global_dof_type_map()( 0, 0 ), -1 ));
            CHECK( equal_to( tCM.get_global_dof_type_map()( 1, 0 ), -1 ));
            CHECK( equal_to( tCM.get_global_dof_type_map()( 2, 0 ), -1 ));
            CHECK( equal_to( tCM.get_global_dof_type_map()( 3, 0 ), 0 ));

            // dof dependency
            //------------------------------------------------------------------------------
            // check dependency on TEMP
            CHECK( tCM.check_dof_dependency({ MSI::Dof_Type::TEMP }) );

            // check dependency on UX
            CHECK( !tCM.check_dof_dependency({ MSI::Dof_Type::UX }) );

            // clean up
            //------------------------------------------------------------------------------

            // delete the field interpolator pointers
            for( Field_Interpolator* tFI : tDofFIs )
            {
                delete tFI;
            }
            tDofFIs.clear();

        }/* TEST_CASE */


    }/* namespace fem */
}/* namespace moris */
