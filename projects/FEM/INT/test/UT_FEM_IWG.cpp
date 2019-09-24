
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
                                                   moris::Cell< moris::fem::Field_Interpolator* > & aFieldInterpolator,
                                                   moris::fem::Geometry_Interpolator              * aGeometryInterpolator )
{
    moris::Matrix< moris::DDRMat > tPropertyVal( 1, 1, 1.0);
    return tPropertyVal;
}

namespace moris
{
    namespace fem
    {

        TEST_CASE( "Property_with_IWG", "[moris],[fem],[Property_with_IWG]" )
        {
            //create a quad4 space element
            Matrix< DDRMat > tXHat( 4, 2 );
            tXHat( 0, 0 ) = 0.0; tXHat( 0, 1 ) = 0.0;
            tXHat( 1, 0 ) = 3.0; tXHat( 1, 1 ) = 1.25;
            tXHat( 2, 0 ) = 4.5; tXHat( 2, 1 ) = 4.0;
            tXHat( 3, 0 ) = 1.0; tXHat( 3, 1 ) = 3.25;

            //create a line time element
            Matrix< DDRMat > tTHat( 2, 1 );
            tTHat( 0 ) = 0.0;
            tTHat( 1 ) = 5.0;

            //create a space geometry interpolation rule
            Interpolation_Rule tGeomInterpRule( mtk::Geometry_Type::QUAD,
                                                Interpolation_Type::LAGRANGE,
                                                mtk::Interpolation_Order::LINEAR,
                                                Interpolation_Type::LAGRANGE,
                                                mtk::Interpolation_Order::LINEAR );

            //create a space and a time geometry interpolator
            Geometry_Interpolator* tGeomInterpolator = new Geometry_Interpolator( tGeomInterpRule );

            //set the coefficients xHat, tHat
            tGeomInterpolator->set_coeff( tXHat, tTHat );

            // create a property object
            Cell< fem::Property* > tProperties( 2, nullptr );
            tProperties( 0 ) = new Property( fem::Property_Type::TEMP_NEUMANN,
                                             {{ MSI::Dof_Type::TEMP }, { MSI::Dof_Type::VX }},
                                             Cell< Matrix< DDRMat > >( 0 ),
                                             tValFunction_UTIWG,
                                             Cell< fem::PropertyFunc >( 0 ),
                                             tGeomInterpolator );

            tProperties( 1 ) = new Property( fem::Property_Type::CONDUCTIVITY,
                                             {{ MSI::Dof_Type::TEMP }, { MSI::Dof_Type::LS1 }},
                                             Cell< Matrix< DDRMat > >( 0 ),
                                             tValFunction_UTIWG,
                                             Cell< fem::PropertyFunc >( 0 ),
                                             tGeomInterpolator );

            // create a constitutive model
            CM_Factory tCMFactory;

            // create a constitutive model
            Cell< fem::Constitutive_Model* > tCMs( 1 );
            tCMs( 0 ) = tCMFactory.create_CM( Constitutive_Type::DIFF_LIN_ISO );

            // set space dim, dof types, property type
            tCMs( 0 )->set_space_dim( 3 );
            tCMs( 0 )->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tCMs( 0 )->set_property_type_list( { fem::Property_Type::CONDUCTIVITY } );
            Cell< fem::Property* > tCMProp( 1, tProperties( 1 ) );
            tCMs( 0 )->set_properties( tCMProp );

            // create an interpolation rule
            Interpolation_Rule tInterpolationRule ( mtk::Geometry_Type::QUAD,
                                                    Interpolation_Type::LAGRANGE,
                                                    mtk::Interpolation_Order::LINEAR,
                                                    Interpolation_Type::LAGRANGE,
                                                    mtk::Interpolation_Order::LINEAR );

            // create a TEMP and a VX field interpolator
            uint tNumberOfFields = 1;
            Cell< Field_Interpolator* > tFieldInterpolators( 3, nullptr );
            tFieldInterpolators( 0 ) = new Field_Interpolator ( tNumberOfFields,
                                                               tInterpolationRule,
                                                               tGeomInterpolator,
                                                               { MSI::Dof_Type::TEMP } );
            tFieldInterpolators( 1 ) = new Field_Interpolator ( tNumberOfFields,
                                                                tInterpolationRule,
                                                                tGeomInterpolator,
                                                                { MSI::Dof_Type::VX } );
            tFieldInterpolators( 2 ) = new Field_Interpolator ( tNumberOfFields,
                                                                tInterpolationRule,
                                                                tGeomInterpolator,
                                                                { MSI::Dof_Type::LS1 } );

            // set coefficients for field interpolators
            Matrix< DDRMat > tUHat0( 8, 1, 2.0 );
            tFieldInterpolators( 0 )->set_coeff( tUHat0 );
            Matrix< DDRMat > tUHat1( 8, 1, 3.0 );
            tFieldInterpolators( 1 )->set_coeff( tUHat1 );

            // a factory to create the IWGs
            fem::IWG_Factory tIWGFactory;

            // create an IWG with the factory for the ith IWG type
            fem::IWG* tIWG = tIWGFactory.create_IWGs( fem::IWG_Type::SPATIALDIFF_NEUMANN );

            // set space dim
            tIWG->set_space_dim( 2 );

            // set residual dof type
            tIWG->set_residual_dof_type( { MSI::Dof_Type::TEMP } );

            // set active dof types
            tIWG->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );

            // set active constitutive types
            tIWG->set_constitutive_type_list( { fem::Constitutive_Type::DIFF_LIN_ISO } );

            // set active property type
            tIWG->set_property_type_list( { fem::Property_Type::TEMP_NEUMANN } );

            // set constitutive models
            tIWG->set_constitutive_models( tCMs );

            // set IWG properties
            tIWG->set_properties( tProperties );

            // set IWG field interpolators
            tIWG->set_field_interpolators( tFieldInterpolators );

            //mMasterGlobalPropTypes check--------------------------------------------------
            //check global property list size
            CHECK( equal_to( tIWG->mMasterGlobalPropTypes.size(), 2 ));

            //check global dof list content
            CHECK( equal_to( static_cast< uint >( tIWG->mMasterGlobalPropTypes( 0 ) ), 2 ) );
            CHECK( equal_to( static_cast< uint >( tIWG->mMasterGlobalPropTypes( 1 ) ), 3 ) );

            //mMasterGlobalDofTypes check---------------------------------------------------
            //check global dof list size
            CHECK( equal_to( tIWG->mMasterGlobalDofTypes.size(), 3 ));

            //check global dof list content
            CHECK( equal_to( static_cast< uint >( tIWG->mMasterGlobalDofTypes( 0 )( 0 ) ), 3 ) );
            CHECK( equal_to( static_cast< uint >( tIWG->mMasterGlobalDofTypes( 1 )( 0 ) ), 11 ) );
            CHECK( equal_to( static_cast< uint >( tIWG->mMasterGlobalDofTypes( 2 )( 0 ) ), 6 ) );

            // clean up

            // delete the property pointers
            for( uint iProp = 0; iProp < tProperties.size(); iProp++ )
            {
                delete tProperties( iProp );
            }

            // delete the field interpolator pointers
            for( uint iFI = 0; iFI < tFieldInterpolators.size(); iFI++ )
            {
                delete tFieldInterpolators( iFI );
            }

            // delete the constitutive model pointers
            for( uint iCM = 0; iCM < tCMs.size(); iCM++ )
            {
                delete tCMs( iCM );
            }

        }/* TEST_CASE */

    }/* namespace fem */
}/* namespace moris */
