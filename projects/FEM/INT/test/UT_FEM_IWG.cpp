
#include "catch.hpp"
#include "fn_equal_to.hpp"
#include "cl_FEM_Property.hpp"              //FEM/INT/src
#include "cl_FEM_Geometry_Interpolator.hpp"         //FEM/INT/src
#include "cl_FEM_Field_Interpolator.hpp"         //FEM/INT/src

#define protected public
#define private   public
#include "cl_FEM_IWG.hpp"         //FEM/INT/src
#undef protected
#undef private

#include "cl_FEM_IWG_Factory.hpp"         //FEM/INT/src

namespace moris
{
    namespace fem
    {

        Matrix< DDRMat > tValFunction( moris::Cell< Matrix< DDRMat > >    & aCoeff,
                                       moris::Cell< Field_Interpolator* > & aFieldInterpolator )
        {
            Matrix< DDRMat > tPropertyVal( 1, 1, 1.0);
            return tPropertyVal;
        }

        Matrix< DDRMat > tDerFunction( moris::Cell< Matrix< DDRMat > >    & aCoeff,
                                       moris::Cell< Field_Interpolator* > & aFieldInterpolator )
        {
            Matrix< DDRMat > tPropertyDer( 1, 1, 2.0);
            return tPropertyDer;
        }

        Matrix< DDRMat > tValFunction2( moris::Cell< Matrix< DDRMat > >    & aCoeff,
                                        moris::Cell< Field_Interpolator* > & aFieldInterpolator )
        {
            Matrix< DDRMat > tPropertyVal( 1, 1, 0.0);
            tPropertyVal = aCoeff( 0 ) + aCoeff( 1 ) * aFieldInterpolator( 0 )->val() + aCoeff( 2 ) * aFieldInterpolator( 1 )->val();
            return tPropertyVal;
        }

        Matrix< DDRMat > tDerFunction2( moris::Cell< Matrix< DDRMat > >    & aCoeff,
                                        moris::Cell< Field_Interpolator* > & aFieldInterpolator )
        {
            Matrix< DDRMat > tPropertyDer;
            tPropertyDer = aCoeff( 1 ) * aFieldInterpolator( 0 )->N();
            return tPropertyDer;
        }

        Matrix< DDRMat > tDerFunction3( moris::Cell< Matrix< DDRMat > >    & aCoeff,
                                        moris::Cell< Field_Interpolator* > & aFieldInterpolator )
        {
            Matrix< DDRMat > tPropertyDer;
            tPropertyDer = aCoeff( 2 ) * aFieldInterpolator( 1 )->N();
            return tPropertyDer;
        }

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

            // create a list of dof dependencies for the property
            Cell< Cell< MSI::Dof_Type > >  tActiveDofTypes = {{ MSI::Dof_Type::TEMP },
                                                              { MSI::Dof_Type::VX }};

            // create the function pointers for the value
            fem::PropertyFunc tValFunction0 = tValFunction;
            Cell< fem::PropertyFunc > tDerFunctions( 0 );

            // create a property object
            Cell< fem::Property* > tProperties( 1, nullptr );
            tProperties( 0 ) = new Property( fem::Property_Type::TEMP_NEUMANN,
                                             tActiveDofTypes,
                                             tValFunction0, tDerFunctions );

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

            // create an interpolation rule
            Interpolation_Rule tInterpolationRule ( mtk::Geometry_Type::QUAD,
                                                    Interpolation_Type::LAGRANGE,
                                                    mtk::Interpolation_Order::LINEAR,
                                                    Interpolation_Type::LAGRANGE,
                                                    mtk::Interpolation_Order::LINEAR );

            // create a TEMP and a VX field interpolator
            uint tNumberOfFields = 1;
            Cell< Field_Interpolator* > tFieldInterpolators( 2, nullptr );
            tFieldInterpolators( 0 ) = new Field_Interpolator ( tNumberOfFields,
                                                               tInterpolationRule,
                                                               tGeomInterpolator,
                                                               { MSI::Dof_Type::TEMP } );
            tFieldInterpolators( 1 ) = new Field_Interpolator ( tNumberOfFields,
                                                                tInterpolationRule,
                                                                tGeomInterpolator,
                                                                { MSI::Dof_Type::VX } );

            // set coefficients for field interpolators
            Matrix< DDRMat > tUHat0( 8, 1, 2.0 );
            tFieldInterpolators( 0 )->set_coeff( tUHat0 );
            Matrix< DDRMat > tUHat1( 8, 1, 3.0 );
            tFieldInterpolators( 1 )->set_coeff( tUHat1 );

            // a factory to create the IWGs
            fem::IWG_Factory tIWGFactory;

            // create an IWG with the factory for the ith IWG type
            fem::IWG* tIWG = tIWGFactory.create_IWGs( fem::IWG_Type::SPATIALDIFF_NEUMANN );

            // set IWG properties
            tIWG->set_properties( tProperties );

            // set IWG field interpolators
            tIWG->set_field_interpolators( tFieldInterpolators );

            //check global dof list size
            CHECK( equal_to( tIWG->mMasterGlobalDofTypes.size(), 2 ));

            //check global dof list content
            CHECK( equal_to( static_cast< uint >( tIWG->mMasterGlobalDofTypes( 0 )( 0 ) ), 3 ) );
            CHECK( equal_to( static_cast< uint >( tIWG->mMasterGlobalDofTypes( 1 )( 0 ) ), 11 ) );

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

        }/* TEST_CASE */


    }/* namespace fem */
}/* namespace moris */
