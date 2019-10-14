
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
                                               moris::Cell< moris::fem::Field_Interpolator* > & aFieldInterpolators,
                                               moris::fem::Geometry_Interpolator              * aGeometryInterpolator )
{
    return aParameters( 0 ) + aParameters( 1 ) * aFieldInterpolators( 0 )->val();
}

moris::Matrix< moris::DDRMat > tDerFunctionCM( moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                                               moris::Cell< moris::fem::Field_Interpolator* > & aFieldInterpolators,
                                               moris::fem::Geometry_Interpolator              * aGeometryInterpolator )
{
    return aParameters( 1 ) * aFieldInterpolators( 0 )->N();
}

namespace moris
{
    namespace fem
    {
        TEST_CASE( "Constitutive_Model", "[moris],[fem],[CM]" )
        {
            // real for check
            real tEpsilon = 1E-6;

            // create a CM factory
            CM_Factory tCMFactory;

            // create a constitutive model
            Constitutive_Model* tCM = tCMFactory.create_CM( Constitutive_Type::DIFF_LIN_ISO );

            // set space dim
            tCM->set_space_dim( 2 );

            // set dof types
            tCM->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );

            // set property type
            tCM->set_property_type_list( { fem::Property_Type::CONDUCTIVITY } );

            //create a quad4 space element
            Matrix< DDRMat > tXHat( 4, 2 );
//            tXHat( 0, 0 ) = 0.0; tXHat( 0, 1 ) = 0.0;
//            tXHat( 1, 0 ) = 3.0; tXHat( 1, 1 ) = 1.25;
//            tXHat( 2, 0 ) = 4.5; tXHat( 2, 1 ) = 4.0;
//            tXHat( 3, 0 ) = 1.0; tXHat( 3, 1 ) = 3.25;
            tXHat( 0, 0 ) = 0.0; tXHat( 0, 1 ) = 0.0;
            tXHat( 1, 0 ) = 1.0; tXHat( 1, 1 ) = 0.0;
            tXHat( 2, 0 ) = 1.0; tXHat( 2, 1 ) = 1.0;
            tXHat( 3, 0 ) = 0.0; tXHat( 3, 1 ) = 1.0;

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
            Geometry_Interpolator* tGI = new Geometry_Interpolator( tGeomInterpRule );

            //set the coefficients xHat, tHat
            tGI->set_coeff( tXHat, tTHat );
            tGI->set_space_time({{0.0}, {0.0}, {-1.0}});

            // create a property object
            Cell< fem::Property* > tProps( 1, nullptr );
            tProps( 0 ) = new Property( fem::Property_Type::CONDUCTIVITY ,
                                        {{ MSI::Dof_Type::TEMP }},
                                        {{{ 1.0}}, {{1.0 }}},
                                        tValFunctionCM,
                                        { tDerFunctionCM },
                                        tGI );

            // create an interpolation rule
            Interpolation_Rule tIPRule ( mtk::Geometry_Type::QUAD,
                                         Interpolation_Type::LAGRANGE,
                                         mtk::Interpolation_Order::LINEAR,
                                         Interpolation_Type::LAGRANGE,
                                         mtk::Interpolation_Order::LINEAR );

            // create a TEMP field interpolator
            Cell< Field_Interpolator* > tFIs( 1, nullptr );
            tFIs( 0 ) = new Field_Interpolator ( 1, tIPRule, tGI, { MSI::Dof_Type::TEMP } );

            // set coefficients for field interpolators
            Matrix< DDRMat > tUHat0 = {{1.0},{1.0},{2.0},{2.0},{1.0},{1.0},{2.0},{2.0}};
            tFIs( 0 )->set_coeff( tUHat0 );
            tFIs( 0 )->set_space_time({{0.0}, {0.0}, {-1.0}});

            // set FI for properties
            tProps( 0 )->set_field_interpolators( tFIs );

            // set properties
            tCM->set_properties( tProps );

            // set field interpolators
            tCM->set_field_interpolators( tFIs );

            // checks
            //------------------------------------------------------------------------------

            // check global dof list size
            CHECK( equal_to( tCM->get_global_dof_type_list().size(), 1 ));

            // check global dof list content
            CHECK( equal_to( static_cast< uint >( tCM->get_global_dof_type_list()( 0 )( 0 ) ), 3 ) ); //TEMP

            // evaluate the constitutive model stress
            Matrix< DDRMat > tFlux;
            tCM->eval_flux( tFlux );
            //print( tStress, "tStress");

            // evaluate the constitutive model stress derivative
            Matrix< DDRMat > tdFluxdDOF;
            tCM->eval_dFluxdDOF( { MSI::Dof_Type::TEMP }, tdFluxdDOF );

            // evaluate the constitutive model stress derivative by FD
            Matrix< DDRMat > tdFluxdDOF_FD;
            tCM->eval_dFluxdDOF_FD( { MSI::Dof_Type::TEMP }, tdFluxdDOF_FD, 1E-6 );

            //check stress derivative
            bool tCheckdStress = true;
            for ( uint iStress = 0; iStress < tdFluxdDOF.n_rows(); iStress++ )
            {
                for( uint jStress = 0; jStress < tdFluxdDOF.n_cols(); jStress++ )
                {
                    tCheckdStress = tCheckdStress && ( tdFluxdDOF( iStress, jStress ) - tdFluxdDOF_FD( iStress, jStress ) < tEpsilon );
                }
            }
            REQUIRE( tCheckdStress );

            // evaluate the constitutive model strain
            Matrix< DDRMat > tStrain;
            tCM->eval_strain( tStrain );
            //print( tStrain, "tStrain");

            // evaluate the constitutive model strain derivative
            Matrix< DDRMat > tdStraindDOF;
            tCM->eval_dStraindDOF( { MSI::Dof_Type::TEMP }, tdStraindDOF );
            //print( tdStraindDOF, "tdStraindDOF" );

            // evaluate the constitutive model strain derivative by FD
            Matrix< DDRMat > tdStraindDOF_FD;
            tCM->eval_dStraindDOF_FD( { MSI::Dof_Type::TEMP }, tdStraindDOF_FD, 1E-6 );

            //check strain derivative
            bool tCheckdStrain = true;
            for ( uint iStress = 0; iStress < tdStraindDOF.n_rows(); iStress++ )
            {
                for( uint jStress = 0; jStress < tdStraindDOF.n_cols(); jStress++ )
                {
                    tCheckdStrain = tCheckdStrain && ( tdStraindDOF( iStress, jStress ) - tdStraindDOF_FD( iStress, jStress ) < tEpsilon );
                }
            }
            REQUIRE( tCheckdStrain );

            // evaluate the constitutive model constitutive matrix
            Matrix< DDRMat > tConst;
            tCM->eval_const( tConst );
            print( tConst, "tConst");

            // evaluate the constitutive model constitutive matrix derivative
            Matrix< DDRMat > tdConstdDOF;
            tCM->eval_dConstdDOF( { MSI::Dof_Type::TEMP }, tdConstdDOF );
            print( tdConstdDOF, "tdConstdDOF" );

            // clean up
            //------------------------------------------------------------------------------

            // delete the property pointers
            for( Property* tProp : tProps )
            {
                delete tProp;
            }
            tProps.clear();

            // delete the field interpolator pointers
            for( Field_Interpolator* tFI : tFIs )
            {
                delete tFI;
            }
            tFIs.clear();

        }/* TEST_CASE */


    }/* namespace fem */
}/* namespace moris */
