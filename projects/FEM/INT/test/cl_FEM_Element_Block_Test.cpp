
#include "catch.hpp"
#include "fn_equal_to.hpp"
#include "cl_FEM_IWG_Factory.hpp"              //FEM/INT/src
#include "cl_FEM_Property.hpp" //FEM/INT/src
#include "cl_FEM_Property_Factory.hpp" //FEM/INT/src
#include "linalg_typedefs.hpp"

#define protected public
#define private   public
#include "cl_FEM_Set.hpp"              //FEM/INT/src
#include "cl_FEM_Field_Interpolator.hpp"         //FEM/INT/src
#undef protected
#undef private

namespace moris
{
    namespace fem
    {
        TEST_CASE( "Property", "[moris],[fem],[Property]" )
        {
            // create a factory
            Property_Factory tPropFactory;

            // create a property object
            Property * tProperty = tPropFactory.create_property( fem::Property_Type::TEMP_DIRICHLET );

            // check property type
            fem::Property_Type tPropertyType = tProperty->get_property_type();
            CHECK( equal_to( static_cast< uint >( tPropertyType ), 1 ) );

            // check val_coeff
            Matrix< DDRMat > tCoeff;
            tProperty->val_coeff( tCoeff );
            print( tCoeff, "tCoeff" );

            // create a field interpolator for the property
            uint tNumberOfFields = 1;

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

            Interpolation_Rule tInterpolationRule ( mtk::Geometry_Type::QUAD,
                                                    Interpolation_Type::LAGRANGE,
                                                    mtk::Interpolation_Order::LINEAR,
                                                    Interpolation_Type::LAGRANGE,
                                                    mtk::Interpolation_Order::LINEAR );

            Field_Interpolator* tFI = new Field_Interpolator ( tNumberOfFields,
                                                               tInterpolationRule,
                                                               tGeomInterpolator,
                                                               tProperty,
                                                               tPropertyType );

            // clean up
            delete tProperty;
            delete tGeomInterpolator;
            delete tFI;
        }

        TEST_CASE( "Set", "[moris],[fem],[ElementBlock]" )
        {
            //MSI::Equation_Set* tElementBlock = new Set();
            Set tElementBlock;

            // create a cell of IWG pointers------------------------------------------------
            // list of IWG types
            moris::Cell< fem::IWG_Type > tIWGTypeList= { fem::IWG_Type::SPATIALDIFF_BULK,
                                                         fem::IWG_Type::SPATIALDIFF_DIRICHLET,
                                                         fem::IWG_Type::HELMHOLTZ,
                                                         fem::IWG_Type::LSNORMAL };
            // number of IWGs to be created
            uint tNumOfIWGs = tIWGTypeList.size();

            // a factory to create the IWGs
            fem::IWG_Factory tIWGFactory;

           // create a cell of IWGs for the problem considered
           moris::Cell< fem::IWG* > tIWGs( tNumOfIWGs , nullptr );

           // loop over the IWG types
           for( uint i = 0; i < tNumOfIWGs; i++)
           {
               // create an IWG with the factory for the ith IWG type
               tIWGs( i ) = tIWGFactory.create_IWGs( tIWGTypeList( i ) );
           }

           // pass in the cell of IWG pointers to the element block
            tElementBlock.mIWGs = tIWGs;

            // create a cell of field interpolator pointers---------------------------------
            moris::Cell< Field_Interpolator* > tFieldInterpolators( 3, nullptr );

            // set the number of coefficients for each field interpolator
            tFieldInterpolators( 0 ) = new Field_Interpolator( 1 );
            tFieldInterpolators( 1 ) = new Field_Interpolator( 3 );
            tFieldInterpolators( 2 ) = new Field_Interpolator( 1 );

            //tFieldInterpolators( 0 )->mNFieldCoeff = 4;
            //tFieldInterpolators( 1 )->mNFieldCoeff = 8;
            //tFieldInterpolators( 2 )->mNFieldCoeff = 4;

            // pass in the cell of FI pointers to the element block
            tElementBlock.mMasterFieldInterpolators = tFieldInterpolators;

            SECTION("Test create_unique_dof_type_lists")
            {
                // call create_uniaue_dof_type_lists
                tElementBlock.create_unique_dof_type_lists();

                // check the size of the list
                CHECK( equal_to( tElementBlock.mEqnObjDofTypeList.size(), 5 ) );

                // check the content of the list
                CHECK( equal_to( static_cast< uint >( tElementBlock.mEqnObjDofTypeList( 0 ) ), 3 ) );
                CHECK( equal_to( static_cast< uint >( tElementBlock.mEqnObjDofTypeList( 1 ) ), 11 ) );
                CHECK( equal_to( static_cast< uint >( tElementBlock.mEqnObjDofTypeList( 2 ) ), 8 ) );
                CHECK( equal_to( static_cast< uint >( tElementBlock.mEqnObjDofTypeList( 3 ) ), 9 ) );
                CHECK( equal_to( static_cast< uint >( tElementBlock.mEqnObjDofTypeList( 4 ) ), 10 ) );
            };


            SECTION("Test create_unique_dof_type_lists")
            {
                // call create_uniaue_dof_type_lists
                tElementBlock.create_dof_type_lists();

                // check the size of mInterpDofTypeList
                CHECK( equal_to( tElementBlock.mInterpDofTypeList.size(), 3 ) );

                // check the content of mInterpDofTypeList
                CHECK( equal_to( static_cast< uint >( tElementBlock.mInterpDofTypeList( 0 )( 0 ) ), 3 ) );
                CHECK( equal_to( static_cast< uint >( tElementBlock.mInterpDofTypeList( 1 )( 0 ) ), 8 ) );
                CHECK( equal_to( static_cast< uint >( tElementBlock.mInterpDofTypeList( 1 )( 1 ) ), 9 ) );
                CHECK( equal_to( static_cast< uint >( tElementBlock.mInterpDofTypeList( 1 )( 2 ) ), 10 ) );
                CHECK( equal_to( static_cast< uint >( tElementBlock.mInterpDofTypeList( 2 )( 0 ) ), 11 ) );

                // check mNumOfInterp
                CHECK( equal_to( tElementBlock.mNumOfInterp, 3 ) );

                // check mInterpDofTypeMap size
                CHECK( equal_to( tElementBlock.mInterpDofTypeMap.n_cols(), 1 ) );
                CHECK( equal_to( tElementBlock.mInterpDofTypeMap.n_rows(), 12 ) );

                // check the content of mInterpDofTypeMap
                CHECK( equal_to( static_cast< uint >( tElementBlock.mInterpDofTypeMap(  0, 0 ) ), -1 ) );
                CHECK( equal_to( static_cast< uint >( tElementBlock.mInterpDofTypeMap(  1, 0 ) ), -1 ) );
                CHECK( equal_to( static_cast< uint >( tElementBlock.mInterpDofTypeMap(  2, 0 ) ), -1 ) );
                CHECK( equal_to( static_cast< uint >( tElementBlock.mInterpDofTypeMap(  3, 0 ) ), 0 ) );
                CHECK( equal_to( static_cast< uint >( tElementBlock.mInterpDofTypeMap(  4, 0 ) ), -1 ) );
                CHECK( equal_to( static_cast< uint >( tElementBlock.mInterpDofTypeMap(  5, 0 ) ), -1 ) );
                CHECK( equal_to( static_cast< uint >( tElementBlock.mInterpDofTypeMap(  6, 0 ) ), -1 ) );
                CHECK( equal_to( static_cast< uint >( tElementBlock.mInterpDofTypeMap(  7, 0 ) ), -1 ) );
                CHECK( equal_to( static_cast< uint >( tElementBlock.mInterpDofTypeMap(  8, 0 ) ), 1 ) );
                CHECK( equal_to( static_cast< uint >( tElementBlock.mInterpDofTypeMap(  9, 0 ) ), -1 ) );
                CHECK( equal_to( static_cast< uint >( tElementBlock.mInterpDofTypeMap( 10, 0 ) ), -1 ) );
                CHECK( equal_to( static_cast< uint >( tElementBlock.mInterpDofTypeMap( 11, 0 ) ), 2 ) );
            };

            SECTION("Test create_dof_assembly_map")
            {
//                // set the number of field interpolators
//                tElementBlock.mNumOfInterp = tElementBlock.mFieldInterpolators.size();
//
//                // check mInterpDofAssemblyMap size
//                CHECK( equal_to( tElementBlock.mInterpDofAssemblyMap.n_cols(), 2 ) );
//                CHECK( equal_to( tElementBlock.mInterpDofAssemblyMap.n_rows(), 3 ) );
//
//                // check mInterpDofAssemblyMap content
//                CHECK( equal_to( static_cast< uint >( tElementBlock.mInterpDofAssemblyMap( 0, 0 ) ), 0 ) );
//                CHECK( equal_to( static_cast< uint >( tElementBlock.mInterpDofAssemblyMap( 0, 1 ) ), 3 ) );
//                CHECK( equal_to( static_cast< uint >( tElementBlock.mInterpDofAssemblyMap( 1, 0 ) ), 4 ) );
//                CHECK( equal_to( static_cast< uint >( tElementBlock.mInterpDofAssemblyMap( 1, 1 ) ), 11 ) );
//                CHECK( equal_to( static_cast< uint >( tElementBlock.mInterpDofAssemblyMap( 2, 0 ) ), 12 ) );
//                CHECK( equal_to( static_cast< uint >( tElementBlock.mInterpDofAssemblyMap( 2, 1 ) ), 15 ) );
            };

            // clean up
            for ( uint i = 0; i < tIWGs.size(); i++ )
            {
                delete tIWGs( i );
            }
//            for ( uint i = 0; i < tFieldInterpolators.size(); i++ )
//            {
//                delete tFieldInterpolators( i );
//            }

        }/* TEST_CASE */
    }/* namespace fem */
}/* namespace moris */
