
#include "catch.hpp"
#include "fn_equal_to.hpp"
#include "cl_FEM_IWG_Factory.hpp"              //FEM/INT/src
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

        TEST_CASE( "Set", "[moris],[fem],[FEMSet]" )
        {
            //create a set
            Set tSet;

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
            tSet.mIWGs = tIWGs;

            std::cout<<"Test create_unique_dof_type_lists"<<std::endl;
            //------------------------------------------------------------------------------

                // call create_uniaue_dof_type_lists
                tSet.create_unique_dof_type_lists();

                // check the size of the list
                CHECK( equal_to( tSet.mEqnObjDofTypeList.size(), 6 ) );

                // check the content of the list
                CHECK( equal_to( static_cast< uint >( tSet.mEqnObjDofTypeList( 0 ) ), 3 ) ); //TEMP
                CHECK( equal_to( static_cast< uint >( tSet.mEqnObjDofTypeList( 1 ) ), 11 ) ); //VX
                CHECK( equal_to( static_cast< uint >( tSet.mEqnObjDofTypeList( 2 ) ), 8 ) ); //NLSX
                CHECK( equal_to( static_cast< uint >( tSet.mEqnObjDofTypeList( 3 ) ), 9 ) ); //NLSY
                CHECK( equal_to( static_cast< uint >( tSet.mEqnObjDofTypeList( 4 ) ), 10 ) ); //NLSZ
                CHECK( equal_to( static_cast< uint >( tSet.mEqnObjDofTypeList( 5 ) ), 6 ) ); //LS1


            std::cout<<"Test create_dof_type_lists"<<std::endl;
            //------------------------------------------------------------------------------

                // call create_uniaue_dof_type_lists
                tSet.create_dof_type_lists();

                // check the size of mInterpDofTypeList
                CHECK( equal_to( tSet.mInterpDofTypeList.size(), 4 ) );

                // check the content of mInterpDofTypeList
                CHECK( equal_to( static_cast< uint >( tSet.mInterpDofTypeList( 0 )( 0 ) ), 3 ) );
                CHECK( equal_to( static_cast< uint >( tSet.mInterpDofTypeList( 1 )( 0 ) ), 11 ) );
                CHECK( equal_to( static_cast< uint >( tSet.mInterpDofTypeList( 2 )( 0 ) ), 8 ) );
                CHECK( equal_to( static_cast< uint >( tSet.mInterpDofTypeList( 2 )( 1 ) ), 9 ) );
                CHECK( equal_to( static_cast< uint >( tSet.mInterpDofTypeList( 2 )( 2 ) ), 10 ) );
                CHECK( equal_to( static_cast< uint >( tSet.mInterpDofTypeList( 3 )( 0 ) ), 6 ) );

                // check mNumOfInterp
                CHECK( equal_to( tSet.mNumOfInterp, 4 ) );

                // check mInterpDofTypeMap size
                CHECK( equal_to( tSet.mInterpDofTypeMap.n_cols(), 1 ) );
                CHECK( equal_to( tSet.mInterpDofTypeMap.n_rows(), 12 ) );

                // check the content of mInterpDofTypeMap
                CHECK( equal_to( static_cast< uint >( tSet.mInterpDofTypeMap(  0, 0 ) ), -1 ) );
                CHECK( equal_to( static_cast< uint >( tSet.mInterpDofTypeMap(  1, 0 ) ), -1 ) );
                CHECK( equal_to( static_cast< uint >( tSet.mInterpDofTypeMap(  2, 0 ) ), -1 ) );
                CHECK( equal_to( static_cast< uint >( tSet.mInterpDofTypeMap(  3, 0 ) ), 0 ) );
                CHECK( equal_to( static_cast< uint >( tSet.mInterpDofTypeMap(  4, 0 ) ), -1 ) );
                CHECK( equal_to( static_cast< uint >( tSet.mInterpDofTypeMap(  5, 0 ) ), -1 ) );
                CHECK( equal_to( static_cast< uint >( tSet.mInterpDofTypeMap(  6, 0 ) ), 3 ) );
                CHECK( equal_to( static_cast< uint >( tSet.mInterpDofTypeMap(  7, 0 ) ), -1 ) );
                CHECK( equal_to( static_cast< uint >( tSet.mInterpDofTypeMap(  8, 0 ) ), 2 ) );
                CHECK( equal_to( static_cast< uint >( tSet.mInterpDofTypeMap(  9, 0 ) ), -1 ) );
                CHECK( equal_to( static_cast< uint >( tSet.mInterpDofTypeMap( 10, 0 ) ), -1 ) );
                CHECK( equal_to( static_cast< uint >( tSet.mInterpDofTypeMap( 11, 0 ) ), 1 ) );


            std::cout<<"Test create_dof_assembly_map"<<std::endl;
            //------------------------------------------------------------------------------

                // create a cell of field interpolator pointers---------------------------------
                moris::Cell< Field_Interpolator* > tFieldInterpolators( tSet.mNumOfInterp, nullptr );

                // set the number of coefficients for each field interpolator
                tFieldInterpolators( 0 ) = new Field_Interpolator( 1 );
                tFieldInterpolators( 1 ) = new Field_Interpolator( 1 );
                tFieldInterpolators( 2 ) = new Field_Interpolator( 3 );
                tFieldInterpolators( 3 ) = new Field_Interpolator( 1 );

                // pass in the cell of FI pointers to the element block
                tSet.mMasterFieldInterpolators = tFieldInterpolators;

                // create dof assembly map
                tSet.create_dof_assembly_map();

                // check mInterpDofAssemblyMap size
                CHECK( equal_to( tSet.mInterpDofAssemblyMap.n_cols(), 2 ) );
                CHECK( equal_to( tSet.mInterpDofAssemblyMap.n_rows(), 4 ) );

                // check mInterpDofAssemblyMap content
                CHECK( equal_to( static_cast< uint >( tSet.mInterpDofAssemblyMap( 0, 0 ) ), 0 ) );
                CHECK( equal_to( static_cast< uint >( tSet.mInterpDofAssemblyMap( 0, 1 ) ), 0 ) );
                CHECK( equal_to( static_cast< uint >( tSet.mInterpDofAssemblyMap( 1, 0 ) ), 1 ) );
                CHECK( equal_to( static_cast< uint >( tSet.mInterpDofAssemblyMap( 1, 1 ) ), 1 ) );
                CHECK( equal_to( static_cast< uint >( tSet.mInterpDofAssemblyMap( 2, 0 ) ), 2 ) );
                CHECK( equal_to( static_cast< uint >( tSet.mInterpDofAssemblyMap( 2, 1 ) ), 4 ) );
                CHECK( equal_to( static_cast< uint >( tSet.mInterpDofAssemblyMap( 3, 0 ) ), 5 ) );
                CHECK( equal_to( static_cast< uint >( tSet.mInterpDofAssemblyMap( 3, 1 ) ), 5 ) );

            std::cout<<"Test create_unique_property_type_list"<<std::endl;
            //------------------------------------------------------------------------------

                // create unique list of property type
                tSet.create_unique_property_type_list();

                // check mPropertyTypeList size
                CHECK( equal_to( tSet.mPropertyTypeList.size(), 2 ) );

                // check mInterpDofAssemblyMap content
                CHECK( equal_to( static_cast< uint >( tSet.mPropertyTypeList( 0 ) ), 3 ) );
                CHECK( equal_to( static_cast< uint >( tSet.mPropertyTypeList( 1 ) ), 1 ) );

            std::cout<<"Test create_IWG_set_info"<<std::endl;
            //------------------------------------------------------------------------------

                // create IWG set info
                tSet.create_IWG_set_info();

                // check number of IWGs
                CHECK( equal_to( tSet.mNumOfIWG, 4 ) );

                // check number of active dof per IWG
                CHECK( equal_to( tSet.mIWGNumActiveDof( 0 ), 1 ) );
                CHECK( equal_to( tSet.mIWGNumActiveDof( 1 ), 1 ) );
                CHECK( equal_to( tSet.mIWGNumActiveDof( 2 ), 1 ) );
                CHECK( equal_to( tSet.mIWGNumActiveDof( 3 ), 2 ) );

                // check number of FI per IWG
                CHECK( equal_to( tSet.mIWGMasterFieldInterpolators( 0 ).size(), 1 ) );
                CHECK( equal_to( tSet.mIWGMasterFieldInterpolators( 1 ).size(), 1 ) );
                CHECK( equal_to( tSet.mIWGMasterFieldInterpolators( 2 ).size(), 1 ) );
                CHECK( equal_to( tSet.mIWGMasterFieldInterpolators( 3 ).size(), 2 ) );

                // check dofAssembly matrix per IWG
                CHECK( equal_to( tSet.mIWGDofAssemblyMap( 0 )( 0, 0 ), 0 ) );
                CHECK( equal_to( tSet.mIWGDofAssemblyMap( 0 )( 0, 1 ), 0 ) );
                CHECK( equal_to( tSet.mIWGDofAssemblyMap( 0 )( 0, 2 ), 0 ) );
                CHECK( equal_to( tSet.mIWGDofAssemblyMap( 0 )( 0, 3 ), 0 ) );

                CHECK( equal_to( tSet.mIWGDofAssemblyMap( 1 )( 0, 0 ), 0 ) );
                CHECK( equal_to( tSet.mIWGDofAssemblyMap( 1 )( 0, 1 ), 0 ) );
                CHECK( equal_to( tSet.mIWGDofAssemblyMap( 1 )( 0, 2 ), 0 ) );
                CHECK( equal_to( tSet.mIWGDofAssemblyMap( 1 )( 0, 3 ), 0 ) );

                CHECK( equal_to( tSet.mIWGDofAssemblyMap( 2 )( 0, 0 ), 1 ) );
                CHECK( equal_to( tSet.mIWGDofAssemblyMap( 2 )( 0, 1 ), 1 ) );
                CHECK( equal_to( tSet.mIWGDofAssemblyMap( 2 )( 0, 2 ), 1 ) );
                CHECK( equal_to( tSet.mIWGDofAssemblyMap( 2 )( 0, 3 ), 1 ) );

                CHECK( equal_to( tSet.mIWGDofAssemblyMap( 3 )( 0, 0 ), 2 ) );
                CHECK( equal_to( tSet.mIWGDofAssemblyMap( 3 )( 0, 1 ), 4 ) );
                CHECK( equal_to( tSet.mIWGDofAssemblyMap( 3 )( 0, 2 ), 2 ) );
                CHECK( equal_to( tSet.mIWGDofAssemblyMap( 3 )( 0, 3 ), 4 ) );
                CHECK( equal_to( tSet.mIWGDofAssemblyMap( 3 )( 1, 0 ), 2 ) );
                CHECK( equal_to( tSet.mIWGDofAssemblyMap( 3 )( 1, 1 ), 4 ) );
                CHECK( equal_to( tSet.mIWGDofAssemblyMap( 3 )( 1, 2 ), 5 ) );
                CHECK( equal_to( tSet.mIWGDofAssemblyMap( 3 )( 1, 3 ), 5 ) );

            // clean up
            //------------------------------------------------------------------------------
            for ( uint i = 0; i < tIWGs.size(); i++ )
            {
                delete tIWGs( i );
            }

        }/* TEST_CASE */
    }/* namespace fem */
}/* namespace moris */
