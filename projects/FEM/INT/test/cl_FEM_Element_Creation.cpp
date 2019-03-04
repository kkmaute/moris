
#include "catch.hpp"

#include "cl_MTK_Vertex.hpp"    //MTK
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Enums.hpp"

#include "cl_Matrix.hpp"        //LINALG
#include "linalg_typedefs.hpp"

#include "cl_FEM_NodeProxy.hpp"                //FEM/INT/src
#include "cl_FEM_ElementProxy.hpp"             //FEM/INT/src
#include "cl_FEM_Node_Base.hpp"                //FEM/INT/src

#include "cl_FEM_Element_Factory.hpp"          //FEM/INT/src
#include "cl_FEM_IWG_Factory.hpp"              //FEM/INT/src

namespace moris
{
    namespace fem
    {
        TEST_CASE( "Element_Creation", "[moris],[fem],[ElementCreate]" )
        {
            //------------------------------------------------------------------------------
            // create nodes
            mtk::Vertex* tVertex1_1 = new NodeProxy( 0.0, -1.0, 0 );
            mtk::Vertex* tVertex1_2 = new NodeProxy( 1.0, -1.0, 1 );
            mtk::Vertex* tVertex1_3 = new NodeProxy( 1.0,  0.0, 2 );
            mtk::Vertex* tVertex1_4 = new NodeProxy( 0.0,  0.0, 3 );

            mtk::Vertex* tVertex2_3 = new NodeProxy( 1.0,  1.0, 4 );
            mtk::Vertex* tVertex2_4 = new NodeProxy( 0.0,  1.0, 5 );

            mtk::Vertex* tVertex3_3 = new NodeProxy( 1.0,  2.0, 6 );
            mtk::Vertex* tVertex3_4 = new NodeProxy( 0.0,  2.0, 7 );

            mtk::Vertex* tVertex4_3 = new NodeProxy( 1.0,  3.0, 8 );
            mtk::Vertex* tVertex4_4 = new NodeProxy( 0.0,  3.0, 9 );

            moris::Cell< mtk::Vertex* > allNodes( 10 );
            allNodes( 0 ) = tVertex1_1;
            allNodes( 1 ) = tVertex1_2;
            allNodes( 2 ) = tVertex1_3;
            allNodes( 3 ) = tVertex1_4;
            allNodes( 4 ) = tVertex2_3;
            allNodes( 5 ) = tVertex2_4;
            allNodes( 6 ) = tVertex3_3;
            allNodes( 7 ) = tVertex3_4;
            allNodes( 8 ) = tVertex4_3;
            allNodes( 9 ) = tVertex4_4;

            //------------------------------------------------------------------------------
            // create cells of nodes for element

            // cell containing nodes for element 1
            moris::Cell< mtk::Vertex* > Name1(4);
            Name1(0) = tVertex1_1;
            Name1(1) = tVertex1_2;
            Name1(2) = tVertex1_3;
            Name1(3) = tVertex1_4;

            // cell containing nodes for element 2
            moris::Cell< mtk::Vertex* > Name2(4);
            Name2(0) = tVertex1_4;
            Name2(1) = tVertex1_3;
            Name2(2) = tVertex2_3;
            Name2(3) = tVertex2_4;

            // cell containing nodes for element 1
            moris::Cell< mtk::Vertex* > Name3(4);
            Name3(0) = tVertex2_4;
            Name3(1) = tVertex2_3;
            Name3(2) = tVertex3_3;
            Name3(3) = tVertex3_4;

            // cell containing nodes for element 1
            moris::Cell< mtk::Vertex* > Name4(4);
            Name4(0) = tVertex3_4;
            Name4(1) = tVertex3_4;
            Name4(2) = tVertex4_3;
            Name4(3) = tVertex4_4;

            //------------------------------------------------------------------------------
            // create element as mtk::Cell
            mtk::Cell* tElement1 = new ElementProxy( Name1, mtk::Geometry_Type::QUAD );
            mtk::Cell* tElement2 = new ElementProxy( Name2, mtk::Geometry_Type::QUAD );
            mtk::Cell* tElement3 = new ElementProxy( Name3, mtk::Geometry_Type::QUAD );
            mtk::Cell* tElement4 = new ElementProxy( Name4, mtk::Geometry_Type::QUAD );

            //------------------------------------------------------------------------------
            // create a mesh
            moris::Cell< mtk::Cell* > Elems(4);
            Elems(0) = tElement1;
            Elems(1) = tElement2;
            Elems(2) = tElement3;
            Elems(3) = tElement4;

            // 1) Create the fem nodes -----------------------------------------------------
            //------------------------------------------------------------------------------
            // number of mesh nodes
            uint tNumOfNodes = allNodes.size();

            //create a celle of fem nodes
            moris::Cell< Node_Base* > tNodes( tNumOfNodes, nullptr );

            // loop obver the mesh nodes
            for( uint k = 0; k < tNumOfNodes; k++ )
            {
                // create a fem node for each mesh node
                tNodes( k ) = new fem::Node( allNodes( k ) );
            }

            // 2) create the IWGs for the element 1------------------------------------------
            //-------------------------------------------------------------------------------
            // input a cell of IWG types to be created
            Cell< IWG_Type > tIWGTypeList = { IWG_Type::HELMHOLTZ, IWG_Type::HJ };

            // number of IWGs to be created
            uint tNumOfIWGs = tIWGTypeList.size();

            // a factory to create the IWGs
            IWG_Factory tIWGFactory;

            // create a cell of IWGs for the problem considered
            Cell< IWG* > tIWGs( tNumOfIWGs , nullptr );

            // loop over the IWG types
            for( uint i = 0; i < tNumOfIWGs; i++)
            {
                // create an IWG with the factory for the ith IWG type
                tIWGs( i ) = tIWGFactory.create_IWGs( tIWGTypeList( i ) );
            }

//            std::cout<<static_cast< uint >(tIWGs2( 0 )->get_residual_dof_type())<<std::endl;
//            std::cout<<static_cast< uint >(tIWGs2( 1 )->get_residual_dof_type())<<std::endl;
//            std::cout<<"---------"<<std::endl;
//
//            Cell< MSI::Dof_Type > tIWGActivedofTypes02 = tIWGs2( 0 )->get_active_dof_types();
//            std::cout<<static_cast< uint >(tIWGActivedofTypes02(0))<<std::endl;
//            std::cout<<static_cast< uint >(tIWGActivedofTypes02(1))<<std::endl;
//            Cell< MSI::Dof_Type > tIWGActivedofTypes12 = tIWGs2( 1 )->get_active_dof_types();
//            std::cout<<static_cast< uint >(tIWGActivedofTypes12(0))<<std::endl;
//            std::cout<<"---------"<<std::endl;

            // 3) create elements------------------------------------------------------------
            //-------------------------------------------------------------------------------
            // a factory to create the elements
            Element_Factory tElementFactory;

            // get the number of mesh elements
            uint tNumOfElements = Elems.size();

            // create a cell of fem elements
            Cell< MSI::Equation_Object * > tListOfElements( tNumOfElements, nullptr );

            // loop over the mesh elements
            for ( uint i = 0; i < tNumOfElements; i++ )
            {
                // create a fem element for the ith mesh element
                tListOfElements( i ) = tElementFactory.create_element( Element_Type::UNDEFINED,
                                                                       Elems( i ),
                                                                       tIWGs,
                                                                       tNodes );

                // evaluate the residual and jacobian for the ith fem element
                tListOfElements( i )->compute_residual();
                tListOfElements( i )->compute_jacobian();
                tListOfElements( i )->compute_jacobian_and_residual();
            }

            //-------------------------------------------------------------------------------

//            // create cell of Dof_Type
//            Cell< Cell < MSI::Dof_Type > > tCellDofType( 2 );
//            tCellDofType( 0 ) = {MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ};
//            tCellDofType( 1 ) = {MSI::Dof_Type::TEMP};
//
//            //std::cout<<static_cast< int >( tCellDofType( 0 )( 0 ) )<<std::endl;
//            //std::cout<<static_cast< int >( tCellDofType( 0 )( 1 ) )<<std::endl;
//            //std::cout<<static_cast< int >( tCellDofType( 0 )( 2 ) )<<std::endl;
//            //std::cout<<static_cast< int >( tCellDofType( 1 )( 0 ) )<<std::endl;
//
//            Cell< MSI::Dof_Type > tCellDofType0 = tCellDofType( 0 );
//            //std::cout<<tCellDofType0.size()<<std::endl;
//
//            Cell< MSI::Dof_Type > tCellDofType1 = tCellDofType( 1 );
//            //std::cout<<tCellDofType1.size()<<std::endl;

            //clean up ----------------------------------------------------------------------
            //-------------------------------------------------------------------------------
            for( uint i = 0; i < tNumOfNodes; i++ )
            {
                delete allNodes( i );
                delete tNodes( i );
            }

            for( uint i = 0; i < tNumOfIWGs; i++)
            {
                delete tIWGs( i );
            }

            for( uint i = 0; i < tNumOfElements; i++ )
            {
                delete Elems( i );
                delete tListOfElements( i );
            }


        }/* TEST_CASE */
    }/* namespace fem */
}/* namespace moris */
