
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
#include "cl_FEM_Element.hpp"                  //FEM/INT/src

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

            //------------------------------------------------------------------------------

            // 1) Create the fem nodes for element 1----------------------------------------
            uint tNumberOfNodes = 10;
            moris::Cell< Node_Base* > tNodes( tNumberOfNodes, nullptr );
            for( uint k = 0; k < tNumberOfNodes; k++ )
            {
                tNodes( k ) = new fem::Node( allNodes( k ) );
            }

            // 2) create the IWGs for the element 1------------------------------------------
            // create IWGs with a factory
            IWG_Factory tIWGFactory;
            Cell< IWG* > tIWGs = tIWGFactory.create_IWGs( Element_Type::HJ );

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
            // create fem elements with a factory
            Element_Factory tElementFactory;

            // get the number of mesh elements
            uint tNumberOfElements = Elems.size();

            // loop over the mesh elements to create the fem elements
            moris::Cell< MSI::Equation_Object * > tListOfElements( tNumberOfElements, nullptr );
            for ( uint i = 0; i < tNumberOfElements; i++ )
            {
                tListOfElements( i ) = tElementFactory.create_element( Element_Type::UNDEFINED,
                                                                       Elems( i ),
                                                                       tIWGs,
                                                                       tNodes );
            }

            // evaluate the residual and jacobian
            tListOfElements( 0 )->compute_residual();
            tListOfElements( 0 )->compute_jacobian();

            tListOfElements( 1 )->compute_residual();
            tListOfElements( 1 )->compute_jacobian();

            tListOfElements( 2 )->compute_residual();
            tListOfElements( 2 )->compute_jacobian();

            tListOfElements( 3 )->compute_residual();
            tListOfElements( 3 )->compute_jacobian();

        }/* TEST_CASE */
    }/* namespace fem */
}/* namespace moris */
