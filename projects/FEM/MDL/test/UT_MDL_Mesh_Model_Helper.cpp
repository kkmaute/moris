
#include "catch.hpp"

#include "cl_MTK_Vertex.hpp"    //MTK
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh.hpp"

#include "cl_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Tools.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Scalar_Field_Info.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Mesh_Data_STK.hpp"
#include "cl_MTK_Mesh_Core_STK.hpp"
#include "cl_MTK_Interpolation_Mesh_STK.hpp"
#include "cl_MTK_Integration_Mesh_STK.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"

#include "cl_Matrix.hpp"        //LINALG
#include "linalg_typedefs.hpp"
#include "fn_equal_to.hpp" // ALG/src

#include "fn_norm.hpp"

//#include "cl_MDL_Node_Proxy.hpp"

namespace moris
{
    namespace mdl
    {
        TEST_CASE( "Mesh_Model_Helper_Test", "[moris],[mdl],[Mesh_Model_Helper]" )
        {
            uint p_rank = moris::par_rank();
            uint p_size = moris::par_size();

            if( p_size == 1 ) // specify it is a serial test only
            {
                /* create 2D 4-element mesh, add circle function, check for intersection */
//                mtk::Vertex* tVertex1_1 = new Node( 0, 0 );
//                mtk::Vertex* tVertex1_2 = new Node( 1, 1 );
//                mtk::Vertex* tVertex1_3 = new Node( 2, 2 );
//                mtk::Vertex* tVertex1_4 = new Node( 3, 3 );
//
//                mtk::Vertex* tVertex2_3 = new Node(1.0,  1.0);
//                mtk::Vertex* tVertex2_4 = new Node(0.0,  1.0);
//
//                mtk::Vertex* tVertex3_3 = new Node(1.0,  2.0);
//                mtk::Vertex* tVertex3_4 = new Node(0.0,  2.0);
//
//                mtk::Vertex* tVertex4_3 = new Node(1.0,  3.0);
//                mtk::Vertex* tVertex4_4 = new Node(0.0,  3.0);
//                //------------------------------------------------------------------------------
//
//                moris::Cell< mtk::Vertex* > Name1(4); //element 1 nodes
//                Name1(0) = tVertex1_1;
//                Name1(1) = tVertex1_2;
//                Name1(2) = tVertex1_3;
//                Name1(3) = tVertex1_4;
//
//                moris::Cell< mtk::Vertex* > Name2(4);
//                Name2(0) = tVertex1_4;
//                Name2(1) = tVertex1_3;
//                Name2(2) = tVertex2_3;
//                Name2(3) = tVertex2_4;
//
//                moris::Cell< mtk::Vertex* > Name3(4);
//                Name3(0) = tVertex2_4;
//                Name3(1) = tVertex2_3;
//                Name3(2) = tVertex3_3;
//                Name3(3) = tVertex3_4;
//
//                moris::Cell< mtk::Vertex* > Name4(4);
//                Name4(0) = tVertex3_4;
//                Name4(1) = tVertex3_4;
//                Name4(2) = tVertex4_3;
//                Name4(3) = tVertex4_4;
//                //------------------------------------------------------------------------------
//
//                mtk::Cell* tElement1 = new Element(Name1);
//                mtk::Cell* tElement2 = new Element(Name2);
//                mtk::Cell* tElement3 = new Element(Name3);
//                mtk::Cell* tElement4 = new Element(Name4);
//
//                moris::Cell< mtk::Cell* > Elems(4);
//                Elems(0) = tElement1;
//                Elems(1) = tElement2;
//                Elems(2) = tElement3;
//                Elems(3) = tElement4;

            }
        }

    }/* namespace mdl */
}/* namespace moris */
