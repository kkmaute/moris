
#include "cl_Stopwatch.hpp" //CHR/src

#include "MTK_Tools.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh_Manager.hpp"                    //MTK/src

#include "cl_FEM_Node_Base.hpp"               //FEM/INT/src
#include "cl_FEM_Node.hpp"               //FEM/INT/src
#include "cl_FEM_Enums.hpp"               //FEM/INT/src

#include "cl_FEM_Element_Bulk.hpp"               //FEM/INT/src
#include "cl_FEM_IWG_Factory.hpp"
#include "cl_FEM_Element_Factory.hpp"
#include "cl_FEM_Set.hpp"

#include "cl_MDL_Mesh_Model_Helper.hpp"

namespace moris
{
    namespace mdl
    {
//------------------------------------------------------------------------------

    Mesh_Model_Helper::Mesh_Model_Helper( mtk::Mesh_Manager * aMesh ) : mMesh( aMesh )
    {
        // get mesh

        // get colors for blocks
        // colormap build on first entry of list
        // create fem::blocks -- input
        //


        // get blocks
        // get vertices on blocks
        // make them unique
        // build nodes based on unique list
        //

    }

//------------------------------------------------------------------------------

    Mesh_Model_Helper::~Mesh_Model_Helper()
    {}

//------------------------------------------------------------------------------

    void Mesh_Model_Helper::compute_unique_vertex_lists()
    {}

//------------------------------------------------------------------------------

    void Mesh_Model_Helper::create_nodes()
    {}


    } /* namespace mdl */
} /* namespace moris */
