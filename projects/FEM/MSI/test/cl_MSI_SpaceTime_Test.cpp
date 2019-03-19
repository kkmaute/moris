/*
 * cl_MSI_Multigrid.cpp
 *
 *  Created on: Jul 14, 2018
 *      Author: schmidt
 */

#include "catch.hpp"
#include "typedefs.hpp"
#include "cl_Map.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_equal_to.hpp"

#include "cl_Communication_Tools.hpp"
#include "cl_Communication_Manager.hpp"

#define protected public
#define private   public
#include "cl_MSI_Equation_Object.hpp"
#include "cl_MSI_Model_Solver_Interface.hpp"
#undef protected
#undef private

#include "cl_MTK_Vertex.hpp"    //MTK
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh.hpp"

#include "cl_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Tools.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Scalar_Field_Info.hpp"
#include "cl_MTK_Mesh.hpp"

#include "cl_FEM_NodeProxy.hpp"                //FEM/INT/src
#include "cl_FEM_ElementProxy.hpp"             //FEM/INT/src
#include "cl_FEM_Node_Base.hpp"                //FEM/INT/src
#include "cl_FEM_Element_Factory.hpp"          //FEM/INT/src
#include "cl_FEM_IWG_Factory.hpp"              //FEM/INT/src

namespace moris
{
namespace MSI
{
TEST_CASE( "MSI_SPace_Time", "[moris],[MSI],[MSI_Space_Time]" )
{

}/* TEST_CASE */

}
}


