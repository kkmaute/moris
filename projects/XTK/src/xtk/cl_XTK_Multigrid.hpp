/*
 * cl_XTK_Multigrid.hpp
 *
 *  Created on: Mar 04, 2020
 *      Author: schmidt
 */

#ifndef XTK_SRC_XTK_CL_XTK_MULTIGRID_HPP_
#define XTK_SRC_XTK_CL_XTK_MULTIGRID_HPP_

// XTKL: Linalg Includes
#include "cl_Matrix.hpp"
#include "cl_XTK_Matrix_Base_Utilities.hpp"

// Std includes
#include <limits>

// XTKL: XTK Includes
#include "cl_Cell.hpp"
#include "cl_XTK_Model.hpp"
#include "cl_XTK_Child_Mesh.hpp"
#include "cl_XTK_Cut_Mesh.hpp"
#include "fn_Pairing.hpp"
#include "fn_equal_to.hpp"

// Mesh includes
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_Mesh_Enums.hpp"
#include "cl_XTK_Background_Mesh.hpp"
#include "cl_Mesh_Enums.hpp"

#include "fn_unique.hpp"

namespace xtk
{

class Model;

class Multigrid
{
private:

    xtk::Model * mXTKModelPtr = nullptr;

    moris::Cell< moris::Matrix< DDRMat > > mFineBasisToCoarseBasis;

public:
    Multigrid(){};

//------------------------------------------------------------------------------

    Multigrid( xtk::Model * aModel );


//------------------------------------------------------------------------------

    ~Multigrid(){};

//------------------------------------------------------------------------------

    void create_fine_to_coarse_relationship();

//------------------------------------------------------------------------------


};
}
#endif /* XTK_SRC_XTK_CL_XTK_MULTIGRID_HPP_ */
