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

    moris::Cell< moris_index > mEnrichedBasisToBackgroundBasis;

    moris::Cell< moris::Matrix< DDSMat > > mFineBasisToCoarseBasis;
    moris::Cell< moris::Matrix< DDSMat > > mCoarseBasisToFineBasis;
    moris::Cell< moris::Matrix< DDRMat > > mCoarseBasisToFineBasisWeights;

public:
    Multigrid(){};

//------------------------------------------------------------------------------

    Multigrid( xtk::Model * aModel );


//------------------------------------------------------------------------------

    ~Multigrid(){};

//------------------------------------------------------------------------------

    void create_fine_to_coarse_relationship();

//------------------------------------------------------------------------------

    void create_coarse_to_fine_relationship();

//------------------------------------------------------------------------------

    void create_coarse_to_fine_weights();

//------------------------------------------------------------------------------

    void build_enriched_coeff_to_background_coeff_map();

    void save_to_vtk( const std::string & aFilePath );

//------------------------------------------------------------------------------
#ifdef DEBUG
    void build_basis_exodus_information();
#endif
//------------------------------------------------------------------------------

#ifdef DEBUG
private:

    moris::Cell< moris::Matrix< DDRMat > > mEnrichedBasisCoords;

    moris::Matrix< DDRMat >             mEnrichedBasisLevel;
    moris::Matrix< DDRMat >             mEnrichedBasisStatus;
#endif




};
}
#endif /* XTK_SRC_XTK_CL_XTK_MULTIGRID_HPP_ */
