/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Multigrid.hpp
 *
 */

#ifndef XTK_SRC_XTK_CL_XTK_MULTIGRID_HPP_
#define XTK_SRC_XTK_CL_XTK_MULTIGRID_HPP_

// XTKL: Linalg Includes
#include "cl_Matrix.hpp"
#include "cl_XTK_Matrix_Base_Utilities.hpp"

// Std includes
#include <limits>

// XTKL: XTK Includes
#include "cl_Vector.hpp"
#include "cl_XTK_Model.hpp"
#include "cl_XTK_Child_Mesh.hpp"
#include "cl_XTK_Cut_Mesh.hpp"
#include "fn_Pairing.hpp"
#include "fn_equal_to.hpp"

// Mesh includes
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_XTK_Background_Mesh.hpp"

#include "fn_unique.hpp"

namespace moris::xtk
{

    class Model;

    class Multigrid
    {
      private:
        xtk::Model* mXTKModelPtr = nullptr;

        moris_index mMeshIndex;

        Vector< moris_index > mEnrichedBasisToBackgroundBasis;

        Vector< Matrix< DDSMat > > mFineBasisToCoarseBasis;
        Vector< Matrix< DDSMat > > mCoarseBasisToFineBasis;
        Vector< Matrix< DDRMat > > mCoarseBasisToFineBasisWeights;

        Matrix< DDRMat > mEnrichedBasisLevel;

        uint mNumBasis = MORIS_UINT_MAX;

      public:
        Multigrid(){};

        //------------------------------------------------------------------------------

        Multigrid( xtk::Model* aModel );

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

        //------------------------------------------------------------------------------

        uint
        get_num_interpolations()
        {
            return 1;
        };

        //-------------------------------------------------------------------------------

        uint
        get_max_level( const moris_index aInterpolationIndex )
        {
            moris::mtk::Interpolation_Mesh& tInterpolationMesh = mXTKModelPtr->get_background_mesh();

            return tInterpolationMesh.get_max_level( aInterpolationIndex );
        };

        //-------------------------------------------------------------------------------

        uint
        get_num_basis( const moris_index aInterpolationIndex )
        {
            return mNumBasis;
        }

        //-------------------------------------------------------------------------------

        uint
        get_basis_level( const moris_index aInterpolationIndex,
                const moris_index          aBasisIndex )
        {
            return mEnrichedBasisLevel( aBasisIndex );
        }

        //-------------------------------------------------------------------------------

        uint
        get_num_coarse_basis_of_basis( const moris_index aInterpolationIndex,
                const moris_index                        aBasisIndex )
        {
            return mFineBasisToCoarseBasis( aBasisIndex ).numel();
        }

        //-------------------------------------------------------------------------------

        uint
        get_coarse_basis_index_of_basis( const moris_index aInterpolationIndex,
                const moris_index                          aBasisIndex,
                const moris_index                          aCoarseParentIndexForBasis )
        {
            return mFineBasisToCoarseBasis( aBasisIndex )( aCoarseParentIndexForBasis );
        }

        //-------------------------------------------------------------------------------

        Matrix< DDSMat >
        get_fine_basis_inds_of_basis( const moris_index aInterpolationIndex,
                const moris_index                       aBasisIndex )
        {
            return mCoarseBasisToFineBasis( aBasisIndex );
        }

        //-------------------------------------------------------------------------------

        Matrix< DDRMat >
        get_fine_basis_weights_of_basis( const moris_index aInterpolationIndex,
                const moris_index                          aBasisIndex )
        {
            return mCoarseBasisToFineBasisWeights( aBasisIndex );
        }

        //-------------------------------------------------------------------------------

#ifdef MORIS_HAVE_DEBUG
        void save_to_vtk( const std::string& aFilePath );
#endif
        //------------------------------------------------------------------------------

        void build_basis_exodus_information( const std::string& aName );

        //------------------------------------------------------------------------------

#ifdef MORIS_HAVE_DEBUG

      private:
        Matrix< DDRMat > mEnrichedBasisCoords;

        Matrix< DDRMat > mEnrichedBasisStatus;
#endif
    };
}    // namespace moris::xtk
#endif /* XTK_SRC_XTK_CL_XTK_MULTIGRID_HPP_ */
