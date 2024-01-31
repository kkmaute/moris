/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_HMR_Helper.hpp
 *
 */

#pragma once

// includes from MORIS
#include "moris_typedefs.hpp"    // COR/src
#include "cl_Cell.hpp"
#include "cl_HMR_Background_Element_Base.hpp"


// forward declaration of the hmr mesh
namespace moris::hmr
{
    class Lagrange_Mesh_Base;
}

using namespace moris;

namespace moris::mtk
{
    class Cell;
}
namespace xtk
{
    // forward declaration of the class model
    class Model;
    class Enrichment_Data;

    class HMR_Helper
    {
      private:
        // model pointer to access data
        xtk::Model* mXTKModelPtr = nullptr;

        // pointer to the hmr mesh
        hmr::Lagrange_Mesh_Base* mHMRLagrangeMesh = nullptr;

        // Cell of enriched basis IDs/Indices, it gets overwritten for each call
        moris::Cell< moris_id > mEnrichedBasisIDs = {};

        // Cell of enriched basis IDs for each Bspline cell and refinement level
        moris::Cell< moris_id > mEnrichedBasisOwners = {};

        // number of basis per element
        uint mNumberOfBasis = 0;

        // is equal to polynomial degree +
        uint mNumberOfNodesPerDimension = 0;

        // Bpsline mesh index
        uint mBsplineMeshIndex = 0;

        Matrix< DDRMat >                mL2ProjectionMatrix = { {} };
        moris::Cell< Matrix< DDRMat > > mMatrices1D;

        // local order of the basis in the bspline element
        Matrix< DDUMat > mBasisIndex = { {} };

        // order if the Bspline basis functions
        uint mBSplineOrder = 0;

        // pointer to enrichment data object ( not owned by this class )
        const Enrichment_Data* mEnrichmentData;

        uint mSpatialDimension = 0;

      public:
        //------------------------------------------------------------------------------------

        /**
         * @brief Construct a new hmr helper object
         *
         * @param aXTKModel
         */

        HMR_Helper( xtk::Model* aXTKModel, moris_index aMeshIndex );

        //------------------------------------------------------------------------------------

        /**
         * @brief Destroy the hmr helper object
         *
         */

        ~HMR_Helper();

        //------------------------------------------------------------------------------------

        /**
         * @brief Get the global domain id of bspline cell ( global domain id is a unique id based on the ijk and refinement level of the cell)
         *
         * @param aCell parent pointer to the bspline cell
         * @return luint global domain id
         */

        luint
        get_global_domain_id_of_cell( const mtk::Cell* aCell );


        //------------------------------------------------------------------------------------

        /**
         * @brief Get the enriched basis id of cell object
         *
         * @param aCell a parent pointer to the hmr element
         * @param aSPGIndex subphase group index of the cell
         * @return moris::Cell< moris_id >& enirched basis ids ordered in the order of hmr nodes
         */

        moris::Cell< moris_id > const&
        get_enriched_basis_id_of_cell( const mtk::Cell* aCell, moris_index aSPGIndex );

        //------------------------------------------------------------------------------------

        /**
         * @brief Get the enriched basis indicies of cell ob
         *
         * @param aCell a parent pointer to the hmr element
         * @param aSPGIndex subphase group index of the cell
         * @return moris::Cell< moris_id >& enirched basis indices ordered in the order of hmr nodes
         */

        moris::Cell< moris_id >&
        get_enriched_basis_indicies_of_cell( const mtk::Cell* aCell, moris_index aSPGIndex );

        //------------------------------------------------------------------------------------

        /**
         * @brief Get the number of basis per element in hmr
         *
         * @return uint number of bassis per element
         */

        uint
        get_number_of_bases_per_element();

        //------------------------------------------------------------------------------------

        /**
         * @brief Get the bg basis indices of cell object
         *
         * @param aCell  parent pointer to the hmr element
         * @return moris::Cell< moris_id >& get the background ( unenriched ) basis indices of this cell
         */

        moris::Cell< moris_id > const &
        get_bg_basis_indices_of_cell( const mtk::Cell* aCell );

        //------------------------------------------------------------------------------------

        /**
         * @brief initiliaze the ordereding of the basis based on the order and dimension
         *
         */
        void
        init_basis_index();

        //------------------------------------------------------------------------------------

        /**
         * @brief Create a background element object
         *
         * @return hmr::Background_Element_Base*
         */

        hmr::Background_Element_Base*
        create_background_element();

        //------------------------------------------------------------------------------------

        /**
         * @brief Get the extention matrix 1d object
         *
         * @param aShift i,j,k shift of the elements w.r.t the root element
         * @param aExtentionMatrix Extnesion matrix relating the exneted and original basis
         */
        void
        get_extention_matrix_1d( real aShift, Matrix< DDRMat >& aExtentionMatrix );

        //------------------------------------------------------------------------------------

        /**
         * @brief Get the l2 projection matrix object
         *
         * @param aExtendedCell parent pointer to the extended cell
         * @param aRootBsplineId global domain id of the root cell
         * @return Matrix< DDRMat >& extension matrix realting basis of the root and extended cell
         * NOTE: this function assumes that both the extended cell and the root cell are in the same refinement level
         */

        Matrix< DDRMat > const &
        get_l2_projection_matrix( const mtk::Cell* aExtendedCell, moris_id aRootBsplineId );

        //------------------------------------------------------------------------------------

        /**
         * @brief Get the l2 projection matrix object, this function is used when both cells are within the same processor domain
         *
         * @param aExtendedCell parent pointer to the extended cell
         * @param aRootCell parent pointer to the root cell
         * @return Matrix< DDRMat >& extension matrix realting basis of the root and extended cell
         */

        Matrix< DDRMat > const &
        get_l2_projection_matrix( const mtk::Cell* aExtendedCell, const mtk::Cell* aRootCell );

        //------------------------------------------------------------------------------------

        /**
         * @brief Get the enriched basis owners of cell object this function is used in conjunction with the get_enriched_basis_id_of_cell function
         *  
         * @return moris::Cell<moris_id>& the owners of the basis in the order of hmr nodes
         * 
         * This is done this way to save memory, since the size enriched basis owners are the same size for all the subphase groups, so this cell gets overwritten
         */

        moris::Cell< moris_id > const &
        get_enriched_basis_owners_of_cell()
        {
            return mEnrichedBasisOwners;
        }

        //------------------------------------------------------------------------------------

        /**
         * @brief Get the ijk bspline cell object, the ijk is local
         *
         * @param aCell parent pointer to the bspline cell
         * @return const luint* array of ijk of the cell( size of 2 or 3 depending on the dimension )
         */

        const luint*
        get_ijk_bspline_cell( const mtk::Cell* aCell );
    };
}    // namespace xtk