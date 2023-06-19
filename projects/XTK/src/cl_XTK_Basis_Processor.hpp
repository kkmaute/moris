/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Basis_Processor.hpp
 *
 */

#include <unordered_set>
#include "cl_Cell.hpp"
#include "cl_Matrix.hpp"
#include "cl_Param_List.hpp"


#pragma once

using namespace moris;

namespace xtk
{
    class Model;

    enum Dependency_Criteria
    {
        VOLUME,
        BASIS_SUPPORT,
        BASIS_INCLUSION,
        END_ENUM,
    };

    struct Basis_Data
    {
        // basis agglomeration information
        // input: enriched BF index || output: 1 if it is a follower basis , 0: leader basis
        moris::Cell< moris::moris_index > mFollowerBasis;

        // The following maps are required to construct the follower basis in terms of leader basis
        moris::Cell< moris::Cell< moris::moris_index > > mFollowerToLeaderBasis;           // input: enriched BF index || output: list of leader basis corresponding to the follower basis
        moris::Cell< moris::Cell< moris::real > >        mFollowerToLeaderBasisWeights;    // input: enriched BF index || output: list of leader basis corresponding to the follower basis

        // averaging weights for each enriched BF index
        // input: enriched BF index |  output: the list of weights corresponding to each SPGs that the basis is active
        moris::Cell< moris::Cell< real > > mAveragingWeights;

        // Subspaces that will be used for blocking in Schwartz method
        // input: # of block || output: list of enriched BFs  on the block
        moris::Cell< moris::Cell< moris::moris_index > > mBasisInSubspace;

        // input: enriched BF index || output: Root SPG of that basis ( used for cell agglomeration )
        moris::Cell< moris::moris_index > mFollowerBasisOwningCell;
    };

    class Basis_Processor
    {
      private:
        // model pointer to access data
        xtk::Model* mXTKModelPtr = nullptr;

        // B-spline mesh indices
        moris::Matrix< moris::IndexMat > mMeshIndices;

        moris::ParameterList* mParameterList;

        // spatial dimension
        uint mSpatialDim;

        // struct containing basis extension info
        // input: DMI
        moris::Cell< Basis_Data > mBasisData;

        // input: SPG index || output: Root SPG index of the SPG
        moris::Cell< moris_index > mRootSPGIndex;


      public:
        // ----------------------------------------------------------------------------------

        /**
         * @brief Construct a new Basis_Processor object
         *
         * @param aXTKModelPtr
         */
        explicit Basis_Processor( xtk::Model* aXTKModelPtr );

        // ----------------------------------------------------------------------------------

        /**
         * @brief Destroy the Basis_Processor object
         *
         */

        ~Basis_Processor();

        // ----------------------------------------------------------------------------------

        /**
         * @brief This function generates basis block for schwartz smoothing methods using volume of the cells
         *
         */

        void
        construct_volumetric_clustering_of_basis();

        // ----------------------------------------------------------------------------------

        /**
         * @brief This function generates basis block for schwartz smoothing methods using the assigned criteria
         *
         * @param aCriteria
         */

        void
        construct_basis_in_subspace( enum Dependency_Criteria aCriteria );

        // ----------------------------------------------------------------------------------


        /**
         * @brief This function generates basis block for schwartz smoothing methods using the assigned criteria
         *
         * @param aCriteria
         */

        void
        construct_follower_basis_list( enum Dependency_Criteria aCriteria, moris_index aMeshIndex );

        // ----------------------------------------------------------------------------------

        /**
         * @brief this function uses basis support information(without computing volume) to construct follower basis
         * This is an optimized version of the volume since
         *
         */


        void
        construct_follower_basis_using_basis_support( moris_index aMeshIndex );

        // ----------------------------------------------------------------------------------

        /**
         * @brief  this function uses volume info to construct follower basis / root cells
         *
         */
        void
        construct_follower_basis_using_volume( moris_index aMeshIndex );


        // void
        // construct_cell_aggregates();

        // ----------------------------------------------------------------------------------

        /**
         * @brief adds the root spg index as a field to xtk enriched ip and ig meshes
         *
         */

        void
        visualize_cell_aggregates( moris_index aMeshIndex );

        // ----------------------------------------------------------------------------------

        /**
         * @brief performer function ( see the function definition for more info)
         *
         */
        void
        perform_basis_extention();

        // ----------------------------------------------------------------------------------

        /**
         * @brief performer function ( see the function definition for more info )
         *
         */
        void
        perform_cell_agglomeration();

        // ----------------------------------------------------------------------------------

        /**
         * @brief construct ownership between follower basis to the owning root spg
         *
         */
        void
        construct_dependent_basis_to_root_map( moris_index aMeshIndex );

        // ----------------------------------------------------------------------------------

        /**
         * @brief  populate the member data mAveragingWeights
         *
         */

        void
        compute_averaging_weights( moris_index aMeshIndex );

        // ----------------------------------------------------------------------------------

        /**
         * @brief This function computes the volume of the treshhold volume that will be used to determine
         *  follower and leader b-spline elements
         *
         */
        real
        compute_threshold_volume( moris_index aMeshIndex );

        // ----------------------------------------------------------------------------------

        /**
         * @brief build the data structure mRootSPGIndex and associated root spgs the cut spgs
         *
         */

        void construct_cell_association( moris_index aMeshIndex );

        // ----------------------------------------------------------------------------------

        /**
         * @brief build the data structure mFollowerToLeaderBasis and mFollowerToLeaderBasisWeights
         *
         */

        void
        construct_follower_to_leader_basis_weights_indices( moris_index aMeshIndex );

        // ----------------------------------------------------------------------------------

        /**
         * @brief replaces the vertex enrichment data(basis ids and wights and indices with the given data)
         *
         */

        void
        replace_t_matrices( moris_index aMeshIndex );

        // ----------------------------------------------------------------------------------

        /**
         * @brief construct cell aggregation method
         *
         */
        void
        construct_cell_aggregates( moris_index aMeshIndex );

        // ----------------------------------------------------------------------------------

        /**
         * @brief over-ride the t matrices, this function is exclusive for the cell agglomeration method
         *
         * @param aMeshIndex
         */

        void
        override_t_matrices( moris_index aMeshIndex );
    };

}    // namespace xtk