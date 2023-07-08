/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Basis_Processor.hpp
 *
 */

#pragma once

#include <unordered_set>
#include "cl_Cell.hpp"
#include "cl_Matrix.hpp"
#include "cl_Param_List.hpp"


using namespace moris;

namespace xtk
{
    class Model;
    class HMR_Helper;

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
        moris::Cell< moris::Cell< moris::real > >        mFollowerToLeaderBasisWeights;    // input: enriched BF index || output: list of leader basis corresponding to the follower basis          // input: enriched BF index || output: list of leader basis corresponding to the follower basis
        moris::Cell< moris::Cell< moris::real > >        mFollowerToLeaderBasisOwners;     // input: enriched BF index || output: list of leader basis corresponding to the follower basis

        // averaging weights for each enriched BF index
        // input: enriched BF index |  output: the list of weights corresponding to each SPGs that the basis is active
        moris::Cell< moris::Cell< real > > mAveragingWeights;

        // averaging weights for each enriched SPG indices
        // input: SPG index |  output: the list of weights corresponding to each enriched BF active on the SPG
        std::unordered_map< moris_index, moris::Cell< real > > mAveragingWeightsSPGBased;

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
        moris::Cell< moris_index > mRootSPGIds;
        moris::Cell< moris_index > mRootSPGOwner;

        // communication table
        moris::Cell< moris_id > mCommTable;
        
        // hmr accesor object for each b-spline mesh
        moris::Cell<HMR_Helper*> mHMRHelper;

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
        construct_follower_cells_using_volume( moris_index aMeshIndex );

        // ----------------------------------------------------------------------------------

        /**
         * @brief  this function uses volume info to construct follower basis / root cells
         *
         */
        void
        construct_follower_basis_using_volume( moris_index aMeshIndex );


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

        // ----------------------------------------------------------------------------------

        /**
         * @brief This function performs the nearest neighbour exchange for the root SPGs
         *
         * @param aMeshIndex
         */

        void
        perform_nearest_neighbour_exchange( moris_index aMeshIndex );

        // ----------------------------------------------------------------------------------

        /**
         * @brief prepare the identifier for the root SPG ids and list of non owned subphase groups
         * based on the processors
         *
         * @param aMeshIndex
         * @param aNotOwnedSpgsToProcs
         * @param aSPGIDs
         */

        void
        prepare_requests_for_not_owned_subphase_groups( moris_index const & aMeshIndex,
                Cell< Cell< moris_index > >&                                aNotOwnedSpgsToProcs,
                Cell< moris::Cell< moris_id > >&                            aSPGIDs );

        // ----------------------------------------------------------------------------------

        /**
         * @brief  uses the information from prepare_requests_for_not_owned_subphase_groups to prepare
         * the data the needs to be sent to other processors
         *
         * @param aMeshIndex
         * @param aSendSubphaseGroupRootIds
         * @param aReceivedSPGIds
         */

        void
        prepare_answers_for_owned_subphase_groups( moris_index const & aMeshIndex,
                moris::Cell< moris::Cell< moris_id > >&                aSendSubphaseGroupRootIds,
                moris::Cell< moris::Cell< moris_id > >&                aSendSubphaseGroupRootOwners,
                moris::Cell< moris::Cell< moris_id > > const &         aReceivedSPGIds );

        // ----------------------------------------------------------------------------------

        /**
         * @brief sets the SPG root IDs received from neighboring processors
         *
         * @param aReceivedSubphaseGroupRootIds
         * @param aReceivedSubphaseGroupRootOwners
         * @param aNotOwnedSpgsToProcs
         */

        void
        handle_requested_subphase_groups_answers(
                moris::Cell< moris::Cell< moris_id > >const &        aReceivedSubphaseGroupRootIds,
                moris::Cell< moris::Cell< moris_id > > const &        aReceivedSubphaseGroupRootOwners,
                moris::Cell< moris::Cell< moris_id > > const & aNotOwnedSpgsToProcs );

        //-----------------------------------------------------------------------------------

        /**
         * @brief Determines if the cell aggregation process should be stopped
         * It stops when all SPGs are assigned a root SPG Id
         *
         * @return true
         * @return false
         */

        bool
        determine_stopping_criteria();

        // ----------------------------------------------------------------------------------
        
        /**
         * @brief constructs a mCommTable based on the layout of the extended and root cells
         * 
         */
        void
        construct_comm_table();

        // ----------------------------------------------------------------------------------

        /**
         * @brief paralle call to detrmine the root spg ids for the not owned spgs
         * 
         * @param aMeshIndex  aDiscretization mesh index
         * @param aNotOwnedSpgsToProcs list of SPGS indices for each neighbouring processor ( identifier that needs to be process the data upon recieving from neighbouring processors )
         * @param aSPGIDs SPG ids that data is needed for to be able to extende the basis ( identifier that needs to be genereated to request from neighbouring processors )
         */


        void
        prepare_requests_for_not_owned_root_spg_projections( moris_index aMeshIndex,
                Cell< Cell< moris_index > >&                             aNotOwnedSpgsToProcs,
                Cell< moris::Cell< moris_id > >&                         aSPGIDs );

        // ----------------------------------------------------------------------------------

        /**
         * @brief the function that after recieving the SPG ids from the neighbouring processors prepares answers to be sent
         * 
         * @param aMeshIndex aDiscretization mesh index
         * @param aSendSubphaseGroupRootBsplineElementIds input: procesor index || output: Bspline global element ids from hmr that was requested
         * @param aSendSubphaseGroupRootBsplineBasisIds input: procesor index || output: enriched bspline basis ids present on the previous b-spline element
         * @param aSendSubphaseGroupRootBsplineBasisOwners input: procesor index || output: enriched bspline basis owner present on the previous b-spline element
         * @param aReceivedSPGIds list of SPG ids that the Bspline element IDs corropspond to ( identifier that will be used to prepare the data)
         */

        void
        prepare_answers_for_owned_root_spg_projections( moris_index const & aMeshIndex,
                moris::Cell< moris::Cell< moris_id > >&                     aSendSubphaseGroupRootBsplineElementIds,
                moris::Cell< moris::Cell< moris_id > >&                     aSendSubphaseGroupRootBsplineBasisIds,
                moris::Cell< moris::Cell< moris_id > >&                     aSendSubphaseGroupRootBsplineBasisOwners,
                moris::Cell< moris::Cell< moris_id > > const &              aReceivedSPGIds );

        // ----------------------------------------------------------------------------------
        
        /**
         * @brief after recieving the data from the neighbouring processors this function sets the data
         * 
         * @param aMeshIndex aDiscretization mesh index
         * @param aReceivedSubphaseGroupRootBsplineElementIds input: procesor index || output: Bspline global element ids from hmr that was requested
         * @param aReceivedSubphaseGroupRootBsplineBasisIds input: procesor index || output: enriched bspline basis ids present on the previous b-spline element
         * @param aReceivedSubphaseGroupRootBsplineBasisOwners input: procesor index || output: enriched bspline basis owner present on the previous b-spline element
         * @param aNotOwnedSpgsToProcs list of SPG indices that the Bspline element IDs corropspond to ( identifier that was generated for in-processor data)
         */
        void
        handle_requested_root_spg_projections( moris_index aMeshIndex,
                moris::Cell< moris::Cell< moris_id > >const &        aReceivedSubphaseGroupRootBsplineElementIds,
                moris::Cell< moris::Cell< moris_id > >const &        aReceivedSubphaseGroupRootBsplineBasisIds,
                moris::Cell< moris::Cell< moris_id > >const &        aReceivedSubphaseGroupRootBsplineBasisOwners,
                moris::Cell< moris::Cell< moris_id > > const & aNotOwnedSpgsToProcs );

        // ----------------------------------------------------------------------------------
        
        /**
         * @brief the parallel communication function that does cell to cell extension for the neighbours
         * 
         * @param aMeshIndex 
         */
        void
        construct_follower_to_leader_basis_weights_indices_neighbours( moris_index aMeshIndex );

        //-----------------------------------------------------------------------------------

        /**
         * @brief constructs cell extension for the extended cell based on the root cell information tat was recieved from the neighbouring processors
         * 
         * @param aSubphaseGroupIndex the subphase group index of the extdended cell
         * @param aRootBsplineId the bpsline global domain id of the root cell
         * @param aRootBsplineBasisId  the enriched basis that is present on the root cell
         * @param aRootBsplineBasisOwner the enriched basis owner that is present on the root cell
         * @param aStartRange location of the first basis id that corresponds to the root cell
         * @param aMeshIndex  aDiscretization mesh index
         */
        void
        construct_follower_to_leader_basis_relationship_for_spg(
                moris_index              aSubphaseGroupIndex,
                moris_id                 aRootBsplineId,
                moris::Cell< moris_id >const & aRootBsplineBasisId,
                moris::Cell< moris_id >const & aRootBsplineBasisOwner,
                uint                     aStartRange,
                uint                     aMeshIndex);

        //-----------------------------------------------------------------------------------

        /**
         * @brief constuctrs cell to cell extension when both the extended and root cell is owned by the current processor
         * 
         * @param aMeshIndex 
         */
        void
        construct_follower_to_leader_basis_weights_indices_mine( moris_index aMeshIndex );

        //-----------------------------------------------------------------------------------

        /**
         * @brief function that ensures that the replaced basis on each processor is communicated to the neighbouring processors
         * such that the replacement is consistent
         * 
         * @param aMeshIndex 
         */

        void
        communicate_shared_basis( moris_index aMeshIndex );

        //-----------------------------------------------------------------------------------

        /**
         * @brief function that prepares the basis data to be sent to the neighbouring processors
         * 
         * @param aMeshIndex aDiscretization mesh index
         * @param aBasisIndexToProcs identifies for the basis index that the current processor need to modify
         * @param tSendBasisIds identifier for the basis that needs to be communicated to this processor
         */
        void
        prepare_requests_for_shared_follower_basis( moris_index aMeshIndex,
                Cell< Cell< moris_index > >&                    aBasisIndexToProcs,
                Cell< moris::Cell< moris_id > >&                tSendBasisIds );

        //-----------------------------------------------------------------------------------

        /**
         * @brief function that fills in the basis data based on the tSendBasisIds identifer that was recieved from the neighbouring processors
         * 
         * @param aMeshIndex discretization mesh index
         * @param tSendFollowerToLeaderBasisIds 
         * @param tSendFollowerToLeaderBasisOwners 
         * @param tSendFollowerToLeaderBasisWeights 
         * @param tSendFollowerToLeaderOffset 
         * @param aReceivedBasisIds 
         */
        void
        prepare_answers_for_follower_shared_basis(
                moris_index const &                            aMeshIndex,
                moris::Cell< moris::Cell< moris_id > >&        tSendFollowerToLeaderBasisIds,
                moris::Cell< moris::Cell< moris_id > >&        tSendFollowerToLeaderBasisOwners,
                moris::Cell< moris::Cell< real > >&            tSendFollowerToLeaderBasisWeights,
                Cell< moris::Cell< moris_id > >&               tSendFollowerToLeaderOffset,
                moris::Cell< moris::Cell< moris_id > > const & aReceivedBasisIds );

        //-----------------------------------------------------------------------------------

        /**
         * @brief function that handles the recieved basis data and changes the wights and owner and thier replacement based on the recieved data
         * 
         * @param aMeshIndex 
         * @param tReceivedFollwerToLeaderBasisIds 
         * @param tReceivedFollwerToLeaderBasisOwners 
         * @param tReceivedFollwerToLeaderBasisWeights 
         * @param tSendFollowerToLeaderOffset 
         * @param tBasisIndexToProcs 
         */
        void
        handle_requested_shared_follower_basis(
                moris_index                                    aMeshIndex,
                moris::Cell< moris::Cell< moris_id > > const & tReceivedFollwerToLeaderBasisIds,
                moris::Cell< moris::Cell< moris_id > > const & tReceivedFollwerToLeaderBasisOwners,
                moris::Cell< moris::Cell< real > > const &     tReceivedFollwerToLeaderBasisWeights,
                Cell< moris::Cell< moris_id > > const &        tSendFollowerToLeaderOffset,
                moris::Cell< moris::Cell< moris_id > > const & tBasisIndexToProcs );
        
        //-----------------------------------------------------------------------------------

        /**
         * @brief updates the commtable based on the mCommTable and passes it to the cut integration meesh
         * 
         */
        void
        update_comm_table(); 
    };

}    // namespace xtk