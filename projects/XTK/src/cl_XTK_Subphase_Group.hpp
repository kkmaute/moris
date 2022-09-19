/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Subphase_Group.hpp  
 * 
 */
#ifndef SRC_cl_XTK_Subphase_Group
#define SRC_cl_XTK_Subphase_Group

#include "cl_Cell.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MPI_Tools.hpp"

using namespace moris;
namespace xtk
{
    // ----------------------------------------------------------------------------------

    class Subphase_Group
    {
        // ----------------------------------------------------------------------------------
      private:
        // index for Subphase_Group
        moris_index mSubphaseGroupIndex;
        moris_id    mSubphaseGroupId = MORIS_ID_MAX;

        // owning proc
        moris_index mOwningProc = MORIS_INDEX_MAX;

        // associated B-spline Element
        moris_index mBsplineCellIndex;

        // index of SPG on associated B-spline Element
        moris_index mLocalIndex;

        // index of bulk phase the SPG lives on
        moris_index mBulkPhaseIndex = -1;

        // List of subphases in group
        moris::Cell< moris_index > mSubphaseIndicesInGroup;
        moris::Cell< moris_id > mSubphaseIdsInGroup;

        // List of ig cells in group
        moris::Cell< moris_index > mIgCellIndicesInGroup;
        bool mIgCellIndicesSet = false;

        // side ordinals of the basis (B-spline) cell through which the SPG is connected to neighboring SPGs
        moris::Cell< moris_index > mLigamentSideOrdinals;
        bool mLigamentSideOrdinalsSet = false;

        // ----------------------------------------------------------------------------------

      public:
        Subphase_Group(
                moris_index                aSubphaseGroupIndex,
                moris_index                aBsplineCellIndex,
                moris_index                aLocalSpgIndex,
                moris::Cell< moris_index > aSubphaseIndicesInGroup )
                // moris::Cell< moris_index > aSubphaseIdsInGroup )
        {
            mSubphaseGroupIndex =     aSubphaseGroupIndex;
            mBsplineCellIndex =       aBsplineCellIndex;
            mLocalIndex =             aLocalSpgIndex;
            mSubphaseIndicesInGroup = aSubphaseIndicesInGroup;
            //mSubphaseIdsInGroup =     aSubphaseIdsInGroup;
        }

        ~Subphase_Group(){}

        moris_index
        get_index() const
        {
            return mSubphaseGroupIndex;
        }

        void
        set_id( moris_id aSubphaseGroupId )
        {
            mSubphaseGroupId = aSubphaseGroupId;
        }

        moris_id
        get_id() const
        {
            return mSubphaseGroupId;
        }

        moris_index
        get_owner() const
        {
            return mOwningProc;
        }

        void
        set_owning_proc( moris_index aOwningProc )
        {
            mOwningProc = aOwningProc;
        }

        moris_index
        get_local_index() const
        {
            return mLocalIndex;
        }
        

        moris_index
        get_bspline_cell_index() const
        {
            return mBsplineCellIndex;
        }

        void 
        set_ligament_side_ordinals( moris::Cell< moris_index > aLigamentSideOrdinals )
        { 
            mLigamentSideOrdinals = aLigamentSideOrdinals;
            mLigamentSideOrdinalsSet = true;
        }

        const moris::Cell< moris_index > &
        get_ligament_side_ordinals() const
        {
            MORIS_ASSERT( mLigamentSideOrdinalsSet, "Subphase_Group::get_ligament_side_ordinals() - Side ordinals have not been set yet." );
            return mLigamentSideOrdinals;
        }

        void 
        set_ig_cell_indices( moris::Cell< moris_index > aIgCellIndicesInGroup )
        {
            MORIS_ASSERT( aIgCellIndicesInGroup.size() > 0, "Subphase_Group::set_ig_cell_indices() - passing empty list of IG cells" );
            mIgCellIndicesInGroup = aIgCellIndicesInGroup;
            mIgCellIndicesSet = true;
        }

        const moris::Cell< moris_index > &
        get_ig_cell_indices_in_group() const
        {
            MORIS_ASSERT( mIgCellIndicesSet, "Subphase_Group::get_ig_cell_indices_in_group() - IG cell indices have not been set yet." );
            return mIgCellIndicesInGroup;
        }

        const moris::Cell< moris_index > &
        get_SP_indices_in_group() const
        {
            return mSubphaseIndicesInGroup;
        }

        uint
        get_num_SPs_in_group() const
        {
            return mSubphaseIndicesInGroup.size();
        }

        void 
        set_bulk_phase( const moris_index aBulkPhaseIndex )
        { 
            mBulkPhaseIndex = aBulkPhaseIndex;
        }

        moris_index
        get_bulk_phase() const
        { 
            return mBulkPhaseIndex;
        }

        bool 
        is_subphase_ID_in_group( const moris_id aSubphaseId )
        {
            // go through Subphases in this group ...
            for( uint iSP = 0; iSP < mSubphaseIdsInGroup.size(); iSP++ )
            {
                // ... check if current SP is one looked for ...
                if( mSubphaseIdsInGroup( iSP ) == aSubphaseId )
                {
                    // ... and return that it has been found if that is the case
                    return true;
                }
            }

            // if the SP has not been found, return a false
            return false; 
        }

    }; // end: class definition

    // ----------------------------------------------------------------------------------

    struct Bspline_Mesh_Info
    {
        // counter storing the maximum number of B-spline cells and subphase groups in mesh
        moris_index mMaxBsplineCellIndex = -1;
        moris_index mMaxSpgIndex = -1;

        // store for all (Lagrange) extraction cells in which (B-spline) basis cell they live
        moris::Cell< moris_index > mExtractionCellToBsplineCell;

        // input: extraction cell index, SPG index local to extraction cell || output: list of Subphases on IP cell associated with SPG
        moris::Cell< moris::Cell< moris::Cell< moris_index > > > mExtractionCellToSubPhase; //TODO: this map remains unused

        // store which (Lagrange) extraction cells sit in a given (B-spline) basis cell
        // input: index of (B-spline) basis cell || output: list of (Lagrange) extraction cells (or their indices)
        moris::Cell< moris::Cell< mtk::Cell* > >  mExtractionCellsInBsplineCells;
        moris::Cell< moris::Cell< moris_index > > mExtractionCellsIndicesInBsplineCells;

        // store which refinement level the B-spline elements are on,
        // Note: this is needed for Ghost to ensure the whole length of a coarser B-spline cell gets penalized at a refinement boundary
        moris::Cell< uint > mBsplineCellLevels; // input: index of (B-spline) basis cell || output: refinement level of the (B-spline) basis cell

        // store which Subphase groups sit in a given (B-spline) basis cells
        // input: index of (B-spline) basis cell || output: list of indices of SPGs living on it
        moris::Cell< moris::Cell< moris_index > > mSpgIndicesInBsplineCells;

        // Subphase Groups
        moris::Cell< Subphase_Group* > mSubphaseGroups;
        moris::Cell< moris_index > mOwnedSubphaseGroupIndices; // list of SPG indices owned by the current proc
        moris::Cell< moris_index > mNotOwnedSubphaseGroupIndices; // list of SPG indices NOT owned by the current proc
        moris_id mAllocatedSpgIds = 1; // tracker for which IDs have already been taken (NOTE: this information only gets updated on proc 0)
        std::unordered_map< moris::moris_id, moris::moris_index > mSpgIdtoIndexMap; // to get the SPG index for a given SPG ID

        // SP to SPG map
        // input: SP index || output: index of SPG the SP belongs to
        moris::Cell< moris_index > mSpToSpgMap;

        // ----------------------------------------------------------------------------------

        ~Bspline_Mesh_Info()
        {
            this->delete_subphase_groups();
        }

        void
        delete_subphase_groups()
        {
            // delete subphase groups in Cell
            for ( auto iSPG : mSubphaseGroups )
            {
                delete iSPG;
            }

            // delete memory for cell
            mSubphaseGroups.clear();
        }

        // ----------------------------------------------------------------------------------

        /**
         * @brief free unused memory
         * 
         */
        void 
        finalize()
        {
            mSubphaseGroups.shrink_to_fit();
        }

        // ----------------------------------------------------------------------------------

        uint
        get_num_Bspline_cells() const
        {
            return mExtractionCellsInBsplineCells.size();
        }

        // ----------------------------------------------------------------------------------

        uint
        get_num_SPGs() const
        {
            // check that the number of SPGs returned matches the list of SPGs
            MORIS_ASSERT( (uint) mMaxSpgIndex + 1 == mSubphaseGroups.size(), "Bspline_Mesh_Info::get_num_SPGs() - mismatch between size of array of SPGs and maximum index stored" );
            
            // return value
            return (uint) mMaxSpgIndex + 1;
        }

        // ----------------------------------------------------------------------------------

        moris_id
        get_id_for_spg_index( moris_index aSubphaseGroupIndex ) const
        {
            return mSubphaseGroups( aSubphaseGroupIndex )->get_id();
        }

        // ----------------------------------------------------------------------------------

        moris_id
        get_index_for_spg_id( moris_id aSubphaseGroupId ) const
        {
            // find ID in map
            auto tIter = mSpgIdtoIndexMap.find( aSubphaseGroupId );

            // make sure the map entry makes sense
            MORIS_ASSERT( tIter != mSpgIdtoIndexMap.end(), 
                    "Bspline_Mesh_Info::get_index_for_spg_id() - Subphase Group ID not in map: %i.", aSubphaseGroupId );

            // return index associated with SPG ID
            return tIter->second;
        }

        // ----------------------------------------------------------------------------------

        moris_index 
        get_bspline_cell_index_for_extraction_cell( moris_index aExtractionCellIndex ) const
        {
            return mExtractionCellToBsplineCell( aExtractionCellIndex );
        }
        
        // ----------------------------------------------------------------------------------

        moris::Cell< moris_index > const&
        get_SPG_indices_in_bspline_cell( moris_index aBsplineCellIndex ) const
        {
            MORIS_ASSERT( (uint) aBsplineCellIndex < mExtractionCellsIndicesInBsplineCells.size(), 
                "Bspline_Mesh_Info::get_SPG_indices_in_bspline_cell() - aBsplineCellIndex out of bounds" );
            return mSpgIndicesInBsplineCells( aBsplineCellIndex );
        }

        // ----------------------------------------------------------------------------------

        moris_index
        get_local_SPG_index( moris_index aSpgIndex ) const
        {
            return mSubphaseGroups( aSpgIndex )->get_local_index();
        }

        // ----------------------------------------------------------------------------------

        uint
        get_num_SPGs_associated_with_extraction_cell( moris_index aExtractionCellIndex ) const 
        {
            // get the underlying B-spline cell's index
            moris_index tBsplineCellIndex = mExtractionCellToBsplineCell( aExtractionCellIndex );

            // get and return the number of SPGs in the B-spline cell
            return mSpgIndicesInBsplineCells( tBsplineCellIndex ).size();
        }

        // ----------------------------------------------------------------------------------

        const moris::Cell< moris_index > &
        get_SPG_indices_associated_with_extraction_cell( moris_index aExtractionCellIndex ) const
        {
            // get the underlying B-spline cell's index
            moris_index tBsplineCellIndex = mExtractionCellToBsplineCell( aExtractionCellIndex );

            // get and return the number of SPGs in the B-spline cell
            return mSpgIndicesInBsplineCells( tBsplineCellIndex );
        }

        // ----------------------------------------------------------------------------------

        void
        add_subphase_group_to_bspline_cell( 
                moris::Cell< moris_index > aSPsInGroup,           // TODO: is it a problem to pass this Cell by reference?
                moris_index                aBsplineElementIndex )
        {
            // track SPG indices and get new one
            mMaxSpgIndex++;

            // get local SPG index
            uint tLocalSpgIndex = mSpgIndicesInBsplineCells( aBsplineElementIndex ).size();

            // create a new SPG and commit it to the B-spline mesh info
            Subphase_Group* tNewSPG = new Subphase_Group( mMaxSpgIndex, aBsplineElementIndex, (moris_index) tLocalSpgIndex, aSPsInGroup );
            mSubphaseGroups.push_back( tNewSPG );

            // store that this new SPG index in the B-spline cell with the given index
            mSpgIndicesInBsplineCells( aBsplineElementIndex ).push_back( mMaxSpgIndex );
        }

        // ----------------------------------------------------------------------------------

        void
        add_ig_cell_indices_to_last_admitted_subphase_group( moris::Cell< moris_index > aIgCellIndicesInGroup ) // TODO: is it a problem to pass this Cell by reference?
        {
            // add ig cells to last admitted SPG
            mSubphaseGroups( mMaxSpgIndex )->set_ig_cell_indices( aIgCellIndicesInGroup );
        }

        // ----------------------------------------------------------------------------------

        bool
        admit_extraction_cell_group( moris::Cell< mtk::Cell * > & tExtractionCellsInBsplineCell )
        {
            // check if list L-to-B-map is initialized
            MORIS_ASSERT( &mExtractionCellToBsplineCell != nullptr && mExtractionCellToBsplineCell.size() > 0, 
                "Bspline_Mesh_Info::admit_extraction_cell_group() -mExtractionCellToBsplineCell has not been initialized." );

            // get size of extraction cell group
            moris::size_t tNumCellsInGroup = tExtractionCellsInBsplineCell.size();
            MORIS_ASSERT( tNumCellsInGroup > 0, 
                "Bspline_Mesh_Info::admit_extraction_cell_group() - empty group cannot be admitted." );

            // get first extraction Cell's index
            moris_index tFirstCellIndex = tExtractionCellsInBsplineCell( 0 )->get_index();
            MORIS_ASSERT( tFirstCellIndex >= 0, 
                "Bspline_Mesh_Info::admit_extraction_cell_group() - trying extraction to admit cell with negative index." );

            // check if Cell index is unknown
            bool tCreateNewBsplineCell = ( mExtractionCellToBsplineCell( (uint) tFirstCellIndex ) == -1 );

            if ( tCreateNewBsplineCell )
            {
                // allocate new B-spline Cell Index
                mMaxBsplineCellIndex++;

                // create List 
                moris::Cell< moris_index > tExtractionCellIndices( tNumCellsInGroup );

                // go over all cells to be admitted
                for ( moris::size_t iCell = 0; iCell < tNumCellsInGroup; iCell++)
                {              
                    // get extraction Cell's index
                    moris_index tCellIndex = tExtractionCellsInBsplineCell( iCell )->get_index();
                    MORIS_ASSERT( tCellIndex >= 0, 
                        "Bspline_Mesh_Info::admit_extraction_cell_group() - trying to admit extraction cell with negative index." );

                    // convert moris_index to uint
                    uint tPosCellIndex = (uint) tCellIndex;

                    // give B-spline cell index to extraction cell
                    mExtractionCellToBsplineCell( tPosCellIndex ) = mMaxBsplineCellIndex;

                    // give extraction cell index to B-spline cell
                    tExtractionCellIndices( iCell ) = tCellIndex;
                }

                // add to B-spline cell to Lagrange cell map
                mExtractionCellsIndicesInBsplineCells.push_back( tExtractionCellIndices );
            }

            // return whether new B-spline Cell was created
            return tCreateNewBsplineCell;
        }

        // ----------------------------------------------------------------------------------

        void
        set_ligament_side_ordinals_of_last_admitted_subphase_group( moris::Cell< bool > aActiveLigamentSideOrdinals )
        {
            // initialize list of side ordinals with correct size
            moris::Cell< moris_index > tLigamentSideOrdinals( 6 );

            // intialize counter for number of side ordinals
            uint tNumSideOrds = 0;

            // go through side ordinals, see whats active and add to list of ordinals
            for ( moris::size_t iOrd = 0; iOrd < aActiveLigamentSideOrdinals.size(); iOrd++ )
            {
                if( aActiveLigamentSideOrdinals( iOrd ) )
                {
                    tLigamentSideOrdinals( tNumSideOrds ) = iOrd;
                    tNumSideOrds++;
                }
            }

            // check that SPGs have been comitted
            MORIS_ASSERT( mMaxSpgIndex > -1, 
                "Bspline_Mesh_Info::set_ligament_side_ordinals_of_last_admitted_subphase_group() - No SPGs have been admitted yet. Unable to add ligaments." );

            // commit information to corresponding subphase group
            mSubphaseGroups( mMaxSpgIndex )->set_ligament_side_ordinals( tLigamentSideOrdinals );
        }

        // ----------------------------------------------------------------------------------

        void
        set_bulk_phase_of_last_admitted_subphase_group( const moris_index aBulkPhaseIndex )
        {
            // access last SPG admitted and set its bulk phase
            mSubphaseGroups( mMaxSpgIndex )->set_bulk_phase( aBulkPhaseIndex );
        }

        // ----------------------------------------------------------------------------------

        moris_index
        get_bulk_phase_for_subphase_group( const moris_index aSpgIndex ) const
        {
            // access last SPG admitted and set its bulk phase
            moris_index tBulkPhaseIndex = mSubphaseGroups( aSpgIndex )->get_bulk_phase();

            // check that the bulk phase for the SPG has been set
            MORIS_ASSERT( tBulkPhaseIndex > -1, 
                "Bspline_Mesh_Info::get_bulk_phase_for_subphase_group() - Bulk phase for SPG has not been set." );

            // return the index
            return tBulkPhaseIndex;
        }

        // ----------------------------------------------------------------------------------

        void
        create_SP_to_SPG_map( luint tTotNumSPs )
        {
            // initialize map with correct size
            mSpToSpgMap.resize( tTotNumSPs, -1 );

            // go over SPGs and get their SP indices, then list them in the map
            for ( luint iSPG = 0; iSPG < mSubphaseGroups.size(); iSPG++ )
            {
                // get list of SP indices in SPG
                const moris::Cell< moris_index > & tSpIndicesInGroup = mSubphaseGroups( iSPG )->get_SP_indices_in_group();

                // get the SPG's index
                const moris_index tSpgIndex = mSubphaseGroups( iSPG )->get_index();

                // loop over the SPs and set their SPG indices in the map
                for ( moris::size_t iSP = 0; iSP < tSpIndicesInGroup.size(); iSP++ )
                {
                    mSpToSpgMap( tSpIndicesInGroup( iSP ) ) = tSpgIndex;
                }
            }
        }

        // ----------------------------------------------------------------------------------

        moris_id
        allocate_subphase_group_ids( moris::size_t aNumIdstoAllocate )
        {
            // get rank of current proc and how big the MPI communicator is
            int tProcRank = moris::par_rank();
            int tProcSize = moris::par_size();

            // size_t is defined as uint here because of aNumRequested
            // Initialize gathered information outputs (information which will be scattered across processors)
            moris::Cell< moris::moris_id > aGatheredInfo;
            moris::Cell< moris::moris_id > tFirstId( 1 );
            moris::Cell< moris::moris_id > tNumIdsRequested( 1 );

            // put current processors ID request size into the Cell that will be shared across procs
            tNumIdsRequested( 0 ) = (moris::moris_id) aNumIdstoAllocate;

            // hand ID range size request to root processor
            moris::gather( tNumIdsRequested, aGatheredInfo );

            // initialize list holding the first ID in range for each processor
            moris::Cell< moris::moris_id > tProcFirstID( tProcSize );

            // Manage information only on the root processor
            if ( tProcRank == 0 )
            {
                // Loop over entities print the number of entities requested by each processor
                for ( int iProc = 0; iProc < tProcSize; ++iProc )
                {
                    // Give each processor their desired amount of IDs
                    tProcFirstID( iProc ) = mAllocatedSpgIds;

                    // update the number of already allocated SPG IDs 
                    mAllocatedSpgIds = mAllocatedSpgIds + aGatheredInfo( iProc );
                }
            }

            // on proc 0: split up the list of first IDs for every proc and send it to every other proc
            // on all procs: receive the assigned first SP ID as tFirstId
            moris::scatter( tProcFirstID, tFirstId );

            // return the first SP ID assigned 
            return tFirstId( 0 );
        }

        // ----------------------------------------------------------------------------------

        void
        assign_owned_subphase_group_ids( moris_id aFirstSpgId )
        {
            // get the number of owned 
            uint tNumOwnedSPGs = mOwnedSubphaseGroupIndices.size();

            // check that the number of owned and non-owned SPGs add up
            MORIS_ASSERT( tNumOwnedSPGs + mNotOwnedSubphaseGroupIndices.size() == this->get_num_SPGs(), 
                "Bspline_Mesh_Info::assign_owned_subphase_group_ids() - number of owned and non-owned SPGs don't add up" );

            // assign IDs to every owned SPG
            for( uint iOwnedSPG = 0; iOwnedSPG < tNumOwnedSPGs; iOwnedSPG++ )
            {
                // get the owned SPG's index
                moris_index mOwnedSpgIndex = mOwnedSubphaseGroupIndices( iOwnedSPG );

                // assign ID to SPG
                mSubphaseGroups( mOwnedSpgIndex )->set_id( aFirstSpgId );

                // increment the ID for the next SPG
                aFirstSpgId++;
            }
        }

        // ----------------------------------------------------------------------------------
    };

    // ----------------------------------------------------------------------------------

}// namespace xtk

#endif /* cl_XTK_Subphase_Group.hpp */