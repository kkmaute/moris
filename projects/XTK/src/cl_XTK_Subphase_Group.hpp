/**
 * cl_XTK_Subphase_Group.hpp  
 * 
 *  Created on: Jan  18, 2022 
 *      Author: Nils Wunsch
 */
#ifndef SRC_cl_XTK_Subphase_Group
#define SRC_cl_XTK_Subphase_Group

#include "cl_Cell.hpp"
#include "cl_MTK_Cell.hpp"

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

        // associated Bspline Element
        moris_index mBsplineCellIndex;

        // index of bulk phase the SPG lives on
        moris_index mBulkPhaseIndex;

        // List of subphases in group
        moris::Cell< moris_index > mSubphaseIndicesInGroup;

        // side ordinals of the basis (B-spline) cell through which the SPG is connected to neighboring SPGs
        moris::Cell< moris_index > mLigamentSideOrdinals;
        bool mLigamentSideOrdinalsSet = false;

        // ----------------------------------------------------------------------------------

      public:
        Subphase_Group(
                moris_index                aSubphaseGroupIndex,
                moris_index                aBsplineCellIndex,
                moris::Cell< moris_index > aSubphaseIndicesInGroup )
        {
            mSubphaseGroupIndex =     aSubphaseGroupIndex;
            mBsplineCellIndex =       aBsplineCellIndex;
            mSubphaseIndicesInGroup = aSubphaseIndicesInGroup;
        }

        ~Subphase_Group();

        moris_index
        get_index() const
        {
            return mSubphaseGroupIndex;
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
    };

    // ----------------------------------------------------------------------------------

    struct Bspline_Mesh_Info
    {
        // counter storing the maximum number of B-spline cells and subphase groups in mesh
        moris_index mMaxBsplineCellIndex = -1;
        moris_index mMaxSpgIndex = -1;

        // store for all (Lagrange) extraction cells in which (B-spline) basis cell they live
        moris::Cell< moris_index > mExtractionCellToBsplineCell;

        // store which (Lagrange) extraction cells sit in a given (B-spline) basis cell
        // input: index of (B-spline) basis cell || output: list of (Lagrange) extraction cells (or their indices)
        moris::Cell< moris::Cell< mtk::Cell* > >  mExtractionCellsInBsplineCells; // TODO: needed?
        moris::Cell< moris::Cell< moris_index > > mExtractionCellsIndicesInBsplineCells;

        // store which Subphase groups sit in a given (B-spline) basis cells
        // input: index of (B-spline) basis cell || output: list of SPGs living in it (or their indices)
        moris::Cell< moris::Cell< moris_index > > mSPGsIndicesInBsplineCells;

        // Subphase Groups
        moris::Cell< Subphase_Group * > mSubphaseGroups;

        // // store which bulk-phases are present in each basis (B-spline) cell
        // // moris::Cell< moris::Cell< moris_index > > mBulkPhasesInBsplineCells; // TODO: needed?

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

        moris_index
        get_num_Bspline_cells()
        {
            return mMaxBsplineCellIndex + 1;
        }

        moris_index
        get_num_SPGs()
        {
            return mMaxSpgIndex + 1;
        }

        // ----------------------------------------------------------------------------------

        void
        add_subphase_group_to_last_admitted_bspline_cell( moris::Cell< moris_index > & aSPsInGroup )
        {
            // check  
            MORIS_ASSERT( mMaxBsplineCellIndex > -1, 
                "Bspline_Mesh_Info::add_subphase_group_to_last_admitted_bspline_cell() - no B-spline cells have been admitted yet. Unable to add SPG." );

            // track SPG indices and get new one
            mMaxSpgIndex++;

            // create and commit a new SPG to the B-spline mesh info
            // FIXME: this is way of allocating SPGs is painfully expensive for big problems, needs to be made more efficient
            Subphase_Group * tNewSPG = new Subphase_Group( mMaxSpgIndex, mMaxBsplineCellIndex, aSPsInGroup );
            mSubphaseGroups.push_back( tNewSPG );
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
    };

    // ----------------------------------------------------------------------------------

}// namespace xtk

#endif /* cl_XTK_Subphase_Group.hpp */