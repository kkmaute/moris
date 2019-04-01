/*
 * cl_XTK_Ghost_Penalization.hpp
 *
 *  Created on: Mar 26, 2019
 *      Author: doble
 */

#ifndef PROJECTS_XTK_SRC_XTK_CL_XTK_GHOST_STABILIZATION_HPP_
#define PROJECTS_XTK_SRC_XTK_CL_XTK_GHOST_STABILIZATION_HPP_


#include "cl_XTK_Model.hpp"
#include "cl_MTK_Mesh.hpp"

using namespace moris;

namespace xtk
{
class Ghost_Cell
{
public:
    Ghost_Cell():
        mCellId(MORIS_ID_MAX),
        mCellIndex(MORIS_INDEX_MAX),
        mLeftCell(0),
        mLeftCellSideOrdinal(0),
        mRightCell(0),
        mRightCellSideOrdinal(0)
    {

    }

    Ghost_Cell(
            moris::moris_index       aBackgroundFacetIndex, /*Background facet index*/
            moris::mtk::Cell const * aLeftCell,             /*Left Cell*/
            moris_index              aLeftCellSideOrdinal,  /*Side ordinal of left Cell on ghost facet*/
            moris::mtk::Cell const * aRightCell,            /*Right Cell*/
            moris_index              aRightCellSideOrdinal): /*Side ordinal of right Cell on ghost facet*/
                mCellId(MORIS_ID_MAX),
                mCellIndex(MORIS_ID_MAX),
                mBackgroundFacetIndex(aBackgroundFacetIndex),
                mLeftCell(aLeftCell),
                mLeftCellSideOrdinal(aLeftCellSideOrdinal),
                mRightCell(aRightCell),
                mRightCellSideOrdinal(aRightCellSideOrdinal)
    {

    }


    /*!
     * Set the Cell identifiers (index and id)
     *
     */
    void
    set_cell_identifiers(moris::moris_index aCellIndex,
                         moris::moris_id    aCellId)
    {
        MORIS_ASSERT(mCellId == MORIS_ID_MAX,"Cell Id has already been set");
        MORIS_ASSERT(mCellIndex == MORIS_INDEX_MAX,"Cell index has already been set");

        mCellId = aCellId;
        mCellIndex = aCellIndex;
    }

    moris::moris_id
    get_cell_id()
    {
        return mCellId;
    }

    moris::moris_index
    get_cell_index()
    {
        return mCellIndex;
    }

    /*
     * Returns the cell pointers in this interface Cell
     */
    moris::mtk::Cell const &
    get_right_cell() const
    {
        return *mRightCell;
    }

    moris::mtk::Cell const &
    get_left_cell() const
    {
        return *mLeftCell;
    }

    /*!
     * returns the side ordinals of the Cell pair on the interface
     */
    moris_index
    get_right_cell_side_ordinal() const
    {
        return mRightCellSideOrdinal;
    }

    /*!
     * returns the side ordinals of the Cell pair on the interface
     */
    moris_index
    get_left_cell_side_ordinal() const
    {
        return mLeftCellSideOrdinal;
    }


    /*!
     * return the outward facing normal relative to the provided pair index
     * @param[in] - aPairIndex
     */
    moris::Matrix<moris::F31RMat>
    get_outward_normal(moris::uint aPairIndex) const
    {
        MORIS_ERROR(0,"get_outward_normal not implemented");
        return moris::Matrix<moris::F31RMat>(0,0);
    }

    /*
     * Extract the interface Cell as a standard Cell rather than a double-sided side set type structure
     * If Cells are
     */
    moris::Matrix<moris::IndexMat>
    extract_as_standard_cell_loc_inds()
    {
        MORIS_ERROR(0,"extract_as_standard_Cell_loc_inds not implemented");
        return moris::Matrix<moris::IndexMat>(0,0);
    };

    /*!
     * Return the extracted interface Cell topology
     */
    CellTopology
    get_extracted_cell_topology()
    {
        MORIS_ERROR(0,"get_extracted_cell_topology not implemented");
        return CellTopology::INVALID;
    }



private:
    moris::moris_id                mCellId;               /*Cell id of this ghost cell*/
    moris::moris_index             mCellIndex;            /*Cell index of this ghost cell*/
    moris::moris_index             mBackgroundFacetIndex; /*Background facet index*/
    moris::mtk::Cell const *       mLeftCell;            /*Left Cells*/
    moris_index                    mLeftCellSideOrdinal; /*Side ordinal of left Cell on ghost facet*/
    moris::mtk::Cell const *       mRightCell;           /*Right Cells multiple allowed*/
    moris_index                    mRightCellSideOrdinal; /**/

};

/*
 * Data used while setting up the ghost stabilization (to not have to pass everything around in function inputs)
 */
struct Ghost_Setup_Data
{
    Cell<moris_index> mCandidateFacets;
};

class Model;

class Ghost_Stabilization
{
public:
    Ghost_Stabilization():
        mXTKModel(nullptr)
    {}

    Ghost_Stabilization( Model* aXTKModel ):
        mXTKModel(aXTKModel)
    {}

    /*!
     * Get the ghost cells
     */
    Cell<Cell<Ghost_Cell>> const &
    get_ghost_cells()
    {
        return mGhostCells;
    }

    /*!
     * Get the ghost cells
     */
    Cell<Ghost_Cell> const &
    get_ghost_cells_by_bulk_phase( moris_index aBulkPhaseIndex)
    {
        MORIS_ASSERT(aBulkPhaseIndex < (moris_index)mGhostCells.size(),"Bulk phase index out of bounds");
        return mGhostCells(aBulkPhaseIndex);
    }


    /*!
     *
     */
    void
    setup_ghost_stabilization_facets();

    /*!
     * Figure out loosely all the candidate facets which may be subject to ghost penalization
     * This basically means collect all facets attached to an intersected element
     */
    void
    identify_candidate_facets_for_ghost_penalization(Ghost_Setup_Data & aGhostSetupData);

    void
    construct_ghost_cells( moris_index        aBulkPhase,
                           Ghost_Setup_Data & aGhostSetupData);


    /*!
     * Choose whether to penalize the facet, and if we do need to penalize it
     * return the relevent information of cells on facet.
     */
    bool
    decide_whether_to_penalize_facet( moris_index aBulkPhase,
                                      moris_index aFacetIndex,
                                      Cell<moris_index>      & aCellsWithChildren,
                                      Cell<moris_index>      & aCellsWithNoChildren,
                                      Cell<Matrix<IndexMat>> & aChildMeshCellIndsOnFacet,
                                      Cell<Matrix<IndexMat>> & aChildMeshCellOnFacetOrdinals);

    /*
     * For a given bulk phase and facet index, construct the ghost cell
     */
     void
     construct_facet_ghost_cell(moris_index              aBulkPhase,
                                moris_index              aFacetIndex,
                                Cell<moris_index>      & aCellsWithChildren,
                                Cell<moris_index>      & aCellsWithNoChildren,
                                Cell<Matrix<IndexMat>> & aChildMeshCellIndsOnFacet,
                                Cell<Matrix<IndexMat>> & aChildMeshCellOnFacetOrdinals);


private:
    Model* mXTKModel; /*Pointer to the model*/

    /*!
     * Outer cell - bulk phase
     * Inner cell - ghost cells
     */
    Cell<Cell<Ghost_Cell>>     mGhostCells;
};

}



#endif /* PROJECTS_XTK_SRC_XTK_CL_XTK_GHOST_STABILIZATION_HPP_ */
