/*
 * cl_MTK_Integration_Mesh.hpp
 *
 *  Created on: Apr 15, 2019
 *      Author: doble
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_INTEGRATION_MESH_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_INTEGRATION_MESH_HPP_

#include "cl_MTK_Mesh_Core.hpp"
#include "cl_MTK_Cell_Cluster.hpp"
#include "cl_MTK_Side_Cluster.hpp"
#include "cl_MTK_Double_Side_Cluster.hpp"
#include "assert.hpp"
#include "cl_Matrix.hpp"
#include "cl_MTK_Block.hpp"

namespace moris
{
namespace mtk
{
class Integration_Mesh : public virtual Mesh
{
protected:
    moris::Cell< moris::mtk::Set * > mListofBlocks;
    moris::Cell< moris::mtk::Set * > mListofSideSets;
    moris::Cell< moris::mtk::Set * > mListofDoubleSideSets;

public:
    Integration_Mesh(){};
    // Functions only valid for integration meshes

    //##############################################
    // Cell Cluster Access
    //##############################################

    /*
     * Get a cell cluster related to an interpolation
     * cell
     */
    virtual
    Cell_Cluster const &
    get_cell_cluster(Cell const & aInterpCell) const = 0;

    // ----------------------------------------------------------------------------

    /*
     * Get block set names
     */
    virtual
    moris::Cell<std::string>
    get_block_set_names() const = 0;

    /*!
     * Returns the label
     */
    virtual
    std::string
    get_block_set_label(moris_index aBlockSetOrdinal) const
    {
        MORIS_ERROR(0, "get_block_set_label has no default implementation");
        return "ERROR";
    }

    // ----------------------------------------------------------------------------

    /*!
     * Returns the index given a label
     */
    virtual
    moris_index
    get_block_set_index(std::string aBlockSetLabel) const
    {
        MORIS_ERROR(0, "get_block_set_index has no default implementation");
        return MORIS_INDEX_MAX;
    }

    // ----------------------------------------------------------------------------
    /*
     * Get number of blocks
     */
    moris::uint
    get_num_blocks() const
    {
        return mListofBlocks.size();
    };

    // ----------------------------------------------------------------------------
    /*
     * Get block by index
     */
    moris::mtk::Set *
    get_block_by_index( moris::uint aBlockIndex) const
    {
        MORIS_ASSERT(aBlockIndex<mListofBlocks.size(),"Block index out of bounds");
        return mListofBlocks(aBlockIndex);
    };

    // ----------------------------------------------------------------------------
    /*
     * Get number of blocks
     * Sometimes num side set * 2. Ask Keenan
     */
    moris::uint
    get_num_side_set() const
    {
        return mListofSideSets.size();
    };

    // ----------------------------------------------------------------------------
    /*
     * Get block by index
     */
    moris::mtk::Set *
    get_side_set_by_index( moris::uint aSideSetIndex) const
    {
        MORIS_ASSERT(aSideSetIndex<mListofSideSets.size(),"Side set index out of bounds");
        return mListofSideSets(aSideSetIndex);
    };

    // ----------------------------------------------------------------------------
    /*
     * Get number of blocks
     * Sometimes num side set * 2. Ask Keenan
     */
    moris::uint
    get_num_double_side_set() const
    {
        return mListofDoubleSideSets.size();
    };

    // ----------------------------------------------------------------------------
    /*
     * Get block by index
     */
    moris::mtk::Set *
    get_double_side_set_by_index( moris::uint aSideSetIndex) const
    {
        MORIS_ASSERT(aSideSetIndex<mListofDoubleSideSets.size(),"Double side set index out of bounds");
        return mListofDoubleSideSets(aSideSetIndex);
    };

    // ----------------------------------------------------------------------------

    /*
     * Get cell clusters within a block set
     */
    virtual
    moris::Cell<Cluster const *>
    get_cell_clusters_in_set(moris_index aBlockSetOrdinal) const = 0;

    //##############################################
    // Side Set Cluster Access
    //##############################################
    /*!
     * Get side clusters within a side set
     */
    virtual
    moris::Cell<Cluster const *>
    get_side_set_cluster(moris_index aSideSetOrdinal) const = 0;

    // ----------------------------------------------------------------------------

    /*!
     * get number of side sets
     */
    virtual
    uint
    get_num_side_sets() const = 0;

    // ----------------------------------------------------------------------------

    /*!
     * Returns the label
     */
    virtual
    std::string
    get_side_set_label(moris_index aSideSetOrdinal) const = 0;

    // ----------------------------------------------------------------------------

    /*!
     * Returns the index given a label
     */
    virtual
    moris_index
    get_side_set_index(std::string aSideSetLabel) const  = 0;

    // ----------------------------------------------------------------------------

    //##############################################
    // Double Side Set Cluster Access
    //##############################################

    /*!
     * Returns the number of double sided side sets in the mesh
     */
    virtual
    uint
    get_num_double_sided_sets() const  = 0;

    // ----------------------------------------------------------------------------

    /*!
     * Returns the label
     */
    virtual
    std::string
    get_double_sided_set_label(moris_index aSideSetOrdinal) const = 0;

    // ----------------------------------------------------------------------------

    /*!
     * Returns the index given a label
     */
    virtual
    moris_index
    get_double_sided_set_index(std::string aDoubleSideSetLabel) const  = 0;

    // ----------------------------------------------------------------------------

    /*!
     * Returns the double side clusters in the side set
     */
    virtual
    moris::Cell<Cluster const*>
    get_double_side_set_cluster(moris_index aSideSetOrdinal) const  = 0;

    // ----------------------------------------------------------------------------


};
}
}


#endif /* PROJECTS_MTK_SRC_CL_MTK_INTEGRATION_MESH_HPP_ */
