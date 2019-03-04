/*
 * cl_XTK_Enrichment.hpp
 *
 *  Created on: Feb 23, 2018
 *      Author: ktdoble
 */

#ifndef XTK_SRC_XTK_CL_XTK_ENRICHMENT_HPP_
#define XTK_SRC_XTK_CL_XTK_ENRICHMENT_HPP_

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
#include "fn_mesh_flood_fill.hpp"
#include "fn_prune_element_to_element.hpp"
#include "fn_generate_element_to_element.hpp"
#include "fn_local_child_mesh_flood_fill.hpp"
#include "fn_generate_shared_face_element_graph.hpp"
#include "fn_assemble_boundary_subphase_constraint.hpp"
#include "fn_mesh_flood_fill.hpp"
#include "fn_Pairing.hpp"
#include "fn_equal_to.hpp"


// Mesh includes
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_Mesh_Enums.hpp"
#include "cl_XTK_Background_Mesh.hpp"
#include "cl_Mesh_Enums.hpp"




/*
 * This class provides all the functions to perform the enrichment strategy on a child mesh
 */
namespace xtk
{

class Enrichment_Parameters
{
public:
    Enrichment_Parameters(){};

    enum moris::EntityRank mBasisToEnrich = EntityRank::NODE ; /*For lagrange mesh this is node, for HMR this may be bsplines*/
};

/*
 * Enrichment of a vertex
 */
class Vertex_Enrichment
{
public:

    Vertex_Enrichment():
     mNodeIndex(MORIS_INDEX_MAX)
    {}

    void
    set_node_index(moris::moris_index aNodeIndex)
    {
        MORIS_ASSERT(mNodeIndex==MORIS_INDEX_MAX,"Node index already set for Vertex Enrichment");
        mNodeIndex = aNodeIndex;
    }

    /*
     * Add the basis information which includes weights, enrichment level, and basis index.
     * There is no "smartness" in this function. Duplicates should have been removed prior to call
     * An assertion will catch duplicates in debug mode
     */
    void
    add_basis_information( moris::Cell<moris::moris_index> aBasisIndices,
                           moris::Cell<moris::moris_index> aBasisEnrichmentLevel )
    {
        MORIS_ASSERT((aBasisIndices.size() == aBasisEnrichmentLevel.size()),"dimension mismatch in add basis information call");
#ifdef DEBUG
        // since I can't write these functions in one line, need to have ifdef
        uint tSize = aBasisIndices.size();
        moris::unique(aBasisIndices);
        uint tUniqueSize = aBasisIndices.size();

        MORIS_ASSERT(tSize == tUniqueSize, "duplicate basis indices detected" );
#endif


        // num basis
        moris::uint tNumBasis = aBasisIndices.size();

        // allocate space
        mBasisIndices.resize(tNumBasis,1);
        mBasisEnrichmentLevel.resize(tNumBasis,1);
        mBasisWeights.resize(tNumBasis,1);

        // iterate to store data
        for(moris::uint i = 0; i<aBasisIndices.size(); i++ )
        {
            moris::uint tBasisLocInd = this->local_basis_index(aBasisIndices(i));
            mBasisIndices(tBasisLocInd)         = aBasisIndices(i);
            mBasisEnrichmentLevel(tBasisLocInd) = aBasisEnrichmentLevel(i);
        }
    }

    void
    add_basis_weights(moris::Matrix<moris::IndexMat> const & aBasisIndices,
                      moris::Matrix<moris::DDRMat>   const & aBasisWeight)
    {
        for(moris::uint i = 0; i <aBasisIndices.numel(); i++)
        {
            moris::uint tBasisLocInd = this->local_basis_index(aBasisIndices(i));
            mBasisWeights(tBasisLocInd) = aBasisWeight(i);
        }
    }


    std::unordered_map<moris::moris_index, moris::moris_index> &
    get_basis_map()
    {
        return mBasisMap;
    }

    moris::uint
    local_basis_index(moris::uint aBasisIndex)
    {
        auto tIter = mBasisMap.find(aBasisIndex);
        MORIS_ASSERT(tIter!=mBasisMap.end(),"Provided basis index not found in map");

        return tIter->second;
    }

    void
    condense_out_basis_with_0_weight()
    {
        moris::uint tCount = 0;

        for( moris::uint i = 0;  i<mBasisIndices.numel(); i++)
        {
            if(moris::equal_to( mBasisWeights(i), 0))
            {
                mBasisMap.erase(mBasisIndices(i));
            }

            else
            {
                mBasisIndices(tCount) = mBasisIndices(i);
                mBasisWeights(tCount) = mBasisWeights(i);
                mBasisEnrichmentLevel(tCount) = mBasisEnrichmentLevel(i);

                // change map index
                mBasisMap[mBasisIndices(i)] = tCount;
                tCount++;
            }
        }

        // remove excess space
        mBasisIndices.resize(tCount,1);
        mBasisWeights.resize(tCount,1);
        mBasisEnrichmentLevel.resize(tCount,1);

    }

    moris::Matrix< moris::IndexMat > const &
    get_basis_basis_indices() const
    {
        return mBasisIndices;
    }

    moris::Matrix< moris::IndexMat > const &
    get_basis_basis_enrichment_level() const
    {
        return mBasisEnrichmentLevel;
    }

    moris::Matrix< moris::DDRMat > const &
    get_basis_weights() const
    {
        return mBasisWeights;
    }

    moris::Matrix< moris::DDRMat > &
    get_basis_weights()
    {
        return mBasisWeights;
    }


private:
    moris::moris_index               mNodeIndex;
    moris::Matrix< moris::IndexMat > mBasisIndices;
    moris::Matrix< moris::IndexMat > mBasisEnrichmentLevel;
    moris::Matrix< moris::DDRMat >   mBasisWeights;
    std::unordered_map<moris::moris_index, moris::moris_index> mBasisMap; /*From basis to local index*/

};



class Model;

class Enrichment
{
public:
    Enrichment(){};

    Enrichment(moris::size_t          aNumBulkPhases,
               xtk::Model*            aXTKModelPtr,
               xtk::Cut_Mesh*         aCutMeshPtr,
               xtk::Background_Mesh*  aBackgroundMeshPtr);

    bool mVerbose = false;
    moris::moris_index INDEX_MAX = std::numeric_limits<moris::moris_index>::max();


    /*!
     * Performs basis function enrichment so that each element in connected regions of a given bulk phase are
     * assigned a unique enrichment level in the basis support.
     * @param[in] aCutMesh - Mesh containing elements around the interface
     * @param[in] aBackgroundMesh - Background mesh (Lagrangian Mesh)
     * @param[in] aMatrixFactory - Means of creating matrix objects
     *
     */
    void
    perform_enrichment();

    /*!
     * Returns the element inds in a basis support constructed in call to perform_enrichment. These are indexed by basis function index.
     */
    Cell<moris::Matrix< moris::IdMat >> const &
    get_element_inds_in_basis_support() const;

    /*!
    * Returns the element enrichment levels in a basis support constructed in call to perform_enrichment. These are indexed by basis function index.
    * Correspond to the element inds found at the same index in mElementIndsInBasis.
    */
    Cell<moris::Matrix< moris::IndexMat >> const &
    get_element_enrichment_levels_in_basis_support() const;


    // ----------------------------------------------------------------------------------
    // Accessing enrichment data
    // ----------------------------------------------------------------------------------

    /*!
     * Returns the vertex enrichment for the provided vertex index
     */
    Vertex_Enrichment const &
    get_vertex_enrichment(moris::moris_index aVertexIndex) const
    {
        MORIS_ASSERT(aVertexIndex<(moris::moris_index)mVertexEnrichments.size(), "Specifed vertex index is out of bounds");
        return mVertexEnrichments(aVertexIndex);
    }

    /*!
     * Get element to basis connectivity
     */
    moris::Cell<moris::moris_index> const &
    get_element_to_basis_connectivity(moris::moris_index aElementIndex) const
    {
        MORIS_ASSERT(aElementIndex<(moris::moris_index)mElementToBasis.size(),"Element index out of bounds");
        return mElementToBasis(aElementIndex);
    }

    /*!
     * Get element to basis enrichment level connectivity
     */
    moris::Cell<moris::moris_index> const &
    get_element_to_basis_enrichment_level(moris::moris_index aElementIndex) const
    {
        MORIS_ASSERT(aElementIndex<(moris::moris_index)mElementToBasis.size(),"Element index out of bounds");
        return mElementToBasisEnrichmentLevel(aElementIndex);
    }

    /*!
     * Get the enriched basis index
     */
    moris::moris_index
    get_enriched_basis_index(moris::moris_index aElementIndex,
                             moris::moris_index aBasisIndex) const
    {

        moris::Cell<moris::moris_index> const & tElementToBasis         = this->get_element_to_basis_connectivity(aElementIndex);
        moris::Cell<moris::moris_index> const & tElementToBasisEnrLevel = this->get_element_to_basis_enrichment_level(aElementIndex);
        moris::print(mElementIndsInBasis(aBasisIndex),"basis elements");
        std::cout<<"aElementIndex = "<<aElementIndex<<std::endl;
        std::cout<<"aBasisIndex = "<<aBasisIndex<<std::endl;
        moris::print(tElementToBasis,"tElementToBasis");


        for(moris::uint i = 0; i <tElementToBasis.size(); i++)
        {
            std::cout<< " i = "<<i<<" | "<<tElementToBasis(i)<<"  |  "<<aBasisIndex<<std::endl;

            if( tElementToBasis(i) == aBasisIndex)
            {
                return mBasisEnrichmentIndices(aBasisIndex)(tElementToBasisEnrLevel(i));
            }
        }

        MORIS_ASSERT(0,"Basis index not found");
        return MORIS_INDEX_MAX;
    }

private:
    moris::size_t mNumBulkPhases;

    // Pointers to Model, Cut and Background meshes (since they are used in most functions)
    Model*    mXTKModelPtr              = nullptr;
    Cut_Mesh* mCutMeshPtr               = nullptr;
    Background_Mesh* mBackgroundMeshPtr = nullptr;
    Enrichment_Parameters mParameters;

    // Enrichment Data ordered by basis function indices
    // For each basis function, the element indices and elemental subphases
    Cell<moris::Matrix< moris::IndexMat >> mElementEnrichmentLevel;
    Cell<moris::Matrix< moris::IndexMat >> mElementIndsInBasis;
    Cell<xtk::Vertex_Enrichment>           mVertexEnrichments;

    // element to basis and enrichment level connectivity
    // for a given element, the basis function and enrichment level of that basis function
    // (transpose of mElementIndsInBasis)
    //TODO: Only store parent element to basis and element index to enrichment level
    moris::Cell<moris::Cell<moris::moris_index>> mElementToBasis;
    moris::Cell<moris::Cell<moris::moris_index>> mElementToBasisEnrichmentLevel;

    // Basis enrichment level indics
    moris::Cell<moris::Matrix<moris::IndexMat>> mBasisEnrichmentIndices;

    // total number of basis enrichment levels (all basis functions)
    moris::uint mNumEnrichmentLevels;

    /*
     * performs local enrichment on all child meshes in the cut mesh. The subphase data (result of floodfill)
     * is stored as a member variable in each child mesh as sub-phase bins
     * @param[in] aCutMesh - Mesh containing elements around the interface
     * @param[in] aMatrixFactory - Means of creating matrix objects
     */
    void
    perform_local_subphase_identification();

    /*
     * Performs enrichment on elements in support of full basis cluster. This enrichment includes all children elements of parents in
     * the basis cluster and parent elements with no children
     * @param[in] aCutMesh - Mesh containing elements around the interface
     * @param[in] aBackgroundMesh - Background mesh (Lagrangian Mesh)
     * @param[in] aMatrixFactory - Means of creating matrix objects
     */
    void
    perform_basis_cluster_enrichment();


    /*
     * Starting from the full element to element graph of elements in the support of basis aBasisIndex,  removes elements from
     *  this graph that are not in the support, creating a pruned element to element graph
     * @param[in] aNumFacePerElement - Number of faces per given element
     * @param[in] aElementsInSupport - Elements in support of basis
     * @param[in] aCutMesh - Mesh containing elements around the interface
     * @param[in] aBackgroundMesh - Background mesh (Lagrangian Mesh)
     * @param[in] aMatrixFactory - Means of creating matrix objects
     * @param[out] Cell(0) Pruned element to element graph only including the elements found in the aElementsInPrunedGraph vector
     *             Cell(1) Pruned shared faces
     */
    Cell<moris::Matrix< moris::IndexMat >>
    generate_pruned_element_graph_in_basis_support(moris::size_t const &                     aNumFacePerElement,
                                                    moris::Matrix< moris::IndexMat > const & aElementsInSupport);


    /*
     * Constructs the subphase bin neighboorhood graph within the support of a basis function (think of element to element graph but for sub-phase bins).
     * This is accomplished by tieing subphase bins across parent element borders.
     * @param[in] aParentElementsInSupport      - Parent elements in support of basis
     * @param[in] aNumSubPhaseBins              - Number of subphase bins in support of basis
     * @param[in] aSubphaseBinIndexToCMBinIndex - Map from child mesh index and child mesh bin index to basis bin index
     * @param[in] aPrunedElementGraph           - Pruned parent element graph only including parent elements in basis support
     * @param[in] aPrunedSharedFaces            - Shared faces corresponding to aPrunedElementGraph
     * @param[in] aCutMesh                      - Mesh containing elements around the interface
     * @param[in] aBackgroundMesh               - Background mesh (Lagrangian Mesh)
     * @param[in] aMatrixFactory                - Means of creating matrix objects
     */
    moris::Matrix< moris::IndexMat >
    construct_subphase_bin_neighborhood(moris::Matrix< moris::IndexMat > const &                    aParentElementsInSupport,
                                        moris::size_t                    const &                    aNumSubPhaseBins,
                                        std::unordered_map<moris::moris_index,moris::moris_index> & aSubphaseBinIndexToCMBinIndex,
                                        moris::Matrix< moris::IndexMat > const &                    aPrunedElementGraph,
                                        moris::Matrix< moris::IndexMat > const &                    aPrunedSharedFaces);


    void
    construct_subphase_bin_to_subphase_bin_2_child_interface(
            moris::moris_index const &                                  aParentElementIndex0,
            moris::moris_index const &                                  aParentElementIndex1,
            moris::moris_index const &                                  aSharedFaceIndex,
            std::unordered_map<moris::moris_index,moris::moris_index> & aSubphaseBinIndexToCMBinIndex,
            moris::Matrix< moris::IndexMat > &                          aSubphaseBinToSubphaseBin,
            moris::Matrix< moris::DDSTMat >  &                          aSubphaseBinCounter);


    void
    construct_subphase_bin_to_subphase_bin_mixed_interface(
            moris::moris_index const &                                   aParentElementWithChildren,
            moris::moris_index const &                                   aParentElementWithoutChildren,
            moris::moris_index const &                                   aSharedFaceIndex,
            std::unordered_map<moris::moris_index,moris::moris_index> &  aSubphaseBinIndexToCMBinIndex,
            moris::Matrix< moris::IndexMat > &                           aSubphaseBinToSubphaseBin,
            moris::Matrix< moris::DDSTMat > &                            aSubphaseBinCounter);


    void
    construct_subphase_bin_to_subphase_bin_2_parent_interface(
            moris::moris_index const &                                 aParentElementIndex0,
            moris::moris_index const &                                 aParentElementIndex1,
            std::unordered_map<moris::moris_index,moris::moris_index>  aSubphaseBinIndexToCMBinIndex,
            moris::Matrix< moris::IndexMat > &                         aSubphaseBinToSubphaseBin,
            moris::Matrix< moris::DDSTMat >  &                         aSubphaseBinCounter);

    /*
     *  Assigns enrichment levels to each of the sub-phase bins in support of a basis function by considering local sub-phase bins and
     *  cross parent element boundary bins tieing.
     * @param[in] aParentElementsInSupport      - Parent elements in support of basis
     * @param[in] aNumSubPhaseBins              - Number of subphase bins in support of basis
     * @param[in] aSubPhaseBinBulkPhase         - Bulk phase index of each of the subphase bins
     * @param[in] aSubphaseBinIndexToCMBinIndex - Map from child mesh index and child mesh bin index to basis bin index
     * @param[in] aPrunedElementGraph           - Pruned parent element graph only including parent elements in basis support
     * @param[in] aPrunedSharedFaces            - Shared faces corresponding to aPrunedElementGraph
     * @param[in] aCutMesh                      - Mesh containing elements around the interface
     * @param[in] aBackgroundMesh               - Background mesh (Lagrangian Mesh)
     * @param[in] aMatrixFactory                - Means of creating matrix objects
     */
    moris::Matrix< moris::IndexMat >
    assign_subphase_bin_enrichment_levels_in_basis_support(moris::Matrix< moris::IndexMat > const &                     aParentElementsInSupport,
                                                           moris::size_t const &                                              aNumSubPhaseBins,
                                                           moris::Matrix< moris::IndexMat > const &                     aSubPhaseBinBulkPhase,
                                                           std::unordered_map<moris::moris_index, moris::moris_index> & aSubphaseBinIndexToCMBinIndex,
                                                           moris::Matrix< moris::IndexMat > const &                     aPrunedElementGraph,
                                                           moris::Matrix< moris::IndexMat > const &                     aPrunedSharedFaces);

    moris::size_t
    setup_all_subphase_bins_in_basis_support(moris::Matrix< moris::IndexMat > const &               aParentElementsInSupport,
                                             std::unordered_map<moris::moris_index,moris::moris_index> & aSubphaseBinIndexToCMBinIndex,
                                             moris::Matrix< moris::IndexMat > &                  aSubPhaseBinBulkPhase);


    void
    unzip_subphase_bin_enrichment_into_element_enrichment(moris::Matrix< moris::IndexMat > const &                     aParentElementsInSupport,
                                                          std::unordered_map<moris::moris_index, moris::moris_index> & aSubphaseBinIndexToCMBinIndex,
                                                          moris::Matrix< moris::IndexMat > const &                      aSubPhaseBinEnrichmentLevel,
                                                          moris::Matrix< moris::IndexMat > &       aElementIndInBasisSupport,
                                                          moris::Matrix< moris::IndexMat > &       aElementEnrichmentLevel);


    void
    setup_vertex_enrichment_data();


    moris::size_t
    count_subphase_bins_in_support(moris::Matrix< moris::IndexMat > const & aParentElementsInSupport);

    /*
     * Given parent elements from the background mesh, counts all children of these parent elements.
     * Used for sizing of matrices for next step
     * @param[in] aParentElementsInSupport - Parent element index in support of the basis function
     * @param[in] aCutMesh                 - Mesh containing elements around the interface
     * @param[in] aBackgroundMesh          - Background mesh (Lagrangian Mesh)
     */
    moris::size_t
    count_elements_in_support(moris::Matrix< moris::IndexMat > const &            aParentElementsInSupport);


    /*!
     * Construct element to basis connectivity and corresponding enrichment level
     */
    void
    construct_element_to_basis_connectivity(moris::Cell<moris::Cell<moris::moris_index>> & aElementToBasis,
                                            moris::Cell<moris::Cell<moris::moris_index>> & aElementToBasisEnrichmentLevel);

    void
    construct_vertex_enrichment_with_element_to_basis(moris::Cell<moris::Cell<moris::moris_index>> const & aElementToBasis,
                                                      moris::Cell<moris::Cell<moris::moris_index>> const & aElementToBasisEnrichmentLevel);


    /*
     * This function sets up the mVertexEnrichments and element to basis and enrichment level connectivity
     */
    void
    compute_vertex_basis_weights();

    /*
     * Combines tmatrices between lagrange nodes and their respective basis counterpart
     */
    void
    combine_t_matrices();

    /*
     * Assigns each enrichment level at the basis function a local proc index and globally unique id
     */
    void
    assign_enrichment_level_identifiers();





};
}
#endif /* XTK_SRC_XTK_CL_XTK_ENRICHMENT_HPP_ */
