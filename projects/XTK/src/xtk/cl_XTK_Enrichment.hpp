/*
 * cl_XTK_Enrichment.hpp
 *
 *  Created on: Feb 23, 2018
 *      Author: ktdoble
 */

#ifndef XTK_SRC_XTK_CL_XTK_ENRICHMENT_HPP_
#define XTK_SRC_XTK_CL_XTK_ENRICHMENT_HPP_

// XTKL: Linalg Includes
#include "linalg/cl_XTK_Matrix.hpp"
#include "linalg/cl_XTK_Matrix_Base_Utilities.hpp"


// Std includes
#include <limits>

// XTKL: XTK Includes
#include "cl_XTK_Child_Mesh.hpp"
#include "xtk/cl_XTK_Cut_Mesh.hpp"
#include "xtk/fn_mesh_flood_fill.hpp"
#include "xtk/fn_prune_element_to_element.hpp"
#include "xtk/fn_generate_element_to_element.hpp"
#include "xtk/fn_local_child_mesh_flood_fill.hpp"
#include "xtk/fn_generate_shared_face_element_graph.hpp"
#include "xtk/fn_assemble_boundary_subphase_constraint.hpp"
#include "xtk/fn_mesh_flood_fill.hpp"
#include "tools/fn_Pairing.hpp"


// Mesh includes
#include "cl_MTK_Mesh.hpp"
#include "cl_XTK_Background_Mesh.hpp"
#include "mesh/cl_Mesh_Enums.hpp"




/*
 * This class provides all the functions to perform the enrichment strategy on a child mesh
 */
namespace xtk
{

class Enrichment
{
public:
    Enrichment(moris::size_t aNumBulkPhases,
               Cut_Mesh* aCutMesh,
               Background_Mesh* aXTKMesh):
        mNumBulkPhases(aNumBulkPhases),
        mCutMesh(aCutMesh),
        mXTKMesh(aXTKMesh)
    {

    };

    bool mOutputFlag = false;
    moris::moris_index INDEX_MAX = std::numeric_limits<moris::moris_index>::max();


    /*
     * Performs basis function enrichment so that each element in connected regions of a given bulk phase are
     * assigned a unique enrichment level in the basis support.
     * @param[in] aCutMesh - Mesh containing elements around the interface
     * @param[in] aBackgroundMesh - Background mesh (Lagrangian Mesh)
     * @param[in] aMatrixFactory - Means of creating matrix objects
     *
     */
    void
    perform_enrichment()
    {
        // Start clock
        std::clock_t start = std::clock();

        // Perform local enrichment for each child mesh (commits local floodfill data to child mesh)
        perform_local_enrichment();

        // Perform enrichment over basis clusters
        perform_basis_cluster_enrichment();

        // Output time
        if(get_rank(get_comm()) == 0 && mOutputFlag)
        {
            std::cout<<"Enrichment completed in "<< (std::clock() - start) / (double)(CLOCKS_PER_SEC)<<" s."<<std::endl;
        }
    }


    /*
     * Returns the element ids in a basis support constructed in call to perform_enrichment. These are indexed by basis function index.
     */
    Cell<moris::Matrix< moris::IdMat >> const &
    get_element_ids_in_basis_support() const
    {
        return mElementIdsInBasis;
    }
    /*
    * Returns the element enrichment levels in a basis support constructed in call to perform_enrichment. These are indexed by basis function index.
    * Correspond to the element ids found at the same index in mElementIdsInBasis.
    */
    Cell<moris::Matrix< moris::IndexMat >> const &
    get_element_enrichment_levels_in_basis_support() const
    {
        return mElementEnrichmentLevel;
    }


private:
    moris::size_t mNumBulkPhases;

    // Pointers to Cut and XTK meshes (since they are used in most functions)
    Cut_Mesh* mCutMesh;
    Background_Mesh* mXTKMesh;

    // Enrichment Data ordered by basis function indices
    // For each basis function, the element ids and elemental subphases
    Cell<moris::Matrix< moris::IndexMat >> mElementEnrichmentLevel;
    Cell<moris::Matrix< moris::IdMat    >> mElementIdsInBasis;

    /*
     * performs local enrichment on all child meshes in the cut mesh. The subphase data (result of floodfill)
     * is stored as a member variable in each child mesh as sub-phase bins
     * @param[in] aCutMesh - Mesh containing elements around the interface
     * @param[in] aMatrixFactory - Means of creating matrix objects
     */
    void
    perform_local_enrichment()
    {

        // get the number of children meshes
        moris::size_t tNumChildMeshes = mCutMesh->get_num_simple_meshes();

        // iterate over children meshes and perform local flood-fill
        for(moris::size_t i = 0; i<tNumChildMeshes; i++)
        {
            // Get child mesh index
            Child_Mesh_Test & tChildMesh = mCutMesh->get_child_mesh(i);

            // Perform local flood-fill on child mesh to identify subphase
            moris::Matrix< moris::IndexMat > tLocalFloodFill = local_child_mesh_flood_fill(tChildMesh);

            // Set the local floodfill data as the elemental subphase values in the child mesh
            // The child mesh then sorts the elements into bins
            tChildMesh.set_elemental_subphase(tLocalFloodFill);
        }
    }

    /*
     * Performs enrichment on elements in support of full basis cluster. This enrichment includes all children elements of parents in
     * the basis cluster and parent elements with no children
     * @param[in] aCutMesh - Mesh containing elements around the interface
     * @param[in] aBackgroundMesh - Background mesh (Lagrangian Mesh)
     * @param[in] aMatrixFactory - Means of creating matrix objects
     */
    void
    perform_basis_cluster_enrichment()
    {
        // Get underlying matrix data to access function
        moris::mtk::Mesh & tXTKMeshData = mXTKMesh->get_mesh_data();

        // Number of basis functions
        moris::size_t tNumBasis          = tXTKMeshData.get_num_basis_functions();
        moris::size_t tNumFacePerElement = tXTKMeshData.get_entity_connected_to_entity_loc_inds(0, moris::EntityRank::ELEMENT, moris::EntityRank::FACE).n_cols();

        // Allocate member variables
        mElementEnrichmentLevel = Cell<moris::Matrix< moris::IndexMat >>(tNumBasis);
        mElementIdsInBasis      = Cell<moris::Matrix< moris::IndexMat >>(tNumBasis);

        for(moris::size_t i = 0; i<tNumBasis; i++)
        {

            // Initialize first available enrichment level for this basis cluster
            moris::size_t tFirstAvailableEnrich = 0;

            // Get elements in support of basis
            moris::Matrix< moris::IndexMat > tParentElementsInSupport = tXTKMeshData.get_elements_in_support_of_basis(i);

            // Cell 0 pruned element to element graph Cell 1 pruned shared face
            Cell<moris::Matrix< moris::IndexMat >> tPrunedData =
            generate_pruned_element_graph_in_basis_support( tNumFacePerElement,
                                                            tParentElementsInSupport);

            // Map from index in basis cluster to child mesh index and child mesh bin index also get the number of bins in basis
            // A cantor pairing of the child mesh index and child mesh bin index is used as the map key (hash value)
            std::unordered_map<moris::moris_index,moris::moris_index>  tSubphaseBinIndexToCMBinIndex;
            moris::Matrix< moris::IndexMat > tSubPhaseBinBulkPhase(1,1);
            moris::size_t tNumBinsInBasis = setup_all_subphase_bins_in_basis_support(tParentElementsInSupport,
                                                                               tSubphaseBinIndexToCMBinIndex,
                                                                               tSubPhaseBinBulkPhase);


            // Assign enrichment levels
            moris::Matrix< moris::IndexMat > tSubPhaseBinEnrichment =
            assign_subphase_bin_enrichment_levels_in_basis_support(tParentElementsInSupport, tNumBinsInBasis,
                                                                   tSubPhaseBinBulkPhase,
                                                                   tSubphaseBinIndexToCMBinIndex, tPrunedData(0),
                                                                   tPrunedData(1));

            // Extract element enrichment levels from assigned sub-phase bin enrichment levels and store these as a member variable
            unzip_subphase_bin_enrichment_into_element_enrichment(tParentElementsInSupport,
                                                                  tSubphaseBinIndexToCMBinIndex,
                                                                  tSubPhaseBinEnrichment,
                                                                  mElementIdsInBasis(i),
                                                                  mElementEnrichmentLevel(i));

        }

    }


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
    generate_pruned_element_graph_in_basis_support(moris::size_t const &                                         aNumFacePerElement,
                                                    moris::Matrix< moris::IndexMat > const &               aElementsInSupport)
     {
        moris::mtk::Mesh & tXTKMeshData = mXTKMesh->get_mesh_data();
        moris::size_t tMax = std::numeric_limits<moris::moris_index>::max();

        // Construct full element neighbor graph in support and the corresponding shared faces
        moris::Matrix< moris::IndexMat > tElementGraph(aElementsInSupport.n_cols(), aNumFacePerElement,tMax);
        moris::Matrix< moris::IndexMat > tSharedFaces(aElementsInSupport.n_cols(), aNumFacePerElement,tMax);

        for(moris::size_t iE = 0; iE<aElementsInSupport.n_cols(); iE++)
        {
            // Get elements connected to element and the corresponding face
            moris::Matrix< moris::IndexMat > tSingleElementToElement = tXTKMeshData.get_elements_connected_to_element_and_face_ind_loc_inds(aElementsInSupport(0,iE));
            replace_row(0,tSingleElementToElement,iE,tElementGraph,false);
            replace_row(1,tSingleElementToElement,iE,tSharedFaces,false);
        }


        // prune the graph to only include elements in the support
        // Note cell(0) is the element graph and cell(1) is the faces
        Cell<moris::Matrix< moris::IndexMat >> tPrunedNeighbors = prune_element_to_element(tElementGraph,aElementsInSupport,tSharedFaces,tMax);
        return tPrunedNeighbors;
     }



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
                                        moris::size_t const &                                             aNumSubPhaseBins,
                                        std::unordered_map<moris::moris_index,moris::moris_index> & aSubphaseBinIndexToCMBinIndex,
                                        moris::Matrix< moris::IndexMat > const &                    aPrunedElementGraph,
                                        moris::Matrix< moris::IndexMat > const &                    aPrunedSharedFaces)
    {
        moris::Matrix< moris::IndexMat > tSubphaseBinToSubphaseBin(aNumSubPhaseBins,aNumSubPhaseBins-1,INDEX_MAX);
        moris::Matrix< moris::DDSTMat > tSubphaseBinCounter(1,aNumSubPhaseBins,0);

        // Iterate through elements in support, constructing shared element graph using function generate shared_face_element_graph
        for(moris::size_t i = 0; i<aParentElementsInSupport.n_cols(); i++)
        {
            // Iterate through neighbors
            moris::moris_index tElementIndex0 = aParentElementsInSupport(0,i);

            for(moris::size_t j = 0; j<aPrunedElementGraph.n_cols(); j++)
            {
                // If there is an moris::size_t max at this entry, this all neighbors have been iterated over
                if( aPrunedElementGraph(i,j) == std::numeric_limits<moris::moris_index>::max())
                {
                    break;
                }

                // Second background  element on boundary
                moris::moris_index tElementIndex1 = aPrunedElementGraph(i,j);

                // The shared parent face index
                moris::moris_index tSharedFaceIndex = aPrunedSharedFaces(i,j);

                // Since both element are neighbors of each other, only create this relationship from the lowest to highest element
                // rather than repeat this twice
                if(tElementIndex0<tElementIndex1)
                {
                    bool tElement0HasChildren = mXTKMesh->entity_has_children(tElementIndex0, EntityRank::ELEMENT);
                    bool tElement1HasChildren = mXTKMesh->entity_has_children(tElementIndex1, EntityRank::ELEMENT);

                    // If both parent elements have children, then use the 2 child on the interface routine
                    if(tElement0HasChildren && tElement1HasChildren)
                    {
                        construct_subphase_bin_to_subphase_bin_2_child_interface(tElementIndex0,            tElementIndex1,
                                                                                 tSharedFaceIndex,          aSubphaseBinIndexToCMBinIndex,
                                                                                 tSubphaseBinToSubphaseBin, tSubphaseBinCounter);
                    }

                    // If element 0 does not have children but element 1 does have children
                    // use the mixed child parent on interface routine
                    else if (!tElement0HasChildren && tElement1HasChildren)
                    {
                        construct_subphase_bin_to_subphase_bin_mixed_interface(tElementIndex1,            tElementIndex0,
                                                                               tSharedFaceIndex,          aSubphaseBinIndexToCMBinIndex,
                                                                               tSubphaseBinToSubphaseBin, tSubphaseBinCounter);
                    }


                    // If element 1 does not have children but element 0 does have children
                    // use the mixed child parent on interface routine
                    else if (tElement0HasChildren && !tElement1HasChildren)
                    {
                        construct_subphase_bin_to_subphase_bin_mixed_interface(tElementIndex0,            tElementIndex1,
                                                                               tSharedFaceIndex,          aSubphaseBinIndexToCMBinIndex,
                                                                               tSubphaseBinToSubphaseBin, tSubphaseBinCounter);
                    }

                    // Both elements do not have children
                    // use the parent interface routing
                    else
                    {
                        construct_subphase_bin_to_subphase_bin_2_parent_interface(tElementIndex0,
                                                                                  tElementIndex1,
                                                                                  aSubphaseBinIndexToCMBinIndex,
                                                                                  tSubphaseBinToSubphaseBin,
                                                                                  tSubphaseBinCounter);
                    }
                }
            }


        }

        return tSubphaseBinToSubphaseBin;
    }


    void
    construct_subphase_bin_to_subphase_bin_2_child_interface(
            moris::moris_index const &                                  aParentElementIndex0,
            moris::moris_index const &                                  aParentElementIndex1,
            moris::moris_index const &                                  aSharedFaceIndex,
            std::unordered_map<moris::moris_index,moris::moris_index> & aSubphaseBinIndexToCMBinIndex,
            moris::Matrix< moris::IndexMat > &                          aSubphaseBinToSubphaseBin,
            moris::Matrix< moris::DDSTMat >  &                          aSubphaseBinCounter)
    {
        // Get the child mesh index
        moris::moris_index tChildMeshIndex0 = mXTKMesh->child_mesh_index(aParentElementIndex0,EntityRank::ELEMENT);
        moris::moris_index tChildMeshIndex1 = mXTKMesh->child_mesh_index(aParentElementIndex1,EntityRank::ELEMENT);

        // Get the child me shes on this boundary
        Child_Mesh_Test tChildMesh0 = mCutMesh->get_child_mesh(tChildMeshIndex0);
        Child_Mesh_Test tChildMesh1 = mCutMesh->get_child_mesh(tChildMeshIndex1);

        // Get child element subphase bin membership
        moris::Matrix< moris::IndexMat > const & tChildElements0BinMembership = tChildMesh0.get_elemental_subphase_bin_membership();
        moris::Matrix< moris::IndexMat > const & tChildElements1BinMembership = tChildMesh1.get_elemental_subphase_bin_membership();

        // Construct element pairs across shared parent face
        moris::Matrix< moris::IndexMat > tBoundaryElementPairs =
            generate_shared_face_element_pairs(aSharedFaceIndex,tChildMeshIndex0,tChildMeshIndex1,*mCutMesh);

        // iterate over pairs and create a relationship between their elements buckets
        for(moris::size_t k = 0; k<tBoundaryElementPairs.n_cols(); k++)
        {
            // Get pairs children element bin membership
            moris::moris_index tChildElemBin0 = tChildElements0BinMembership(0,(tBoundaryElementPairs)(0,k));
            moris::moris_index tChildElemBin1 = tChildElements1BinMembership(0,(tBoundaryElementPairs)(1,k));

            // get these bins basis cluster index
            moris::moris_index tChildElemBinBasisIndex0 = aSubphaseBinIndexToCMBinIndex[cantor_pairing(tChildMeshIndex0+1,tChildElemBin0)];
            moris::moris_index tChildElemBinBasisIndex1 = aSubphaseBinIndexToCMBinIndex[cantor_pairing(tChildMeshIndex1+1,tChildElemBin1)];

            //Check whether this bin relationship already exists
            bool tRelationshipExists = false;
            for(moris::size_t l = 0; l <= aSubphaseBinCounter(0,tChildElemBinBasisIndex0); l++)
            {
                if(aSubphaseBinToSubphaseBin(tChildElemBinBasisIndex0,l) == tChildElemBinBasisIndex1)
                {
                    tRelationshipExists = true;
                    break;
                }
            }

            // If this bin neighbor relationship does not exist, add it to each basis function
            if(tRelationshipExists == false)
            {
                aSubphaseBinToSubphaseBin(tChildElemBinBasisIndex0,aSubphaseBinCounter(0,tChildElemBinBasisIndex0)) = tChildElemBinBasisIndex1;
                aSubphaseBinCounter(0,tChildElemBinBasisIndex0)++;
                aSubphaseBinToSubphaseBin(tChildElemBinBasisIndex1,aSubphaseBinCounter(0,tChildElemBinBasisIndex1)) = tChildElemBinBasisIndex0;
                aSubphaseBinCounter(0,tChildElemBinBasisIndex1)++;
            }
        }
    }


    void
    construct_subphase_bin_to_subphase_bin_mixed_interface(
            moris::moris_index const &                                   aParentElementWithChildren,
            moris::moris_index const &                                   aParentElementWithoutChildren,
            moris::moris_index const &                                   aSharedFaceIndex,
            std::unordered_map<moris::moris_index,moris::moris_index> &  aSubphaseBinIndexToCMBinIndex,
            moris::Matrix< moris::IndexMat > &                           aSubphaseBinToSubphaseBin,
            moris::Matrix< moris::DDSTMat > &                            aSubphaseBinCounter)
    {
        // Get the child mesh index
        moris::moris_index tChildMeshIndex = mXTKMesh->child_mesh_index(aParentElementWithChildren,EntityRank::ELEMENT);

        // Get the child mesh on this shared face
        Child_Mesh_Test tChildMesh = mCutMesh->get_child_mesh(tChildMeshIndex);

        // Get child element subphase bin membership
        moris::Matrix< moris::IndexMat > const & tChildElementsBinMembership = tChildMesh.get_elemental_subphase_bin_membership();

        // Allocate Matrixes
        moris::Matrix< moris::IndexMat > tFaceOrdinals(1,1);
        moris::Matrix< moris::IndexMat > tChildrenElementCMInds(1,1);
        moris::Matrix< moris::IdMat > tChildrenElementIds(1,1);

        // Get children elements attached to aFaceIndex on the side of child mesh index 0
        mCutMesh->get_child_elements_connected_to_parent_face(tChildMeshIndex,
                                                             aSharedFaceIndex,
                                                             tChildrenElementIds,
                                                             tChildrenElementCMInds,
                                                             tFaceOrdinals);


        // iterate over child elements on boundary and construct the subphase bin relationship
        for(moris::size_t k = 0; k<tChildrenElementCMInds.n_cols(); k++)
        {
            // Get pairs children element bin membership
            moris::moris_index tChildElemBin = tChildElementsBinMembership(0,tChildrenElementCMInds(0,k));

            // Get the child  element bin basis index
            moris::moris_index tChildElemBinBasisIndex = aSubphaseBinIndexToCMBinIndex[cantor_pairing(tChildMeshIndex+1,tChildElemBin)];

            // Get the child  element bin basis index
            moris::moris_index tParentElementDummyIndex = 0;
            moris::moris_index tParentElemBinBasisIndex = aSubphaseBinIndexToCMBinIndex[cantor_pairing(tParentElementDummyIndex,aParentElementWithoutChildren)];

            //Check whether this bin relationship already exists
            bool tRelationshipExists = false;
            for(moris::size_t l = 0; l <= aSubphaseBinCounter(0,tChildElemBinBasisIndex); l++)
            {
                if(aSubphaseBinToSubphaseBin(tChildElemBinBasisIndex,l) == tParentElemBinBasisIndex)
                {
                    tRelationshipExists = true;
                    break;
                }
            }

            // If this bin neighbor relationship does not exist, add it to each basis function
            if(tRelationshipExists == false)
            {
                aSubphaseBinToSubphaseBin(tChildElemBinBasisIndex,aSubphaseBinCounter(0,tChildElemBinBasisIndex)) = tParentElemBinBasisIndex;
                aSubphaseBinCounter(0,tChildElemBinBasisIndex)++;
                aSubphaseBinToSubphaseBin(tParentElemBinBasisIndex,aSubphaseBinCounter(0,tParentElemBinBasisIndex)) = tChildElemBinBasisIndex;
                aSubphaseBinCounter(0,tParentElemBinBasisIndex)++;
           }
        }
    }


    void
    construct_subphase_bin_to_subphase_bin_2_parent_interface(
            moris::moris_index const &                                 aParentElementIndex0,
            moris::moris_index const &                                 aParentElementIndex1,
            std::unordered_map<moris::moris_index,moris::moris_index>  aSubphaseBinIndexToCMBinIndex,
            moris::Matrix< moris::IndexMat > &                         aSubphaseBinToSubphaseBin,
            moris::Matrix< moris::DDSTMat >  &                         aSubphaseBinCounter)
    {

            // get these bins basis cluster index
        moris::moris_index tParentElementDummyIndex = 0;
        moris::moris_index tParentElemBinBasisIndex0 = aSubphaseBinIndexToCMBinIndex[cantor_pairing(tParentElementDummyIndex,aParentElementIndex0)];
        moris::moris_index tParentElemBinBasisIndex1 = aSubphaseBinIndexToCMBinIndex[cantor_pairing(tParentElementDummyIndex,aParentElementIndex1)];

        //Check whether this bin relationship already exists
        bool tRelationshipExists = false;

        for(moris::size_t l = 0; l <= aSubphaseBinCounter(0,tParentElemBinBasisIndex0); l++)
        {
            if(aSubphaseBinToSubphaseBin(tParentElemBinBasisIndex0,l) == tParentElemBinBasisIndex1)
            {
                tRelationshipExists = true;
                break;
            }
        }

        // If this bin neighbor relationship does not exist, add it to each basis function
        if(tRelationshipExists == false)
        {
            aSubphaseBinToSubphaseBin(tParentElemBinBasisIndex0,aSubphaseBinCounter(0,tParentElemBinBasisIndex0)) = tParentElemBinBasisIndex1;
            aSubphaseBinCounter(0,tParentElemBinBasisIndex0)++;
            aSubphaseBinToSubphaseBin(tParentElemBinBasisIndex1,aSubphaseBinCounter(0,tParentElemBinBasisIndex1)) = tParentElemBinBasisIndex0;
            aSubphaseBinCounter(0,tParentElemBinBasisIndex1)++;
        }

    }

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
                                                           moris::Matrix< moris::IndexMat > const &                     aPrunedSharedFaces)
    {

        // Initialize Enrichment Vector now that the count includes all children
        moris::Matrix< moris::IndexMat > tSubPhaseBinEnrichmentVals(1, aNumSubPhaseBins);

        // Consider cutting down the number of arguements by moving aBackgroundMesh,aCutMesh,aMatrixFactory into reference member variables
        moris::Matrix< moris::IndexMat > tSubPhaseBinNeighborhood =
        construct_subphase_bin_neighborhood(aParentElementsInSupport,
                                            aNumSubPhaseBins,
                                            aSubphaseBinIndexToCMBinIndex,
                                            aPrunedElementGraph,
                                            aPrunedSharedFaces);

        // Variables needed for floodfill, consider removing these.
        // Active bins to include in floodfill (We include all bins)
        moris::Matrix< moris::IndexMat > tActiveBins(1,tSubPhaseBinNeighborhood.n_rows());
        for(moris::size_t i = 0; i< tSubPhaseBinNeighborhood.n_rows(); i++)
        {
            (tActiveBins)(0,i) = i;
        }

        // Mark all as included
        moris::Matrix< moris::IndexMat > tIncludedBins(1,tSubPhaseBinNeighborhood.n_rows(),1);

        tSubPhaseBinEnrichmentVals = flood_fill(tSubPhaseBinNeighborhood,
                                                aSubPhaseBinBulkPhase,
                                                tActiveBins,
                                                tIncludedBins,
                                                mNumBulkPhases,
                                                INDEX_MAX,
                                                true);


        return tSubPhaseBinEnrichmentVals;
    }

    moris::size_t
    setup_all_subphase_bins_in_basis_support(moris::Matrix< moris::IndexMat > const &               aParentElementsInSupport,
                                             std::unordered_map<moris::moris_index,moris::moris_index> & aSubphaseBinIndexToCMBinIndex,
                                             moris::Matrix< moris::IndexMat > &                  aSubPhaseBinBulkPhase)
    {
        // Count all elements including children of the elements in support of the basis
        moris::size_t tNumSubPhaseBinsInSupport = count_subphase_bins_in_support(aParentElementsInSupport);

        aSubPhaseBinBulkPhase.resize(1,tNumSubPhaseBinsInSupport);

        // Counter
        moris::size_t tCount = 0;

        // Child Mesh Index
        moris::size_t tChildMeshIndex = 0;

        // Number of children elements
        moris::size_t tNumChildrenElements = 0;

        // Bin Index
        moris::size_t tBinIndex = 0;

        for(moris::size_t i = 0; i<aParentElementsInSupport.n_cols(); i++)
        {
            // Add children elements to all elements in support vector but do not add parent
            if(mXTKMesh->entity_has_children(aParentElementsInSupport(0,i),EntityRank::ELEMENT))
            {
                tChildMeshIndex = mXTKMesh->child_mesh_index(aParentElementsInSupport(0,i),EntityRank::ELEMENT);

                // Get the child mesh
                Child_Mesh_Test & tChildMesh = mCutMesh->get_child_mesh(tChildMeshIndex);

                moris::size_t tNumLocSubPhaseBins = tChildMesh.get_num_subphase_bins();

                Cell<moris::moris_index> const & tSubPhaseBinBulkPhase = tChildMesh.get_subphase_bin_bulk_phase();

                for(moris::size_t j = 0; j<tNumLocSubPhaseBins; j++)
                {

                    // Create a hash for child mesh index and bin index
                    moris::size_t tCMIndexBinHash = cantor_pairing(tChildMeshIndex+1,j);

                    // Add to map
                    aSubphaseBinIndexToCMBinIndex[tCMIndexBinHash] = tBinIndex;

                    //
                    aSubPhaseBinBulkPhase(0,tBinIndex) = tSubPhaseBinBulkPhase(j);
                    tBinIndex++;

                }


            }

            // Parent element does not have children  case
            else
            {
                // Assign a child mesh index of zero to parent elements without children elements
                moris::size_t tChildMeshIndex = 0;

                // hash a child mesh index of 0 and parent element index
                moris::size_t tHash = cantor_pairing(tChildMeshIndex,(moris::size_t)aParentElementsInSupport(0,i));

                aSubphaseBinIndexToCMBinIndex[tHash] = tBinIndex;

                aSubPhaseBinBulkPhase(0,tBinIndex) = mXTKMesh->get_element_phase_index(aParentElementsInSupport(0,i));

                tBinIndex++;
            }
        }

        return tNumSubPhaseBinsInSupport;
    }


    void
    unzip_subphase_bin_enrichment_into_element_enrichment(moris::Matrix< moris::IndexMat > const &                     aParentElementsInSupport,
                                                          std::unordered_map<moris::moris_index, moris::moris_index> & aSubphaseBinIndexToCMBinIndex,
                                                          moris::Matrix< moris::IndexMat > const &                      aSubPhaseBinEnrichmentLevel,
                                                          moris::Matrix< moris::IndexMat > &       aElementIndInBasisSupport,
                                                          moris::Matrix< moris::IndexMat > &       aElementEnrichmentLevel)
    {

        // Count all elements including children of the elements in support of the basis
        moris::size_t tNumAllElementsInSupport = count_elements_in_support(aParentElementsInSupport);

        // Background mesh underlying meshd ata
        moris::mtk::Mesh const & tBackgroundMeshData = mXTKMesh->get_mesh_data();

        aElementIndInBasisSupport = moris::Matrix< moris::IndexMat >(1,tNumAllElementsInSupport);
        aElementEnrichmentLevel   = moris::Matrix< moris::IndexMat >(1,tNumAllElementsInSupport);

        // Counter
        moris::size_t tCount = 0;

        // Child Mesh Index
        moris::moris_index tChildMeshIndex = 0;

        // Number of children elements
        moris::size_t tNumChildrenElements = 0;

        for(moris::size_t i = 0; i<aParentElementsInSupport.n_cols(); i++)
        {
            // Add children elements to all elements in support vector but do not add parent
            if(mXTKMesh->entity_has_children(aParentElementsInSupport(0,i),EntityRank::ELEMENT))
            {
                moris::moris_index tChildMeshIndex = mXTKMesh->child_mesh_index(aParentElementsInSupport(0,i),EntityRank::ELEMENT);

                Child_Mesh_Test & tChildMesh = mCutMesh->get_child_mesh(tChildMeshIndex);

                moris::Matrix< moris::IndexMat > const & tChildElementSubphaseBin = tChildMesh.get_elemental_subphase_bin_membership();
                moris::Matrix< moris::IdMat > const & tChildElementIds = tChildMesh.get_element_ids();
                tNumChildrenElements = tChildElementIds.n_cols();

                for(moris::size_t j = 0; j<tNumChildrenElements; j++)
                {
                    aElementIndInBasisSupport(0,tCount) = tChildElementIds(0,j);

                    moris::moris_index tChildElementBinHash = cantor_pairing(tChildMeshIndex+1,tChildElementSubphaseBin(0,j));

                    moris::moris_index tBinBasisIndex = aSubphaseBinIndexToCMBinIndex[tChildElementBinHash];

                    aElementEnrichmentLevel(0,tCount) = aSubPhaseBinEnrichmentLevel(0,tBinBasisIndex);

                    tCount++;
                }

            }

            else
            {
                moris::moris_index tParentElementDummyIndex = 0;
                moris::moris_index tParentHash = cantor_pairing(tParentElementDummyIndex,aParentElementsInSupport(0,i));
                moris::moris_index tBinBasisIndex = aSubphaseBinIndexToCMBinIndex[tParentHash];


                aElementIndInBasisSupport(0,tCount) = tBackgroundMeshData.get_glb_entity_id_from_entity_loc_index(aParentElementsInSupport(0,i),moris::EntityRank::ELEMENT);
                aElementEnrichmentLevel(0,tCount) = aSubPhaseBinEnrichmentLevel(0,tBinBasisIndex);

                tCount++;
            }

        }
    }



    moris::size_t
    count_subphase_bins_in_support(moris::Matrix< moris::IndexMat > const &               aParentElementsInSupport)
    {
        // Number of elements in this support (need both parent and total)
        moris::size_t tNumParentElementsInSupport = aParentElementsInSupport.n_cols();

        // Initialize number of sub phase bins
        moris::size_t tNumSubPhaseBins  = 0;

        // initialize variable for child mesh index if an element has children
        moris::size_t tChildMeshIndex = 0;

        // initialize variable for number of children elements in a mesh
        moris::size_t tNumChildElements = 0;

        // Count children elements in support
        for(moris::size_t i = 0; i<tNumParentElementsInSupport; i++)
        {
            // Check if this element has children and if it does add them to the count
            if(mXTKMesh->entity_has_children(aParentElementsInSupport(0,i),EntityRank::ELEMENT))
            {
                // The child mesh index
                tChildMeshIndex = mXTKMesh->child_mesh_index(aParentElementsInSupport(0,i),EntityRank::ELEMENT);

                // Get the child mesh
                Child_Mesh_Test & tChildMesh = mCutMesh->get_child_mesh(tChildMeshIndex);

                // Get the number of subphase bins
                tNumSubPhaseBins = tNumSubPhaseBins + tChildMesh.get_num_subphase_bins();
            }

            else
            {
                tNumSubPhaseBins++;
            }
        }

        return tNumSubPhaseBins;
    }

    /*
     * Given parent elements from the background mesh, counts all children of these parent elements.
     * Used for sizing of matrices for next step
     * @param[in] aParentElementsInSupport - Parent element index in support of the basis function
     * @param[in] aCutMesh                 - Mesh containing elements around the interface
     * @param[in] aBackgroundMesh          - Background mesh (Lagrangian Mesh)
     */
    moris::size_t
    count_elements_in_support(moris::Matrix< moris::IndexMat > const &            aParentElementsInSupport)
    {

        // Number of elements in this support (need both parent and total)
        moris::size_t tNumParentElementsInSupport = aParentElementsInSupport.n_cols();
        moris::size_t tNumElementsInSupport = tNumParentElementsInSupport;

        // initialize variable for child mesh index if an element has children
        moris::size_t tChildMeshIndex = 0;

        // initialize variable for number of children elements in a mesh
        moris::size_t tNumChildElements = 0;

        // Count children elements in support
        for(moris::size_t i = 0; i<tNumParentElementsInSupport; i++)
        {
            // Check if this element has children and if it does add them to the count
            if(mXTKMesh->entity_has_children(aParentElementsInSupport(0,i),EntityRank::ELEMENT))
            {
                // The child mesh index
                tChildMeshIndex = mXTKMesh->child_mesh_index(aParentElementsInSupport(0,i),EntityRank::ELEMENT);

                // Number of child elements in this mesh
                tNumChildElements = mCutMesh->get_num_entities(tChildMeshIndex,EntityRank::ELEMENT);

                // Add the number of elements in the child mesh to the account (-1 to remove parent)
                tNumElementsInSupport = tNumElementsInSupport + tNumChildElements - 1;
            }
        }

        return tNumElementsInSupport;
    }





};
}
#endif /* XTK_SRC_XTK_CL_XTK_ENRICHMENT_HPP_ */
