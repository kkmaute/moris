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
#include "xtk/cl_XTK_Mesh.hpp"
#include "xtk/cl_XTK_Face_Registry.hpp"
#include "xtk/fn_mesh_flood_fill.hpp"
#include "xtk/fn_prune_element_to_element.hpp"
#include "xtk/fn_generate_element_to_element.hpp"
#include "xtk/fn_local_child_mesh_flood_fill.hpp"
#include "xtk/fn_generate_shared_face_element_graph.hpp"
#include "xtk/fn_assemble_boundary_subphase_constraint.hpp"
#include "xtk/fn_mesh_flood_fill.hpp"
#include "tools/fn_Pairing.hpp"


// Mesh includes
#include "mesh/cl_Mesh_Data.hpp"
#include "mesh/cl_Mesh_Enums.hpp"




/*
 * This class provides all the functions to perform the enrichment strategy on a child mesh
 */
namespace xtk
{

template<typename Real, typename Integer, typename Real_Matrix, typename Integer_Matrix>
class Enrichment
{
public:
    Enrichment(Integer aNumBulkPhases):
        mNumBulkPhases(aNumBulkPhases)
    {

    };

    Integer INTEGER_MAX = std::numeric_limits<Integer>::max();


    /*
     * Performs basis function enrichment so that each element in connected regions of a given bulk phase are
     * assigned a unique enrichment level in the basis support.
     * @param[in] aCutMesh - Mesh containing elements around the interface
     * @param[in] aBackgroundMesh - Background mesh (Lagrangian Mesh)
     * @param[in] aMatrixFactory - Means of creating matrix objects
     *
     */
    void
    perform_enrichment(Cut_Mesh<Real, Integer, Real_Matrix, Integer_Matrix> &       aCutMesh,
                       XTK_Mesh<Real, Integer, Real_Matrix, Integer_Matrix> &       aBackgroundMesh)
    {
        // Start clock
        std::clock_t start = std::clock();

        // Perform local enrichment for each child mesh (commits local floodfill data to child mesh)
        perform_local_enrichment(aCutMesh);

        // Perform enrichment over basis clusters
        perform_basis_cluster_enrichment(aCutMesh, aBackgroundMesh);

        // Output time
        if(get_rank(get_comm()) == 0)
        {
            std::cout<<"Enrichment completed in "<< (std::clock() - start) / (double)(CLOCKS_PER_SEC)<<" s."<<std::endl;
        }
    }


    /*
     * Returns the element ids in a basis support constructed in call to perform_enrichment. These are indexed by basis function index.
     */
    Cell<moris::Mat_New<Integer, Integer_Matrix>> const &
    get_element_ids_in_basis_support() const
    {
        return mElementIdsInBasis;
    }
    /*
    * Returns the element enrichment levels in a basis support constructed in call to perform_enrichment. These are indexed by basis function index.
    * Correspond to the element ids found at the same index in mElementIdsInBasis.
    */
    Cell<moris::Mat_New<Integer, Integer_Matrix>> const &
    get_element_enrichment_levels_in_basis_support() const
    {
        return mElementEnrichmentLevel;
    }


private:
    Integer mNumBulkPhases;

    // Enrichment Data ordered by basis function indices
    // For each basis function, the element ids and elemental subphases
    Cell<moris::Mat_New<Integer, Integer_Matrix>> mElementEnrichmentLevel;
    Cell<moris::Mat_New<Integer, Integer_Matrix>> mElementIdsInBasis;

    /*
     * performs local enrichment on all child meshes in the cut mesh. The subphase data (result of floodfill)
     * is stored as a member variable in each child mesh as sub-phase bins
     * @param[in] aCutMesh - Mesh containing elements around the interface
     * @param[in] aMatrixFactory - Means of creating matrix objects
     */
    void
    perform_local_enrichment(Cut_Mesh<Real, Integer, Real_Matrix, Integer_Matrix> &       aCutMesh)
    {

        // get the number of children meshes
        Integer tNumChildMeshes = aCutMesh.get_num_simple_meshes();

        // iterate over children meshes and perform local flood-fill
        for(Integer i = 0; i<tNumChildMeshes; i++)
        {
            // Get child mesh index
            Child_Mesh_Test<Real, Integer, Real_Matrix, Integer_Matrix> & tChildMesh = aCutMesh.get_child_mesh(i);

            // Perform local flood-fill on child mesh to identify subphase
            moris::Mat_New<Integer, Integer_Matrix> tLocalFloodFill = local_child_mesh_flood_fill(tChildMesh);

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
    perform_basis_cluster_enrichment(Cut_Mesh<Real, Integer, Real_Matrix, Integer_Matrix> & aCutMesh,
                                     XTK_Mesh<Real, Integer, Real_Matrix, Integer_Matrix> & aBackgroundMesh)
    {
        // Get underlying matrix data to access function
        mesh::Mesh_Data<Real, Integer, Real_Matrix, Integer_Matrix> & tXTKMeshData = aBackgroundMesh.get_mesh_data();

        // Number of basis functions
        Integer tNumBasis          = tXTKMeshData.get_num_basis_functions();
        Integer tNumFacePerElement = tXTKMeshData.get_entity_connected_to_entity_loc_inds(0, EntityRank::ELEMENT, EntityRank::FACE).n_cols();

        // Allocate member variables
        mElementEnrichmentLevel = Cell<moris::Mat_New<Integer, Integer_Matrix>>(tNumBasis);
        mElementIdsInBasis      = Cell<moris::Mat_New<Integer, Integer_Matrix>>(tNumBasis);

        for(Integer i = 0; i<tNumBasis; i++)
        {

            // Initialize first available enrichment level for this basis cluster
            Integer tFirstAvailableEnrich = 0;

            // Get elements in support of basis
            moris::Mat_New<Integer, Integer_Matrix> tParentElementsInSupport = tXTKMeshData.get_elements_in_basis_support(i);

            // Cell 0 pruned element to element graph Cell 1 pruned shared face
            Cell<moris::Mat_New<Integer, Integer_Matrix>> tPrunedData =
            generate_pruned_element_graph_in_basis_support( tNumFacePerElement,
                                                            tParentElementsInSupport,
                                                            aCutMesh,
                                                            aBackgroundMesh);

            // Map from index in basis cluster to child mesh index and child mesh bin index also get the number of bins in basis
            // A cantor pairing of the child mesh index and child mesh bin index is used as the map key (hash value)
            std::unordered_map<Integer,Integer>  tSubphaseBinIndexToCMBinIndex;
            moris::Mat_New<Integer, Integer_Matrix> tSubPhaseBinBulkPhase(1,1);
            Integer tNumBinsInBasis = setup_all_subphase_bins_in_basis_support(tParentElementsInSupport,
                                                                               aCutMesh,
                                                                               aBackgroundMesh,
                                                                               tSubphaseBinIndexToCMBinIndex,
                                                                               tSubPhaseBinBulkPhase);


            // Assign enrichment levels
            moris::Mat_New<Integer, Integer_Matrix> tSubPhaseBinEnrichment =
            assign_subphase_bin_enrichment_levels_in_basis_support(tParentElementsInSupport, tNumBinsInBasis,
                                                                   tSubPhaseBinBulkPhase,
                                                                   tSubphaseBinIndexToCMBinIndex, tPrunedData(0),
                                                                   tPrunedData(1), aCutMesh,
                                                                   aBackgroundMesh);

            // Extract element enrichment levels from assigned sub-phase bin enrichment levels and store these as a member variable
            unzip_subphase_bin_enrichment_into_element_enrichment(tParentElementsInSupport,
                                                                  tSubphaseBinIndexToCMBinIndex,
                                                                  tSubPhaseBinEnrichment,
                                                                  aCutMesh,
                                                                  aBackgroundMesh,
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
    Cell<moris::Mat_New<Integer, Integer_Matrix>>
    generate_pruned_element_graph_in_basis_support(Integer const &                                              aNumFacePerElement,
                                                    moris::Mat_New<Integer, Integer_Matrix> const &                  aElementsInSupport,
                                                    Cut_Mesh<Real, Integer, Real_Matrix, Integer_Matrix> &       aCutMesh,
                                                    XTK_Mesh<Real, Integer, Real_Matrix, Integer_Matrix> &       aBackgroundMesh)
     {
        mesh::Mesh_Data<Real, Integer, Real_Matrix, Integer_Matrix> & tXTKMeshData = aBackgroundMesh.get_mesh_data();
        Integer tMax = std::numeric_limits<Integer>::max();

        // Construct full element neighbor graph in support and the corresponding shared faces
        moris::Mat_New<Integer, Integer_Matrix> tElementGraph(aElementsInSupport.n_cols(), aNumFacePerElement,tMax);
        moris::Mat_New<Integer, Integer_Matrix> tSharedFaces(aElementsInSupport.n_cols(), aNumFacePerElement,tMax);

        for(Integer iE = 0; iE<aElementsInSupport.n_cols(); iE++)
        {
            // Get elements connected to element and the corresponding face
            moris::Mat_New<Integer, Integer_Matrix> tSingleElementToElement = tXTKMeshData.get_element_connected_to_element_loc_inds(aElementsInSupport(0,iE));
            replace_row(0,tSingleElementToElement,iE,tElementGraph,false);
            replace_row(1,tSingleElementToElement,iE,tSharedFaces,false);
        }


        // prune the graph to only include elements in the support
        // Note cell(0) is the element graph and cell(1) is the faces
        Cell<moris::Mat_New<Integer, Integer_Matrix>> tPrunedNeighbors = prune_element_to_element(tElementGraph,aElementsInSupport,tSharedFaces,tMax);
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
    moris::Mat_New<Integer, Integer_Matrix>
    construct_subphase_bin_neighborhood(moris::Mat_New<Integer, Integer_Matrix> const &                  aParentElementsInSupport,
                                        Integer const &                                              aNumSubPhaseBins,
                                        std::unordered_map<Integer,Integer> &                        aSubphaseBinIndexToCMBinIndex,
                                        moris::Mat_New<Integer, Integer_Matrix> const &                  aPrunedElementGraph,
                                        moris::Mat_New<Integer, Integer_Matrix> const &                  aPrunedSharedFaces,
                                        Cut_Mesh<Real, Integer, Real_Matrix, Integer_Matrix> &       aCutMesh,
                                        XTK_Mesh<Real, Integer, Real_Matrix, Integer_Matrix> &       aBackgroundMesh)
    {
        moris::Mat_New<Integer, Integer_Matrix> tSubphaseBinToSubphaseBin(aNumSubPhaseBins,aNumSubPhaseBins-1,INTEGER_MAX);
        moris::Mat_New<Integer, Integer_Matrix> tSubphaseBinCounter(1,aNumSubPhaseBins,0);

        // Iterate through elements in support, constructing shared element graph using function generate shared_face_element_graph
        for(Integer i = 0; i<aParentElementsInSupport.n_cols(); i++)
        {
            // Iterate through neighbors
            Integer tElementIndex0 = aParentElementsInSupport(0,i);

            for(Integer j = 0; j<aPrunedElementGraph.n_cols(); j++)
            {
                // If there is an integer max at this entry, this all neighbors have been iterated over
                if( aPrunedElementGraph(i,j) == INTEGER_MAX)
                {
                    break;
                }

                // Second background  element on boundary
                Integer tElementIndex1 = aPrunedElementGraph(i,j);

                // The shared parent face index
                Integer tSharedFaceIndex = aPrunedSharedFaces(i,j);

                // Since both element are neighbors of each other, only create this relationship from the lowest to highest element
                // rather than repeat this twice
                if(tElementIndex0<tElementIndex1)
                {
                    bool tElement0HasChildren = aBackgroundMesh.entity_has_children(tElementIndex0, EntityRank::ELEMENT);
                    bool tElement1HasChildren = aBackgroundMesh.entity_has_children(tElementIndex1, EntityRank::ELEMENT);

                    // If both parent elements have children, then use the 2 child on the interface routine
                    if(tElement0HasChildren && tElement1HasChildren)
                    {
                        construct_subphase_bin_to_subphase_bin_2_child_interface(tElementIndex0,            tElementIndex1,
                                                                                 tSharedFaceIndex,          aSubphaseBinIndexToCMBinIndex,
                                                                                 tSubphaseBinToSubphaseBin, tSubphaseBinCounter,
                                                                                 aCutMesh,                  aBackgroundMesh);
                    }

                    // If element 0 does not have children but element 1 does have children
                    // use the mixed child parent on interface routine
                    else if (!tElement0HasChildren && tElement1HasChildren)
                    {
                        construct_subphase_bin_to_subphase_bin_mixed_interface(tElementIndex1,            tElementIndex0,
                                                                               tSharedFaceIndex,          aSubphaseBinIndexToCMBinIndex,
                                                                               tSubphaseBinToSubphaseBin, tSubphaseBinCounter,
                                                                               aCutMesh,                  aBackgroundMesh);
                    }


                    // If element 1 does not have children but element 0 does have children
                    // use the mixed child parent on interface routine
                    else if (tElement0HasChildren && !tElement1HasChildren)
                    {
                        construct_subphase_bin_to_subphase_bin_mixed_interface(tElementIndex0,            tElementIndex1,
                                                                               tSharedFaceIndex,          aSubphaseBinIndexToCMBinIndex,
                                                                               tSubphaseBinToSubphaseBin, tSubphaseBinCounter,
                                                                               aCutMesh,                  aBackgroundMesh);
                    }

                    // Both elements do not have children
                    // use the parent interface routing
                    else
                    {
                        construct_subphase_bin_to_subphase_bin_2_parent_interface(tElementIndex0,
                                                                                  tElementIndex1,
                                                                                  aSubphaseBinIndexToCMBinIndex,
                                                                                  tSubphaseBinToSubphaseBin,
                                                                                  tSubphaseBinCounter,
                                                                                  aCutMesh,
                                                                                  aBackgroundMesh);
                    }
                }
            }


        }

        return tSubphaseBinToSubphaseBin;
    }


    void
    construct_subphase_bin_to_subphase_bin_2_child_interface(
            Integer const &                                              aParentElementIndex0,
            Integer const &                                              aParentElementIndex1,
            Integer const &                                              aSharedFaceIndex,
            std::unordered_map<Integer,Integer> &                        aSubphaseBinIndexToCMBinIndex,
            moris::Mat_New<Integer, Integer_Matrix> &                        aSubphaseBinToSubphaseBin,
            moris::Mat_New<Integer, Integer_Matrix> &                        aSubphaseBinCounter,
            Cut_Mesh<Real, Integer, Real_Matrix, Integer_Matrix> &       aCutMesh,
            XTK_Mesh<Real, Integer, Real_Matrix, Integer_Matrix> &       aBackgroundMesh)
    {
        // Get the child mesh index
        Integer tChildMeshIndex0 = aBackgroundMesh.child_mesh_index(aParentElementIndex0,EntityRank::ELEMENT);
        Integer tChildMeshIndex1 = aBackgroundMesh.child_mesh_index(aParentElementIndex1,EntityRank::ELEMENT);

        // Get the child me shes on this boundary
        Child_Mesh_Test<Real, Integer, Real_Matrix, Integer_Matrix> tChildMesh0 = aCutMesh.get_child_mesh(tChildMeshIndex0);
        Child_Mesh_Test<Real, Integer, Real_Matrix, Integer_Matrix> tChildMesh1 = aCutMesh.get_child_mesh(tChildMeshIndex1);

        // Get child element subphase bin membership
        moris::Mat_New<Integer, Integer_Matrix> const & tChildElements0BinMembership = tChildMesh0.get_elemental_subphase_bin_membership();
        moris::Mat_New<Integer, Integer_Matrix> const & tChildElements1BinMembership = tChildMesh1.get_elemental_subphase_bin_membership();

        // Construct element pairs across shared parent face
        moris::Mat_New<Integer, Integer_Matrix> tBoundaryElementPairs =
            generate_shared_face_element_pairs(aSharedFaceIndex,tChildMeshIndex0,tChildMeshIndex1,aCutMesh);

        // iterate over pairs and create a relationship between their elements buckets
        for(Integer k = 0; k<tBoundaryElementPairs.n_cols(); k++)
        {
            // Get pairs children element bin membership
            Integer tChildElemBin0 = tChildElements0BinMembership(0,(tBoundaryElementPairs)(0,k));
            Integer tChildElemBin1 = tChildElements1BinMembership(0,(tBoundaryElementPairs)(1,k));

            // get these bins basis cluster index
            Integer tChildElemBinBasisIndex0 = aSubphaseBinIndexToCMBinIndex[cantor_pairing(tChildMeshIndex0+1,tChildElemBin0)];
            Integer tChildElemBinBasisIndex1 = aSubphaseBinIndexToCMBinIndex[cantor_pairing(tChildMeshIndex1+1,tChildElemBin1)];

            //Check whether this bin relationship already exists
            bool tRelationshipExists = false;
            for(Integer l = 0; l <= aSubphaseBinCounter(0,tChildElemBinBasisIndex0); l++)
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
            Integer const &                                              aParentElementWithChildren,
            Integer const &                                              aParentElementWithoutChildren,
            Integer const &                                              aSharedFaceIndex,
            std::unordered_map<Integer,Integer> &                        aSubphaseBinIndexToCMBinIndex,
            moris::Mat_New<Integer, Integer_Matrix> &                        aSubphaseBinToSubphaseBin,
            moris::Mat_New<Integer, Integer_Matrix> &                        aSubphaseBinCounter,
            Cut_Mesh<Real, Integer, Real_Matrix, Integer_Matrix> &       aCutMesh,
            XTK_Mesh<Real, Integer, Real_Matrix, Integer_Matrix> &       aBackgroundMesh)
    {
        // Get the child mesh index
        Integer tChildMeshIndex = aBackgroundMesh.child_mesh_index(aParentElementWithChildren,EntityRank::ELEMENT);

        // Get the child mesh on this shared face
        Child_Mesh_Test<Real, Integer, Real_Matrix, Integer_Matrix> tChildMesh = aCutMesh.get_child_mesh(tChildMeshIndex);

        // Get child element subphase bin membership
        moris::Mat_New<Integer, Integer_Matrix> const & tChildElementsBinMembership = tChildMesh.get_elemental_subphase_bin_membership();

        // Allocate Matrixes
        moris::Mat_New<Integer, Integer_Matrix> tFaceOrdinals(1,1);
        moris::Mat_New<Integer, Integer_Matrix> tChildrenElementCMInds(1,1);
        moris::Mat_New<Integer, Integer_Matrix> tChildrenElementIds(1,1);

        // Get children elements attached to aFaceIndex on the side of child mesh index 0
        aCutMesh.get_child_elements_connected_to_parent_face(tChildMeshIndex,
                                                             aSharedFaceIndex,
                                                             tChildrenElementIds,
                                                             tChildrenElementCMInds,
                                                             tFaceOrdinals);


        // iterate over child elements on boundary and construct the subphase bin relationship
        for(Integer k = 0; k<tChildrenElementCMInds.n_cols(); k++)
        {
            // Get pairs children element bin membership
            Integer tChildElemBin = tChildElementsBinMembership(0,tChildrenElementCMInds(0,k));

            // Get the child  element bin basis index
            Integer tChildElemBinBasisIndex = aSubphaseBinIndexToCMBinIndex[cantor_pairing(tChildMeshIndex+1,tChildElemBin)];

            // Get the child  element bin basis index
            Integer tParentElementDummyIndex = 0;
            Integer tParentElemBinBasisIndex = aSubphaseBinIndexToCMBinIndex[cantor_pairing(tParentElementDummyIndex,aParentElementWithoutChildren)];

            //Check whether this bin relationship already exists
            bool tRelationshipExists = false;
            for(Integer l = 0; l <= aSubphaseBinCounter(0,tChildElemBinBasisIndex); l++)
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
            Integer const &                                              aParentElementIndex0,
            Integer const &                                              aParentElementIndex1,
            std::unordered_map<Integer,Integer> &                        aSubphaseBinIndexToCMBinIndex,
            moris::Mat_New<Integer, Integer_Matrix> &                        aSubphaseBinToSubphaseBin,
            moris::Mat_New<Integer, Integer_Matrix> &                        aSubphaseBinCounter,
            Cut_Mesh<Real, Integer, Real_Matrix, Integer_Matrix> &       aCutMesh,
            XTK_Mesh<Real, Integer, Real_Matrix, Integer_Matrix> &       aBackgroundMesh)
    {

            // get these bins basis cluster index
        Integer tParentElementDummyIndex = 0;
        Integer tParentElemBinBasisIndex0 = aSubphaseBinIndexToCMBinIndex[cantor_pairing(tParentElementDummyIndex,aParentElementIndex0)];
        Integer tParentElemBinBasisIndex1 = aSubphaseBinIndexToCMBinIndex[cantor_pairing(tParentElementDummyIndex,aParentElementIndex1)];

        //Check whether this bin relationship already exists
        bool tRelationshipExists = false;

        for(Integer l = 0; l <= aSubphaseBinCounter(0,tParentElemBinBasisIndex0); l++)
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
    moris::Mat_New<Integer, Integer_Matrix>
    assign_subphase_bin_enrichment_levels_in_basis_support(moris::Mat_New<Integer, Integer_Matrix> const &                  aParentElementsInSupport,
                                                           Integer const &                                              aNumSubPhaseBins,
                                                           moris::Mat_New<Integer, Integer_Matrix> const &                  aSubPhaseBinBulkPhase,
                                                           std::unordered_map<Integer,Integer> &                        aSubphaseBinIndexToCMBinIndex,
                                                           moris::Mat_New<Integer, Integer_Matrix> const &                  aPrunedElementGraph,
                                                           moris::Mat_New<Integer, Integer_Matrix> const &                  aPrunedSharedFaces,
                                                           Cut_Mesh<Real, Integer, Real_Matrix, Integer_Matrix> &       aCutMesh,
                                                           XTK_Mesh<Real, Integer, Real_Matrix, Integer_Matrix> &       aBackgroundMesh)
    {

        // Initialize Enrichment Vector now that the count includes all children
        moris::Mat_New<Integer, Integer_Matrix> tSubPhaseBinEnrichmentVals(1, aNumSubPhaseBins);

        // Consider cutting down the number of arguements by moving aBackgroundMesh,aCutMesh,aMatrixFactory into reference member variables
        moris::Mat_New<Integer, Integer_Matrix> tSubPhaseBinNeighborhood =
        construct_subphase_bin_neighborhood(aParentElementsInSupport,
                                            aNumSubPhaseBins,
                                            aSubphaseBinIndexToCMBinIndex,
                                            aPrunedElementGraph,
                                            aPrunedSharedFaces,
                                            aCutMesh,
                                            aBackgroundMesh);

        // Variables needed for floodfill, consider removing these.
        // Active bins to include in floodfill (We include all bins)
        moris::Mat_New<Integer, Integer_Matrix> tActiveBins(1,tSubPhaseBinNeighborhood.n_rows());
        for(Integer i = 0; i< tSubPhaseBinNeighborhood.n_rows(); i++)
        {
            (tActiveBins)(0,i) = i;
        }

        // Mark all as included
        moris::Mat_New<Integer, Integer_Matrix> tIncludedBins(1,tSubPhaseBinNeighborhood.n_rows(),1);

        // Dummy value
        Integer tMax = INTEGER_MAX;

        tSubPhaseBinEnrichmentVals = flood_fill(tSubPhaseBinNeighborhood,
                                                aSubPhaseBinBulkPhase,
                                                tActiveBins,
                                                tIncludedBins,
                                                mNumBulkPhases,
                                                tMax,
                                                true);


        return tSubPhaseBinEnrichmentVals;
    }

    Integer
    setup_all_subphase_bins_in_basis_support(moris::Mat_New<Integer, Integer_Matrix> const &           aParentElementsInSupport,
                                            Cut_Mesh<Real, Integer, Real_Matrix, Integer_Matrix> & aCutMesh,
                                            XTK_Mesh<Real, Integer, Real_Matrix, Integer_Matrix> & aBackgroundMesh,
                                            std::unordered_map<Integer,Integer> &                  aSubphaseBinIndexToCMBinIndex,
                                            moris::Mat_New<Integer, Integer_Matrix> &                  aSubPhaseBinBulkPhase)
    {
        // Count all elements including children of the elements in support of the basis
        Integer tNumSubPhaseBinsInSupport = count_subphase_bins_in_support(aParentElementsInSupport,aBackgroundMesh,aCutMesh);

        aSubPhaseBinBulkPhase.resize(1,tNumSubPhaseBinsInSupport);

        // Counter
        Integer tCount = 0;

        // Child Mesh Index
        Integer tChildMeshIndex = 0;

        // Number of children elements
        Integer tNumChildrenElements = 0;

        // Bin Index
        Integer tBinIndex = 0;

        for(Integer i = 0; i<aParentElementsInSupport.n_cols(); i++)
        {
            // Add children elements to all elements in support vector but do not add parent
            if(aBackgroundMesh.entity_has_children(aParentElementsInSupport(0,i),EntityRank::ELEMENT))
            {
                tChildMeshIndex = aBackgroundMesh.child_mesh_index(aParentElementsInSupport(0,i),EntityRank::ELEMENT);

                // Get the child mesh
                Child_Mesh_Test<Real, Integer, Real_Matrix, Integer_Matrix> & tChildMesh = aCutMesh.get_child_mesh(tChildMeshIndex);

                Integer tNumLocSubPhaseBins = tChildMesh.get_num_subphase_bins();

                Cell<Integer> const & tSubPhaseBinBulkPhase = tChildMesh.get_subphase_bin_bulk_phase();

                for(Integer j = 0; j<tNumLocSubPhaseBins; j++)
                {

                    // Create a hash for child mesh index and bin index
                    Integer tCMIndexBinHash = cantor_pairing(tChildMeshIndex+1,j);

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
                Integer tChildMeshIndex = 0;

                // hash a child mesh index of 0 and parent element index
                Integer tHash = cantor_pairing(tChildMeshIndex,aParentElementsInSupport(0,i));

                aSubphaseBinIndexToCMBinIndex[tHash] = tBinIndex;

                aSubPhaseBinBulkPhase(0,tBinIndex) = aBackgroundMesh.get_element_phase_index(aParentElementsInSupport(0,i));

                tBinIndex++;
            }
        }

        return tNumSubPhaseBinsInSupport;
    }


    void
    unzip_subphase_bin_enrichment_into_element_enrichment(moris::Mat_New<Integer, Integer_Matrix> const &                  aParentElementsInSupport,
                                                          std::unordered_map<Integer,Integer> &                        aSubphaseBinIndexToCMBinIndex,
                                                          moris::Mat_New<Integer, Integer_Matrix> const &                  aSubPhaseBinEnrichmentLevel,
                                                          Cut_Mesh<Real, Integer, Real_Matrix, Integer_Matrix> &       aCutMesh,
                                                          XTK_Mesh<Real, Integer, Real_Matrix, Integer_Matrix> &       aBackgroundMesh,
                                                          moris::Mat_New<Integer, Integer_Matrix> &       aElementIndInBasisSupport,
                                                          moris::Mat_New<Integer, Integer_Matrix> &       aElementEnrichmentLevel)
    {

        // Count all elements including children of the elements in support of the basis
        Integer tNumAllElementsInSupport = count_elements_in_support(aParentElementsInSupport,aBackgroundMesh,aCutMesh);

        // Background mesh underlying meshd ata
        mesh::Mesh_Data<Real, Integer, Real_Matrix, Integer_Matrix> const & tBackgroundMeshData = aBackgroundMesh.get_mesh_data();

        aElementIndInBasisSupport = moris::Mat_New<Integer, Integer_Matrix>(1,tNumAllElementsInSupport);
        aElementEnrichmentLevel   = moris::Mat_New<Integer, Integer_Matrix>(1,tNumAllElementsInSupport);

        // Counter
        Integer tCount = 0;

        // Child Mesh Index
        Integer tChildMeshIndex = 0;

        // Number of children elements
        Integer tNumChildrenElements = 0;

        for(Integer i = 0; i<aParentElementsInSupport.n_cols(); i++)
        {
            // Add children elements to all elements in support vector but do not add parent
            if(aBackgroundMesh.entity_has_children(aParentElementsInSupport(0,i),EntityRank::ELEMENT))
            {
                Integer tChildMeshIndex = aBackgroundMesh.child_mesh_index(aParentElementsInSupport(0,i),EntityRank::ELEMENT);

                Child_Mesh_Test<Real, Integer, Real_Matrix, Integer_Matrix> & tChildMesh = aCutMesh.get_child_mesh(tChildMeshIndex);

                moris::Mat_New<Integer, Integer_Matrix> const & tChildElementSubphaseBin = tChildMesh.get_elemental_subphase_bin_membership();
                moris::Mat_New<Integer, Integer_Matrix> const & tChildElementIds = tChildMesh.get_element_ids();
                tNumChildrenElements = tChildElementIds.n_cols();

                for(Integer j = 0; j<tNumChildrenElements; j++)
                {
                    aElementIndInBasisSupport(0,tCount) = tChildElementIds(0,j);

                    Integer tChildElementBinHash = cantor_pairing(tChildMeshIndex+1,tChildElementSubphaseBin(0,j));

                    Integer tBinBasisIndex = aSubphaseBinIndexToCMBinIndex[tChildElementBinHash];

                    aElementEnrichmentLevel(0,tCount) = aSubPhaseBinEnrichmentLevel(0,tBinBasisIndex);

                    tCount++;
                }

            }

            else
            {
                Integer tParentElementDummyIndex = 0;
                Integer tParentHash = cantor_pairing(tParentElementDummyIndex,aParentElementsInSupport(0,i));
                Integer tBinBasisIndex = aSubphaseBinIndexToCMBinIndex[tParentHash];


                aElementIndInBasisSupport(0,tCount) = tBackgroundMeshData.get_glb_entity_id_from_entity_loc_index(aParentElementsInSupport(0,i),EntityRank::ELEMENT);
                aElementEnrichmentLevel(0,tCount) = aSubPhaseBinEnrichmentLevel(0,tBinBasisIndex);

                tCount++;
            }

        }
    }



    Integer
    count_subphase_bins_in_support(moris::Mat_New<Integer, Integer_Matrix> const &            aParentElementsInSupport,
                                   XTK_Mesh<Real, Integer, Real_Matrix, Integer_Matrix> & aBackgroundMesh,
                                   Cut_Mesh<Real, Integer, Real_Matrix, Integer_Matrix> & aCutMesh)
    {
        // Number of elements in this support (need both parent and total)
        Integer tNumParentElementsInSupport = aParentElementsInSupport.n_cols();

        // Initialize number of sub phase bins
        Integer tNumSubPhaseBins  = 0;

        // initialize variable for child mesh index if an element has children
        Integer tChildMeshIndex = 0;

        // initialize variable for number of children elements in a mesh
        Integer tNumChildElements = 0;

        // Count children elements in support
        for(Integer i = 0; i<tNumParentElementsInSupport; i++)
        {
            // Check if this element has children and if it does add them to the count
            if(aBackgroundMesh.entity_has_children(aParentElementsInSupport(0,i),EntityRank::ELEMENT))
            {
                // The child mesh index
                tChildMeshIndex = aBackgroundMesh.child_mesh_index(aParentElementsInSupport(0,i),EntityRank::ELEMENT);

                // Get the child mesh
                Child_Mesh_Test<Real, Integer, Real_Matrix, Integer_Matrix> & tChildMesh = aCutMesh.get_child_mesh(tChildMeshIndex);

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
    Integer
    count_elements_in_support(moris::Mat_New<Integer, Integer_Matrix> const &            aParentElementsInSupport,
                              XTK_Mesh<Real, Integer, Real_Matrix, Integer_Matrix> & aBackgroundMesh,
                              Cut_Mesh<Real, Integer, Real_Matrix, Integer_Matrix> & aCutMesh)
    {

        // Number of elements in this support (need both parent and total)
        Integer tNumParentElementsInSupport = aParentElementsInSupport.n_cols();
        Integer tNumElementsInSupport = tNumParentElementsInSupport;

        // initialize variable for child mesh index if an element has children
        Integer tChildMeshIndex = 0;

        // initialize variable for number of children elements in a mesh
        Integer tNumChildElements = 0;

        // Count children elements in support
        for(Integer i = 0; i<tNumParentElementsInSupport; i++)
        {
            // Check if this element has children and if it does add them to the count
            if(aBackgroundMesh.entity_has_children(aParentElementsInSupport(0,i),EntityRank::ELEMENT))
            {
                // The child mesh index
                tChildMeshIndex = aBackgroundMesh.child_mesh_index(aParentElementsInSupport(0,i),EntityRank::ELEMENT);

                // Number of child elements in this mesh
                tNumChildElements = aCutMesh.get_num_entities(tChildMeshIndex,EntityRank::ELEMENT);

                // Add the number of elements in the child mesh to the account (-1 to remove parent)
                tNumElementsInSupport = tNumElementsInSupport + tNumChildElements - 1;
            }
        }

        return tNumElementsInSupport;
    }





};
}
#endif /* XTK_SRC_XTK_CL_XTK_ENRICHMENT_HPP_ */
