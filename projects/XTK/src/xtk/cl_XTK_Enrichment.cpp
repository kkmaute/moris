/*
 * cl_XTK_Enrichment.cpp
 *
 *  Created on: Feb 18, 2019
 *      Author: doble
 */

#include "cl_XTK_Enrichment.hpp"
#include "cl_XTK_Hexahedron_8_Basis_Function.hpp"
#include "fn_sum.hpp"
#include "fn_equal_to.hpp"
namespace xtk
{

Enrichment::Enrichment(moris::size_t         aNumBulkPhases,
                       xtk::Model*           aXTKModelPtr,
                       xtk::Cut_Mesh*        aCutMeshPtr,
                       xtk::Background_Mesh* aBackgroundMeshPtr):
    mNumBulkPhases(aNumBulkPhases),
    mXTKModelPtr(aXTKModelPtr),
    mCutMeshPtr(aCutMeshPtr),
    mBackgroundMeshPtr(aBackgroundMeshPtr)
{

}



void
Enrichment::perform_enrichment()
{
    // Start clock
    std::clock_t start = std::clock();

    // Verify initialized properly
    MORIS_ERROR(mCutMeshPtr!=nullptr,"mCutMesh nullptr detected, this is probably because the enrichment has not been initialized properly");
    MORIS_ERROR(mBackgroundMeshPtr!=nullptr,"mBackgroundMesh nullptr detected, this is probably because the enrichment has not been initialized properly");

    // Perform local enrichment for each child mesh (commits local floodfill data to child mesh)
    perform_local_subphase_identification();

    // Perform enrichment over basis clusters
    perform_basis_cluster_enrichment();


    // Output time
    if(moris::par_rank() == 0 && mVerbose)
    {
        std::cout<<"XTK: Enrichment completed in "<< (std::clock() - start) / (double)(CLOCKS_PER_SEC)<<" s."<<std::endl;
    }
}


Cell<moris::Matrix< moris::IdMat >> const &
Enrichment::get_element_inds_in_basis_support() const
{
    return mElementIndsInBasis;
}


Cell<moris::Matrix< moris::IndexMat >> const &
Enrichment::get_element_enrichment_levels_in_basis_support() const
{
    return mElementEnrichmentLevel;
}


void
Enrichment::perform_local_subphase_identification()
{

    // get the number of children meshes
    moris::size_t tNumChildMeshes = mCutMeshPtr->get_num_child_meshes();

    // iterate over children meshes and perform local flood-fill
    for(moris::size_t i = 0; i<tNumChildMeshes; i++)
    {
        // Get child mesh index
        Child_Mesh & tChildMesh = mCutMeshPtr->get_child_mesh(i);

        // Perform local flood-fill on child mesh to identify subphase
        moris::Matrix< moris::IndexMat > tLocalFloodFill = local_child_mesh_flood_fill(tChildMesh);

        // Set the local floodfill data as the elemental subphase values in the child mesh
        // The child mesh then sorts the elements into bins
        tChildMesh.set_elemental_subphase(tLocalFloodFill);
    }
}

void
Enrichment::perform_basis_cluster_enrichment()
{
    // Get underlying matrix data to access function
    moris::mtk::Mesh & tXTKMeshData = mBackgroundMeshPtr->get_mesh_data();

    // Number of basis functions
    moris::size_t tNumBasis          = tXTKMeshData.get_num_basis_functions();
    moris::size_t tNumFacePerElement = tXTKMeshData.get_entity_connected_to_entity_loc_inds(0, moris::EntityRank::ELEMENT, moris::EntityRank::FACE).n_cols();

    // Allocate member variables
    mElementEnrichmentLevel = Cell<moris::Matrix< moris::IndexMat >>(tNumBasis);
    mElementIndsInBasis      = Cell<moris::Matrix< moris::IndexMat >>(tNumBasis);

    for(moris::size_t i = 0; i<tNumBasis; i++)
    {

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
                                                              mElementIndsInBasis(i),
                                                              mElementEnrichmentLevel(i));

    }

    setup_vertex_enrichment_data();

    // assign enrichment level indices
    assign_enrichment_level_identifiers();

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
Enrichment::generate_pruned_element_graph_in_basis_support(moris::size_t const &                                         aNumFacePerElement,
                                                moris::Matrix< moris::IndexMat > const &               aElementsInSupport)
 {
    moris::mtk::Mesh & tXTKMeshData = mBackgroundMeshPtr->get_mesh_data();
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
Enrichment::construct_subphase_bin_neighborhood(moris::Matrix< moris::IndexMat > const &        aParentElementsInSupport,
                                    moris::size_t const &                                       aNumSubPhaseBins,
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
                bool tElement0HasChildren = mBackgroundMeshPtr->entity_has_children(tElementIndex0, EntityRank::ELEMENT);
                bool tElement1HasChildren = mBackgroundMeshPtr->entity_has_children(tElementIndex1, EntityRank::ELEMENT);

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
Enrichment::construct_subphase_bin_to_subphase_bin_2_child_interface(
        moris::moris_index const &                                  aParentElementIndex0,
        moris::moris_index const &                                  aParentElementIndex1,
        moris::moris_index const &                                  aSharedFaceIndex,
        std::unordered_map<moris::moris_index,moris::moris_index> & aSubphaseBinIndexToCMBinIndex,
        moris::Matrix< moris::IndexMat > &                          aSubphaseBinToSubphaseBin,
        moris::Matrix< moris::DDSTMat >  &                          aSubphaseBinCounter)
{
    // Get the child mesh index
    moris::moris_index tChildMeshIndex0 = mBackgroundMeshPtr->child_mesh_index(aParentElementIndex0,EntityRank::ELEMENT);
    moris::moris_index tChildMeshIndex1 = mBackgroundMeshPtr->child_mesh_index(aParentElementIndex1,EntityRank::ELEMENT);

    // Get the child meshes on this boundary
    Child_Mesh tChildMesh0 = mCutMeshPtr->get_child_mesh(tChildMeshIndex0);
    Child_Mesh tChildMesh1 = mCutMeshPtr->get_child_mesh(tChildMeshIndex1);

    // Get child element subphase bin membership
    moris::Matrix< moris::IndexMat > const & tChildElements0BinMembership = tChildMesh0.get_elemental_subphase_bin_membership();
    moris::Matrix< moris::IndexMat > const & tChildElements1BinMembership = tChildMesh1.get_elemental_subphase_bin_membership();

    // Construct element pairs across shared parent face
    moris::Matrix< moris::IndexMat > tBoundaryElementPairs =
        generate_shared_face_element_pairs(aSharedFaceIndex,tChildMeshIndex0,tChildMeshIndex1,*mCutMeshPtr);

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
Enrichment::construct_subphase_bin_to_subphase_bin_mixed_interface(
        moris::moris_index const &                                   aParentElementWithChildren,
        moris::moris_index const &                                   aParentElementWithoutChildren,
        moris::moris_index const &                                   aSharedFaceIndex,
        std::unordered_map<moris::moris_index,moris::moris_index> &  aSubphaseBinIndexToCMBinIndex,
        moris::Matrix< moris::IndexMat > &                           aSubphaseBinToSubphaseBin,
        moris::Matrix< moris::DDSTMat > &                            aSubphaseBinCounter)
{

    // Get the child mesh index
    moris::moris_index tChildMeshIndex = mBackgroundMeshPtr->child_mesh_index(aParentElementWithChildren,EntityRank::ELEMENT);

    // Get the child mesh on this shared face
    Child_Mesh tChildMesh = mCutMeshPtr->get_child_mesh(tChildMeshIndex);

    // Get child element subphase bin membership
    moris::Matrix< moris::IndexMat > const & tChildElementsBinMembership = tChildMesh.get_elemental_subphase_bin_membership();

    // Allocate Matrixes
    moris::Matrix< moris::IndexMat > tFaceOrdinals(1,1);
    moris::Matrix< moris::IndexMat > tChildrenElementCMInds(1,1);
    moris::Matrix< moris::IdMat > tChildrenElementIds(1,1);

    // Get children elements attached to aFaceIndex on the side of child mesh index 0
    mCutMeshPtr->get_child_elements_connected_to_parent_face(tChildMeshIndex,
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
Enrichment::construct_subphase_bin_to_subphase_bin_2_parent_interface(
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
Enrichment::assign_subphase_bin_enrichment_levels_in_basis_support(moris::Matrix< moris::IndexMat > const &                     aParentElementsInSupport,
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
Enrichment::setup_all_subphase_bins_in_basis_support(moris::Matrix< moris::IndexMat > const &               aParentElementsInSupport,
                                         std::unordered_map<moris::moris_index,moris::moris_index> & aSubphaseBinIndexToCMBinIndex,
                                         moris::Matrix< moris::IndexMat > &                  aSubPhaseBinBulkPhase)
{
    // Count all elements including children of the elements in support of the basis
    moris::size_t tNumSubPhaseBinsInSupport = count_subphase_bins_in_support(aParentElementsInSupport);

    aSubPhaseBinBulkPhase.resize(1,tNumSubPhaseBinsInSupport);

    // Child Mesh Index
    moris::size_t tChildMeshIndex = 0;

    // Bin Index
    moris::size_t tBinIndex = 0;

    for(moris::size_t i = 0; i<aParentElementsInSupport.n_cols(); i++)
    {
        // Add children elements to all elements in support vector but do not add parent
        if(mBackgroundMeshPtr->entity_has_children(aParentElementsInSupport(0,i),EntityRank::ELEMENT))
        {
            tChildMeshIndex = mBackgroundMeshPtr->child_mesh_index(aParentElementsInSupport(0,i),EntityRank::ELEMENT);

            // Get the child mesh
            Child_Mesh & tChildMesh = mCutMeshPtr->get_child_mesh(tChildMeshIndex);

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

            aSubPhaseBinBulkPhase(0,tBinIndex) = mBackgroundMeshPtr->get_element_phase_index(aParentElementsInSupport(0,i));

            tBinIndex++;
        }
    }

    return tNumSubPhaseBinsInSupport;
}


void
Enrichment::unzip_subphase_bin_enrichment_into_element_enrichment(moris::Matrix< moris::IndexMat > const &                     aParentElementsInSupport,
                                                      std::unordered_map<moris::moris_index, moris::moris_index> & aSubphaseBinIndexToCMBinIndex,
                                                      moris::Matrix< moris::IndexMat > const &                      aSubPhaseBinEnrichmentLevel,
                                                      moris::Matrix< moris::IndexMat > &       aElementIndInBasisSupport,
                                                      moris::Matrix< moris::IndexMat > &       aElementEnrichmentLevel)
{

    // Count all elements including children of the elements in support of the basis
    moris::size_t tNumAllElementsInSupport = count_elements_in_support(aParentElementsInSupport);

    aElementIndInBasisSupport = moris::Matrix< moris::IndexMat >(1,tNumAllElementsInSupport);
    aElementEnrichmentLevel   = moris::Matrix< moris::IndexMat >(1,tNumAllElementsInSupport);

    // Counter
    moris::size_t tCount = 0;

    // Number of children elements
    moris::size_t tNumChildrenElements = 0;

    for(moris::size_t i = 0; i<aParentElementsInSupport.n_cols(); i++)
    {
        // Add children elements to all elements in support vector but do not add parent
        if(mBackgroundMeshPtr->entity_has_children(aParentElementsInSupport(0,i),EntityRank::ELEMENT))
        {
            moris::moris_index tChildMeshIndex = mBackgroundMeshPtr->child_mesh_index(aParentElementsInSupport(0,i),EntityRank::ELEMENT);

            Child_Mesh & tChildMesh = mCutMeshPtr->get_child_mesh(tChildMeshIndex);

            moris::Matrix< moris::IndexMat > const & tChildElementSubphaseBin = tChildMesh.get_elemental_subphase_bin_membership();
            moris::Matrix< moris::IdMat > const & tChildElementInds = tChildMesh.get_element_inds();
            tNumChildrenElements = tChildElementInds.n_cols();

            for(moris::size_t j = 0; j<tNumChildrenElements; j++)
            {
                aElementIndInBasisSupport(0,tCount) = tChildElementInds(0,j);

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


            aElementIndInBasisSupport(0,tCount) = aParentElementsInSupport(0,i);
            aElementEnrichmentLevel(0,tCount) = aSubPhaseBinEnrichmentLevel(0,tBinBasisIndex);

            tCount++;
        }

    }
}


void
Enrichment::setup_vertex_enrichment_data()
{
    // Allocate member data created by this routine
    moris::uint tNumNodes = mBackgroundMeshPtr->get_num_entities(EntityRank::NODE);
    mVertexEnrichments = moris::Cell<xtk::Vertex_Enrichment>(tNumNodes);

    // get element to basis connectivity
    construct_element_to_basis_connectivity(mElementToBasis,mElementToBasisEnrichmentLevel);

    // get node to basis connectivity
    construct_vertex_enrichment_with_element_to_basis(mElementToBasis, mElementToBasisEnrichmentLevel);

    //compute the basis weights and store them
    compute_vertex_basis_weights();

}


moris::size_t
Enrichment::count_subphase_bins_in_support(moris::Matrix< moris::IndexMat > const & aParentElementsInSupport)
{
    // Number of elements in this support (need both parent and total)
    moris::size_t tNumParentElementsInSupport = aParentElementsInSupport.n_cols();

    // Initialize number of sub phase bins
    moris::size_t tNumSubPhaseBins  = 0;

    // initialize variable for child mesh index if an element has children
    moris::size_t tChildMeshIndex = 0;

    // Count children elements in support
    for(moris::size_t i = 0; i<tNumParentElementsInSupport; i++)
    {
        // Check if this element has children and if it does add them to the count
        if(mBackgroundMeshPtr->entity_has_children(aParentElementsInSupport(0,i),EntityRank::ELEMENT))
        {
            // The child mesh index
            tChildMeshIndex = mBackgroundMeshPtr->child_mesh_index(aParentElementsInSupport(0,i),EntityRank::ELEMENT);

            // Get the child mesh
            Child_Mesh & tChildMesh = mCutMeshPtr->get_child_mesh(tChildMeshIndex);

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

moris::size_t
Enrichment::count_elements_in_support(moris::Matrix< moris::IndexMat > const &            aParentElementsInSupport)
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
        if(mBackgroundMeshPtr->entity_has_children(aParentElementsInSupport(0,i),EntityRank::ELEMENT))
        {
            // The child mesh index
            tChildMeshIndex = mBackgroundMeshPtr->child_mesh_index(aParentElementsInSupport(0,i),EntityRank::ELEMENT);

            // Number of child elements in this mesh
            tNumChildElements = mCutMeshPtr->get_num_entities(tChildMeshIndex,EntityRank::ELEMENT);

            // Add the number of elements in the child mesh to the account (-1 to remove parent)
            tNumElementsInSupport = tNumElementsInSupport + tNumChildElements - 1;
        }
    }

    return tNumElementsInSupport;
}


void
Enrichment::construct_element_to_basis_connectivity(moris::Cell<moris::Cell<moris::moris_index>> & aElementToBasis,
                                                    moris::Cell<moris::Cell<moris::moris_index>> & aElementToBasisEnrichmentLevel)
{
    // mesh data
    moris::mtk::Mesh & tXTKMeshData = mBackgroundMeshPtr->get_mesh_data();

    // Number of basis functions
    moris::size_t tNumBasis = tXTKMeshData.get_num_basis_functions();

    // member data access
    Cell<moris::Matrix< moris::IdMat >> const & tElementInds        = this->get_element_inds_in_basis_support();
    Cell<moris::Matrix< moris::IdMat >> const & tElementEnrichments = this->get_element_enrichment_levels_in_basis_support();

    // allocate outputs
    aElementToBasis = moris::Cell<moris::Cell<moris::moris_index>>(mXTKModelPtr->get_num_elements());
    aElementToBasisEnrichmentLevel = moris::Cell<moris::Cell<moris::moris_index>>(mXTKModelPtr->get_num_elements());

    // Iterate through basis
    for(moris::moris_index  iB = 0; iB<(moris::moris_index)tNumBasis; iB++)
    {
        // iterate through elements (child and not child in basis support)
        for(moris::uint iEl = 0; iEl < tElementInds(iB).numel(); iEl++ )
        {
            moris::moris_index tElemIndex = tElementInds(iB)(iEl);
            moris::moris_index tEnrichmentLevel = tElementEnrichments(iB)(iEl);

            // add basis index to cell for the element and also store enrichment level
            aElementToBasis(tElemIndex).push_back(iB);
            aElementToBasisEnrichmentLevel(tElemIndex).push_back(tEnrichmentLevel);
        }
    }
}

void
Enrichment::construct_vertex_to_basis_connectivity(moris::Cell<moris::Cell<moris::moris_index>> const & aElementToBasis,
                                                 moris::Cell<moris::Cell<moris::moris_index>> const & aElementToBasisEnrichmentLevel)
{
    moris::uint tNumVerts = mBackgroundMeshPtr->get_num_entities(EntityRank::NODE);
    moris::Cell<moris::Cell<moris::moris_index>> tVertexToBasis(tNumVerts);

    for(moris::uint iEl = 0; iEl <aElementToBasis.size(); iEl++)
    {

    }

}

void
Enrichment::construct_vertex_enrichment_with_element_to_basis(moris::Cell<moris::Cell<moris::moris_index>> const & aElementToBasis,
                                                              moris::Cell<moris::Cell<moris::moris_index>> const & aElementToBasisEnrichmentLevel)
{
    // allocate node to basis data
    moris::Cell<moris::Cell<moris::moris_index>> tVertexToBasis(mVertexEnrichments.size());
    moris::Cell<moris::Cell<moris::moris_index>> tNodeToBasisEnrLvl(mVertexEnrichments.size());

    Cell<moris::moris_index> tNodeToBasisCount(mVertexEnrichments.size(),0);

    // iterate through elements
    for(moris::uint iEl = 0; iEl <aElementToBasis.size(); iEl++)
    {
        // number of basis functions interp into element
        moris::uint tNumBasis = aElementToBasis(iEl).size();

        if(tNumBasis >0)
        {
            // get information about the cell
            moris::mtk::Cell & tCell = mBackgroundMeshPtr->get_mtk_cell(iEl);

            // We construct the element to basis first such that we don't have to make this
            // matrix a bunch of times
            moris::Matrix< moris::IndexMat > tCellVertInds = tCell.get_vertex_inds();

            // iterate through basis associated with element
            for(moris::uint iB = 0; iB<tNumBasis; iB++)
            {
                moris::moris_index tBasisIndex = aElementToBasis(iEl)(iB);

                // iterate through nodes attached to element
                for(moris::uint iN = 0; iN < tCellVertInds.numel(); iN++)
                {
                    moris::moris_index tVertInd = tCellVertInds(iN);

                    std::unordered_map<moris::moris_index,moris::moris_index> & tVertEnrichMap = mVertexEnrichments(tVertInd).get_basis_map();

                    auto tIter = tVertEnrichMap.find(tBasisIndex);
                    if(tIter == tVertEnrichMap.end())
                    {
                        // place this basis in the map for the vertex (ensures we don't add it twice
                        // to the node to basis connectivity)
                        tVertEnrichMap[tBasisIndex] = tNodeToBasisCount(tVertInd);
                        tVertexToBasis(tVertInd).push_back(tBasisIndex);
                        tNodeToBasisEnrLvl(tVertInd).push_back(aElementToBasisEnrichmentLevel(iEl)(iB));
                        tNodeToBasisCount(tVertInd) ++;

                    }
                    // verify the enrichment levels match
                    else
                    {
//                        MORIS_ASSERT(aElementToBasisEnrichmentLevel(iEl)(iB) == tNodeToBasisEnrLvl(tVertInd)(tIter->second),
//                                     "Enrichment levels for the same basis do not coincide: New EL = %u, EL in map = %u",
//                                     aElementToBasisEnrichmentLevel(iEl)(iB),
//                                     tNodeToBasisEnrLvl(tVertInd)(tIter->second));
                    }
                }
            }
        }
    }

    // iterate over vertex and add to Vertex enrichment data
    for(moris::uint i = 0 ; i <mVertexEnrichments.size(); i++)
    {
        // add basis information
        mVertexEnrichments(i).add_basis_information(tVertexToBasis(i),
                                                    tNodeToBasisEnrLvl(i));
    }

}

void
Enrichment::compute_vertex_basis_weights()
{
    //TODO: either uses INT or Create these in a factory
    // using hexa 8 as the basis function
    Hexahedron_8_Basis_Function tHex8Basis;
    moris::Matrix< moris::DDRMat > tBasisWeights(1,8);

    // Get underlying matrix data to access function
    moris::mtk::Mesh & tXTKMeshData = mBackgroundMeshPtr->get_mesh_data();

    // Iterate through children meshes and compute basis weights
    moris::uint tNumChildMeshes = mCutMeshPtr->get_num_child_meshes();

    for(moris::uint iCM = 0; iCM < tNumChildMeshes; iCM++)
    {
       // reference to child mesh
       Child_Mesh & tChildMesh = mCutMeshPtr->get_child_mesh(iCM);

       // Parent element index
       moris::moris_index tParentIndex = tChildMesh.get_parent_element_index();

       // Nodes attached to parent element
       moris::Matrix< moris::IndexMat > tNodesAttachedToParentElem = tXTKMeshData.get_entity_connected_to_entity_loc_inds(tParentIndex,moris::EntityRank::ELEMENT,moris::EntityRank::NODE);

       // get the node indices in this child mesh
       moris::Matrix<moris::IndexMat> const & tNodeIndicesOfCM = tChildMesh.get_node_indices();

       // reference to the parametric coordinates
       moris::Matrix< moris::DDRMat > const & tParamCoords = tChildMesh.get_parametric_coordinates();

       for(moris::uint iN = 0; iN < tChildMesh.get_num_entities(EntityRank::NODE); iN++)
       {
           tHex8Basis.evaluate_basis_function(tParamCoords.get_row(iN),tBasisWeights);

           MORIS_ASSERT(moris::equal_to(moris::sum(tBasisWeights),1.0),"Basis weights do not sum to 1");
           mVertexEnrichments(tNodeIndicesOfCM(iN)).add_basis_weights(tNodesAttachedToParentElem,tBasisWeights);

       }
    }

    // iterate over all vertex enrichments and condense out 0s
    for(moris::uint i = 0; i<mVertexEnrichments.size(); i++)
    {
        mVertexEnrichments(i).condense_out_basis_with_0_weight();
    }

    // perform multiplication of tmatrices is the basis is not a node
    if(mParameters.mBasisToEnrich !=EntityRank::NODE)
    {
        combine_t_matrices();
    }



}

void
Enrichment::combine_t_matrices()
{
    // iterate through background nodes and construct
//    moris::uint tNumBGNodes = mBackgroundMeshPtr->get_num_entities_background(EntityRank::NODE);
//    moris::uint tNumNodes   = mBackgroundMeshPtr->get_num_entities(EntityRank::NODE);
//
//    for(moris::uint iN = 0; iN<tNumBGNodes; iN++)
//    {
//        moris::Matrix< moris::DDRMat > & tNodeWeight = mVertexEnrichments(iN).get_basis_weights();
//
//        mBackgroundMeshPtr->
//
//    }

}

void
Enrichment::assign_enrichment_level_identifiers()
{

    mBasisEnrichmentIndices.resize(mElementIndsInBasis.size());

    moris::uint mNumEnrichmentLevels = 0;
    for(moris::uint i = 0; i <mElementIndsInBasis.size(); i++)
    {
        moris::moris_index tMaxEnrLev = mElementEnrichmentLevel(i).max() + 1;
        mNumEnrichmentLevels = mNumEnrichmentLevels + tMaxEnrLev;
        mBasisEnrichmentIndices(i) = moris::Matrix<moris::IndexMat>(tMaxEnrLev,1);

    }
    //TODO: Parallel strategy (change this to basis)
//    moris::moris_id    tIDOffset = mBackgroundMeshPtr->allocate_entity_ids(mNumEnrichmentLevels,EntityRank::NODE);
    moris::moris_index tIndOffset = mBackgroundMeshPtr->get_first_available_index(EntityRank::ELEMENT);

    for(moris::uint  i = 0; i < mBasisEnrichmentIndices.size(); i++)
    {
        moris::Matrix<moris::IndexMat> &  tBasisEnrichmentInds = mBasisEnrichmentIndices(i);
        tBasisEnrichmentInds(0) = i;
        for(moris::uint j = 1 ; j < tBasisEnrichmentInds.numel(); j++)
        {
            tBasisEnrichmentInds(j) = tIndOffset;
            tIndOffset++;
        }

    }

    mBackgroundMeshPtr->update_first_available_index(tIndOffset,EntityRank::ELEMENT);

}



}


