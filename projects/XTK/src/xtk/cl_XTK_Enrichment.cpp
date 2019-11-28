/*
 * cl_XTK_Enrichment.cpp
 *
 *  Created on: Feb 18, 2019
 *      Author: doble
 */

#include "cl_XTK_Enrichment.hpp"
#include "cl_XTK_Cut_Mesh.hpp"
#include "cl_XTK_Hexahedron_8_Basis_Function.hpp"
#include "fn_sum.hpp"
#include "fn_equal_to.hpp"
#include "fn_iscol.hpp"
#include "fn_trans.hpp"
#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"
//#include "cl_HMR_Database.hpp"
namespace xtk
{

Enrichment::Enrichment(enum Enrichment_Method aMethod,
                       enum EntityRank        aBasisRank,
                       moris::moris_index     aInterpIndex,
                       moris::moris_index     aNumBulkPhases,
                       xtk::Model*            aXTKModelPtr,
                       xtk::Cut_Mesh*         aCutMeshPtr,
                       xtk::Background_Mesh*  aBackgroundMeshPtr):
        mEnrichmentMethod(aMethod),
        mBasisRank(aBasisRank),
        mInterpIndex(aInterpIndex),
        mNumBulkPhases(aNumBulkPhases),
        mXTKModelPtr(aXTKModelPtr),
        mCutMeshPtr(aCutMeshPtr),
        mBackgroundMeshPtr(aBackgroundMeshPtr),
        mInterpCellBasis(aBackgroundMeshPtr->get_mesh_data().get_num_elems()),
        mInterpCellBasisEnrLev(aBackgroundMeshPtr->get_mesh_data().get_num_elems())
{

}



void
Enrichment::perform_enrichment()
{
    // Verify initialized properly
    MORIS_ERROR(mCutMeshPtr!=nullptr,"mCutMesh nullptr detected, this is probably because the enrichment has not been initialized properly");
    MORIS_ERROR(mBackgroundMeshPtr!=nullptr,"mBackgroundMesh nullptr detected, this is probably because the enrichment has not been initialized properly");

    // Perform enrichment over basis clusters
    perform_basis_cluster_enrichment();

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
Enrichment::perform_basis_cluster_enrichment()
{
    // Get underlying matrix data to access function
    moris::mtk::Interpolation_Mesh & tXTKMeshData = mBackgroundMeshPtr->get_mesh_data();

    // Number of basis functions
    moris::size_t tNumBasis = tXTKMeshData.get_num_basis_functions();

    MORIS_ASSERT(tXTKMeshData.get_num_elems() > 0,"0 cell interpolation mesh passed");

    // construct cell in xtk conformal model neighborhood connectivity
    this->construct_neighborhoods();

    // Allocate member variables
    mElementEnrichmentLevel = moris::Cell<moris::Matrix< moris::IndexMat >>(tNumBasis);
    mElementIndsInBasis     = moris::Cell<moris::Matrix< moris::IndexMat >>(tNumBasis);

    for(moris::size_t i = 0; i<tNumBasis; i++)
    {
        // Get elements in support of basis (these are interpolation cells)
        moris::Matrix< moris::IndexMat > tParentElementsInSupport;
        tXTKMeshData.get_elements_in_support_of_basis(mInterpIndex,i,tParentElementsInSupport);

        // get subphase clusters in support (seperated by phase)
        Matrix<moris::IndexMat> tSubphaseClusterIndicesInSupport = this->get_subphase_clusters_in_support(tParentElementsInSupport);

        // construct subphase in support map
        IndexMap  tSubPhaseIndexToSupportIndex;
        this->construct_subphase_in_support_map(tSubphaseClusterIndicesInSupport,tSubPhaseIndexToSupportIndex);

        // prune the subphase to remove subphases outside of basis support
        moris::Matrix< moris::IndexMat > tPrunedSubphaseNeighborhood;
        this->generate_pruned_subphase_graph_in_basis_support(tSubphaseClusterIndicesInSupport, tSubPhaseIndexToSupportIndex, tPrunedSubphaseNeighborhood);

        // Assign enrichment levels to subphases
        moris::Matrix< moris::IndexMat > tSubPhaseBinEnrichment;
        this->assign_subphase_bin_enrichment_levels_in_basis_support(tSubphaseClusterIndicesInSupport, tSubPhaseIndexToSupportIndex, tPrunedSubphaseNeighborhood,tSubPhaseBinEnrichment);

        // Extract element enrichment levels from assigned sub-phase bin enrichment levels and store these as a member variable
        this->unzip_subphase_bin_enrichment_into_element_enrichment(i, tParentElementsInSupport, tSubphaseClusterIndicesInSupport, tSubPhaseIndexToSupportIndex, tPrunedSubphaseNeighborhood,tSubPhaseBinEnrichment);
    }

    // assign enrichment level indices
    assign_enrichment_level_identifiers();

    // set element to basis connectivity
    this->construct_element_to_basis_connectivity(mElementToBasis,mElementToBasisEnrichmentLevel);

    this->construct_enriched_interpolation_mesh();

    this->construct_enriched_integration_mesh();
}

void
Enrichment::construct_neighborhoods()
{
    // construct full mesh neighbhoord
    mXTKModelPtr->construct_neighborhood();

    // construct subphase neighborhood
    mXTKModelPtr->construct_subphase_neighborhood();
}


Matrix<IndexMat>
Enrichment::get_subphase_clusters_in_support(moris::Matrix< moris::IndexMat > const & aElementsInSupport)
{
    // count the number of subphase cluster in support
    moris::uint tCount = 0;
    for(moris::size_t iE = 0; iE<aElementsInSupport.numel(); iE++)
    {
        if(mBackgroundMeshPtr->entity_has_children(aElementsInSupport(iE),EntityRank::ELEMENT))
        {
            moris::moris_index tCMIndex = mBackgroundMeshPtr->child_mesh_index(aElementsInSupport(iE),EntityRank::ELEMENT);
            Child_Mesh const * tCM      = & mCutMeshPtr->get_child_mesh(tCMIndex);
            tCount = tCount + tCM->get_num_subphase_bins();
        }

        else
        {
            tCount++;
        }
    }
    Matrix<IndexMat> tSubPhaseClusters(1,tCount);
    tCount = 0;
    for(moris::size_t iE = 0; iE<aElementsInSupport.numel(); iE++)
    {
        if(mBackgroundMeshPtr->entity_has_children(aElementsInSupport(iE),EntityRank::ELEMENT))
        {
            moris::moris_index tCMIndex = mBackgroundMeshPtr->child_mesh_index(aElementsInSupport(iE),EntityRank::ELEMENT);
            Child_Mesh const * tCM      = & mCutMeshPtr->get_child_mesh(tCMIndex);
            Cell<moris_index> const & tCMSubPhaseIndices = tCM->get_subphase_indices();

            // add subphase indices to output matrix
            for(moris::uint iSP = 0; iSP < tCMSubPhaseIndices.size(); iSP++)
            {
                tSubPhaseClusters(tCount) = tCMSubPhaseIndices(iSP);
                tCount++;
            }
        }

        else
        {
            tSubPhaseClusters(tCount) = aElementsInSupport(iE);
            tCount++;
        }
    }
    return tSubPhaseClusters;
}

void
Enrichment::construct_subphase_in_support_map(moris::Matrix< moris::IndexMat > const & aSubphaseClusterIndicesInSupport,
                                              IndexMap & aSubPhaseIndexToSupportIndex)
{
    for(moris::moris_index i = 0; i <(moris::moris_index) aSubphaseClusterIndicesInSupport.numel(); i++)
    {
        aSubPhaseIndexToSupportIndex[aSubphaseClusterIndicesInSupport(i)] = i;
    }
}


void
Enrichment::generate_pruned_subphase_graph_in_basis_support(moris::Matrix< moris::IndexMat > const & aSubphasesInSupport,
                                                            IndexMap &                               aSubPhaseIndexToSupportIndex,
                                                            moris::Matrix< moris::IndexMat >       & aPrunedSubPhaseToSubphase)
 {
    // Construct full element neighbor graph in support and the corresponding shared faces
    aPrunedSubPhaseToSubphase.resize(aSubphasesInSupport.numel(), 50);
    aPrunedSubPhaseToSubphase.fill(MORIS_INDEX_MAX);

    // get subphase neighborhood information
    moris::Cell<moris::Cell<moris_index>>  const & tSubphasetoSubphase = mXTKModelPtr->get_subphase_to_subphase();

    for(moris::size_t i = 0; i<aSubphasesInSupport.numel(); i++)
    {
        moris::Cell<moris_index> const & tSingleSubPhaseNeighbors = tSubphasetoSubphase(aSubphasesInSupport(i));

        // iterate through and prune subphases not in support
        moris::uint tCount = 0;
        for(moris::size_t j = 0; j < tSingleSubPhaseNeighbors.size(); j++)
        {
            moris_index tNeighborSubphaseIndex = tSingleSubPhaseNeighbors(j);

            auto tNeighborIter = aSubPhaseIndexToSupportIndex.find(tNeighborSubphaseIndex) ;

            if(tNeighborIter != aSubPhaseIndexToSupportIndex.end())
            {
                aPrunedSubPhaseToSubphase(i,tCount) = tNeighborIter->second;
                tCount++;
            }
        }
    }

 }


void
Enrichment::assign_subphase_bin_enrichment_levels_in_basis_support(moris::Matrix< moris::IndexMat > const & aSubphasesInSupport,
                                                                   IndexMap &                               aSubPhaseIndexToSupportIndex,
                                                                   moris::Matrix< moris::IndexMat > const & aPrunedSubPhaseToSubphase,
                                                                   moris::Matrix< moris::IndexMat >       & aSubPhaseBinEnrichmentVals)
{
    // Variables needed for floodfill, consider removing these.
    // Active bins to include in floodfill (We include all bins)
    moris::Matrix< moris::IndexMat > tActiveBins(1,aPrunedSubPhaseToSubphase.n_rows());
    for(moris::size_t i = 0; i< aPrunedSubPhaseToSubphase.n_rows(); i++)
    {
        (tActiveBins)(0,i) = i;
    }

    // Mark all as included
    moris::Matrix< moris::IndexMat > tIncludedBins(1,aSubphasesInSupport.numel(),1);

    // Flood fill metric value (since all the subphases do not connect to dissimilar phases)
    moris::Matrix< moris::IndexMat > tDummyPhase(1,aSubphasesInSupport.numel(),1);

    aSubPhaseBinEnrichmentVals = flood_fill(aPrunedSubPhaseToSubphase,
                                            tDummyPhase,
                                            tActiveBins,
                                            tIncludedBins,
                                            mNumBulkPhases,
                                            INDEX_MAX,
                                            true);

}
void
Enrichment::unzip_subphase_bin_enrichment_into_element_enrichment(moris_index aBasisIndex,
                                                                  moris::Matrix< moris::IndexMat > const & aParentElementsInSupport,
                                                                  moris::Matrix< moris::IndexMat > const & aSubphasesInSupport,
                                                                  IndexMap &                               aSubPhaseIndexToSupportIndex,
                                                                  moris::Matrix< moris::IndexMat > const & aPrunedSubPhaseToSubphase,
                                                                  moris::Matrix< moris::IndexMat >       & aSubPhaseBinEnrichmentVals)
{

    // access the subphase index to child mesh index connectivity in cut mesh
    Matrix<IndexMat> const & tSubphaseToCM = mCutMeshPtr->get_subphase_to_child_mesh_connectivity();

    // resize member data
    moris::size_t tNumAllElementsInSupport = count_elements_in_support(aParentElementsInSupport);
    mElementIndsInBasis(aBasisIndex)       = moris::Matrix< moris::IndexMat >(1,tNumAllElementsInSupport);
    mElementEnrichmentLevel(aBasisIndex)   = moris::Matrix< moris::IndexMat >(1,tNumAllElementsInSupport);

    moris::uint tCount = 0;

    for(moris::size_t i = 0; i<aSubphasesInSupport.numel(); i++)
    {
        bool tSubphaseInChildMesh = mXTKModelPtr->subphase_is_in_child_mesh(aSubphasesInSupport(i));

        moris_index tSubphaseIndex = aSubphasesInSupport(i);

        if(tSubphaseInChildMesh)
        {
            // get the child mesh
            Child_Mesh & tChildMesh = mCutMeshPtr->get_child_mesh(tSubphaseToCM(tSubphaseIndex));

            // get subphase child mesh
            tChildMesh.add_basis_and_enrichment_to_subphase_group(tSubphaseIndex,aBasisIndex,aSubPhaseBinEnrichmentVals(i));

            // get the child elements in subphas
            Cell<moris::Matrix< moris::IndexMat >> const & tCellToSubphase = tChildMesh.get_subphase_groups();

            // element indices of child mesh
            Matrix<IndexMat> const & tCellInds = tChildMesh.get_element_inds();

            // get the subphase index local to child mesh
            moris_index tCMLocalSubphaseIndex = tChildMesh.get_subphase_loc_index(tSubphaseIndex);

            for(moris::uint iCE = 0; iCE < tCellToSubphase(tCMLocalSubphaseIndex).numel(); iCE++)
            {
                mElementIndsInBasis(aBasisIndex)(tCount)     = tCellInds(tCellToSubphase(tCMLocalSubphaseIndex)(iCE));
                mElementEnrichmentLevel(aBasisIndex)(tCount) = aSubPhaseBinEnrichmentVals(i);
                tCount++;
            }
        }
        else
        {
            mElementIndsInBasis(aBasisIndex)(tCount)     = aSubphasesInSupport(i);
            mElementEnrichmentLevel(aBasisIndex)(tCount) = aSubPhaseBinEnrichmentVals(i);
            tCount++;


            // add information to interp cells about which basis/enrichment level interpolates in it
            mInterpCellBasis( aSubphasesInSupport(i) ).push_back(aBasisIndex);
            mInterpCellBasisEnrLev(aSubphasesInSupport(i)).push_back(aSubPhaseBinEnrichmentVals(i));
        }
    }
}

void
Enrichment::assign_enrichment_level_identifiers()
{

    mBasisEnrichmentIndices.resize(mElementIndsInBasis.size());

     mNumEnrichmentLevels = 0;
    //moris::uint mNumEnrichmentLevels = 0;
    for(moris::uint i = 0; i <mElementIndsInBasis.size(); i++)
    {
        moris::moris_index tMaxEnrLev = mElementEnrichmentLevel(i).max() + 1;
        mNumEnrichmentLevels = mNumEnrichmentLevels + tMaxEnrLev;
        mBasisEnrichmentIndices(i) = moris::Matrix<moris::IndexMat>(tMaxEnrLev,1);
    }

    //TODO: Parallel strategy (change this to basis)
//    moris::moris_id    tIDOffset = mBackgroundMeshPtr->allocate_entity_ids(mNumEnrichmentLevels,EntityRank::NODE);
    moris::moris_index tIndOffset = mBackgroundMeshPtr->get_num_entities_background(mBasisRank);


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
    //mBackgroundMeshPtr->update_first_available_index(tIndOffset,EntityRank::ELEMENT);
}


moris::size_t
Enrichment::count_elements_in_support(moris::Matrix< moris::IndexMat > const & aParentElementsInSupport)
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
    aElementToBasis = moris::Cell<moris::Cell<moris::moris_index>>(mXTKModelPtr->get_num_elements_total());
    aElementToBasisEnrichmentLevel = moris::Cell<moris::Cell<moris::moris_index>>(mXTKModelPtr->get_num_elements_total());

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
Enrichment::print_basis_support_debug(moris_index aBasisIndex,
                          moris::Matrix< moris::IndexMat > const & aParentElementsInSupport,
                          moris::Matrix< moris::IndexMat > const & aSubphasesInSupport,
                          IndexMap &                               aSubPhaseIndexToSupportIndex,
                          moris::Matrix< moris::IndexMat > const & aPrunedSubPhaseToSubphase,
                          moris::Matrix< moris::IndexMat >       & aSubPhaseBinEnrichmentVals)
{
    std::cout<<"--------------------------------------------------"<<std::endl;
    std::cout<<"Basis Index: "<<aBasisIndex<<std::endl;
    std::cout<<"Parent Cell In Support:";
    for(moris::uint i = 0; i < aParentElementsInSupport.numel(); i++)
    {
        std::cout<<std::setw(8)<<aParentElementsInSupport(i);
    }
    std::cout<<"\nSubphases In Support:";
    for(moris::uint i = 0; i < aSubphasesInSupport.numel(); i++)
    {
        std::cout<<std::setw(8)<<aSubphasesInSupport(i);
    }

    std::cout<<"\nSubphase Neighborhood In Support:"<<std::endl;
    for(moris::uint i = 0; i<aPrunedSubPhaseToSubphase.n_rows(); i++ )
    {
        std::cout<<std::setw(6)<<aSubphasesInSupport(i)<<" | ";


        for(moris::uint j = 0; j< aPrunedSubPhaseToSubphase.n_cols(); j++)
        {
            if(aPrunedSubPhaseToSubphase(i,j) != MORIS_INDEX_MAX)
            {
                std::cout<<std::setw(6)<<aSubphasesInSupport(aPrunedSubPhaseToSubphase(i,j));
            }
        }
        std::cout<<std::endl;
    }

    std::cout<<"\nSubphase Enrichment Level: \n";
    for(moris::uint i = 0; i < aSubPhaseBinEnrichmentVals.numel(); i++)
    {
        std::cout<<std::setw(8)<<aSubPhaseBinEnrichmentVals(i);
    }
    std::cout<<"\n--------------------------------------------------"<<std::endl;
}



void
Enrichment::construct_enriched_interpolation_mesh()
{
    // initialize a new enriched interpolation mesh
    mXTKModelPtr->mEnrichedInterpMesh(mInterpIndex) = new Enriched_Interpolation_Mesh(mXTKModelPtr);

    // allocate space for interpolation cells
    this->allocate_interpolation_cells();

    this->construct_enriched_interpolation_vertices_and_cells();

    // add the coeff to enriched coefs to enriched interpolation mesh
    mXTKModelPtr->mEnrichedInterpMesh(mInterpIndex)->mCoeffToEnrichCoeffs = mBasisEnrichmentIndices;

    mXTKModelPtr->mEnrichedInterpMesh(mInterpIndex)->finalize_setup();
}

void
Enrichment::construct_enriched_integration_mesh()
{
    MORIS_ASSERT(mXTKModelPtr->mEnrichedInterpMesh(mInterpIndex) != nullptr,"No enriched interpolation mesh to link enriched integration mesh to");

    mXTKModelPtr->mEnrichedIntegMesh(mInterpIndex) = new Enriched_Integration_Mesh(mXTKModelPtr,mInterpIndex);
}

void
Enrichment::allocate_interpolation_cells()
{
    // enriched interpolation mesh pointer
    Enriched_Interpolation_Mesh* tEnrInterpMesh = mXTKModelPtr->mEnrichedInterpMesh(mInterpIndex);

    // iterate through child meshes
    uint tNumChildMeshes = mCutMeshPtr->get_num_child_meshes();

    // count how many subphases there are total in the entire mesh
    uint tNumSubphases = 0;
    for(moris::uint iCM = 0;  iCM <tNumChildMeshes; iCM++)
    {
        Child_Mesh const & tCM = mCutMeshPtr->get_child_mesh(iCM);

        tNumSubphases = tNumSubphases + tCM.get_num_subphase_bins();
    }

    // there is one enrichment cell per subphase in children meshes and 1 for each unintersected elements
    moris::uint tNumEnrInterpCells = tNumSubphases + mBackgroundMeshPtr->get_num_entities(EntityRank::ELEMENT) - tNumChildMeshes;

    tEnrInterpMesh->mEnrichedInterpCells.resize(tNumEnrInterpCells);

    // assuming all interpolation cells are the same
    // figure out how many vertices there are per interpolation cell
    uint tNumVertsPerCell = 0;
    if(mBackgroundMeshPtr->get_mesh_data().get_num_elems()> 0)
    {
        tNumVertsPerCell = mBackgroundMeshPtr->get_mesh_data().get_mtk_cell(0).get_number_of_vertices();
    }

    tEnrInterpMesh->mNumVertsPerInterpCell = tNumVertsPerCell;

    // allocate maximum number of enriched vertices
    tEnrInterpMesh->mEnrichedInterpVerts.resize(tNumVertsPerCell*tNumEnrInterpCells);

    // allocate the base vertices to vertex enrichment data
    tEnrInterpMesh->mBaseInterpVertToVertEnrichmentIndex.resize(mBackgroundMeshPtr->get_mesh_data().get_num_nodes());

    // allocate base cell to enriched cell data
    tEnrInterpMesh->mBaseCelltoEnrichedCell.resize(mBackgroundMeshPtr->get_mesh_data().get_num_elems());
}

void
Enrichment::construct_enriched_interpolation_vertices_and_cells()
{
    // allocate indices and ids
    moris_index tCellIndex = 0;

    // fixme: allocate ids to each processor
    moris_id    tCellId    = 1;

    // allocate vertex indices and ids
    moris_index tVertId    = 1;
    uint        tVertexCount = 0;

    // iterate through children meshes and create the unzipped interpolation cells
    uint tNumChildMeshes = mCutMeshPtr->get_num_child_meshes();

    // not enriched interpolation mesh data
    moris::mtk::Mesh & tMesh = mBackgroundMeshPtr->get_mesh_data();

    // enriched interpolation mesh pointer
    Enriched_Interpolation_Mesh* tEnrInterpMesh = mXTKModelPtr->mEnrichedInterpMesh(mInterpIndex);

    // Enriched Interpolation Cell to Vertex Index
    Matrix<IndexMat> tEnrInterpCellToVertex(tEnrInterpMesh->get_num_elements(),tEnrInterpMesh->mNumVertsPerInterpCell);

    // construct a connectivity for the enriched interpolation cells
    if(tMesh.get_num_elems() > 0)
    {
        mtk::Cell const & tDummyCell = tMesh.get_mtk_cell(0);
        mtk::Cell_Info_Factory tFactory;
        tEnrInterpMesh->mCellInfo = tFactory.create_cell_info(tDummyCell.get_geometry_type(),tDummyCell.get_interpolation_order());
    }


    // construct enriched cells/vertices for the unintersected background cells
    for(moris::uint iEl = 0; iEl < tMesh.get_num_elems(); iEl++ )
    {
        if(!mBackgroundMeshPtr->entity_has_children((moris_index)iEl,EntityRank::ELEMENT))
        {
            // information about this cell
            moris::mtk::Cell & tParentCell = tMesh.get_mtk_cell((moris_index)iEl);

            // owner
            moris_id tOwner = tParentCell.get_owner();

            // vertexes of cell
            moris::Cell<mtk::Vertex*> tVertices = tParentCell.get_vertex_pointers();

            // vertex interpolations of the parent cell
            moris::Cell<mtk::Vertex_Interpolation*> tVertexInterpolations = tParentCell.get_vertex_interpolations( mInterpIndex );

            // number of vertices
            uint tNumVertices = tParentCell.get_number_of_vertices();

            // basis interpolating into the cell
            Cell<moris_index> const & tBasisInCell = mInterpCellBasis((moris_index)iEl);

            // bulk phase
            moris_index tBulkPhase = mXTKModelPtr->mBackgroundMesh.get_element_phase_index(tParentCell.get_index());

            // construct a map between basis index and index relative to the subphase cluster
            std::unordered_map<moris_id,moris_id> tCellBasisMap = construct_subphase_basis_to_basis_map(tBasisInCell);

            // Subphase basis enrichment level of the current subphase in child mesh
            Cell<moris_index> const & tEnrLevOfBasis = mInterpCellBasisEnrLev((moris_index)iEl);

            // construct unzipped enriched vertices
            for(uint iEV = 0; iEV < tNumVertices; iEV++)
            {
                // construct vertex enrichment
                Vertex_Enrichment tVertEnrichment;
                this->construct_enriched_vertex_interpolation(tVertexInterpolations(iEV),tEnrLevOfBasis,tCellBasisMap,tVertEnrichment);

                // add vertex enrichment to enriched interpolation mesh
                bool tNewVertFlag = false;
                moris_index tVertEnrichIndex = tEnrInterpMesh->add_vertex_enrichment(tVertices(iEV),tVertEnrichment, tNewVertFlag);

                // if this vertex interpolation is new, create a vertex
                if(tNewVertFlag)
                {
                    // Create interpolation vertex
                    tEnrInterpMesh->mEnrichedInterpVerts(tVertEnrichIndex)
                                = Interpolation_Vertex_Unzipped(tVertices(iEV), tVertId, tVertEnrichIndex,
                                                                tVertices(iEV)->get_owner(), mInterpIndex,
                                                                tEnrInterpMesh->get_vertex_enrichment(tVertEnrichIndex));

                    tVertId++;
                    tVertexCount++;
                }

                tEnrInterpCellToVertex(tCellIndex,iEV) = tVertEnrichIndex;
            }

            // create new enriched interpolation cell
            tEnrInterpMesh->mEnrichedInterpCells(tCellIndex) = Interpolation_Cell_Unzipped(&tParentCell, tParentCell.get_index(),
                                                                                           tBulkPhase, tCellId, tCellIndex, tOwner, tEnrInterpMesh->mCellInfo);

            // add enriched interpolation cell to base cell to enriched cell data
            tEnrInterpMesh->mBaseCelltoEnrichedCell(tParentCell.get_index()).push_back(&tEnrInterpMesh->mEnrichedInterpCells(tCellIndex));

            // increment the cell index/id
            tCellIndex++;
            tCellId++;

        }
    }

    // iterate through child meshes and construct the unzipped vertices if necessary and the cells
    for(moris::uint iCM = 0;  iCM <tNumChildMeshes; iCM++)
    {
        // get the child mesh
        Child_Mesh & tCM = mCutMeshPtr->get_child_mesh(iCM);

        // get the parent cell of this child mesh
        moris_index tParentCellIndex = tCM.get_parent_element_index();
        moris::mtk::Cell & tParentCell = tMesh.get_mtk_cell(tParentCellIndex);

        // owner
        moris_id tOwner = tParentCell.get_owner();

        // vertexes of cell
        moris::Cell<mtk::Vertex*> tVertices = tParentCell.get_vertex_pointers();

        // vertex interpolations of the parent cell
        moris::Cell<mtk::Vertex_Interpolation*> tVertexInterpolations = tParentCell.get_vertex_interpolations( mInterpIndex );

        // number of vertices
        uint tNumVertices = tParentCell.get_number_of_vertices();

        // get bulkphase of subphases
        Cell<moris::moris_index> const & tSubPhaseBulkPhase =tCM.get_subphase_bin_bulk_phase();

        //
        Cell<moris_index> const & tSubphaseIndices = tCM.get_subphase_indices();

        // create a interpolation cell per subphase
        for(moris::uint iSP = 0; iSP<tCM.get_num_subphase_bins(); iSP++)
        {
            // Subphase basis of the current subphase in child mesh
            Cell<moris_index> const & tSubPhaseBasis = tCM.get_subphase_basis_indices((moris_index)iSP);

            // construct a map between basis index and index relative to the subphase cluster
            std::unordered_map<moris_id,moris_id> tSubPhaseBasisMap = construct_subphase_basis_to_basis_map(tSubPhaseBasis);

            // Subphase basis enrichment level of the current subphase in child mesh
            Cell<moris_index> const & tSubPhaseBasisEnrLev = tCM.get_subphase_basis_enrichment_levels((moris_index)iSP);

            // construct unzipped enriched vertices
            for(uint iEV = 0; iEV < tNumVertices; iEV++)
            {
                // construct vertex enrichment
                Vertex_Enrichment tVertEnrichment;
                this->construct_enriched_vertex_interpolation(tVertexInterpolations(iEV),tSubPhaseBasisEnrLev,tSubPhaseBasisMap,tVertEnrichment);

                // add vertex enrichment to enriched interpolation mesh
                bool tNewVertFlag = false;
                moris_index tVertEnrichIndex = tEnrInterpMesh->add_vertex_enrichment(tVertices(iEV),tVertEnrichment, tNewVertFlag);

                // if this vertex interpolation is new, create a vertex
                if(tNewVertFlag)
                {
                    // Create interpolation vertex
                    tEnrInterpMesh->mEnrichedInterpVerts(tVertEnrichIndex) = Interpolation_Vertex_Unzipped(tVertices(iEV),tVertId,tVertEnrichIndex,tVertices(iEV)->get_owner(),mInterpIndex,tEnrInterpMesh->get_vertex_enrichment(tVertEnrichIndex));

                    tVertId ++;
                    tVertexCount++;
                }

                tEnrInterpCellToVertex(tCellIndex,iEV) = tVertEnrichIndex;
            }

            // create new enriched interpolation cell
            tEnrInterpMesh->mEnrichedInterpCells(tCellIndex) = Interpolation_Cell_Unzipped(&tParentCell, tSubphaseIndices(iSP), tSubPhaseBulkPhase(iSP), tCellId, tCellIndex, tOwner, tEnrInterpMesh->mCellInfo );

            // add enriched interpolation cell to base cell to enriched cell data
            tEnrInterpMesh->mBaseCelltoEnrichedCell(tParentCell.get_index()).push_back(&tEnrInterpMesh->mEnrichedInterpCells(tCellIndex));

            // increment the cell index/id
            tCellIndex++;
            tCellId++;
        }
    }

    // with the cell to vertex data fully setup, add the vertex pointers to the cell

    for(moris::uint iE = 0; iE< tEnrInterpCellToVertex.n_rows(); iE++)
    {
        moris::Cell<Interpolation_Vertex_Unzipped*> tVertices(tEnrInterpCellToVertex.n_cols());
        // iterate and get vertices
        for(moris::uint iV =0; iV<tEnrInterpCellToVertex.n_cols(); iV++)
        {

         tVertices(iV) = tEnrInterpMesh->get_unzipped_vertex_pointer(tEnrInterpCellToVertex(iE,iV));
        }

        // set vertices in cell
        tEnrInterpMesh->mEnrichedInterpCells(iE).set_vertices(tVertices);
    }

    tEnrInterpMesh->mEnrichedInterpVerts.resize(tVertexCount);

}


void
Enrichment::construct_enriched_vertex_interpolation(mtk::Vertex_Interpolation*              aBaseVertexInterp,
                                                    Cell<moris_index> const &               aSubPhaseBasisEnrLev,
                                                    std::unordered_map<moris_id,moris_id> & aMapBasisIndexToLocInSubPhase,
                                                    Vertex_Enrichment &                     aVertexEnrichment)
{
    // allocate a new vertex enrichment
    aVertexEnrichment = Vertex_Enrichment();

    // iterate through basis in the base vertex interpolation
    moris::uint tNumCoeffs = aBaseVertexInterp->get_number_of_coefficients();

    // indices of the coefficients
    Matrix< IndexMat > tBaseVertCoeffInds = aBaseVertexInterp->get_indices();

    // weights of the coefficients
    const Matrix< DDRMat > * tBaseVertWeights = aBaseVertexInterp->get_weights();

    // enriched Basis Coefficient indices
    Matrix< IndexMat > tEnrichCoeffInds(tBaseVertCoeffInds.n_rows(),tBaseVertCoeffInds.n_cols());

    for(moris::uint iBC = 0; iBC<tNumCoeffs; iBC++)
    {
        // coefficient of the base index
        moris_index tCoeffIndex = tBaseVertCoeffInds(iBC);

        // find the coefficients index within this subphase cluster
        auto tIter = aMapBasisIndexToLocInSubPhase.find(tCoeffIndex);
        MORIS_ASSERT(tIter != aMapBasisIndexToLocInSubPhase.end(),"Basis not found in vertex map");

        // The basis local index relative to the subphase
        moris_index tSubphaseBasisIndex = tIter->second;

        // enrichment level of the basis
        moris_index tEnrLev = aSubPhaseBasisEnrLev(tSubphaseBasisIndex);

        moris_index tEnrichedCoeffIndex = mBasisEnrichmentIndices(tCoeffIndex)(tEnrLev);

        tEnrichCoeffInds(iBC) = tEnrichedCoeffIndex;

    }

    std::unordered_map<moris::moris_index,moris::moris_index> & tVertEnrichMap = aVertexEnrichment.get_basis_map();

    for(moris::uint iB = 0; iB<tEnrichCoeffInds.numel(); iB++)
    {
        moris::moris_index tBasisIndex = tEnrichCoeffInds(iB);

        tVertEnrichMap[tBasisIndex] = iB;
    }

    aVertexEnrichment.add_basis_information(tEnrichCoeffInds);
    aVertexEnrichment.add_basis_weights(tEnrichCoeffInds,*tBaseVertWeights);
}

std::unordered_map<moris_id,moris_id>
Enrichment::construct_subphase_basis_to_basis_map(Cell<moris_id> const & aSubPhaseBasisIndex)
{
    uint tNumBasisOfSubphase = aSubPhaseBasisIndex.size();

    std::unordered_map<moris_id,moris_id> tSubphaseBasisMap;
    for(moris::uint iB =0; iB < tNumBasisOfSubphase; iB++)
    {
        tSubphaseBasisMap[aSubPhaseBasisIndex(iB)] = iB;
    }

    return tSubphaseBasisMap;
}



Cell<std::string>
Enrichment::get_cell_enrichment_field_names() const
{
    // background mesh data
    moris::mtk::Mesh & tXTKMeshData = mBackgroundMeshPtr->get_mesh_data();

    // number of basis
    moris::uint tNumBasis = tXTKMeshData.get_num_basis_functions();

    // declare  enrichment fields
    Cell<std::string> tEnrichmentFields(tNumBasis);
    std::string tBaseEnrich = "enrlev_basis_";
    for(size_t i = 0; i<tNumBasis; i++)
    {
        tEnrichmentFields(i) = tBaseEnrich + std::to_string(mBackgroundMeshPtr->get_glb_entity_id_from_entity_loc_index(i,EntityRank::NODE));
    }

    // Add local floodfill field to the output mesh
    std::string tLocalFFStr = "child_ff";
    tEnrichmentFields.push_back(tLocalFFStr);

    std::string tSubPhaseStr = "subphase";
    tEnrichmentFields.push_back(tSubPhaseStr);

    return tEnrichmentFields;
}

void
Enrichment::write_cell_enrichment_to_fields(Cell<std::string>  & aEnrichmentFieldStrs,
                                            mtk::Mesh*           aMeshWithEnrFields) const
{
    // background mesh data
    moris::mtk::Mesh & tXTKMeshData = mBackgroundMeshPtr->get_mesh_data();

    // number of basis
    moris::uint tNumBasis = tXTKMeshData.get_num_basis_functions();

    // Local subphase bins
    moris::Matrix<moris::DDRMat> tLocalSubphaseVal(aMeshWithEnrFields->get_num_entities(moris::EntityRank::ELEMENT),1);
    for(size_t i = 0; i<mCutMeshPtr->get_num_child_meshes(); i++)
    {

        Child_Mesh & tChildMesh = mCutMeshPtr->get_child_mesh(i);

        moris::Matrix< moris::IndexMat > const & tElementSubphases = tChildMesh.get_elemental_subphase_bin_membership();

        moris::Matrix< moris::IdMat > const & tChildCellIds = tChildMesh.get_element_ids();

        for(size_t j = 0; j<tChildCellIds.n_cols(); j++)
        {
            moris_index tNewMeshInd = aMeshWithEnrFields->get_loc_entity_ind_from_entity_glb_id(tChildCellIds(j),EntityRank::ELEMENT);
            tLocalSubphaseVal(tNewMeshInd) = (real)(tElementSubphases(0,j));
        }
    }
    std::string tLocalFFStr = "child_ff";
    aMeshWithEnrFields->add_mesh_field_real_scalar_data_loc_inds(tLocalFFStr, moris::EntityRank::ELEMENT, tLocalSubphaseVal);

    // subphase field
    moris::Matrix<moris::DDRMat> tSubPhase(aMeshWithEnrFields->get_num_entities(moris::EntityRank::ELEMENT),1);

    moris::Matrix<moris::IndexMat> tXTKSubphases = mXTKModelPtr->get_element_to_subphase();


    for(moris::uint i = 0; i < aMeshWithEnrFields->get_num_entities(EntityRank::ELEMENT); i++)
    {
        moris_id    tGlbId  = aMeshWithEnrFields->get_glb_entity_id_from_entity_loc_index((moris_index)i,EntityRank::ELEMENT);
        moris_index tXTKInd = mXTKModelPtr->get_cell_xtk_index(tGlbId);
        tSubPhase(i) = tXTKSubphases(tXTKInd);
    }
    std::string tSubPhaseStr = "subphase";
    aMeshWithEnrFields->add_mesh_field_real_scalar_data_loc_inds(tSubPhaseStr, moris::EntityRank::ELEMENT, tSubPhase);

    // Enrichment values
    Cell<moris::Matrix< moris::IndexMat >> const & tElementIndsInBasis = this->get_element_inds_in_basis_support();
    Cell<moris::Matrix< moris::IndexMat >> const & tElementEnrichmentInBasis = this->get_element_enrichment_levels_in_basis_support();


    for(size_t i = 0; i<tNumBasis; i++)
    {
        moris::Matrix<moris::DDRMat> tEnrichmentLevels(aMeshWithEnrFields->get_num_entities(moris::EntityRank::ELEMENT),1,10);

        for(size_t j = 0; j<tElementIndsInBasis(i).numel(); j++)
        {
            moris_index tXTKMeshInd = (tElementIndsInBasis(i))(j);
            moris_id    tXTKMeshId  = mBackgroundMeshPtr->get_glb_entity_id_from_entity_loc_index(tXTKMeshInd,EntityRank::ELEMENT);
            moris_index tNewMeshInd = aMeshWithEnrFields->get_loc_entity_ind_from_entity_glb_id(tXTKMeshId,EntityRank::ELEMENT);

            tEnrichmentLevels(tNewMeshInd) = (real)(((tElementEnrichmentInBasis(i)))(j));

        }

        aMeshWithEnrFields->add_mesh_field_real_scalar_data_loc_inds(aEnrichmentFieldStrs(i), moris::EntityRank::ELEMENT, tEnrichmentLevels);
    }
}



//
//void
//Enrichment::create_multilevel_enrichments()
//{
//    //-----------------------------------------------------------
//    //---------------------Hardcoded stuff-----------------------
//
//    //moris::sint tNumLevels = 3;
//
//    moris::sint tMeshIndex = 1;
//
//    //-----------------------------------------------------------
//
//
//    std::shared_ptr< moris::mtk::Mesh > tMeshPointer = mXTKModelPtr->mHMRMesh;
//
////    moris::uint tMaxMeshLevel = tMeshPointer->get_HMR_database()
////                                            ->get_bspline_mesh_by_index( tMeshIndex )
////                                            ->get_max_level();
//
//    moris::sint tMaxIndexForOrder = tMeshPointer->get_HMR_database()->get_bspline_mesh_by_index( tMeshIndex )
//                                                                    ->get_number_of_indexed_basis();
//
//    moris::uint tMaxNumBasis = mNumEnrichmentLevels * tMaxIndexForOrder;
//
//    moris::uint tOriginalBasisEnridgmentsize = mBasisEnrichmentIndices.size();
//
//    moris::uint tCounter = 1;
//
//    //print(mBasisEnrichmentIndices,"mBasisEnrichmentIndices");
//
//    mEnrichedMultilevelBasis       .set_size( tMaxNumBasis, 1, -1 );
//    mLevelOfEnrichedMultilevelBasis.set_size( tMaxNumBasis, 1, -1 );
//
//    moris::sint tMaxEnrichmentIndex = -1;
//    // Loop over
//    for ( moris::uint Ii = 0; Ii < mBasisEnrichmentIndices.size(); Ii++ )
//    {
//        moris::sint tMaxIndex = mBasisEnrichmentIndices( Ii ).max();
//
//        tMaxEnrichmentIndex = std::max( tMaxEnrichmentIndex, tMaxIndex );
//    }
//
//    mBasisEnrichmentIndices  .resize( tMaxIndexForOrder );
//    mBasisEnrichmentBulkPhase.resize( tMaxIndexForOrder );
//
//    for ( moris::uint Ii = 0; Ii < (uint)tMaxIndexForOrder; Ii++ )
//    {
//        moris::uint tBasisLevel = tMeshPointer->get_HMR_database()->get_bspline_mesh_by_index( tMeshIndex )
//                                                                  ->get_basis_by_index( Ii )
//                                                                  ->get_level();
//        if ( Ii < tOriginalBasisEnridgmentsize )
//        {
//            for ( moris::uint Ik = 0; Ik < mBasisEnrichmentIndices( Ii ).n_rows(); Ik++ )
//            {
//                moris::sint tEnrichedBasisIndex = mBasisEnrichmentIndices( Ii )( Ik, 0 );
//
//                MORIS_ASSERT( mLevelOfEnrichedMultilevelBasis( tEnrichedBasisIndex, 0 ) == -1,
//                        "Enrichment::create_multilevel_enrichments(), More than one HMR basis related to enriched basis %-5i", tEnrichedBasisIndex);
//
//                mLevelOfEnrichedMultilevelBasis( tEnrichedBasisIndex, 0 ) = tBasisLevel;
//            }
//        }
//        else
//        {
//            moris::Matrix< DDSMat > tIndices = tMeshPointer->get_HMR_database()
//                                                           ->get_bspline_mesh_by_index( tMeshIndex )
//                                                           ->get_children_ind_for_basis( Ii );
//
//            moris::uint tNumMaxEnrichLevel = 0;
//
//            for ( moris::uint Ik = 0; Ik < tIndices.n_rows(); Ik++ )
//            {
//                moris::sint tIndex = tIndices( Ik );
//
//                moris::uint tNumEnrichLevel = mBasisEnrichmentIndices( tIndex ).numel();
//
//                tNumMaxEnrichLevel = std::max( tNumMaxEnrichLevel, tNumEnrichLevel );
//            }
//
//            mBasisEnrichmentIndices( Ii )  .resize( tNumMaxEnrichLevel, 1 );
//            mBasisEnrichmentBulkPhase( Ii ).resize( tNumMaxEnrichLevel, 1 );
//
//            moris::uint tBasisLevel = tMeshPointer->get_HMR_database()->get_bspline_mesh_by_index( tMeshIndex )
//                                                                      ->get_basis_by_index( Ii )
//                                                                      ->get_level();
//
//            for ( moris::uint Ik = 0; Ik < tNumMaxEnrichLevel; Ik++ )
//            {
//                mBasisEnrichmentIndices( Ii )( Ik )   = tMaxEnrichmentIndex + (tCounter++);
//                mBasisEnrichmentBulkPhase( Ii )( Ik ) = Ik;
//
//                mLevelOfEnrichedMultilevelBasis( mBasisEnrichmentIndices( Ii )( Ik ), 0 ) = tBasisLevel;
//            }
//        }
//
//    }
//
//    mEnrichedMultilevelBasis       .resize( tCounter + tMaxEnrichmentIndex, 1);
//    mLevelOfEnrichedMultilevelBasis.resize( tCounter + tMaxEnrichmentIndex, 1);
//
//    mEnrichmentToBasisIndex.set_size(tCounter + tMaxEnrichmentIndex, 1, -1 );
//    mEnrichmentToBulk      .set_size(tCounter + tMaxEnrichmentIndex, 1, -1 );
//
//    for ( moris::uint Ii = 0; Ii < mBasisEnrichmentIndices.size(); Ii++ )
//    {
//        for ( moris::uint Ia = 0; Ia < mBasisEnrichmentIndices( Ii ).numel(); Ia++ )
//        {
//             mEnrichmentToBasisIndex( mBasisEnrichmentIndices( Ii )( Ia ) ) = Ii;
//             mEnrichmentToBulk( mBasisEnrichmentIndices( Ii )( Ia ) ) = mBasisEnrichmentBulkPhase( Ii )( Ia );
//        }
//    }
//
//    print(mLevelOfEnrichedMultilevelBasis, "mLevelOfEnrichedMultilevelBasis");
//
//    print(mBasisEnrichmentIndices,"mBasisEnrichmentIndices");
//
//    print(mEnrichmentToBasisIndex,"mEnrichmentToBasisIndex");
//
//    this->create_multilevel_children_to_parent_relations();
//
//    this->create_multilevel_parent_to_children_relations();
//
//    print(mChildrenToParents,"mChildrenToParents");
//
//    print(mParentsToChildren,"mParentsToChildren");
//
//
//    //---------------------------------------------------------------------------------
//
////
////    tMeshPointer->get_HMR_database()->get_bspline_mesh_by_index( tAdofOrderHack )->get_children_ind_for_basis( tExtDofInd );
////
////
////
////    moris::uint tDofLevel = mMesh->get_HMR_database()->get_bspline_mesh_by_index( tAdofOrderHack )
////                                                     ->get_basis_by_index( tExtDofInd )
////                                                     ->get_level();
////
////    moris::sint tMaxIndexForOrder = mMesh->get_HMR_database()->get_bspline_mesh_by_index( tMeshIndex )
////                                                             ->get_number_of_indexed_basis();
//
//}
//
//void
//Enrichment::create_multilevel_children_to_parent_relations()
//{
//    moris::sint tMeshIndex = 1;
//
//    std::shared_ptr< moris::mtk::Mesh > tMeshPointer = mXTKModelPtr->mHMRMesh;
//
//    moris::uint tNumEnrichmentIndices= mEnrichmentToBasisIndex.numel();
//
//    mChildrenToParents.resize( tNumEnrichmentIndices );
//
//    for ( moris::uint Ii = 0; Ii < tNumEnrichmentIndices; Ii++ )
//    {
//        if( mLevelOfEnrichedMultilevelBasis( Ii ) != 0 )
//        {
//            moris::sint tBulkPhase = mEnrichmentToBulk( Ii );
//
//            moris::sint tBasisIndex = mEnrichmentToBasisIndex( Ii );
//
//            moris::sint tNumOfParents = tMeshPointer->get_HMR_database()
//                                                     ->get_bspline_mesh_by_index( tMeshIndex )
//                                                     ->get_basis_by_index( tBasisIndex )
//                                                     ->get_number_of_parents();
//
//            mChildrenToParents( Ii ).set_size( tNumOfParents, 1, -1 );
//
//            for ( moris::sint Ia = 0; Ia < tNumOfParents; Ia++ )
//            {
//                moris::sint tParentBasisIndex = tMeshPointer->get_HMR_database()->get_bspline_mesh_by_index( tMeshIndex )
//                                                            ->get_basis_by_index( tBasisIndex )
//                                                            ->get_parent( Ia )->get_index();
//
//                for ( moris::uint Ik = 0; Ik < mBasisEnrichmentIndices( tParentBasisIndex ).numel(); Ik++ )
//                {
//                    if ( mBasisEnrichmentBulkPhase( tParentBasisIndex )( Ik ) == tBulkPhase )
//                    {
//                        mChildrenToParents( Ii )( Ia ) = mBasisEnrichmentIndices( tParentBasisIndex )( Ik );
//                    }
//                }
//            }
//        }
//    }
//}
//
//void
//Enrichment::create_multilevel_parent_to_children_relations()
//{
//    moris::sint tMeshIndex = 1;
//
//    std::shared_ptr< moris::mtk::Mesh > tMeshPointer = mXTKModelPtr->mHMRMesh;
//
//    moris::sint tMaxMeshLevel = tMeshPointer->get_HMR_database()
//                                            ->get_bspline_mesh_by_index( tMeshIndex )
//                                            ->get_max_level();
//
//    moris::uint tNumEnrichmentIndices= mEnrichmentToBasisIndex.numel();
//
//    mParentsToChildren.resize( tNumEnrichmentIndices );
//
//    for ( moris::uint Ii = 0; Ii < tNumEnrichmentIndices; Ii++ )
//    {
//        if( mLevelOfEnrichedMultilevelBasis( Ii ) != tMaxMeshLevel )
//        {
//            moris::sint tBulkPhase = mEnrichmentToBulk( Ii );
//
//            moris::sint tBasisIndex = mEnrichmentToBasisIndex( Ii );
//
//            // Get vector with external fine indices
//            moris::Matrix< DDSMat > tIndices = tMeshPointer->get_HMR_database()
//                                                           ->get_bspline_mesh_by_index( tMeshIndex )
//                                                           ->get_children_ind_for_basis( tBasisIndex );
//
//            moris::sint tNumChildren = tIndices.numel();
//
//            mParentsToChildren( Ii ).set_size( tNumChildren, 1, -1 );
//
//            for ( moris::sint Ia = 0; Ia < tNumChildren; Ia++ )
//            {
//                moris::sint tChildBasisIndex = tIndices( Ia );
//
//                for ( moris::uint Ik = 0; Ik < mBasisEnrichmentIndices( tChildBasisIndex ).numel(); Ik++ )
//                {
//                    if ( mBasisEnrichmentBulkPhase( tChildBasisIndex )( Ik ) == tBulkPhase )
//                    {
//                        //mVertexEnrichments( Ii )
//                        mParentsToChildren( Ii )( Ia ) = mBasisEnrichmentIndices( tChildBasisIndex )( Ik );
//                    }
//                }
//            }
//        }
//    }
//}


}


