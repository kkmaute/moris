/*
 * cl_XTK_Enrichment.cpp
 *
 *  Created on: Feb 18, 2019
 *      Author: doble
 */

#include "cl_XTK_Enrichment.hpp"

#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_map>

#include "cl_Communication_Tools.hpp"
#include "linalg_typedefs.hpp"
#include "fn_assert.hpp"
#include "fn_sort.hpp"
#include "typedefs.hpp"
#include "cl_MTK_Cell_Info_Factory.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh_Core.hpp"
#include "cl_MTK_Vertex.hpp"
#include "xtk_typedefs.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
#include "cl_XTK_Interpolation_Cell_Unzipped.hpp"
#include "cl_XTK_Interpolation_Vertex_Unzipped.hpp"
#include "cl_XTK_Cut_Integration_Mesh.hpp"
#include "cl_XTK_Integration_Mesh_Generator.hpp"

#include "cl_Tracer.hpp"
//#include "cl_HMR_Database.hpp"
namespace xtk
{
//-------------------------------------------------------------------------------------

Enrichment::Enrichment(
    enum Enrichment_Method const& aMethod,
    enum EntityRank const&        aBasisRank,
    Matrix<IndexMat> const&       aInterpIndex,
    moris::moris_index const&     aNumBulkPhases,
    xtk::Model*                   aXTKModelPtr,
    moris::mtk::Mesh*         aBackgroundMeshPtr)
    : mEnrichmentMethod(aMethod)
    , mBasisRank(aBasisRank)
    , mMeshIndices(aInterpIndex)
    , mNumBulkPhases(aNumBulkPhases)
    , mXTKModelPtr(aXTKModelPtr)
    , mBackgroundMeshPtr(aBackgroundMeshPtr)
    , mEnrichmentData(aInterpIndex.max() + 1, mXTKModelPtr->get_cut_integration_mesh()->get_num_subphases())
{
    mCutIgMesh   = mXTKModelPtr->get_cut_integration_mesh();
    mIgMeshTools = new Integration_Mesh_Generator();
}

//-------------------------------------------------------------------------------------

void
Enrichment::perform_enrichment()
{
    MORIS_ERROR(mBackgroundMeshPtr != nullptr,
        "mBackgroundMesh nullptr detected, this is probably because the enrichment has not been initialized properly");

    // Perform enrichment over basis clusters
    perform_basis_cluster_enrichment();
}

//-------------------------------------------------------------------------------------

Cell<moris::Matrix<moris::IdMat>> const&
Enrichment::get_element_inds_in_basis_support(moris_index const& aEnrichmentDataIndex) const
{
    return mEnrichmentData(aEnrichmentDataIndex).mElementIndsInBasis;
}

//-------------------------------------------------------------------------------------

Cell<moris::Matrix<moris::IndexMat>> const&
Enrichment::get_element_enrichment_levels_in_basis_support(moris_index const& aEnrichmentDataIndex) const
{
    return mEnrichmentData(aEnrichmentDataIndex).mElementEnrichmentLevel;
}

//-------------------------------------------------------------------------------------

Cell<moris::Matrix<moris::IndexMat>> const&
Enrichment::get_subphases_loc_inds_in_enriched_basis(moris_index const& aEnrichmentDataIndex) const
{
    return mEnrichmentData(aEnrichmentDataIndex).mSubphaseIndsInEnrichedBasis;
}

//-------------------------------------------------------------------------------------

moris::Memory_Map
Enrichment::get_memory_usage()
{
    Memory_Map tMemoryMap;

    tMemoryMap.mMemoryMapData["mEnrichmentMethod"] = sizeof(mEnrichmentMethod);
    tMemoryMap.mMemoryMapData["mBasisRank"]        = sizeof(mBasisRank);
    tMemoryMap.mMemoryMapData["mMeshIndices"]      = mMeshIndices.capacity();
    tMemoryMap.mMemoryMapData["mNumBulkPhases"]    = sizeof(mNumBulkPhases);
    tMemoryMap.mMemoryMapData["mMeshIndices"]      = mMeshIndices.capacity();

    tMemoryMap.mMemoryMapData["mElementEnrichmentLevel"]      = moris::internal_capacity(mEnrichmentData(0).mElementEnrichmentLevel);
    tMemoryMap.mMemoryMapData["mElementIndsInBasis"]          = moris::internal_capacity(mEnrichmentData(0).mElementIndsInBasis);
    tMemoryMap.mMemoryMapData["mSubphaseIndsInEnrichedBasis"] = moris::internal_capacity(mEnrichmentData(0).mSubphaseIndsInEnrichedBasis);
    tMemoryMap.mMemoryMapData["mBasisEnrichmentIndices"]      = moris::internal_capacity(mEnrichmentData(0).mBasisEnrichmentIndices);
    tMemoryMap.mMemoryMapData["mEnrichedBasisIndexToId"]      = mEnrichmentData(0).mEnrichedBasisIndexToId.capacity();
    tMemoryMap.mMemoryMapData["mSubphaseBGBasisIndices"]      = moris::internal_capacity(mEnrichmentData(0).mSubphaseBGBasisIndices);
    tMemoryMap.mMemoryMapData["mSubphaseBGBasisEnrLev"]       = moris::internal_capacity(mEnrichmentData(0).mSubphaseBGBasisEnrLev);
    tMemoryMap.mMemoryMapData["mNumEnrichmentLevels"]         = sizeof(mEnrichmentData(0).mNumEnrichmentLevels);
    tMemoryMap.mMemoryMapData["mBGVertexInterpolations ptrs"] = mEnrichmentData(0).mBGVertexInterpolations.capacity();

    return tMemoryMap;
}

void
Enrichment::write_diagnostics()
{
    for (moris::size_t iBasisType = 0; iBasisType < mMeshIndices.numel(); iBasisType++)
    {
        // get the mesh index
        moris_index tMeshIndex = mMeshIndices(iBasisType);

        std::string tEnrBasisIdToSubphaseId = mXTKModelPtr->get_diagnostic_file_name( std::string( "Enr_Basis_To_Subphase_" + std::to_string(tMeshIndex)) );
    
        this->print_enriched_basis_to_subphase_id(tMeshIndex,tEnrBasisIdToSubphaseId);
    }

    mXTKModelPtr->mEnrichedInterpMesh(0)->write_diagnostics();
}
//-------------------------------------------------------------------------------------

void
Enrichment::print_enriched_basis_to_subphase_id(
    const moris_index& aMeshIndex,
    std::string        aFileName )
{
    
    Cell< moris::Matrix< moris::IndexMat > > const& tSubphasesInEnrBasis = mEnrichmentData( aMeshIndex ).mSubphaseIndsInEnrichedBasis;

    std::ostringstream tStringStream;
    tStringStream << "Enriched_Basis_Id,";
    tStringStream << "Mesh_Index,";
    tStringStream << "Owner,";
    tStringStream << "PRank,";

    // global max size of
    moris_index tLocalMaxIGCellGroupSize = 0;
    for ( moris::uint iEB = 0; iEB < tSubphasesInEnrBasis.size(); iEB++ )
    {
        if ( (moris_index) tSubphasesInEnrBasis( iEB ).numel() > tLocalMaxIGCellGroupSize )
        {
            tLocalMaxIGCellGroupSize = tSubphasesInEnrBasis( iEB ).numel();
        }
    }

    moris_index tGlbMaxIgCellGroupSize = moris::max_all( tLocalMaxIGCellGroupSize );

    for ( moris_index iCH = 0; iCH < tGlbMaxIgCellGroupSize; iCH++ )
    {
        tStringStream << "Subphase_ID" + std::to_string( iCH );

        if ( iCH != tGlbMaxIgCellGroupSize - 1 )
        {
            tStringStream << ",";
        }
    }

    tStringStream << "\n";


    Matrix< IndexMat > const& tLocToGlbEnrBasisId = mXTKModelPtr->mEnrichedInterpMesh( 0 )->get_enriched_coefficient_local_to_global_map( aMeshIndex );

    for ( moris::uint iEB = 0; iEB < tSubphasesInEnrBasis.size(); iEB++ )
    {
        tStringStream << std::to_string( tLocToGlbEnrBasisId( iEB ) ) << ",";
        tStringStream << std::to_string( aMeshIndex ) << ",";
        tStringStream << std::to_string( mXTKModelPtr->mEnrichedInterpMesh( 0 )->get_basis_owner( (moris_index)iEB, aMeshIndex ) ) << ",";
        tStringStream << std::to_string( par_rank() ) << ",";
        for ( size_t iSp = 0; iSp < tSubphasesInEnrBasis( iEB ).numel(); iSp++ )
        {
            tStringStream << std::to_string( mXTKModelPtr->get_cut_integration_mesh()->get_subphase_id( tSubphasesInEnrBasis( iEB )( iSp ) ) );

            if ( iSp != tSubphasesInEnrBasis( iEB ).numel() - 1 )
            {
                tStringStream << ",";
            }
        }
        tStringStream << "\n";
    }

    if ( aFileName.empty() == false )
    {
        std::ofstream tOutputFile( aFileName );
        tOutputFile << tStringStream.str() << std::endl;
        tOutputFile.close();
    }
}

//-------------------------------------------------------------------------------------

void
Enrichment::perform_basis_cluster_enrichment()
{
    MORIS_ASSERT( mBackgroundMeshPtr->get_num_elems() > 0, "0 cell interpolation mesh passed" );

    // construct cell in xtk conformal model neighborhood connectivity
    this->construct_neighborhoods();

    // make sure we have access to all the vertex interpolation
    this->setup_background_vertex_interpolations();

    // iterate through basis types (i.e. linear and quadratic)
    for ( moris::size_t iBasisType = 0; iBasisType < mMeshIndices.numel(); iBasisType++ )
    {
        // get the mesh index
        moris_index tMeshIndex = mMeshIndices( iBasisType );

        Tracer tTracer( "XTK", "Enrichment", "Mesh Index " + std::to_string( tMeshIndex ) );

        // Number of basis functions (= number of B-Splines on the A-mesh ?)
        moris::size_t tNumBasis = mBackgroundMeshPtr->get_num_basis_functions( tMeshIndex );

        // allocate member variables
        mEnrichmentData( tMeshIndex ).mElementEnrichmentLevel = moris::Cell< moris::Matrix< moris::IndexMat > >( tNumBasis );
        mEnrichmentData( tMeshIndex ).mElementIndsInBasis     = moris::Cell< moris::Matrix< moris::IndexMat > >( tNumBasis );

        // allocate data used after basis loop
        moris::Cell< moris::Matrix< moris::IndexMat > > tSubPhaseBinEnrichment( tNumBasis );
        moris::Cell< moris::Matrix< moris::IndexMat > > tSubphaseClusterIndicesInSupport( tNumBasis );
        moris::Cell< moris_index >                      tMaxEnrichmentLevel( tNumBasis, 0 );

        for ( moris::size_t iBasis = 0; iBasis < tNumBasis; iBasis++ )
        {
            // Get elements in support of basis (these are interpolation cells)
            moris::Matrix< moris::IndexMat > tParentElementsInSupport;

            mBackgroundMeshPtr->get_elements_in_support_of_basis( tMeshIndex, iBasis, tParentElementsInSupport );

            // get subphase clusters in support (separated by phase)
            tSubphaseClusterIndicesInSupport( iBasis ) = this->get_subphase_clusters_in_support( tParentElementsInSupport );

            // construct subphase in support map
            IndexMap tSubPhaseIndexToSupportIndex;

            this->construct_subphase_in_support_map( tSubphaseClusterIndicesInSupport( iBasis ), tSubPhaseIndexToSupportIndex );

            // prune the subphase to remove subphases outside of basis support
            moris::Matrix< moris::IndexMat > tPrunedSubphaseNeighborhood;

            this->generate_pruned_subphase_graph_in_basis_support(
                tSubphaseClusterIndicesInSupport( iBasis ),
                tSubPhaseIndexToSupportIndex,
                tPrunedSubphaseNeighborhood );

            // Assign enrichment levels to subphases
            this->assign_subphase_bin_enrichment_levels_in_basis_support(
                tSubphaseClusterIndicesInSupport( iBasis ),
                tSubPhaseIndexToSupportIndex,
                tPrunedSubphaseNeighborhood,
                tSubPhaseBinEnrichment( iBasis ),
                tMaxEnrichmentLevel( iBasis ) );

            // Extract element enrichment levels from assigned sub-phase bin enrichment levels and store these as a member variable
            this->unzip_subphase_bin_enrichment_into_element_enrichment(
                tMeshIndex,
                i,
                tParentElementsInSupport,
                tSubphaseClusterIndicesInSupport( iBasis ),
                tSubPhaseIndexToSupportIndex,
                tPrunedSubphaseNeighborhood,
                tSubPhaseBinEnrichment( iBasis ) );
        }

        // construct subphase to enriched index
        this->construct_enriched_basis_to_subphase_connectivity(
            tMeshIndex,
            tSubPhaseBinEnrichment,
            tSubphaseClusterIndicesInSupport,
            tMaxEnrichmentLevel );

        // Assign enriched basis indices Only indices here because interpolation cells needs basis
        // indices and basis ids are created using the interpolation cell ids
        this->assign_enriched_coefficients_identifiers(
            tMeshIndex,
            tMaxEnrichmentLevel );

        MORIS_LOG_SPEC( "Num Enriched Basis", mEnrichmentData( tMeshIndex ).mNumEnrichmentLevels );
    }

    // create the enriched interpolation mesh
    this->construct_enriched_interpolation_mesh();

    // create the integration mesh
    this->construct_enriched_integration_mesh();
}

//-------------------------------------------------------------------------------------

void
Enrichment::construct_neighborhoods()
{
    // // construct full mesh neighborhood
    // mXTKModelPtr->construct_neighborhood();

    // // construct subphase neighborhood
    // mXTKModelPtr->construct_subphase_neighborhood();
}

//-------------------------------------------------------------------------------------

void
Enrichment::setup_background_vertex_interpolations()
{
    // access the background interpolation mesh
    mtk::Interpolation_Mesh& tBGIPMesh = mXTKModelPtr->get_background_mesh();

    // size the data in the enrichment data
    moris::uint tNumVertices = tBGIPMesh.get_num_nodes();

    // iterate through meshes
    for (moris::uint iM = 0; iM < mMeshIndices.numel(); iM++)
    {
        // size the vertex interpolations (allocating everything as null)
        mEnrichmentData(mMeshIndices(iM)).mBGVertexInterpolations.resize(tNumVertices, nullptr);

        // iterate through vertices
        for (moris::uint iV = 0; iV < tNumVertices; iV++)
        {
            // access the vertex
            moris::mtk::Vertex& tVert = tBGIPMesh.get_mtk_vertex((moris_index)iV);

            // if the vertex interpolation exits this vertex is own by the current processor
            // and the vertex interpolation is stored; otherwise the nullptr is kept indicating
            // that the vertex is not owned the by current processor
            if (tVert.has_interpolation(mMeshIndices(iM)))
            {
                mEnrichmentData(mMeshIndices(iM)).mBGVertexInterpolations(iV) = tVert.get_interpolation(mMeshIndices(iM));
            }
        }
    }
}

//-------------------------------------------------------------------------------------

Matrix<IndexMat>
Enrichment::get_subphase_clusters_in_support(moris::Matrix<moris::IndexMat> const& aElementsInSupport)
{

    // count the number of subphase cluster in support
    moris::uint tCount = 0;

    for (moris::size_t iE = 0; iE < aElementsInSupport.numel(); iE++)
    {
        moris::Cell<moris_index> const& tSubphaseIndices = mCutIgMesh->get_parent_cell_subphases(aElementsInSupport(iE));

        tCount = tCount + tSubphaseIndices.size();
    }

    Matrix<IndexMat> tSubPhaseClusters(1, tCount);
    tCount = 0;


    for (moris::size_t iE = 0; iE < aElementsInSupport.numel(); iE++)
    {

        moris::Cell<moris_index> const& tSubphaseIndices = mCutIgMesh->get_parent_cell_subphases(aElementsInSupport(iE));

        for (moris::uint iSP = 0; iSP < tSubphaseIndices.size(); iSP++)
        {
            tSubPhaseClusters(tCount++) = tSubphaseIndices(iSP);
        }
    }

    return tSubPhaseClusters;
}

//-------------------------------------------------------------------------------------

void
Enrichment::construct_subphase_in_support_map(
    moris::Matrix<moris::IndexMat> const& aSubphaseClusterIndicesInSupport,
    IndexMap&                             aSubPhaseIndexToSupportIndex)
{
    for (moris::moris_index i = 0; i < (moris::moris_index)aSubphaseClusterIndicesInSupport.numel(); i++)
    {
        aSubPhaseIndexToSupportIndex[aSubphaseClusterIndicesInSupport(i)] = i;
    }
}

//-------------------------------------------------------------------------------------

void
Enrichment::generate_pruned_subphase_graph_in_basis_support(
    moris::Matrix< moris::IndexMat > const& aSubphasesInSupport,
    IndexMap&                               aSubPhaseIndexToSupportIndex,
    moris::Matrix< moris::IndexMat >&       aPrunedSubPhaseToSubphase )
{
    std::shared_ptr<Subphase_Neighborhood_Connectivity> tSubphaseNeighborhood = mCutIgMesh->get_subphase_neighborhood();

    // Construct full element neighbor graph in support and the corresponding shared faces
    aPrunedSubPhaseToSubphase.resize(aSubphasesInSupport.numel(), 50);// FIXME: Keenan this allocation needs to done smarter
    aPrunedSubPhaseToSubphase.fill(MORIS_INDEX_MAX);

    // get subphase neighborhood information
    moris::Cell<std::shared_ptr<moris::Cell<moris_index>>> const& tSubphasetoSubphase = tSubphaseNeighborhood->mSubphaseToSubPhase;

    for (moris::size_t i = 0; i < aSubphasesInSupport.numel(); i++)
    {
        moris::Cell<moris_index> const& tSingleSubPhaseNeighbors = *tSubphasetoSubphase(aSubphasesInSupport(i));

        // iterate through and prune subphases not in support
        moris::uint tCount = 0;
        for (moris::size_t j = 0; j < tSingleSubPhaseNeighbors.size(); j++)
        {
            moris_index tNeighborSubphaseIndex = tSingleSubPhaseNeighbors(j);

            auto tNeighborIter = aSubPhaseIndexToSupportIndex.find(tNeighborSubphaseIndex);

            if (tNeighborIter != aSubPhaseIndexToSupportIndex.end())
            {
                aPrunedSubPhaseToSubphase(i, tCount) = tNeighborIter->second;
                tCount++;
            }
        }
    }
}

//-------------------------------------------------------------------------------------

void
Enrichment::assign_subphase_bin_enrichment_levels_in_basis_support(
    moris::Matrix<moris::IndexMat> const& aSubphasesInSupport,
    IndexMap&                             aSubPhaseIndexToSupportIndex,
    moris::Matrix<moris::IndexMat> const& aPrunedSubPhaseToSubphase,
    moris::Matrix<moris::IndexMat>&       aSubPhaseBinEnrichmentVals,
    moris_index&                          aMaxEnrichmentLevel)
{
    // Variables needed for floodfill, consider removing these.
    // Active bins to include in floodfill (We include all bins)
    moris::Matrix<moris::IndexMat> tActiveBins(1, aPrunedSubPhaseToSubphase.n_rows());

    for (moris::size_t i = 0; i < aPrunedSubPhaseToSubphase.n_rows(); i++)
    {
        (tActiveBins)(0, i) = i;
    }

    // Mark all as included
    moris::Matrix<moris::IndexMat> tIncludedBins(1, aSubphasesInSupport.numel(), 1);

    // Flood fill metric value (since all the subphases do not connect to dissimilar phases)
    moris::Matrix<moris::IndexMat> tDummyPhase(1, aSubphasesInSupport.numel(), 1);

    aSubPhaseBinEnrichmentVals = flood_fill(
        aPrunedSubPhaseToSubphase,
        tDummyPhase,
        tActiveBins,
        tIncludedBins,
        mNumBulkPhases,
        MORIS_INDEX_MAX,
        aMaxEnrichmentLevel,
        true);
}

//-------------------------------------------------------------------------------------

void
Enrichment::unzip_subphase_bin_enrichment_into_element_enrichment(
    moris_index const&                    aEnrichmentDataIndex,
    moris_index const&                    aBasisIndex,
    moris::Matrix<moris::IndexMat> const& aParentElementsInSupport,
    moris::Matrix<moris::IndexMat> const& aSubphasesInSupport,
    IndexMap&                             aSubPhaseIndexToSupportIndex,
    moris::Matrix<moris::IndexMat> const& aPrunedSubPhaseToSubphase,
    moris::Matrix<moris::IndexMat>&       aSubPhaseBinEnrichmentVals)
{
    // resize member data
    moris::size_t tNumAllElementsInSupport                                     = this->count_elements_in_support(aParentElementsInSupport);
    mEnrichmentData(aEnrichmentDataIndex).mElementIndsInBasis(aBasisIndex)     = moris::Matrix<moris::IndexMat>(1, tNumAllElementsInSupport);
    mEnrichmentData(aEnrichmentDataIndex).mElementEnrichmentLevel(aBasisIndex) = moris::Matrix<moris::IndexMat>(1, tNumAllElementsInSupport);

    moris::uint tCount = 0;

    for (moris::size_t i = 0; i < aSubphasesInSupport.numel(); i++)
    {
        moris_index tSubphaseIndex = aSubphasesInSupport(i);

        // iterate through cell sin the subphase
        std::shared_ptr<IG_Cell_Group> tSubphaseIgCells = mCutIgMesh->get_subphase_ig_cells(tSubphaseIndex);

        // iterate through cell sin subphase
        for (moris::uint iSPCell = 0; iSPCell < tSubphaseIgCells->mIgCellGroup.size(); iSPCell++)
        {
            mEnrichmentData(aEnrichmentDataIndex).mElementIndsInBasis(aBasisIndex)(tCount)     = tSubphaseIgCells->mIgCellGroup(iSPCell)->get_index();
            mEnrichmentData(aEnrichmentDataIndex).mElementEnrichmentLevel(aBasisIndex)(tCount) = aSubPhaseBinEnrichmentVals(i);
            tCount++;
        }

        // add information to interp cells about which basis/enrichment level interpolates in it
        mEnrichmentData(aEnrichmentDataIndex).mSubphaseBGBasisIndices(tSubphaseIndex).push_back(aBasisIndex);
        mEnrichmentData(aEnrichmentDataIndex).mSubphaseBGBasisEnrLev(tSubphaseIndex).push_back(aSubPhaseBinEnrichmentVals(i));
    }
}

//-------------------------------------------------------------------------------------

void
Enrichment::construct_enriched_basis_to_subphase_connectivity(
    moris_index const&                                 aEnrichmentDataIndex,
    moris::Cell<moris::Matrix<moris::IndexMat>> const& aSubPhaseBinEnrichment,
    moris::Cell<moris::Matrix<moris::IndexMat>> const& aSubphaseClusterIndicesInSupport,
    moris::Cell<moris_index> const&                    aMaxEnrichmentLevel)
{
    moris::moris_index tNumEnrichmentLevels = 0;

    for (moris::uint i = 0; i < aMaxEnrichmentLevel.size(); i++)
    {
        tNumEnrichmentLevels += aMaxEnrichmentLevel(i) + 1;
    }

    // size data
    mEnrichmentData(aEnrichmentDataIndex).mSubphaseIndsInEnrichedBasis.resize(tNumEnrichmentLevels);

    moris_index tBaseIndex = 0;

    for (moris::uint i = 0; i < aSubPhaseBinEnrichment.size(); i++)
    {
        // get the maximum enrichment level in this basis support
        moris::moris_index tMaxEnrLev = aMaxEnrichmentLevel(i);

        // counter
        Cell<moris_index> tCounter(tMaxEnrLev + 1, 0);

        // allocate member data for these basis functions
        for (moris::moris_index iEnr = tBaseIndex; iEnr < tBaseIndex + tMaxEnrLev + 1; iEnr++)
        {
            mEnrichmentData(aEnrichmentDataIndex).mSubphaseIndsInEnrichedBasis(iEnr).resize(1, aSubphaseClusterIndicesInSupport(i).numel());
        }

        // iterate through subphases in support and add them to appropriate location in mSubphaseIndsInEnrichedBasis
        for (moris::uint iSp = 0; iSp < aSubphaseClusterIndicesInSupport(i).numel(); iSp++)
        {
            // get cluster enrichment level
            moris_index tClusterEnrLev = aSubPhaseBinEnrichment(i)(iSp);

            // add to the member data
            mEnrichmentData(aEnrichmentDataIndex).mSubphaseIndsInEnrichedBasis(tBaseIndex + tClusterEnrLev)(tCounter(tClusterEnrLev)) =
                aSubphaseClusterIndicesInSupport(i)(iSp);

            // increment count
            tCounter(tClusterEnrLev)++;
        }

        // size out unused space
        for (moris::moris_index iEnr = 0; iEnr < tMaxEnrLev + 1; iEnr++)
        {
            moris_index tIndex = tBaseIndex + iEnr;
            mEnrichmentData(aEnrichmentDataIndex).mSubphaseIndsInEnrichedBasis(tIndex).resize(1, tCounter(iEnr));

            // sort in ascending order (easier to find in MPI)
            // if this sort is removed the function  subphase_is_in_support needs to be updated
            moris::sort(mEnrichmentData(aEnrichmentDataIndex).mSubphaseIndsInEnrichedBasis(tIndex),
                mEnrichmentData(aEnrichmentDataIndex).mSubphaseIndsInEnrichedBasis(tIndex),
                "ascend",
                1);
        }

        // update starting index
        tBaseIndex = tBaseIndex + tMaxEnrLev + 1;
    }
}

//-------------------------------------------------------------------------------------

void
Enrichment::assign_enriched_coefficients_identifiers(
    moris_index const&              aEnrichmentDataIndex,
    moris::Cell<moris_index> const& aMaxEnrichmentLevel)
{
    mEnrichmentData(aEnrichmentDataIndex).mBasisEnrichmentIndices.resize(mEnrichmentData(aEnrichmentDataIndex).mElementIndsInBasis.size());

    mEnrichmentData(aEnrichmentDataIndex).mNumEnrichmentLevels = 0;

    for (moris::uint i = 0; i < mEnrichmentData(aEnrichmentDataIndex).mElementIndsInBasis.size(); i++)
    {
        moris::moris_index tMaxEnrLev                                    = aMaxEnrichmentLevel(i) + 1;
        mEnrichmentData(aEnrichmentDataIndex).mNumEnrichmentLevels       = mEnrichmentData(aEnrichmentDataIndex).mNumEnrichmentLevels + tMaxEnrLev;
        mEnrichmentData(aEnrichmentDataIndex).mBasisEnrichmentIndices(i) = moris::Matrix<moris::IndexMat>(1, tMaxEnrLev);
    }

    // allocate enriched basis index to id data
    mEnrichmentData(aEnrichmentDataIndex).mEnrichedBasisIndexToId.resize(1, mEnrichmentData(aEnrichmentDataIndex).mNumEnrichmentLevels);
    mEnrichmentData(aEnrichmentDataIndex).mEnrichedBasisIndexToId.fill(MORIS_INDEX_MAX);

    moris_index tParRank = par_rank();

    // get the comm table
    Matrix<IndexMat> tCommTable = mXTKModelPtr->get_communication_table();

    // Procs CEll
    Cell<moris_index> tProcRanks(tCommTable.numel());

    std::unordered_map<moris_id, moris_id> tProcRankToIndexInData;

    // resize proc ranks and setup map to comm table
    for (moris::uint i = 0; i < tCommTable.numel(); i++)
    {
        tProcRankToIndexInData[tCommTable(i)] = i;
        tProcRanks(i)                         = (tCommTable(i));
    }

    // first index and id ( not first gets background basis information)
    moris::moris_index tIndOffset     = 0;
    moris::moris_id    tBasisIdOffset = this->allocate_basis_ids(aEnrichmentDataIndex, mEnrichmentData(aEnrichmentDataIndex).mNumEnrichmentLevels);

    Cell<Cell<moris_index>> tBasisIdToBasisOwner(tCommTable.numel());
    Cell<Cell<moris_index>> tSubphaseIdInSupport(tCommTable.numel());
    Cell<Cell<moris_index>> tBasisIndexToBasisOwner(tCommTable.numel());

    for (moris::uint i = 0; i < mEnrichmentData(aEnrichmentDataIndex).mBasisEnrichmentIndices.size(); i++)
    {
        moris_index tOwner = mBackgroundMeshPtr->get_entity_owner(
            i,
            mBasisRank,
            aEnrichmentDataIndex);

        moris_id tBackBasisId =
            mBackgroundMeshPtr->get_glb_entity_id_from_entity_loc_index(
                i,
                mBasisRank,
                aEnrichmentDataIndex);

        moris_index tProcDataIndex = tProcRankToIndexInData[tOwner];

        // always set indices
        moris::Matrix<moris::IndexMat>& tBasisEnrichmentInds = mEnrichmentData(aEnrichmentDataIndex).mBasisEnrichmentIndices(i);

        for (moris::uint j = 0; j < tBasisEnrichmentInds.numel(); j++)
        {
            tBasisEnrichmentInds(j) = tIndOffset;
            tIndOffset++;
        }

        // only set id if we own it and package data for communication if shared
        if (tOwner == tParRank)
        {
            mEnrichmentData(aEnrichmentDataIndex).mEnrichedBasisIndexToId(tBasisEnrichmentInds(0)) = tBackBasisId;

            for (moris::uint j = 1; j < tBasisEnrichmentInds.numel(); j++)
            {
                MORIS_ASSERT(mEnrichmentData(aEnrichmentDataIndex).mEnrichedBasisIndexToId(tBasisEnrichmentInds(j)) == MORIS_INDEX_MAX, "Already set enriched basis id");

                mEnrichmentData(aEnrichmentDataIndex).mEnrichedBasisIndexToId(tBasisEnrichmentInds(j)) = tBasisIdOffset;

                tBasisIdOffset++;
            }
        }
        // if we don't own the basis setup the communication to get the basis
        else
        {
            for (moris::uint j = 0; j < tBasisEnrichmentInds.numel(); j++)
            {
                moris_index tEnrichedBasisIndex = tBasisEnrichmentInds(j);

                moris_index tFirstSubphaseInSupportIndex = mEnrichmentData(aEnrichmentDataIndex).mSubphaseIndsInEnrichedBasis(tEnrichedBasisIndex)(0);

                moris_index tFirstSubphaseInSupportId = mXTKModelPtr->get_subphase_id(tFirstSubphaseInSupportIndex);

                tBasisIdToBasisOwner(tProcDataIndex).push_back(tBackBasisId);
                tSubphaseIdInSupport(tProcDataIndex).push_back(tFirstSubphaseInSupportId);
                tBasisIndexToBasisOwner(tProcDataIndex).push_back(tEnrichedBasisIndex);
            }
        }
    }

    // send information about not owned enriched basis to owner processor
    Cell<moris::Matrix<moris::IndexMat>> tNotOwnedEnrichedBasisId;

    communicate_basis_information_with_owner(
        aEnrichmentDataIndex,
        tBasisIdToBasisOwner,
        tSubphaseIdInSupport,
        tProcRanks,
        tProcRankToIndexInData,
        tNotOwnedEnrichedBasisId);

    // set the received information in my data
    set_received_enriched_basis_ids(
        aEnrichmentDataIndex,
        tNotOwnedEnrichedBasisId,
        tBasisIndexToBasisOwner,
        tSubphaseIdInSupport);
}


//-------------------------------------------------------------------------------------

void
Enrichment::communicate_basis_information_with_owner(
    moris_index const&                      aEnrichmentDataIndex,
    Cell<Cell<moris_index>> const&          aBasisIdToBasisOwner,
    Cell<Cell<moris_index>> const&          aSubphaseIdInSupport,
    Cell<moris_index> const&                aProcRanks,
    std::unordered_map<moris_id, moris_id>& aProcRankToIndexInData,
    Cell<moris::Matrix<moris::IndexMat>>&   aEnrichedBasisId)
{
    // copy into a matrix
    Cell<moris::Matrix<moris::IndexMat>> tBasisIdToBasisOwnerMat(aProcRanks.size());
    Cell<moris::Matrix<moris::IndexMat>> tSubphaseIdInSupport(aProcRanks.size());

    for (moris::uint i = 0; i < aBasisIdToBasisOwner.size(); i++)
    {
        tBasisIdToBasisOwnerMat(i) = Matrix<IndexMat>(1, aBasisIdToBasisOwner(i).size());
        tSubphaseIdInSupport(i)    = Matrix<IndexMat>(1, aBasisIdToBasisOwner(i).size());

        for (moris::uint j = 0; j < aBasisIdToBasisOwner(i).size(); j++)
        {
            tBasisIdToBasisOwnerMat(i)(j) = aBasisIdToBasisOwner(i)(j);
            tSubphaseIdInSupport(i)(j)    = aSubphaseIdInSupport(i)(j);
        }

        if (aBasisIdToBasisOwner(i).size() == 0)
        {
            tBasisIdToBasisOwnerMat(i).resize(1, 1);
            tBasisIdToBasisOwnerMat(i)(0) = MORIS_INDEX_MAX;

            tSubphaseIdInSupport(i).resize(1, 1);
            tSubphaseIdInSupport(i)(0) = MORIS_INDEX_MAX;
        }
    }

    // send to processor (no particular reason for these tag number choices)
    moris_index tBasisIdsTag    = 1119;
    moris_index tMaxSubphasetag = 1120;
    moris_index tReturnIdTag    = 1121;

    for (moris::uint i = 0; i < tBasisIdToBasisOwnerMat.size(); i++)
    {
        // send child element ids
        nonblocking_send(
            tBasisIdToBasisOwnerMat(i),
            tBasisIdToBasisOwnerMat(i).n_rows(),
            tBasisIdToBasisOwnerMat(i).n_cols(),
            aProcRanks(i),
            tBasisIdsTag);

        // send element parent ids
        nonblocking_send(
            tSubphaseIdInSupport(i),
            tSubphaseIdInSupport(i).n_rows(),
            tSubphaseIdInSupport(i).n_cols(),
            aProcRanks(i),
            tMaxSubphasetag);
    }

    // make sure that all sends have been completed
    barrier();

    // receive
    moris::moris_index tNumRow = 1;

    Cell<moris::Matrix<IdMat>> tReceiveInfoBasisId(aProcRanks.size());
    Cell<moris::Matrix<IdMat>> tReceivedSubphaseId(aProcRanks.size());
    Cell<moris::Matrix<IdMat>> tEnrichedBasisIds(aProcRanks.size());
    aEnrichedBasisId.resize(aProcRanks.size());

    for (moris::uint i = 0; i < aProcRanks.size(); i++)
    {
        tReceiveInfoBasisId(i).resize(1, 1);
        tReceivedSubphaseId(i).resize(1, 1);

        receive(tReceiveInfoBasisId(i), tNumRow, aProcRanks(i), tBasisIdsTag);
        receive(tReceivedSubphaseId(i), tNumRow, aProcRanks(i), tMaxSubphasetag);

        tEnrichedBasisIds(i).resize(tReceiveInfoBasisId(i).n_rows(), tReceiveInfoBasisId(i).n_cols());
    }

    // iterate through received information and setup response information
    for (moris::uint i = 0; i < aProcRanks.size(); i++)
    {
        moris_id    tBasisId    = 0;
        moris_id    tBasisIndex = 0;
        moris_index tCount      = 0;
        aEnrichedBasisId(i).resize(tReceiveInfoBasisId(i).n_rows(), tReceiveInfoBasisId(i).n_cols());

        if (tReceiveInfoBasisId(i)(0) != MORIS_INDEX_MAX)
        {
            for (moris::uint j = 0; j < tReceiveInfoBasisId(i).n_cols(); j++)
            {
                if (tBasisId != tReceiveInfoBasisId(i)(j))
                {
                    tBasisId    = tReceiveInfoBasisId(i)(j);
                    tBasisIndex = mBackgroundMeshPtr->get_loc_entity_ind_from_entity_glb_id(
                        tBasisId,
                        mBasisRank,
                        aEnrichmentDataIndex);
                }

                moris_id    tSubphaseId    = tReceivedSubphaseId(i)(j);
                moris_index tSubphaseIndex = mXTKModelPtr->get_subphase_index(tSubphaseId);

                // iterate through the max integer ids and match

                bool tFound = false;
                for (moris::uint k = 0; k < mEnrichmentData(aEnrichmentDataIndex).mBasisEnrichmentIndices(tBasisIndex).numel(); k++)
                {
                    if (this->subphase_is_in_support(aEnrichmentDataIndex, tSubphaseIndex, mEnrichmentData(aEnrichmentDataIndex).mBasisEnrichmentIndices(tBasisIndex)(k)))
                    {
                        moris_index tEnrichedBasisIndex = mEnrichmentData(aEnrichmentDataIndex).mBasisEnrichmentIndices(tBasisIndex)(k);
                        moris_index tEnrichedBasisId    = mEnrichmentData(aEnrichmentDataIndex).mEnrichedBasisIndexToId(tEnrichedBasisIndex);

                        tEnrichedBasisIds(i)(tCount) = tEnrichedBasisId;
                        tCount++;
                        tFound = true;
                    }
                }
                MORIS_ERROR(tFound, "Basis not found");
            }
        }
        else
        {
            tEnrichedBasisIds(i)(0) = MORIS_INDEX_MAX;
        }

        // send enriched ids
        nonblocking_send(
            tEnrichedBasisIds(i),
            tEnrichedBasisIds(i).n_rows(),
            tEnrichedBasisIds(i).n_cols(),
            aProcRanks(i),
            tReturnIdTag);
    }

    barrier();

    // receive information
    for (moris::uint i = 0; i < aProcRanks.size(); i++)
    {
        aEnrichedBasisId(i).resize(1, 1);
        receive(aEnrichedBasisId(i), tNumRow, aProcRanks(i), tReturnIdTag);
    }

    barrier();
}

//-------------------------------------------------------------------------------------

void
Enrichment::set_received_enriched_basis_ids(
    moris_index const&                          aEnrichmentDataIndex,
    Cell<moris::Matrix<moris::IndexMat>> const& aReceivedEnrichedIds,
    Cell<Cell<moris_index>> const&              aBasisIndexToBasisOwner,
    Cell<Cell<moris_index>> const&              aSubphaseIdInSupport)
{
    for (moris::uint i = 0; i < aReceivedEnrichedIds.size(); i++)
    {
        if (aReceivedEnrichedIds(i)(0) != MORIS_INDEX_MAX)
        {
            MORIS_ASSERT(aReceivedEnrichedIds(i).numel() == aBasisIndexToBasisOwner(i).size(),
                "Dimension mismatch between received information and expected information");

            for (moris::uint j = 0; j < aReceivedEnrichedIds(i).numel(); j++)
            {
                moris_index tLocalBasisIndex = aBasisIndexToBasisOwner(i)(j);
                moris_index tGlobaId         = aReceivedEnrichedIds(i)(j);


                MORIS_ASSERT(mEnrichmentData(aEnrichmentDataIndex).mEnrichedBasisIndexToId(tLocalBasisIndex) == MORIS_INDEX_MAX,
                    "Id already set for this basis function");

                mEnrichmentData(aEnrichmentDataIndex).mEnrichedBasisIndexToId(tLocalBasisIndex) = tGlobaId;
            }
        }
    }
}

//-------------------------------------------------------------------------------------

moris::size_t
Enrichment::count_elements_in_support(moris::Matrix<moris::IndexMat> const& aParentElementsInSupport)
{

    // Number of elements in this support (need both parent and total)
    moris::size_t tNumParentElementsInSupport = aParentElementsInSupport.n_cols();
    moris::size_t tNumElementsInSupport       = 0;

    // Count children elements in support
    for (moris::size_t i = 0; i < tNumParentElementsInSupport; i++)
    {
        moris::Cell<moris_index> const& tSubphaseIndices = mCutIgMesh->get_parent_cell_subphases(aParentElementsInSupport(i));

        // iterate through subphases
        for (moris::uint iSP = 0; iSP < tSubphaseIndices.size(); iSP++)
        {
            tNumElementsInSupport = tNumElementsInSupport + mCutIgMesh->get_subphase_ig_cells(tSubphaseIndices(iSP))->mIgCellGroup.size();
        }
    }

    return tNumElementsInSupport;
}

//-------------------------------------------------------------------------------------


bool
Enrichment::subphase_is_in_support(
    moris_index const& aEnrichmentDataIndex,
    moris_index        aSubphaseIndex,
    moris_index        aEnrichedBasisIndex)
{
    for (moris::uint i = 0; i < mEnrichmentData(aEnrichmentDataIndex).mSubphaseIndsInEnrichedBasis(aEnrichedBasisIndex).numel(); i++)
    {
        if (mEnrichmentData(aEnrichmentDataIndex).mSubphaseIndsInEnrichedBasis(aEnrichedBasisIndex)(i) == aSubphaseIndex)
        {
            return true;
        }
        else if (mEnrichmentData(aEnrichmentDataIndex).mSubphaseIndsInEnrichedBasis(aEnrichedBasisIndex)(i) > aSubphaseIndex)
        {
            return false;
        }
    }

    return false;
}

//-------------------------------------------------------------------------------------

void
Enrichment::print_basis_support_debug(
    moris_index                           aBasisIndex,
    moris::Matrix<moris::IndexMat> const& aParentElementsInSupport,
    moris::Matrix<moris::IndexMat> const& aSubphasesInSupport,
    IndexMap&                             aSubPhaseIndexToSupportIndex,
    moris::Matrix<moris::IndexMat> const& aPrunedSubPhaseToSubphase,
    moris::Matrix<moris::IndexMat>&       aSubPhaseBinEnrichmentVals)
{
    std::cout << "--------------------------------------------------" << std::endl;
    std::cout << "Basis Index: " << aBasisIndex << std::endl;
    std::cout << "Parent Cell In Support:";
    for (moris::uint i = 0; i < aParentElementsInSupport.numel(); i++)
    {
        std::cout << std::setw(8) << aParentElementsInSupport(i);
    }
    std::cout << "\nSubphases In Support:";
    for (moris::uint i = 0; i < aSubphasesInSupport.numel(); i++)
    {
        std::cout << std::setw(8) << aSubphasesInSupport(i);
    }

    std::cout << "\nSubphase Neighborhood In Support:" << std::endl;
    for (moris::uint i = 0; i < aPrunedSubPhaseToSubphase.n_rows(); i++)
    {
        std::cout << std::setw(6) << aSubphasesInSupport(i) << " | ";

        for (moris::uint j = 0; j < aPrunedSubPhaseToSubphase.n_cols(); j++)
        {
            if (aPrunedSubPhaseToSubphase(i, j) != MORIS_INDEX_MAX)
            {
                std::cout << std::setw(6) << aSubphasesInSupport(aPrunedSubPhaseToSubphase(i, j));
            }
        }
        std::cout << std::endl;
    }

    std::cout << "\nSubphase Enrichment Level: \n";
    for (moris::uint i = 0; i < aSubPhaseBinEnrichmentVals.numel(); i++)
    {
        std::cout << std::setw(8) << aSubPhaseBinEnrichmentVals(i);
    }
    std::cout << "\n--------------------------------------------------" << std::endl;
}

//-------------------------------------------------------------------------------------

void
Enrichment::construct_enriched_interpolation_mesh()
{
    Tracer tTracer("XTK", "Enrich", "Construct Enriched Interpolation Mesh");
    // initialize a new enriched interpolation mesh
    mXTKModelPtr->mEnrichedInterpMesh(0) = new Enriched_Interpolation_Mesh(mXTKModelPtr);

    // set enriched basis rank
    mXTKModelPtr->mEnrichedInterpMesh(0)->mBasisRank = mBasisRank;

    // set mesh index
    mXTKModelPtr->mEnrichedInterpMesh(0)->mMeshIndices = mMeshIndices;
    mXTKModelPtr->mEnrichedInterpMesh(0)->setup_mesh_index_map();

    // allocate space for interpolation cells
    this->allocate_interpolation_cells();

    // construct all cells with interpolation vertices being attached to a single cell
    // this handles the case of multiple enrichments where the number of interpolation vertices vary
    this->construct_enriched_interpolation_vertices_and_cells();



    mXTKModelPtr->mEnrichedInterpMesh(0)->mCoeffToEnrichCoeffs.resize(mMeshIndices.max() + 1);
    mXTKModelPtr->mEnrichedInterpMesh(0)->mEnrichCoeffLocToGlob.resize(mMeshIndices.max() + 1);

    for (moris::uint i = 0; i < mMeshIndices.numel(); i++)
    {
        // Mesh index
        moris_index tMeshIndex = mMeshIndices(i);

        // add the coeff to enriched coefs to enriched interpolation mesh
        mXTKModelPtr->mEnrichedInterpMesh(0)->mCoeffToEnrichCoeffs(tMeshIndex) =
            mEnrichmentData(tMeshIndex).mBasisEnrichmentIndices;

        // add the local to global map
        mXTKModelPtr->mEnrichedInterpMesh(0)->mEnrichCoeffLocToGlob(tMeshIndex) =
            mEnrichmentData(tMeshIndex).mEnrichedBasisIndexToId;
    }

    
    // tell the mesh to finish setting itself up
    mXTKModelPtr->mEnrichedInterpMesh(0)->finalize_setup();

    // initialize local to global maps
    mXTKModelPtr->mEnrichedInterpMesh(0)->mLocalToGlobalMaps = Cell< Matrix< IdMat > >( 4 );
    mXTKModelPtr->mEnrichedInterpMesh(0)->mGlobaltoLobalMaps = Cell< std::unordered_map< moris_id, moris_index > >( 4 );
    mXTKModelPtr->mEnrichedInterpMesh(0)->setup_cell_maps();
    mXTKModelPtr->mEnrichedInterpMesh(0)->setup_basis_maps();
    mXTKModelPtr->mEnrichedInterpMesh(0)->assign_ip_vertex_ids();
    mXTKModelPtr->mEnrichedInterpMesh(0)->setup_vertex_maps();

    moris::Cell< mtk::Vertex* > tVerticesToCommunicate;
    for(auto & iVert: mXTKModelPtr->mEnrichedInterpMesh(0)->mEnrichedInterpVerts)
    {
        if(iVert->get_owner() != moris::par_rank() )
        {
            tVerticesToCommunicate.push_back(iVert);
        }
    }

    mXTKModelPtr->mEnrichedInterpMesh(0)->communicate_select_vertex_interpolation(tVerticesToCommunicate);

    // in most cases all the interpolation vertices are the same. We merge them back together with this call
    // post-processing to construct_enriched_interpolation_vertices_and_cells in an effort to not add complexity to the function (already too)
    mXTKModelPtr->mEnrichedInterpMesh(0)->merge_duplicate_interpolation_vertices();

    // tell the mesh to finish setting itself up
    // mXTKModelPtr->mEnrichedInterpMesh(0)->finalize_setup();
    mXTKModelPtr->mEnrichedInterpMesh(0)->mGlobaltoLobalMaps(0).clear();
    mXTKModelPtr->mEnrichedInterpMesh(0)->setup_vertex_maps();
}

//-------------------------------------------------------------------------------------

void
Enrichment::construct_enriched_integration_mesh()
{
    MORIS_ASSERT(mXTKModelPtr->mEnrichedInterpMesh(0) != nullptr,
        "No enriched interpolation mesh to link enriched integration mesh to");

    mXTKModelPtr->mEnrichedIntegMesh(0) = new Enriched_Integration_Mesh(mXTKModelPtr, 0);
}

//-------------------------------------------------------------------------------------

void
Enrichment::allocate_interpolation_cells()
{
    Tracer tTracer("XTK", "Enrich", "Allocate Interpolation Cells");
    // enriched interpolation mesh pointer
    Enriched_Interpolation_Mesh* tEnrInterpMesh = mXTKModelPtr->mEnrichedInterpMesh(0);

    // set mesh indices
    tEnrInterpMesh->mMeshIndices = mMeshIndices;

    // count how many subphases there are total in the entire mesh
    uint tNumSubphases = mCutIgMesh->get_num_subphases();

    // there is one interpolation cell per subphase
    tEnrInterpMesh->mEnrichedInterpCells.resize(tNumSubphases);

    // assuming all interpolation cells are the same
    // figure out how many vertices there are per interpolation cell
    uint tNumVertsPerCell = 0;

    if (mBackgroundMeshPtr->get_num_elems() > 0)
    {
        tNumVertsPerCell = mBackgroundMeshPtr->get_mtk_cell(0).get_number_of_vertices();
    }

    tEnrInterpMesh->mNumVertsPerInterpCell = tNumVertsPerCell;

    // allocate maximum number of enriched vertices
    tEnrInterpMesh->mEnrichedInterpVerts.resize(tNumVertsPerCell * tNumSubphases);

    // allocate the base vertices to vertex enrichment data
    tEnrInterpMesh->mBaseInterpVertToVertEnrichmentIndex.resize(mMeshIndices.max() + 1, mBackgroundMeshPtr->get_num_nodes());

    // allocate space in the vertex enrichment index to parent vertex enrichment data
    tEnrInterpMesh->mVertexEnrichmentParentVertexIndex.resize(mMeshIndices.max() + 1);

    tEnrInterpMesh->mInterpVertEnrichment.resize(mMeshIndices.max() + 1);

    // allocate base cell to enriched cell data
    tEnrInterpMesh->mBaseCelltoEnrichedCell.resize(mBackgroundMeshPtr->get_num_elems());
}

//-------------------------------------------------------------------------------------

void
Enrichment::construct_enriched_interpolation_vertices_and_cells()
{
    // enriched interpolation mesh pointer, this one is constructed here
    Enriched_Interpolation_Mesh* tEnrInterpMesh = mXTKModelPtr->mEnrichedInterpMesh(0);

    // geometry and interpolation order, limited to a single interpolation order mesh
    mtk::Cell const& tFirstCell = mBackgroundMeshPtr->get_mtk_cell(0);

    // set the interpolation mesh cell info
    mtk::Cell_Info_Factory tFactory;
    tEnrInterpMesh->mCellInfo = tFactory.create_cell_info_sp(tFirstCell.get_geometry_type(), tFirstCell.get_interpolation_order());

    // allocate indices and ids
    moris_index tCellIndex = 0;

    // maximum mesh index
    moris_index tMaxMeshIndex = mMeshIndices.max();

    // Enriched Interpolation Cell to Vertex Index
    Matrix<IndexMat> tEnrInterpCellToVertex(
        tEnrInterpMesh->get_num_elements(),
        tEnrInterpMesh->mNumVertsPerInterpCell);

    // allocate vertex indices and ids
    // NOTE: THESE ARE NOT PARALLEL IDS
    moris_index tVertId      = 1;
    uint        tVertexCount = 0;

    // iterate through subphases and construct an interpolation cell in the interpolation mesh for each one
    for (moris::uint iSP = 0; iSP < mCutIgMesh->get_num_subphases(); iSP++)
    {
        // information about this cell
        moris::mtk::Cell* tParentCell = mCutIgMesh->get_subphase_parent_cell(iSP);

        // owner
        moris_id tOwner = tParentCell->get_owner();

        // vertexes of cell
        moris::Cell<mtk::Vertex*> tVertices = tParentCell->get_vertex_pointers();

        // bulk phase
        moris_index tBulkPhase = mCutIgMesh->get_subphase_bulk_phase((moris_index)iSP);

        for (moris::uint iMT = 0; iMT < mMeshIndices.numel(); iMT++)
        {
            // Mesh Index
            moris_index tMeshIndex = mMeshIndices(iMT);

            // vertex interpolations of the parent cell
            moris::Cell<mtk::Vertex_Interpolation*> tVertexInterpolations = this->get_vertex_interpolations(*tParentCell, tMeshIndex);

            // basis interpolating into the subphase
            Cell<moris_index> const& tBasisInCell = mEnrichmentData(tMeshIndex).mSubphaseBGBasisIndices((moris_index)iSP);

            // construct a map between basis index and index relative to the subphase cluster
            std::unordered_map<moris_id, moris_id> tCellBasisMap = construct_subphase_basis_to_basis_map(tBasisInCell);

            // Subphase basis enrichment level of the current subphase in child mesh
            Cell<moris_index> const& tEnrLevOfBasis = mEnrichmentData(tMeshIndex).mSubphaseBGBasisEnrLev((moris_index)iSP);

            uint tNumVertices = tParentCell->get_number_of_vertices();

            // construct unzipped enriched vertices
            for (uint iEV = 0; iEV < tNumVertices; iEV++)
            {
                // construct vertex enrichment
                Vertex_Enrichment tVertEnrichment;

                this->construct_enriched_vertex_interpolation(
                    tMeshIndex,
                    tVertexInterpolations(iEV),
                    tEnrLevOfBasis,
                    tCellBasisMap,
                    tVertEnrichment);

                // add vertex enrichment to enriched interpolation mesh
                bool tNewVertFlag = false;

                moris_index tVertEnrichIndex = tEnrInterpMesh->add_vertex_enrichment(
                    tMeshIndex,
                    tVertices(iEV),
                    tVertEnrichment,
                    tNewVertFlag);

                // create this vertex on the first go around
                if (iMT == 0)
                {

                    // Create interpolation vertex
                    tEnrInterpMesh->mEnrichedInterpVerts(tVertexCount) =
                        new Interpolation_Vertex_Unzipped(
                            tVertices(iEV),
                            tVertId,
                            tVertexCount,
                            tVertices(iEV)->get_owner(),
                            tMeshIndex,
                            tEnrInterpMesh->get_vertex_enrichment(tMeshIndex, tVertEnrichIndex),
                            tMaxMeshIndex);

                    tEnrInterpCellToVertex(tCellIndex, iEV) = tVertexCount;
                    tVertId++;
                    tVertexCount++;
                }
                else
                {
                    // interpolation vertex index
                    moris_index tVertexIndexInIp = tEnrInterpCellToVertex(tCellIndex - 1, iEV);

                    // add the vertex interpolation
                    tEnrInterpMesh->mEnrichedInterpVerts(tVertexIndexInIp)->add_vertex_interpolation(tMeshIndex, tEnrInterpMesh->get_vertex_enrichment(tMeshIndex, tVertEnrichIndex));
                }
            }

            if (iMT == 0)
            {
                // create new enriched interpolation cell
                tEnrInterpMesh->mEnrichedInterpCells(tCellIndex) =
                    new Interpolation_Cell_Unzipped(
                        tParentCell,
                        iSP,
                        tBulkPhase,
                        mCutIgMesh->get_subphase_id(iSP),
                        tCellIndex,
                        tOwner,
                        tEnrInterpMesh->mCellInfo);

                // add enriched interpolation cell to base cell to enriched cell data
                tEnrInterpMesh->mBaseCelltoEnrichedCell(
                                  tParentCell->get_index())
                    .push_back(tEnrInterpMesh->mEnrichedInterpCells(tCellIndex));

                // increment the cell index/id
                tCellIndex++;
            }
        }
    }

    //resize out aura cells
    tEnrInterpCellToVertex.resize(tCellIndex, tEnrInterpMesh->mNumVertsPerInterpCell);
    tEnrInterpMesh->mEnrichedInterpCells.resize(tCellIndex);

    // with the cell to vertex data fully setup, add the vertex pointers to the cell

    for (moris::uint iE = 0; iE < tEnrInterpCellToVertex.n_rows(); iE++)
    {
        moris::Cell<Interpolation_Vertex_Unzipped*> tVertices(tEnrInterpCellToVertex.n_cols());

        // iterate and get vertices
        for (moris::uint iV = 0; iV < tEnrInterpCellToVertex.n_cols(); iV++)
        {
            tVertices(iV) = tEnrInterpMesh->get_unzipped_vertex_pointer(tEnrInterpCellToVertex(iE, iV));
        }

        // set vertices in cell
        tEnrInterpMesh->mEnrichedInterpCells(iE)->set_vertices(tVertices);
    }

    tEnrInterpMesh->mEnrichedInterpVerts.resize(tVertexCount);
}
//-------------------------------------------------------------------------------------

void
Enrichment::construct_enriched_vertex_interpolation(
    moris_index const&                      aEnrichmentDataIndex,
    mtk::Vertex_Interpolation*              aBaseVertexInterp,
    Cell<moris_index> const&                aSubPhaseBasisEnrLev,
    std::unordered_map<moris_id, moris_id>& aMapBasisIndexToLocInSubPhase,
    Vertex_Enrichment&                      aVertexEnrichment)
{
    // allocate a new vertex enrichment
    aVertexEnrichment = Vertex_Enrichment();

    // a nullptr here would indicate an aura node without a t-matrix
    if (aBaseVertexInterp != nullptr)
    {
        // iterate through basis in the base vertex interpolation
        moris::uint tNumCoeffs = aBaseVertexInterp->get_number_of_coefficients();

        // indices of the coefficients
        Matrix<IndexMat> tBaseVertCoeffInds = aBaseVertexInterp->get_indices();

        // get owners
        Matrix<IndexMat> tBaseVertOwners = aBaseVertexInterp->get_owners();

        // weights of the coefficients
        const Matrix<DDRMat>* tBaseVertWeights = aBaseVertexInterp->get_weights();

        // enriched Basis Coefficient indices
        Matrix<IndexMat> tEnrichCoeffInds(tBaseVertCoeffInds.n_rows(), tBaseVertCoeffInds.n_cols());
        Matrix<IndexMat> tEnrichCoeffIds(tBaseVertCoeffInds.n_rows(), tBaseVertCoeffInds.n_cols());

        for (moris::uint iBC = 0; iBC < tNumCoeffs; iBC++)
        {
            // coefficient of the base index
            moris_index tCoeffIndex = tBaseVertCoeffInds(iBC);

            // find the coefficients index within this subphase cluster
            auto tIter = aMapBasisIndexToLocInSubPhase.find(tCoeffIndex);

            MORIS_ASSERT(tIter != aMapBasisIndexToLocInSubPhase.end(), "Basis not found in vertex map");

            // The basis local index relative to the subphase
            moris_index tSubphaseBasisIndex = tIter->second;

            // enrichment level of the basis
            moris_index tEnrLev = aSubPhaseBasisEnrLev(tSubphaseBasisIndex);

            moris_index tEnrichedCoeffIndex = mEnrichmentData(aEnrichmentDataIndex).mBasisEnrichmentIndices(tCoeffIndex)(tEnrLev);

            tEnrichCoeffInds(iBC) = tEnrichedCoeffIndex;

            tEnrichCoeffIds(iBC) = mEnrichmentData(aEnrichmentDataIndex).mEnrichedBasisIndexToId(tEnrichedCoeffIndex);
        }

        std::unordered_map<moris::moris_index, moris::moris_index>& tVertEnrichMap = aVertexEnrichment.get_basis_map();

        for (moris::uint iB = 0; iB < tEnrichCoeffInds.numel(); iB++)
        {
            moris::moris_index tBasisIndex = tEnrichCoeffInds(iB);

            tVertEnrichMap[tBasisIndex] = iB;
        }

        aVertexEnrichment.add_basis_information(tEnrichCoeffInds, tEnrichCoeffIds);
        aVertexEnrichment.add_basis_owners(tEnrichCoeffInds, tBaseVertOwners);
        aVertexEnrichment.add_basis_weights(tEnrichCoeffInds, *tBaseVertWeights);
        aVertexEnrichment.add_base_vertex_interpolation(aBaseVertexInterp);
    }
}

//-------------------------------------------------------------------------------------

std::unordered_map<moris_id, moris_id>
Enrichment::construct_subphase_basis_to_basis_map(Cell<moris_id> const& aSubPhaseBasisIndex)
{
    uint tNumBasisOfSubphase = aSubPhaseBasisIndex.size();

    std::unordered_map<moris_id, moris_id> tSubphaseBasisMap;

    for (moris::uint iB = 0; iB < tNumBasisOfSubphase; iB++)
    {
        tSubphaseBasisMap[aSubPhaseBasisIndex(iB)] = iB;
    }

    return tSubphaseBasisMap;
}

//-------------------------------------------------------------------------------------

moris::Cell<mtk::Vertex_Interpolation*>
Enrichment::get_vertex_interpolations(
    moris::mtk::Cell& aParentCell,
    const uint        aMeshIndex) const
{
    uint tNumVerts = aParentCell.get_number_of_vertices();

    moris::Cell<mtk::Vertex*> tVertexPointers = aParentCell.get_vertex_pointers();

    moris::Cell<mtk::Vertex_Interpolation*> tVertexInterp(tNumVerts);

    for (moris::uint i = 0; i < tNumVerts; i++)
    {
        moris_index tVertexIndex = tVertexPointers(i)->get_index();

        tVertexInterp(i) = mEnrichmentData(aMeshIndex).mBGVertexInterpolations(tVertexIndex);
    }

    return tVertexInterp;
}

//-------------------------------------------------------------------------------------

Cell<std::string>
Enrichment::get_cell_enrichment_field_names() const
{
    // number of basis types
    moris::uint tNumBasisTypes = mEnrichmentData.size();

    // declare  enrichment fields
    Cell<std::string> tEnrichmentFields;

    for (moris::uint iBT = 0; iBT < tNumBasisTypes; iBT++)
    {
        // number of basis
        moris::size_t tNumBasis = mBackgroundMeshPtr->get_num_basis_functions(mMeshIndices(iBT));

        std::string tBaseEnrich = "el_bt_" + std::to_string(mMeshIndices(iBT)) + "b_";
        for (size_t i = 0; i < tNumBasis; i++)
        {
            tEnrichmentFields.push_back(
                tBaseEnrich + std::to_string(mBackgroundMeshPtr->get_glb_entity_id_from_entity_loc_index(i, mBasisRank, mMeshIndices(iBT))));
        }
    }

    // Add local floodfill field to the output mesh
    std::string tLocalFFStr = "child_ff";
    tEnrichmentFields.push_back(tLocalFFStr);

    std::string tSubPhaseStr = "subphase";
    tEnrichmentFields.push_back(tSubPhaseStr);

    return tEnrichmentFields;
}

//-------------------------------------------------------------------------------------

void
Enrichment::write_cell_enrichment_to_fields(
    Cell<std::string>& aEnrichmentFieldStrs,
    mtk::Mesh*         aMeshWithEnrFields) const
{
    MORIS_ERROR(0, "Deprecated. (Removal in progress)");
}

//-------------------------------------------------------------------------------------

moris_index
Enrichment::allocate_basis_ids(moris_index const& aMeshIndex,
    moris_index const&                            aNumIdsToAllocate)
{
    int tProcRank = par_rank();
    int tProcSize = par_size();

    // size_t is defined as uint here because of aNumRequested
    //Initialize gathered information outputs (information which will be scattered across processors)
    moris::Cell<moris::moris_id> aGatheredInfo;
    moris::Cell<moris::moris_id> tFirstId(1);
    moris::Cell<moris::moris_id> tNumIdsRequested(1);

    tNumIdsRequested(0) = (moris::moris_id)aNumIdsToAllocate;

    moris::gather(tNumIdsRequested, aGatheredInfo);

    moris::Cell<moris::moris_id> tProcFirstID(tProcSize);

    // local maximum
    moris_index tLocMaxId = this->get_max_basis_id(aMeshIndex);

    // global maximum
    moris_index tFirstAvailGlbMaxId = moris::max_all(tLocMaxId) + 1;
    if (tProcRank == 0)
    {
        // Loop over entities print the number of entities requested by each processor
        for (int iProc = 0; iProc < tProcSize; ++iProc)
        {
            // Give each processor their desired amount of IDs
            tProcFirstID(iProc) = tFirstAvailGlbMaxId;

            // Increment the first available node ID
            tFirstAvailGlbMaxId = tFirstAvailGlbMaxId + aGatheredInfo(iProc);
        }
    }

    moris::scatter(tProcFirstID, tFirstId);


    return tFirstId(0);
}

moris_index
Enrichment::get_max_basis_id(moris_index const& aMeshIndex)
{
    // Number of basis functions
    moris::size_t tNumBasis = mXTKModelPtr->get_background_mesh().get_num_basis_functions(aMeshIndex);

    //maximum id
    moris_id tMaxId = 0;

    for (moris::uint i = 0; i < tNumBasis; i++)
    {
        // get the basis id
        moris_id tBasisId = mBackgroundMeshPtr->get_glb_entity_id_from_entity_loc_index(i, mBasisRank, aMeshIndex);

        if (tMaxId < tBasisId)
        {
            tMaxId = tBasisId;
        }
    }
    return tMaxId;
}
}// namespace xtk
