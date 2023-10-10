/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Enrichment.cpp
 *
 */

#include "cl_XTK_Enrichment.hpp"

#include <iomanip>
#include <iostream>
#include <string>
#include <set>

#include "cl_Communication_Tools.hpp"
#include "linalg_typedefs.hpp"
#include "fn_assert.hpp"
#include "fn_sort.hpp"
#include "typedefs.hpp"
#include "fn_stringify_matrix.hpp"
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

#include "fn_sort_points_by_coordinates.hpp"
#include "fn_XTK_convert_cell_to_map.hpp"
#include "fn_XTK_punch_card_sum.hpp"
#include "fn_XTK_find_in_cell.hpp"

#include "cl_Tracer.hpp"

using namespace moris;

namespace xtk
{
    //-------------------------------------------------------------------------------------

    Enrichment::Enrichment(
            enum Enrichment_Method const & aMethod,
            mtk::EntityRank const &        aBasisRank,
            Matrix< IndexMat > const &     aInterpIndex,
            moris_index const &            aNumBulkPhases,
            xtk::Model*                    aXTKModelPtr,
            moris::mtk::Mesh*              aBackgroundMeshPtr,
            bool                           aSortBasisEnrichmentLevels,
            bool                           aUseSpgBasedEnrichment )
            : mEnrichmentMethod( aMethod )
            , mBasisRank( aBasisRank )
            , mMeshIndices( aInterpIndex )
            , mNumBulkPhases( aNumBulkPhases )
            , mXTKModelPtr( aXTKModelPtr )
            , mBackgroundMeshPtr( aBackgroundMeshPtr )
            , mEnrichmentData( aInterpIndex.max() + 1, mXTKModelPtr->get_cut_integration_mesh()->get_num_subphases() )
            , mSortBasisEnrichmentLevels( aSortBasisEnrichmentLevels )
    {
        mCutIgMesh        = mXTKModelPtr->get_cut_integration_mesh();
        mIgMeshTools      = new Integration_Mesh_Generator();
        mBsplineMeshInfos = mCutIgMesh->get_bspline_mesh_info();

        // FIXME: this needs to go once the SPG based enrichment is validated
        // reinitialize the enrichment data for SPGs if needed
        if ( aUseSpgBasedEnrichment )
        {
            // initialize the enrichment data for every discretization mesh index, as each mesh will result in a different number of SPGs
            for ( uint iMeshIndex = 0; iMeshIndex < aInterpIndex.numel(); iMeshIndex++ )
            {
                moris_index tMeshIndex = aInterpIndex( iMeshIndex );
                uint        tNumSPGs   = aXTKModelPtr->get_cut_integration_mesh()->get_num_subphase_groups( iMeshIndex );
                mEnrichmentData( tMeshIndex ).reinitialize_for_SPGs( tNumSPGs );
            }
        }
    }

    //-------------------------------------------------------------------------------------

    Enrichment::~Enrichment()
    {
        delete mIgMeshTools;
    }

    //-------------------------------------------------------------------------------------

    void
    Enrichment::perform_enrichment()
    {
        MORIS_ERROR( mBackgroundMeshPtr != nullptr,
                "mBackgroundMesh nullptr detected, this is probably because the enrichment has not been initialized properly" );

        // Perform enrichment over basis clusters
        perform_basis_cluster_enrichment();
    }

    //-------------------------------------------------------------------------------------

    void
    Enrichment::perform_enrichment_new()
    {
        MORIS_ERROR( mBackgroundMeshPtr != nullptr,
                "mBackgroundMesh nullptr detected, this is probably because the enrichment has not been initialized properly" );

        // Perform enrichment over basis clusters
        perform_basis_cluster_enrichment_new();
    }

    //-------------------------------------------------------------------------------------

    Cell< Matrix< IdMat > > const &
    Enrichment::get_element_inds_in_basis_support( moris_index const & aEnrichmentDataIndex ) const
    {
        return mEnrichmentData( aEnrichmentDataIndex ).mElementIndsInBasis;
    }

    //-------------------------------------------------------------------------------------

    Cell< Matrix< IndexMat > > const &
    Enrichment::get_element_enrichment_levels_in_basis_support( moris_index const & aEnrichmentDataIndex ) const
    {
        return mEnrichmentData( aEnrichmentDataIndex ).mElementEnrichmentLevel;
    }

    //-------------------------------------------------------------------------------------

    Cell< Matrix< IndexMat > > const &
    Enrichment::get_subphases_loc_inds_in_enriched_basis( moris_index const & aEnrichmentDataIndex ) const
    {
        return mEnrichmentData( aEnrichmentDataIndex ).mSubphaseIndsInEnrichedBasis;
    }

    //-------------------------------------------------------------------------------------

    moris::Memory_Map
    Enrichment::get_memory_usage()
    {
        Memory_Map tMemoryMap;

        tMemoryMap.mMemoryMapData[ "mEnrichmentMethod" ] = sizeof( mEnrichmentMethod );
        tMemoryMap.mMemoryMapData[ "mBasisRank" ]        = sizeof( mBasisRank );
        tMemoryMap.mMemoryMapData[ "mMeshIndices" ]      = mMeshIndices.capacity();
        tMemoryMap.mMemoryMapData[ "mNumBulkPhases" ]    = sizeof( mNumBulkPhases );
        tMemoryMap.mMemoryMapData[ "mMeshIndices" ]      = mMeshIndices.capacity();

        tMemoryMap.mMemoryMapData[ "mElementEnrichmentLevel" ]      = moris::internal_capacity( mEnrichmentData( 0 ).mElementEnrichmentLevel );
        tMemoryMap.mMemoryMapData[ "mElementIndsInBasis" ]          = moris::internal_capacity( mEnrichmentData( 0 ).mElementIndsInBasis );
        tMemoryMap.mMemoryMapData[ "mSubphaseIndsInEnrichedBasis" ] = moris::internal_capacity( mEnrichmentData( 0 ).mSubphaseIndsInEnrichedBasis );
        tMemoryMap.mMemoryMapData[ "mBasisEnrichmentIndices" ]      = moris::internal_capacity( mEnrichmentData( 0 ).mBasisEnrichmentIndices );
        tMemoryMap.mMemoryMapData[ "mEnrichedBasisIndexToId" ]      = mEnrichmentData( 0 ).mEnrichedBasisIndexToId.capacity();
        tMemoryMap.mMemoryMapData[ "mSubphaseBGBasisIndices" ]      = moris::internal_capacity( mEnrichmentData( 0 ).mSubphaseBGBasisIndices );
        tMemoryMap.mMemoryMapData[ "mSubphaseBGBasisEnrLev" ]       = moris::internal_capacity( mEnrichmentData( 0 ).mSubphaseBGBasisEnrLev );
        tMemoryMap.mMemoryMapData[ "mNumEnrichedBasisFunctions" ]   = sizeof( mEnrichmentData( 0 ).mNumEnrichedBasisFunctions );
        tMemoryMap.mMemoryMapData[ "mBGVertexInterpolations ptrs" ] = mEnrichmentData( 0 ).mBGVertexInterpolations.capacity();

        return tMemoryMap;
    }

    //-------------------------------------------------------------------------------------

    void
    Enrichment::write_diagnostics()
    {
        for ( moris::size_t iBasisType = 0; iBasisType < mMeshIndices.numel(); iBasisType++ )
        {
            // get the mesh index
            moris_index tMeshIndex = mMeshIndices( iBasisType );

            std::string tEnrBasisIdToSubphaseId = mXTKModelPtr->get_diagnostic_file_name( std::string( "Enr_Basis_To_Subphase_" + std::to_string( tMeshIndex ) ) );

            this->print_enriched_basis_to_subphase_id( tMeshIndex, tEnrBasisIdToSubphaseId );
        }

        mXTKModelPtr->mEnrichedInterpMesh( 0 )->write_diagnostics();
    }

    //-------------------------------------------------------------------------------------

    moris_index
    Enrichment::get_list_index_for_mesh_index( moris_index aMeshIndex )
    {
        // compute the mesh index to List index map if it has not been set yet
        if ( !mMeshIndexMapIsSet )
        {
            // loop over the mesh indices and fill the map
            for ( uint iMeshListIndex = 0; iMeshListIndex < mMeshIndices.numel(); iMeshListIndex++ )
            {
                mMeshIndexToListIndexMap[ mMeshIndices( iMeshListIndex ) ] = iMeshListIndex;
            }
        }

        // find the mesh index in the map
        auto tIter = mMeshIndexToListIndexMap.find( aMeshIndex );

        // return its position in list
        return tIter->second;
    }

    //-------------------------------------------------------------------------------------

    void
    Enrichment::print_enriched_basis_to_subphase_id(
            const moris_index& aMeshIndex,
            std::string        aFileName )
    {
        Cell< Matrix< IndexMat > > const & tSubphasesInEnrBasis = mEnrichmentData( aMeshIndex ).mSubphaseIndsInEnrichedBasis;

        std::ostringstream tStringStream;
        tStringStream << "Enriched_Basis_Id,";
        tStringStream << "Mesh_Index,";
        tStringStream << "Owner,";
        tStringStream << "PRank,";

        // global max size of
        moris_index tLocalMaxIGCellGroupSize = 0;
        for ( uint iEB = 0; iEB < tSubphasesInEnrBasis.size(); iEB++ )
        {
            if ( (moris_index)tSubphasesInEnrBasis( iEB ).numel() > tLocalMaxIGCellGroupSize )
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

        Matrix< IndexMat > const & tLocToGlbEnrBasisId = mXTKModelPtr->mEnrichedInterpMesh( 0 )->get_enriched_coefficient_local_to_global_map( aMeshIndex );

        for ( uint iEB = 0; iEB < tSubphasesInEnrBasis.size(); iEB++ )
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
        // make sure mesh is not empty
        MORIS_ASSERT( mBackgroundMeshPtr->get_num_elems() > 0, "Enrichment::perform_basis_cluster_enrichment() - IP mesh without cells passed" );

        // construct cell in xtk conformal model neighborhood connectivity
        // todo: this does nothing right now, still needed?
        this->construct_neighborhoods();

        // make sure we have access to all the vertex interpolation
        this->setup_background_vertex_interpolations();

        // construct data needed for enrichment for every B-spline mesh the Lagrange mesh is related to
        for ( moris::size_t iMeshIndex = 0; iMeshIndex < mMeshIndices.numel(); iMeshIndex++ )
        {
            // get the mesh index
            moris_index tMeshIndex = mMeshIndices( iMeshIndex );

            // log/trace enrichment for every discretization mesh index
            Tracer tTracer( "XTK", "Enrichment", "Mesh Index " + std::to_string( tMeshIndex ) );

            // Number of basis functions (= number of B-Splines on the A-mesh ?)
            moris::size_t tNumBasisFunctions = mBackgroundMeshPtr->get_num_basis_functions( tMeshIndex );

            // allocate member variables
            mEnrichmentData( tMeshIndex ).mElementEnrichmentLevel = Cell< Matrix< IndexMat > >( tNumBasisFunctions );
            mEnrichmentData( tMeshIndex ).mElementIndsInBasis     = Cell< Matrix< IndexMat > >( tNumBasisFunctions );

            // allocate data used after basis loop
            Cell< Matrix< IndexMat > > tSubPhaseBinEnrichment( tNumBasisFunctions );
            Cell< Matrix< IndexMat > > tSubphaseClusterIndicesInSupport( tNumBasisFunctions );
            Cell< moris_index >        tMaxEnrichmentLevel( tNumBasisFunctions, 0 );

            for ( moris::size_t iBasisFunction = 0; iBasisFunction < tNumBasisFunctions; iBasisFunction++ )
            {
                // Get elements in support of basis (these are interpolation cells)
                Matrix< IndexMat > tParentElementsInSupport;

                mBackgroundMeshPtr->get_elements_in_support_of_basis( tMeshIndex, iBasisFunction, tParentElementsInSupport );

                // get subphase clusters in support (separated by phase)
                tSubphaseClusterIndicesInSupport( iBasisFunction ) = this->get_subphase_clusters_in_support( tParentElementsInSupport );

                // construct subphase in support map
                IndexMap tSubPhaseIndexToSupportIndex;

                this->construct_subphase_in_support_map( tSubphaseClusterIndicesInSupport( iBasisFunction ), tSubPhaseIndexToSupportIndex );

                // prune the subphase to remove subphases outside of basis support
                Matrix< IndexMat > tPrunedSubphaseNeighborhood;

                this->generate_pruned_subphase_graph_in_basis_support(
                        tSubphaseClusterIndicesInSupport( iBasisFunction ),
                        tSubPhaseIndexToSupportIndex,
                        tPrunedSubphaseNeighborhood );

                // Assign enrichment levels to subphases
                this->assign_subphase_bin_enrichment_levels_in_basis_support(
                        tSubphaseClusterIndicesInSupport( iBasisFunction ),
                        tSubPhaseIndexToSupportIndex,
                        tPrunedSubphaseNeighborhood,
                        tSubPhaseBinEnrichment( iBasisFunction ),
                        tMaxEnrichmentLevel( iBasisFunction ) );

                // Sort enrichment levels
                if ( mSortBasisEnrichmentLevels )
                {
                    this->sort_enrichment_levels_in_basis_support(
                            tSubphaseClusterIndicesInSupport( iBasisFunction ),
                            tSubPhaseBinEnrichment( iBasisFunction ),
                            tMaxEnrichmentLevel( iBasisFunction ) );
                }

                // Extract element enrichment levels from assigned sub-phase bin enrichment levels and store these as a member variable
                this->unzip_subphase_bin_enrichment_into_element_enrichment(
                        tMeshIndex,
                        iBasisFunction,
                        tParentElementsInSupport,
                        tSubphaseClusterIndicesInSupport( iBasisFunction ),
                        tSubPhaseIndexToSupportIndex,
                        tPrunedSubphaseNeighborhood,
                        tSubPhaseBinEnrichment( iBasisFunction ) );
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

            MORIS_LOG_SPEC( "Num Non-enriched Bases", mEnrichmentData( tMeshIndex ).mElementIndsInBasis.size() );
            MORIS_LOG_SPEC( "Num Enriched Bases", mEnrichmentData( tMeshIndex ).mNumEnrichedBasisFunctions );
        }

        // create the enriched interpolation mesh
        this->construct_enriched_interpolation_mesh();

        // create the integration mesh
        this->construct_enriched_integration_mesh();
    }

    //-------------------------------------------------------------------------------------

    void
    Enrichment::perform_basis_cluster_enrichment_new()
    {
        // make sure mesh is not empty
        MORIS_ASSERT( mBackgroundMeshPtr->get_num_elems() > 0, "Enrichment::perform_basis_cluster_enrichment_new() - 0 cell interpolation mesh passed" );

        // make sure we have access to all the vertex interpolation
        this->setup_background_vertex_interpolations();

        // iterate through B-spline meshes
        for ( moris::size_t iMeshIndex = 0; iMeshIndex < mMeshIndices.numel(); iMeshIndex++ )
        {
            // get the mesh index
            moris_index tMeshIndex = mMeshIndices( iMeshIndex );

            Tracer tTracer( "XTK", "Enrichment", "Mesh Index " + std::to_string( tMeshIndex ) );

            // Number of basis functions (= number of B-Splines on the A-mesh )
            moris::size_t tNumBasisFunctions = mBackgroundMeshPtr->get_num_basis_functions( tMeshIndex );

            // allocate member variables
            mEnrichmentData( tMeshIndex ).mElementEnrichmentLevel = Cell< Matrix< IndexMat > >( tNumBasisFunctions );
            mEnrichmentData( tMeshIndex ).mElementIndsInBasis     = Cell< Matrix< IndexMat > >( tNumBasisFunctions );

            // allocate data used after basis loop
            Cell< Matrix< IndexMat > > tSpgBinEnrichment( tNumBasisFunctions );
            Cell< Matrix< IndexMat > > tSpgIndicesInSupport( tNumBasisFunctions );
            Cell< moris_index >        tMaxEnrichmentLevel( tNumBasisFunctions, 0 );

            for ( moris::size_t iBasisFunction = 0; iBasisFunction < tNumBasisFunctions; iBasisFunction++ )
            {
                // Get elements in support of basis (these are interpolation cells)
                Matrix< IndexMat > tParentElementsInSupport;

                mBackgroundMeshPtr->get_elements_in_support_of_basis( tMeshIndex, iBasisFunction, tParentElementsInSupport );

                // collect Subphase group indices in support
                IndexMap tSpgIndexToSupportIndex;
                this->get_subphase_groups_in_support(
                        iMeshIndex,
                        tParentElementsInSupport,
                        tSpgIndicesInSupport( iBasisFunction ),
                        tSpgIndexToSupportIndex );

                // prune the subphase to remove subphases outside of basis support
                Matrix< IndexMat > tPrunedSpgNeighborhood;
                this->generate_pruned_subphase_group_graph_in_basis_support(
                        iMeshIndex,
                        tSpgIndicesInSupport( iBasisFunction ),
                        tSpgIndexToSupportIndex,
                        tPrunedSpgNeighborhood );

                // Assign enrichment levels to subphases
                this->assign_subphase_group_bin_enrichment_levels_in_basis_support(
                        tSpgIndicesInSupport( iBasisFunction ),
                        tPrunedSpgNeighborhood,
                        tSpgBinEnrichment( iBasisFunction ),
                        tMaxEnrichmentLevel( iBasisFunction ) );

                // Sort enrichment levels
                if ( mSortBasisEnrichmentLevels )
                {
                    MORIS_ERROR( false, "Enrichment::perform_basis_cluster_enrichment_new() - function: sort_enrichment_levels_in_basis_support() not supported yet with new SPG based enrichment" );
                    // FIXME: this function doesn't work for SPGs yet
                    this->sort_enrichment_levels_in_basis_support(
                            tSpgIndicesInSupport( iBasisFunction ),
                            tSpgBinEnrichment( iBasisFunction ),
                            tMaxEnrichmentLevel( iBasisFunction ) );
                }

                // Extract element enrichment levels from assigned sub-phase bin enrichment levels and store these as a member variable
                this->unzip_subphase_group_bin_enrichment_into_element_enrichment(
                        tMeshIndex,
                        iBasisFunction,
                        tParentElementsInSupport,
                        tSpgIndicesInSupport( iBasisFunction ),
                        tSpgBinEnrichment( iBasisFunction ) );
            }

            // construct subphase to enriched index
            this->construct_enriched_basis_to_subphase_group_connectivity(
                    tMeshIndex,
                    tSpgBinEnrichment,
                    tSpgIndicesInSupport,
                    tMaxEnrichmentLevel );

            // Assign enriched basis indices Only indices here because interpolation cells needs basis
            // indices and basis ids are created using the interpolation cell ids
            this->assign_enriched_coefficients_identifiers_new(
                    tMeshIndex,
                    tMaxEnrichmentLevel );

            // log the total number of enriched Basis functions across all procs
            uint        tNumEnrBFs     = mEnrichmentData( tMeshIndex ).mNumEnrichedBasisFunctions;
            std::string tTempLogString = "Number of enriched basis functions on B-spline mesh #" + std::to_string( tMeshIndex );
            MORIS_LOG_SPEC( tTempLogString, sum_all( tNumEnrBFs ) );

            // log how big the function basis is for each material (this is useful information for the EXHUME project)

            // get the number of bulk phases in mesh
            size_t tNumBulkPhases = mXTKModelPtr->get_geom_engine()->get_num_bulk_phase();

            // initialize list that collects number of enriched basis functions in each bulk phase
            Cell< uint > tNumEnrBfsInBulkPhase( tNumBulkPhases, 0 );

            // go through enriched basis functions and count them towards the bulk phase they interpolate into
            for ( uint iEnrBF = 0; iEnrBF < tNumEnrBFs; iEnrBF++ )
            {
                // get the bulk phase for the current enriched basis function
                moris_index tBulkPhaseBfInterpolatesInto = mEnrichmentData( tMeshIndex ).mBulkPhaseInEnrichedBasis( iEnrBF );

                // count BF towards BFs on bulk-phase
                tNumEnrBfsInBulkPhase( tBulkPhaseBfInterpolatesInto )++;
            }

            // log/output how may enr. basis functions there are in each bulk phase
            for ( uint iBulkPhase = 0; iBulkPhase < tNumBulkPhases; iBulkPhase++ )
            {
                MORIS_LOG(
                        "B-spline mesh #%i, Bulk phase #%i, Number of enriched basis functions on material domain: %i",
                        tMeshIndex,
                        iBulkPhase,
                        sum_all( tNumEnrBfsInBulkPhase( iBulkPhase ) ) );
            }

            // find the SPs within one IP cell that are also within the same SPG
            // NOTE: this is needed for cluster construction later; also needs to be repeated for every Mesh index, hence it's here
            // this->establish_IP_SPG_SP_relationship( tMeshIndex );
        }

        // create the enriched interpolation mesh
        this->construct_enriched_interpolation_mesh_new();

        // FIXME: the IdMat below should eventually be created and passed as a shared-ptr to avoid the copy operation
        // pass Enr. BF to Bulkphase map to enriched IP mesh
        for ( uint iMeshIndex = 0; iMeshIndex < mMeshIndices.numel(); iMeshIndex++ )
        {
            moris_index tMeshIndex = mMeshIndices( iMeshIndex );

            Matrix< IdMat > tBulkPhaseInEnrichedBasis = mEnrichmentData( tMeshIndex ).mBulkPhaseInEnrichedBasis;
            mXTKModelPtr->mEnrichedInterpMesh( 0 )->set_enriched_basis_to_bulkphase_map( tMeshIndex, tBulkPhaseInEnrichedBasis );
        }

        // construct the enriched IG mesh (i.e. clusters for all B-spline meshes)
        this->construct_enriched_integration_mesh( mMeshIndices );
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
        mtk::Interpolation_Mesh& tBgIpMesh = mXTKModelPtr->get_background_mesh();

        // size the data in the enrichment data
        uint tNumVertices = tBgIpMesh.get_num_nodes();

        // iterate through meshes
        for ( uint iM = 0; iM < mMeshIndices.numel(); iM++ )
        {
            // get the mesh index
            moris_index tMeshIndex = mMeshIndices( iM );

            // size the vertex interpolations (allocating everything as null)
            mEnrichmentData( tMeshIndex ).mBGVertexInterpolations.resize( tNumVertices, nullptr );

            // iterate through vertices
            for ( uint iV = 0; iV < tNumVertices; iV++ )
            {
                // access the vertex
                moris::mtk::Vertex& tVert = tBgIpMesh.get_mtk_vertex( (moris_index)iV );

                // if the vertex interpolation exits this vertex is own by the current processor
                // and the vertex interpolation is stored; otherwise the nullptr is kept indicating
                // that the vertex is not owned the by current processor
                if ( tVert.has_interpolation( tMeshIndex ) )
                {
                    mEnrichmentData( tMeshIndex ).mBGVertexInterpolations( iV ) = tVert.get_interpolation( tMeshIndex );
                }
            }
        }
    }

    //-------------------------------------------------------------------------------------

    Matrix< IndexMat >
    Enrichment::get_subphase_clusters_in_support( Matrix< IndexMat > const & aElementsInSupport )
    {
        // count the number of subphase cluster in support
        uint tCount = 0;

        for ( moris::size_t iE = 0; iE < aElementsInSupport.numel(); iE++ )
        {
            Cell< moris_index > const & tSubphaseIndices = mCutIgMesh->get_parent_cell_subphases( aElementsInSupport( iE ) );

            tCount = tCount + tSubphaseIndices.size();
        }

        Matrix< IndexMat > tSubPhaseClusters( 1, tCount );
        tCount = 0;

        for ( moris::size_t iE = 0; iE < aElementsInSupport.numel(); iE++ )
        {

            Cell< moris_index > const & tSubphaseIndices = mCutIgMesh->get_parent_cell_subphases( aElementsInSupport( iE ) );

            for ( uint iSP = 0; iSP < tSubphaseIndices.size(); iSP++ )
            {
                tSubPhaseClusters( tCount++ ) = tSubphaseIndices( iSP );
            }
        }

        return tSubPhaseClusters;
    }

    //-------------------------------------------------------------------------------------

    void
    Enrichment::get_subphase_groups_in_support(
            moris_index                aMeshIndexPosition,
            Matrix< IndexMat > const & aLagElementsInSupport,
            Matrix< IndexMat >&        aSubphaseGroupIndicesInSupport,
            IndexMap&                  aSubphaseGroupIndexToSupportIndex )
    {
        // initialize map between Lagrange an Bspline elements
        IndexMap tBspElemIndexMap;

        // get number of lagrange elements in support
        uint tNumLagElemsInSupport = aLagElementsInSupport.numel();

        // initialize counter tracking the number of B-spline elements in support
        uint tBspElemCounter = 0;

        // go through Lag elements and find their corresponding B-Spline elements
        for ( uint iLagElem = 0; iLagElem < tNumLagElemsInSupport; iLagElem++ )
        {
            // temporarily store current Lagrange element's index
            moris_index tLagElemIndex = aLagElementsInSupport( iLagElem );

            // get the corresponding B-spline element index
            moris_index tBsplineElemIndex = mBsplineMeshInfos( aMeshIndexPosition )->get_bspline_cell_index_for_extraction_cell( tLagElemIndex );

            // look for
            auto tBsplineMapPosition = tBspElemIndexMap.find( tBsplineElemIndex );

            // check if the B-spline element has already been found in basis support, if not ...
            if ( tBsplineMapPosition == tBspElemIndexMap.end() )
            {
                // ... register the B-spline element in the map
                tBspElemIndexMap[ tBsplineElemIndex ] = tBspElemCounter;

                // increment the B-spline element counter
                tBspElemCounter++;
            }
        }

        // copy the found elements into an array
        Cell< moris_index > tBspElemIndicesInSupport( tBspElemCounter );
        for ( uint iLagElem = 0; iLagElem < tNumLagElemsInSupport; iLagElem++ )
        {
            // get the Lagrange element's index
            moris_index tLagElemIndex = aLagElementsInSupport( iLagElem );

            // get the corresponding B-spline element index
            moris_index tBsplineElemIndex = mBsplineMeshInfos( aMeshIndexPosition )->get_bspline_cell_index_for_extraction_cell( tLagElemIndex );

            // look for where in the local list of B-spline elements this element is
            moris_index tPos = tBspElemIndexMap[ tBsplineElemIndex ];

            // copy and store
            tBspElemIndicesInSupport( tPos ) = tBsplineElemIndex;
        }

        // initialize counter for subphase groups
        uint tSpgCounter = 0;

        // for each B-spline element get the SPGs, count them, and store the local ordering
        for ( auto& iBspElemIndex : tBspElemIndicesInSupport )
        {
            // get the SPGs associated with the current B-spline element
            Cell< moris_index > const & tSPGsInBsplineElem =
                    mBsplineMeshInfos( aMeshIndexPosition )->get_SPG_indices_in_bspline_cell( iBspElemIndex );

            // get the number of SPGs in the current B-spline element
            uint tNumSpgsInBspElem = tSPGsInBsplineElem.size();

            // store all SPGs in map
            for ( uint iSPG = 0; iSPG < tNumSpgsInBspElem; iSPG++ )
            {
                // get the locally
                moris_index tSpgIndex = tSPGsInBsplineElem( iSPG );

                // save SPG index in map
                aSubphaseGroupIndexToSupportIndex[ tSpgIndex ] = tSpgCounter;

                // count SPGs in support
                tSpgCounter++;
            }
        }

        // initialize linear list storing all SPGs in support with correct size
        aSubphaseGroupIndicesInSupport.resize( tSpgCounter, 1 );

        // store all SPGs found in support in linear list
        for ( auto& iSPG : aSubphaseGroupIndexToSupportIndex )
        {
            // temporarily store current SPG's index
            moris_index tSpgIndex      = iSPG.first;
            moris_index tLocalSpgIndex = iSPG.second;

            // save SPG index in list
            aSubphaseGroupIndicesInSupport( tLocalSpgIndex ) = tSpgIndex;
        }
    }

    //-------------------------------------------------------------------------------------

    void
    Enrichment::construct_subphase_in_support_map(
            Matrix< IndexMat > const & aSubphaseClusterIndicesInSupport,
            IndexMap&                  aSubPhaseIndexToSupportIndex )
    {
        for ( moris_index i = 0; i < (moris_index)aSubphaseClusterIndicesInSupport.numel(); i++ )
        {
            aSubPhaseIndexToSupportIndex[ aSubphaseClusterIndicesInSupport( i ) ] = i;
        }
    }

    //-------------------------------------------------------------------------------------

    void
    Enrichment::generate_pruned_subphase_graph_in_basis_support(
            Matrix< IndexMat > const & aSubphasesInSupport,
            IndexMap&                  aSubPhaseIndexToSupportIndex,
            Matrix< IndexMat >&        aPrunedSubPhaseToSubphase )
    {
        std::shared_ptr< Subphase_Neighborhood_Connectivity > tSubphaseNeighborhood = mCutIgMesh->get_subphase_neighborhood();

        // get the number of subphases in the current BF's support
        uint tNumSPsInSupport = aSubphasesInSupport.numel();

        // Construct full element neighbor graph in support and the corresponding shared faces
        aPrunedSubPhaseToSubphase.resize( tNumSPsInSupport, 50 );    // FIXME: Keenan this allocation needs to done smarter
        aPrunedSubPhaseToSubphase.fill( MORIS_INDEX_MAX );

        // get subphase neighborhood information
        Cell< std::shared_ptr< Cell< moris_index > > > const & tSubphaseToSubphase = tSubphaseNeighborhood->mSubphaseToSubPhase;

        for ( moris::size_t iSP = 0; iSP < tNumSPsInSupport; iSP++ )
        {
            Cell< moris_index > const & tSingleSubPhaseNeighbors = *tSubphaseToSubphase( aSubphasesInSupport( iSP ) );

            // iterate through and prune subphases not in support
            uint tCount = 0;
            for ( moris::size_t iNeighborSP = 0; iNeighborSP < tSingleSubPhaseNeighbors.size(); iNeighborSP++ )
            {
                moris_index tNeighborSubphaseIndex = tSingleSubPhaseNeighbors( iNeighborSP );

                auto tNeighborIter = aSubPhaseIndexToSupportIndex.find( tNeighborSubphaseIndex );

                if ( tNeighborIter != aSubPhaseIndexToSupportIndex.end() )
                {
                    aPrunedSubPhaseToSubphase( iSP, tCount ) = tNeighborIter->second;
                    tCount++;
                }
            }
        }
    }

    //-------------------------------------------------------------------------------------

    void
    Enrichment::generate_pruned_subphase_group_graph_in_basis_support(
            const moris_index          aMeshIndex,
            Matrix< IndexMat > const & aSubphaseGroupIndicesInSupport,
            IndexMap&                  aSubphaseGroupIndexToSupportIndex,
            Matrix< IndexMat >&        aPrunedSpgToSpg )
    {
        // get the global SPG connectivity graph
        std::shared_ptr< Subphase_Neighborhood_Connectivity > tSubphaseGroupNeighborhood = mCutIgMesh->get_subphase_group_neighborhood( aMeshIndex );

        // Construct full SPG neighbor graph in support and the corresponding shared faces
        aPrunedSpgToSpg.resize( aSubphaseGroupIndicesInSupport.numel(), 50 );    // FIXME: this allocation needs to be done smarter
        aPrunedSpgToSpg.fill( MORIS_INDEX_MAX );

        // get subphase group neighborhood information
        Cell< std::shared_ptr< Cell< moris_index > > > const & tSpgToSpg = tSubphaseGroupNeighborhood->mSubphaseToSubPhase;

        for ( uint iSPG = 0; iSPG < aSubphaseGroupIndicesInSupport.numel(); iSPG++ )
        {
            // get the index of the current SPG
            moris_index tCenterSpgIndex = aSubphaseGroupIndicesInSupport( iSPG );

            // get temporarily list of Spg Neighbors
            Cell< moris_index > const & tSubphaseGroupNeighbors = *tSpgToSpg( tCenterSpgIndex );

            // iterate through and prune/remove subphases groups not in support
            uint tCount = 0;
            for ( uint iNeighborSpg = 0; iNeighborSpg < tSubphaseGroupNeighbors.size(); iNeighborSpg++ )
            {
                moris_index tNeighborSubphaseGroupIndex = tSubphaseGroupNeighbors( iNeighborSpg );

                auto tNeighborIter = aSubphaseGroupIndexToSupportIndex.find( tNeighborSubphaseGroupIndex );

                if ( tNeighborIter != aSubphaseGroupIndexToSupportIndex.end() )
                {
                    aPrunedSpgToSpg( iSPG, tCount ) = tNeighborIter->second;
                    tCount++;
                }
            }
        }
    }

    //-------------------------------------------------------------------------------------

    void
    Enrichment::assign_subphase_bin_enrichment_levels_in_basis_support(
            Matrix< IndexMat > const & aSubphasesInSupport,
            IndexMap&                  aSubPhaseIndexToSupportIndex,
            Matrix< IndexMat > const & aPrunedSubPhaseToSubphase,
            Matrix< IndexMat >&        aSubPhaseBinEnrichmentVals,
            moris_index&               aMaxEnrichmentLevel )
    {
        // Variables needed for flood-fill, consider removing these.
        // Active bins to include in flood-fill (We include all bins)
        Matrix< IndexMat > tActiveBins( 1, aPrunedSubPhaseToSubphase.n_rows() );

        for ( moris::size_t i = 0; i < aPrunedSubPhaseToSubphase.n_rows(); i++ )
        {
            ( tActiveBins )( 0, i ) = i;
        }

        // Mark all as included
        Matrix< IndexMat > tIncludedBins( 1, aSubphasesInSupport.numel(), 1 );

        // Flood fill metric value (since all the subphases do not connect to dissimilar phases)
        Matrix< IndexMat > tDummyPhase( 1, aSubphasesInSupport.numel(), 1 );

        aSubPhaseBinEnrichmentVals = flood_fill(
                aPrunedSubPhaseToSubphase,
                tDummyPhase,
                tActiveBins,
                tIncludedBins,
                mNumBulkPhases,
                MORIS_INDEX_MAX,
                aMaxEnrichmentLevel,
                true );
    }

    //-------------------------------------------------------------------------------------

    void
    Enrichment::assign_subphase_group_bin_enrichment_levels_in_basis_support(
            Matrix< IndexMat > const & aSpgsInSupport,
            Matrix< IndexMat > const & aPrunedSpgToSpg,
            Matrix< IndexMat >&        aSpgBinEnrichmentVals,
            moris_index&               aMaxEnrichmentLevel )
    {
        // Variables needed for flood-fill, consider removing these.
        // Active bins to include in flood-fill (We include all bins)
        Matrix< IndexMat > tActiveBins( 1, aPrunedSpgToSpg.n_rows() );

        for ( moris::size_t i = 0; i < aPrunedSpgToSpg.n_rows(); i++ )
        {
            ( tActiveBins )( 0, i ) = i;
        }

        // Mark all as included
        Matrix< IndexMat > tIncludedBins( 1, aSpgsInSupport.numel(), 1 );

        // Flood fill metric value (since all the subphases do not connect to dissimilar phases)
        Matrix< IndexMat > tDummyPhase( 1, aSpgsInSupport.numel(), 1 );

        aSpgBinEnrichmentVals = flood_fill(
                aPrunedSpgToSpg,
                tDummyPhase,
                tActiveBins,
                tIncludedBins,
                mNumBulkPhases,
                MORIS_INDEX_MAX,
                aMaxEnrichmentLevel,
                true );
    }

    //-------------------------------------------------------------------------------------

    void
    Enrichment::sort_enrichment_levels_in_basis_support(
            Matrix< IndexMat > const & aSubphasesInSupport,
            Matrix< IndexMat >&        aSubPhaseBinEnrichmentVals,
            moris_index const          aMaxEnrichmentLevel )
    {
        // number of enrichment levels
        uint tNumEnrichLevel = aMaxEnrichmentLevel + 1;

        // if only one or less enrichment levels used; no need to sort them
        if ( tNumEnrichLevel < 2 )
        {
            return;
        }

        // check that maximum enrichment level index is smaller than number of enrichment levels
        MORIS_ASSERT( (uint)aSubPhaseBinEnrichmentVals.max() < tNumEnrichLevel,
                "Enrichment::sort_enrichment_levels_in_basis_support - incorrect enrichment index." );

        // allocate list of centroids and volume for each enrichment level
        Cell< Matrix< DDRMat > > tCentroids( tNumEnrichLevel );
        Cell< moris::real >      tTotalVolumes( tNumEnrichLevel );

        for ( moris::size_t i = 0; i < aSubphasesInSupport.numel(); i++ )
        {
            // get enrichment level
            moris_index tEnrichLevel = aSubPhaseBinEnrichmentVals( i );

            // get subphase index
            moris_index tSubphaseIndex = aSubphasesInSupport( i );

            // get cell group in subphase
            std::shared_ptr< IG_Cell_Group > tSubphaseIgCells = mCutIgMesh->get_subphase_ig_cells( tSubphaseIndex );

            // compute volume of first cell in IG cell group
            moris::real tVolume = tSubphaseIgCells->mIgCellGroup( 0 )->compute_cell_measure();

            // initialize or add volume weighted centroid of first cell in IG cell group
            if ( tCentroids( tEnrichLevel ).numel() == 0 )
            {
                tCentroids( tEnrichLevel )    = tVolume * tSubphaseIgCells->mIgCellGroup( 0 )->compute_cell_centroid();
                tTotalVolumes( tEnrichLevel ) = tVolume;
            }
            else
            {
                tCentroids( tEnrichLevel ) += tVolume * tSubphaseIgCells->mIgCellGroup( 0 )->compute_cell_centroid();
                tTotalVolumes( tEnrichLevel ) += tVolume;
            }

            // iterate through remaining IG cells
            for ( uint iSPCell = 1; iSPCell < tSubphaseIgCells->mIgCellGroup.size(); iSPCell++ )
            {
                tVolume = tSubphaseIgCells->mIgCellGroup( iSPCell )->compute_cell_measure();

                tCentroids( tEnrichLevel ) += tVolume * tSubphaseIgCells->mIgCellGroup( iSPCell )->compute_cell_centroid();
                tTotalVolumes( tEnrichLevel ) += tVolume;
            }
        }

        for ( uint iL = 0; iL < (uint)tNumEnrichLevel; ++iL )
        {
            // check that centroid for current enrichment level has been computed
            MORIS_ERROR( tCentroids( iL ).numel() > 0,
                    "Enrichment::sort_enrichment_levels - enrichment levels not consecutive." );

            // compute centroid
            tCentroids( iL ) = 1.0 / ( tTotalVolumes( iL ) + MORIS_REAL_EPS ) * tCentroids( iL );
        }

        // sort centroids by coordinates
        Cell< moris_index > tSortingIndex = sort_points_by_coordinates( tCentroids );

        // build enrichment level old to new map
        Cell< moris_index > tEnrichmentMap( tNumEnrichLevel );

        for ( uint i = 0; i < (uint)tNumEnrichLevel; ++i )
        {
            tEnrichmentMap( tSortingIndex( i ) ) = i;
        }

        // resort enrichment levels of subphases
        Matrix< IndexMat > tSortedEnrichmentVals( 1, aSubphasesInSupport.numel(), MORIS_SINT_MAX );

        for ( uint i = 0; i < aSubphasesInSupport.numel(); ++i )
        {
            // old enrichment level
            moris_index tOldEnrichLevel = aSubPhaseBinEnrichmentVals( i );

            // new enrichment level
            tSortedEnrichmentVals( i ) = tEnrichmentMap( tOldEnrichLevel );
        }

        // check that for proper new enrichment levels
        MORIS_ERROR( tSortedEnrichmentVals.max() < MORIS_SINT_MAX,
                "Enrichment::sort_enrichment_levels - error in sorted enrichment levels." );

        // overwrite old enrichment levels
        aSubPhaseBinEnrichmentVals = tSortedEnrichmentVals;
    }

    //-------------------------------------------------------------------------------------

    void
    Enrichment::unzip_subphase_bin_enrichment_into_element_enrichment(
            moris_index const &        aEnrichmentDataIndex,
            moris_index const &        aBasisIndex,
            Matrix< IndexMat > const & aParentElementsInSupport,
            Matrix< IndexMat > const & aSubphasesInSupport,
            IndexMap&                  aSubPhaseIndexToSupportIndex,
            Matrix< IndexMat > const & aPrunedSubPhaseToSubphase,
            Matrix< IndexMat >&        aSubPhaseBinEnrichmentVals )
    {
        // resize member data
        moris::size_t tNumAllElementsInSupport                                         = this->count_elements_in_support( aParentElementsInSupport );
        mEnrichmentData( aEnrichmentDataIndex ).mElementIndsInBasis( aBasisIndex )     = Matrix< IndexMat >( 1, tNumAllElementsInSupport );
        mEnrichmentData( aEnrichmentDataIndex ).mElementEnrichmentLevel( aBasisIndex ) = Matrix< IndexMat >( 1, tNumAllElementsInSupport );

        uint tCount = 0;

        // go through all subphases in support and save the current basis function's index to them,
        // and which enrichment level of this basis function is active on them
        for ( moris::size_t iSP = 0; iSP < aSubphasesInSupport.numel(); iSP++ )
        {
            moris_index tSubphaseIndex = aSubphasesInSupport( iSP );

            // iterate through cells in the subphase
            std::shared_ptr< IG_Cell_Group > tSubphaseIgCells = mCutIgMesh->get_subphase_ig_cells( tSubphaseIndex );

            // iterate through cells in subphase
            for ( uint iIgCell = 0; iIgCell < tSubphaseIgCells->mIgCellGroup.size(); iIgCell++ )
            {
                mEnrichmentData( aEnrichmentDataIndex ).mElementIndsInBasis( aBasisIndex )( tCount )     = tSubphaseIgCells->mIgCellGroup( iIgCell )->get_index();
                mEnrichmentData( aEnrichmentDataIndex ).mElementEnrichmentLevel( aBasisIndex )( tCount ) = aSubPhaseBinEnrichmentVals( iSP );
                tCount++;
            }

            // add information to interp cells about which basis/enrichment level interpolates in it
            mEnrichmentData( aEnrichmentDataIndex ).mSubphaseBGBasisIndices( tSubphaseIndex ).push_back( aBasisIndex );
            mEnrichmentData( aEnrichmentDataIndex ).mSubphaseBGBasisEnrLev( tSubphaseIndex ).push_back( aSubPhaseBinEnrichmentVals( iSP ) );
        }
    }

    //-------------------------------------------------------------------------------------

    void
    Enrichment::unzip_subphase_group_bin_enrichment_into_element_enrichment(
            moris_index const &        aEnrichmentDataIndex,
            moris_index const &        aBasisIndex,
            Matrix< IndexMat > const & aParentElementsInSupport,
            Matrix< IndexMat > const & aSpgsInSupport,
            Matrix< IndexMat >&        aSpgBinEnrichmentVals )
    {
        // resize member data
        moris::size_t tNumAllElementsInSupport                                         = this->count_elements_in_support( aParentElementsInSupport );
        mEnrichmentData( aEnrichmentDataIndex ).mElementIndsInBasis( aBasisIndex )     = Matrix< IndexMat >( 1, tNumAllElementsInSupport );
        mEnrichmentData( aEnrichmentDataIndex ).mElementEnrichmentLevel( aBasisIndex ) = Matrix< IndexMat >( 1, tNumAllElementsInSupport );

        // initialize counter tracking number of IG cells in support of basis function
        uint tIgCellCount = 0;

        // go through all subphases in support and save the current basis function's index to them,
        // and which enrichment level of this basis function is active on them
        for ( moris::size_t iSPG = 0; iSPG < aSpgsInSupport.numel(); iSPG++ )
        {
            // get the index of the SPG currently treated
            moris_index tSpgIndex = aSpgsInSupport( iSPG );

            // iterate through cells in the subphase
            Cell< moris_index > tIgCellsInSPG = mCutIgMesh->get_ig_cells_in_SPG( this->get_list_index_for_mesh_index( aEnrichmentDataIndex ), tSpgIndex );

            // iterate through cells in subphase
            for ( uint iIgCell = 0; iIgCell < tIgCellsInSPG.size(); iIgCell++ )
            {
                mEnrichmentData( aEnrichmentDataIndex ).mElementIndsInBasis( aBasisIndex )( tIgCellCount )     = tIgCellsInSPG( iIgCell );
                mEnrichmentData( aEnrichmentDataIndex ).mElementEnrichmentLevel( aBasisIndex )( tIgCellCount ) = aSpgBinEnrichmentVals( iSPG );
                tIgCellCount++;
            }

            // add information to interp cells about which basis/enrichment level interpolates in it
            mEnrichmentData( aEnrichmentDataIndex ).mSubphaseGroupBGBasisIndices( tSpgIndex ).push_back( aBasisIndex );
            mEnrichmentData( aEnrichmentDataIndex ).mSubphaseGroupBGBasisEnrLev( tSpgIndex ).push_back( aSpgBinEnrichmentVals( iSPG ) );
        }
    }

    //-------------------------------------------------------------------------------------

    void
    Enrichment::construct_enriched_basis_to_subphase_connectivity(
            moris_index const &                aEnrichmentDataIndex,
            Cell< Matrix< IndexMat > > const & aSubPhaseBinEnrichment,
            Cell< Matrix< IndexMat > > const & aSubphaseClusterIndicesInSupport,
            Cell< moris_index > const &        aMaxEnrichmentLevel )
    {
        // initialize counter for total number of enriched BFs
        moris_index tNumEnrichmentBasis = 0;

        // count up the number of enriched basis functions
        for ( uint i = 0; i < aMaxEnrichmentLevel.size(); i++ )
        {
            tNumEnrichmentBasis += aMaxEnrichmentLevel( i ) + 1;
        }

        // size data
        mEnrichmentData( aEnrichmentDataIndex ).mSubphaseIndsInEnrichedBasis.resize( tNumEnrichmentBasis );

        moris_index tBaseIndex = 0;

        for ( uint iEnrLvl = 0; iEnrLvl < aSubPhaseBinEnrichment.size(); iEnrLvl++ )
        {
            // get the maximum enrichment level in this basis support
            moris_index tMaxEnrLev = aMaxEnrichmentLevel( iEnrLvl );

            // counter
            Cell< moris_index > tCounter( tMaxEnrLev + 1, 0 );

            // allocate member data for these basis functions
            for ( moris_index iEnr = tBaseIndex; iEnr < tBaseIndex + tMaxEnrLev + 1; iEnr++ )
            {
                mEnrichmentData( aEnrichmentDataIndex ).mSubphaseIndsInEnrichedBasis( iEnr ).resize( 1, aSubphaseClusterIndicesInSupport( iEnrLvl ).numel() );
            }

            // iterate through subphases in support and add them to appropriate location in mSubphaseIndsInEnrichedBasis
            for ( uint iSp = 0; iSp < aSubphaseClusterIndicesInSupport( iEnrLvl ).numel(); iSp++ )
            {
                // get cluster enrichment level
                moris_index tClusterEnrLev = aSubPhaseBinEnrichment( iEnrLvl )( iSp );

                // add to the member data
                mEnrichmentData( aEnrichmentDataIndex ).mSubphaseIndsInEnrichedBasis( tBaseIndex + tClusterEnrLev )( tCounter( tClusterEnrLev ) ) =
                        aSubphaseClusterIndicesInSupport( iEnrLvl )( iSp );

                // increment count
                tCounter( tClusterEnrLev )++;
            }

            // size out unused space
            for ( moris_index iEnr = 0; iEnr < tMaxEnrLev + 1; iEnr++ )
            {
                moris_index tIndex = tBaseIndex + iEnr;
                mEnrichmentData( aEnrichmentDataIndex ).mSubphaseIndsInEnrichedBasis( tIndex ).resize( 1, tCounter( iEnr ) );

                // sort in ascending order (easier to find in MPI)
                // if this sort is removed the function  subphase_is_in_support needs to be updated
                moris::sort( mEnrichmentData( aEnrichmentDataIndex ).mSubphaseIndsInEnrichedBasis( tIndex ),
                        mEnrichmentData( aEnrichmentDataIndex ).mSubphaseIndsInEnrichedBasis( tIndex ),
                        "ascend",
                        1 );
            }

            // update starting index
            tBaseIndex = tBaseIndex + tMaxEnrLev + 1;
        }
    }

    //-------------------------------------------------------------------------------------

    void
    Enrichment::construct_enriched_basis_to_subphase_group_connectivity(
            moris_index const &                aEnrichmentDataIndex,
            Cell< Matrix< IndexMat > > const & aSpgBinEnrichment,
            Cell< Matrix< IndexMat > > const & aSpgIndicesInSupport,
            Cell< moris_index > const &        aMaxEnrichmentLevel )
    {
        // get easy access to the enrichment data
        xtk::Enrichment_Data* tEnrichmentData = &mEnrichmentData( aEnrichmentDataIndex );

        // initialize counter for total number of enriched BFs
        moris_index tNumEnrichmentBasis = 0;

        // count up the number of enriched basis functions
        for ( uint i = 0; i < aMaxEnrichmentLevel.size(); i++ )
        {
            tNumEnrichmentBasis += aMaxEnrichmentLevel( i ) + 1;
        }

        // get the position of the DMI in the list of B-spline meshes
        moris_index tMeshListIndex = this->get_list_index_for_mesh_index( aEnrichmentDataIndex );

        // get access to the Bspline mesh info for the current mesh
        Bspline_Mesh_Info* tBsplineMeshInfo = mBsplineMeshInfos( tMeshListIndex );

        // size data
        tEnrichmentData->mSubphaseGroupIndsInEnrichedBasis.resize( tNumEnrichmentBasis );
        tEnrichmentData->mSubphaseIndsInEnrichedBasis.resize( tNumEnrichmentBasis );
        tEnrichmentData->mBulkPhaseInEnrichedBasis.set_size( 1, tNumEnrichmentBasis, gNoID );

        moris_index tBaseIndex = 0;

        for ( uint iBaseBF = 0; iBaseBF < aSpgBinEnrichment.size(); iBaseBF++ )
        {
            // get the maximum enrichment level in this basis support
            moris_index tMaxEnrLev = aMaxEnrichmentLevel( iBaseBF );

            // initialize counters counting SPGs and SPs in the basis support
            Cell< moris_index > tSpgCounter( tMaxEnrLev + 1, 0 );
            Cell< moris_index > tSpCounter( tMaxEnrLev + 1, 0 );

            // get the number of SPGs in the support of the enr. BF
            uint tNumSpgIndicesInSupport = aSpgIndicesInSupport( iBaseBF ).numel();

            // get the number of SPs in the support of the enr. BF
            uint tNumSpsInSupport = 0;
            for ( uint iSPG = 0; iSPG < tNumSpgIndicesInSupport; iSPG++ )
            {
                // get the current SPG index
                moris_index tSpgIndex = aSpgIndicesInSupport( iBaseBF )( iSPG );

                // count up number of SPs in support
                tNumSpsInSupport += tBsplineMeshInfo->mSubphaseGroups( tSpgIndex )->get_num_SPs_in_group();
            }

            // allocate member data for these basis functions
            for ( moris_index iEnr = tBaseIndex; iEnr < tBaseIndex + tMaxEnrLev + 1; iEnr++ )
            {
                tEnrichmentData->mSubphaseGroupIndsInEnrichedBasis( iEnr ).resize( 1, tNumSpgIndicesInSupport );
                tEnrichmentData->mSubphaseIndsInEnrichedBasis( iEnr ).resize( 1, tNumSpsInSupport );
            }

            // iterate through SPGs in support and add them to appropriate location in mSubphaseGroupIndsInEnrichedBasis
            for ( uint iSPG = 0; iSPG < tNumSpgIndicesInSupport; iSPG++ )
            {
                // get cluster enrichment level
                moris_index tClusterEnrLev = aSpgBinEnrichment( iBaseBF )( iSPG );

                // get the current SPGs index
                moris_index tSpgIndex = aSpgIndicesInSupport( iBaseBF )( iSPG );

                // count up index of the enriched Basis function
                moris_index tEnrBfIndex = tBaseIndex + tClusterEnrLev;

                // add to the member data
                tEnrichmentData->mSubphaseGroupIndsInEnrichedBasis( tEnrBfIndex )( tSpgCounter( tClusterEnrLev ) ) = tSpgIndex;

                // increment count of SPGs in the enriched basis' support
                tSpgCounter( tClusterEnrLev )++;

                // get SPs on current SPG
                Cell< moris_index > const & tSpIndicesInSpg = tBsplineMeshInfo->mSubphaseGroups( tSpgIndex )->get_SP_indices_in_group();
                uint                        tNumSpsOnSpg    = tSpIndicesInSpg.size();

                // fill list of SPs in basis support
                for ( uint iSP = 0; iSP < tNumSpsOnSpg; iSP++ )
                {
                    // add SP to list of SPs in basis support
                    tEnrichmentData->mSubphaseIndsInEnrichedBasis( tEnrBfIndex )( tSpCounter( tClusterEnrLev ) ) = tSpIndicesInSpg( iSP );

                    // count up number of SPs in the enriched basis' support
                    tSpCounter( tClusterEnrLev )++;
                }

                // find the bulk-phase index corresponding to the SPG and put it in map
                moris_index tPrevSetBpIndex = tEnrichmentData->mBulkPhaseInEnrichedBasis( tEnrBfIndex );
                moris_index tBulkPhaseIndex = mBsplineMeshInfos( tMeshListIndex )->get_bulk_phase_for_subphase_group( tSpgIndex );

                // check if the bulk phase has not already been set for this enr. BF
                if ( tPrevSetBpIndex == -1 )
                {
                    tEnrichmentData->mBulkPhaseInEnrichedBasis( tEnrBfIndex ) = tBulkPhaseIndex;
                }
#ifdef MORIS_HAVE_DEBUG
                else    // if the bulk phase has already been set for this enr. BF  ...
                {
                    // ... check if the newly found one is still the same as the one previously found
                    MORIS_ASSERT(
                            tPrevSetBpIndex == tBulkPhaseIndex,
                            "Enrichment::construct_enriched_basis_to_subphase_group_connectivity() - "
                            "SPGs of different Bulk-phases associated with the same enriched BF. Something is wrong" );
                }
#endif
            }    // end: loop over all SPGs in support of current Base BF

            // size out unused space
            for ( moris_index iEnrLvl = 0; iEnrLvl < tMaxEnrLev + 1; iEnrLvl++ )
            {
                moris_index tIndex = tBaseIndex + iEnrLvl;
                tEnrichmentData->mSubphaseGroupIndsInEnrichedBasis( tIndex ).resize( 1, tSpgCounter( iEnrLvl ) );

                // sort in ascending order (easier to find in MPI)
                // if this sort is removed the function subphase_is_in_support needs to be updated
                moris::sort(
                        tEnrichmentData->mSubphaseGroupIndsInEnrichedBasis( tIndex ),
                        tEnrichmentData->mSubphaseGroupIndsInEnrichedBasis( tIndex ),
                        "ascend",
                        1 );
            }

            // update starting index
            tBaseIndex = tBaseIndex + tMaxEnrLev + 1;

        }    // end for: each non-enriched Basis function

    }        // end function:

    //-------------------------------------------------------------------------------------

    void
    Enrichment::assign_enriched_coefficients_identifiers(
            moris_index const &         aEnrichmentDataIndex,
            Cell< moris_index > const & aMaxEnrichmentLevel )
    {
        // get number of non-enriched BFs
        uint tNumNonEnrichedBFs = mEnrichmentData( aEnrichmentDataIndex ).mElementIndsInBasis.size();

        // initialize array holding the enriched basis function function indices living on any given non-enriched BF
        mEnrichmentData( aEnrichmentDataIndex ).mBasisEnrichmentIndices.resize( tNumNonEnrichedBFs );

        // initialize counter for number of enriched BFs
        mEnrichmentData( aEnrichmentDataIndex ).mNumEnrichedBasisFunctions = 0;

        // continue to initialize array holding the enriched basis function indices living on any given non-enriched BF
        uint tNumEnrichedBasisFunctionCount = 0;
        for ( uint iBF = 0; iBF < tNumNonEnrichedBFs; iBF++ )
        {
            moris_index tMaxEnrLev = aMaxEnrichmentLevel( iBF ) + 1;
            tNumEnrichedBasisFunctionCount += (uint)tMaxEnrLev;
            mEnrichmentData( aEnrichmentDataIndex ).mBasisEnrichmentIndices( iBF ) = Matrix< IndexMat >( 1, tMaxEnrLev );
        }

        // allocate enriched basis index to id data
        mEnrichmentData( aEnrichmentDataIndex ).mNumEnrichedBasisFunctions = tNumEnrichedBasisFunctionCount;
        mEnrichmentData( aEnrichmentDataIndex ).mEnrichedBasisIndexToId.resize( 1, tNumEnrichedBasisFunctionCount );
        mEnrichmentData( aEnrichmentDataIndex ).mEnrichedBasisIndexToId.fill( MORIS_INDEX_MAX );

        // get current processor's ID
        moris_index tParRank = par_rank();

        // get the XTK comm table
        Matrix< IndexMat > tCommTable = mXTKModelPtr->get_communication_table();

        // Initialize map relating the global MPI proc rank to its position in the XTK Comm Table
        Cell< moris_index > tProcRanks( tCommTable.numel() );

        std::unordered_map< moris_id, moris_id > tProcRankToIndexInData;

        // relate the global MPI proc rank to its position in the XTK Comm Table
        for ( uint iProc = 0; iProc < tCommTable.numel(); iProc++ )
        {
            tProcRankToIndexInData[ tCommTable( iProc ) ] = iProc;

            tProcRanks( iProc ) = ( tCommTable( iProc ) );
        }

        // get the first first free global ID (not first gets background basis information)
        moris_index tIndOffset     = 0;
        moris_id    tBasisIdOffset = this->allocate_basis_ids( aEnrichmentDataIndex, mEnrichmentData( aEnrichmentDataIndex ).mNumEnrichedBasisFunctions );

        Cell< Cell< moris_index > > tBasisIdToBasisOwner( tCommTable.numel() );
        Cell< Cell< moris_index > > tSubphaseIdInSupport( tCommTable.numel() );
        Cell< Cell< moris_index > > tBasisIndexToBasisOwner( tCommTable.numel() );

        // for each non-enriched BF ...
        for ( uint iNonEnrichedBF = 0; iNonEnrichedBF < tNumNonEnrichedBFs; iNonEnrichedBF++ )
        {
            // ... get their basis owner and ...
            moris_index tOwner = mBackgroundMeshPtr->get_entity_owner(
                    iNonEnrichedBF,
                    mBasisRank,
                    aEnrichmentDataIndex );

            // .. get their basis ID
            moris_id tBackBasisId = mBackgroundMeshPtr->get_glb_entity_id_from_entity_loc_index(
                    iNonEnrichedBF,
                    mBasisRank,
                    aEnrichmentDataIndex );

            // get the owning processor's position in the communication arrays
            moris_index tProcDataIndex = tProcRankToIndexInData[ tOwner ];

            // get access to the list of enriched BF indices living on the current non-enriched BF
            Matrix< IndexMat >& tBasisEnrichmentInds =    //
                    mEnrichmentData( aEnrichmentDataIndex ).mBasisEnrichmentIndices( iNonEnrichedBF );

            uint tNumEnrichedBFsOnBasis = tBasisEnrichmentInds.numel();

            // assign indices to each enriched Basis function by counting up the enriched bases living on each non-enriched Basis
            for ( uint jEnrBF = 0; jEnrBF < tNumEnrichedBFsOnBasis; jEnrBF++ )
            {
                tBasisEnrichmentInds( jEnrBF ) = tIndOffset;
                tIndOffset++;
            }

            // only set id if we own it and package data for communication if shared
            if ( tOwner == tParRank )
            {
                // get access to the non-enriched basis' ID
                mEnrichmentData( aEnrichmentDataIndex ).mEnrichedBasisIndexToId( tBasisEnrichmentInds( 0 ) ) = tBackBasisId;

                // give all enr. BFs their IDs
                for ( uint jEnrBF = 1; jEnrBF < tNumEnrichedBFsOnBasis; jEnrBF++ )
                {
                    // get the current enriched basis' index
                    moris_index tEnrBfInd = tBasisEnrichmentInds( jEnrBF );

                    // check that the basis doesn't already have an ID attached to it
                    MORIS_ASSERT( mEnrichmentData( aEnrichmentDataIndex ).mEnrichedBasisIndexToId( tEnrBfInd ) == MORIS_INDEX_MAX,
                            "Enrichment::assign_enriched_coefficients_identifiers() - Already set enriched basis id" );

                    // assign ID to enriched BF
                    mEnrichmentData( aEnrichmentDataIndex ).mEnrichedBasisIndexToId( tEnrBfInd ) = tBasisIdOffset;

                    // increment ID for next enr. BF
                    tBasisIdOffset++;
                }
            }

            // if we don't own the basis setup the communication to get the basis
            else
            {
                // prepare ID requests for all enr. BFs
                for ( uint jEnrBF = 0; jEnrBF < tNumEnrichedBFsOnBasis; jEnrBF++ )
                {
                    // get the current enriched basis' index
                    moris_index tEnrichedBasisIndex = tBasisEnrichmentInds( jEnrBF );

                    // TODO: once the SPGs have IDs and parallel consistent, this here needs to be changed to use SPGs
                    moris_index tFirstSubphaseInSupportIndex = mEnrichmentData( aEnrichmentDataIndex ).mSubphaseIndsInEnrichedBasis( tEnrichedBasisIndex )( 0 );
                    // moris_index tFirstSpgInSupportIndex = mEnrichmentData( aEnrichmentDataIndex ).mSubphaseGroupIndsInEnrichedBasis( tEnrichedBasisIndex )( 0 );

                    moris_index tFirstSubphaseInSupportId = mXTKModelPtr->get_subphase_id( tFirstSubphaseInSupportIndex );
                    // moris_index tFirstSubphaseGroupInSupportId = mXTKModelPtr->get_subphase_group_id( tFirstSpgInSupportIndex );

                    tBasisIdToBasisOwner( tProcDataIndex ).push_back( tBackBasisId );
                    tSubphaseIdInSupport( tProcDataIndex ).push_back( tFirstSubphaseInSupportId );
                    // tSubphaseGroupIdInSupport( tProcDataIndex ).push_back( tFirstSubphaseGroupInSupportId );
                    tBasisIndexToBasisOwner( tProcDataIndex ).push_back( tEnrichedBasisIndex );
                }
            }
        }

        // send information about not owned enriched basis to owner processor
        Cell< Matrix< IndexMat > > tNotOwnedEnrichedBasisId;

        this->communicate_basis_information_with_owner(
                aEnrichmentDataIndex,
                tBasisIdToBasisOwner,
                tSubphaseIdInSupport,    // tSubphaseGroupIdInSupport,
                tProcRanks,
                tProcRankToIndexInData,
                tNotOwnedEnrichedBasisId );

        // set the received information in my data
        this->set_received_enriched_basis_ids(
                aEnrichmentDataIndex,
                tNotOwnedEnrichedBasisId,
                tBasisIndexToBasisOwner );
    }

    //-------------------------------------------------------------------------------------

    void
    Enrichment::assign_enriched_coefficients_identifiers_new(
            moris_index const &         aEnrichmentDataIndex,
            Cell< moris_index > const & aMaxEnrichmentLevel )
    {
        // log this function when verbose output is requested
        Tracer tTracer( "XTK", "Enrichment", "assign enriched basis function indices and IDs", mXTKModelPtr->mVerboseLevel, 1 );

        // access current enrichment data
        xtk::Enrichment_Data* tEnrichmentData = &mEnrichmentData( aEnrichmentDataIndex );

        /* ---------------------------------------------------------------------------------------- */
        /* Step 0: assign indices and store the owned and non-owned entities */

        this->sort_enriched_coefficients_into_owned_and_not_owned( aEnrichmentDataIndex, aMaxEnrichmentLevel );

        /* ---------------------------------------------------------------------------------------- */
        /* Step 0.5: access the comm table */

        // get the communication table
        Matrix< IdMat > tCommTable = mCutIgMesh->get_communication_table();

        /* ---------------------------------------------------------------------------------------- */
        /* Step 1: Let each proc decide how many entity IDs it needs & communicate ID ranges */

        // get the maximum BF ID already assigned on this proc
        moris_index tLocMaxId = this->get_max_basis_id( aEnrichmentDataIndex );

        // get the first free BF ID globally
        moris_index tFirstFreeEntityId = moris::max_all( tLocMaxId ) + 1;

        // Get the number of subphases (on the current proc)
        moris_id tNumEntitiesToAllocate = (moris_id)tEnrichmentData->mNumEnrichedBasisFunctions;

        // reserve IDs for his proc
        moris_id tMyFirstId = get_processor_offset( tNumEntitiesToAllocate ) + tFirstFreeEntityId;

        /* ---------------------------------------------------------------------------------------- */
        /* Step 2: Assign IDs to owned entities */

        this->assign_IDs_to_owned_enriched_coefficients( aEnrichmentDataIndex, tMyFirstId );

        /* ---------------------------------------------------------------------------------------- */
        /* The following steps are only necessary if code runs in parallel */

        if ( par_size() == 1 )    // serial
        {
            // check that all entities are owned in serial
            MORIS_ASSERT( tEnrichmentData->mNotOwnedEnrBasisIndices.size() == 0,
                    "Enrichment::assign_enriched_coefficients_identifiers_new() - "
                    "Code running in serial, but not all entities are owned by proc 0." );
        }
        else    // parallel
        {
            /* ---------------------------------------------------------------------------------------- */
            /* Step 3: Prepare requests for non-owned entities */

            // initialize lists of information that identifies the entities (on other procs)
            Cell< Cell< moris_index > > tNotOwnedEnrBfsToProcs;       // enr. BF indices communicated with other procs
            Cell< Matrix< IdMat > >     tNonEnrBasisIDs;              // non-enriched BF IDs which the enr. BF is constructed from
            Cell< Matrix< IdMat > >     tSubphaseGroupIdInSupport;    // Sub-phase group ID the enr. basis function interpolates into

            // fill information
            this->prepare_requests_for_not_owned_enriched_coefficient_IDs(
                    aEnrichmentDataIndex,
                    tNotOwnedEnrBfsToProcs,
                    tNonEnrBasisIDs,
                    tSubphaseGroupIdInSupport );

            /* ---------------------------------------------------------------------------------------- */
            /* Step 4: Send and Receive requests about non-owned entities to and from other procs */

            // initialize arrays for receiving
            Cell< Matrix< IdMat > > tReceivedNonEnrBasisIDs;
            Cell< Matrix< IdMat > > tReceivedSubphaseGroupIdInSupport;

            // communicate information
            moris::communicate_mats( tCommTable, tNonEnrBasisIDs, tReceivedNonEnrBasisIDs );
            moris::communicate_mats( tCommTable, tSubphaseGroupIdInSupport, tReceivedSubphaseGroupIdInSupport );

            // clear memory not needed anymore
            tNonEnrBasisIDs.clear();
            tSubphaseGroupIdInSupport.clear();

            /* ---------------------------------------------------------------------------------------- */
            /* Step 5: Find answers to the requests */

            // initialize lists of ID answers to other procs
            Cell< Matrix< IdMat > > tEnrBfIds;

            // find the IDs for the enr. BFs requested by the other processors
            this->prepare_answers_for_owned_enriched_coefficient_IDs(
                    aEnrichmentDataIndex,
                    tReceivedNonEnrBasisIDs,
                    tReceivedSubphaseGroupIdInSupport,
                    tEnrBfIds );

            // clear memory from requests (the answers to which have been found)
            tReceivedNonEnrBasisIDs.clear();
            tReceivedSubphaseGroupIdInSupport.clear();

            /* ---------------------------------------------------------------------------------------- */
            /* Step 6: Send and receive answers to and from other procs */

            // initialize arrays for receiving
            Cell< Matrix< IdMat > > tReceivedEnrBfIds;

            // communicate answers
            moris::communicate_mats( tCommTable, tEnrBfIds, tReceivedEnrBfIds );

            // clear unused memory
            tEnrBfIds.clear();

            /* ---------------------------------------------------------------------------------------- */
            /* Step 7: Use answers to assign IDs to non-owned entities */

            this->handle_requested_unzipped_enriched_coefficient_answers(
                    aEnrichmentDataIndex,
                    tNotOwnedEnrBfsToProcs,
                    tReceivedEnrBfIds );

        }    // end if: parallel

    }        // end function: assign_enriched_coefficients_identifiers_new()

    //-------------------------------------------------------------------------------------

    void
    Enrichment::sort_enriched_coefficients_into_owned_and_not_owned(
            moris_index const &         aEnrichmentDataIndex,
            Cell< moris_index > const & aMaxEnrichmentLevel )
    {
        // access current enrichment data
        xtk::Enrichment_Data* tEnrichmentData = &mEnrichmentData( aEnrichmentDataIndex );

        // get number of non-enriched BFs
        uint tNumNonEnrichedBFs = tEnrichmentData->mElementIndsInBasis.size();

        // counters for owned and not-owned enr. basis functions
        uint tNumOwnedEnrBFs    = 0;
        uint tNumNotOwnedEnrBFs = 0;

        // counter for total number of enr. BFs
        uint tNumEnrBFs = 0;

        // count number of owned and not owned enr. basis functions
        for ( uint iNonEnrBF = 0; iNonEnrBF < tNumNonEnrichedBFs; iNonEnrBF++ )
        {
            //  get the non-enriched BF owner
            moris_index tOwner = mBackgroundMeshPtr->get_entity_owner(
                    iNonEnrBF,
                    mBasisRank,
                    aEnrichmentDataIndex );

            // get the number of enrichment levels on the current BF obtained from the flood fill
            uint tNumEnrLvlsOnBf = (uint)aMaxEnrichmentLevel( iNonEnrBF ) + 1;

            // count up total number of enr. BFs
            tNumEnrBFs += tNumEnrLvlsOnBf;

            // count towards owned or not-owned basis functions
            if ( tOwner == par_rank() )    // owned
            {
                tNumOwnedEnrBFs += tNumEnrLvlsOnBf;
            }
            else    // not owned
            {
                tNumNotOwnedEnrBFs += tNumEnrLvlsOnBf;
            }
        }

        // store total number of enr. BFs
        tEnrichmentData->mNumEnrichedBasisFunctions = tNumEnrBFs;

        // check that the number of owned and not owned basis functions adds up
        MORIS_ASSERT( tNumOwnedEnrBFs + tNumNotOwnedEnrBFs == tNumEnrBFs,
                "Enrichment::sort_enriched_coefficients_into_owned_and_not_owned() - "
                "Numbers of owned and not-owned enriched basis functions do not add up "
                "to the total number of enr. basis functions: %i + %i =/= %i",
                tNumOwnedEnrBFs,
                tNumNotOwnedEnrBFs,
                tNumEnrBFs );

        // initialize lists of owned and non-owned enr. basis function indices
        tEnrichmentData->mOwnedEnrBasisIndices.resize( tNumOwnedEnrBFs );
        tEnrichmentData->mNotOwnedEnrBasisIndices.resize( tNumNotOwnedEnrBFs );

        // initialize enriched basis index to ID data
        tEnrichmentData->mEnrichedBasisIndexToId.resize( 1, tNumEnrBFs );
        tEnrichmentData->mEnrichedBasisIndexToId.fill( MORIS_ID_MAX );

        // initialize map relating the enr. BF indices back to their non-enriched counterpart
        tEnrichmentData->mNonEnrBfIndForEnrBfInd.resize( tNumEnrBFs, MORIS_INDEX_MAX );
        tEnrichmentData->mEnrLvlOfEnrBf.resize( tNumEnrBFs, MORIS_INDEX_MAX );

        // initialize array holding the enriched basis function indices living on any given non-enriched BF
        tEnrichmentData->mBasisEnrichmentIndices.resize( tNumNonEnrichedBFs );

        // initialize counters for basis functions indices
        moris_index tOwnedEnrBfIndexCounter    = 0;
        moris_index tNotOwnedEnrBfIndexCounter = 0;
        moris_index tEnrBfIndexCounter         = 0;

        // continue to initialize array holding the enriched basis function indices living on any given non-enriched BF
        for ( uint iNonEnrBF = 0; iNonEnrBF < tNumNonEnrichedBFs; iNonEnrBF++ )
        {
            //  get the non-enriched BF owner
            moris_index tOwner = mBackgroundMeshPtr->get_entity_owner(
                    iNonEnrBF,
                    mBasisRank,
                    aEnrichmentDataIndex );

            // get the number of enrichment levels on the current BF obtained from the flood fill
            uint tNumEnrLvlsOnBf = (uint)aMaxEnrichmentLevel( iNonEnrBF ) + 1;

            // get access to the list of enriched BF indices living on the current non-enriched BF
            Matrix< IndexMat >& tBasisEnrichmentInds = tEnrichmentData->mBasisEnrichmentIndices( iNonEnrBF );

            // count up the total number of enriched BFs
            tBasisEnrichmentInds = Matrix< IndexMat >( tNumEnrLvlsOnBf, 1 );

            // assign indices to each enriched basis function by counting up the enriched bases living on each non-enriched Basis,
            // also
            for ( uint iEnrBf = 0; iEnrBf < tNumEnrLvlsOnBf; iEnrBf++ )
            {
                // store as either owned or not owned basis function
                if ( tOwner == par_rank() )    // owned BF
                {
                    tEnrichmentData->mOwnedEnrBasisIndices( tOwnedEnrBfIndexCounter ) = tEnrBfIndexCounter;
                    tOwnedEnrBfIndexCounter++;
                }
                else    // not owned BF
                {
                    tEnrichmentData->mNotOwnedEnrBasisIndices( tNotOwnedEnrBfIndexCounter ) = tEnrBfIndexCounter;
                    tNotOwnedEnrBfIndexCounter++;
                }

                // relate enr. BF back to its non-enriched basis
                tEnrichmentData->mNonEnrBfIndForEnrBfInd( tEnrBfIndexCounter ) = iNonEnrBF;
                tEnrichmentData->mEnrLvlOfEnrBf( tEnrBfIndexCounter )          = iEnrBf;

                // assign enr. basis function indices
                tBasisEnrichmentInds( iEnrBf ) = tEnrBfIndexCounter;
                tEnrBfIndexCounter++;
            }

        }    // end for: each background (non-enriched) basis function

    }        // end function: Enrichment::sort_enriched_coefficients_into_owned_and_not_owned()

    //-------------------------------------------------------------------------------------

    void
    Enrichment::assign_IDs_to_owned_enriched_coefficients(
            moris_index const & aEnrichmentDataIndex,
            moris_id&           aFirstEnrBasisId )
    {
        // access current enrichment data
        xtk::Enrichment_Data* tEnrichmentData = &mEnrichmentData( aEnrichmentDataIndex );

        // assign IDs to owned entities and store non-owned entities
        for ( uint iOwnedEnrBF = 0; iOwnedEnrBF < tEnrichmentData->mOwnedEnrBasisIndices.size(); iOwnedEnrBF++ )
        {
            // get the index of the enriched BF
            moris_index tEnrBfIndex = tEnrichmentData->mOwnedEnrBasisIndices( iOwnedEnrBF );

            // get the enrichment level of this basis wrt the background basis
            moris_index tBfEnrLvl = tEnrichmentData->mEnrLvlOfEnrBf( tEnrBfIndex );

            // if this is the first enr. BF on this background basis, give it the same ID as the background basis
            // this convention is used in other parts of the code
            if ( tBfEnrLvl == 0 )
            {
                // get the background basis index and ID
                moris_index tBaseBfIndex = tEnrichmentData->mNonEnrBfIndForEnrBfInd( tEnrBfIndex );
                moris_id    tBaseBfId    = mBackgroundMeshPtr->get_glb_entity_id_from_entity_loc_index( tBaseBfIndex, mBasisRank, aEnrichmentDataIndex );

                // assign this ID
                mEnrichmentData( aEnrichmentDataIndex ).mEnrichedBasisIndexToId( tEnrBfIndex ) = tBaseBfId;
            }
            // otherwise, count up the enr. BF IDs
            else
            {
                // assign its ID
                mEnrichmentData( aEnrichmentDataIndex ).mEnrichedBasisIndexToId( tEnrBfIndex ) = aFirstEnrBasisId;

                // count number of assigned Basis IDs
                aFirstEnrBasisId++;
            }

        }    // end for: each non-enriched basis function

    }        // end function: Enrichment::assign_IDs_to_owned_enriched_coefficients()

    //-------------------------------------------------------------------------------------

    void
    Enrichment::prepare_requests_for_not_owned_enriched_coefficient_IDs(
            moris_index const &          aEnrichmentDataIndex,
            Cell< Cell< moris_index > >& aNotOwnedEnrBfsToProcs,
            Cell< Matrix< IdMat > >&     aNonEnrBasisIDs,
            Cell< Matrix< IdMat > >&     aSubphaseGroupIdInSupport )
    {
        // access current enrichment data
        xtk::Enrichment_Data* tEnrichmentData = &mEnrichmentData( aEnrichmentDataIndex );

        // get the communication table and map
        Matrix< IdMat >                   tCommTable              = mCutIgMesh->get_communication_table();
        uint                              tCommTableSize          = tCommTable.numel();
        std::map< moris_id, moris_index > tProcIdToCommTableIndex = mCutIgMesh->get_communication_map();

        // initialize lists of identifying information
        aNotOwnedEnrBfsToProcs.resize( tCommTableSize );
        aNonEnrBasisIDs.resize( tCommTableSize );
        aSubphaseGroupIdInSupport.resize( tCommTableSize );

        // initialize map relating enriched BF position in the request arrays with the position in the not-owned arrays
        Cell< Cell< moris_index > > tEnrBfPositionsInNotOwnedList( tCommTableSize );

        // get the non-owned entities on the executing processor
        Cell< moris_index > const & tNotOwnedEnrBFs = tEnrichmentData->mNotOwnedEnrBasisIndices;

        // go through entities that the executing proc knows about, but doesn't own
        // and sort them into bins with other entities that are owned by each proc communicated with
        for ( uint iNotOwnedEnrBF = 0; iNotOwnedEnrBF < tNotOwnedEnrBFs.size(); iNotOwnedEnrBF++ )
        {
            // get the index of the enr. BF
            moris_index tEnrBfIndex = tNotOwnedEnrBFs( iNotOwnedEnrBF );

            // to determine their owner, first determine their non-enriched basis index
            moris_index tNonEnrBfIndex = tEnrichmentData->mNonEnrBfIndForEnrBfInd( tEnrBfIndex );

            // get the respective owners and position in the comm-table ...
            moris_index tOwnerProc = mBackgroundMeshPtr->get_entity_owner( tNonEnrBfIndex, mBasisRank, aEnrichmentDataIndex );
            auto        tIter      = tProcIdToCommTableIndex.find( tOwnerProc );
            MORIS_ASSERT(
                    tIter != tProcIdToCommTableIndex.end(),
                    "Integration_Mesh_Generator::prepare_requests_for_not_owned_enriched_coefficient_IDs() - "
                    "Basis function owner (Proc #%i) not found in communication table of current proc #%i which is: %s",
                    tOwnerProc,
                    par_rank(),
                    ios::stringify_log( tCommTable ).c_str() );
            moris_index tProcDataIndex = tIter->second;

            // ... and finally add the non-owned SP in the list of SPs to be requested from that owning proc
            aNotOwnedEnrBfsToProcs( tProcDataIndex ).push_back( tEnrBfIndex );

            // store where in the not owned list of enr. BFs this enr. BF sits
            tEnrBfPositionsInNotOwnedList( tProcDataIndex ).push_back( iNotOwnedEnrBF );
        }

        // for each processor communicated with assemble the identifying information of the entities communicated
        for ( uint iProc = 0; iProc < tCommTableSize; iProc++ )
        {
            // get the number of non-owned entities to be sent to each processor
            uint tNumEntitiesToSendToProc = aNotOwnedEnrBfsToProcs( iProc ).size();

            // resize arrays to correct size
            aNonEnrBasisIDs( iProc ).resize( tNumEntitiesToSendToProc, 1 );
            aSubphaseGroupIdInSupport( iProc ).resize( tNumEntitiesToSendToProc, 1 );

            // go through the entities for which IDs will need to be answered by the other processor
            for ( uint iEnrBf = 0; iEnrBf < tNumEntitiesToSendToProc; iEnrBf++ )
            {
                // get the index of the entity on the executing proc
                moris_index tPosInNotOwnedList = tEnrBfPositionsInNotOwnedList( iProc )( iEnrBf );
                moris_index tEnrBfIndex        = tNotOwnedEnrBFs( tPosInNotOwnedList );

                // get the ID of the non-enriched BF the enriched BF is constructed from
                moris_index tNonEnrBfIndex = tEnrichmentData->mNonEnrBfIndForEnrBfInd( tEnrBfIndex );
                moris_id    tNonEnrBfId    = mBackgroundMeshPtr->get_glb_entity_id_from_entity_loc_index( tNonEnrBfIndex, mBasisRank, aEnrichmentDataIndex );

                // check that the ID obtained from HMR makes sense
                MORIS_ASSERT( tNonEnrBfId != MORIS_ID_MAX && tNonEnrBfId > -1,
                        "Integration_Mesh_Generator::prepare_requests_for_not_owned_enriched_coefficient_IDs() - "
                        "HMR mesh returned MORIS_ID_MAX for background basis function ID." );

                // get teh ID of the first SPG in the support of the enr. BF
                moris_index tFirstSpgInSupportIndex = tEnrichmentData->mSubphaseGroupIndsInEnrichedBasis( tEnrBfIndex )( 0 );
                moris_id    tFirstSpgInSupportId    = mXTKModelPtr->get_subphase_group_id( tFirstSpgInSupportIndex, aEnrichmentDataIndex );

                // store identifying information
                aNonEnrBasisIDs( iProc )( iEnrBf )           = tNonEnrBfId;
                aSubphaseGroupIdInSupport( iProc )( iEnrBf ) = tFirstSpgInSupportId;
            }

        }    // end for: each proc communicated with

        // size out unused memory
        aNotOwnedEnrBfsToProcs.shrink_to_fit();
        aNonEnrBasisIDs.shrink_to_fit();
        aSubphaseGroupIdInSupport.shrink_to_fit();

    }    // end function Enrichment::prepare_requests_for_not_owned_enriched_coefficient_IDs()

    //-------------------------------------------------------------------------------------

    void
    Enrichment::prepare_answers_for_owned_enriched_coefficient_IDs(
            moris_index const &             aEnrichmentDataIndex,
            Cell< Matrix< IdMat > > const & aReceivedNonEnrBasisIDs,
            Cell< Matrix< IdMat > > const & aReceivedSubphaseGroupIdInSupport,
            Cell< Matrix< IdMat > >&        aEnrBfIds )
    {
        // access current enrichment data
        xtk::Enrichment_Data* tEnrichmentData = &mEnrichmentData( aEnrichmentDataIndex );

        // get the communication table
        Matrix< IdMat > tCommTable     = mCutIgMesh->get_communication_table();
        uint            tCommTableSize = tCommTable.numel();

        // initialize answer array with correct size
        aEnrBfIds.resize( tCommTableSize );

        // check that the received data is complete
        MORIS_ASSERT(
                aReceivedNonEnrBasisIDs.size() == tCommTableSize && aReceivedSubphaseGroupIdInSupport.size() == tCommTableSize,
                "Enrichment::assign_enriched_coefficients_identifiers_new() - Received information incomplete." );

        // go through the list of processors in the array of ID requests
        for ( uint iProc = 0; iProc < tCommTableSize; iProc++ )
        {
            // get the number of SP IDs requested from the current proc position
            uint tNumReceivedReqs = aReceivedNonEnrBasisIDs( iProc ).numel();
            MORIS_ASSERT(
                    tNumReceivedReqs == aReceivedSubphaseGroupIdInSupport( iProc ).numel(),
                    "Enrichment::assign_enriched_coefficients_identifiers_new() - "
                    "Size of received arrays from proc #%i inconsistent.",
                    tCommTable( iProc ) );

            // size the list of answers / IDs accordingly
            aEnrBfIds( iProc ).resize( tNumReceivedReqs, 1 );

            // iterate through entities for which the IDs are requested
            for ( uint iEnrBF = 0; iEnrBF < tNumReceivedReqs; iEnrBF++ )
            {
                // access identifying information
                moris_id tNonEnrBfId = aReceivedNonEnrBasisIDs( iProc )( iEnrBF );
                moris_id tFirstSpgId = aReceivedSubphaseGroupIdInSupport( iProc )( iEnrBF );

                // convert IDs to local indices
                moris_index tNonEnrBfIndex = mBackgroundMeshPtr->get_loc_entity_ind_from_entity_glb_id( tNonEnrBfId, mBasisRank, aEnrichmentDataIndex );
                moris_index tFirstSpgIndex = mXTKModelPtr->get_subphase_group_index( tFirstSpgId, aEnrichmentDataIndex );

                // go through enriched BFs associated with current non-enriched BF and see which one of them contains the SPG
                bool tFound = false;
                for ( uint iEnrLvl = 0; iEnrLvl < tEnrichmentData->mBasisEnrichmentIndices( tNonEnrBfIndex ).numel(); iEnrLvl++ )
                {
                    // check if the SP checked for is associated with the current enriched BF ...
                    bool tSpgIsInSupport = this->subphase_group_is_in_support(
                            aEnrichmentDataIndex,
                            tFirstSpgIndex,
                            tEnrichmentData->mBasisEnrichmentIndices( tNonEnrBfIndex )( iEnrLvl ) );

                    // ... if so, the ID for the enr. BF is known
                    if ( tSpgIsInSupport )
                    {
                        // get the ID for the current enr. BF
                        moris_index tEnrBfIndex = tEnrichmentData->mBasisEnrichmentIndices( tNonEnrBfIndex )( iEnrLvl );
                        moris_index tEnrBfId    = tEnrichmentData->mEnrichedBasisIndexToId( tEnrBfIndex );

                        // make sure that the enr. BF found is unique
                        MORIS_ERROR( !tFound,
                                "Enrichment::assign_enriched_coefficients_identifiers_new() - "
                                "Data indicates that two different enriched BFs constructed from the same base BF have support in the same SPG. "
                                "Something must have gone wrong in the enrichment." );

                        // store ID answer
                        aEnrBfIds( iProc )( iEnrBF ) = tEnrBfId;

                        // mark BF as being found
                        tFound = true;
                    }
                }

                // check that the enriched BF looked for has been found
                MORIS_ERROR( tFound,
                        "Enrichment::assign_enriched_coefficients_identifiers_new() - "
                        "Proc #%i requested the enriched BF ID for the basis interpolating into SPG ID %i from base BF %i. "
                        "But no such enriched BF has been found.",
                        tCommTable( iProc ),
                        tFirstSpgId,
                        tNonEnrBfId );

            }    // end for: communicated entities from current proc

        }        // end for: communication list for each processor

    }            // end function: Enrichment::prepare_answers_for_owned_enriched_coefficient_IDs()

    //-------------------------------------------------------------------------------------

    void
    Enrichment::handle_requested_unzipped_enriched_coefficient_answers(
            moris_index const &                 aEnrichmentDataIndex,
            Cell< Cell< moris_index > > const & aNotOwnedEnrBfsToProcs,
            Cell< Matrix< IdMat > > const &     aReceivedEnrBfIds )
    {
        // access current enrichment data
        xtk::Enrichment_Data* tEnrichmentData = &mEnrichmentData( aEnrichmentDataIndex );

        // process answers from each proc communicated with
        for ( uint iProc = 0; iProc < aReceivedEnrBfIds.size(); iProc++ )
        {
            // get the number of requests and answers from the current proc
            uint tNumReceivedAns = aReceivedEnrBfIds( iProc ).numel();

            // make sure that all requests have been answered
            MORIS_ASSERT( tNumReceivedAns == aNotOwnedEnrBfsToProcs( iProc ).size(),
                    "Enrichment::handle_requested_unzipped_enriched_coefficient_answers() - "
                    "Number of answers received different from number of requests sent." );

            // assign IDs for each communicated entity
            for ( uint iEnrBf = 0; iEnrBf < tNumReceivedAns; iEnrBf++ )
            {
                // get the local enr. basis index and corresponding ID from the communicated information
                moris_index tEnrBfIndex = aNotOwnedEnrBfsToProcs( iProc )( iEnrBf );
                moris_index tEnrBfId    = aReceivedEnrBfIds( iProc )( iEnrBf );

                // check that a basis-ID doesn't get assigned twice
                MORIS_ASSERT( mEnrichmentData( aEnrichmentDataIndex ).mEnrichedBasisIndexToId( tEnrBfIndex ) == MORIS_ID_MAX,
                        "Enrichment::handle_requested_unzipped_enriched_coefficient_answers() - "
                        "ID already set for this basis function. This shouldn't happen." );

                // assign enr. basis ID to proc local enr. basis index
                tEnrichmentData->mEnrichedBasisIndexToId( tEnrBfIndex ) = tEnrBfId;
            }

        }    // end for: each proc communicated with

    }        // end function: Enrichment::handle_requested_unzipped_enriched_coefficient_answers()

    //-------------------------------------------------------------------------------------

    void
    Enrichment::communicate_basis_information_with_owner(
            moris_index const &                       aEnrichmentDataIndex,      // B-spline mesh index
            Cell< Cell< moris_index > > const &       aBasisIdToBasisOwner,      // request lists: outer cell: proc index in comm table to request from, inner cell: list of non-enriched Basis IDs to request
            Cell< Cell< moris_index > > const &       aSubphaseIdInSupport,      // inner cell: list of SP IDs for unique identification purposes
            Cell< moris_index > const &               aProcRanks,                // comm table: proc IDs associated with outer cell indices in above lists
            std::unordered_map< moris_id, moris_id >& aProcRankToIndexInData,    // map corresponding to comm table
            Cell< Matrix< IndexMat > >&               aEnrichedBasisId )                       // output: list of Basis IDs as requested above
    {
        // STEP 0: prepare send information // copy into a matrix
        Cell< Matrix< IdMat > > tBasisIdToBasisOwnerMat( aProcRanks.size() );
        Cell< Matrix< IdMat > > tSubphaseIdInSupport( aProcRanks.size() );

        // copy all information for each processor into matrices, instead of cells (for easier communication later)
        for ( uint iProc = 0; iProc < aBasisIdToBasisOwner.size(); iProc++ )
        {
            tBasisIdToBasisOwnerMat( iProc ) = Matrix< IdMat >( 1, aBasisIdToBasisOwner( iProc ).size() );
            tSubphaseIdInSupport( iProc )    = Matrix< IdMat >( 1, aBasisIdToBasisOwner( iProc ).size() );

            for ( uint jBF = 0; jBF < aBasisIdToBasisOwner( iProc ).size(); jBF++ )
            {
                tBasisIdToBasisOwnerMat( iProc )( jBF ) = aBasisIdToBasisOwner( iProc )( jBF );
                tSubphaseIdInSupport( iProc )( jBF )    = aSubphaseIdInSupport( iProc )( jBF );
            }

            if ( aBasisIdToBasisOwner( iProc ).size() == 0 )
            {
                tBasisIdToBasisOwnerMat( iProc ).resize( 1, 1 );
                tBasisIdToBasisOwnerMat( iProc )( 0 ) = MORIS_INDEX_MAX;

                tSubphaseIdInSupport( iProc ).resize( 1, 1 );
                tSubphaseIdInSupport( iProc )( 0 ) = MORIS_INDEX_MAX;
            }
        }

        // convert proc ranks from cell to matrix
        Matrix< IdMat > tCommTable( 1, aProcRanks.size() );
        for ( uint iProc = 0; iProc < aProcRanks.size(); iProc++ )
        {
            tCommTable( iProc ) = aProcRanks( iProc );
        }

        // initialize arrays in which to receive information
        Cell< Matrix< IdMat > > tReceiveInfoBasisId;
        Cell< Matrix< IdMat > > tReceivedSubphaseId;

        // send and receive basis information
        moris::communicate_mats( tCommTable, tBasisIdToBasisOwnerMat, tReceiveInfoBasisId );
        moris::communicate_mats( tCommTable, tSubphaseIdInSupport, tReceivedSubphaseId );

        // initialize answer array
        Cell< Matrix< IndexMat > > tEnrichedBasisIds( aProcRanks.size() );

        // iterate through received information and setup response information
        for ( uint iProc = 0; iProc < aProcRanks.size(); iProc++ )
        {
            // initialize answer array
            tEnrichedBasisIds( iProc ).resize( tReceiveInfoBasisId( iProc ).n_rows(), tReceiveInfoBasisId( iProc ).n_cols() );

            // initialize counters
            moris_id tBasisId    = 0;
            moris_id tBasisIndex = 0;
            uint     tCount      = 0;

            // check that this isn't just an empty dummy request
            if ( tReceiveInfoBasisId( iProc )( 0 ) != MORIS_INDEX_MAX )
            {
                // loop through the bases whose IDs are to be assigned
                for ( uint jBF = 0; jBF < tReceiveInfoBasisId( iProc ).n_cols(); jBF++ )
                {
                    // get the index of the requested non-enriched BF wrt. to the current (i.e. owning) processor
                    // if the received non-enriched basis ID is not the same as the previous, then ...
                    if ( tBasisId != tReceiveInfoBasisId( iProc )( jBF ) )
                    {
                        // ... set current requested non-enriched basis ID
                        tBasisId = tReceiveInfoBasisId( iProc )( jBF );

                        // ... and get the index on the requested non-enriched basis index on the current proc
                        tBasisIndex = mBackgroundMeshPtr->get_loc_entity_ind_from_entity_glb_id(
                                tBasisId,
                                mBasisRank,
                                aEnrichmentDataIndex );
                    }

                    // get the subphase index of the currently requested Subphase
                    moris_id    tSubphaseId    = tReceivedSubphaseId( iProc )( jBF );
                    moris_index tSubphaseIndex = mXTKModelPtr->get_subphase_index( tSubphaseId );

                    // iterate through the max integer ids and match
                    bool tFound = false;

                    // go through enriched BFs associated with current non-enriched BF
                    for ( uint k = 0; k < mEnrichmentData( aEnrichmentDataIndex ).mBasisEnrichmentIndices( tBasisIndex ).numel(); k++ )
                    {
                        // check if the SP checked for is associated with the current enriched BF ...
                        bool tSpIsInSupport = this->subphase_is_in_support(
                                aEnrichmentDataIndex,
                                tSubphaseIndex,
                                mEnrichmentData( aEnrichmentDataIndex ).mBasisEnrichmentIndices( tBasisIndex )( k ) );

                        // ... if so, the ID for the enr. BF is known
                        if ( tSpIsInSupport )
                        {
                            moris_index tEnrichedBasisIndex = mEnrichmentData( aEnrichmentDataIndex ).mBasisEnrichmentIndices( tBasisIndex )( k );
                            moris_index tEnrichedBasisId    = mEnrichmentData( aEnrichmentDataIndex ).mEnrichedBasisIndexToId( tEnrichedBasisIndex );

                            tEnrichedBasisIds( iProc )( tCount ) = tEnrichedBasisId;
                            tCount++;
                            tFound = true;
                        }
                    }

                    // check that the enr. BF ID has actually been found
                    MORIS_ERROR( tFound,
                            "Enrichment::communicate_basis_information_with_owner() - Basis %i not found on proc %i",
                            tReceiveInfoBasisId( iProc )( jBF ),
                            aProcRanks( iProc ) );

                }    // end for: each request for the current processor

                // check that the number of found Basis IDs matches the number of requested Basis IDs and nothing was left out or answered twice
                MORIS_ASSERT( tReceiveInfoBasisId( iProc ).numel() == tCount,
                        "Enrichment::communicate_basis_information_with_owner() - "
                        "Number of enriched basis ID requests answered by proc %i not equal to number of basis IDs requested.",
                        aProcRanks( iProc ) );
            }

            // otherwise, answer dummy request with dummy info
            else
            {
                tEnrichedBasisIds( iProc )( 0 ) = MORIS_INDEX_MAX;
            }

        }    // end for: each proc communicated with

        // send and receive enr. basis IDs
        moris::communicate_mats( tCommTable, tEnrichedBasisIds, aEnrichedBasisId );
    }

    //-------------------------------------------------------------------------------------

    void
    Enrichment::set_received_enriched_basis_ids(
            moris_index const &                 aEnrichmentDataIndex,
            Cell< Matrix< IndexMat > > const &  aReceivedEnrichedIds,
            Cell< Cell< moris_index > > const & aBasisIndexToBasisOwner )
    {
        for ( uint iProc = 0; iProc < aReceivedEnrichedIds.size(); iProc++ )
        {
            if ( aReceivedEnrichedIds( iProc )( 0 ) != MORIS_INDEX_MAX )
            {
                // check size of received data
                MORIS_ASSERT( aReceivedEnrichedIds( iProc ).numel() == aBasisIndexToBasisOwner( iProc ).size(),
                        "Enrichment::set_received_enriched_basis_ids() - Dimension mismatch between received information and expected information" );

                // loop over the received Basis-Ids for every proc
                for ( uint jBasisId = 0; jBasisId < aReceivedEnrichedIds( iProc ).numel(); jBasisId++ )
                {
                    // get the local enr. basis index and corresponding ID from the communicated information
                    moris_index tLocalBasisIndex = aBasisIndexToBasisOwner( iProc )( jBasisId );
                    moris_index tGlobalId        = aReceivedEnrichedIds( iProc )( jBasisId );

                    // check that a basis-ID doesn't get assigned twice
                    MORIS_ASSERT( mEnrichmentData( aEnrichmentDataIndex ).mEnrichedBasisIndexToId( tLocalBasisIndex ) == MORIS_INDEX_MAX,
                            "Enrichment::set_received_enriched_basis_ids() - Id already set for this basis function" );

                    // assign enr. basis ID to proc local enr. basis index
                    mEnrichmentData( aEnrichmentDataIndex ).mEnrichedBasisIndexToId( tLocalBasisIndex ) = tGlobalId;
                }
            }
        }
    }

    //-------------------------------------------------------------------------------------

    moris::size_t
    Enrichment::count_elements_in_support( Matrix< IndexMat > const & aParentElementsInSupport )
    {

        // Number of elements in this support (need both parent and total)
        moris::size_t tNumParentElementsInSupport = aParentElementsInSupport.n_cols();
        moris::size_t tNumElementsInSupport       = 0;

        // Count children elements in support
        for ( moris::size_t i = 0; i < tNumParentElementsInSupport; i++ )
        {
            Cell< moris_index > const & tSubphaseIndices = mCutIgMesh->get_parent_cell_subphases( aParentElementsInSupport( i ) );

            // iterate through subphases
            for ( uint iSP = 0; iSP < tSubphaseIndices.size(); iSP++ )
            {
                tNumElementsInSupport = tNumElementsInSupport + mCutIgMesh->get_subphase_ig_cells( tSubphaseIndices( iSP ) )->mIgCellGroup.size();
            }
        }

        return tNumElementsInSupport;
    }

    //-------------------------------------------------------------------------------------

    bool
    Enrichment::subphase_is_in_support(
            moris_index const & aEnrichmentDataIndex,
            moris_index         aSubphaseIndex,
            moris_index         aEnrichedBasisIndex )
    {
        // get list of
        Matrix< IndexMat >& tSPsInEnrBasis = mEnrichmentData( aEnrichmentDataIndex ).mSubphaseIndsInEnrichedBasis( aEnrichedBasisIndex );

        // get the number of SPs in support of curent enr. BF
        uint tNumSPsInSupport = tSPsInEnrBasis.numel();

        // go through all subphases in support and see if one matches the one to be checked for
        for ( uint iSpInSupport = 0; iSpInSupport < tNumSPsInSupport; iSpInSupport++ )
        {
            // get the current SP's index
            moris_index tSpIndex = tSPsInEnrBasis( iSpInSupport );

            // check if SP is the one we're looking for
            if ( tSpIndex == aSubphaseIndex )
            {
                return true;
            }

            // NOTE [NW]: I believe the below comment is not true. Why was this case specified in the first place? Leaving it commented out in here for now...
            // SP indices in support are ordered in ascending order, so if the one we're looking at is already greater than the one we're looking for, it will not be found
            else if ( tSpIndex > aSubphaseIndex )
            {
                return false;
            }
        }

        // if SP not found until here, return not found
        return false;
    }

    //-------------------------------------------------------------------------------------

    bool
    Enrichment::subphase_group_is_in_support(
            moris_index const & aBsplineMeshListIndex,
            moris_index         aSubphaseGroupIndex,
            moris_index         aEnrichedBasisIndex )
    {
        // get list of SPGs in the enriched basis' support
        Matrix< IndexMat >& tSPGsInEnrBasis = mEnrichmentData( aBsplineMeshListIndex ).mSubphaseGroupIndsInEnrichedBasis( aEnrichedBasisIndex );

        // get the number of SPGs in support of curent enr. BF
        uint tNumSPGsInSupport = tSPGsInEnrBasis.numel();

        // go through all subphase groups in support and see if one matches the one to be checked for
        for ( uint iSpgInSupport = 0; iSpgInSupport < tNumSPGsInSupport; iSpgInSupport++ )
        {
            // get the current SPG's index
            moris_index tSpgIndex = tSPGsInEnrBasis( iSpgInSupport );

            // check if SPG is the one we're looking for
            if ( tSpgIndex == aSubphaseGroupIndex )
            {
                return true;
            }
        }

        // if SPG not found until here, return not found
        return false;
    }

    //-------------------------------------------------------------------------------------

    void
    Enrichment::print_basis_support_debug(
            moris_index                aBasisIndex,
            Matrix< IndexMat > const & aParentElementsInSupport,
            Matrix< IndexMat > const & aSubphasesInSupport,
            IndexMap&                  aSubPhaseIndexToSupportIndex,
            Matrix< IndexMat > const & aPrunedSubPhaseToSubphase,
            Matrix< IndexMat >&        aSubPhaseBinEnrichmentVals )
    {
        std::cout << "--------------------------------------------------" << std::endl;
        std::cout << "Basis Index: " << aBasisIndex << std::endl;
        std::cout << "Parent Cell In Support:";
        for ( uint i = 0; i < aParentElementsInSupport.numel(); i++ )
        {
            std::cout << std::setw( 8 ) << aParentElementsInSupport( i );
        }
        std::cout << "\nSubphases In Support:";
        for ( uint i = 0; i < aSubphasesInSupport.numel(); i++ )
        {
            std::cout << std::setw( 8 ) << aSubphasesInSupport( i );
        }

        std::cout << "\nSubphase Neighborhood In Support:" << std::endl;
        for ( uint i = 0; i < aPrunedSubPhaseToSubphase.n_rows(); i++ )
        {
            std::cout << std::setw( 6 ) << aSubphasesInSupport( i ) << " | ";

            for ( uint j = 0; j < aPrunedSubPhaseToSubphase.n_cols(); j++ )
            {
                if ( aPrunedSubPhaseToSubphase( i, j ) != MORIS_INDEX_MAX )
                {
                    std::cout << std::setw( 6 ) << aSubphasesInSupport( aPrunedSubPhaseToSubphase( i, j ) );
                }
            }
            std::cout << std::endl;
        }

        std::cout << "\nSubphase Enrichment Level: \n";
        for ( uint i = 0; i < aSubPhaseBinEnrichmentVals.numel(); i++ )
        {
            std::cout << std::setw( 8 ) << aSubPhaseBinEnrichmentVals( i );
        }
        std::cout << "\n--------------------------------------------------" << std::endl;
    }

    //-------------------------------------------------------------------------------------

    void
    Enrichment::construct_enriched_interpolation_mesh()
    {
        // log/trace this function
        Tracer tTracer( "XTK", "Enrichment", "Construct Enriched Interpolation Mesh" );

        // initialize a new enriched interpolation mesh
        mXTKModelPtr->mEnrichedInterpMesh( 0 ) = new Enriched_Interpolation_Mesh( mXTKModelPtr );

        // set enriched basis rank
        mXTKModelPtr->mEnrichedInterpMesh( 0 )->mBasisRank = mBasisRank;

        // set mesh index
        mXTKModelPtr->mEnrichedInterpMesh( 0 )->mMeshIndices = mMeshIndices;
        mXTKModelPtr->mEnrichedInterpMesh( 0 )->setup_mesh_index_map();

        // allocate memory for enriched interpolation cells
        this->allocate_interpolation_cells();

        // unzip all IP cells and vertices and construct the whole enr. IP mesh
        /* Note: This constructs all unzipped (i.e. enriched) IP cells with unzipped interpolation vertices being attached to a single cell
         * this also handles the case of multiple enrichments where the number of interpolation vertices vary */
        this->construct_enriched_interpolation_vertices_and_cells();

        // TODO: check that there are no empty enriched T-matrices
        // MORIS_ASSERT( this->check_for_empty_T_matrices(),
        //         "Enrichment::construct_enriched_interpolation_mesh() - "
        //         "Check for emtpy T-matrices failed. There should be a error statement previous to this." );

        mXTKModelPtr->mEnrichedInterpMesh( 0 )->mCoeffToEnrichCoeffs.resize( mMeshIndices.max() + 1 );
        mXTKModelPtr->mEnrichedInterpMesh( 0 )->mEnrichCoeffLocToGlob.resize( mMeshIndices.max() + 1 );

        // enriched to non-enriched coefficient mapping for all B-spline meshes
        for ( uint iMesh = 0; iMesh < mMeshIndices.numel(); iMesh++ )
        {
            // get the discretization mesh index (DMI)
            moris_index tMeshIndex = mMeshIndices( iMesh );

            // add the coeff to enriched coeffs to enriched interpolation mesh
            mXTKModelPtr->mEnrichedInterpMesh( 0 )->mCoeffToEnrichCoeffs( tMeshIndex ) =
                    mEnrichmentData( tMeshIndex ).mBasisEnrichmentIndices;

            // add the local to global map
            mXTKModelPtr->mEnrichedInterpMesh( 0 )->mEnrichCoeffLocToGlob( tMeshIndex ) =
                    mEnrichmentData( tMeshIndex ).mEnrichedBasisIndexToId;
        }

        // tell the enriched IP mesh to finish setting itself up
        mXTKModelPtr->mEnrichedInterpMesh( 0 )->finalize_setup();

        // in most cases all the interpolation vertices are the same. We merge them back together with this call
        // post-processing to construct_enriched_interpolation_vertices_and_cells in an effort to not add complexity to the function
        // (as that function is already too complex/loaded)
        mXTKModelPtr->mEnrichedInterpMesh( 0 )->merge_duplicate_interpolation_vertices();

        // reset global to local maps (delete & setup again) with deleted duplicate vertices
        mXTKModelPtr->mEnrichedInterpMesh( 0 )->mGlobalToLocalMaps( 0 ).clear();
        mXTKModelPtr->mEnrichedInterpMesh( 0 )->setup_vertex_maps();
    }

    //-------------------------------------------------------------------------------------

    void
    Enrichment::construct_enriched_interpolation_mesh_new()
    {
        // log/trace this function
        Tracer tTracer( "XTK", "Enrichment", "Construct Enriched Interpolation Mesh" );

        // initialize a new enriched interpolation mesh
        mXTKModelPtr->mEnrichedInterpMesh( 0 ) = new Enriched_Interpolation_Mesh( mXTKModelPtr );

        // set enriched basis rank
        mXTKModelPtr->mEnrichedInterpMesh( 0 )->mBasisRank = mBasisRank;

        // set mesh index
        mXTKModelPtr->mEnrichedInterpMesh( 0 )->mMeshIndices = mMeshIndices;
        mXTKModelPtr->mEnrichedInterpMesh( 0 )->setup_mesh_index_map();

        // allocate memory for enriched interpolation cells
        // this->allocate_interpolation_cells_based_on_SPGs();
        this->allocate_interpolation_cells_based_on_SPGs_new();

        // construct inverse of the map relating base IP cells and their unzipping to the resulting UIPC index
        this->construct_UIPC_to_unzipping_index();

        // unzip all IP cells and vertices and construct the whole enr. IP mesh
        // Note: This constructs all unzipped (i.e. enriched) IP cells with unzipped interpolation vertices being attached to
        // note: a single cell this also handles the case of multiple enrichments where the number of interpolation vertices vary
        // this->construct_enriched_interpolation_vertices_and_cells_based_on_SPGs();
        this->construct_enriched_interpolation_vertices_and_cells_based_on_SPGs_new();

        // assign IDs to all UIPCs and communicate them across all procs
        this->communicate_unzipped_ip_cells();

        mXTKModelPtr->mEnrichedInterpMesh( 0 )->mCoeffToEnrichCoeffs.resize( mMeshIndices.max() + 1 );
        mXTKModelPtr->mEnrichedInterpMesh( 0 )->mEnrichCoeffLocToGlob.resize( mMeshIndices.max() + 1 );

        // enriched to non-enriched coefficient mapping for all B-spline meshes
        for ( uint iMesh = 0; iMesh < mMeshIndices.numel(); iMesh++ )
        {
            // get the discretization mesh index (DMI)
            moris_index tMeshIndex = mMeshIndices( iMesh );

            // add the coeff to enriched coeffs to enriched interpolation mesh
            mXTKModelPtr->mEnrichedInterpMesh( 0 )->mCoeffToEnrichCoeffs( tMeshIndex ) =
                    mEnrichmentData( tMeshIndex ).mBasisEnrichmentIndices;

            // add the local to global map
            mXTKModelPtr->mEnrichedInterpMesh( 0 )->mEnrichCoeffLocToGlob( tMeshIndex ) =
                    mEnrichmentData( tMeshIndex ).mEnrichedBasisIndexToId;
        }

        // tell the enriched IP mesh to finish setting itself up
        mXTKModelPtr->mEnrichedInterpMesh( 0 )->finalize_setup_new();

        // in most cases all the interpolation vertices are the same. We merge them back together with this call
        // post-processing to construct_enriched_interpolation_vertices_and_cells in an effort to not add complexity to the function
        // (as that function is already too complex/loaded)
        mXTKModelPtr->mEnrichedInterpMesh( 0 )->merge_duplicate_interpolation_vertices();

        // reset global to local maps (delete & setup again) with deleted duplicate vertices
        mXTKModelPtr->mEnrichedInterpMesh( 0 )->mGlobalToLocalMaps( 0 ).clear();
        mXTKModelPtr->mEnrichedInterpMesh( 0 )->setup_vertex_maps();
    }

    //-------------------------------------------------------------------------------------

    void
    Enrichment::construct_enriched_integration_mesh()
    {
        MORIS_ASSERT( mXTKModelPtr->mEnrichedInterpMesh( 0 ) != nullptr,
                "Enrichment::construct_enriched_integration_mesh() - No enriched interpolation mesh to link enriched integration mesh to" );

        mXTKModelPtr->mEnrichedIntegMesh( 0 ) = new Enriched_Integration_Mesh( mXTKModelPtr );
    }

    //-------------------------------------------------------------------------------------

    void
    Enrichment::construct_enriched_integration_mesh( const Matrix< IndexMat > aBsplineMeshIndices )
    {
        MORIS_ASSERT( mXTKModelPtr->mEnrichedInterpMesh( 0 ) != nullptr,
                "Enrichment::construct_enriched_integration_mesh_new() - No enriched interpolation mesh to link enriched integration mesh to" );

        mXTKModelPtr->mEnrichedIntegMesh( 0 ) = new Enriched_Integration_Mesh( mXTKModelPtr, aBsplineMeshIndices );
    }

    //-------------------------------------------------------------------------------------

    void
    Enrichment::allocate_interpolation_cells()
    {
        // log trace this function
        Tracer tTracer( "XTK", "Enrich", "Allocate Interpolation Cells" );

        // get the enriched interpolation mesh pointer
        Enriched_Interpolation_Mesh* tEnrInterpMesh = mXTKModelPtr->mEnrichedInterpMesh( 0 );

        // set discretization mesh indices
        tEnrInterpMesh->mMeshIndices = mMeshIndices;

        // count how many subphases there are total in the entire mesh
        uint tNumSubphases = mCutIgMesh->get_num_subphases();

        // there is one interpolation cell per subphase
        tEnrInterpMesh->mEnrichedInterpCells.resize( tNumSubphases );

        // assuming all interpolation cells are the same
        // figure out how many vertices there are per interpolation cell
        uint tNumVertsPerCell = 0;

        if ( mBackgroundMeshPtr->get_num_elems() > 0 )
        {
            tNumVertsPerCell = mBackgroundMeshPtr->get_mtk_cell( 0 ).get_number_of_vertices();
        }

        tEnrInterpMesh->mNumVertsPerInterpCell = tNumVertsPerCell;

        // allocate maximum number of enriched vertices
        tEnrInterpMesh->mEnrichedInterpVerts.resize( tNumVertsPerCell * tNumSubphases );

        // allocate the base vertices to vertex enrichment data
        tEnrInterpMesh->mBaseInterpVertToVertEnrichmentIndex.resize( mMeshIndices.max() + 1, mBackgroundMeshPtr->get_num_nodes() );

        // allocate space in the vertex enrichment index to parent vertex enrichment data
        tEnrInterpMesh->mVertexEnrichmentParentVertexIndex.resize( mMeshIndices.max() + 1 );

        tEnrInterpMesh->mInterpVertEnrichment.resize( mMeshIndices.max() + 1 );

        // allocate base cell to enriched cell data
        tEnrInterpMesh->mBaseCellToEnrichedCell.resize( mBackgroundMeshPtr->get_num_elems() );

        // initialize subphase index to UIPC index map, for the SP based enrichment this map is trivial
        uint tNumSPs = mCutIgMesh->get_num_subphases();
        mSubphaseIndexToEnrIpCellIndex.resize( tNumSPs );
        for ( uint iSP = 0; iSP < tNumSPs; iSP++ )
        {
            mSubphaseIndexToEnrIpCellIndex( iSP ) = iSP;
        }
    }

    //-------------------------------------------------------------------------------------

    void
    Enrichment::allocate_interpolation_cells_based_on_SPGs()
    {
        // log trace this function
        Tracer tTracer( "XTK", "Enrich", "Allocate Interpolation Cells" );

        // get the enriched interpolation mesh pointer
        Enriched_Interpolation_Mesh* tEnrInterpMesh = mXTKModelPtr->mEnrichedInterpMesh( 0 );

        // set discretization mesh indices
        tEnrInterpMesh->mMeshIndices = mMeshIndices;

        // initialize counter tracking the number of enriched IP cells that need to be created
        uint tEnrIpCellCounter = 0;

        // get number of IP cells
        uint tNumIpElems = mBackgroundMeshPtr->get_num_elems();

        // initialize list storing number each IP cell needs to be unzipped
        mEnrIpCellIndices.resize( tNumIpElems );

        // initialize lists storing SPGs for unzipping of IP cells
        uint tNumBspMeshes = mMeshIndices.numel();
        mMaterialSpgsUnzippedOnIpCell.resize( tNumBspMeshes );
        mVoidSpgsUnzippedOnIpCell.resize( tNumBspMeshes );

        // initialize map relating SPGs to unzipping of a given IP cell
        mBaseIpCellAndSpgToUnzipping.resize( tNumBspMeshes );

        // initialize map relating SPGs and UIPCs
        mSpgToUipcIndex.resize( tNumBspMeshes );

        // initialize subphase to UIPC map
        mSubphaseIndexToEnrIpCellIndex.resize( mCutIgMesh->get_num_subphases() );

        // initialize lists with correct size for each B-spline mesh
        for ( moris::size_t iBspMesh = 0; iBspMesh < tNumBspMeshes; iBspMesh++ )
        {
            mMaterialSpgsUnzippedOnIpCell( iBspMesh ).resize( tNumIpElems );
            mVoidSpgsUnzippedOnIpCell( iBspMesh ).resize( tNumIpElems );
            mBaseIpCellAndSpgToUnzipping( iBspMesh ).resize( tNumIpElems );

            Bspline_Mesh_Info* tBsplineMeshInfo = mBsplineMeshInfos( iBspMesh );
            mSpgToUipcIndex( iBspMesh ).resize( tBsplineMeshInfo->get_num_SPGs() );
        }

        // loop to obtain the number of max possible number of enr. IP cells
        for ( moris::size_t iIpCell = 0; iIpCell < tNumIpElems; iIpCell++ )
        {
            // get number of times this IP cell needs to be unzipped
            uint tNumUnzippingsOfIpCell = this->maximum_number_of_unzippings_for_IP_cell( iIpCell );

            // correct size for list of enr. IP cell indices on current Ip cell
            mEnrIpCellIndices( iIpCell ).resize( tNumUnzippingsOfIpCell );

            // get the number of SPs on the current IP cell
            Cell< moris_index > const & tSPsOnCell = mCutIgMesh->get_parent_cell_subphases( iIpCell );

            for ( uint iEnrLvl = 0; iEnrLvl < tNumUnzippingsOfIpCell; iEnrLvl++ )
            {
                // fill map relating IP-cell index and local enr. lvl to the index of the resulting enr. IP cell
                mEnrIpCellIndices( iIpCell )( iEnrLvl ) = (moris_index)tEnrIpCellCounter;

                // assemble subphase to UIPC map
                if ( iEnrLvl < tSPsOnCell.size() )
                {
                    // get the index of the subphase associated with the current non-void UIPC
                    moris_index tSpIndex = tSPsOnCell( iEnrLvl );

                    // link the Subphase index to the current enr. IP cell index
                    mSubphaseIndexToEnrIpCellIndex( tSpIndex ) = tEnrIpCellCounter;
                }

                // count total number of enr. IP cells
                tEnrIpCellCounter++;
            }

            // convert the local unzipping index in mBaseIpCellAndSpgToUnzipping to the global enr. IP cell index
            for ( uint iBspMesh = 0; iBspMesh < mBaseIpCellAndSpgToUnzipping.size(); iBspMesh++ )
            {
                for ( uint iSPG = 0; iSPG < mBaseIpCellAndSpgToUnzipping( iBspMesh )( iIpCell ).size(); iSPG++ )
                {
                    // replace unzipping index with enr. IP cell index
                    moris_index tUnzippingIndex                                 = mBaseIpCellAndSpgToUnzipping( iBspMesh )( iIpCell )( iSPG );
                    mBaseIpCellAndSpgToUnzipping( iBspMesh )( iIpCell )( iSPG ) = mEnrIpCellIndices( iIpCell )( tUnzippingIndex );
                }
            }
        }

        // store number of Enr IP cells
        mNumEnrIpCells = (moris_index)tEnrIpCellCounter;

        // there is one interpolation cell per subphase
        tEnrInterpMesh->mEnrichedInterpCells.resize( mNumEnrIpCells );

        // initialize reverse map relating SPGs and UIPCs
        mUipcToSpgIndex.resize( tNumBspMeshes );
        for ( moris::size_t iBspMesh = 0; iBspMesh < tNumBspMeshes; iBspMesh++ )
        {
            mUipcToSpgIndex( iBspMesh ).resize( mNumEnrIpCells );
        }

        // assuming all interpolation cells are the same
        // figure out how many vertices there are per interpolation cell
        uint tNumVertsPerCell = 0;

        // initialize lists relating SPGs and UIPCs
        if ( mBackgroundMeshPtr->get_num_elems() > 0 )
        {
            tNumVertsPerCell = mBackgroundMeshPtr->get_mtk_cell( 0 ).get_number_of_vertices();
        }

        tEnrInterpMesh->mNumVertsPerInterpCell = tNumVertsPerCell;

        // allocate maximum number of enriched vertices
        tEnrInterpMesh->mEnrichedInterpVerts.resize( tNumVertsPerCell * mNumEnrIpCells );

        // allocate the base vertices to vertex enrichment data
        tEnrInterpMesh->mBaseInterpVertToVertEnrichmentIndex.resize( mMeshIndices.max() + 1, mBackgroundMeshPtr->get_num_nodes() );

        // allocate space in the vertex enrichment index to parent vertex enrichment data
        tEnrInterpMesh->mVertexEnrichmentParentVertexIndex.resize( mMeshIndices.max() + 1 );

        tEnrInterpMesh->mInterpVertEnrichment.resize( mMeshIndices.max() + 1 );

        // allocate base cell to enriched cell data
        tEnrInterpMesh->mBaseCellToEnrichedCell.resize( mBackgroundMeshPtr->get_num_elems() );
    }

    //-------------------------------------------------------------------------------------

    void
    Enrichment::allocate_interpolation_cells_based_on_SPGs_new()
    {
        // log trace this function
        Tracer tTracer( "XTK", "Enrich", "Allocate Interpolation Cells" );

        // get the enriched interpolation mesh pointer
        Enriched_Interpolation_Mesh* tEnrInterpMesh = mXTKModelPtr->mEnrichedInterpMesh( 0 );

        // set discretization mesh indices
        tEnrInterpMesh->mMeshIndices = mMeshIndices;

        // initialize counter tracking the number of enriched IP cells that need to be created
        uint tEnrIpCellCounter = 0;

        // get number of base IP cells
        uint tNumIpElems = mBackgroundMeshPtr->get_num_elems();

        // initialize list storing number each IP cell needs to be unzipped
        mEnrIpCellIndices.resize( tNumIpElems );

        // initialize lists storing SPGs for unzipping of IP cells
        uint tNumBspMeshes = mMeshIndices.numel();
        mMaterialSpgsUnzippedOnIpCell.resize( tNumBspMeshes );
        mVoidSpgsUnzippedOnIpCell.resize( tNumBspMeshes );

        // initialize map relating SPGs to unzipping of a given IP cell
        mBaseIpCellAndSpgToUnzipping.resize( tNumBspMeshes );

        // initialize map relating unzippings to SPGs
        mUnzippingToSpgIndex.resize( tNumBspMeshes );

        // initialize map relating SPGs and UIPCs
        mSpgToUipcIndex.resize( tNumBspMeshes );

        // initialize subphase to UIPC map
        mSubphaseIndexToEnrIpCellIndex.resize( mCutIgMesh->get_num_subphases() );

        /* Maps assembled below:
         * - mEnrIpCellIndices               // input: base IP cell index, index of unzipping || output: index of enr. IP cell
         * - mSubphaseIndexToEnrIpCellIndex  // input: subphase index || output: enriched interpolation cell index
         * - mUnzippingToSpgIndex            // input: B-spline mesh list index, base IP cell index, index of unzipping || output: SPG index the UIPC belongs to
         * - mBaseIpCellAndSpgToUnzipping    // input: B-spline mesh list index, base IP cell index, SPG index local to corresponding B-spline element || output: Enr. IP cell index
         * - mSpgToUipcIndex                 // input: B-spline mesh list index, SPG index || output: List of enr. IP cells / cell clusters on SPG
         */

        // initialize lists with correct size for each B-spline mesh
        for ( uint iBspMesh = 0; iBspMesh < tNumBspMeshes; iBspMesh++ )
        {
            mBaseIpCellAndSpgToUnzipping( iBspMesh ).resize( tNumIpElems );
            mUnzippingToSpgIndex( iBspMesh ).resize( tNumIpElems );

            Bspline_Mesh_Info* tBsplineMeshInfo = mBsplineMeshInfos( iBspMesh );
            mSpgToUipcIndex( iBspMesh ).resize( tBsplineMeshInfo->get_num_SPGs() );
        }

        // loop to obtain the number of max possible number of enr. IP cells
        for ( uint iIpCell = 0; iIpCell < tNumIpElems; iIpCell++ )
        {
            // get the number of SPs on the current IP cell
            const Cell< moris_index >& tSPsOnCell = mCutIgMesh->get_parent_cell_subphases( iIpCell );
            const uint                 tNumSPs    = tSPsOnCell.size();

            // get the Union MSD indices for the current IP cell
            const Cell< moris_index >& tUnionMsdIndices = mCutIgMesh->get_union_MSD_indices_for_base_IP_cell( iIpCell );
            const uint                 tNumVoidUIPCs    = tUnionMsdIndices.size();

            // get number of times this IP cell needs to be unzipped
            uint tNumUnzippingsOfIpCell = tNumSPs + tNumVoidUIPCs;

            // correct size for list of enr. IP cell indices on current Ip cell and their SPGs
            mEnrIpCellIndices( iIpCell ).resize( tNumUnzippingsOfIpCell );

            // count up enriched IP cells to be constructed on the current base IP cell and index them
            for ( uint iUnzipping = 0; iUnzipping < tNumUnzippingsOfIpCell; iUnzipping++ )
            {
                // fill map relating IP-cell index and local enr. lvl to the index of the resulting enr. IP cell
                mEnrIpCellIndices( iIpCell )( iUnzipping ) = (moris_index)tEnrIpCellCounter;

                // assemble subphase to UIPC map
                if ( iUnzipping < tSPsOnCell.size() )
                {
                    // get the index of the subphase associated with the current non-void UIPC
                    moris_index tSpIndex = tSPsOnCell( iUnzipping );

                    // link the Subphase index to the current enr. IP cell index
                    mSubphaseIndexToEnrIpCellIndex( tSpIndex ) = tEnrIpCellCounter;
                }

                // count total number of enr. IP cells
                tEnrIpCellCounter++;
            }

            // construct information on the UIPCs related to B-spline meshes
            for ( uint iBspMesh = 0; iBspMesh < tNumBspMeshes; iBspMesh++ )
            {
                // get the current B-spline mesh info
                Bspline_Mesh_Info* tBsplineMeshInfo = mBsplineMeshInfos( iBspMesh );

                // get the number of SPGs associated with the current IP cell
                uint tNumSpgsAssociatedWithIpCell = tBsplineMeshInfo->get_SPG_indices_associated_with_extraction_cell( iIpCell ).size();

                // resize maps relating UIPCs to their SPGs
                mUnzippingToSpgIndex( iBspMesh )( iIpCell ).resize( tNumUnzippingsOfIpCell );
                mBaseIpCellAndSpgToUnzipping( iBspMesh )( iIpCell ).resize( tNumSpgsAssociatedWithIpCell, -1 );

                // assign the SPGs to the material clusters
                for ( uint iUnzipping = 0; iUnzipping < tNumSPs; iUnzipping++ )
                {
                    // get the material SPG index
                    moris_index tSpgIndex = mCutIgMesh->get_material_SPG_indices_for_base_IP_cell( iBspMesh, iIpCell )( iUnzipping );

                    // store the SPG index for this unzipping
                    mUnzippingToSpgIndex( iBspMesh )( iIpCell )( iUnzipping ) = tSpgIndex;

                    // get the correspond UIPC's index
                    moris_index tUipcIndex = mEnrIpCellIndices( iIpCell )( iUnzipping );

                    // get the local index of this SPG
                    moris_index tLocalSpgIndex = tBsplineMeshInfo->mSubphaseGroups( tSpgIndex )->get_local_index();

                    // store that the current UIPC is associated with this SPG
                    mSpgToUipcIndex( iBspMesh )( tSpgIndex ).push_back( tUipcIndex );
                    mBaseIpCellAndSpgToUnzipping( iBspMesh )( iIpCell )( tLocalSpgIndex ) = tUipcIndex;
                }

                // get the list of void MSD indices and convert it to a map
                const Cell< moris_index >& tVoidMsdIndices = mCutIgMesh->get_void_MSD_indices_for_base_IP_cell( iBspMesh, iIpCell );

                // get the list of corresponding SPGs
                const Cell< moris_index >& tVoidSpgIndices = mCutIgMesh->get_void_SPG_indices_for_base_IP_cell( iBspMesh, iIpCell );

                // get the list of free void MSD indices
                const Cell< moris_index >& tFreeVoidMsdIndices = mCutIgMesh->get_free_void_MSD_indices_for_base_IP_cell( iBspMesh, iIpCell );

                // initialize punch-cards indicating which elements of the above lists have been found
                Cell< bool > tVoidMsdIndicesFound( tVoidMsdIndices.size(), false );
                Cell< bool > tFreeVoidMsdIndicesFound( tFreeVoidMsdIndices.size(), false );

                // go over the union MSD indices (i.e. void clusters to be constructed on the current IP cell)
                for ( uint iUnionMSDI = 0; iUnionMSDI < tUnionMsdIndices.size(); iUnionMSDI++ )
                {
                    // get the index of the unzipping
                    uint tUnzipping = tNumSPs + iUnionMSDI;

                    // get the current Union MSD index
                    moris_index tUnionMsdIndex = tUnionMsdIndices( iUnionMSDI );

                    // find the index in the list of void MSD indices
                    moris_index tPositionInListOfVoidMSDIs = xtk::find_in_cell( tVoidMsdIndices, tUnionMsdIndex, tVoidMsdIndicesFound );

                    // if the MSD index is found among the void MSD indices, note this down
                    if ( tPositionInListOfVoidMSDIs != -1 )
                    {
                        // get the corresponding SPG
                        moris_index tSpgIndex = tVoidSpgIndices( tPositionInListOfVoidMSDIs );

                        // assign the SPG used for enrichment
                        mUnzippingToSpgIndex( iBspMesh )( iIpCell )( tUnzipping ) = tSpgIndex;

                        // get the correspond UIPC's index
                        moris_index tUipcIndex = mEnrIpCellIndices( iIpCell )( tUnzipping );

                        // get the local index of this SPG
                        moris_index tLocalSpgIndex = tBsplineMeshInfo->mSubphaseGroups( tSpgIndex )->get_local_index();

                        // store that the current UIPC is associated with this SPG
                        mSpgToUipcIndex( iBspMesh )( tSpgIndex ).push_back( tUipcIndex );
                        mBaseIpCellAndSpgToUnzipping( iBspMesh )( iIpCell )( tLocalSpgIndex ) = tUipcIndex;
                    }
                    else    // if the MSD index doesn't correspond to one of the SPGs related to the current IP cell
                    {
                        // find the MSD index in the free void MSD indices for the current IP cell
                        moris_index tPositionInListOfFreeVoidMSDIs = xtk::find_in_cell( tFreeVoidMsdIndices, tUnionMsdIndex, tFreeVoidMsdIndicesFound );

                        // check that the index is found, if not something must be wrong
                        MORIS_ERROR( tPositionInListOfFreeVoidMSDIs != -1,
                                "Enrichment::allocate_interpolation_cells_based_on_SPGs_new() - "
                                "Union MSD index #%i for IP cell #%i on B-spline mesh #%i neither found "
                                "in void or free void MSD indices for this cell. There must be a bug.",
                                tUnionMsdIndex,
                                iIpCell,
                                mMeshIndices( iBspMesh ) );

                        // assign SPG -1 to mark need for basis extension
                        mUnzippingToSpgIndex( iBspMesh )( iIpCell )( tUnzipping ) = -1;
                    }

                }    // end for: all Union free MSD indices on IP cell for current B-spline mesh

                // check that the UIPCs on the current base IP cell on the current B-spline mesh have been fully assigned
                MORIS_ERROR( xtk::punch_card_sum( tVoidMsdIndicesFound ) == tVoidMsdIndices.size(),
                        "Enrichment::allocate_interpolation_cells_based_on_SPGs_new() - "
                        "Not all  void MSD indices assigned to an enr. IP cell."
                        "IP cell #%i on B-spline mesh #%i",
                        iIpCell,
                        mMeshIndices( iBspMesh ) );

                MORIS_ERROR( xtk::punch_card_sum( tFreeVoidMsdIndicesFound ) == tFreeVoidMsdIndices.size(),
                        "Enrichment::allocate_interpolation_cells_based_on_SPGs_new() - "
                        "Not all free void MSD indices assigned to an enr. IP cell."
                        "IP cell #%i on B-spline mesh #%i",
                        iIpCell,
                        mMeshIndices( iBspMesh ) );

            }    // end for: all B-spline mesh indices

        }        // end for: all base IP cells

        // store the number of enr. IP cells
        mNumEnrIpCells = (moris_index)tEnrIpCellCounter;

        // there is one interpolation cell per subphase
        tEnrInterpMesh->mEnrichedInterpCells.resize( mNumEnrIpCells );

        // initialize reverse map relating SPGs and UIPCs
        mUipcToSpgIndex.resize( tNumBspMeshes );
        for ( moris::size_t iBspMesh = 0; iBspMesh < tNumBspMeshes; iBspMesh++ )
        {
            mUipcToSpgIndex( iBspMesh ).resize( mNumEnrIpCells );
        }

        // assuming all interpolation cells are the same
        // figure out how many vertices there are per interpolation cell
        uint tNumVertsPerCell = 0;

        // initialize lists relating SPGs and UIPCs
        if ( mBackgroundMeshPtr->get_num_elems() > 0 )
        {
            tNumVertsPerCell = mBackgroundMeshPtr->get_mtk_cell( 0 ).get_number_of_vertices();
        }

        tEnrInterpMesh->mNumVertsPerInterpCell = tNumVertsPerCell;

        // allocate maximum number of enriched vertices
        tEnrInterpMesh->mEnrichedInterpVerts.resize( tNumVertsPerCell * mNumEnrIpCells );

        // allocate the base vertices to vertex enrichment data
        tEnrInterpMesh->mBaseInterpVertToVertEnrichmentIndex.resize( mMeshIndices.max() + 1, mBackgroundMeshPtr->get_num_nodes() );

        // allocate space in the vertex enrichment index to parent vertex enrichment data
        tEnrInterpMesh->mVertexEnrichmentParentVertexIndex.resize( mMeshIndices.max() + 1 );

        tEnrInterpMesh->mInterpVertEnrichment.resize( mMeshIndices.max() + 1 );

        // allocate base cell to enriched cell data
        tEnrInterpMesh->mBaseCellToEnrichedCell.resize( mBackgroundMeshPtr->get_num_elems() );
    }

    //-------------------------------------------------------------------------------------

    uint
    Enrichment::maximum_number_of_unzippings_for_IP_cell( moris_index aIpCellIndex )
    {
        // initialize counter for number of unzippings
        uint tMaxNumUnzippings = 0;

        // get the number of SPs on the current IP cell
        Cell< moris_index > const & tSPsOnCell = mCutIgMesh->get_parent_cell_subphases( aIpCellIndex );
        const uint                  tNumSPs    = tSPsOnCell.size();

        // for each DMI get how often the IP cell needs to be unzipped
        for ( moris::size_t iBspMesh = 0; iBspMesh < mMeshIndices.numel(); iBspMesh++ )
        {
            // set size of list of SPGs with material associated with current IP cell for current B-spline mesh
            mMaterialSpgsUnzippedOnIpCell( iBspMesh )( aIpCellIndex ).resize( tNumSPs );

            // get the pointer to the current B-spline mesh info
            Bspline_Mesh_Info* tBsplineMeshInfo = mBsplineMeshInfos( iBspMesh );

            // get the SPGs that are associated with the current IP cell
            Cell< moris_index > const & tSPGsOnCell =
                    tBsplineMeshInfo->get_SPG_indices_associated_with_extraction_cell( aIpCellIndex );
            uint tNumSPGsOnCell = tSPGsOnCell.size();

            // initialize punch-card logging which SPGs have material on the current IP cell
            Cell< bool > tVoidSPGs( tNumSPGsOnCell, true );

            // initialize list relating each SPG to a specific unzipping
            mBaseIpCellAndSpgToUnzipping( iBspMesh )( aIpCellIndex ).resize( tNumSPGsOnCell );

            // convert List of SPGs to map
            IndexMap tSpgIndexToLocalMap;
            convert_cell_to_map( tSPGsOnCell, tSpgIndexToLocalMap );

            for ( uint iSP = 0; iSP < tNumSPs; iSP++ )
            {
                // get the index of the current subphase
                moris_index tSpIndex = tSPsOnCell( iSP );

                // get the index of SPG the currently treated SP belongs to
                moris_index tSpgIndex = tBsplineMeshInfo->mSpToSpgMap( tSpIndex );

                // mSpgToUipcIndex( iBspMesh )( tSpgIndex ) = tEnrIpCellIndex;

                // store SPG index containing material
                mMaterialSpgsUnzippedOnIpCell( iBspMesh )( aIpCellIndex )( iSP ) = tSpgIndex;

                // find where the SPG is in the list of SPGs on the
                auto tIter = tSpgIndexToLocalMap.find( tSpgIndex );
                MORIS_ASSERT( tIter != tSpgIndexToLocalMap.end(), "Enrichment::maximum_number_of_unzippings_for_IP_cell() - SPG for SP on IP cell not found." );
                moris_index tLocalSpgIndex = tIter->second;

                // mark the SPG as having material in the punch card
                tVoidSPGs( tLocalSpgIndex ) = false;

                // relate the unzipping to the local SPG index, store this information in map
                mBaseIpCellAndSpgToUnzipping( iBspMesh )( aIpCellIndex )( tLocalSpgIndex ) = iSP;
            }

            // count the number of void IP cells that need to be constructed
            uint tNumVoidClusters = xtk::punch_card_sum( tVoidSPGs );

            // set size of list of SPGs without material associated with current IP cell for current B-spline mesh
            mVoidSpgsUnzippedOnIpCell( iBspMesh )( aIpCellIndex ).resize( tNumVoidClusters );

            // store SPG indices for void clusters
            uint tVoidSpgCounter = 0;
            for ( uint iSPG = 0; iSPG < tNumSPGsOnCell; iSPG++ )
            {
                // if SPG has noted to not have material in it
                if ( tVoidSPGs( iSPG ) )
                {
                    // relate the unzipping to the local SPG index, store this information in map
                    mBaseIpCellAndSpgToUnzipping( iBspMesh )( aIpCellIndex )( iSPG ) = tNumSPs + tVoidSpgCounter;

                    mVoidSpgsUnzippedOnIpCell( iBspMesh )( aIpCellIndex )( tVoidSpgCounter ) = tSPGsOnCell( iSPG );
                    tVoidSpgCounter++;
                }
            }

            // get the number of valid Enr. IP cells and clusters on the current IP cell wrt. to the current B-spline mesh
            uint tNumValidClusters = tNumSPs + tNumVoidClusters;

            // overwrite max if new number of unzippings is greater
            tMaxNumUnzippings = std::max( tMaxNumUnzippings, tNumValidClusters );
        }

        // return max number of unzippings
        return tMaxNumUnzippings;
    }

    //-------------------------------------------------------------------------------------

    void
    Enrichment::establish_IP_SPG_SP_relationship( const moris_index aMeshIndex )
    {
        // log/trace this function
        Tracer tTracer( "XTK", "Integration_Mesh_Generator", "Establish relationship between IP cells, SPGs and SPs for mesh index " + std::to_string( aMeshIndex ) );

        // get the pointer to the current B-spline mesh info
        Bspline_Mesh_Info* tBsplineMeshInfo = mBsplineMeshInfos( this->get_list_index_for_mesh_index( aMeshIndex ) );

        // initialize IP-cell to SP map with correct size
        tBsplineMeshInfo->mExtractionCellToSubPhase.resize( tBsplineMeshInfo->mExtractionCellToBsplineCell.size() );

        // get the number of active B-spline elements
        uint tNumBspElems = tBsplineMeshInfo->mExtractionCellsIndicesInBsplineCells.size();

        // sorting algorithm for every B-spline element
        for ( uint iBspElem = 0; iBspElem < tNumBspElems; iBspElem++ )
        {
            // initialize map of all SPs on BSp element
            IndexMap tSpToIpElem;

            // get number of IP cells in current B-spline element
            uint tNumIpCellsInBsplineElement = tBsplineMeshInfo->mExtractionCellsIndicesInBsplineCells( iBspElem ).size();

            // get the number of SPGs present on the current B-spline cell
            uint tNumSPGs = tBsplineMeshInfo->mSpgIndicesInBsplineCells( iBspElem ).size();

            // for every IP cell in current B-spline element ...
            for ( uint iIpCell = 0; iIpCell < tNumIpCellsInBsplineElement; iIpCell++ )
            {
                // get current IP cell's index
                moris_index tIpCellIndex = tBsplineMeshInfo->mExtractionCellsIndicesInBsplineCells( iBspElem )( iIpCell );

                // initialize the map storing the subphases in it
                tBsplineMeshInfo->mExtractionCellToSubPhase( tIpCellIndex ).resize( tNumSPGs );

                // get the SPs on the current element
                Cell< moris_index > const & tSubphaseIndicesOnLagElem = mCutIgMesh->get_parent_cell_subphases( tIpCellIndex );

                // populate the SP to IP-cell map
                for ( uint iSpOnIpCell = 0; iSpOnIpCell < tSubphaseIndicesOnLagElem.size(); iSpOnIpCell++ )
                {
                    // add SP to the map
                    tSpToIpElem[ tSubphaseIndicesOnLagElem( iSpOnIpCell ) ] = tIpCellIndex;
                }
            }

            // go through SPGs on B-spline cell
            for ( uint iSPG = 0; iSPG < tNumSPGs; iSPG++ )
            {
                // get the index of the SPG
                moris_index tSpgIndex = tBsplineMeshInfo->mSpgIndicesInBsplineCells( iBspElem )( iSPG );

                // get the list of SPs on the current SPG
                const Cell< moris_index >& tSPsInSPG = tBsplineMeshInfo->mSubphaseGroups( tSpgIndex )->get_SP_indices_in_group();

                // loop over SPs on current SPG and find which element they correspond to
                for ( uint iSpInSPG = 0; iSpInSPG < tSPsInSPG.size(); iSpInSPG++ )
                {
                    // get the SP's index
                    moris_index tSpIndex = tSPsInSPG( iSpInSPG );

                    // find the IP cell the SP lives in
                    auto tIter = tSpToIpElem.find( tSpIndex );

                    // check that the SP is actually found in one of the IP cells
                    MORIS_ASSERT( tIter != tSpToIpElem.end(),
                            "Integration_Mesh_Generator::establish_IP_SPG_SP_relationship() - SP not found in any IP cell within B-spline element" );

                    // get the index of the IP cell SP is found in
                    moris_index tIpCellIndex = tIter->second;

                    // add the IP-cell -- SPG -- SP relationship to the corresponding map
                    tBsplineMeshInfo->mExtractionCellToSubPhase( tIpCellIndex )( iSPG ).push_back( tSpIndex );

                    // TODO: is it better for efficiency to delete the map entry that is not needed anymore, to speed up subsequent find-calls?
                }
            }    // end: loop over all SPGs on B-spline element
        }        // end: loop over all B-spline elements
    }            // end: function

    // ----------------------------------------------------------------------------------

    void
    Enrichment::construct_enriched_interpolation_vertices_and_cells()
    {
        // log trace this function
        Tracer tTracer( "XTK", "Enrichment", "Construct Enriched Interpolation Vertices and Cells" );

        // get the enriched interpolation mesh pointer, this one is constructed here
        Enriched_Interpolation_Mesh* tEnrInterpMesh = mXTKModelPtr->mEnrichedInterpMesh( 0 );

        // geometry and interpolation order, limited to a single interpolation order mesh
        mtk::Cell const & tFirstCell = mBackgroundMeshPtr->get_mtk_cell( 0 );

        // set the interpolation mesh cell info (i.e. let the enriched IP mesh know what element type it uses)
        mtk::Cell_Info_Factory tFactory;
        tEnrInterpMesh->mCellInfo = tFactory.create_cell_info_sp( tFirstCell.get_geometry_type(), tFirstCell.get_interpolation_order() );

        // allocate indices and ids
        moris_index tCellIndex = 0;

        // maximum mesh index
        moris_index tMaxMeshIndex = mMeshIndices.max();

        // Enriched Interpolation Cell Index to Vertex Index map
        Matrix< IndexMat > tEnrInterpCellToVertex( tEnrInterpMesh->get_num_elements(), tEnrInterpMesh->mNumVertsPerInterpCell );

        // allocate vertex indices and ids
        // NOTE: THESE ARE NOT PARALLEL IDS
        moris_index tVertId      = 1;
        uint        tVertexCount = 0;

        // get the number of subphases in the current mesh
        uint tNumSPs = mCutIgMesh->get_num_subphases();

        // iterate through subphases and construct an interpolation cell in the interpolation mesh for each one
        for ( uint iSP = 0; iSP < tNumSPs; iSP++ )
        {
            // information about this cell
            moris::mtk::Cell* tParentCell = mCutIgMesh->get_subphase_parent_cell( iSP );

            // ID of the owning processor ("Owner")
            moris_id tOwner = tParentCell->get_owner();

            // vertices of cell
            Cell< mtk::Vertex* > tVertices = tParentCell->get_vertex_pointers();

            // bulk phase
            moris_index tBulkPhase = mCutIgMesh->get_subphase_bulk_phase( (moris_index)iSP );

            // loop over the B-spline meshes (i.e. discretization meshes)
            // this is done as the T-matrices at the enriched IP vertices need to be constructed wrt each B-spline mesh
            for ( uint iBspMesh = 0; iBspMesh < mMeshIndices.numel(); iBspMesh++ )
            {
                // Mesh Index
                moris_index tMeshIndex = mMeshIndices( iBspMesh );

                // get the T-matrix of the current IP node wrt the current B-spline mesh
                Cell< mtk::Vertex_Interpolation* > tVertexInterpolations = this->get_vertex_interpolations( *tParentCell, tMeshIndex );

                // get list of non-enriched basis indices interpolating into the current subphase
                Cell< moris_index > const & tBasisInCell = mEnrichmentData( tMeshIndex ).mSubphaseBGBasisIndices( (moris_index)iSP );

                // construct a map between non-enriched BF index and index relative to the subphase cluster
                Mini_Map< moris_id, moris_id > tCellBasisMap = construct_subphase_basis_to_basis_map( tBasisInCell );

                // get the enrichment levels of the BFs that are/need to be used for the interpolation in the current subphase
                Cell< moris_index > const & tEnrLevOfBasis = mEnrichmentData( tMeshIndex ).mSubphaseBGBasisEnrLev( (moris_index)iSP );

                // get the number of vertices per IP cell
                uint tNumVertices = tParentCell->get_number_of_vertices();

                // construct unzipped enriched vertices
                for ( uint iParentCellVertex = 0; iParentCellVertex < tNumVertices; iParentCellVertex++ )
                {
                    // construct vertex enrichment
                    Vertex_Enrichment tVertEnrichment;

                    // find out the enriched BF indices and IDs interpolating into the current vertex
                    // and store it in the Vertex_Enrichment object
                    this->construct_enriched_vertex_interpolation(
                            tMeshIndex,
                            tVertexInterpolations( iParentCellVertex ),
                            tEnrLevOfBasis,
                            tCellBasisMap,
                            tVertEnrichment );

                    // add vertex enrichment to enriched interpolation mesh
                    bool tNewVertFlag = false;

                    moris_index tVertEnrichIndex = tEnrInterpMesh->add_vertex_enrichment(
                            tMeshIndex,
                            tVertices( iParentCellVertex ),
                            tVertEnrichment,
                            tNewVertFlag );

                    // create this vertex on the first go around
                    // Note: the Interpolation_Vertex_Unzipped (UIPV) carries a list of vertex enrichments, each VE corresponds to one mesh index
                    // note though, that the UIPV is still created for every subphase (i.e. for every material sub-domain within the IP element)
                    if ( iBspMesh == 0 )
                    {
                        // Create interpolation vertex with only the first Vertex enrichment
                        tEnrInterpMesh->mEnrichedInterpVerts( tVertexCount ) =
                                new Interpolation_Vertex_Unzipped(
                                        tVertices( iParentCellVertex ),
                                        tVertId,
                                        tVertexCount,
                                        tVertices( iParentCellVertex )->get_owner(),
                                        tMeshIndex,
                                        tEnrInterpMesh->get_vertex_enrichment( tMeshIndex, tVertEnrichIndex ),
                                        tMaxMeshIndex );

                        // store enriched vertex's index for given parent cell index and element local node index
                        tEnrInterpCellToVertex( tCellIndex, iParentCellVertex ) = tVertexCount;

                        // update vertex ID to use for nex unzipped vertex
                        tVertId++;

                        // track number of unzipped vertices that have been created
                        tVertexCount++;
                    }
                    else
                    {
                        // the unzipped interpolation vertex' index
                        moris_index tVertexIndexInIp = tEnrInterpCellToVertex( tCellIndex - 1, iParentCellVertex );

                        // add the vertex interpolation for new mesh index
                        tEnrInterpMesh->mEnrichedInterpVerts( tVertexIndexInIp )    //
                                ->add_vertex_interpolation(
                                        tMeshIndex,
                                        tEnrInterpMesh->get_vertex_enrichment( tMeshIndex, tVertEnrichIndex ) );
                    }
                }

                // create the unzipped interpolation cell on first go
                /* Note: the Interpolation_Cell_Unzipped carries a list of Interpolation_Vertex_Unzipped (UIPV) which themselves get updated for every DMI
                 * Hence, the Interpolation_Cell_Unzipped can be left alone after initial creation.
                 * Access to the right UIPVs is given once they're all constructed (see code section with double for-loop just below) */
                if ( iBspMesh == 0 )
                {
                    // create new enriched interpolation cell and put it in list associating it with the underlying parent IP cell
                    tEnrInterpMesh->mEnrichedInterpCells( tCellIndex ) =
                            new Interpolation_Cell_Unzipped(
                                    tParentCell,
                                    iSP,
                                    tBulkPhase,
                                    mCutIgMesh->get_subphase_id( iSP ),
                                    tCellIndex,
                                    tOwner,
                                    tEnrInterpMesh->mCellInfo );

                    // add enriched interpolation cell to base cell to enriched cell data
                    tEnrInterpMesh->mBaseCellToEnrichedCell( tParentCell->get_index() ).push_back(    //
                            tEnrInterpMesh->mEnrichedInterpCells( tCellIndex ) );

                    // increment the cell index/id
                    tCellIndex++;
                }
            }    // end: loop over the Bspline meshes / discretization mesh indices
        }        // end: loop over subphases

        // resize out aura cells
        tEnrInterpCellToVertex.resize( tCellIndex, tEnrInterpMesh->mNumVertsPerInterpCell );
        tEnrInterpMesh->mEnrichedInterpCells.resize( tCellIndex );

        // with the cell to vertex data fully setup, add the vertex pointers to the cell
        // for every (unzipped) IP cell get its unzipped vertices
        for ( uint iIpCell = 0; iIpCell < tEnrInterpCellToVertex.n_rows(); iIpCell++ )
        {
            // initialize list of unzipped vertices on cell
            Cell< Interpolation_Vertex_Unzipped* > tVertices( tEnrInterpCellToVertex.n_cols() );

            // iterate through and get unzipped vertices on cell
            for ( uint iVertex = 0; iVertex < tEnrInterpCellToVertex.n_cols(); iVertex++ )
            {
                // store pointer to unzipped vertices in list
                tVertices( iVertex ) = tEnrInterpMesh->get_unzipped_vertex_pointer( tEnrInterpCellToVertex( iIpCell, iVertex ) );
            }

            // set vertices in cell
            tEnrInterpMesh->mEnrichedInterpCells( iIpCell )->set_vertices( tVertices );
        }

        // make sure list is only as big as it needs to be
        // FIXME: shouldn't this be a shrink-to-fit?
        tEnrInterpMesh->mEnrichedInterpVerts.resize( tVertexCount );
    }

    //-------------------------------------------------------------------------------------

    void
    Enrichment::construct_enriched_interpolation_vertices_and_cells_based_on_SPGs_new()
    {
        // log trace this function
        Tracer tTracer( "XTK", "Enrichment", "Construct Enriched Interpolation Vertices and Cells" );

        // get the enriched interpolation mesh pointer, this one is constructed here
        Enriched_Interpolation_Mesh* tEnrInterpMesh = mXTKModelPtr->mEnrichedInterpMesh( 0 );

        // geometry and interpolation order, limited to a single interpolation order mesh
        mtk::Cell const & tFirstCell = mBackgroundMeshPtr->get_mtk_cell( 0 );

        // set the interpolation mesh cell info (i.e. let the enriched IP mesh know what element type it uses)
        mtk::Cell_Info_Factory tFactory;
        tEnrInterpMesh->mCellInfo = tFactory.create_cell_info_sp( tFirstCell.get_geometry_type(), tFirstCell.get_interpolation_order() );

        // maximum mesh index
        moris_index tMaxMeshIndex = mMeshIndices.max();
        uint        tNumBspMeshes = mMeshIndices.numel();

        // Enriched Interpolation Cell Index to Vertex Index map
        Matrix< IndexMat > tEnrInterpCellToVertex( tEnrInterpMesh->get_num_elements(), tEnrInterpMesh->mNumVertsPerInterpCell );

        // allocate vertex indices and ids
        // NOTE: THESE ARE NOT PARALLEL IDS
        moris_index tVertId      = 1;
        uint        tVertexCount = 0;

        // loop over the B-spline meshes (i.e. discretization mesh indices)
        // this is done as the T-matrices at the enriched IP vertices need to be constructed wrt each B-spline mesh
        for ( uint iBspMesh = 0; iBspMesh < tNumBspMeshes; iBspMesh++ )
        {
            // Mesh Index
            moris_index tMeshIndex = mMeshIndices( iBspMesh );

            // iterate through all IP cells (iterator iIpCell is equal to IP cell's index)
            for ( uint iIpCell = 0; iIpCell < mBackgroundMeshPtr->get_num_elems(); iIpCell++ )
            {
                // convert to index
                moris_index tBaseIpCellIndex = (moris_index)iIpCell;

                // get the MTK cell
                mtk::Cell* tIpCell = &( mCutIgMesh->get_mtk_cell( tBaseIpCellIndex ) );

                // ID of the owning processor ("Owner")
                moris_id tOwner = tIpCell->get_owner();

                // get the number of vertices per IP cell
                uint tNumVertices = tIpCell->get_number_of_vertices();

                // vertices of IP cell
                Cell< mtk::Vertex* > tVertices = tIpCell->get_vertex_pointers();

                // get the non-enriched T-matrices of the current IP cell's nodes wrt the current B-spline mesh
                Cell< mtk::Vertex_Interpolation* > tVertexInterpolations = this->get_vertex_interpolations( *tIpCell, tMeshIndex );

                // sanity check that number of T-matrices matches number of nodes
                MORIS_ASSERT( tVertexInterpolations.size() == tNumVertices,
                        "Enrichment::construct_enriched_interpolation_vertices_and_cells_new() - "
                        "Number of T-matrices on cell and number of nodes don't match up." );

                // get max number of enr IP cells associated with the current IP cell
                uint tNumUnzippings = this->get_num_unzippings_of_base_ip_cell( iIpCell );

                // get the number of clusters with material in them
                Cell< moris_index > const & tMaterialSpgIndices  = mCutIgMesh->get_material_SPG_indices_for_base_IP_cell( iBspMesh, iIpCell );
                uint                        tNumMaterialClusters = tMaterialSpgIndices.size();

                // iterate through the enriched IP cells that are to be constructed on the current base IP cell
                for ( uint iUnzipping = 0; iUnzipping < tNumUnzippings; iUnzipping++ )
                {
                    // get the enr. IP cell's index
                    moris_index tUipcIndex = mEnrIpCellIndices( iIpCell )( iUnzipping );

                    // get the SPG index associated with the current UIPC
                    moris_index tAssociatedSpgIndex = mUnzippingToSpgIndex( iBspMesh )( iIpCell )( iUnzipping );

                    // store the UIPCs to SPG relationship
                    mUipcToSpgIndex( iBspMesh )( tUipcIndex ) = tAssociatedSpgIndex;

                    // check whether this UIPC needs a basis extension wrt. the current B-spline mesh
                    bool tNeedsBasisExtension = ( tAssociatedSpgIndex == -1 );

                    // get the SPG whose B-spline basis coefficients will be used
                    moris_index tBasisSpgIndex = -1;
                    if ( !tNeedsBasisExtension )
                    {
                        // for clusters with associated SPGs the assignment is trivial
                        tBasisSpgIndex = tAssociatedSpgIndex;
                    }
                    else    // extension needed
                    {
                        // if one of the material clusters requests an extension something must have gone wrong
                        bool tIsVoidCluster = ( iUnzipping >= tNumMaterialClusters );
                        MORIS_ERROR( tIsVoidCluster,
                                "Enrichment::construct_enriched_interpolation_vertices_and_cells_new() - "
                                "Unzipping #%i on base IP cell #%i (which should result in a material cluster) "
                                "marked for extension on B-spline mesh #%i. This shouldn't happen.",
                                iUnzipping,
                                iIpCell,
                                iBspMesh );

                        // see if a trivial extension is possible and if so, get the SPG from which the basis will be copied
                        tBasisSpgIndex = this->get_SPG_for_basis_extension( iIpCell, iUnzipping, iBspMesh );
                    }

                    // check whether the basis extension is trivial
                    bool tBasisExtensionIsNonTrivial = ( tBasisSpgIndex == -1 );

                    // initialize the list of active basis
                    Mini_Map< moris_id, moris_id > tBasisInCellMap;

                    // initialize the enriched T-matrix
                    Vertex_Enrichment tVertEnrichment;

                    // if the cluster doesn't need extension, or the extension is trivial,
                    // just collect the basis indices and compute the T-matrices for the individual vertices based on position later
                    if ( !tBasisExtensionIsNonTrivial )
                    {
                        // get list of non-enriched basis indices interpolating into the current subphase group
                        Cell< moris_index > tBasisInCell = mEnrichmentData( tMeshIndex ).mSubphaseGroupBGBasisIndices( tBasisSpgIndex );

                        // construct a map between non-enriched BF index and index relative to the list of basis interpolating into the current SPG
                        convert_cell_to_map( tBasisInCell, tBasisInCellMap );
                    }
                    else    // if a non-trivial basis extension is needed ...
                    {
                        // ... collect the basis indices and their weights in a larger B-spline element
                        // and construct a T-matrix that produces represents an average of the field
                        // this T-matrix will be copied into the enriched vertices where needed
                        this->construct_averaged_T_matrix_for_extension_on_enriched_IP_cell(
                                iIpCell,
                                iUnzipping,
                                iBspMesh,
                                tVertEnrichment );
                    }

                    // construction of unzipped enriched IP vertices
                    // i.e. only construct T-matrices if underlying enriched IP cell is valid
                    for ( uint iIpCellVertex = 0; iIpCellVertex < tNumVertices; iIpCellVertex++ )
                    {
                        // if the cluster doesn't need extension, or the extension is trivial ...
                        if ( !tBasisExtensionIsNonTrivial )
                        {
                            // ... create the enriched T-matrices in the usual manner
                            this->construct_enriched_vertex_interpolation(
                                    tMeshIndex,
                                    tVertexInterpolations( iIpCellVertex ),
                                    mEnrichmentData( tMeshIndex ).mSubphaseGroupBGBasisEnrLev( tBasisSpgIndex ),
                                    tBasisInCellMap,
                                    tVertEnrichment );
                        }
                        else    // else: the T-matrix to be copied is already filled with the average
                        {
                            // only set the related base vertex to the already constructed T-matrix
                            tVertEnrichment.add_base_vertex_interpolation( tVertexInterpolations( iIpCellVertex ) );
                        }

                        // add vertex enrichment to enriched interpolation mesh
                        bool         tNewVertFlag     = false;
                        mtk::Vertex* tBaseVertexPtr   = tVertices( iIpCellVertex );
                        moris_index  tVertEnrichIndex = tEnrInterpMesh->add_vertex_enrichment(
                                tMeshIndex,
                                tBaseVertexPtr,    // feed non-enriched Vertex
                                tVertEnrichment,
                                tNewVertFlag );

                        // create this vertex on the first go around
                        // NOTE: the Interpolation_Vertex_Unzipped (UIPV) carries a list of vertex enrichments, each VE corresponds to one mesh index
                        // note though, that the UIPV is still created for all subphase groups associated with every IP cell
                        if ( iBspMesh == 0 )
                        {
                            // Create interpolation vertex with only the first Vertex enrichment
                            tEnrInterpMesh->mEnrichedInterpVerts( tVertexCount ) =
                                    new Interpolation_Vertex_Unzipped(
                                            tBaseVertexPtr,                 // non-enriched IP vertices
                                            tVertId,                        // current IP vertex' ID
                                            tVertexCount,                   // index for enriched IP vertex
                                            tBaseVertexPtr->get_owner(),    // owning proc
                                            // NOTE: the two inputs below are only for initialization, not final information
                                            tMeshIndex,                                                               // current DMI
                                            tEnrInterpMesh->get_vertex_enrichment( tMeshIndex, tVertEnrichIndex ),    // T-matrix
                                            tMaxMeshIndex );                                                          // maximum DMI

                            // store enriched vertex's index for given parent cell index and element local node index
                            tEnrInterpCellToVertex( tUipcIndex, iIpCellVertex ) = tVertexCount;

                            // update vertex ID to use for next unzipped vertex
                            tVertId++;

                            // track number of unzipped vertices that have been created
                            tVertexCount++;
                        }
                        else
                        {
                            // the unzipped interpolation vertex' index
                            moris_index tUipvIndex = tEnrInterpCellToVertex( tUipcIndex, iIpCellVertex );

                            // get the vertex interpolation to be added
                            Vertex_Enrichment* tVertEnrichmentToAddToUipv = tEnrInterpMesh->get_vertex_enrichment( tMeshIndex, tVertEnrichIndex );

                            // add the vertex interpolation for new mesh index
                            tEnrInterpMesh->mEnrichedInterpVerts( tUipvIndex )    //
                                    ->add_vertex_interpolation( tMeshIndex, tVertEnrichmentToAddToUipv );
                        }
                    }    // end for: loop over IP vertices

                    // create the unzipped interpolation cell on first go
                    /* Note: the Interpolation_Cell_Unzipped carries a list of Interpolation_Vertex_Unzipped (UIPV) which themselves get updated for every DMI
                     * Hence, the Interpolation_Cell_Unzipped can be left alone after initial creation.
                     * Access to the right UIPVs is given once they're all constructed (see code section with double for-loop just below) */
                    if ( iBspMesh == 0 )
                    {
                        // get the bulk- and sub-phase indices for the primary phase
                        moris_index tPrimaryBulkPhase;
                        moris_index tPrimarySubPhase;
                        if ( iUnzipping < tNumMaterialClusters )    // material UIPC/cluster
                        {
                            // get the bulk-phase directly through the associated SPG
                            tPrimarySubPhase  = mCutIgMesh->get_parent_cell_subphases( tBaseIpCellIndex )( iUnzipping );
                            tPrimaryBulkPhase = mBsplineMeshInfos( iBspMesh )->get_bulk_phase_for_subphase_group( tAssociatedSpgIndex );
                        }
                        else    // UIPC/cluster is void, i.e. has no material
                        {
                            // set the subphase to -1 to mark this UIPC as void
                            tPrimarySubPhase = -1;

                            // get the bulk-phase of the sub-domain which the basis is extended from
                            tPrimaryBulkPhase = mCutIgMesh->get_bulk_phases_for_union_MSD_indices_for_base_IP_cell( tBaseIpCellIndex )( iUnzipping - tNumMaterialClusters );
                        }

                        // create new enriched interpolation cell and put it in list associating it with the underlying parent IP cell
                        tEnrInterpMesh->mEnrichedInterpCells( tUipcIndex ) =
                                new Interpolation_Cell_Unzipped(
                                        tIpCell,                      // Base IP cell
                                        tPrimarySubPhase,             // sub-phase index of the primary cells
                                        tPrimaryBulkPhase,            // Bulk-phase index of the primary cells
                                        tUipcIndex,                   // Index of Enr. IP cell
                                        tOwner,                       // Owning Proc of Enr. IP cell
                                        tNumBspMeshes,                // number of B-spline meshes being treated (for initializing internal lists)
                                        tEnrInterpMesh->mCellInfo,    // Cell info for the IP Cell type (e.g. Quad4, Hex27, etc.)
                                        true );                       // signal that the SPG base construction is used

                        // add enriched interpolation cell to base cell to enriched cell data
                        tEnrInterpMesh->mBaseCellToEnrichedCell( tIpCell->get_index() ).push_back(    //
                                tEnrInterpMesh->mEnrichedInterpCells( tUipcIndex ) );
                    }

                    // get the bulk-phase index, even for void clusters their hypothetical primary bulk-phase index
                    moris_index tBulkPhaseIndex = -1;
                    if ( !tNeedsBasisExtension )
                    {
                        tBulkPhaseIndex = mBsplineMeshInfos( iBspMesh )->get_bulk_phase_for_subphase_group( tAssociatedSpgIndex );
                    }

                    // store SPG and Bulk-phase indices for the current mesh index
                    tEnrInterpMesh->mEnrichedInterpCells( tUipcIndex )    //
                            ->set_SPG_and_BP_indices_for_DM_list_index( iBspMesh, tAssociatedSpgIndex, tBulkPhaseIndex );

                }    // end: loop over enr. IP cells on base IP cell
            }        // end: loop over IP cells on mesh
        }            // end: loop over DMIs

        // resize out aura cells
        tEnrInterpCellToVertex.resize( mNumEnrIpCells, tEnrInterpMesh->mNumVertsPerInterpCell );
        tEnrInterpMesh->mEnrichedInterpCells.resize( mNumEnrIpCells );

        // with the cell to vertex data fully setup, add the vertex pointers to the cell
        // for every (unzipped) IP cell get its unzipped vertices
        for ( uint iEnrIpCell = 0; iEnrIpCell < tEnrInterpCellToVertex.n_rows(); iEnrIpCell++ )
        {
            // initialize list of unzipped vertices on cell
            Cell< Interpolation_Vertex_Unzipped* > tVertices( tEnrInterpCellToVertex.n_cols() );

            // iterate through and get unzipped vertices on cell
            for ( uint iVertex = 0; iVertex < tEnrInterpCellToVertex.n_cols(); iVertex++ )
            {
                // store pointer to unzipped vertices in list
                tVertices( iVertex ) = tEnrInterpMesh->get_unzipped_vertex_pointer( tEnrInterpCellToVertex( iEnrIpCell, iVertex ) );
            }

            // set vertices in cell
            tEnrInterpMesh->mEnrichedInterpCells( iEnrIpCell )->set_vertices( tVertices );
        }

        // make sure list is only as big as it needs to be
        // FIXME: shouldn't this be a shrink-to-fit?
        tEnrInterpMesh->mEnrichedInterpVerts.resize( tVertexCount );
    }

    //-------------------------------------------------------------------------------------

    moris_index
    Enrichment::get_SPG_for_basis_extension(
            const moris_index aIpCellIndex,
            const moris_index aUnzippingOnIpCell,
            const moris_index aBsplineMeshIndex ) const
    {

        // get the number of clusters with material in them for the current IP cell
        uint tNumMaterialClusters = mCutIgMesh->get_material_SPG_indices_for_base_IP_cell( aBsplineMeshIndex, aIpCellIndex ).size();

        // get the unzipping index for the void UIPCs
        moris_index tVoidUnzipping = aUnzippingOnIpCell - (moris_index)tNumMaterialClusters;

        // make sure the input for the unzipping makes sense
        MORIS_ASSERT( tVoidUnzipping >= 0,
                "Enrichment::get_SPG_for_basis_extension() - "
                "UIPC provided is material cluster. Basis extensions don't exist on material clusters." );

        // get the MSD index associated with the current unzipped IP cell
        moris_index tMsdIndex = mCutIgMesh->get_union_MSD_indices_for_base_IP_cell( aIpCellIndex )( tVoidUnzipping );

        // get the list of void Msd Indices
        Cell< moris_index > const & tVoidMsdIndices = mCutIgMesh->get_void_MSD_indices_for_base_IP_cell( aBsplineMeshIndex, aIpCellIndex );

        // check if the the MSD index is found in that list
        moris_index tPositionFound = xtk::find_in_cell( tVoidMsdIndices, tMsdIndex );

        // if found ...
        if ( tPositionFound > -1 )
        {
            // ... return the corresponding SPG
            return mCutIgMesh->get_void_SPG_indices_for_base_IP_cell( aBsplineMeshIndex, aIpCellIndex )( tPositionFound );
        }
        else    // if not found among the void MSD indices...
        {
            // ... check if the MSD index is found in the list of material MSD indices
            Cell< moris_index > const & tMaterialMsdIndices    = mCutIgMesh->get_material_MSD_indices_for_base_IP_cell( aBsplineMeshIndex, aIpCellIndex );
            moris_index                 tPositionFoundMaterial = xtk::find_in_cell( tMaterialMsdIndices, tMsdIndex );

            // if it is found in the list of material MSD indices ...
            if ( tPositionFoundMaterial > -1 )
            {
                // ... return the corresponding SPG
                return mCutIgMesh->get_material_SPG_indices_for_base_IP_cell( aBsplineMeshIndex, aIpCellIndex )( tPositionFoundMaterial );
            }
        }

        // if the MSD index has not been found in the previous two lists, mark for non-trivial basis extension
        return -1;
    }

    //-------------------------------------------------------------------------------------

    void
    Enrichment::construct_averaged_T_matrix_for_extension_on_enriched_IP_cell(
            const moris_index  aIpCellIndex,
            const moris_index  aUnzippingIndex,
            const moris_index  aBsplineMeshIndex,
            Vertex_Enrichment& aAverageEnrichedTmatrix ) const
    {
        // check that the B-spline mesh index makes sense
        MORIS_ERROR( aBsplineMeshIndex != mCutIgMesh->get_coarsest_bspline_mesh_index_on_base_ip_cell( aIpCellIndex ),
                "Enrichment::collect_enriched_basis_indices_in_coarsest_Bspline_element() - "
                "A non-trivial basis extension should never be requested on the coarsest B-spline mesh, as a trivial one always exists." );

        // get the SPs in the material subdomain whose field needs to be extended
        Cell< moris_index > const & tSpsInMaterialSubDomain =
                this->collect_SPs_in_coarsest_element_in_same_material_subdomain_as_unzipped_cell( aIpCellIndex, aUnzippingIndex );
        uint tNumSpsInCoarsestElement = tSpsInMaterialSubDomain.size();

        // initialize lists storing cell-wise averaged T-matrices
        Cell< Vertex_Enrichment* > tAvgTmatrices( tNumSpsInCoarsestElement );
        Cell< real >               tRelativeIpCellVolumes( tNumSpsInCoarsestElement );

        // for each SP, and therefore material cluster, that exists in this material sub-domain on the coarsest B-spline element
        // collect an averaged T-matrix
        for ( uint iMaterialCluster = 0; iMaterialCluster < tNumSpsInCoarsestElement; iMaterialCluster++ )
        {
            // create a new T-matrix
            tAvgTmatrices( iMaterialCluster ) = new Vertex_Enrichment();

            // compute the averaged enriched T-matrix for this cluster
            this->get_averaged_vertex_enrichment_for_subphase(
                    tSpsInMaterialSubDomain( iMaterialCluster ),
                    aBsplineMeshIndex,
                    tAvgTmatrices( iMaterialCluster ) );

            // get the relative volume of the Lagrange element
            moris_index tSpIndex                       = tSpsInMaterialSubDomain( iMaterialCluster );
            mtk::Cell*  tBaseIpCell                    = mCutIgMesh->get_subphase_parent_cell( tSpIndex );
            real        tRefinementLevel               = (real)tBaseIpCell->get_base_cell()->get_level();
            real        tNumSpatialDims                = (real)mXTKModelPtr->get_spatial_dim();
            tRelativeIpCellVolumes( iMaterialCluster ) = std::pow( 0.5, tNumSpatialDims * tRefinementLevel );

        }    // end for: loop over material clusters for a given material sub-domain

        // create and compute average T-matrix
        aAverageEnrichedTmatrix = Vertex_Enrichment();
        this->average_T_matrices( tAvgTmatrices, tRelativeIpCellVolumes, aAverageEnrichedTmatrix );

        // delete the temporarily created T-matrices
        for ( auto iIter : tAvgTmatrices )
        {
            delete iIter;
        }
    }

    //-------------------------------------------------------------------------------------

    void
    Enrichment::fill_T_matrix_dummy( Vertex_Enrichment& aDummyEnrichedTmatrix ) const
    {
        // create dummy information
        Matrix< IndexMat > tBasisIndices( 1, 1, MORIS_INDEX_MAX - 1 );
        Matrix< IdMat >    tBasisIds( 1, 1, MORIS_ID_MAX - 1 );
        Matrix< IndexMat > tBasisOwners( 1, 1, 0 );
        Matrix< DDRMat >   tBasisWeights( 1, 1, std::numeric_limits< real >::signaling_NaN() );

        // create map storing relationship between proc-local basis index and the position in the map
        IndexMap& tBasisIndexToPositionMap              = aDummyEnrichedTmatrix.get_basis_map();
        tBasisIndexToPositionMap[ MORIS_INDEX_MAX - 1 ] = 0;

        // add the dummy information to the T-matrix
        aDummyEnrichedTmatrix.add_basis_information( tBasisIndices, tBasisIds );
        aDummyEnrichedTmatrix.add_basis_owners( tBasisIndices, tBasisOwners );
        aDummyEnrichedTmatrix.add_basis_weights( tBasisIndices, tBasisWeights );
        aDummyEnrichedTmatrix.add_base_vertex_interpolation( nullptr );
    }

    //-------------------------------------------------------------------------------------

    void
    Enrichment::average_T_matrices(
            Cell< Vertex_Enrichment* > const & aAverageTmatrices,
            Cell< real > const &               aWeights,
            Vertex_Enrichment&                 aAverageEnrichedTmatrix ) const
    {
        // check that the inputs make sense
        MORIS_ASSERT( aAverageTmatrices.size() == aWeights.size(),
                "Enrichment::average_T_matrices() - Number of T-matrices and weights provided do not match." );

        // initialize map finding position of basis indices in complete list of indices
        IndexMap& tBasisIndexToPositionMap = aAverageEnrichedTmatrix.get_basis_map();
        MORIS_ASSERT( tBasisIndexToPositionMap.size() == 0,
                "Enrichment::average_T_matrices() - Non-empty T-matrix passed for output." );

        // initialize counter indexing the positions of the basis functions
        uint tNumBfIndices = 0;

        // go over all basis indices and register them in the map
        for ( uint iTmatrix = 0; iTmatrix < aAverageTmatrices.size(); iTmatrix++ )
        {
            // get all the basis indices
            Matrix< IndexMat > const & tBasisIndices = aAverageTmatrices( iTmatrix )->get_basis_indices();

            // loop  over the basis indices in the current T-matrix and register them
            for ( uint iBF = 0; iBF < tBasisIndices.numel(); iBF++ )
            {
                // get the current enr. basis index
                moris_index tBfIndex = tBasisIndices( iBF );

                // find them in map
                auto tIter = tBasisIndexToPositionMap.find( tBfIndex );

                // if the index is not yet registered in the map add it and give it a position index
                if ( tIter == tBasisIndexToPositionMap.end() )
                {
                    tBasisIndexToPositionMap[ tBfIndex ] = tNumBfIndices;
                    tNumBfIndices++;
                }
            }
        }

        // initialize T-matrix data
        Matrix< IndexMat > tAvgBasisIndices( tNumBfIndices, 1, MORIS_INDEX_MAX );
        Matrix< IdMat >    tAvgBasisIds( tNumBfIndices, 1, MORIS_ID_MAX );
        Matrix< IndexMat > tAvgBasisOwners( tNumBfIndices, 1, MORIS_INDEX_MAX );
        Matrix< DDRMat >   tAvgBasisWeights( tNumBfIndices, 1, 0.0 );

        // loop over all T-matrices and add their information into the average T-matrix information
        for ( uint iTmatrix = 0; iTmatrix < aAverageTmatrices.size(); iTmatrix++ )
        {
            // get the information
            Matrix< IndexMat > const & tBasisIndices = aAverageTmatrices( iTmatrix )->get_basis_indices();
            Matrix< IdMat > const &    tBasisIds     = aAverageTmatrices( iTmatrix )->get_basis_ids();
            Matrix< IdMat >            tBasisOwners  = aAverageTmatrices( iTmatrix )->get_owners();
            Matrix< DDRMat > const &   tBasisWeights = aAverageTmatrices( iTmatrix )->get_basis_weights();

            // loop  over the basis indices in the current T-matrix and register them
            for ( uint iBF = 0; iBF < tBasisIndices.numel(); iBF++ )
            {
                // get the current enr. basis index and position in map
                moris_index tBfIndex = tBasisIndices( iBF );
                moris_index tPos     = tBasisIndexToPositionMap.find( tBfIndex )->second;

                // fill information
                tAvgBasisIndices( tPos ) = tBasisIndices( iBF );
                tAvgBasisIds( tPos )     = tBasisIds( iBF );
                tAvgBasisOwners( tPos )  = tBasisOwners( iBF );

                // sum weights
                tAvgBasisWeights( tPos ) += tBasisWeights( iBF );
            }
        }

        // get the total sum of weights
        real tSumOfWeights = xtk::sum_over_cell( aWeights );

        // average the weighted sum of weights
        tAvgBasisWeights = tAvgBasisWeights / tSumOfWeights;

        // copy the information into the
        aAverageEnrichedTmatrix.add_basis_information( tAvgBasisIndices, tAvgBasisIds );
        aAverageEnrichedTmatrix.add_basis_owners( tAvgBasisIndices, tAvgBasisOwners );
        aAverageEnrichedTmatrix.add_basis_weights( tAvgBasisIndices, tAvgBasisWeights );
        aAverageEnrichedTmatrix.add_base_vertex_interpolation( nullptr );

        // return: filled aAverageEnrichedTmatrix
    }

    //-------------------------------------------------------------------------------------

    Cell< moris_index > const &
    Enrichment::collect_SPs_in_coarsest_element_in_same_material_subdomain_as_unzipped_cell(
            const moris_index aIpCellIndex,
            const moris_index aUnzippingIndex ) const
    {
        // get the number of clusters with material in them for the current IP cell
        uint tNumMaterialClusters = mCutIgMesh->get_parent_cell_subphases( aIpCellIndex ).size();

        // get the unzipping index for the void UIPCs
        moris_index tVoidUnzipping = aUnzippingIndex - (moris_index)tNumMaterialClusters;

        // make sure the input for the unzipping makes sense
        MORIS_ASSERT( tVoidUnzipping >= 0,
                "Enrichment::collect_SPs_in_coarsest_element_in_same_material_subdomain_as_unzipped_cell() - "
                "UIPC provided is a material cluster. This function shouldn't be called here because "
                "basis extensions don't exist on material clusters." );

        // get the MSD index associated with the current unzipped IP cell
        moris_index tMsdIndex = mCutIgMesh->get_union_MSD_indices_for_base_IP_cell( aIpCellIndex )( tVoidUnzipping );

        // get access to the B-spline mesh information for the coarsest and current B-spline meshes
        moris_index        tCoarsestBspElemMeshIndex = mCutIgMesh->get_coarsest_bspline_mesh_index_on_base_ip_cell( aIpCellIndex );
        Bspline_Mesh_Info* tCoarsestBspMeshInfo      = mBsplineMeshInfos( tCoarsestBspElemMeshIndex );

        // get the SPs in the material subdomain whose field needs to be extended and return them
        moris_index tCoarseSpgIndex = tCoarsestBspMeshInfo->get_SPG_indices_associated_with_extraction_cell( aIpCellIndex )( tMsdIndex );
        return tCoarsestBspMeshInfo->mSubphaseGroups( tCoarseSpgIndex )->get_SP_indices_in_group();
    }

    //-------------------------------------------------------------------------------------

    void
    Enrichment::get_averaged_vertex_enrichment_for_subphase(
            const moris_index  aSubphaseIndex,
            const moris_index  aBsplineMeshIndex,
            Vertex_Enrichment* aAveragedTmatrix ) const
    {
        // get access to the B-spline mesh information for the B-spline mesh currently considered
        Bspline_Mesh_Info* tBspMeshInfo = mBsplineMeshInfos( aBsplineMeshIndex );

        // get the SPG index on the B-spline mesh that this SP is associated with
        moris_index tSpgIndex = tBspMeshInfo->mSpToSpgMap( aSubphaseIndex );

        // get the pointer to the base IP cell
        mtk::Cell* tIpCell = mCutIgMesh->get_subphase_parent_cell( aSubphaseIndex );

        // get access to the T-matrices of the non-enriched vertices
        Cell< mtk::Vertex_Interpolation* > tVertexInterpolations = this->get_vertex_interpolations( *tIpCell, aBsplineMeshIndex );
        uint                               tNumVertices          = tVertexInterpolations.size();

        // check that there are no empty aura vertices
        for ( uint iVert = 0; iVert < tVertexInterpolations.size(); iVert++ )
        {
            if ( tVertexInterpolations( iVert ) == nullptr )
            {
                // if there is one, return a nullptr for the averaged T-matrix
                aAveragedTmatrix = nullptr;
                return;
            }
        }

        // get the list of un-enriched basis indices and the enrichment levels for the current SPG
        Cell< moris_index > tBasisInSpg          = mEnrichmentData( aBsplineMeshIndex ).mSubphaseGroupBGBasisIndices( tSpgIndex );
        Cell< moris_index > tEnrLvlsOfBasisInSpg = mEnrichmentData( aBsplineMeshIndex ).mSubphaseGroupBGBasisEnrLev( tSpgIndex );

        // convert the list of basis into a searchable map
        IndexMap tBasisInSpgMap;
        xtk::convert_cell_to_map( tBasisInSpg, tBasisInSpgMap );

        // initialize map tracking where in the T-matrices the basis indices go
        IndexMap& tEnrBasisIndexToPositionMap = aAveragedTmatrix->get_basis_map();

        // initialize tracker for number of basis functions active on current enr. IP cell
        uint tNumActiveBFs = 0;

        // count the number of active enriched basis functions and establish map
        for ( uint iVert = 0; iVert < tNumVertices; iVert++ )
        {
            // get the non-enriched basis coefficients for the current vertex
            Matrix< IndexMat > tBaseVertCoeffInds = tVertexInterpolations( iVert )->get_indices();

            // treat each basis index from previous list
            for ( uint iBF = 0; iBF < tBaseVertCoeffInds.numel(); iBF++ )
            {
                // get the index of the current base BF
                moris_index tBfIndex = tBaseVertCoeffInds( iBF );

                // find it in list of basis indices in the current SPG
                auto tIter = tBasisInSpgMap.find( tBfIndex );
                MORIS_ASSERT( tIter != tBasisInSpgMap.end(),
                        "Enrichment::get_averaged_vertex_enrichment_for_subphase() - "
                        "Interpolation vertex contains contribution from basis that is not in the "
                        "list of basis indices in SPG #%i in the enrichment data for B-spline mesh #%i",
                        tSpgIndex,
                        aBsplineMeshIndex );

                // get the position of the basis index in the list of basis in the SPG
                moris_index tPosOfBfInSpg = tIter->second;

                // get the correct enrichment level
                moris_index tEnrLev = tEnrLvlsOfBasisInSpg( tPosOfBfInSpg );

                // get the enriched BF's index
                moris_index tEnrBfIndex = mEnrichmentData( aBsplineMeshIndex ).mBasisEnrichmentIndices( tBfIndex )( tEnrLev );

                // check if this enriched basis function has already been found, if not ...
                if ( tEnrBasisIndexToPositionMap.find( tEnrBfIndex ) == tEnrBasisIndexToPositionMap.end() )
                {
                    // ... store it and count up the number of active basis functions
                    tEnrBasisIndexToPositionMap[ tEnrBfIndex ] = tNumActiveBFs;
                    tNumActiveBFs++;
                }
            }
        }

        // initialize list of averaged weights
        Matrix< IndexMat > tEnrBasisIndices( tNumActiveBFs, 1, MORIS_INDEX_MAX );
        Matrix< IdMat >    tEnrBasisIds( tNumActiveBFs, 1, MORIS_ID_MAX );
        Matrix< IndexMat > tEnrBasisOwners( tNumActiveBFs, 1, MORIS_INDEX_MAX );
        Matrix< DDRMat >   tAveragedWeights( tNumActiveBFs, 1, 0.0 );

        // find enriched BF indices and IDs for every basis coefficient interpolating into the current subphase
        for ( uint iVert = 0; iVert < tNumVertices; iVert++ )
        {
            // collect information from non-enriched T-matrices
            uint                    tNumCoeffs         = tVertexInterpolations( iVert )->get_number_of_coefficients();
            Matrix< IndexMat >      tBaseVertCoeffInds = tVertexInterpolations( iVert )->get_indices();
            Matrix< IndexMat >      tBaseVertOwners    = tVertexInterpolations( iVert )->get_owners();
            const Matrix< DDRMat >* tBaseVertWeights   = tVertexInterpolations( iVert )->get_weights();

            // find enriched BF indices and IDs for every basis coefficient interpolating into the current element
            for ( uint iBF = 0; iBF < tNumCoeffs; iBF++ )
            {
                // get the current basis' index
                moris_index tBfIndex = tBaseVertCoeffInds( iBF );

                // find it in list of basis indices in the current SPG
                auto tIterPosInSpg = tBasisInSpgMap.find( tBfIndex );
                MORIS_ASSERT( tIterPosInSpg != tBasisInSpgMap.end(),
                        "Enrichment::get_averaged_vertex_enrichment_for_subphase() - "
                        "Interpolation vertex contains contribution from basis that is not in the "
                        "list of basis indices in SPG #%i in the enrichment data for B-spline mesh #%i",
                        tSpgIndex,
                        aBsplineMeshIndex );

                // get the position of the basis index in the list of basis in the SPG
                moris_index tPosOfBfInSpg = tIterPosInSpg->second;

                // get the correct enrichment level
                moris_index tEnrLev = tEnrLvlsOfBasisInSpg( tPosOfBfInSpg );

                // get the enriched BF's index
                moris_index tEnrBfIndex = mEnrichmentData( aBsplineMeshIndex ).mBasisEnrichmentIndices( tBfIndex )( tEnrLev );

                // get the position in the T-matrix
                auto tIterPosInTmat = tEnrBasisIndexToPositionMap.find( tEnrBfIndex );
                MORIS_ASSERT( tIterPosInTmat != tEnrBasisIndexToPositionMap.end(),
                        "Enrichment::get_averaged_vertex_enrichment_for_subphase() - "
                        "Map giving position to enriched basis indices in averaged T-matrix doesn't contain index. "
                        "Something must have gone wrong in the construction of that map." );
                moris_index tPosOfBfInTmat = tIterPosInTmat->second;

                // store the indices and IDs
                tEnrBasisIndices( tPosOfBfInTmat ) = tEnrBfIndex;
                tEnrBasisIds( tPosOfBfInTmat )     = mEnrichmentData( aBsplineMeshIndex ).mEnrichedBasisIndexToId( tEnrBfIndex );
                tEnrBasisOwners( tPosOfBfInTmat )  = tBaseVertOwners( iBF );

                // sum up the total basis weights
                tAveragedWeights( tPosOfBfInTmat ) += ( *tBaseVertWeights )( iBF );

            }    // end for: basis interpolating into vertex

        }        // end for: vertices on IP cell

        // average the weights
        tAveragedWeights = tAveragedWeights / (real)tNumVertices;

        // Fill the information of the enriched T-matrix
        aAveragedTmatrix->add_basis_information( tEnrBasisIndices, tEnrBasisIds );
        aAveragedTmatrix->add_basis_owners( tEnrBasisIndices, tEnrBasisOwners );
        aAveragedTmatrix->add_basis_weights( tEnrBasisIndices, tAveragedWeights );
        aAveragedTmatrix->add_base_vertex_interpolation( nullptr );
    }

    //-------------------------------------------------------------------------------------

    void
    Enrichment::communicate_unzipped_ip_cells()
    {
        // get a pointer to the enriched interpolation mesh which is constructed
        Enriched_Interpolation_Mesh* tEnrInterpMesh = mXTKModelPtr->mEnrichedInterpMesh( 0 );

        /* ---------------------------------------------------------------------------------------- */
        /* Step 0: Sort entities into owned and non-owned*/

        this->sort_unzipped_IP_cells_into_owned_and_not_owned();

        /* ---------------------------------------------------------------------------------------- */
        /* Step 0.5: get the parallel communication information */

        // get the communication table and equivalent map
        Matrix< IdMat > tCommTable = mXTKModelPtr->get_communication_table();

        /* ---------------------------------------------------------------------------------------- */
        /* Step 1: Let each proc decide how many UIPC-IDs it needs & communicate ID ranges */

        moris_id tMyFirstUipcId = (moris_id)get_processor_offset( mNumEnrIpCells );

        /* ---------------------------------------------------------------------------------------- */
        /* Step 2: Assign IDs to owned entities*/

        // iterate through UIPCs that the current proc owns and assign IDs to them
        for ( uint iCell = 0; iCell < tEnrInterpMesh->mOwnedEnrichedInterpCells.size(); iCell++ )
        {
            moris_index tUipcIndex = tEnrInterpMesh->mOwnedEnrichedInterpCells( iCell );
            tEnrInterpMesh->mEnrichedInterpCells( tUipcIndex )->set_id( tMyFirstUipcId );
            tMyFirstUipcId++;
        }

        /* ---------------------------------------------------------------------------------------- */
        /* The following steps are only necessary if code runs in parallel */

        if ( par_size() == 1 )    // serial
        {
            // check that all UIPCs are owned in serial
            MORIS_ASSERT( tEnrInterpMesh->mOwnedEnrichedInterpCells.size() == tEnrInterpMesh->mEnrichedInterpCells.size(),
                    "Enrichment::communicate_unzipped_ip_cells() - Code running in serial; not all UIPCs are owned by proc 0." );
        }
        else    // parallel
        {

            /* ---------------------------------------------------------------------------------------- */
            /* Step 3: Prepare requests for non-owned entities */

            // initialize lists of information that identifies UIPCs on other procs
            Cell< Cell< moris_index > > tUnzippedCellIndices;
            Cell< Matrix< IdMat > >     tBaseCellIds;
            Cell< Matrix< IndexMat > >  tUnzippingOnCells;
            Cell< Matrix< IndexMat > >  tBulkPhaseIndices;

            // fill the identifying information
            this->prepare_requests_for_not_owned_unzipped_IP_cell_IDs(
                    tUnzippedCellIndices,
                    tBaseCellIds,
                    tUnzippingOnCells,
                    tBulkPhaseIndices );

            /* ---------------------------------------------------------------------------------------- */
            /* Step 4: Send and Receive requests about non-owned entities to and from other procs */

            // initialize arrays for receiving
            Cell< Matrix< IdMat > >    tReceivedBaseCellIds;
            Cell< Matrix< IndexMat > > tReceivedUnzippingOnCells;
            Cell< Matrix< IndexMat > > tReceivedBulkPhaseIndices;

            // communicate information
            moris::communicate_mats( tCommTable, tBaseCellIds, tReceivedBaseCellIds );
            moris::communicate_mats( tCommTable, tUnzippingOnCells, tReceivedUnzippingOnCells );
            moris::communicate_mats( tCommTable, tBulkPhaseIndices, tReceivedBulkPhaseIndices );

            // clear unused memory
            tBaseCellIds.clear();
            tUnzippingOnCells.clear();
            tBulkPhaseIndices.clear();

            /* ---------------------------------------------------------------------------------------- */
            /* Step 5: Find answers to the requests */

            // initialize lists of ID answers to other procs
            Cell< Matrix< IdMat > > tAnswerUipcIds;

            this->prepare_answers_for_owned_unzipped_IP_cell_IDs(
                    tAnswerUipcIds,
                    tReceivedBaseCellIds,
                    tReceivedUnzippingOnCells,
                    tReceivedBulkPhaseIndices );

            // clear memory from now obsolete requests
            tReceivedBaseCellIds.clear();
            tReceivedUnzippingOnCells.clear();
            tReceivedBulkPhaseIndices.clear();

            /* ---------------------------------------------------------------------------------------- */
            /* Step 6: Send and receive answers to and from other procs */

            // initialize arrays for receiving
            Cell< Matrix< IdMat > > tReceivedAnswerBaseCellIds;

            // communicate answers
            communicate_mats( tCommTable, tAnswerUipcIds, tReceivedAnswerBaseCellIds );

            // clear unused memory
            tAnswerUipcIds.clear();

            /* ---------------------------------------------------------------------------------------- */
            /* Step 7: Use answers to assign IDs to non-owned UIPCs */

            this->handle_requested_unzipped_unzipped_IP_cell_answers( tUnzippedCellIndices, tReceivedAnswerBaseCellIds );

        }    // end if: parallel

    }        // end function: Enrichment::communicate_unzipped_ip_cells()

    //-------------------------------------------------------------------------------------

    void
    Enrichment::sort_unzipped_IP_cells_into_owned_and_not_owned()
    {
        // get a pointer to the enriched interpolation mesh which is constructed
        Enriched_Interpolation_Mesh* tEnrInterpMesh = mXTKModelPtr->mEnrichedInterpMesh( 0 );

        // allocate memory for lists of owned and non-owned UIPCs
        tEnrInterpMesh->mOwnedEnrichedInterpCells.reserve( mNumEnrIpCells );
        tEnrInterpMesh->mNotOwnedEnrichedInterpCells.reserve( mNumEnrIpCells );

        // go through all UIPCs and check whether they're owned or not
        for ( moris_index iUIPC = 0; iUIPC < mNumEnrIpCells; iUIPC++ )    // case: UIPC is owned by current proc
        {
            // get the underlying background element's owner
            moris_id tUipcOwner = tEnrInterpMesh->mEnrichedInterpCells( iUIPC )->get_base_cell()->get_owner();

            // sort into list according to owned and non-owned
            if ( tUipcOwner == par_rank() )    // case: owned UIPC
            {
                // add index to list of owned UIPCs
                tEnrInterpMesh->mOwnedEnrichedInterpCells.push_back( iUIPC );
            }
            else    // case: non-owned UIPC
            {
                // add index to list of non-owned UIPCs (do not assign an ID to these)
                tEnrInterpMesh->mNotOwnedEnrichedInterpCells.push_back( iUIPC );
            }
        }

        // free unused memory
        tEnrInterpMesh->mOwnedEnrichedInterpCells.shrink_to_fit();
        tEnrInterpMesh->mNotOwnedEnrichedInterpCells.shrink_to_fit();

    }    // end function: Enrichment::sort_unzipped_IP_cells_into_owned_and_not_owned()

    //-------------------------------------------------------------------------------------

    void
    Enrichment::prepare_requests_for_not_owned_unzipped_IP_cell_IDs(
            Cell< Cell< moris_index > >& aUnzippedCellIndices,
            Cell< Matrix< IdMat > >&     aBaseCellIds,
            Cell< Matrix< IndexMat > >&  aUnzippingOnCells,
            Cell< Matrix< IndexMat > >&  aBulkPhaseIndices )
    {
        // get a pointer to the enriched interpolation mesh which is constructed
        Enriched_Interpolation_Mesh* tEnrInterpMesh = mXTKModelPtr->mEnrichedInterpMesh( 0 );

        // get the communication table and equivalent map
        Matrix< IdMat >                   tCommTable              = mXTKModelPtr->get_communication_table();
        uint                              tCommTableSize          = tCommTable.numel();
        std::map< moris_id, moris_index > tProcIdToCommTableIndex = mXTKModelPtr->get_communication_map();

        // get the number of non-owned entities on the executing processor
        uint tNumNotOwnedUIPCs = tEnrInterpMesh->mNotOwnedEnrichedInterpCells.size();

        // initialize request arrays with correct size
        aUnzippedCellIndices.resize( tCommTableSize );
        aUnzippedCellIndices.reserve( tNumNotOwnedUIPCs );
        aBaseCellIds.resize( tCommTableSize );
        aUnzippingOnCells.resize( tCommTableSize );
        aBulkPhaseIndices.resize( tCommTableSize );

        // prepare list that give the position of the requested UIPV in the array of non-owned UIPVs
        Cell< Cell< moris_index > > tUipcPositionInNotOwnedList( tCommTableSize );

        // go through SPGs that executing proc knows about, but doesn't own, ...
        for ( uint iNotOwnedUIPC = 0; iNotOwnedUIPC < tNumNotOwnedUIPCs; iNotOwnedUIPC++ )
        {
            // ... get their index ...
            moris_index tUipcIndex = tEnrInterpMesh->mNotOwnedEnrichedInterpCells( iNotOwnedUIPC );

            // ... get their respective owners, and position in the comm table ...
            moris_index tOwnerProc = tEnrInterpMesh->mEnrichedInterpCells( tUipcIndex )->get_owner();
            auto        tIter      = tProcIdToCommTableIndex.find( tOwnerProc );
            MORIS_ASSERT(
                    tIter != tProcIdToCommTableIndex.end(),
                    "Enriched_Interpolation_Mesh::prepare_requests_for_not_owned_unzipped_vertex_IDs() - "
                    "Entity owner (Proc #%i) not found in communication table of current proc #%i which is: %s",
                    tOwnerProc,
                    par_rank(),
                    ios::stringify_log( tCommTable ).c_str() );
            moris_index tProcDataIndex = tIter->second;

            // ... and finally add the non-owned SPGs in the list of SPs to be requested from that owning proc
            aUnzippedCellIndices( tProcDataIndex ).push_back( tUipcIndex );

            // store the where the current UIPV can be found in the list of not owned UIPVs
            tUipcPositionInNotOwnedList( tProcDataIndex ).push_back( iNotOwnedUIPC );
        }

        // size out unused memory
        aUnzippedCellIndices.shrink_to_fit();
        tUipcPositionInNotOwnedList.shrink_to_fit();

        // assemble identifying information for every processor communicated with
        for ( uint iProc = 0; iProc < tCommTableSize; iProc++ )
        {
            // get the number of non-owned entities to be sent to each processor processor
            uint tNumNotOwnedEntitiesOnProc = aUnzippedCellIndices( iProc ).size();

            // allocate matrix
            aBaseCellIds( iProc ).resize( tNumNotOwnedEntitiesOnProc, 1 );
            aUnzippingOnCells( iProc ).resize( tNumNotOwnedEntitiesOnProc, 1 );
            aBulkPhaseIndices( iProc ).resize( tNumNotOwnedEntitiesOnProc, 1 );

            // go through the Subphase groups for which IDs will be requested by the other processor
            for ( uint iCell = 0; iCell < tNumNotOwnedEntitiesOnProc; iCell++ )
            {
                // get the index of the UIPC on the executing proc and its ID
                moris_index tUipcIndex = aUnzippedCellIndices( iProc )( iCell );

                // get the UIPC
                Interpolation_Cell_Unzipped* tUIPC = tEnrInterpMesh->mEnrichedInterpCells( tUipcIndex );

                // get the base cell ID
                moris_id tBaseCellId = tUIPC->get_base_cell()->get_id();

                // get the Unzipping index of the current UIPC
                moris_index tUnzipping = mUipcUnzippingIndices( tUipcIndex );

                // check that the unzipping makes sense
                MORIS_ASSERT( tUnzipping < MORIS_INDEX_MAX,
                        "Enrichment::assign_IDs_to_owned_unzipped_IP_cells() - "
                        "UIPC to Unzipping map returned MORIS_INDEX_MAX, incomplete map construction" );

                // get the bulk-phase index
                moris_index tBulkPhase = tUIPC->get_bulkphase_index();

                // store the identifying information in the output arrays
                aBaseCellIds( iProc )( iCell )      = tBaseCellId;
                aUnzippingOnCells( iProc )( iCell ) = tUnzipping;
                aBulkPhaseIndices( iProc )( iCell ) = tBulkPhase;
            }

        }    // end for: each proc communicated with

    }        // end function: Enrichment::assign_IDs_to_owned_unzipped_IP_cells()

    //-------------------------------------------------------------------------------------

    void
    Enrichment::prepare_answers_for_owned_unzipped_IP_cell_IDs(
            Cell< Matrix< IdMat > >&           aAnswerUipcIds,
            Cell< Matrix< IdMat > > const &    aReceivedBaseCellIds,
            Cell< Matrix< IndexMat > > const & aReceivedUnzippingOnCells,
            Cell< Matrix< IndexMat > > const & aReceivedBulkPhaseIndices )
    {
        // get a pointer to the enriched interpolation mesh which is constructed
        Enriched_Interpolation_Mesh* tEnrInterpMesh = mXTKModelPtr->mEnrichedInterpMesh( 0 );

        // get the size of the communication table
        moris_index tCommTableSize = mXTKModelPtr->get_communication_table().numel();

        // initialize answer array with correct size
        aAnswerUipcIds.resize( tCommTableSize );

        // answer requests from each proc
        for ( moris_id iProcInCommTable = 0; iProcInCommTable < tCommTableSize; iProcInCommTable++ )
        {
            // check for how many UIPCs IDs have been requested by the current Proc
            uint tNumRequestedUipcIds = aReceivedBaseCellIds( iProcInCommTable ).numel();

            // initialize lists of ID answers to other procs
            aAnswerUipcIds( iProcInCommTable ).set_size( tNumRequestedUipcIds, 1 );

            // go through all UIPCs for which IDs have been requested
            for ( uint iReqUIPC = 0; iReqUIPC < tNumRequestedUipcIds; iReqUIPC++ )
            {
                // make temporary copy of the received information
                moris_id    tBaseCellId = aReceivedBaseCellIds( iProcInCommTable )( iReqUIPC );
                moris_index tUnzipping  = aReceivedUnzippingOnCells( iProcInCommTable )( iReqUIPC );

                // get the index of the base cell
                moris_index tBaseCellIndex =
                        mXTKModelPtr->mBackgroundMesh->get_loc_entity_ind_from_entity_glb_id( tBaseCellId, mtk::EntityRank::ELEMENT );

                // get the index of the requested UIPC on the current proc
                moris_index tUipcIndex = mEnrIpCellIndices( tBaseCellIndex )( tUnzipping );

                // get the pointer to the UIPC
                Interpolation_Cell_Unzipped* tUIPC = tEnrInterpMesh->mEnrichedInterpCells( tUipcIndex );

                // get the ID of the UIPC
                moris_id tUipcID = tUIPC->get_id();

                // check if the output makes sense
                MORIS_ASSERT( tUipcID != MORIS_ID_MAX,
                        "Enrichment::communicate_unzipped_ip_cells() - UIPC request from other proc cannot be answered by current proc. "
                        "This UIPC has ID = MORIS_ID_MAX on current proc, i.e. it is not owned by current proc" );

                MORIS_ASSERT( tUIPC->get_bulkphase_index() == aReceivedBulkPhaseIndices( iProcInCommTable )( iReqUIPC ),
                        "Enrichment::communicate_unzipped_ip_cells() - "
                        "UIPC has different bulk phases on requesting and receiving proc." );

                // write the id into the answer array
                aAnswerUipcIds( iProcInCommTable )( iReqUIPC ) = tUipcID;
            }

        }    // end for: each proc communicated with

    }        // end function: Enrichment::prepare_answers_for_owned_unzipped_IP_cell_IDs()

    //-------------------------------------------------------------------------------------

    void
    Enrichment::handle_requested_unzipped_unzipped_IP_cell_answers(
            Cell< Cell< moris_index > > const & aUnzippedCellIndices,
            Cell< Matrix< IdMat > > const &     aReceivedAnswerBaseCellIds )
    {
        // get a pointer to the enriched interpolation mesh which is constructed
        Enriched_Interpolation_Mesh* tEnrInterpMesh = mXTKModelPtr->mEnrichedInterpMesh( 0 );

        // answers received from each proc
        for ( uint iProc = 0; iProc < aUnzippedCellIndices.size(); iProc++ )
        {
            // check for how many UIPCs IDs have been requested by the current Proc
            uint tNumAnsweredUipcIds = aReceivedAnswerBaseCellIds( iProc ).numel();

            // go through all UIPCs for which IDs have been requested
            for ( uint iAnsUIPC = 0; iAnsUIPC < tNumAnsweredUipcIds; iAnsUIPC++ )
            {
                // get the ID of the UIPC
                moris_id tUipcId = aReceivedAnswerBaseCellIds( iProc )( iAnsUIPC );

                // get the index of the UIPC
                moris_index tUipcIndex = aUnzippedCellIndices( iProc )( iAnsUIPC );

                // check that UIPC doesn't have an index attached to it already
                MORIS_ASSERT( tEnrInterpMesh->mEnrichedInterpCells( tUipcIndex )->get_id() == MORIS_ID_MAX,
                        "Enrichment::communicate_unzipped_ip_cells() - "
                        "Trying to assign parallel communicated UIPC ID to an UIPC that already has an ID" );

                // assign ID to UIPC
                tEnrInterpMesh->mEnrichedInterpCells( tUipcIndex )->set_id( tUipcId );
            }

        }    // end for: each proc communicated with

    }        // end function: Enrichment::handle_requested_unzipped_unzipped_IP_cell_answers()

    //-------------------------------------------------------------------------------------

    void
    Enrichment::construct_UIPC_to_unzipping_index()
    {
        // allocate the size of the map
        mUipcUnzippingIndices.resize( mNumEnrIpCells, MORIS_INDEX_MAX );

        // go through all base IP cells and their unzippings to construct the reverse map
        for ( uint iBaseIpCell = 0; iBaseIpCell < mEnrIpCellIndices.size(); iBaseIpCell++ )
        {
            for ( uint iUnzipping = 0; iUnzipping < mEnrIpCellIndices( iBaseIpCell ).size(); iUnzipping++ )
            {
                // get the UIPC index
                moris_index tUipcIndex = mEnrIpCellIndices( iBaseIpCell )( iUnzipping );

                // assign the UIPCs unzipping index to it
                mUipcUnzippingIndices( tUipcIndex ) = (moris_index)iUnzipping;
            }
        }
    }

    //-------------------------------------------------------------------------------------

    void
    Enrichment::construct_enriched_vertex_interpolation(
            moris_index const &             aEnrichmentDataIndex,
            mtk::Vertex_Interpolation*      aBaseVertexInterp,
            Cell< moris_index > const &     aSubPhaseBasisEnrLev,
            Mini_Map< moris_id, moris_id >& aMapBasisIndexToLocInSubPhase,
            Vertex_Enrichment&              aVertexEnrichment )
    {
        // allocate a new vertex enrichment
        aVertexEnrichment = Vertex_Enrichment();

        // a nullptr here would indicate an aura node without a t-matrix
        if ( aBaseVertexInterp != nullptr )
        {
            // iterate through basis in the base vertex interpolation
            uint tNumCoeffs = aBaseVertexInterp->get_number_of_coefficients();

            // indices of the coefficients
            Matrix< IndexMat > tBaseVertCoeffInds = aBaseVertexInterp->get_indices();

            // get owners
            Matrix< IndexMat > tBaseVertOwners = aBaseVertexInterp->get_owners();

            // weights of the coefficients
            const Matrix< DDRMat >* tBaseVertWeights = aBaseVertexInterp->get_weights();

            // enriched Basis Coefficient indices and IDs
            Matrix< IndexMat > tEnrichCoeffInds( tBaseVertCoeffInds.n_rows(), tBaseVertCoeffInds.n_cols() );
            Matrix< IndexMat > tEnrichCoeffIds( tBaseVertCoeffInds.n_rows(), tBaseVertCoeffInds.n_cols() );

            // find enriched BF indices and IDs for every basis coefficient interpolating into the current subphase
            for ( uint iBC = 0; iBC < tNumCoeffs; iBC++ )
            {
                // coefficient of the base index
                moris_index tCoeffIndex = tBaseVertCoeffInds( iBC );

                // find the coefficients index within this subphase cluster
                auto tIter = aMapBasisIndexToLocInSubPhase.find( tCoeffIndex );

                MORIS_ASSERT( tIter != aMapBasisIndexToLocInSubPhase.end(),
                        "Enrichment::construct_enriched_vertex_interpolation() - Basis not found in vertex map" );

                // The basis local index relative to the subphase
                moris_index tSubphaseBasisIndex = tIter->second;

                // enrichment level of the basis
                moris_index tEnrLev = aSubPhaseBasisEnrLev( tSubphaseBasisIndex );

                // get the enriched BF's index
                moris_index tEnrichedCoeffIndex = mEnrichmentData( aEnrichmentDataIndex ).mBasisEnrichmentIndices( tCoeffIndex )( tEnrLev );

                // store the enriched BF's index associated with current basis coefficient
                tEnrichCoeffInds( iBC ) = tEnrichedCoeffIndex;

                // store the enriched BF's ID associated with current basis coefficient
                tEnrichCoeffIds( iBC ) = mEnrichmentData( aEnrichmentDataIndex ).mEnrichedBasisIndexToId( tEnrichedCoeffIndex );
            }

            // get access to the basis to local index map of the vertex enrichment for modification
            Mini_Map< moris_index, moris_index >& tVertEnrichMap = aVertexEnrichment.get_basis_map();

            // store in the vertex enrichment the list of enriched BF coefficients associated with it (in the form of an index map)
            for ( uint iBC = 0; iBC < tEnrichCoeffInds.numel(); iBC++ )
            {
                moris_index tBasisIndex       = tEnrichCoeffInds( iBC );
                tVertEnrichMap[ tBasisIndex ] = iBC;
            }

            // feed the generated data about the enriched vertex interpolation to the vertex enrichment object
            aVertexEnrichment.add_basis_information( tEnrichCoeffInds, tEnrichCoeffIds );
            aVertexEnrichment.add_basis_owners( tEnrichCoeffInds, tBaseVertOwners );
            aVertexEnrichment.add_basis_weights( tEnrichCoeffInds, *tBaseVertWeights );
            aVertexEnrichment.add_base_vertex_interpolation( aBaseVertexInterp );
        }
    }

    //-------------------------------------------------------------------------------------

    Mini_Map< moris_id, moris_id >
    Enrichment::construct_subphase_basis_to_basis_map( Cell< moris_id > const & aSubPhaseBasisIndex )
    {
        // get number of BFs interpolating into subphase
        uint tNumBasisOfSubphase = aSubPhaseBasisIndex.size();

        // initialize output map
        Mini_Map< moris_id, moris_id > tSubphaseBasisMap;

        // for each basis interpolating into the current Subphase, ...
        for ( uint iB = 0; iB < tNumBasisOfSubphase; iB++ )
        {
            // ... get its Index and relate it to the position in the list of basis for that particular subphase
            tSubphaseBasisMap[ aSubPhaseBasisIndex( iB ) ] = iB;
        }

        // return the map
        return tSubphaseBasisMap;
    }

    //-------------------------------------------------------------------------------------

    Cell< mtk::Vertex_Interpolation* >
    Enrichment::get_vertex_interpolations(
            moris::mtk::Cell& aParentCell,
            const uint        aMeshIndex ) const
    {
        uint tNumVerts = aParentCell.get_number_of_vertices();

        Cell< mtk::Vertex* > tVertexPointers = aParentCell.get_vertex_pointers();

        Cell< mtk::Vertex_Interpolation* > tVertexInterp( tNumVerts );

        for ( uint i = 0; i < tNumVerts; i++ )
        {
            moris_index tVertexIndex = tVertexPointers( i )->get_index();

            tVertexInterp( i ) = mEnrichmentData( aMeshIndex ).mBGVertexInterpolations( tVertexIndex );
        }

        return tVertexInterp;
    }

    //-------------------------------------------------------------------------------------

    Cell< std::string >
    Enrichment::get_cell_enrichment_field_names() const
    {
        // number of basis types
        uint tNumBasisTypes = mEnrichmentData.size();

        // declare  enrichment fields
        Cell< std::string > tEnrichmentFields;

        for ( uint iBT = 0; iBT < tNumBasisTypes; iBT++ )
        {
            // number of basis
            moris::size_t tNumBasis = mBackgroundMeshPtr->get_num_basis_functions( mMeshIndices( iBT ) );

            std::string tBaseEnrich = "el_bt_" + std::to_string( mMeshIndices( iBT ) ) + "b_";
            for ( size_t i = 0; i < tNumBasis; i++ )
            {
                tEnrichmentFields.push_back(
                        tBaseEnrich + std::to_string( mBackgroundMeshPtr->get_glb_entity_id_from_entity_loc_index( i, mBasisRank, mMeshIndices( iBT ) ) ) );
            }
        }

        // Add local flood-fill field to the output mesh
        std::string tLocalFFStr = "child_ff";
        tEnrichmentFields.push_back( tLocalFFStr );

        std::string tSubPhaseStr = "subphase";
        tEnrichmentFields.push_back( tSubPhaseStr );

        return tEnrichmentFields;
    }

    //-------------------------------------------------------------------------------------

    void
    Enrichment::write_cell_enrichment_to_fields(
            Cell< std::string >& aEnrichmentFieldStrs,
            mtk::Mesh*           aMeshWithEnrFields ) const
    {
        MORIS_ERROR( 0, "Deprecated. (Removal in progress)" );
    }

    //-------------------------------------------------------------------------------------

    moris_index
    Enrichment::allocate_basis_ids(
            moris_index const & aMeshIndex,
            moris_index const & aNumIdsToAllocate )
    {
        int tProcRank = par_rank();
        int tProcSize = par_size();

        // size_t is defined as uint here because of aNumRequested
        // Initialize gathered information outputs (information which will be scattered across processors)
        Cell< moris_id > aGatheredInfo;
        Cell< moris_id > tFirstId( 1 );
        Cell< moris_id > tNumIdsRequested( 1 );

        tNumIdsRequested( 0 ) = (moris_id)aNumIdsToAllocate;

        moris::gather( tNumIdsRequested, aGatheredInfo );

        Cell< moris_id > tProcFirstID( tProcSize );

        // local maximum
        moris_index tLocMaxId = this->get_max_basis_id( aMeshIndex );

        // global maximum
        moris_index tFirstAvailGlbMaxId = moris::max_all( tLocMaxId ) + 1;
        if ( tProcRank == 0 )
        {
            // Loop over entities print the number of entities requested by each processor
            for ( int iProc = 0; iProc < tProcSize; ++iProc )
            {
                // Give each processor their desired amount of IDs
                tProcFirstID( iProc ) = tFirstAvailGlbMaxId;

                // Increment the first available node ID
                tFirstAvailGlbMaxId = tFirstAvailGlbMaxId + aGatheredInfo( iProc );
            }
        }

        moris::scatter( tProcFirstID, tFirstId );

        return tFirstId( 0 );
    }

    //-------------------------------------------------------------------------------------

    moris_index
    Enrichment::get_max_basis_id( moris_index const & aMeshIndex )
    {
        // Number of basis functions
        moris::size_t tNumBasis = mXTKModelPtr->get_background_mesh().get_num_basis_functions( aMeshIndex );

        // maximum id
        moris_id tMaxId = 0;

        for ( uint i = 0; i < tNumBasis; i++ )
        {
            // get the basis id
            moris_id tBasisId = mBackgroundMeshPtr->get_glb_entity_id_from_entity_loc_index( i, mBasisRank, aMeshIndex );

            if ( tMaxId < tBasisId )
            {
                tMaxId = tBasisId;
            }
        }
        return tMaxId;
    }

    //-------------------------------------------------------------------------------------

    void
    Enrichment::construct_enriched_basis_in_subphase_group_map()
    {
        // loop over DMI
        for ( const auto& iMeshIndex : mMeshIndices )
        {
            // get a reference of the map and its transpose
            Cell< Matrix< IndexMat > >&  tSubphaseGroupIndsInEnrichedBasis = mEnrichmentData( iMeshIndex ).mSubphaseGroupIndsInEnrichedBasis;
            Cell< Cell< moris_index > >& tEnrichedBasisInSubphaseGroup     = mEnrichmentData( iMeshIndex ).mEnrichedBasisInSubphaseGroup;

            // Get the number of subphase groups (on the current proc) and resize
            moris_id tNumSubphaseGroups = mBsplineMeshInfos( iMeshIndex )->get_num_SPGs();
            tEnrichedBasisInSubphaseGroup.resize( tNumSubphaseGroups );

            // estimate number of inner cells and use the estimate to resize
            uint tOrder = mXTKModelPtr->mBackgroundMesh->get_discretization_order( iMeshIndex );
            uint tDim   = mXTKModelPtr->mBackgroundMesh->get_spatial_dim();

            // resize the inner cells based on the order estimate
            for ( auto& iEnrichedBasisInSubphaseGroup : tEnrichedBasisInSubphaseGroup )
            {
                iEnrichedBasisInSubphaseGroup.reserve( tOrder + tDim );
            }

            // loop over the BF indices to populate the map
            for ( size_t iEnrichedBFIndex = 0; iEnrichedBFIndex < tSubphaseGroupIndsInEnrichedBasis.size(); iEnrichedBFIndex++ )
            {
                // get the subphases on which BF is active
                Matrix< IndexMat >& iSPGsInSupport = tSubphaseGroupIndsInEnrichedBasis( iEnrichedBFIndex );

                // loop over the SPG indices and put the enriched BF index in the corresponding SPG index
                std::for_each( iSPGsInSupport.begin(), iSPGsInSupport.end(), [ &tEnrichedBasisInSubphaseGroup, &iEnrichedBFIndex ]    //
                        ( moris_index const & aSPGIndex ) {
                            // construct the inverse map
                            tEnrichedBasisInSubphaseGroup( aSPGIndex ).push_back( iEnrichedBFIndex );
                        } );
            }
            // trim the extra spaces generated
            shrink_to_fit_all( tEnrichedBasisInSubphaseGroup );
        }
    }

    //-------------------------------------------------------------------------------------

    moris::Cell< Enrichment_Data > const &
    Enrichment::get_enrichment_data() const
    {
        return mEnrichmentData;
    }

}    // namespace xtk
