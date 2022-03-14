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

#include "fn_sort_points_by_coordinates.hpp"

#include "cl_Tracer.hpp"

namespace xtk
{
    //-------------------------------------------------------------------------------------

    Enrichment::Enrichment(
            enum Enrichment_Method const & aMethod,
            enum EntityRank const &        aBasisRank,
            Matrix< IndexMat > const &     aInterpIndex,
            moris::moris_index const &     aNumBulkPhases,
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
        mCutIgMesh   = mXTKModelPtr->get_cut_integration_mesh();
        mIgMeshTools = new Integration_Mesh_Generator();
        mBsplineMeshInfos = mCutIgMesh->get_bspline_mesh_info();

        // FIXME: this needs to go once the SPG based enrichment is validated
        // reinitialize the enrichment data for SPGs if needed
        if( aUseSpgBasedEnrichment )
        {
            // initialize the enrichment data for every discretization mesh index, as each mesh will result in a different number of SPGs
            for ( uint iMeshIndex = 0; iMeshIndex < aInterpIndex.numel(); iMeshIndex++)
            {
                moris_index tMeshIndex = aInterpIndex( iMeshIndex );
                uint tNumSPGs = aXTKModelPtr->get_cut_integration_mesh()->get_num_subphase_groups( iMeshIndex );
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

    Cell< moris::Matrix< moris::IdMat > > const &
    Enrichment::get_element_inds_in_basis_support( moris_index const & aEnrichmentDataIndex ) const
    {
        return mEnrichmentData( aEnrichmentDataIndex ).mElementIndsInBasis;
    }

    //-------------------------------------------------------------------------------------

    Cell< moris::Matrix< moris::IndexMat > > const &
    Enrichment::get_element_enrichment_levels_in_basis_support( moris_index const & aEnrichmentDataIndex ) const
    {
        return mEnrichmentData( aEnrichmentDataIndex ).mElementEnrichmentLevel;
    }

    //-------------------------------------------------------------------------------------

    Cell< moris::Matrix< moris::IndexMat > > const &
    Enrichment::get_subphases_loc_inds_in_enriched_basis( moris_index const & aEnrichmentDataIndex ) const
    {
        return mEnrichmentData( aEnrichmentDataIndex ).mSubphaseIndsInEnrichedBasis;
    }

    //-------------------------------------------------------------------------------------

    moris::Memory_Map
    Enrichment::get_memory_usage()
    {
        Memory_Map tMemoryMap;

        tMemoryMap.mMemoryMapData["mEnrichmentMethod"] = sizeof( mEnrichmentMethod );
        tMemoryMap.mMemoryMapData["mBasisRank"]        = sizeof( mBasisRank );
        tMemoryMap.mMemoryMapData["mMeshIndices"]      = mMeshIndices.capacity();
        tMemoryMap.mMemoryMapData["mNumBulkPhases"]    = sizeof( mNumBulkPhases );
        tMemoryMap.mMemoryMapData["mMeshIndices"]      = mMeshIndices.capacity();

        tMemoryMap.mMemoryMapData["mElementEnrichmentLevel"]      = moris::internal_capacity( mEnrichmentData( 0 ).mElementEnrichmentLevel );
        tMemoryMap.mMemoryMapData["mElementIndsInBasis"]          = moris::internal_capacity( mEnrichmentData( 0 ).mElementIndsInBasis );
        tMemoryMap.mMemoryMapData["mSubphaseIndsInEnrichedBasis"] = moris::internal_capacity( mEnrichmentData( 0 ).mSubphaseIndsInEnrichedBasis );
        tMemoryMap.mMemoryMapData["mBasisEnrichmentIndices"]      = moris::internal_capacity( mEnrichmentData( 0 ).mBasisEnrichmentIndices );
        tMemoryMap.mMemoryMapData["mEnrichedBasisIndexToId"]      = mEnrichmentData( 0 ).mEnrichedBasisIndexToId.capacity();
        tMemoryMap.mMemoryMapData["mSubphaseBGBasisIndices"]      = moris::internal_capacity( mEnrichmentData( 0 ).mSubphaseBGBasisIndices );
        tMemoryMap.mMemoryMapData["mSubphaseBGBasisEnrLev"]       = moris::internal_capacity( mEnrichmentData( 0 ).mSubphaseBGBasisEnrLev );
        tMemoryMap.mMemoryMapData["mNumEnrichmentLevels"]         = sizeof( mEnrichmentData( 0 ).mNumEnrichmentLevels );
        tMemoryMap.mMemoryMapData["mBGVertexInterpolations ptrs"] = mEnrichmentData( 0 ).mBGVertexInterpolations.capacity();

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
        if( !mMeshIndexMapIsSet )
        {
            // loop over the mesh indices and fill the map
            for ( uint iMeshListIndex = 0; iMeshListIndex < mMeshIndices.numel(); iMeshListIndex++)
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
        Cell< moris::Matrix< moris::IndexMat > > const & tSubphasesInEnrBasis = mEnrichmentData( aMeshIndex ).mSubphaseIndsInEnrichedBasis;

        std::ostringstream tStringStream;
        tStringStream << "Enriched_Basis_Id,";
        tStringStream << "Mesh_Index,";
        tStringStream << "Owner,";
        tStringStream << "PRank,";

        // global max size of
        moris_index tLocalMaxIGCellGroupSize = 0;
        for ( moris::uint iEB = 0; iEB < tSubphasesInEnrBasis.size(); iEB++ )
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
        // make sure mesh is not empty
        MORIS_ASSERT( mBackgroundMeshPtr->get_num_elems() > 0, "Enrichment::perform_basis_cluster_enrichment() - IP mesh without cells passed" );

        // construct cell in xtk conformal model neighborhood connectivity
        // todo: this does nothing right now, still needed?
        this->construct_neighborhoods();

        // make sure we have access to all the vertex interpolation
        this->setup_background_vertex_interpolations();

        // iterate through basis types (i.e. linear and quadratic)
        for ( moris::size_t iMeshIndex = 0; iMeshIndex < mMeshIndices.numel(); iMeshIndex++ )
        {
            // get the mesh index
            moris_index tMeshIndex = mMeshIndices( iMeshIndex );

            // log/trace enrichment for every discretization mesh index
            Tracer tTracer( "XTK", "Enrichment", "Mesh Index " + std::to_string( tMeshIndex ) );

            // Number of basis functions (= number of B-Splines on the A-mesh ?)
            moris::size_t tNumBasisFunctions = mBackgroundMeshPtr->get_num_basis_functions( tMeshIndex );

            // allocate member variables
            mEnrichmentData( tMeshIndex ).mElementEnrichmentLevel = moris::Cell< moris::Matrix< moris::IndexMat > >( tNumBasisFunctions );
            mEnrichmentData( tMeshIndex ).mElementIndsInBasis     = moris::Cell< moris::Matrix< moris::IndexMat > >( tNumBasisFunctions );

            // allocate data used after basis loop
            moris::Cell< moris::Matrix< moris::IndexMat > > tSubPhaseBinEnrichment( tNumBasisFunctions );
            moris::Cell< moris::Matrix< moris::IndexMat > > tSubphaseClusterIndicesInSupport( tNumBasisFunctions );
            moris::Cell< moris_index >                      tMaxEnrichmentLevel( tNumBasisFunctions, 0 );

            for ( moris::size_t iBasisFunction = 0; iBasisFunction < tNumBasisFunctions; iBasisFunction++ )
            {
                // Get elements in support of basis (these are interpolation cells)
                moris::Matrix< moris::IndexMat > tParentElementsInSupport;

                mBackgroundMeshPtr->get_elements_in_support_of_basis( tMeshIndex, iBasisFunction, tParentElementsInSupport );

                // get subphase clusters in support (separated by phase)
                tSubphaseClusterIndicesInSupport( iBasisFunction ) = this->get_subphase_clusters_in_support( tParentElementsInSupport );

                // construct subphase in support map
                IndexMap tSubPhaseIndexToSupportIndex;

                this->construct_subphase_in_support_map( tSubphaseClusterIndicesInSupport( iBasisFunction ), tSubPhaseIndexToSupportIndex );

                // prune the subphase to remove subphases outside of basis support
                moris::Matrix< moris::IndexMat > tPrunedSubphaseNeighborhood;

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

            MORIS_LOG_SPEC( "Num Enriched Basis", mEnrichmentData( tMeshIndex ).mNumEnrichmentLevels );
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

            // Number of basis functions (= number of B-Splines on the A-mesh ?)
            moris::size_t tNumBasisFunctions = mBackgroundMeshPtr->get_num_basis_functions( tMeshIndex );

            // allocate member variables
            mEnrichmentData( tMeshIndex ).mElementEnrichmentLevel = moris::Cell< moris::Matrix< moris::IndexMat > >( tNumBasisFunctions );
            mEnrichmentData( tMeshIndex ).mElementIndsInBasis     = moris::Cell< moris::Matrix< moris::IndexMat > >( tNumBasisFunctions );

            // allocate data used after basis loop
            moris::Cell< moris::Matrix< moris::IndexMat > > tSpgBinEnrichment( tNumBasisFunctions );
            moris::Cell< moris::Matrix< moris::IndexMat > > tSpgIndicesInSupport( tNumBasisFunctions );
            moris::Cell< moris_index >                      tMaxEnrichmentLevel( tNumBasisFunctions, 0 );

            for ( moris::size_t iBasisFunction = 0; iBasisFunction < tNumBasisFunctions; iBasisFunction++ )
            {
                // Get elements in support of basis (these are interpolation cells)
                moris::Matrix< moris::IndexMat > tParentElementsInSupport;

                mBackgroundMeshPtr->get_elements_in_support_of_basis( tMeshIndex, iBasisFunction, tParentElementsInSupport );

                // collect Subphase grounp indices in support
                IndexMap tSpgIndexToSupportIndex;
                this->get_subphase_groups_in_support( 
                    iMeshIndex, 
                    tParentElementsInSupport, 
                    tSpgIndicesInSupport( iBasisFunction ),
                    tSpgIndexToSupportIndex );

                // prune the subphase to remove subphases outside of basis support
                moris::Matrix< moris::IndexMat > tPrunedSpgNeighborhood;
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
            this->construct_enriched_basis_to_subphase_connectivity(
                    tMeshIndex,
                    tSpgBinEnrichment,
                    tSpgIndicesInSupport,
                    tMaxEnrichmentLevel );

            // FIXME: move communication to SPGs once everything is parallel consistent
            // Assign enriched basis indices Only indices here because interpolation cells needs basis
            // indices and basis ids are created using the interpolation cell ids
            this->assign_enriched_coefficients_identifiers(
                    tMeshIndex,
                    tMaxEnrichmentLevel );

            MORIS_LOG_SPEC( "Num Enriched Basis", mEnrichmentData( tMeshIndex ).mNumEnrichmentLevels );

            // find the SPs within one IP cell that are also within the same SPG 
            // NOTE: this is needed for cluster construction later; also needs to be repeated for every Mesh index, hence it's here
            this->establish_IP_SPG_SP_relationship( tMeshIndex );
        }

        // FIXME: enriched IP cells not parallel consistent yet
        // create the enriched interpolation mesh
        this->construct_enriched_interpolation_mesh_new();

        // create the enriched integration meshes with the clusters for every B-spline mesh
        for ( moris::uint iMeshIndexInList = 0; iMeshIndexInList < mMeshIndices.numel(); iMeshIndexInList++ )
        {
            this->construct_enriched_integration_mesh( (moris_index) iMeshIndexInList );
        }
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
        for ( moris::uint iM = 0; iM < mMeshIndices.numel(); iM++ )
        {
            // size the vertex interpolations (allocating everything as null)
            mEnrichmentData( mMeshIndices( iM ) ).mBGVertexInterpolations.resize( tNumVertices, nullptr );

            // iterate through vertices
            for ( moris::uint iV = 0; iV < tNumVertices; iV++ )
            {
                // access the vertex
                moris::mtk::Vertex& tVert = tBGIPMesh.get_mtk_vertex( (moris_index)iV );

                // if the vertex interpolation exits this vertex is own by the current processor
                // and the vertex interpolation is stored; otherwise the nullptr is kept indicating
                // that the vertex is not owned the by current processor
                if ( tVert.has_interpolation( mMeshIndices( iM ) ) )
                {
                    mEnrichmentData( mMeshIndices( iM ) ).mBGVertexInterpolations( iV ) = tVert.get_interpolation( mMeshIndices( iM ) );
                }
            }
        }
    }

    //-------------------------------------------------------------------------------------

    Matrix< IndexMat >
    Enrichment::get_subphase_clusters_in_support( moris::Matrix< moris::IndexMat > const & aElementsInSupport )
    {
        // count the number of subphase cluster in support
        moris::uint tCount = 0;

        for ( moris::size_t iE = 0; iE < aElementsInSupport.numel(); iE++ )
        {
            moris::Cell< moris_index > const & tSubphaseIndices = mCutIgMesh->get_parent_cell_subphases( aElementsInSupport( iE ) );

            tCount = tCount + tSubphaseIndices.size();
        }

        Matrix< IndexMat > tSubPhaseClusters( 1, tCount );
        tCount = 0;

        for ( moris::size_t iE = 0; iE < aElementsInSupport.numel(); iE++ )
        {

            moris::Cell< moris_index > const & tSubphaseIndices = mCutIgMesh->get_parent_cell_subphases( aElementsInSupport( iE ) );

            for ( moris::uint iSP = 0; iSP < tSubphaseIndices.size(); iSP++ )
            {
                tSubPhaseClusters( tCount++ ) = tSubphaseIndices( iSP );
            }
        }

        return tSubPhaseClusters;
    }

    //-------------------------------------------------------------------------------------

    void
    Enrichment::get_subphase_groups_in_support( 
        moris_index                             aMeshIndexPosition, 
        moris::Matrix< moris::IndexMat > const& aLagElementsInSupport,
        moris::Matrix< moris::IndexMat >&       aSubphaseGroupIndicesInSupport,
        IndexMap&                               aSubphaseGroupIndexToSupportIndex )
    {
        // initialize map between Lagrange an Bspline elements
        IndexMap tBspElemIndexMap;

        // intialize counter tracking the number of B-spline elements in support
        uint tBspElemCounter = 0;

        // go through Lag elements and find their corresponding B-Spline element
        for ( moris::size_t iLagElem = 0; iLagElem < aLagElementsInSupport.numel(); iLagElem++ )
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

        // initialize counter for subphase groups
        uint tSpgCounter = 0;

        // for each B-spline element get the SPGs and count them
        for( auto& iBspElem : tBspElemIndexMap )
        {
            // temporarily store current B-spline element's index
            moris_index tBspElemIndex = iBspElem.first;

            // get the SPGs associated with the current B-spline element
            moris::Cell< moris_index > const& tSPGsInBsplineElem = 
                mBsplineMeshInfos( aMeshIndexPosition )->get_SPG_indices_in_bspline_cell( tBspElemIndex );

            // store all SPGs in map
            for ( moris::size_t iSPG = 0; iSPG < tSPGsInBsplineElem.size(); iSPG++)
            {
                // save SPG index in map
                aSubphaseGroupIndexToSupportIndex[ tSPGsInBsplineElem( iSPG ) ] = tSpgCounter;

                // count SPGs in support
                tSpgCounter++;
            }           
        }

        // initialize linear list storing all SPGs in support with correct size
        aSubphaseGroupIndicesInSupport.resize( 1, tSpgCounter );

        // store all SPGs found in support in linear list
        for( auto& iSPG : aSubphaseGroupIndexToSupportIndex )
        {
            // temporarily store current SPG's index
            moris_index tSpgIndex = iSPG.first;
            moris_index tLocalSpgIndex = iSPG.second;
        
            // save SPG index in list
            aSubphaseGroupIndicesInSupport( tLocalSpgIndex ) = tSpgIndex;
        }
    }

    //-------------------------------------------------------------------------------------

    void
    Enrichment::construct_subphase_in_support_map(
            moris::Matrix< moris::IndexMat > const & aSubphaseClusterIndicesInSupport,
            IndexMap&                                aSubPhaseIndexToSupportIndex )
    {
        for ( moris::moris_index i = 0; i < (moris::moris_index)aSubphaseClusterIndicesInSupport.numel(); i++ )
        {
            aSubPhaseIndexToSupportIndex[aSubphaseClusterIndicesInSupport( i )] = i;
        }
    }

    //-------------------------------------------------------------------------------------

    void
    Enrichment::generate_pruned_subphase_graph_in_basis_support(
            moris::Matrix< moris::IndexMat > const & aSubphasesInSupport,
            IndexMap&                                aSubPhaseIndexToSupportIndex,
            moris::Matrix< moris::IndexMat >&        aPrunedSubPhaseToSubphase )
    {
        std::shared_ptr< Subphase_Neighborhood_Connectivity > tSubphaseNeighborhood = mCutIgMesh->get_subphase_neighborhood();

        // Construct full element neighbor graph in support and the corresponding shared faces
        aPrunedSubPhaseToSubphase.resize( aSubphasesInSupport.numel(), 50 );// FIXME: Keenan this allocation needs to done smarter
        aPrunedSubPhaseToSubphase.fill( MORIS_INDEX_MAX );

        // get subphase neighborhood information
        moris::Cell< std::shared_ptr< moris::Cell< moris_index > > > const & tSubphasetoSubphase = tSubphaseNeighborhood->mSubphaseToSubPhase;

        for ( moris::size_t i = 0; i < aSubphasesInSupport.numel(); i++ )
        {
            moris::Cell< moris_index > const & tSingleSubPhaseNeighbors = *tSubphasetoSubphase( aSubphasesInSupport( i ) );

            // iterate through and prune subphases not in support
            moris::uint tCount = 0;
            for ( moris::size_t j = 0; j < tSingleSubPhaseNeighbors.size(); j++ )
            {
                moris_index tNeighborSubphaseIndex = tSingleSubPhaseNeighbors( j );

                auto tNeighborIter = aSubPhaseIndexToSupportIndex.find( tNeighborSubphaseIndex );

                if ( tNeighborIter != aSubPhaseIndexToSupportIndex.end() )
                {
                    aPrunedSubPhaseToSubphase( i, tCount ) = tNeighborIter->second;
                    tCount++;
                }
            }
        }
    }

    //-------------------------------------------------------------------------------------

    void
    Enrichment::generate_pruned_subphase_group_graph_in_basis_support(
            const moris_index                       aMeshIndex,
            moris::Matrix< moris::IndexMat > const& aSubphaseGroupIndicesInSupport,
            IndexMap&                               aSubphaseGroupIndexToSupportIndex,
            moris::Matrix< moris::IndexMat >&       aPrunedSpgToSpg )
    {
        // get the global SPG connectivity graph 
        std::shared_ptr< Subphase_Neighborhood_Connectivity > tSubphaseGroupNeighborhood = mCutIgMesh->get_subphase_group_neighborhood( aMeshIndex );

        // Construct full SPG neighbor graph in support and the corresponding shared faces
        aPrunedSpgToSpg.resize( aSubphaseGroupIndicesInSupport.numel(), 50 ); // FIXME: this allocation needs to be done smarter
        aPrunedSpgToSpg.fill( MORIS_INDEX_MAX );

        // get subphase group neighborhood information
        moris::Cell< std::shared_ptr< moris::Cell< moris_index > > > const & tSpgToSpg = tSubphaseGroupNeighborhood->mSubphaseToSubPhase;

        for ( moris::size_t iSPG = 0; iSPG < aSubphaseGroupIndicesInSupport.numel(); iSPG++ )
        {
            // get temporariy list of Spg Neighbors 
            moris::Cell< moris_index > const & tSubphaseGroupNeighbors = *tSpgToSpg( aSubphaseGroupIndicesInSupport( iSPG ) );

            // iterate through and prune/remove subphases groups not in support
            moris::uint tCount = 0;
            for ( moris::size_t j = 0; j < tSubphaseGroupNeighbors.size(); j++ )
            {
                moris_index tNeighborSubphaseGroupIndex = tSubphaseGroupNeighbors( j );

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
            moris::Matrix< moris::IndexMat > const & aSubphasesInSupport,
            IndexMap&                                aSubPhaseIndexToSupportIndex,
            moris::Matrix< moris::IndexMat > const & aPrunedSubPhaseToSubphase,
            moris::Matrix< moris::IndexMat >&        aSubPhaseBinEnrichmentVals,
            moris_index&                             aMaxEnrichmentLevel )
    {
        // Variables needed for floodfill, consider removing these.
        // Active bins to include in floodfill (We include all bins)
        moris::Matrix< moris::IndexMat > tActiveBins( 1, aPrunedSubPhaseToSubphase.n_rows() );

        for ( moris::size_t i = 0; i < aPrunedSubPhaseToSubphase.n_rows(); i++ )
        {
            ( tActiveBins )( 0, i ) = i;
        }

        // Mark all as included
        moris::Matrix< moris::IndexMat > tIncludedBins( 1, aSubphasesInSupport.numel(), 1 );

        // Flood fill metric value (since all the subphases do not connect to dissimilar phases)
        moris::Matrix< moris::IndexMat > tDummyPhase( 1, aSubphasesInSupport.numel(), 1 );

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
            moris::Matrix< moris::IndexMat > const & aSpgsInSupport,
            moris::Matrix< moris::IndexMat > const & aPrunedSpgToSpg,
            moris::Matrix< moris::IndexMat >&        aSpgBinEnrichmentVals,
            moris_index&                             aMaxEnrichmentLevel )
    {
        // Variables needed for floodfill, consider removing these.
        // Active bins to include in floodfill (We include all bins)
        moris::Matrix< moris::IndexMat > tActiveBins( 1, aPrunedSpgToSpg.n_rows() );

        for ( moris::size_t i = 0; i < aPrunedSpgToSpg.n_rows(); i++ )
        {
            ( tActiveBins )( 0, i ) = i;
        }

        // Mark all as included
        moris::Matrix< moris::IndexMat > tIncludedBins( 1, aSpgsInSupport.numel(), 1 );

        // Flood fill metric value (since all the subphases do not connect to dissimilar phases)
        moris::Matrix< moris::IndexMat > tDummyPhase( 1, aSpgsInSupport.numel(), 1 );

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
            moris::Matrix< moris::IndexMat > const & aSubphasesInSupport,
            moris::Matrix< moris::IndexMat >&        aSubPhaseBinEnrichmentVals,
            moris_index const                        aMaxEnrichmentLevel )
    {
        // number of enrichment levels
        moris::uint tNumEnrichLevel = aMaxEnrichmentLevel + 1;

        // if only one or less enrichment levels used; no need to sort them
        if ( tNumEnrichLevel < 2 )
        {
            return;
        }

        // check that maximum enrichment level index is smaller than number of enrichment levels
        MORIS_ASSERT( (uint)aSubPhaseBinEnrichmentVals.max() < tNumEnrichLevel,
                "Enrichment::sort_enrichment_levels_in_basis_support - incorrect enrichment index." );

        // allocate list of centroids and volume for each enrichment level
        moris::Cell< moris::Matrix< DDRMat > > tCentroids( tNumEnrichLevel );
        moris::Cell< moris::real >             tTotalVolumes( tNumEnrichLevel );

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
            for ( moris::uint iSPCell = 1; iSPCell < tSubphaseIgCells->mIgCellGroup.size(); iSPCell++ )
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

        // sort centroids by cooridinates
        moris::Cell< moris_index > tSortingIndex = sort_points_by_coordinates( tCentroids );

        // build enrichment level old to new map
        moris::Cell< moris_index > tEnrichmentMap( tNumEnrichLevel );

        for ( uint i = 0; i < (uint)tNumEnrichLevel; ++i )
        {
            tEnrichmentMap( tSortingIndex( i ) ) = i;
        }

        // resort enrichment levels of subphases
        moris::Matrix< moris::IndexMat > tSortedEnrichmentVals( 1, aSubphasesInSupport.numel(), MORIS_SINT_MAX );

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
            moris_index const &                      aEnrichmentDataIndex,
            moris_index const &                      aBasisIndex,
            moris::Matrix< moris::IndexMat > const & aParentElementsInSupport,
            moris::Matrix< moris::IndexMat > const & aSubphasesInSupport,
            IndexMap&                                aSubPhaseIndexToSupportIndex,
            moris::Matrix< moris::IndexMat > const & aPrunedSubPhaseToSubphase,
            moris::Matrix< moris::IndexMat >&        aSubPhaseBinEnrichmentVals )
    {
        // resize member data
        moris::size_t tNumAllElementsInSupport                                         = this->count_elements_in_support( aParentElementsInSupport );
        mEnrichmentData( aEnrichmentDataIndex ).mElementIndsInBasis( aBasisIndex )     = moris::Matrix< moris::IndexMat >( 1, tNumAllElementsInSupport );
        mEnrichmentData( aEnrichmentDataIndex ).mElementEnrichmentLevel( aBasisIndex ) = moris::Matrix< moris::IndexMat >( 1, tNumAllElementsInSupport );

        moris::uint tCount = 0;

        // go through all subphases in support and save the curretn basis function's index to them, 
        // and which enrichment level of this basis function is active on them
        for ( moris::size_t iSP = 0; iSP < aSubphasesInSupport.numel(); iSP++ )
        {
            moris_index tSubphaseIndex = aSubphasesInSupport( iSP );

            // iterate through cells in the subphase
            std::shared_ptr< IG_Cell_Group > tSubphaseIgCells = mCutIgMesh->get_subphase_ig_cells( tSubphaseIndex );

            // iterate through cells in subphase
            for ( moris::uint iIgCell = 0; iIgCell < tSubphaseIgCells->mIgCellGroup.size(); iIgCell++ )
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
            moris_index const &                      aEnrichmentDataIndex,
            moris_index const &                      aBasisIndex,
            moris::Matrix< moris::IndexMat > const & aParentElementsInSupport,
            moris::Matrix< moris::IndexMat > const & aSpgsInSupport,
            moris::Matrix< moris::IndexMat >&        aSpgBinEnrichmentVals )
    {
        // resize member data
        moris::size_t tNumAllElementsInSupport                                         = this->count_elements_in_support( aParentElementsInSupport );
        mEnrichmentData( aEnrichmentDataIndex ).mElementIndsInBasis( aBasisIndex )     = moris::Matrix< moris::IndexMat >( 1, tNumAllElementsInSupport );
        mEnrichmentData( aEnrichmentDataIndex ).mElementEnrichmentLevel( aBasisIndex ) = moris::Matrix< moris::IndexMat >( 1, tNumAllElementsInSupport );

        // initialize counter tracking number of IG cells in support of basis function
        moris::uint tIgCellCount = 0;

        // go through all subphases in support and save the current basis function's index to them, 
        // and which enrichment level of this basis function is active on them
        for ( moris::size_t iSPG = 0; iSPG < aSpgsInSupport.numel(); iSPG++ )
        {
            // get the index of the SPG currently treated
            moris_index tSpgIndex = aSpgsInSupport( iSPG );

            // iterate through cells in the subphase
            moris::Cell< moris_index > tIgCellsInSPG = mCutIgMesh->get_ig_cells_in_SPG( this->get_list_index_for_mesh_index( aEnrichmentDataIndex ), tSpgIndex );

            // iterate through cells in subphase
            for ( moris::uint iIgCell = 0; iIgCell < tIgCellsInSPG.size(); iIgCell++ )
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
            moris_index const &                                     aEnrichmentDataIndex,
            moris::Cell< moris::Matrix< moris::IndexMat > > const & aSubPhaseBinEnrichment,
            moris::Cell< moris::Matrix< moris::IndexMat > > const & aSubphaseClusterIndicesInSupport,
            moris::Cell< moris_index > const &                      aMaxEnrichmentLevel )
    {
        // intialize counter for total number of enriched BFs
        moris::moris_index tNumEnrichmentBasis = 0;

        // count up the number of enriched basis functions
        for ( moris::uint i = 0; i < aMaxEnrichmentLevel.size(); i++ )
        {
            tNumEnrichmentBasis += aMaxEnrichmentLevel( i ) + 1;
        }

        // size data
        mEnrichmentData( aEnrichmentDataIndex ).mSubphaseIndsInEnrichedBasis.resize( tNumEnrichmentBasis );

        moris_index tBaseIndex = 0;

        for ( moris::uint iEnrLvl = 0; iEnrLvl < aSubPhaseBinEnrichment.size(); iEnrLvl++ )
        {
            // get the maximum enrichment level in this basis support
            moris::moris_index tMaxEnrLev = aMaxEnrichmentLevel( iEnrLvl );

            // counter
            Cell< moris_index > tCounter( tMaxEnrLev + 1, 0 );

            // allocate member data for these basis functions
            for ( moris::moris_index iEnr = tBaseIndex; iEnr < tBaseIndex + tMaxEnrLev + 1; iEnr++ )
            {
                mEnrichmentData( aEnrichmentDataIndex ).mSubphaseIndsInEnrichedBasis( iEnr ).resize( 1, aSubphaseClusterIndicesInSupport( iEnrLvl ).numel() );
            }

            // iterate through subphases in support and add them to appropriate location in mSubphaseIndsInEnrichedBasis
            for ( moris::uint iSp = 0; iSp < aSubphaseClusterIndicesInSupport( iEnrLvl ).numel(); iSp++ )
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
            for ( moris::moris_index iEnr = 0; iEnr < tMaxEnrLev + 1; iEnr++ )
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
            moris_index const &                                     aEnrichmentDataIndex,
            moris::Cell< moris::Matrix< moris::IndexMat > > const & aSpgBinEnrichment,
            moris::Cell< moris::Matrix< moris::IndexMat > > const & aSpgIndicesInSupport,
            moris::Cell< moris_index > const &                      aMaxEnrichmentLevel )
    {
        // intialize counter for total number of enriched BFs
        moris::moris_index tNumEnrichmentBasis = 0;

        // count up the number of enriched basis functions
        for ( moris::uint i = 0; i < aMaxEnrichmentLevel.size(); i++ )
        {
            tNumEnrichmentBasis += aMaxEnrichmentLevel( i ) + 1;
        }

        // size data
        mEnrichmentData( aEnrichmentDataIndex ).mSubphaseGroupIndsInEnrichedBasis.resize( tNumEnrichmentBasis );

        moris_index tBaseIndex = 0;

        for ( moris::uint iEnrLvl = 0; iEnrLvl < aSpgBinEnrichment.size(); iEnrLvl++ )
        {
            // get the maximum enrichment level in this basis support
            moris::moris_index tMaxEnrLev = aMaxEnrichmentLevel( iEnrLvl );

            // counter
            Cell< moris_index > tCounter( tMaxEnrLev + 1, 0 );

            // allocate member data for these basis functions
            for ( moris::moris_index iEnr = tBaseIndex; iEnr < tBaseIndex + tMaxEnrLev + 1; iEnr++ )
            {
                mEnrichmentData( aEnrichmentDataIndex ).mSubphaseGroupIndsInEnrichedBasis( iEnr ).resize( 1, aSpgIndicesInSupport( iEnrLvl ).numel() );
            }

            // iterate through SPGs in support and add them to appropriate location in mSubphaseGroupIndsInEnrichedBasis
            for ( moris::uint iSPG = 0; iSPG < aSpgIndicesInSupport( iEnrLvl ).numel(); iSPG++ )
            {
                // get cluster enrichment level
                moris_index tClusterEnrLev = aSpgBinEnrichment( iEnrLvl )( iSPG );

                // add to the member data
                mEnrichmentData( aEnrichmentDataIndex ).mSubphaseGroupIndsInEnrichedBasis( tBaseIndex + tClusterEnrLev )( tCounter( tClusterEnrLev ) ) =
                        aSpgIndicesInSupport( iEnrLvl )( iSPG );

                // increment count
                tCounter( tClusterEnrLev )++;
            }

            // size out unused space
            for ( moris::moris_index iEnr = 0; iEnr < tMaxEnrLev + 1; iEnr++ )
            {
                moris_index tIndex = tBaseIndex + iEnr;
                mEnrichmentData( aEnrichmentDataIndex ).mSubphaseGroupIndsInEnrichedBasis( tIndex ).resize( 1, tCounter( iEnr ) );

                // sort in ascending order (easier to find in MPI)
                // if this sort is removed the function  subphase_is_in_support needs to be updated
                // FIXME: does this function need to be updated for SPGs?
                moris::sort( mEnrichmentData( aEnrichmentDataIndex ).mSubphaseGroupIndsInEnrichedBasis( tIndex ),
                        mEnrichmentData( aEnrichmentDataIndex ).mSubphaseGroupIndsInEnrichedBasis( tIndex ),
                        "ascend",
                        1 );
            }

            // update starting index
            tBaseIndex = tBaseIndex + tMaxEnrLev + 1;
        }
    }

    //-------------------------------------------------------------------------------------

    void
    Enrichment::assign_enriched_coefficients_identifiers(
            moris_index const &                aEnrichmentDataIndex,
            moris::Cell< moris_index > const & aMaxEnrichmentLevel )
    {
        // get number of non-enriched BFs
        uint tNumBasis = mEnrichmentData( aEnrichmentDataIndex ).mElementIndsInBasis.size();
        
        mEnrichmentData( aEnrichmentDataIndex ).mBasisEnrichmentIndices.resize( tNumBasis );

        mEnrichmentData( aEnrichmentDataIndex ).mNumEnrichmentLevels = 0;

        for ( moris::uint i = 0; i < tNumBasis; i++ )
        {
            moris::moris_index tMaxEnrLev                                        = aMaxEnrichmentLevel( i ) + 1;
            mEnrichmentData( aEnrichmentDataIndex ).mNumEnrichmentLevels         = mEnrichmentData( aEnrichmentDataIndex ).mNumEnrichmentLevels + tMaxEnrLev;
            mEnrichmentData( aEnrichmentDataIndex ).mBasisEnrichmentIndices( i ) = moris::Matrix< moris::IndexMat >( 1, tMaxEnrLev );
        }

        // allocate enriched basis index to id data
        mEnrichmentData( aEnrichmentDataIndex ).mEnrichedBasisIndexToId.resize( 1, mEnrichmentData( aEnrichmentDataIndex ).mNumEnrichmentLevels );
        mEnrichmentData( aEnrichmentDataIndex ).mEnrichedBasisIndexToId.fill( MORIS_INDEX_MAX );

        moris_index tParRank = par_rank();

        // get the comm table
        Matrix< IndexMat > tCommTable = mXTKModelPtr->get_communication_table();

        // Procs CEll
        Cell< moris_index > tProcRanks( tCommTable.numel() );

        std::unordered_map< moris_id, moris_id > tProcRankToIndexInData;

        // resize proc ranks and setup map to comm table
        for ( moris::uint i = 0; i < tCommTable.numel(); i++ )
        {
            tProcRankToIndexInData[tCommTable( i )] = i;
            tProcRanks( i )                         = ( tCommTable( i ) );
        }

        // first index and id ( not first gets background basis information)
        moris::moris_index tIndOffset     = 0;
        moris::moris_id    tBasisIdOffset = this->allocate_basis_ids( aEnrichmentDataIndex, mEnrichmentData( aEnrichmentDataIndex ).mNumEnrichmentLevels );

// TODO: once the SPGs have IDs and parallel consistent, this here needs to be changed to use SPG IDs
        Cell< Cell< moris_index > > tBasisIdToBasisOwner( tCommTable.numel() );
        Cell< Cell< moris_index > > tSubphaseIdInSupport( tCommTable.numel() );
        //Cell< Cell< moris_index > > tSubphaseGroupIdInSupport( tCommTable.numel() );
        Cell< Cell< moris_index > > tBasisIndexToBasisOwner( tCommTable.numel() );

        for ( moris::uint i = 0; i < mEnrichmentData( aEnrichmentDataIndex ).mBasisEnrichmentIndices.size(); i++ )
        {
            moris_index tOwner = mBackgroundMeshPtr->get_entity_owner(
                    i,
                    mBasisRank,
                    aEnrichmentDataIndex );

            moris_id tBackBasisId =
                    mBackgroundMeshPtr->get_glb_entity_id_from_entity_loc_index(
                            i,
                            mBasisRank,
                            aEnrichmentDataIndex );

            moris_index tProcDataIndex = tProcRankToIndexInData[tOwner];

            // always set indices
            moris::Matrix< moris::IndexMat >& tBasisEnrichmentInds = mEnrichmentData( aEnrichmentDataIndex ).mBasisEnrichmentIndices( i );

            for ( moris::uint j = 0; j < tBasisEnrichmentInds.numel(); j++ )
            {
                tBasisEnrichmentInds( j ) = tIndOffset;
                tIndOffset++;
            }

            // only set id if we own it and package data for communication if shared
            if ( tOwner == tParRank )
            {
                mEnrichmentData( aEnrichmentDataIndex ).mEnrichedBasisIndexToId( tBasisEnrichmentInds( 0 ) ) = tBackBasisId;

                for ( moris::uint j = 1; j < tBasisEnrichmentInds.numel(); j++ )
                {
                    MORIS_ASSERT( mEnrichmentData( aEnrichmentDataIndex ).mEnrichedBasisIndexToId( tBasisEnrichmentInds( j ) ) == MORIS_INDEX_MAX, "Already set enriched basis id" );

                    mEnrichmentData( aEnrichmentDataIndex ).mEnrichedBasisIndexToId( tBasisEnrichmentInds( j ) ) = tBasisIdOffset;

                    tBasisIdOffset++;
                }
            }

            // if we don't own the basis setup the communication to get the basis
            else
            {
                for ( moris::uint j = 0; j < tBasisEnrichmentInds.numel(); j++ )
                {
                    moris_index tEnrichedBasisIndex = tBasisEnrichmentInds( j );

// TODO: once the SPGs have IDs and parallel consistent, this here needs to be changed to use SPGs
                    moris_index tFirstSubphaseInSupportIndex = mEnrichmentData( aEnrichmentDataIndex ).mSubphaseIndsInEnrichedBasis( tEnrichedBasisIndex )( 0 );
                    //moris_index tFirstSpgInSupportIndex = mEnrichmentData( aEnrichmentDataIndex ).mSubphaseGroupIndsInEnrichedBasis( tEnrichedBasisIndex )( 0 );

                    moris_index tFirstSubphaseInSupportId = mXTKModelPtr->get_subphase_id( tFirstSubphaseInSupportIndex );
                    //moris_index tFirstSubphaseGroupInSupportId = mXTKModelPtr->get_subphase_group_id( tFirstSpgInSupportIndex );

                    tBasisIdToBasisOwner( tProcDataIndex ).push_back( tBackBasisId );
                    tSubphaseIdInSupport( tProcDataIndex ).push_back( tFirstSubphaseInSupportId );
                    // tSubphaseGroupIdInSupport( tProcDataIndex ).push_back( tFirstSubphaseGroupInSupportId );
                    tBasisIndexToBasisOwner( tProcDataIndex ).push_back( tEnrichedBasisIndex );
                }
            }
        }

        // send information about not owned enriched basis to owner processor
        Cell< moris::Matrix< moris::IndexMat > > tNotOwnedEnrichedBasisId;

        this->communicate_basis_information_with_owner(
                aEnrichmentDataIndex,
                tBasisIdToBasisOwner,
                tSubphaseIdInSupport, // tSubphaseGroupIdInSupport,
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
    Enrichment::communicate_basis_information_with_owner(
            moris_index const &                       aEnrichmentDataIndex,
            Cell< Cell< moris_index > > const &       aBasisIdToBasisOwner,
            Cell< Cell< moris_index > > const &       aSubphaseIdInSupport,
            Cell< moris_index > const &               aProcRanks,
            std::unordered_map< moris_id, moris_id >& aProcRankToIndexInData,
            Cell< moris::Matrix< moris::IndexMat > >& aEnrichedBasisId )
    {
        // copy into a matrix
        Cell< moris::Matrix< moris::IndexMat > > tBasisIdToBasisOwnerMat( aProcRanks.size() );
        Cell< moris::Matrix< moris::IndexMat > > tSubphaseIdInSupport( aProcRanks.size() );

        for ( moris::uint i = 0; i < aBasisIdToBasisOwner.size(); i++ )
        {
            tBasisIdToBasisOwnerMat( i ) = Matrix< IndexMat >( 1, aBasisIdToBasisOwner( i ).size() );
            tSubphaseIdInSupport( i )    = Matrix< IndexMat >( 1, aBasisIdToBasisOwner( i ).size() );

            for ( moris::uint j = 0; j < aBasisIdToBasisOwner( i ).size(); j++ )
            {
                tBasisIdToBasisOwnerMat( i )( j ) = aBasisIdToBasisOwner( i )( j );
                tSubphaseIdInSupport( i )( j )    = aSubphaseIdInSupport( i )( j );
            }

            if ( aBasisIdToBasisOwner( i ).size() == 0 )
            {
                tBasisIdToBasisOwnerMat( i ).resize( 1, 1 );
                tBasisIdToBasisOwnerMat( i )( 0 ) = MORIS_INDEX_MAX;

                tSubphaseIdInSupport( i ).resize( 1, 1 );
                tSubphaseIdInSupport( i )( 0 ) = MORIS_INDEX_MAX;
            }
        }

        // send to processor (no particular reason for these tag number choices)
        moris_index tBasisIdsTag    = 1119;
        moris_index tMaxSubphasetag = 1120;
        moris_index tReturnIdTag    = 1121;

        for ( moris::uint i = 0; i < tBasisIdToBasisOwnerMat.size(); i++ )
        {
            // send child element ids
            nonblocking_send(
                    tBasisIdToBasisOwnerMat( i ),
                    tBasisIdToBasisOwnerMat( i ).n_rows(),
                    tBasisIdToBasisOwnerMat( i ).n_cols(),
                    aProcRanks( i ),
                    tBasisIdsTag );

            // send element parent ids
            nonblocking_send(
                    tSubphaseIdInSupport( i ),
                    tSubphaseIdInSupport( i ).n_rows(),
                    tSubphaseIdInSupport( i ).n_cols(),
                    aProcRanks( i ),
                    tMaxSubphasetag );
        }

        // make sure that all sends have been completed
        barrier();

        // receive
        moris::moris_index tNumRow = 1;

        Cell< moris::Matrix< IdMat > > tReceiveInfoBasisId( aProcRanks.size() );
        Cell< moris::Matrix< IdMat > > tReceivedSubphaseId( aProcRanks.size() );
        Cell< moris::Matrix< IdMat > > tEnrichedBasisIds( aProcRanks.size() );
        aEnrichedBasisId.resize( aProcRanks.size() );

        for ( moris::uint i = 0; i < aProcRanks.size(); i++ )
        {
            tReceiveInfoBasisId( i ).resize( 1, 1 );
            tReceivedSubphaseId( i ).resize( 1, 1 );

            receive( tReceiveInfoBasisId( i ), tNumRow, aProcRanks( i ), tBasisIdsTag );
            receive( tReceivedSubphaseId( i ), tNumRow, aProcRanks( i ), tMaxSubphasetag );

            tEnrichedBasisIds( i ).resize( tReceiveInfoBasisId( i ).n_rows(), tReceiveInfoBasisId( i ).n_cols() );
        }

        // iterate through received information and setup response information
        for ( moris::uint i = 0; i < aProcRanks.size(); i++ )
        {
            moris_id    tBasisId    = 0;
            moris_id    tBasisIndex = 0;
            moris_index tCount      = 0;
            aEnrichedBasisId( i ).resize( tReceiveInfoBasisId( i ).n_rows(), tReceiveInfoBasisId( i ).n_cols() );

            if ( tReceiveInfoBasisId( i )( 0 ) != MORIS_INDEX_MAX )
            {
                for ( moris::uint j = 0; j < tReceiveInfoBasisId( i ).n_cols(); j++ )
                {
                    if ( tBasisId != tReceiveInfoBasisId( i )( j ) )
                    {
                        tBasisId    = tReceiveInfoBasisId( i )( j );
                        tBasisIndex = mBackgroundMeshPtr->get_loc_entity_ind_from_entity_glb_id(
                                tBasisId,
                                mBasisRank,
                                aEnrichmentDataIndex );
                    }

                    moris_id    tSubphaseId    = tReceivedSubphaseId( i )( j );
                    moris_index tSubphaseIndex = mXTKModelPtr->get_subphase_index( tSubphaseId );

                    // iterate through the max integer ids and match

                    bool tFound = false;
                    for ( moris::uint k = 0; k < mEnrichmentData( aEnrichmentDataIndex ).mBasisEnrichmentIndices( tBasisIndex ).numel(); k++ )
                    {
                        if ( this->subphase_is_in_support( aEnrichmentDataIndex, tSubphaseIndex, mEnrichmentData( aEnrichmentDataIndex ).mBasisEnrichmentIndices( tBasisIndex )( k ) ) )
                        {
                            moris_index tEnrichedBasisIndex = mEnrichmentData( aEnrichmentDataIndex ).mBasisEnrichmentIndices( tBasisIndex )( k );
                            moris_index tEnrichedBasisId    = mEnrichmentData( aEnrichmentDataIndex ).mEnrichedBasisIndexToId( tEnrichedBasisIndex );

                            tEnrichedBasisIds( i )( tCount ) = tEnrichedBasisId;
                            tCount++;
                            tFound = true;
                        }
                    }
                    MORIS_ERROR( tFound, "Basis not found" );
                }
            }
            else
            {
                tEnrichedBasisIds( i )( 0 ) = MORIS_INDEX_MAX;
            }

            // send enriched ids
            nonblocking_send(
                    tEnrichedBasisIds( i ),
                    tEnrichedBasisIds( i ).n_rows(),
                    tEnrichedBasisIds( i ).n_cols(),
                    aProcRanks( i ),
                    tReturnIdTag );
        }

        barrier();

        // receive information
        for ( moris::uint i = 0; i < aProcRanks.size(); i++ )
        {
            aEnrichedBasisId( i ).resize( 1, 1 );
            receive( aEnrichedBasisId( i ), tNumRow, aProcRanks( i ), tReturnIdTag );
        }

        barrier();
    }

    //-------------------------------------------------------------------------------------

    void
    Enrichment::set_received_enriched_basis_ids(
            moris_index const &                              aEnrichmentDataIndex,
            Cell< moris::Matrix< moris::IndexMat > > const & aReceivedEnrichedIds,
            Cell< Cell< moris_index > > const &              aBasisIndexToBasisOwner )
    {
        for ( moris::uint i = 0; i < aReceivedEnrichedIds.size(); i++ )
        {
            if ( aReceivedEnrichedIds( i )( 0 ) != MORIS_INDEX_MAX )
            {
                MORIS_ASSERT( aReceivedEnrichedIds( i ).numel() == aBasisIndexToBasisOwner( i ).size(),
                        "Dimension mismatch between received information and expected information" );

                for ( moris::uint j = 0; j < aReceivedEnrichedIds( i ).numel(); j++ )
                {
                    moris_index tLocalBasisIndex = aBasisIndexToBasisOwner( i )( j );
                    moris_index tGlobaId         = aReceivedEnrichedIds( i )( j );

                    MORIS_ASSERT( mEnrichmentData( aEnrichmentDataIndex ).mEnrichedBasisIndexToId( tLocalBasisIndex ) == MORIS_INDEX_MAX,
                            "Id already set for this basis function" );

                    mEnrichmentData( aEnrichmentDataIndex ).mEnrichedBasisIndexToId( tLocalBasisIndex ) = tGlobaId;
                }
            }
        }
    }

    //-------------------------------------------------------------------------------------

    moris::size_t
    Enrichment::count_elements_in_support( moris::Matrix< moris::IndexMat > const & aParentElementsInSupport )
    {

        // Number of elements in this support (need both parent and total)
        moris::size_t tNumParentElementsInSupport = aParentElementsInSupport.n_cols();
        moris::size_t tNumElementsInSupport       = 0;

        // Count children elements in support
        for ( moris::size_t i = 0; i < tNumParentElementsInSupport; i++ )
        {
            moris::Cell< moris_index > const & tSubphaseIndices = mCutIgMesh->get_parent_cell_subphases( aParentElementsInSupport( i ) );

            // iterate through subphases
            for ( moris::uint iSP = 0; iSP < tSubphaseIndices.size(); iSP++ )
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
        for ( moris::uint i = 0; i < mEnrichmentData( aEnrichmentDataIndex ).mSubphaseIndsInEnrichedBasis( aEnrichedBasisIndex ).numel(); i++ )
        {
            if ( mEnrichmentData( aEnrichmentDataIndex ).mSubphaseIndsInEnrichedBasis( aEnrichedBasisIndex )( i ) == aSubphaseIndex )
            {
                return true;
            }
            else if ( mEnrichmentData( aEnrichmentDataIndex ).mSubphaseIndsInEnrichedBasis( aEnrichedBasisIndex )( i ) > aSubphaseIndex )
            {
                return false;
            }
        }

        return false;
    }

    //-------------------------------------------------------------------------------------

    void
    Enrichment::print_basis_support_debug(
            moris_index                              aBasisIndex,
            moris::Matrix< moris::IndexMat > const & aParentElementsInSupport,
            moris::Matrix< moris::IndexMat > const & aSubphasesInSupport,
            IndexMap&                                aSubPhaseIndexToSupportIndex,
            moris::Matrix< moris::IndexMat > const & aPrunedSubPhaseToSubphase,
            moris::Matrix< moris::IndexMat >&        aSubPhaseBinEnrichmentVals )
    {
        std::cout << "--------------------------------------------------" << std::endl;
        std::cout << "Basis Index: " << aBasisIndex << std::endl;
        std::cout << "Parent Cell In Support:";
        for ( moris::uint i = 0; i < aParentElementsInSupport.numel(); i++ )
        {
            std::cout << std::setw( 8 ) << aParentElementsInSupport( i );
        }
        std::cout << "\nSubphases In Support:";
        for ( moris::uint i = 0; i < aSubphasesInSupport.numel(); i++ )
        {
            std::cout << std::setw( 8 ) << aSubphasesInSupport( i );
        }

        std::cout << "\nSubphase Neighborhood In Support:" << std::endl;
        for ( moris::uint i = 0; i < aPrunedSubPhaseToSubphase.n_rows(); i++ )
        {
            std::cout << std::setw( 6 ) << aSubphasesInSupport( i ) << " | ";

            for ( moris::uint j = 0; j < aPrunedSubPhaseToSubphase.n_cols(); j++ )
            {
                if ( aPrunedSubPhaseToSubphase( i, j ) != MORIS_INDEX_MAX )
                {
                    std::cout << std::setw( 6 ) << aSubphasesInSupport( aPrunedSubPhaseToSubphase( i, j ) );
                }
            }
            std::cout << std::endl;
        }

        std::cout << "\nSubphase Enrichment Level: \n";
        for ( moris::uint i = 0; i < aSubPhaseBinEnrichmentVals.numel(); i++ )
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
        Tracer tTracer( "XTK", "Enrich", "Construct Enriched Interpolation Mesh" );

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

        mXTKModelPtr->mEnrichedInterpMesh( 0 )->mCoeffToEnrichCoeffs.resize( mMeshIndices.max() + 1 );
        mXTKModelPtr->mEnrichedInterpMesh( 0 )->mEnrichCoeffLocToGlob.resize( mMeshIndices.max() + 1 );

        // enriched to non-enriched coefficient mapping for all B-spline meshes
        for ( moris::uint iMesh = 0; iMesh < mMeshIndices.numel(); iMesh++ )
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

        // initialize local to global maps
        mXTKModelPtr->mEnrichedInterpMesh( 0 )->mLocalToGlobalMaps = Cell< Matrix< IdMat > >( 4 );
        mXTKModelPtr->mEnrichedInterpMesh( 0 )->mGlobaltoLobalMaps = Cell< std::unordered_map< moris_id, moris_index > >( 4 );
        mXTKModelPtr->mEnrichedInterpMesh( 0 )->setup_cell_maps();
        mXTKModelPtr->mEnrichedInterpMesh( 0 )->setup_basis_maps();
        mXTKModelPtr->mEnrichedInterpMesh( 0 )->assign_ip_vertex_ids();
        mXTKModelPtr->mEnrichedInterpMesh( 0 )->setup_vertex_maps();
        // FIXME: aren't the above steps already done in the finalize_setup() call in Enriched_Interpolation_Mesh::setup_local_to_global_maps()?

        // moris::Cell< mtk::Vertex* > tVerticesToCommunicate;
        // for(auto & iVert: mXTKModelPtr->mEnrichedInterpMesh(0)->mEnrichedInterpVerts)
        // {
        //     if(iVert->get_owner() != moris::par_rank() )
        //     {
        //         tVerticesToCommunicate.push_back(iVert);
        //     }
        // }
        // 
        // mXTKModelPtr->mEnrichedInterpMesh(0)->communicate_select_vertex_interpolation(tVerticesToCommunicate);

        // in most cases all the interpolation vertices are the same. We merge them back together with this call
        // post-processing to construct_enriched_interpolation_vertices_and_cells in an effort to not add complexity to the function 
        // (as that function is already too complex/loaded)
        mXTKModelPtr->mEnrichedInterpMesh( 0 )->merge_duplicate_interpolation_vertices();

        // reset global to local maps (delete & setup again) with deleted duplicate vertices
        mXTKModelPtr->mEnrichedInterpMesh( 0 )->mGlobaltoLobalMaps( 0 ).clear();
        mXTKModelPtr->mEnrichedInterpMesh( 0 )->setup_vertex_maps();
    }

    //-------------------------------------------------------------------------------------

    void
    Enrichment::construct_enriched_interpolation_mesh_new()
    {
        // log/trace this function
        Tracer tTracer( "XTK", "Enrich", "Construct Enriched Interpolation Mesh" );

        // initialize a new enriched interpolation mesh
        mXTKModelPtr->mEnrichedInterpMesh( 0 ) = new Enriched_Interpolation_Mesh( mXTKModelPtr );

        // set enriched basis rank
        mXTKModelPtr->mEnrichedInterpMesh( 0 )->mBasisRank = mBasisRank;

        // set mesh index
        mXTKModelPtr->mEnrichedInterpMesh( 0 )->mMeshIndices = mMeshIndices;
        mXTKModelPtr->mEnrichedInterpMesh( 0 )->setup_mesh_index_map();

        // allocate memory for enriched interpolation cells
        this->allocate_interpolation_cells_new();

        // unzip all IP cells and vertices and construct the whole enr. IP mesh
        /* Note: This constructs all unzipped (i.e. enriched) IP cells with unzipped interpolation vertices being attached to a single cell
         * this also handles the case of multiple enrichments where the number of interpolation vertices vary */
        this->construct_enriched_interpolation_vertices_and_cells_new();

        mXTKModelPtr->mEnrichedInterpMesh( 0 )->mCoeffToEnrichCoeffs.resize( mMeshIndices.max() + 1 );
        mXTKModelPtr->mEnrichedInterpMesh( 0 )->mEnrichCoeffLocToGlob.resize( mMeshIndices.max() + 1 );

        // enriched to non-enriched coefficient mapping for all B-spline meshes
        for ( moris::uint iMesh = 0; iMesh < mMeshIndices.numel(); iMesh++ )
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

        // initialize local to global maps
        mXTKModelPtr->mEnrichedInterpMesh( 0 )->mLocalToGlobalMaps = Cell< Matrix< IdMat > >( 4 );
        mXTKModelPtr->mEnrichedInterpMesh( 0 )->mGlobaltoLobalMaps = Cell< std::unordered_map< moris_id, moris_index > >( 4 );
        mXTKModelPtr->mEnrichedInterpMesh( 0 )->setup_cell_maps();
        mXTKModelPtr->mEnrichedInterpMesh( 0 )->setup_basis_maps();
        mXTKModelPtr->mEnrichedInterpMesh( 0 )->assign_ip_vertex_ids();
        mXTKModelPtr->mEnrichedInterpMesh( 0 )->setup_vertex_maps();
        // FIXME: aren't the above steps already done in the finalize_setup() call in Enriched_Interpolation_Mesh::setup_local_to_global_maps()?

        // moris::Cell< mtk::Vertex* > tVerticesToCommunicate;
        // for(auto & iVert: mXTKModelPtr->mEnrichedInterpMesh(0)->mEnrichedInterpVerts)
        // {
        //     if(iVert->get_owner() != moris::par_rank() )
        //     {
        //         tVerticesToCommunicate.push_back(iVert);
        //     }
        // }
        // 
        // mXTKModelPtr->mEnrichedInterpMesh(0)->communicate_select_vertex_interpolation(tVerticesToCommunicate);

        // in most cases all the interpolation vertices are the same. We merge them back together with this call
        // post-processing to construct_enriched_interpolation_vertices_and_cells in an effort to not add complexity to the function 
        // (as that function is already too complex/loaded)
        mXTKModelPtr->mEnrichedInterpMesh( 0 )->merge_duplicate_interpolation_vertices();

        // reset global to local maps (delete & setup again) with deleted duplicate vertices
        mXTKModelPtr->mEnrichedInterpMesh( 0 )->mGlobaltoLobalMaps( 0 ).clear();
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
    Enrichment::construct_enriched_integration_mesh( const moris_index aMeshIndexInList )
    {
        MORIS_ASSERT( mXTKModelPtr->mEnrichedInterpMesh( 0 ) != nullptr,
                "Enrichment::construct_enriched_integration_mesh_new() - No enriched interpolation mesh to link enriched integration mesh to" );

        MORIS_ASSERT( mXTKModelPtr->mEnrichedIntegMesh.size() > (uint) aMeshIndexInList, 
                "Enrichment::construct_enriched_integration_mesh_new() - Place for Enriched IG mesh doesn't exist in list of enriched IG meshes." );

        mXTKModelPtr->mEnrichedIntegMesh( aMeshIndexInList ) = new Enriched_Integration_Mesh( mXTKModelPtr, aMeshIndexInList );
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
        tEnrInterpMesh->mBaseCelltoEnrichedCell.resize( mBackgroundMeshPtr->get_num_elems() );
    }

    //-------------------------------------------------------------------------------------

    void
    Enrichment::allocate_interpolation_cells_new()
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
        mNumUnzippingsOnIpCell.resize( tNumIpElems );
        mEnrIpCellIndices.resize( tNumIpElems );

        // loop to obtain the number of max possible number of enr. IP cells
        for ( moris::size_t iIpCell = 0; iIpCell < tNumIpElems; iIpCell++)
        {
            // get number of times this IP cell needs to be unzipped
            uint tNumEnrLvlsOnIpCell = this->maximum_number_of_unzippings_for_IP_cell( iIpCell );

            // store it in list
            mNumUnzippingsOnIpCell( iIpCell ) = tNumEnrLvlsOnIpCell;

            // correct size for list of enr. IP cell indices on current Ip cell
            mEnrIpCellIndices( iIpCell ).resize( tNumEnrLvlsOnIpCell );

            for ( moris::uint iEnrLvl = 0; iEnrLvl < tNumEnrLvlsOnIpCell; iEnrLvl++)
            {
                // fill map relating IP-cell index and local enr. lvl to the index of the resulting enr. IP cell
                mEnrIpCellIndices( iIpCell )( iEnrLvl ) = (moris_index) tEnrIpCellCounter;

                // count total number of enr. IP cells
                tEnrIpCellCounter ++;
            } 
        }
        
        // store number of Enr IP cells
        mNumEnrIpCells = tEnrIpCellCounter;

        // there is one interpolation cell per subphase
        tEnrInterpMesh->mEnrichedInterpCells.resize( mNumEnrIpCells );

        // assuming all interpolation cells are the same
        // figure out how many vertices there are per interpolation cell
        uint tNumVertsPerCell = 0;

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
        tEnrInterpMesh->mBaseCelltoEnrichedCell.resize( mBackgroundMeshPtr->get_num_elems() );
    }

    //-------------------------------------------------------------------------------------

    uint
    Enrichment::maximum_number_of_unzippings_for_IP_cell( moris_index aIpCellIndex )
    {
        // initialize counter for number of unzippings
        uint tMaxNumUnzippings = 0;

        // for each DMI get how often the IP cell needs to be unzipped
        for ( moris::size_t iMeshIndex = 0; iMeshIndex < mMeshIndices.numel(); iMeshIndex++)
        {
            // get number of times IP cell would need to get unzipped for current mesh index
            uint tNumUnzippingsOnMeshIndex = mBsplineMeshInfos( iMeshIndex )->get_num_SPGs_associated_with_extraction_cell( aIpCellIndex );
            
            // overwrite max if new number of unzippings is greater
            tMaxNumUnzippings = std::max( tMaxNumUnzippings, tNumUnzippingsOnMeshIndex );
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
        for( uint iBspElem = 0; iBspElem < tNumBspElems; iBspElem++ )
        {
            // initialize map of all SPs on BSp element
            IndexMap tSpToIpElem;

            // get number of IP cells in current B-spline element
            uint tNumIpCellsInBsplineElement = tBsplineMeshInfo->mExtractionCellsIndicesInBsplineCells( iBspElem ).size();

            // get the number of SPGs present on the current B-spline cell
            uint tNumSPGs = tBsplineMeshInfo->mSpgIndicesInBsplineCells( iBspElem ).size();

            // for every IP cell in current B-spline element ...
            for ( uint iIpCell = 0; iIpCell < tNumIpCellsInBsplineElement; iIpCell++)
            {
                // get current IP cell's index
                moris_index tIpCellIndex = tBsplineMeshInfo->mExtractionCellsIndicesInBsplineCells( iBspElem )( iIpCell );

                // initialize the map storing the subphases in it
                tBsplineMeshInfo->mExtractionCellToSubPhase( tIpCellIndex ).resize( tNumSPGs );

                // get the SPs on the current element
                moris::Cell< moris_index > const & tSubphaseIndicesOnLagElem = mCutIgMesh->get_parent_cell_subphases( tIpCellIndex );

                // populate the SP to IP-cell map
                for ( uint iSpOnIpCell = 0; iSpOnIpCell < tSubphaseIndicesOnLagElem.size(); iSpOnIpCell++ )
                {
                    // add SP to the map
                    tSpToIpElem[tSubphaseIndicesOnLagElem( iSpOnIpCell )] = tIpCellIndex;
                }
            }

            // go through SPGs on B-spline cell
            for ( uint iSPG = 0; iSPG < tNumSPGs; iSPG++ )
            {
                // get the index of the SPG
                moris_index tSpgIndex = tBsplineMeshInfo->mSpgIndicesInBsplineCells( iBspElem )( iSPG );

                // get the list of SPs on the current SPG
                const moris::Cell< moris_index > & tSPsInSPG = tBsplineMeshInfo->mSubphaseGroups( tSpgIndex )->get_SP_indices_in_group(); 

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
            } // end: loop over all SPGs on B-spline element
        } // end: loop over all B-spline elements
    } // end: function

    // ----------------------------------------------------------------------------------

    void
    Enrichment::construct_enriched_interpolation_vertices_and_cells()
    {
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

        // iterate through subphases and construct an interpolation cell in the interpolation mesh for each one
        for ( moris::uint iSP = 0; iSP < mCutIgMesh->get_num_subphases(); iSP++ )
        {
            // information about this cell
            moris::mtk::Cell* tParentCell = mCutIgMesh->get_subphase_parent_cell( iSP );

            // ID of the owning processor ("Owner")
            moris_id tOwner = tParentCell->get_owner();

            // vertices of cell
            moris::Cell< mtk::Vertex* > tVertices = tParentCell->get_vertex_pointers();

            // bulk phase
            moris_index tBulkPhase = mCutIgMesh->get_subphase_bulk_phase( (moris_index)iSP );

            // loop over the B-spline meshes (i.e. discretization meshes)
            // this is done as the T-matrices at the enriched IP vertices need to be constructed wrt each B-spline mesh
            for ( moris::uint iBspMesh = 0; iBspMesh < mMeshIndices.numel(); iBspMesh++ )
            {
                // Mesh Index
                moris_index tMeshIndex = mMeshIndices( iBspMesh );

                // get the T-matrix of the current IP node wrt the current B-spline mesh
                moris::Cell< mtk::Vertex_Interpolation* > tVertexInterpolations = this->get_vertex_interpolations( *tParentCell, tMeshIndex );

                // get list of non-enriched basis indices interpolating into the current subphase
                Cell< moris_index > const & tBasisInCell = mEnrichmentData( tMeshIndex ).mSubphaseBGBasisIndices( (moris_index)iSP );

                // construct a map between non-enriched BF index and index relative to the subphase cluster
                std::unordered_map< moris_id, moris_id > tCellBasisMap = construct_subphase_basis_to_basis_map( tBasisInCell );

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
                    // Note: the Interpolation_Vertex_Unzipped (IVU) carries a list of vertex enrichments, each VE corresponds to one mesh index
                    // note though, that the IVU is still created for every subphase (i.e. for every material sub-domain within the IP element)
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
                        tEnrInterpMesh->mEnrichedInterpVerts( tVertexIndexInIp )->add_vertex_interpolation( 
                            tMeshIndex, tEnrInterpMesh->get_vertex_enrichment( tMeshIndex, tVertEnrichIndex ) );
                    }
                }

                // create the unzipped interpolation cell on first go
                /* Note: the Interpolation_Cell_Unzipped carries a list of Interpolation_Vertex_Unzipped (IVU) which themselves get updated for every DMI
                 * Hence, the Interpolation_Cell_Unzipped can be left alone after initial creation.
                 * Access to the right IVUs is given once they're all constructed (see code section with double for-loop just below) */
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
                    tEnrInterpMesh->mBaseCelltoEnrichedCell( tParentCell->get_index() ).push_back( 
                        tEnrInterpMesh->mEnrichedInterpCells( tCellIndex ) );

                    // increment the cell index/id
                    tCellIndex++;
                }
            } // end: loop over the Bspline meshes / discretization mesh indices
        } // end: loop over subphases

        // resize out aura cells
        tEnrInterpCellToVertex.resize( tCellIndex, tEnrInterpMesh->mNumVertsPerInterpCell );
        tEnrInterpMesh->mEnrichedInterpCells.resize( tCellIndex );

        // with the cell to vertex data fully setup, add the vertex pointers to the cell
        // for every (unzipped) IP cell get its unzipped vertices
        for ( moris::uint iIpCell = 0; iIpCell < tEnrInterpCellToVertex.n_rows(); iIpCell++ )
        {
            // initialize list of unzipped vertices on cell
            moris::Cell< Interpolation_Vertex_Unzipped* > tVertices( tEnrInterpCellToVertex.n_cols() );

            // iterate through and get unzipped vertices on cell
            for ( moris::uint iVertex = 0; iVertex < tEnrInterpCellToVertex.n_cols(); iVertex++ )
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
    Enrichment::construct_enriched_interpolation_vertices_and_cells_new()
    {
        // get the enriched interpolation mesh pointer, this one is constructed here
        Enriched_Interpolation_Mesh* tEnrInterpMesh = mXTKModelPtr->mEnrichedInterpMesh( 0 );

        // geometry and interpolation order, limited to a single interpolation order mesh
        mtk::Cell const & tFirstCell = mBackgroundMeshPtr->get_mtk_cell( 0 );

        // set the interpolation mesh cell info (i.e. let the enriched IP mesh know what element type it uses)
        mtk::Cell_Info_Factory tFactory;
        tEnrInterpMesh->mCellInfo = tFactory.create_cell_info_sp( tFirstCell.get_geometry_type(), tFirstCell.get_interpolation_order() );

        // maximum mesh index
        moris_index tMaxMeshIndex = mMeshIndices.max();

        // Enriched Interpolation Cell Index to Vertex Index map
        Matrix< IndexMat > tEnrInterpCellToVertex( tEnrInterpMesh->get_num_elements(), tEnrInterpMesh->mNumVertsPerInterpCell );

        // allocate vertex indices and ids
        // NOTE: THESE ARE NOT PARALLEL IDS
        moris_index tVertId      = 1;
        uint        tVertexCount = 0;

        // loop over the B-spline meshes (i.e. discretization mesh indices)
        // this is done as the T-matrices at the enriched IP vertices need to be constructed wrt each B-spline mesh
        for ( moris::uint iBspMesh = 0; iBspMesh < mMeshIndices.numel(); iBspMesh++ )
        {
            // Mesh Index
            moris_index tMeshIndex = mMeshIndices( iBspMesh );

            // iterate through all IP cells (iterator iIpCell is equal to IP cell's index)
            for ( moris::uint iIpCell = 0; iIpCell < mBackgroundMeshPtr->get_num_elems(); iIpCell++ )
            {    
                // get the MTK cell 
                mtk::Cell* tIpCell = &( mCutIgMesh->get_mtk_cell( (moris_index) iIpCell ) );

                // ID of the owning processor ("Owner")
                moris_id tOwner = tIpCell->get_owner();

                // get the number of vertices per IP cell
                uint tNumVertices = tIpCell->get_number_of_vertices();

                // vertices of IP cell
                moris::Cell< mtk::Vertex* > tVertices = tIpCell->get_vertex_pointers();

                // get the T-matrix of the current IP node wrt the current B-spline mesh
                moris::Cell< mtk::Vertex_Interpolation* > tVertexInterpolations = this->get_vertex_interpolations( *tIpCell, tMeshIndex );

                // get the number of SPGs that are present in the associated B-spline element
                uint tNumSpgsAssociatedWithIpCell = mBsplineMeshInfos( iBspMesh )->get_num_SPGs_associated_with_extraction_cell( iIpCell );

                // get max number of enr IP cells associated with the current IP cell
                uint tMaxNumUnzippings = mNumUnzippingsOnIpCell( iIpCell );

                // iterate through the enriched IP cells that are to be constructed on the current base IP cell
                for ( moris::uint iEnrIpCell = 0; iEnrIpCell < tMaxNumUnzippings; iEnrIpCell++ )
                {
                    // get the enr. IP cell's index
                    moris_index tCellIndex = mEnrIpCellIndices( iIpCell )( iEnrIpCell );

                    // get the T-matrix of the current IP node wrt the current B-spline mesh
                    moris::Cell< mtk::Vertex_Interpolation* > tVertexInterpolations = this->get_vertex_interpolations( *tIpCell, tMeshIndex );

                    // sanity check that number of T-matrices matches number of nodes
                    MORIS_ASSERT( tVertexInterpolations.size() == tNumVertices, 
                        "Enrichment::construct_enriched_interpolation_vertices_and_cells_new() -  number of T-matrices on cell and number of Nodes don't match up." );

                    // initialize variables that are only used when constructiong valid enr. IP cell
                    std::unordered_map< moris_id, moris_id > tCellBasisMap;
                    moris_index tSpgIndex = -1;
                    moris_index tVertEnrichIndex = -1;

                    // check if current Enr. IP cell is invalid wrt to current B-spline mesh
                    // if so, get the information related to the active SPG
                    if( iEnrIpCell < tNumSpgsAssociatedWithIpCell )
                    {
                        // get the index of the current SPG, set to -1 if for Enr. IP cell is invalid wrt to current B-spline mesh
                        tSpgIndex = mBsplineMeshInfos( iBspMesh )->get_SPG_indices_associated_with_extraction_cell( iIpCell )( iEnrIpCell );

                        // get list of non-enriched basis indices interpolating into the current subphase
                        Cell< moris_index > const & tBasisInSpg = mEnrichmentData( tMeshIndex ).mSubphaseGroupBGBasisIndices( tSpgIndex );

                        // construct a map between non-enriched BF index and index relative to the subphase cluster
                        // NOTE: the function is indifferent to whether it's operating on SPs or SPGs
                        tCellBasisMap = construct_subphase_basis_to_basis_map( tBasisInSpg );
                    }

                    // construction of unzipped enriched IP vertices
                    // i.e. only construct T-matrices if underlying enriched IP cell is valid
                    for ( uint iIpCellVertex = 0; iIpCellVertex < tNumVertices; iIpCellVertex++ )
                    {
                        // check if current Enr. IP cell is invalid wrt to current B-spline mesh
                        // if so, compute and add T-matrix 
                        if( iEnrIpCell < tNumSpgsAssociatedWithIpCell )
                        {
                            // initialize vertex enrichment object
                            Vertex_Enrichment tVertEnrichment;

                            // compute T-matrix (aka. Vertex Enrichment)                       
                            this->construct_enriched_vertex_interpolation(
                                    tMeshIndex,
                                    tVertexInterpolations( iIpCellVertex ),
                                    mEnrichmentData( tMeshIndex ).mSubphaseGroupBGBasisEnrLev( tSpgIndex ),
                                    tCellBasisMap,
                                    tVertEnrichment );

                            // add vertex enrichment to enriched interpolation mesh
                            bool tNewVertFlag = false;
                            tVertEnrichIndex = tEnrInterpMesh->add_vertex_enrichment(
                                    tMeshIndex,
                                    tVertices( iIpCellVertex ), // feed non-enriched Vertex 
                                    tVertEnrichment,
                                    tNewVertFlag );
                        }
                        else // Enr. IP cell is invalid
                        {
                            // add empty T-matrix to enriched interpolation mesh
                            bool tNewVertFlag = false;
                            tVertEnrichIndex = tEnrInterpMesh->add_empty_vertex_enrichment(
                                    tMeshIndex,
                                    tVertices( iIpCellVertex ),
                                    tNewVertFlag );
                        }

                        // create this vertex on the first go around
                        // Note: the Interpolation_Vertex_Unzipped (IVU) carries a list of vertex enrichments, each VE corresponds to one mesh index
                        // note though, that the IVU is still created for all subphase groups associdated with every IP cell
                        if ( iBspMesh == 0 )
                        {
                            // Create interpolation vertex with only the first Vertex enrichment
                            tEnrInterpMesh->mEnrichedInterpVerts( tVertexCount ) =
                                    new Interpolation_Vertex_Unzipped(
                                            tVertices( iIpCellVertex ),              // non-enriched IP vertices
                                            tVertId,                                 // current IP vertex' ID  
                                            tVertexCount,                            // index for enriched IP vertex
                                            tVertices( iIpCellVertex )->get_owner(), // owning proc                // NOTE: the two inputs below are only for initialization, not final information
                                            tMeshIndex,                                                            // current DMI
                                            tEnrInterpMesh->get_vertex_enrichment( tMeshIndex, tVertEnrichIndex ), // T-matrix
                                            tMaxMeshIndex );                         // maximum DMI 

                            // store enriched vertex's index for given parent cell index and element local node index 
                            tEnrInterpCellToVertex( tCellIndex, iIpCellVertex ) = tVertexCount;

                            // update vertex ID to use for nex unzipped vertex
                            tVertId++;
                            
                            // track number of unzipped vertices that have been created
                            tVertexCount++;
                        }
                        else
                        {
                            // the unzipped interpolation vertex' index
                            moris_index tVertexIndexInIp = tEnrInterpCellToVertex( tCellIndex, iIpCellVertex );

                            // add the vertex interpolation for new mesh index
                            tEnrInterpMesh->mEnrichedInterpVerts( tVertexIndexInIp )->add_vertex_interpolation( 
                                tMeshIndex, tEnrInterpMesh->get_vertex_enrichment( tMeshIndex, tVertEnrichIndex ) );
                        }
                    } // end: loop over IP vertices

                    // create the unzipped interpolation cell on first go
                    /* Note: the Interpolation_Cell_Unzipped carries a list of Interpolation_Vertex_Unzipped (IVU) which themselves get updated for every DMI
                     * Hence, the Interpolation_Cell_Unzipped can be left alone after initial creation.
                     * Access to the right IVUs is given once they're all constructed (see code section with double for-loop just below) */
                    if ( iBspMesh == 0 )
                    {
                        // FIXME: use a communicated proc global ID once everything is parallel consistent
                        moris_index tCellID = tCellIndex;
                        
                        // create new enriched interpolation cell and put it in list associating it with the underlying parent IP cell
                        tEnrInterpMesh->mEnrichedInterpCells( tCellIndex ) =
                                new Interpolation_Cell_Unzipped(
                                        tIpCell,
                                        (moris_index) iEnrIpCell,
                                        tCellID,
                                        tCellIndex,
                                        tOwner,
                                        tEnrInterpMesh->mCellInfo,
                                        true );

                        // add enriched interpolation cell to base cell to enriched cell data
                        tEnrInterpMesh->mBaseCelltoEnrichedCell( tIpCell->get_index() ).push_back( 
                            tEnrInterpMesh->mEnrichedInterpCells( tCellIndex ) );
                    }

                } // end: loop over SPGs associated with IP cell
            } // end: loop over IP cells on mesh
        } // end: loop over DMIs

        // resize out aura cells
        tEnrInterpCellToVertex.resize( mNumEnrIpCells, tEnrInterpMesh->mNumVertsPerInterpCell );
        tEnrInterpMesh->mEnrichedInterpCells.resize( mNumEnrIpCells );

        // with the cell to vertex data fully setup, add the vertex pointers to the cell
        // for every (unzipped) IP cell get its unzipped vertices
        for ( moris::uint iEnrIpCell = 0; iEnrIpCell < tEnrInterpCellToVertex.n_rows(); iEnrIpCell++ )
        {
            // initialize list of unzipped vertices on cell
            moris::Cell< Interpolation_Vertex_Unzipped* > tVertices( tEnrInterpCellToVertex.n_cols() );

            // iterate through and get unzipped vertices on cell
            for ( moris::uint iVertex = 0; iVertex < tEnrInterpCellToVertex.n_cols(); iVertex++ )
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

    void
    Enrichment::construct_enriched_vertex_interpolation(
            moris_index const &                       aEnrichmentDataIndex,
            mtk::Vertex_Interpolation*                aBaseVertexInterp,
            Cell< moris_index > const &               aSubPhaseBasisEnrLev,
            std::unordered_map< moris_id, moris_id >& aMapBasisIndexToLocInSubPhase,
            Vertex_Enrichment&                        aVertexEnrichment )
    {
        // allocate a new vertex enrichment
        aVertexEnrichment = Vertex_Enrichment();

        // a nullptr here would indicate an aura node without a t-matrix
        if ( aBaseVertexInterp != nullptr )
        {
            // iterate through basis in the base vertex interpolation
            moris::uint tNumCoeffs = aBaseVertexInterp->get_number_of_coefficients();

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
            for ( moris::uint iBC = 0; iBC < tNumCoeffs; iBC++ )
            {
                // coefficient of the base index
                moris_index tCoeffIndex = tBaseVertCoeffInds( iBC );

                // find the coefficients index within this subphase cluster
                auto tIter = aMapBasisIndexToLocInSubPhase.find( tCoeffIndex );

                MORIS_ASSERT( tIter != aMapBasisIndexToLocInSubPhase.end(), "Basis not found in vertex map" );

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
            std::unordered_map< moris::moris_index, moris::moris_index >& tVertEnrichMap = aVertexEnrichment.get_basis_map();

            // store in the vertex enrichment the list of enriched BF coefficients associated with it (in the form of an index map)
            for ( moris::uint iBC = 0; iBC < tEnrichCoeffInds.numel(); iBC++ )
            {
                moris::moris_index tBasisIndex = tEnrichCoeffInds( iBC );
                tVertEnrichMap[tBasisIndex] = iBC;
            }

            // feed the generated data about the enriched vertex interpolation to the vertex enrichment object
            aVertexEnrichment.add_basis_information( tEnrichCoeffInds, tEnrichCoeffIds );
            aVertexEnrichment.add_basis_owners( tEnrichCoeffInds, tBaseVertOwners );
            aVertexEnrichment.add_basis_weights( tEnrichCoeffInds, *tBaseVertWeights );
            aVertexEnrichment.add_base_vertex_interpolation( aBaseVertexInterp );
        }
    }

    //-------------------------------------------------------------------------------------

    std::unordered_map< moris_id, moris_id >
    Enrichment::construct_subphase_basis_to_basis_map( Cell< moris_id > const & aSubPhaseBasisIndex )
    {
        // get number of BFs interpolating into subphase
        uint tNumBasisOfSubphase = aSubPhaseBasisIndex.size();

        // initialize output map
        std::unordered_map< moris_id, moris_id > tSubphaseBasisMap;

        // for each basis interpolating into the current Subphase, ...
        for ( moris::uint iB = 0; iB < tNumBasisOfSubphase; iB++ )
        {
            // ... get its Index and relate it to the position in the list of basis for that particular subphase
            tSubphaseBasisMap[aSubPhaseBasisIndex( iB )] = iB;
        }

        // return the map
        return tSubphaseBasisMap;
    }

    //-------------------------------------------------------------------------------------

    moris::Cell< mtk::Vertex_Interpolation* >
    Enrichment::get_vertex_interpolations(
            moris::mtk::Cell& aParentCell,
            const uint        aMeshIndex ) const
    {
        uint tNumVerts = aParentCell.get_number_of_vertices();

        moris::Cell< mtk::Vertex* > tVertexPointers = aParentCell.get_vertex_pointers();

        moris::Cell< mtk::Vertex_Interpolation* > tVertexInterp( tNumVerts );

        for ( moris::uint i = 0; i < tNumVerts; i++ )
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
        moris::uint tNumBasisTypes = mEnrichmentData.size();

        // declare  enrichment fields
        Cell< std::string > tEnrichmentFields;

        for ( moris::uint iBT = 0; iBT < tNumBasisTypes; iBT++ )
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

        // Add local floodfill field to the output mesh
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
        moris::Cell< moris::moris_id > aGatheredInfo;
        moris::Cell< moris::moris_id > tFirstId( 1 );
        moris::Cell< moris::moris_id > tNumIdsRequested( 1 );

        tNumIdsRequested( 0 ) = (moris::moris_id)aNumIdsToAllocate;

        moris::gather( tNumIdsRequested, aGatheredInfo );

        moris::Cell< moris::moris_id > tProcFirstID( tProcSize );

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

        for ( moris::uint i = 0; i < tNumBasis; i++ )
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
}// namespace xtk
