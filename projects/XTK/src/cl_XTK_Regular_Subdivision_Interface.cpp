/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Regular_Subdivision_Interface.cpp
 *
 */

#include "cl_XTK_Decomposition_Algorithm.hpp"
#include "cl_XTK_Regular_Subdivision_Interface.hpp"
#include "cl_XTK_Integration_Mesh_Generator.hpp"
#include "cl_XTK_Cut_Integration_Mesh.hpp"
#include "cl_MTK_Cell_Info_Factory.hpp"
#include "cl_MTK_Cell_Info.hpp"
#include "cl_XTK_Octree_Interface.hpp"

namespace moris::xtk
{
    Regular_Subdivision_Interface::Regular_Subdivision_Interface( Parameter_List& aParameterList, mtk::CellTopology aCellTopology )
    {
        // get number of spatial dimensions and decide on subdivision template
        if ( aCellTopology == mtk::CellTopology::QUAD4 )
        {
            mRegularSubdivisionTemplate = std::make_shared< Regular_Subdivision_4_TRIS >();
        }
        else if ( aCellTopology == mtk::CellTopology::HEX8 )
        {
            mRegularSubdivisionTemplate = std::make_shared< Regular_Subdivision_24_TETS >();
        }
        else
        {
            MORIS_ERROR( false, "Regular_Subdivision_Interface::Regular_Subdivision_Interface() - subdivision routine for BG element type unknown" );
        }
    }

    //--------------------------------------------------------------------------------------------------

    bool
    Regular_Subdivision_Interface::is_eligible( std::pair< mtk::Cell*, Vector< Decomposition_Algorithm_Type > >& aElementContext,
            Cut_Integration_Mesh*                                                                                aCutIntegrationMesh,
            Integration_Mesh_Generator*                                                                          aMeshGenerator ) const
    {
        // Get the existing decomposition types
        Vector< Decomposition_Algorithm_Type > tDecompTypes = aElementContext.second;

        // Check if this element is doing Delaunay
        bool tIsDelaunay = std::find( tDecompTypes.begin(), tDecompTypes.end(), Decomposition_Algorithm_Type::DELAUNAY ) != tDecompTypes.end();

        // Decompose if it the element has not been Delaunay triangulated and is intersected
        // return not tIsDelaunay and aMeshGenerator->is_intersected( aElementContext.first );
        return not tIsDelaunay;
    }

    //--------------------------------------------------------------------------------------------------

    Vector< moris_index >
    Regular_Subdivision_Interface::get_decomposed_cell_indices()
    {
        return mMeshGenerationData->mAllIntersectedBgCellInds;
    }

    //--------------------------------------------------------------------------------------------------

    bool
    Regular_Subdivision_Interface::has_geometric_independent_vertices() const
    {
        return true;
    }

    //--------------------------------------------------------------------------------------------------

    void
    Regular_Subdivision_Interface::perform_impl_vertex_requests(
            Integration_Mesh_Generation_Data* aMeshGenerationData,
            Decomposition_Data*               aDecompositionData,
            Cut_Integration_Mesh*             aCutIntegrationMesh,
            moris::mtk::Mesh*                 aBackgroundMesh,
            Integration_Mesh_Generator*       aMeshGenerator )
    {
        // populate member data
        mMeshGenerationData = aMeshGenerationData;
        mDecompositionData  = aDecompositionData;
        mCutIntegrationMesh = aCutIntegrationMesh;
        mBackgroundMesh     = aBackgroundMesh;
        mGenerator          = aMeshGenerator;

        // get number of spatial dimensions and decide on subdivision template
        if ( aMeshGenerator->get_spatial_dim() == 2 )
        {
            mRegularSubdivisionTemplate = std::make_shared< Regular_Subdivision_4_TRIS >();
        }
        else if ( aMeshGenerator->get_spatial_dim() == 3 )
        {
            mRegularSubdivisionTemplate = std::make_shared< Regular_Subdivision_24_TETS >();
        }
        else
        {
            MORIS_ERROR( false, "Regular_Subdivision_Interface() - spatial dimension not 2 or 3" );
        }

        aDecompositionData->mDecompId = this->get_signature();

        aDecompositionData->mHasSecondaryIdentifier = false;

        moris_index tNumChildMeshes = mCutIntegrationMesh->get_num_ig_cell_groups();

        // allocate data in decomposition data
        aDecompositionData->tCMNewNodeLoc        = Vector< Vector< moris_index > >( tNumChildMeshes );
        aDecompositionData->tCMNewNodeParamCoord = Vector< Vector< Matrix< DDRMat > > >(
                tNumChildMeshes, Vector< Matrix< DDRMat > >( this->get_num_new_nodes(), Matrix< DDRMat >( 1, this->get_parametric_dimension() ) ) );

        // allocate a data struct to pass into functions - stores the matrix but doesnt require the method to store data
        Regular_Subdivision_Interface_Data tRegSubInterfaceData;

        // get the new nodes on each entity rank
        tRegSubInterfaceData.mNewNodesOnEdges = this->get_new_node_on_parent_edge();
        tRegSubInterfaceData.mNewNodesOnFaces = this->get_new_node_on_parent_face();
        tRegSubInterfaceData.mNewNodesOnCells = this->get_new_node_in_parent_cell();

        // get the new node edge,face,cell ordinal
        tRegSubInterfaceData.mNewNodesOnEdgesOrd = this->get_new_node_on_parent_edge_edge_ordinal();
        tRegSubInterfaceData.mNewNodesOnFacesOrd = this->get_new_node_on_parent_face_face_ordinal();
        tRegSubInterfaceData.mNewNodesOnCellsOrd = this->get_new_node_on_parent_cell_cell_ordinal();
        tRegSubInterfaceData.mVertexAncestry     = mRegularSubdivisionTemplate->get_vertex_ancestry();

        tRegSubInterfaceData.mNewNodeXi = this->get_new_vertex_parametric_coordinates_wrt_parent();

        // iterate through all intersected background cells and make vertex requests
        for ( auto& iCell : aMeshGenerationData->mAllIntersectedBgCellInds )
        {
            std::shared_ptr< Child_Mesh_Experimental > tChildMesh = aCutIntegrationMesh->get_child_mesh( iCell );

            // make new vertex requests
            this->make_new_vertex_requests( tChildMesh.get(), aCutIntegrationMesh, this, &tRegSubInterfaceData, aBackgroundMesh, aDecompositionData );
        }

        // handle the requests
        aDecompositionData->tSecondaryIdentifiers = Vector< moris_index >( aDecompositionData->tNewNodeParentIndex.size(), MORIS_INDEX_MAX );
    }

    //--------------------------------------------------------------------------------------------------

    void
    Regular_Subdivision_Interface::perform_impl_generate_mesh(
            Integration_Mesh_Generation_Data* aMeshGenerationData,
            Decomposition_Data*               aDecompositionData,
            Cut_Integration_Mesh*             aCutIntegrationMesh,
            moris::mtk::Mesh*                 aBackgroundMesh,
            Integration_Mesh_Generator*       aMeshGenerator )
    {
        // new cell info
        moris::mtk::Cell_Info_Factory            tFactory;
        std::shared_ptr< moris::mtk::Cell_Info > tIgCellInfo = tFactory.create_cell_info_sp( this->get_ig_cell_topology() );

        // this is going to be the cell info pointer for all new cells

        // number of vertices per cell
        moris::uint tVerticesPerCell = tIgCellInfo->get_num_verts();

        // allocate data in the new ig cell data
        mNewCellToVertexConnectivity = Vector< Vector< moris::moris_index > >( mNumTotalCells );
        mNewCellChildMeshIndex       = Vector< moris::moris_index >( mNumTotalCells );
        mNewCellCellIndexToReplace   = Vector< moris::moris_index >( mNumTotalCells, MORIS_INDEX_MAX );
        mNewCellCellInfo             = Vector< std::shared_ptr< moris::mtk::Cell_Info > >( mNumTotalCells, tIgCellInfo );

        // get the cell to vertex template
        Vector< Vector< moris::moris_index > > tIgCellToVertexTemplate = this->get_ig_cell_to_vertex_connectivity();

        // populate new cell data
        moris::moris_index tCurrentCellIndex = 0;
        for ( auto& iCM : aMeshGenerationData->mRegularSubdivisionBgCellInds )
        {
            std::shared_ptr< Child_Mesh_Experimental > tChildMesh = aCutIntegrationMesh->get_child_mesh( iCM );

            if ( tChildMesh->mIgCells->mIgCellGroup.size() <= 1 )
            {
                for ( moris::moris_index iNewCell = 0; iNewCell < this->get_num_ig_cells(); iNewCell++ )
                {
                    mNewCellChildMeshIndex( tCurrentCellIndex ) = iCM;

                    for ( moris::uint iV = 0; iV < tVerticesPerCell; iV++ )
                    {
                        moris_index tNewVertexCMOrdinal = tIgCellToVertexTemplate( iNewCell )( iV );
                        MORIS_ERROR( tNewVertexCMOrdinal < (moris::moris_index)tChildMesh->mIgVerts->size(), "Template ordinal out of bounds" );
                        mNewCellToVertexConnectivity( tCurrentCellIndex ).push_back( tChildMesh->mIgVerts->get_vertex( tNewVertexCMOrdinal )->get_index() );
                    }

                    tCurrentCellIndex++;
                }
            }
            // if the child mesh has more than one cell, i.e. octree/quadtree refinement has been performed beforehand, the template needs to be applied repeatedly
            else
            {
                std::shared_ptr< Generated_Regular_Subdivision_Template > tGeneratedTemplate = mGeneratedTemplate( tChildMesh->mIgCells->mIgCellGroup.size() );

                // change the replace data
                for ( moris::moris_index iReplace = 0; iReplace < (moris::moris_index)tChildMesh->mIgCells->mIgCellGroup.size(); iReplace++ )
                {
                    mNewCellCellIndexToReplace( tCurrentCellIndex + iReplace ) = tChildMesh->mIgCells->mIgCellGroup( iReplace )->get_index();
                }

                for ( moris::moris_index iNewCell = 0; iNewCell < tGeneratedTemplate->get_num_ig_cells(); iNewCell++ )
                {
                    mNewCellChildMeshIndex( tCurrentCellIndex ) = iCM;

                    for ( moris::uint iV = 0; iV < tVerticesPerCell; iV++ )
                    {
                        moris_index tNewVertexCMOrdinal = tGeneratedTemplate->mIgCellToVertOrd( iNewCell )( iV );
                        MORIS_ERROR( tNewVertexCMOrdinal < (moris::moris_index)tChildMesh->mIgVerts->size(), "Template ordinal out of bounds" );
                        mNewCellToVertexConnectivity( tCurrentCellIndex )( iV ) = tChildMesh->mIgVerts->get_vertex( tNewVertexCMOrdinal )->get_index();
                    }
                    tCurrentCellIndex++;
                }
            }
        }    // end for: each intersected background element
    }    // end function: Regular_Subdivision_Interface::perform_impl_generate_mesh()

    //--------------------------------------------------------------------------------------------------

    enum Decomposition_Algorithm_Type
    Regular_Subdivision_Interface::get_algorithm_type() const
    {
        return Decomposition_Algorithm_Type::REGULAR_TEMPLATE_NONCONFORMING;
    }

    //--------------------------------------------------------------------------------------------------

    void
    Regular_Subdivision_Interface::make_new_vertex_requests(
            Child_Mesh_Experimental*            aChildMesh,
            Cut_Integration_Mesh*               aCutIntegrationMesh,
            Regular_Subdivision_Interface*      aRegularSubdivisionInterface,
            Regular_Subdivision_Interface_Data* aRegularSubdivisionInterfaceData,
            moris::mtk::Mesh*                   aBackgroundMesh,
            Decomposition_Data*                 aDecompositionData )
    {

        // standard case where ancestry is directly known
        if ( aChildMesh->mIgCells->mIgCellGroup.size() == 0 || aChildMesh->mIgCells->mIgCellGroup.size() == 1 )
        {
            this->make_new_vertex_requests_trivial(
                    aChildMesh,
                    aCutIntegrationMesh,
                    aRegularSubdivisionInterface,
                    aRegularSubdivisionInterfaceData,
                    aBackgroundMesh,
                    aDecompositionData );

            mNumNewCells   = mNumNewCells + mRegularSubdivisionTemplate->get_num_ig_cells();
            mNumTotalCells = mNumNewCells;
        }
        else
        {
            this->make_new_vertex_requests_octree(
                    aChildMesh,
                    aCutIntegrationMesh,
                    aRegularSubdivisionInterface,
                    aRegularSubdivisionInterfaceData,
                    aBackgroundMesh,
                    aDecompositionData );

            mNumNewCells =
                    mNumNewCells + mGeneratedTemplate( aChildMesh->mIgCells->mIgCellGroup.size() )->get_num_ig_cells() - aChildMesh->mIgCells->mIgCellGroup.size();

            mNumTotalCells =
                    mNumTotalCells + mGeneratedTemplate( aChildMesh->mIgCells->mIgCellGroup.size() )->get_num_ig_cells();
        }
    }

    //--------------------------------------------------------------------------------------------------

    void
    Regular_Subdivision_Interface::make_new_vertex_requests_trivial(
            Child_Mesh_Experimental*            aChildMesh,
            Cut_Integration_Mesh*               aCutIntegrationMesh,
            Regular_Subdivision_Interface*      aRegularSubdivisionInterface,
            Regular_Subdivision_Interface_Data* aRegularSubdivisionInterfaceData,
            moris::mtk::Mesh*                   aBackgroundMesh,
            Decomposition_Data*                 aDecompositionData )
    {
        // collect information from parent BG cell
        moris::mtk::Cell*             tParentCell     = aChildMesh->get_parent_cell();
        moris::mtk::Cell_Info const * tParentCellInfo = tParentCell->get_cell_info();
        Cell_Connectivity             tCellConn       = aCutIntegrationMesh->get_background_cell_connectivity( tParentCell->get_index() );
        moris_index                   tChildMeshIndex = aChildMesh->get_child_mesh_index();

        // process every node requested on a faces
        for ( moris::uint iVertsOnFaces = 0; iVertsOnFaces < aRegularSubdivisionInterfaceData->mNewNodesOnFaces.numel(); iVertsOnFaces++ )
        {
            // gather template data
            moris_index       tNewNodeFaceOrdinal = aRegularSubdivisionInterfaceData->mNewNodesOnFacesOrd( iVertsOnFaces );
            moris_index       tNewNodeTemplateOrd = aRegularSubdivisionInterfaceData->mNewNodesOnFaces( iVertsOnFaces );
            Matrix< DDRMat >& tNewNodeXi          = aRegularSubdivisionInterfaceData->mNewNodeXi( tNewNodeTemplateOrd );

            // check whether a request on this entity (face) already exists
            moris_index tRequestLoc    = MORIS_INDEX_MAX;
            moris_id    tBgFacetIndex  = tCellConn.mCellFacesInds( tNewNodeFaceOrdinal );
            bool        tRequestExists = aDecompositionData->request_exists( tBgFacetIndex, mtk::EntityRank::FACE, tRequestLoc );

            // create a new request if it doesn't exist already
            if ( !tRequestExists )
            {
                moris_index tOwner = aBackgroundMesh->get_entity_owner( tBgFacetIndex, mtk::EntityRank::FACE );

                // evaluate the shape functions at this point relative to the background cell
                tParentCellInfo->eval_N( tNewNodeXi, aRegularSubdivisionInterfaceData->mNXi );

                // create a new request
                std::shared_ptr< Matrix< DDRMat > > tNewNodeXiPtr              = std::make_shared< Matrix< DDRMat > >( tNewNodeXi );
                Matrix< DDRMat >                    tNewNodePhysCoords         = aRegularSubdivisionInterfaceData->mNXi * tParentCell->get_vertex_coords();
                moris_index                         tNewNodeIndexInSubdivision = aDecompositionData->register_new_request(
                        tBgFacetIndex,
                        tOwner,
                        mtk::EntityRank::FACE,
                        tNewNodePhysCoords,
                        tParentCell,
                        tNewNodeXiPtr );

                // store the new node in the decomposition data
                aDecompositionData->tCMNewNodeParamCoord( tChildMeshIndex )( tNewNodeTemplateOrd ) = tNewNodeXi;
                aDecompositionData->tCMNewNodeLoc( tChildMeshIndex ).push_back( tNewNodeIndexInSubdivision );
            }

            // if the request already exists, reference the request and add the location of the new node wrt. the current BG cell to the child mesh
            else
            {
                aDecompositionData->tCMNewNodeLoc( tChildMeshIndex ).push_back( tRequestLoc );
                aDecompositionData->tCMNewNodeParamCoord( tChildMeshIndex )( tNewNodeTemplateOrd ) = tNewNodeXi;
            }

        }    // end for: iVertsOnFaces

        // process nodes requested in BG cell interior
        Matrix< IndexMat > tElementIndices( 0, 0 );
        for ( moris::uint iVertsInCell = 0; iVertsInCell < aRegularSubdivisionInterfaceData->mNewNodesOnCells.numel(); iVertsInCell++ )
        {
            moris_index tRequestLoc    = MORIS_INDEX_MAX;
            bool        tRequestExists = aDecompositionData->request_exists( tParentCell->get_index(), mtk::EntityRank::ELEMENT, tRequestLoc );

            moris_index tNewNodeTemplateOrd = aRegularSubdivisionInterfaceData->mNewNodesOnCells( iVertsInCell );
            if ( !tRequestExists )
            {
                moris_index tOwner = aBackgroundMesh->get_entity_owner( tParentCell->get_index(), mtk::EntityRank::ELEMENT );

                // evaluate the shape functions at this point relative to the background cell
                tParentCellInfo->eval_N( aRegularSubdivisionInterfaceData->mNewNodeXi( tNewNodeTemplateOrd ), aRegularSubdivisionInterfaceData->mNXi );

                std::shared_ptr< Matrix< DDRMat > > tNewNodeXi =
                        std::make_shared< Matrix< DDRMat > >( aRegularSubdivisionInterfaceData->mNewNodeXi( tNewNodeTemplateOrd ) );

                Matrix< DDRMat > tNewNodeCoordinates = aRegularSubdivisionInterfaceData->mNXi * tParentCell->get_vertex_coords();

                moris_index tNewNodeIndexInSubdivision = aDecompositionData->register_new_request(
                        tParentCell->get_index(),
                        tOwner,
                        mtk::EntityRank::ELEMENT,
                        tNewNodeCoordinates,
                        tParentCell,
                        tNewNodeXi );

                aDecompositionData->tCMNewNodeLoc( aChildMesh->get_child_mesh_index() ).push_back( tNewNodeIndexInSubdivision );

                aDecompositionData->tCMNewNodeParamCoord( aChildMesh->get_child_mesh_index() )( tNewNodeTemplateOrd ) =
                        aRegularSubdivisionInterfaceData->mNewNodeXi( tNewNodeTemplateOrd );
            }

            else
            {
                aDecompositionData->tCMNewNodeLoc( aChildMesh->get_child_mesh_index() ).push_back( tRequestLoc );

                aDecompositionData->tCMNewNodeParamCoord( aChildMesh->get_child_mesh_index() )( tNewNodeTemplateOrd ) =
                        aRegularSubdivisionInterfaceData->mNewNodeXi( tNewNodeTemplateOrd );
            }
        }    // end for: iVertsInCell
    }    // end function: Regular_Subdivision_Interface::make_new_vertex_requests()

    //--------------------------------------------------------------------------------------------------

    void
    Regular_Subdivision_Interface::make_new_vertex_requests_octree(
            Child_Mesh_Experimental*            aChildMesh,
            Cut_Integration_Mesh*               aCutIntegrationMesh,
            Regular_Subdivision_Interface*      aRegularSubdivisionInterface,
            Regular_Subdivision_Interface_Data* aRegularSubdivisionInterfaceData,
            moris::mtk::Mesh*                   aBackgroundMesh,
            Decomposition_Data*                 aDecompositionData )
    {

        aDecompositionData->mHasSecondaryIdentifier = true;

        aDecompositionData->tCMNewNodeParamCoord( aChildMesh->get_child_mesh_index() ).clear();
        aDecompositionData->tCMNewNodeLoc( aChildMesh->get_child_mesh_index() ).clear();

        // generate a template for this child mesh if one hasn't been created yet
        this->generate_new_node_parent_information_ijk_mesh( aRegularSubdivisionInterfaceData, aChildMesh->mIgCells->mIgCellGroup.size(), aChildMesh );

        std::shared_ptr< Generated_Regular_Subdivision_Template > tGeneratedTemplate = mGeneratedTemplate( aChildMesh->mIgCells->mIgCellGroup.size() );

        moris::mtk::Cell* tCell = aChildMesh->get_parent_cell();

        Cell_Connectivity tCellConnectivity = mCutIntegrationMesh->get_background_cell_connectivity( tCell->get_index() );

        Matrix< DDRMat > tBGCellCoords = tCell->get_vertex_coords();

        Vertex_Ancestry* tVertexAncestry = &tGeneratedTemplate->mNewVertexAncestry;

        // iterate through the vertex requests
        for ( moris::moris_index iV = 0; iV < tGeneratedTemplate->mNumNewNodes; iV++ )
        {

            Matrix< DDRMat > tBasisWeights;

            // local coordinate of this vertex wrt the current cell group
            Matrix< DDRMat > const & tVertexLocalCoords = tGeneratedTemplate->mParamCoords( iV );

            // evaluate the basis function
            tCell->get_cell_info()->eval_N( tVertexLocalCoords, tBasisWeights );

            Matrix< DDRMat > tNewCoordinate = tBasisWeights * tBGCellCoords;

            if ( tVertexAncestry->get_vertex_parent_rank( iV ) == mtk::EntityRank::ELEMENT )
            {
                moris_index tNewNodeIndexInSubdivision = aDecompositionData->register_new_request(
                        tCell->get_index(),
                        tGeneratedTemplate->mVertexHash( iV ),
                        tCell->get_owner(),
                        mtk::EntityRank::ELEMENT,
                        tNewCoordinate,
                        tCell,
                        std::make_shared< Matrix< DDRMat > >( tVertexLocalCoords ) );

                aDecompositionData->tCMNewNodeParamCoord( aChildMesh->get_child_mesh_index() ).push_back( tVertexLocalCoords );
                aDecompositionData->tCMNewNodeLoc( aChildMesh->get_child_mesh_index() ).push_back( tNewNodeIndexInSubdivision );
            }

            else
            {
                moris_index tOwningProc = mBackgroundMesh->get_entity_owner(
                        tVertexAncestry->get_vertex_parent_index( iV ),
                        tVertexAncestry->get_vertex_parent_rank( iV ) );

                moris_index tNewNodeIndexInSubdivision = aDecompositionData->register_new_request(
                        tCellConnectivity.get_entity_index( tVertexAncestry->get_vertex_parent_index( iV ), tVertexAncestry->get_vertex_parent_rank( iV ) ),
                        tGeneratedTemplate->mVertexHash( iV ),
                        tOwningProc,
                        tVertexAncestry->get_vertex_parent_rank( iV ),
                        tNewCoordinate,
                        tCell,
                        std::make_shared< Matrix< DDRMat > >( tVertexLocalCoords ) );

                aDecompositionData->tCMNewNodeParamCoord( aChildMesh->get_child_mesh_index() ).push_back( tVertexLocalCoords );
                aDecompositionData->tCMNewNodeLoc( aChildMesh->get_child_mesh_index() ).push_back( tNewNodeIndexInSubdivision );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------

    void
    Regular_Subdivision_Interface::generate_new_node_parent_information_ijk_mesh(
            Regular_Subdivision_Interface_Data* aRegularSubdivisionInterfaceData,
            moris::uint                         aNumIgCells,
            Child_Mesh_Experimental*            aChildMesh )
    {
        if ( mGeneratedTemplate.size() < aNumIgCells )
        {
            mGeneratedTemplate.resize( aNumIgCells + 1, nullptr );
        }

        if ( mGeneratedTemplate( aNumIgCells ) == nullptr )
        {
            mGeneratedTemplate( aNumIgCells ) = std::make_shared< Generated_Regular_Subdivision_Template >();

            // Start by getting the ancestry of each facet in the current child mesh
            std::shared_ptr< Facet_Based_Connectivity > tFaceConnectivity = std::make_shared< Facet_Based_Connectivity >();

            // FIXME: the below method is deprecated, so is the whole octree interface
            MORIS_ERROR( false, "Regular_Subdivision_Interface::generate_new_node_parent_information_ijk_mesh() - method is deprecated" );
            // mGenerator->create_facet_from_element_to_node( aChildMesh->mIgCells->mIgCellGroup, tFaceConnectivity );

            std::shared_ptr< IG_Vertex_Group > tVertexGroup = aChildMesh->mIgVerts;

            Vector< moris::mtk::Cell* > tBGCellForFacet( tFaceConnectivity->mFacetVertices.size(), aChildMesh->get_parent_cell() );

            std::shared_ptr< Facet_Based_Ancestry > tFacetAncestry = std::make_shared< Facet_Based_Ancestry >();
            mGenerator->deduce_facet_ancestry( mCutIntegrationMesh, mBackgroundMesh, tFaceConnectivity, tBGCellForFacet, tFacetAncestry );

            Cell_Connectivity tCellConn                     = mCutIntegrationMesh->get_background_cell_connectivity( aChildMesh->get_parent_cell()->get_index() );
            const moris::uint tNumNewCells                  = mRegularSubdivisionTemplate->get_num_ig_cells();
            mGeneratedTemplate( aNumIgCells )->mNumIgCells  = tNumNewCells * aNumIgCells;
            mGeneratedTemplate( aNumIgCells )->mNumNewNodes = tFaceConnectivity->mFacetVertices.size() + aNumIgCells;

            // local coord dim
            moris::uint tParamDim = tVertexGroup->get_vertex_local_coords_dim();

            mGeneratedTemplate( aNumIgCells )->mParamCoords.resize( mGeneratedTemplate( aNumIgCells )->mNumNewNodes, Matrix< DDRMat >( 1, tParamDim ) );
            mGeneratedTemplate( aNumIgCells )->mNewVertexAncestry.mVertexParentEntityIndex.resize( mGeneratedTemplate( aNumIgCells )->mNumNewNodes, MORIS_INDEX_MAX );
            mGeneratedTemplate( aNumIgCells )->mNewVertexAncestry.mVertexParentEntityRank.resize( mGeneratedTemplate( aNumIgCells )->mNumNewNodes );
            mGeneratedTemplate( aNumIgCells )->mVertexHash.resize( mGeneratedTemplate( aNumIgCells )->mNumNewNodes );
            mGeneratedTemplate( aNumIgCells )->mIgCellToVertOrd.resize( mGeneratedTemplate( aNumIgCells )->mNumIgCells, Vector< moris_index >( mRegularSubdivisionTemplate->get_num_verts_per_cell(), MORIS_INDEX_MAX ) );

            moris_index tFirstNewFacetVertexOrdinal = tVertexGroup->size();
            moris_index tNewVertexOrdinal           = 0;

            Matrix< DDRMat > tLocalCoords( 1, tParamDim, 0.0 );

            // generate the template
            // place vertices on each new facet and evaluate the local coordinates
            std::cout << "TODO:GENERATE PROPER HASHES" << '\n';
            for ( moris::uint iFacet = 0; iFacet < tFaceConnectivity->mFacetVertices.size(); iFacet++ )
            {
                tLocalCoords.fill( 0.0 );
                // iterate all vertices on the face and sum up local coords
                for ( moris::uint iV = 0; iV < tFaceConnectivity->mFacetVertices( iFacet ).size(); iV++ )
                {
                    tLocalCoords.matrix_data() += tVertexGroup->get_vertex_local_coords( tFaceConnectivity->mFacetVertices( iFacet )( iV )->get_index() )->matrix_data();
                }

                mGeneratedTemplate( aNumIgCells )->mParamCoords( tNewVertexOrdinal ) = tLocalCoords / (moris::real)tFaceConnectivity->mFacetVertices( iFacet ).size();
                if ( tFacetAncestry->mFacetParentEntityRank( iFacet ) == 3 )
                {
                    mGeneratedTemplate( aNumIgCells )->mNewVertexAncestry.mVertexParentEntityIndex( tNewVertexOrdinal ) = 0;
                }
                else
                {
                    mGeneratedTemplate( aNumIgCells )->mNewVertexAncestry.mVertexParentEntityIndex( tNewVertexOrdinal ) = tCellConn.get_entity_ordinal( tFacetAncestry->mFacetParentEntityIndex( iFacet ), moris::mtk::get_entity_rank_from_index( tFacetAncestry->mFacetParentEntityRank( iFacet ) ) );
                }

                mGeneratedTemplate( aNumIgCells )->mNewVertexAncestry.mVertexParentEntityRank( tNewVertexOrdinal ) = moris::mtk::get_entity_rank_from_index( tFacetAncestry->mFacetParentEntityRank( iFacet ) );
                mGeneratedTemplate( aNumIgCells )->mVertexHash( tNewVertexOrdinal )                                = iFacet;
                tNewVertexOrdinal++;
            }

            Vector< Vector< moris::moris_index > > tIgCellToVertexTemplate    = this->get_ig_cell_to_vertex_connectivity();
            moris_index                            tCellIndex                 = 0;
            moris_index                            tFirstNewCellVertexOrdinal = tNewVertexOrdinal + tFirstNewFacetVertexOrdinal;
            for ( moris::uint iCell = 0; iCell < aNumIgCells; iCell++ )
            {
                tLocalCoords.fill( 0.0 );

                moris::mtk::Cell* tCell = aChildMesh->mIgCells->mIgCellGroup( iCell );

                Vector< moris::mtk::Vertex* > tVertices = tCell->get_vertex_pointers();

                // iterate all vertices on the face and sum up local coords
                for ( moris::uint iV = 0; iV < tVertices.size(); iV++ )
                {
                    tLocalCoords.matrix_data() += tVertexGroup->get_vertex_local_coords( tVertices( iV )->get_index() )->matrix_data();
                }

                mGeneratedTemplate( aNumIgCells )->mParamCoords( tNewVertexOrdinal )                                = tLocalCoords / (moris::real)tVertices.size();
                mGeneratedTemplate( aNumIgCells )->mNewVertexAncestry.mVertexParentEntityIndex( tNewVertexOrdinal ) = 0;
                mGeneratedTemplate( aNumIgCells )->mNewVertexAncestry.mVertexParentEntityRank( tNewVertexOrdinal )  = mtk::EntityRank::ELEMENT;
                mGeneratedTemplate( aNumIgCells )->mVertexHash( tNewVertexOrdinal )                                 = tNewVertexOrdinal;    // no special hashing needed here
                tNewVertexOrdinal++;

                // iterate through template and construct the new template cells
                for ( moris::moris_index iTemplateCells = 0; iTemplateCells < this->mRegularSubdivisionTemplate->get_num_ig_cells(); iTemplateCells++ )
                {
                    for ( moris::moris_index iTemplateVert = 0; iTemplateVert < this->mRegularSubdivisionTemplate->get_num_verts_per_cell(); iTemplateVert++ )
                    {
                        // base template index
                        moris_index tBaseTemplateVertexIndex = tIgCellToVertexTemplate( iTemplateCells )( iTemplateVert );

                        // node parent relative to the single ig cell
                        const moris_index     tParentEntityOrd  = aRegularSubdivisionInterfaceData->mVertexAncestry.get_vertex_parent_index( tBaseTemplateVertexIndex );
                        const mtk::EntityRank tParentEntityRank = aRegularSubdivisionInterfaceData->mVertexAncestry.get_vertex_parent_rank( tBaseTemplateVertexIndex );

                        if ( tParentEntityRank == mtk::EntityRank::NODE )
                        {
                            mGeneratedTemplate( aNumIgCells )->mIgCellToVertOrd( tCellIndex )( iTemplateVert ) = tVertexGroup->get_vertex_group_ordinal( tVertices( tParentEntityOrd )->get_index() );
                        }
                        else if ( tParentEntityRank == mtk::EntityRank::FACE )
                        {
                            moris_index tFaceIndex                                                             = tFaceConnectivity->mCellToFacet( iCell )( tParentEntityOrd );
                            mGeneratedTemplate( aNumIgCells )->mIgCellToVertOrd( tCellIndex )( iTemplateVert ) = tFirstNewFacetVertexOrdinal + tFaceIndex;
                        }
                        else
                        {
                            mGeneratedTemplate( aNumIgCells )->mIgCellToVertOrd( tCellIndex )( iTemplateVert ) = tFirstNewCellVertexOrdinal + iCell;
                        }
                    }
                    tCellIndex++;
                }
            }
        }
    }

    //--------------------------------------------------------------------------------------------------

}    // namespace moris::xtk
