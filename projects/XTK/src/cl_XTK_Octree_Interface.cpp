/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Octree_Interface.cpp
 *
 */

#include "cl_XTK_Octree_Interface.hpp"
#include "cl_XTK_Decomposition_Algorithm.hpp"
#include "cl_XTK_Regular_Subdivision_Interface.hpp"
#include "cl_XTK_Integration_Mesh_Generator.hpp"
#include "cl_XTK_Cut_Integration_Mesh.hpp"
#include "cl_MTK_Cell_Info_Factory.hpp"
#include "cl_MTK_Cell_Info.hpp"
#include "fn_norm.hpp"
namespace xtk
{
Vertex_Ancestry
IJK_Mesh::get_vertex_parent_entities() const
{

    Vertex_Ancestry tVertexAncestry;
    tVertexAncestry.mVertexParentEntityIndex.resize( this->num_verts() );
    tVertexAncestry.mVertexParentEntityRank.resize( this->num_verts() );

    moris::uint tI = 0;
    moris::uint tJ = 0;
    moris::uint tK = 0;

    // Vertex 0
    moris_index tIndex                                 = this->get_vertex_index( this->min_vert_i(), this->min_vert_j(), this->min_vert_k() );
    tVertexAncestry.mVertexParentEntityIndex( tIndex ) = 0;
    tVertexAncestry.mVertexParentEntityRank( tIndex )  = mtk::EntityRank::NODE;

    // Vertex 1
    tIndex                                             = this->get_vertex_index( this->max_vert_i(), this->min_vert_j(), this->min_vert_k() );
    tVertexAncestry.mVertexParentEntityIndex( tIndex ) = 1;
    tVertexAncestry.mVertexParentEntityRank( tIndex )  = mtk::EntityRank::NODE;

    // Vertex 2
    tIndex                                             = this->get_vertex_index( this->max_vert_i(), this->max_vert_j(), this->min_vert_k() );
    tVertexAncestry.mVertexParentEntityIndex( tIndex ) = 2;
    tVertexAncestry.mVertexParentEntityRank( tIndex )  = mtk::EntityRank::NODE;

    // Vertex 3
    tIndex                                             = this->get_vertex_index( this->min_vert_i(), this->max_vert_j(), this->min_vert_k() );
    tVertexAncestry.mVertexParentEntityIndex( tIndex ) = 3;
    tVertexAncestry.mVertexParentEntityRank( tIndex )  = mtk::EntityRank::NODE;

    // Vertex 4
    tIndex                                             = this->get_vertex_index( this->min_vert_i(), this->min_vert_j(), this->max_vert_k() );
    tVertexAncestry.mVertexParentEntityIndex( tIndex ) = 4;
    tVertexAncestry.mVertexParentEntityRank( tIndex )  = mtk::EntityRank::NODE;

    // Vertex 5
    tIndex                                             = this->get_vertex_index( this->max_vert_i(), this->min_vert_j(), this->max_vert_k() );
    tVertexAncestry.mVertexParentEntityIndex( tIndex ) = 5;
    tVertexAncestry.mVertexParentEntityRank( tIndex )  = mtk::EntityRank::NODE;

    // Vertex 6
    tIndex                                             = this->get_vertex_index( this->max_vert_i(), this->max_vert_j(), this->max_vert_k() );
    tVertexAncestry.mVertexParentEntityIndex( tIndex ) = 6;
    tVertexAncestry.mVertexParentEntityRank( tIndex )  = mtk::EntityRank::NODE;

    // Vertex 7
    tIndex                                             = this->get_vertex_index( this->min_vert_i(), this->max_vert_j(), this->max_vert_k() );
    tVertexAncestry.mVertexParentEntityIndex( tIndex ) = 7;
    tVertexAncestry.mVertexParentEntityRank( tIndex )  = mtk::EntityRank::NODE;

    // Face 0 Interior Vertices - Variable i and k - fixed j
    tJ = this->min_vert_j();
    for ( tI = 1; tI < this->num_vert_x() - 1; tI++ )
    {
        for ( tK = 1; tK < this->num_vert_z() - 1; tK++ )
        {
            tIndex                                             = this->get_vertex_index( tI, tJ, tK );
            tVertexAncestry.mVertexParentEntityIndex( tIndex ) = 0;
            tVertexAncestry.mVertexParentEntityRank( tIndex )  = mtk::EntityRank::FACE;
        }
    }

    // Face 1 Interior Vertices - Variable j and k - fixed i+1
    tI = this->max_vert_i();
    for ( tJ = 1; tJ < this->num_vert_y() - 1; tJ++ )
    {
        for ( tK = 1; tK < this->num_vert_z() - 1; tK++ )
        {
            tIndex                                             = this->get_vertex_index( tI, tJ, tK );
            tVertexAncestry.mVertexParentEntityIndex( tIndex ) = 1;
            tVertexAncestry.mVertexParentEntityRank( tIndex )  = mtk::EntityRank::FACE;
        }
    }

    // Face 2 Interior Vertices - Variable i and k - fixed j+1
    tJ = this->max_vert_j();
    for ( tI = 1; tI < this->num_vert_x() - 1; tI++ )
    {
        for ( tK = 1; tK < this->num_vert_z() - 1; tK++ )
        {
            tIndex                                             = this->get_vertex_index( tI, tJ, tK );
            tVertexAncestry.mVertexParentEntityIndex( tIndex ) = 2;
            tVertexAncestry.mVertexParentEntityRank( tIndex )  = mtk::EntityRank::FACE;
        }
    }

    // Face 3 Interior Vertices - variable j and k - fixed i
    tI = this->min_vert_i();
    for ( tJ = 1; tJ < this->num_vert_y() - 1; tJ++ )
    {
        for ( tK = 1; tK < this->num_vert_z() - 1; tK++ )
        {
            tIndex                                             = this->get_vertex_index( tI, tJ, tK );
            tVertexAncestry.mVertexParentEntityIndex( tIndex ) = 3;
            tVertexAncestry.mVertexParentEntityRank( tIndex )  = mtk::EntityRank::FACE;
        }
    }
    // Face 4 Interior Vertices - variable i and j - fixed k
    tK = this->min_vert_k();
    for ( tI = 1; tI < this->num_vert_x() - 1; tI++ )
    {
        for ( tJ = 1; tJ < this->num_vert_y() - 1; tJ++ )
        {
            tIndex                                             = this->get_vertex_index( tI, tJ, tK );
            tVertexAncestry.mVertexParentEntityIndex( tIndex ) = 4;
            tVertexAncestry.mVertexParentEntityRank( tIndex )  = mtk::EntityRank::FACE;
        }
    }
    // Face 5 Interior Vertices - variable i and j - fixed k+1
    tK = this->max_vert_k();
    for ( tI = 1; tI < this->num_vert_x() - 1; tI++ )
    {
        for ( tJ = 1; tJ < this->num_vert_y() - 1; tJ++ )
        {
            tIndex                                             = this->get_vertex_index( tI, tJ, tK );
            tVertexAncestry.mVertexParentEntityIndex( tIndex ) = 5;
            tVertexAncestry.mVertexParentEntityRank( tIndex )  = mtk::EntityRank::FACE;
        }
    }

    // EDGES WITH VARIABLE I
    // fixed J and K edges
    tJ = 0;// fixed for edge 0
    tK = 0;// fixed for edge 0
    for ( tI = 1; tI < this->num_vert_x() - 1; tI++ )
    {
        tIndex                                             = this->get_vertex_index( tI, tJ, tK );
        tVertexAncestry.mVertexParentEntityIndex( tIndex ) = 0;
        tVertexAncestry.mVertexParentEntityRank( tIndex )  = mtk::EntityRank::EDGE;
    }

    tJ = this->max_vert_j();
    tK = 0;
    for ( tI = 1; tI < this->num_vert_x() - 1; tI++ )
    {
        tIndex                                             = this->get_vertex_index( tI, tJ, tK );
        tVertexAncestry.mVertexParentEntityIndex( tIndex ) = 2;
        tVertexAncestry.mVertexParentEntityRank( tIndex )  = mtk::EntityRank::EDGE;
    }

    tJ = 0;
    tK = this->max_vert_k();
    for ( tI = 1; tI < this->num_vert_x() - 1; tI++ )
    {
        tIndex                                             = this->get_vertex_index( tI, tJ, tK );
        tVertexAncestry.mVertexParentEntityIndex( tIndex ) = 4;
        tVertexAncestry.mVertexParentEntityRank( tIndex )  = mtk::EntityRank::EDGE;
    }

    tJ = this->max_vert_j();
    tK = this->max_vert_k();
    for ( tI = 1; tI < this->num_vert_x() - 1; tI++ )
    {
        tIndex                                             = this->get_vertex_index( tI, tJ, tK );
        tVertexAncestry.mVertexParentEntityIndex( tIndex ) = 6;
        tVertexAncestry.mVertexParentEntityRank( tIndex )  = mtk::EntityRank::EDGE;
    }

    // EDGES WITH VARIABLE J
    tI = 0;
    tK = 0;
    for ( tJ = 1; tJ < this->num_vert_y() - 1; tJ++ )
    {
        tIndex                                             = this->get_vertex_index( tI, tJ, tK );
        tVertexAncestry.mVertexParentEntityIndex( tIndex ) = 3;
        tVertexAncestry.mVertexParentEntityRank( tIndex )  = mtk::EntityRank::EDGE;
    }

    tI = this->max_vert_i();
    tK = 0;
    for ( tJ = 1; tJ < this->num_vert_y() - 1; tJ++ )
    {
        tIndex                                             = this->get_vertex_index( tI, tJ, tK );
        tVertexAncestry.mVertexParentEntityIndex( tIndex ) = 1;
        tVertexAncestry.mVertexParentEntityRank( tIndex )  = mtk::EntityRank::EDGE;
    }

    tI = 0;
    tK = this->max_vert_k();
    for ( tJ = 1; tJ < this->num_vert_y() - 1; tJ++ )
    {
        tIndex                                             = this->get_vertex_index( tI, tJ, tK );
        tVertexAncestry.mVertexParentEntityIndex( tIndex ) = 7;
        tVertexAncestry.mVertexParentEntityRank( tIndex )  = mtk::EntityRank::EDGE;
    }

    tI = this->max_vert_k();
    tK = this->max_vert_k();
    for ( tJ = 1; tJ < this->num_vert_y() - 1; tJ++ )
    {
        tIndex                                             = this->get_vertex_index( tI, tJ, tK );
        tVertexAncestry.mVertexParentEntityIndex( tIndex ) = 5;
        tVertexAncestry.mVertexParentEntityRank( tIndex )  = mtk::EntityRank::EDGE;
    }

    // VARIABLE K
    tI = 0;
    tJ = 0;
    for ( tK = 1; tK < this->num_vert_z() - 1; tK++ )
    {
        tIndex                                             = this->get_vertex_index( tI, tJ, tK );
        tVertexAncestry.mVertexParentEntityIndex( tIndex ) = 8;
        tVertexAncestry.mVertexParentEntityRank( tIndex )  = mtk::EntityRank::EDGE;
    }

    tI = this->max_vert_i();
    tJ = 0;
    for ( tK = 1; tK < this->num_vert_z() - 1; tK++ )
    {
        tIndex                                             = this->get_vertex_index( tI, tJ, tK );
        tVertexAncestry.mVertexParentEntityIndex( tIndex ) = 9;
        tVertexAncestry.mVertexParentEntityRank( tIndex )  = mtk::EntityRank::EDGE;
    }

    tI = 0;
    tJ = this->max_vert_j();
    for ( tK = 1; tK < this->num_vert_z() - 1; tK++ )
    {
        tIndex                                             = this->get_vertex_index( tI, tJ, tK );
        tVertexAncestry.mVertexParentEntityIndex( tIndex ) = 11;
        tVertexAncestry.mVertexParentEntityRank( tIndex )  = mtk::EntityRank::EDGE;
    }

    tI = this->max_vert_i();
    tJ = this->max_vert_j();
    for ( tK = 1; tK < this->num_vert_z() - 1; tK++ )
    {
        tIndex                                             = this->get_vertex_index( tI, tJ, tK );
        tVertexAncestry.mVertexParentEntityIndex( tIndex ) = 10;
        tVertexAncestry.mVertexParentEntityRank( tIndex )  = mtk::EntityRank::EDGE;
    }

    // flag the interior vertices
    for ( moris::uint k = 1; k < this->num_vert_z() - 1; k++ )
    {
        for ( moris::uint j = 1; j < this->num_vert_y() - 1; j++ )
        {
            for ( moris::uint i = 1; i < this->num_vert_x() - 1; i++ )
            {
                tIndex                                             = this->get_vertex_index( i, j, k );
                tVertexAncestry.mVertexParentEntityIndex( tIndex ) = 0;
                tVertexAncestry.mVertexParentEntityRank( tIndex )  = mtk::EntityRank::ELEMENT;
            }
        }
    }

    return tVertexAncestry;
}

Octree_Interface::Octree_Interface( ParameterList &aParameterList )
{
    mOctreeRefinementLevel = std::stoi( aParameterList.get< std::string >( "octree_refinement_level" ) );
    std::cout << "mOctreeRefinementLevel = " << mOctreeRefinementLevel << std::endl;
}
// ----------------------------------------------------------------------------------

// ----------------------------------------------------------------------------------
void
Octree_Interface::perform(
    Integration_Mesh_Generation_Data *aMeshGenerationData,
    Decomposition_Data *              aDecompositionData,
    Cut_Integration_Mesh *            aCutIntegrationMesh,
    moris::mtk::Mesh *                aBackgroundMesh,
    Integration_Mesh_Generator *      aMeshGenerator )
{
    // store the inputs as member data
    mIgMeshGenData  = aMeshGenerationData;
    mDecompData     = aDecompositionData;
    mCutIgMesh      = aCutIntegrationMesh;
    mBackgroundMesh = aBackgroundMesh;
    mIgMeshGen      = aMeshGenerator;

    Tracer tTracer( "XTK", "Decomposition_Algorithm", "Octree" );

    // generate the local tensor grid templates
    mOctreeTemplates = this->generate_octree_templates();

    mIgVertexGroupIndexToIjkIndex = this->generate_octree_template_vertex_group_to_ijk_map();

    this->perform_impl_vertex_requests( aMeshGenerationData, aDecompositionData, aCutIntegrationMesh, aBackgroundMesh, aMeshGenerator );

    // in parallel give all these nodes ids
    aMeshGenerator->assign_node_requests_identifiers( *aDecompositionData, aCutIntegrationMesh, aBackgroundMesh );

    // commit vertices to the mesh
    aMeshGenerator->commit_new_ig_vertices_to_cut_mesh( aMeshGenerationData, aDecompositionData, aCutIntegrationMesh, aBackgroundMesh, this );

    // aDecompositionData->print_requests(*aBackgroundMesh);
    // aDecompositionData->print(*aBackgroundMesh);

    // generate_decomposed_mesh - perform the mesh generation
    this->perform_impl_generate_mesh( aMeshGenerationData, aDecompositionData, aCutIntegrationMesh, aBackgroundMesh, aMeshGenerator );

    // commit the cells to the mesh
    aMeshGenerator->commit_new_ig_cells_to_cut_mesh( aMeshGenerationData, aDecompositionData, aCutIntegrationMesh, aBackgroundMesh, this );

    // std::cout<<"tOctreeTemplates.size() == "<<tOctreeTemplates.size()<<std::endl;
}

// ----------------------------------------------------------------------------------

enum Decomposition_Algorithm_Type
Octree_Interface::get_algorithm_type() const
{
    return Decomposition_Algorithm_Type::OCTREE;
}
// ----------------------------------------------------------------------------------
moris_index
Octree_Interface::get_signature() const
{
    return 8;
}
// ----------------------------------------------------------------------------------
bool
Octree_Interface::has_geometric_independent_vertices() const
{
    return true;
}
// ----------------------------------------------------------------------------------
void
Octree_Interface::perform_impl_vertex_requests(
    Integration_Mesh_Generation_Data *aMeshGenerationData,
    Decomposition_Data *              aDecompositionData,
    Cut_Integration_Mesh *            aCutIntegrationMesh,
    moris::mtk::Mesh *                aBackgroundMesh,
    Integration_Mesh_Generator *      aMeshGenerator )
{
    mNumNewCells = 0;

    aDecompositionData->mHasSecondaryIdentifier = true;
    aDecompositionData->tDecompId               = this->get_signature();

    // size the CM New Node Loc and param coord data
    moris_index const tNumChildMeshes        = aCutIntegrationMesh->get_num_ig_cell_groups();
    aDecompositionData->tCMNewNodeLoc        = Cell< Cell< moris_index > >( tNumChildMeshes );
    aDecompositionData->tCMNewNodeParamCoord = Cell< Cell< Matrix< DDRMat > > >( tNumChildMeshes );

    // initialize data
    moris::Matrix< moris::DDRMat > tBasisWeights( 1, 8 );

    mDifference.reserve( tNumChildMeshes );

    // iterate through cells
    for ( auto &iCell : aMeshGenerationData->mAllIntersectedBgCellInds )
    {
        // background cell
        moris::mtk::Cell &tCell = aBackgroundMesh->get_mtk_cell( iCell );

        // background cell vertices
        moris::Matrix< moris::DDRMat > tBGCellCoords = tCell.get_vertex_coords();

        // level of this cell
        moris_index tLevel = tCell.get_level();

        mDifference.push_back( mOctreeRefinementLevel - tLevel );

        // specific template for this cell
        std::shared_ptr< Octree_Template const > tCurrentTemplate = mOctreeTemplates( mOctreeRefinementLevel - tLevel );

        mNumNewCells = mNumNewCells + tCurrentTemplate->get_mesh_grid()->num_cells();

        MORIS_ASSERT( tCurrentTemplate != nullptr, "Selected a uninitialized template" );

        // number of vertices in the selected grid
        const moris::uint tNumVerts = tCurrentTemplate->get_mesh_grid()->num_verts();

        // vertex ancestr
        const Vertex_Ancestry *tVertexAncestry = tCurrentTemplate->get_vertex_ancestry();

        // hashes for the current template
        const moris::Cell< moris_index > *tVertexHashes = tCurrentTemplate->get_vertex_hashes();

        Cell_Connectivity const tCellConnectivity = mCutIgMesh->get_background_cell_connectivity( tCell.get_index() );

        // parameteric coordinates for the current template
        const Matrix< DDRMat > &tVertexParamCoords = tCurrentTemplate->get_vertex_param_coords();

        for ( moris::uint iV = 0; iV < tNumVerts; iV++ )
        {

            moris::Matrix< moris::DDRMat > tBasisWeights;

            // local coordinate of this vertex wrt the current cell group
            std::shared_ptr< Matrix< DDRMat > > tVertexLocalCoords = std::make_shared< moris::Matrix< moris::DDRMat > >( tVertexParamCoords.get_row( iV ) );

            // evaluate the basis function
            tCell.get_cell_info()->eval_N( *tVertexLocalCoords, tBasisWeights );

            moris::Matrix< moris::DDRMat > tNewCoordinate = tBasisWeights * tBGCellCoords;

            if ( tVertexAncestry->get_vertex_parent_rank( iV ) != mtk::EntityRank::NODE )
            {
                moris_index tNewNodeIndexInSubdivision = MORIS_INDEX_MAX;

                if ( tVertexAncestry->get_vertex_parent_rank( iV ) == mtk::EntityRank::ELEMENT )
                {
                    if ( !mDecompData->request_exists( tCell.get_index(), ( *tVertexHashes )( iV ), mtk::EntityRank::ELEMENT, tNewNodeIndexInSubdivision ) )
                    {
                        tNewNodeIndexInSubdivision = mDecompData->register_new_request(
                            tCell.get_index(),
                            ( *tVertexHashes )( iV ),
                            tCell.get_owner(),
                            mtk::EntityRank::ELEMENT,
                            tNewCoordinate,
                            &tCell,
                            tVertexLocalCoords );
                    }
                    else
                    {
                        MORIS_ASSERT( moris::norm( tNewCoordinate - mDecompData->tNewNodeCoordinate( tNewNodeIndexInSubdivision ) ) < 1e-12,
                             "Vertex coordinate mismatch, could be a hashing collision, parent rank: ELEMENT , index: %i  hash: %i", tCell.get_index(), ( *tVertexHashes )( iV ) );
                    }
                }

                else
                {
                    if ( !mDecompData->request_exists(
                             tCellConnectivity.get_entity_index( tVertexAncestry->get_vertex_parent_index( iV ), tVertexAncestry->get_vertex_parent_rank( iV ) ),
                             ( *tVertexHashes )( iV ),
                             tVertexAncestry->get_vertex_parent_rank( iV ),
                             tNewNodeIndexInSubdivision ) )
                    {
                        moris_index tOwningProc    = mBackgroundMesh->get_entity_owner( tVertexAncestry->get_vertex_parent_index( iV ), tVertexAncestry->get_vertex_parent_rank( iV ) );
                        tNewNodeIndexInSubdivision = mDecompData->register_new_request(
                            tCellConnectivity.get_entity_index( tVertexAncestry->get_vertex_parent_index( iV ), tVertexAncestry->get_vertex_parent_rank( iV ) ),
                            ( *tVertexHashes )( iV ),
                            tOwningProc,
                            tVertexAncestry->get_vertex_parent_rank( iV ),
                            tNewCoordinate,
                            &tCell,
                            tVertexLocalCoords );
                    }
                    else
                    {
                        MORIS_ASSERT( moris::norm( tNewCoordinate - mDecompData->tNewNodeCoordinate( tNewNodeIndexInSubdivision ) ) < 1e-12,
                            "Vertex coordinate mismatch, could be a hashing collision, parent rank: %s, parent index: %d, hash: %d",
                            get_enum_str( tVertexAncestry->get_vertex_parent_rank( iV ) ).c_str(),
                            tCellConnectivity.get_entity_index( tVertexAncestry->get_vertex_parent_index( iV ), tVertexAncestry->get_vertex_parent_rank( iV ) ),
                            ( *tVertexHashes )( iV ) );
                    }
                }
                aDecompositionData->tCMNewNodeParamCoord( iCell ).push_back( *tVertexLocalCoords );
                aDecompositionData->tCMNewNodeLoc( iCell ).push_back( tNewNodeIndexInSubdivision );
            }

            // make a vertex request
        }
    }
}
// ----------------------------------------------------------------------------------
void
Octree_Interface::perform_impl_generate_mesh(
    Integration_Mesh_Generation_Data *aMeshGenerationData,
    Decomposition_Data *              aDecompositionData,
    Cut_Integration_Mesh *            aCutIntegrationMesh,
    moris::mtk::Mesh *                aBackgroundMesh,
    Integration_Mesh_Generator *      aMeshGenerator )
{

    // new cell info
    moris::mtk::Cell_Info_Factory            tFactory;
    std::shared_ptr< moris::mtk::Cell_Info > tIgCellInfo = tFactory.create_cell_info_sp( this->get_ig_cell_topology() );

    // this is going to be the cell infor pointer for all new cells
    mNewCellCellInfo = moris::Cell< std::shared_ptr< moris::mtk::Cell_Info > >( mNumNewCells, tIgCellInfo );

    // number of vertices per cell
    moris::uint tVerticesPerCell = tIgCellInfo->get_num_verts();

    // allocate data in the new ig cell data
    mNewCellToVertexConnectivity = moris::Cell< moris::Cell< moris::moris_index > >( mNumNewCells, tVerticesPerCell );
    mNewCellChildMeshIndex       = moris::Cell< moris::moris_index >( mNumNewCells );

    // for this method we are not going to replace any cells, this is because the cell we would replace corresponds to the background mesh cell
    mNewCellCellIndexToReplace = moris::Cell< moris::moris_index >( mNumNewCells, MORIS_INDEX_MAX );

    // populate new cell data
    moris::moris_index tCurrentCellIndex = 0;
    for ( moris::moris_index iCM = 0; iCM < aMeshGenerationData->mNumChildMeshes; iCM++ )
    {
        std::shared_ptr< Child_Mesh_Experimental > tChildMesh = aCutIntegrationMesh->get_child_mesh( iCM );

        moris::uint tDiff = mOctreeRefinementLevel - tChildMesh->get_parent_cell()->get_level();

        std::shared_ptr< Octree_Template const > tTemplate = mOctreeTemplates( tDiff );

        const moris::Matrix< IndexMat > *tCellTemplates = tTemplate->get_cells();

        std::cout << "iCM = " << iCM << " | tChildMesh->get_parent_cell()->get_level() = " << tChildMesh->get_parent_cell()->get_level() << std::endl;

        for ( moris::uint iNewCell = 0; iNewCell < tCellTemplates->n_rows(); iNewCell++ )
        {
            mNewCellChildMeshIndex( tCurrentCellIndex ) = iCM;

            for ( moris::uint iV = 0; iV < tCellTemplates->n_cols(); iV++ )
            {
                // Ijk index
                moris_index tIgCellIJKIndex = ( *tCellTemplates )( iNewCell, iV );

                // std::cout << " mDifference( iCM ) = " << mDifference( iCM ) << std::endl;
                // index in child mesh vertex group
                moris_index tNewVertexCMOrdinal = mIgVertexGroupIndexToIjkIndex( tDiff ).find( tIgCellIJKIndex )->second;
                MORIS_ERROR( tNewVertexCMOrdinal < (moris::moris_index)tChildMesh->mIgVerts->size(), "Template ordinal out of bounds" );
                mNewCellToVertexConnectivity( tCurrentCellIndex )( iV ) = tChildMesh->mIgVerts->get_vertex( tNewVertexCMOrdinal )->get_index();
            }

            tCurrentCellIndex++;
        }
    }
}
// ----------------------------------------------------------------------------------
moris::Cell< std::shared_ptr< Octree_Template const > >
Octree_Interface::generate_octree_templates()
{
    // figure out the lowest and highest octree templates to create
    moris::Cell< moris_index > const tOctreeBounds = this->determine_octree_bounds();

    moris::Cell< std::shared_ptr< Octree_Template const > > tTemplate( mOctreeRefinementLevel + 1, nullptr );

    // highest level is used for vertex hashing in other templates
    tTemplate( mOctreeRefinementLevel ) = std::make_shared< Octree_Template const >( mOctreeRefinementLevel );

    // iterate from low to high bounds
    for ( moris_index iBounds = 0; iBounds < mOctreeRefinementLevel; iBounds++ )
    {
        std::cout << "iBounds = " << iBounds << std::endl;

        tTemplate( iBounds ) = std::make_shared< Octree_Template const >( iBounds );
    }

    std::cout << "tMaxLevel = " << tOctreeBounds( 1 ) << std::endl;
    std::cout << "tMinLevel = " << tOctreeBounds( 0 ) << std::endl;

    return tTemplate;
}

moris::Cell< moris_index > const
Octree_Interface::determine_octree_bounds()
{
    moris_index tMaxLevel = 0;
    moris_index tMinLevel = 0;

    // iterate through cells
    for ( moris::uint iGeom = 0; iGeom < mIgMeshGenData->mIntersectedBackgroundCellIndex.size(); iGeom++ )
    {
        for ( moris::uint iCell = 0; iCell < mIgMeshGenData->mIntersectedBackgroundCellIndex( iGeom ).size(); iCell++ )
        {
            const moris_index tCellLevel = mBackgroundMesh->get_mtk_cell( mIgMeshGenData->mIntersectedBackgroundCellIndex( iGeom )( iCell ) ).get_level();

            if ( tCellLevel > tMaxLevel ) { tMaxLevel = tCellLevel; }
            if ( tCellLevel < tMinLevel ) { tMinLevel = tCellLevel; }
        }
    }

    return { tMinLevel, tMaxLevel };
}
// ----------------------------------------------------------------------------------

mtk::CellTopology
Octree_Interface::get_ig_cell_topology() const
{
    std::cout << " WARNING NEED ABSTRACTION " << std::endl;
    return mtk::CellTopology::HEX8;
}

moris::Cell< std::unordered_map< moris_index, moris_index > >
Octree_Interface::generate_octree_template_vertex_group_to_ijk_map()
{

    moris::Cell< std::unordered_map< moris_index, moris_index > > tMaps( mOctreeTemplates.size() );

    // iterate through templates
    for ( moris::uint iTempl = 0; iTempl < mOctreeTemplates.size(); iTempl++ )
    {
        if ( mOctreeTemplates( iTempl ) != nullptr )
        {
            moris_index            tCurrentIndex   = mOctreeTemplates( iTempl )->get_num_corner_points();
            Vertex_Ancestry const *tVertexAncestry = mOctreeTemplates( iTempl )->get_vertex_ancestry();

            for ( moris::uint iV = 0; iV < mOctreeTemplates( iTempl )->get_mesh_grid()->num_verts(); iV++ )
            {
                if ( tVertexAncestry->get_vertex_parent_rank( iV ) == mtk::EntityRank::NODE )
                {
                    tMaps( iTempl )[iV] = tVertexAncestry->get_vertex_parent_index( iV );
                }
                else
                {
                    tMaps( iTempl )[iV] = tCurrentIndex++;
                }
            }
        }
    }

    return tMaps;
}

}// namespace xtk

