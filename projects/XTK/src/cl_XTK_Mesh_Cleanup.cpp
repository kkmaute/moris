/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Mesh_Cleanup.cpp
 *
 */

#include "cl_XTK_Mesh_Cleanup.hpp"
#include "cl_XTK_Model.hpp"
#include "cl_Parameter_List.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_XTK_Child_Mesh.hpp"
#include "cl_Tracer.hpp"
namespace moris::xtk
{
    Mesh_Cleanup::Mesh_Cleanup( Model* aModel,
            moris::Parameter_List*      aParamList )
            : mModel( aModel )
    {
        mMeshCleanupParameters.mDeactivateOneBPChildMeshes = aParamList->get< bool >( "cleanup_cut_mesh" );
    }

    void
    Mesh_Cleanup::perform()
    {
        // remove inactive child meshes ones that do not have an interior interface
        if ( mMeshCleanupParameters.mDeactivateOneBPChildMeshes )
        {
            this->cleanup_cut_mesh();
        }
    }

    void
    Mesh_Cleanup::cleanup_cut_mesh()
    {
        Tracer tTracer( "XTK", "Mesh Cleanup", "Remove Inactive Child Mesh" );

        // Using a map here so that I can remove some from the removal process.
        Vector< moris_index >                          tChildMeshesToDelete;
        std::unordered_map< moris_index, moris_index > tChildMeshesToDeleteMap;

        // select candidate cell groups for deletion
        this->select_candidate_child_meshes_for_cleanup( tChildMeshesToDeleteMap );

        // remove the child meshes with interfaces from candidate list

        // place the map keys in a vector
        this->get_vector_of_child_meshes_for_removal( tChildMeshesToDeleteMap, tChildMeshesToDelete );

        // remove subphases associated with child

        // remove child meshes that have an interface between two bulk phases from the candidates

        // tell the cut integration mesh to remove some of the child meshes

        // for(moris::uint iCM = 0; iCM < tNumCutMeshes; iCM++)
        // {
        //     Child_Mesh & tCM = mModel->get_cut_mesh().get_child_mesh(iCM);
        //     Cell<moris::moris_index> const & tSubphasebinBulkPhase = tCM.get_subphase_bin_bulk_phase();

        //    if(tSubphasebinBulkPhase.size() > 1 or tCM.has_inter_child_mesh_interfaces())
        //    {
        //         tChildMeshesToKeep.push_back(iCM);
        //    }
        //    else
        //    {
        //         tChildMeshesToDelete.push_back(iCM);
        //    }
        // }

        // Vector<moris_index> tCellsToRemoveFromMesh;
        // mModel->get_cut_mesh().remove_all_child_meshes_but_selected(tChildMeshesToKeep,tChildMeshesToDelete, tCellsToRemoveFromMesh );

        // // moris::print(tCellsToRemoveFromMesh,"tCellsToRemoveFromMesh");
        // Vector<moris_index> tNewCellIndices;
        // mModel->get_background_mesh().remove_cells_from_mesh(tCellsToRemoveFromMesh,tNewCellIndices);

        // // reindex the child mesh
        // mModel->get_cut_mesh().reindex_cells(tNewCellIndices);

        MORIS_LOG_SPEC( "Num Child Meshes Removed", tChildMeshesToDelete.size() );
        MORIS_LOG_SPEC( "Num Child Meshes Kept", mModel->get_cut_integration_mesh()->get_num_child_meshes() - tChildMeshesToDelete.size() );

        MORIS_ERROR( 0, "Intentional Bail" );
    }

    void
    Mesh_Cleanup::select_candidate_child_meshes_for_cleanup(
            std::unordered_map< moris_index, moris_index >& aRemoveChildMeshes )
    {

        // cut integration mesh
        Cut_Integration_Mesh* tCutIgMesh = mModel->get_cut_integration_mesh();

        std::cout << "Num Child Meshes = " << tCutIgMesh->get_num_child_meshes() << '\n';
        std::cout << "Num Background Cells = " << mModel->get_background_mesh().get_num_elems() << '\n';
        // iterate through ig cell groups
        for ( moris::uint iCM = 0; iCM < tCutIgMesh->get_num_child_meshes(); iCM++ )
        {
            std::shared_ptr< Child_Mesh_Experimental > tChildMesh = tCutIgMesh->get_child_mesh( iCM );

            // 1 here ignores the child meshes that are empty for non-intersecte background cells
            if ( tChildMesh->get_num_subphase_cell_groups() == 1 )
            {
                std::cout << "(tChildMesh->get_num_subphase_cell_groups() = " << tChildMesh->get_num_subphase_cell_groups() << '\n';
                std::shared_ptr< IG_Cell_Group > tSubphaseIgCells = tChildMesh->get_subphase_cell_group( 0 );

                std::cout << "tSubphaseIgCells size = " << tSubphaseIgCells->mIgCellGroup.size() << '\n';
                aRemoveChildMeshes[ iCM ] = 1;
            }
        }
    }
    void
    Mesh_Cleanup::get_vector_of_child_meshes_for_removal(
            std::unordered_map< moris_index, moris_index > const & aChildMeshesToDeleteMap,
            Vector< moris_index >&                                 aChildMeshesToDelete )
    {
        aChildMeshesToDelete.clear();
        aChildMeshesToDelete.reserve( aChildMeshesToDeleteMap.size() );
        for ( auto const & iIter : aChildMeshesToDeleteMap )
        {
            if ( iIter.second == 1 )
            {
                aChildMeshesToDelete.push_back( iIter.first );
            }
        }
    }

    // void
    // Mesh_Cleanup::save_candidate_child_meshes
}    // namespace moris::xtk
