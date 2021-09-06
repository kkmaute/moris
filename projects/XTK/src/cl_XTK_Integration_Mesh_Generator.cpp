#include "cl_XTK_Integration_Mesh_Generator.hpp"

#include "cl_XTK_Regular_Subdivision_Interface.hpp"
#include "cl_XTK_Cell_No_CM.hpp"

namespace xtk
{
        bool
        Integration_Mesh_Generator::perform()
        {
            Integration_Mesh_Generation_Data tGenerationData;

            moris::mtk::Mesh* tBackgroundMesh = &mXTKModel->get_background_mesh().get_mesh_data();

            Cut_Integration_Mesh tCutIntegrationMesh(tBackgroundMesh);

            bool tDetermineIntersectedCellsByGeom = this->determine_intersected_background_cells(tGenerationData,&tCutIntegrationMesh,tBackgroundMesh);

            if(!tDetermineIntersectedCellsByGeom)
            {
                return tDetermineIntersectedCellsByGeom;
            }

            bool tAllocatedCM = this->allocate_child_meshes(tGenerationData,&tCutIntegrationMesh,tBackgroundMesh);

            if(!tAllocatedCM)
            {
                return tAllocatedCM;
            }
            
            // perform the regular subdivision
            Decomposition_Algorithm* tRegSubMethod = new Regular_Subdivision_Interface();
            Decomposition_Data tDecompositionData;
            tRegSubMethod->perform(&tGenerationData, &tDecompositionData,&tCutIntegrationMesh, tBackgroundMesh, this);
            delete tRegSubMethod;

            // perform the conformal subdivision
            this->determine_intersected_ig_edges(tGenerationData,&tCutIntegrationMesh);


            return true;
        }

        void
        Integration_Mesh_Generator::commit_new_ig_cells_to_cut_mesh(
            Integration_Mesh_Generation_Data*  aMeshGenerationData,
            Decomposition_Data *               aDecompositionData,
            Cut_Integration_Mesh *             aCutIntegrationMesh,
            moris::mtk::Mesh *                 aBackgroundMesh,
            Decomposition_Algorithm*           aDecompositionAlgorithm)
        {
            // iterate through cells that the decomposition constructed
            moris::uint tNumNewCells = aDecompositionAlgorithm->mNumNewCells;

            MORIS_ERROR(aDecompositionAlgorithm->mNewCellToVertexConnectivity.size() == aDecompositionAlgorithm->mNewCellChildMeshIndex.size()  &&
                        aDecompositionAlgorithm->mNewCellChildMeshIndex.size() == aDecompositionAlgorithm->mNewCellCellIndexToReplace.size()
                        ,"Inconsistent size from decomposition algorithm");

            // add space to the mesh
            moris::uint tNumStartingCellsControlled = aCutIntegrationMesh->mControlledIgCells.size();
            moris::uint tNumStartingTotalIgCells    = aCutIntegrationMesh->mIntegrationCells.size();

            aCutIntegrationMesh->mControlledIgCells.resize(tNumStartingCellsControlled + tNumNewCells);
            aCutIntegrationMesh->mIntegrationCells.resize(tNumStartingTotalIgCells + tNumNewCells);

            // current index
            moris_index tCellIndex = tNumStartingTotalIgCells;

            // iterate through new and add to the mesh
            for(moris::uint iCell = 0;  iCell < aDecompositionAlgorithm->mNewCellToVertexConnectivity.size(); iCell++ )
            {
                // child mesh pointer
                std::shared_ptr<Child_Mesh_Experimental> tCM = aCutIntegrationMesh->get_child_mesh(aDecompositionAlgorithm->mNewCellChildMeshIndex(iCell));

                // parent cell owner
                moris_index tOwner = tCM->get_parent_cell()->get_owner();

                // collect the vertex pointers for the cell
                moris::Cell<moris::mtk::Vertex*> tVertexPointers(aDecompositionAlgorithm->mNewCellToVertexConnectivity(iCell).size());
                for(moris::uint iV = 0; iV < aDecompositionAlgorithm->mNewCellToVertexConnectivity(iCell).size(); iV++)
                {
                    tVertexPointers(iV) = aCutIntegrationMesh->get_mtk_vertex_pointer( aDecompositionAlgorithm->mNewCellToVertexConnectivity(iCell)(iV));
                }

                bool tReplaceExistingCell = aDecompositionAlgorithm->mNewCellCellIndexToReplace(iCell) != MORIS_INDEX_MAX;

                // cell index (if I replace one that is the index of this cell)
                moris_index tNewCellIndex = tReplaceExistingCell? aDecompositionAlgorithm->mNewCellCellIndexToReplace(iCell) : tCellIndex++;

                // create the new cell no id
                std::shared_ptr<moris::mtk::Cell> tNewCell = std::make_shared<xtk::Cell_XTK_No_CM>(
                    MORIS_ID_MAX,tNewCellIndex,tOwner,aDecompositionAlgorithm->mNewCellCellInfo(iCell),tVertexPointers);

                // add the cell to the mesh
                aCutIntegrationMesh->set_integration_cell(tNewCellIndex,tNewCell,!tReplaceExistingCell);

                // add the cell to a child mesh
                aCutIntegrationMesh->add_cell_to_integration_mesh(tNewCellIndex,tCM->mChildMeshIndex);
                
            }        

        }
}