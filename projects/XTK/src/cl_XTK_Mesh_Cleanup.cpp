#include "cl_XTK_Mesh_Cleanup.hpp"
#include "cl_XTK_Model.hpp"
#include "cl_Param_List.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_XTK_Child_Mesh.hpp"
#include "cl_Tracer.hpp"
namespace xtk
{
    Mesh_Cleanup::Mesh_Cleanup(Model*                aModel,
                               moris::ParameterList* aParamList):
    mModel(aModel)
    {
        mMeshCleanupParameters.mDeactivateOneBPChildMeshes = aParamList->get<bool>("cleanup_cut_mesh");
        // mMeshCleanupParameters.mMethod                     = aParamList->get<std::string>("method");
        // mMeshCleanupParameters.mSnapTolerance              = aParamList->get<moris::real>("snap_tolerance");
        
    }

    void
    Mesh_Cleanup::perform()
    {
        // remove inactive child meshes ones that do not have an interior interface
        if(mMeshCleanupParameters.mDeactivateOneBPChildMeshes)
        {
            this->cleanup_cut_mesh();
        }
        // // get the mesh quality data
        // Mesh_Quality_Data tMeshQualityData;
        // this->compute_mesh_quality(tMeshQualityData);
    }

    void
    Mesh_Cleanup::compute_mesh_quality(Mesh_Quality_Data & aMeshQuality)
    {
        // compute volumes
        this->compute_mesh_volumes(aMeshQuality.mVolumes);

        moris::print(aMeshQuality.mVolumes,"aMeshQuality.mVolumes");

    }

    void
    Mesh_Cleanup::compute_mesh_volumes(moris::Matrix<moris::DDRMat> & aVolumes)
    {
        // number of cells in background mesh
        uint tNumCells = mModel->get_background_mesh().get_num_entities(EntityRank::ELEMENT);

        // allocate volumes
        aVolumes.resize( tNumCells , 1 );

        // iterate through background cells and compute the volume
        for(moris::uint i = 0; i<tNumCells; i++)
        {
            // cell
            moris::mtk::Cell & tCell = mModel->get_background_mesh().get_mtk_cell((moris_index) i);
            
            // compute the volume
            aVolumes(i) = tCell.compute_cell_measure();
        }
    }

    void 
    Mesh_Cleanup::compute_mesh_side_lengths(moris::Cell<moris::Matrix<moris::DDRMat>> & aEdgeLengths)
    {
        // number of cells in background mesh
        uint tNumCells = mModel->get_background_mesh().get_num_entities(EntityRank::ELEMENT);
        
        // allocate volumes
        aEdgeLengths.resize(tNumCells );

        // iterate through background cells and compute the volume
        for(moris::uint i = 0; i<tNumCells; i++)
        {
            // cell
            moris::mtk::Cell & tCell = mModel->get_background_mesh().get_mtk_cell((moris_index) i);

            std::cout<<"tCell->get_num_facets() = "<<tCell.get_number_of_facets();

            // allocate lengths
            aEdgeLengths(tCell.get_index()).resize(tCell.get_number_of_facets(),1);
            
            // compute the volume
        //     tVolumes(i) = tCell.compute_cell_measure();
        }
        // return tVolumes;
    }


    void
    Mesh_Cleanup::cleanup_cut_mesh()
    {   
        Tracer tTracer( "XTK", "Mesh Cleanup", "Remove Inactive Child Mesh" );
        moris::uint tNumCutMeshes = mModel->get_cut_mesh().get_num_child_meshes();
        moris::Cell<moris::uint> tChildMeshesToKeep;
        moris::Cell<moris::uint> tChildMeshesToDelete;

        for(moris::uint iCM = 0; iCM < tNumCutMeshes; iCM++)
        {
            Child_Mesh & tCM = mModel->get_cut_mesh().get_child_mesh(iCM);
            Cell<moris::moris_index> const & tSubphasebinBulkPhase = tCM.get_subphase_bin_bulk_phase();

           if(tSubphasebinBulkPhase.size() > 1 or tCM.has_inter_child_mesh_interfaces())
           {
                tChildMeshesToKeep.push_back(iCM);
           }
           else
           {
                tChildMeshesToDelete.push_back(iCM);
           }
        }

        mModel->get_cut_mesh().remove_all_child_meshes_but_selected(tChildMeshesToKeep,tChildMeshesToDelete);

        mModel->get_background_mesh().setup_downward_inheritance(mModel->get_cut_mesh());

        MORIS_LOG_SPEC("Num Child Meshes Removed",tChildMeshesToDelete.size());
        MORIS_LOG_SPEC("Num Child Meshes Kept",tChildMeshesToKeep.size());
    }
}