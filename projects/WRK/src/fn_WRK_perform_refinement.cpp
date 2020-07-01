#include "fn_WRK_perform_refinement.hpp"
#include "cl_WRK_Performer.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Mesh.hpp"

namespace moris
{
    namespace wrk
    {

        //--------------------------------------------------------------------------------------------------------------

        void perform_refinement(std::shared_ptr<hmr::HMR>          aHMR,
                                Cell< std::shared_ptr<Performer> > aPerformers,
                                bool                               aSimultaneous)
        {
            // Create mesh
            std::shared_ptr<hmr::Mesh> tMesh = aHMR->create_mesh(0);

            // Set refinement index/flag
            bool tPerformRefinement = true;
            sint tRefinementNumber = 0;

            // Set loop bounds
            uint tFirstPerformer = 0;
            uint tLastPerformer = (aSimultaneous ? aPerformers.size() : 0);

            // Loop over set number of refinement levels
            while (tPerformRefinement)
            {
                // Reset flag
                tPerformRefinement = false;

                for (uint tPerformerIndex = tFirstPerformer; tPerformerIndex < tLastPerformer; tPerformerIndex++)
                {
                    // Queue refinement
                    tPerformRefinement = (tPerformRefinement or
                            queue_single_refinement(aHMR, tMesh, aPerformers(tPerformerIndex), tRefinementNumber));
                }

                // Perform refinement and update index
                if (tPerformRefinement)
                {
                    aHMR->perform_refinement_based_on_working_pattern( 0, false );
                }

                // Update performers
                tFirstPerformer++;
                if (aSimultaneous or tFirstPerformer == aPerformers.size())
                {
                    tRefinementNumber++;
                    tFirstPerformer = 0;
                }
                if (!aSimultaneous)
                {
                    tLastPerformer = tFirstPerformer;
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        bool queue_single_refinement(std::shared_ptr<hmr::HMR>  aHMR,
                                     std::shared_ptr<hmr::Mesh> aMesh,
                                     std::shared_ptr<Performer> aPerformer,
                                     uint                       aRefinementNumber)
        {
            // Set refinement flag
            bool tPerformRefinement = false;

            // Loop over fields
            for (uint tFieldIndex = 0; tFieldIndex < aPerformer->get_num_refinement_fields(); tFieldIndex++)
            {
                // Determine if refinement is needed
                if (aPerformer->refinement_needed(tFieldIndex, aRefinementNumber))
                {
                    // Set flag
                    tPerformRefinement = true;

                    // Loop over nodes and get field values
                    Matrix<DDRMat> tFieldValues(aMesh->get_num_nodes(), 1);
                    for (uint tNodeIndex = 0; tNodeIndex < aMesh->get_num_nodes(); tNodeIndex++)
                    {
                        tFieldValues(tNodeIndex) = aPerformer->get_field_value(tFieldIndex, tNodeIndex, aMesh->get_node_coordinate(tNodeIndex));
                    }

                    // Put elements on queue and set flag for refinement
                    aHMR->based_on_field_put_elements_on_queue(tFieldValues, 0, aPerformer->get_refinement_function_index(tFieldIndex, aRefinementNumber));
                }
            }

            return tPerformRefinement;
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
