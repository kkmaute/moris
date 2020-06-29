#ifndef MORIS_FN_WRK_PERFORM_REFINEMENT_HPP
#define MORIS_FN_WRK_PERFORM_REFINEMENT_HPP

#include "cl_Matrix.hpp"

namespace moris
{
    namespace hmr
    {
        class HMR;
        class Mesh;
    }

    namespace wrk
    {
        class Performer;

        void perform_refinement(std::shared_ptr<hmr::HMR>          aHMR,
                                Cell< std::shared_ptr<Performer> > aPerformers,
                                bool                               aSimultaneous = true);

        bool queue_single_refinement(std::shared_ptr<hmr::HMR>  aHMR,
                                     std::shared_ptr<hmr::Mesh> aMesh,
                                     std::shared_ptr<Performer> aPerformer,
                                     uint                       aRefinementNumber);
    }
}

#endif //MORIS_FN_WRK_PERFORM_REFINEMENT_HPP
