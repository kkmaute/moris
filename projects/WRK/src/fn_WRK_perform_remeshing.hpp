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
    namespace mtk
    {
        class Field;
    }


    namespace wrk
    {
        class Performer;

        void perform_remeshing(
                mtk::Field                                 * aInputField,
                moris::Cell< std::shared_ptr< hmr::HMR > >   aHMRPerformers );
    }
}

#endif //MORIS_FN_WRK_PERFORM_REFINEMENT_HPP
