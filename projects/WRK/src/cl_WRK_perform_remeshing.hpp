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
        struct Remeshing_Parameters
        {
                //! Mode name. Mainly for output purposes
                std::string mMode;

                //! Mode index. Transleted based in mode name
                uint        mModeIndex;
        };

        class Remeshing_Mini_Performer
        {
            private:
                wrk::Remeshing_Parameters mParameters;

                std::shared_ptr< Library_IO > mLibrary = nullptr;

            public:

                Remeshing_Mini_Performer(
                        ParameterList                 & aParameterlist,
                        std::shared_ptr< Library_IO >   aLibrary );

                //------------------------------------------------------------------------------

                ~Remeshing_Mini_Performer(){};

                //------------------------------------------------------------------------------

                void perform_remeshing(
                        mtk::Field                                 * aInputField,
                        moris::Cell< std::shared_ptr< hmr::HMR > > & aHMRPerformers );

                //------------------------------------------------------------------------------

                void perform_refinement(
                        std::shared_ptr< hmr::HMR >   aHMRPerformer,
                        mtk::Field                  * aSourceField );

                //------------------------------------------------------------------------------

                void perform_refinement_mode_0(
                        std::shared_ptr< hmr::HMR >   aHMRPerformer,
                        mtk::Field                  * aSourceField );

                //------------------------------------------------------------------------------

                void perform_refinement_mode_1(
                        std::shared_ptr< hmr::HMR >   aHMRPerformer,
                        mtk::Field                  * aSourceField );

                //------------------------------------------------------------------------------
        };
    }
}

#endif //MORIS_FN_WRK_PERFORM_REFINEMENT_HPP
