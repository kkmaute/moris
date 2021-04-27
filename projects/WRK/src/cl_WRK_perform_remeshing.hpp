#ifndef MORIS_FN_WRK_PERFORM_REFINEMENT_HPP
#define MORIS_FN_WRK_PERFORM_REMESHING_HPP

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

                //mode ab_initio
                Matrix< DDUMat > mRefinementsMode_0;
                Matrix< DDUMat > mRefinementPatternMode_0;
        };

        class Remeshing_Mini_Performer
        {
            private:
                wrk::Remeshing_Parameters mParameters;

                std::shared_ptr< Library_IO > mLibrary = nullptr;

            public:

                Remeshing_Mini_Performer(
                        ParameterList                 & aParameterlist,
                        std::shared_ptr< Library_IO >   aLibrary = nullptr );

                //------------------------------------------------------------------------------

                ~Remeshing_Mini_Performer(){};

                //------------------------------------------------------------------------------

                void perform_remeshing(
                        moris::Cell< std::shared_ptr< mtk::Field > > const  & aSourceFields,
                        moris::Cell< std::shared_ptr< hmr::HMR > >          & aHMRPerformers,
                        moris::Cell< std::shared_ptr< mtk::Mesh_Manager > > & mMTKPerformer,
                        moris::Cell< std::shared_ptr< mtk::Field > >        & aNewFields);

                //------------------------------------------------------------------------------

                void perform_refinement(
                        std::shared_ptr< hmr::HMR >           aHMRPerformer,
                        Cell< std::shared_ptr< mtk::Field > > aSourceFields );

                //------------------------------------------------------------------------------

                /**
                 * Refinement before creating the new mesh.
                 * This mode assumes a uniformly, refined starting mesh for all used pattern.
                 * Pattern are refined separately up to the specific level.
                 *
                 */
                void perform_refinement_mode_0(
                        std::shared_ptr< hmr::HMR >           aHMRPerformer,
                        Cell< std::shared_ptr< mtk::Field > > aSourceFields );

                //------------------------------------------------------------------------------

                void perform_refinement_mode_1(
                        std::shared_ptr< hmr::HMR >           aHMRPerformer,
                        Cell< std::shared_ptr< mtk::Field > > aSourceFields);

                //------------------------------------------------------------------------------

                void map_field(
                std::shared_ptr< hmr::HMR >     aHMRPerformer,
                mtk::Field                  *   aSourceField,
                mtk::Field                  * & aTargetField );

                //------------------------------------------------------------------------------
        };
    }
}

#endif //MORIS_FN_WRK_PERFORM_REFINEMENT_HPP
