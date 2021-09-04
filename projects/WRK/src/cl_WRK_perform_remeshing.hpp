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
                moris::Cell< std::string > mRefinementsFieldNames_0;
                Cell< Matrix< DDSMat > >   mRefinementsMode_0;
                Cell< Matrix< DDSMat > >   mRefinementPatternMode_0;
                std::string                mRefinementFunction;

                moris::Cell< moris::Cell < std::string > > mRefinementsFieldNames_1;
                Matrix< DDSMat >                           mMaxRefinementsMode_1;
                Matrix< DDSMat >                           mMinRefinementsMode_1;
                Matrix< DDSMat >                           mRefinementPatternMode_1;
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
                        moris::Cell< std::shared_ptr< mtk::Field > >          aSourceFields,
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
                        std::shared_ptr< hmr::HMR >             aHMRPerformer,
                        Cell< std::shared_ptr< mtk::Field > > & aSourceFields );

                //------------------------------------------------------------------------------

                void perform_refinement_mode_1(
                        std::shared_ptr< hmr::HMR >           aHMRPerformer,
                        Cell< std::shared_ptr< mtk::Field > > aSourceFields);

                //------------------------------------------------------------------------------

                void map_fields(
                        Cell< std::shared_ptr< mtk::Field > > & aSourceFields,
                        Cell< std::shared_ptr< mtk::Field > > & aTargetFields,
                        mtk::Mesh_Pair                        & aMeshPair,
                        uint                                    aDiscretizationMeshIndex,
                        bool                                    aMapFields);

                //------------------------------------------------------------------------------

                void prepare_input_for_refinement(
                        Cell< moris_index >                       & aPatternForRefinement,
                        moris::Cell< moris::Cell< std::string > > & aFieldsForRefinement,
                        moris::Cell< moris::Cell< uint > >        & aRefinements,
                        moris::Cell< sint >                       & aMaxRefinementPerPattern );

                //------------------------------------------------------------------------------

                void create_refinement_input_list(
                        moris::ParameterList & aRefinementParameterlist,
                        uint                   aPattern );

                //------------------------------------------------------------------------------
        };
    }
}

#endif //MORIS_FN_WRK_PERFORM_REFINEMENT_HPP
