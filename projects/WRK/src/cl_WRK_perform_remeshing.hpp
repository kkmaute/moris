/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_WRK_perform_remeshing.hpp
 *
 */

#pragma once

#include "cl_Matrix.hpp"

namespace moris
{
    namespace hmr
    {
        class HMR;
        class Mesh;
    }    // namespace hmr
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

            //! Mode index. Translated based in mode name
            uint mModeIndex;

            //! Remeshing frequency
            sint mRemeshingFrequency;

            bool mOutputMeshes = false;

            // mode ab_initio
            Vector< Matrix< DDSMat > > mRefinementsMode_0;
            Vector< Matrix< DDSMat > > mRefinementPatternMode_0;
            std::string                mRefinementFunction;

            Matrix< DDSMat > mMaxRefinementsMode_1;
            Matrix< DDSMat > mMinRefinementsMode_1;
            Matrix< DDSMat > mRefinementPatternMode_1;

            Vector< Matrix< DDSMat > > mRefinementPatternMode_3;
            Vector< Matrix< DDSMat > > mRefinemenCopytPatternToPattern_3;
            Vector< Matrix< DDSMat > > mMinimumRefinementLevel;
            Vector< std::string >      mRefinementFunctionForField;
        };

        class Remeshing_Mini_Performer
        {
          private:
            wrk::Remeshing_Parameters mParameters;

            std::shared_ptr< Library_IO > mLibrary = nullptr;

          public:
            Remeshing_Mini_Performer(
                    ParameterList&                aParameterlist,
                    std::shared_ptr< Library_IO > aLibrary = nullptr );

            //------------------------------------------------------------------------------

            ~Remeshing_Mini_Performer(){};

            //------------------------------------------------------------------------------

            sint
            get_reinitialization_frequency()
            {
                return mParameters.mRemeshingFrequency;
            }

            //------------------------------------------------------------------------------

            void perform_remeshing(
                    Vector< std::shared_ptr< mtk::Field > >         aSourceFields,
                    Vector< std::shared_ptr< hmr::HMR > >&          aHMRPerformers,
                    Vector< std::shared_ptr< mtk::Mesh_Manager > >& mMTKPerformer,
                    Vector< std::shared_ptr< mtk::Field > >&        aNewFields );

            //------------------------------------------------------------------------------

            void perform_refinement(
                    std::shared_ptr< hmr::HMR >             aHMRPerformer,
                    Vector< std::shared_ptr< mtk::Field > > aSourceFields );

            //------------------------------------------------------------------------------

            /**
             * Refinement before creating the new mesh.
             * This mode assumes a uniformly, refined starting mesh for all used pattern.
             * Pattern are refined separately up to the specific level.
             *
             */
            void perform_refinement_mode_0(
                    std::shared_ptr< hmr::HMR >              aHMRPerformer,
                    Vector< std::shared_ptr< mtk::Field > >& aSourceFields );

            //------------------------------------------------------------------------------

            void perform_refinement_mode_1(
                    std::shared_ptr< hmr::HMR >             aHMRPerformer,
                    Vector< std::shared_ptr< mtk::Field > > aSourceFields );

            //------------------------------------------------------------------------------

            void perform_refinement_mode_2(
                    std::shared_ptr< hmr::HMR >             aHMRPerformer,
                    Vector< std::shared_ptr< mtk::Field > > aSourceFields );

            //------------------------------------------------------------------------------

            void map_fields(
                    Vector< std::shared_ptr< mtk::Field > >& aSourceFields,
                    Vector< std::shared_ptr< mtk::Field > >& aTargetFields,
                    mtk::Mesh_Pair&                          aMeshPair,
                    uint                                     aDiscretizationMeshIndex,
                    bool                                     aMapFields );

            //------------------------------------------------------------------------------

            void prepare_input_for_refinement(
                    Vector< moris_index >&    aPatternForRefinement,
                    Vector< Vector< uint > >& aRefinements,
                    Vector< sint >&           aMaxRefinementPerPattern );

            //------------------------------------------------------------------------------

            void create_refinement_input_list(
                    moris::ParameterList&                  aRefinementParameterlist,
                    Cell< std::shared_ptr< mtk::Field > >& aFields,
                    uint                                   aPattern );

            //------------------------------------------------------------------------------

            void create_refinement_input_list_2(
                    moris::ParameterList&                  aRefinementParameterlist,
                    Cell< std::shared_ptr< mtk::Field > >& aFields );

            //------------------------------------------------------------------------------

            void output_meshes( std::shared_ptr< hmr::HMR > aHMRPerformer );

            //------------------------------------------------------------------------------

            void unite_all_pattern_for_lagrange( std::shared_ptr< hmr::HMR > aHMRPerformer );

            //------------------------------------------------------------------------------
        };
    }    // namespace wrk
}    // namespace moris
