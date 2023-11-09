/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_WRK_perform_refinement.hpp
 *
 */

#pragma once

#include "cl_Matrix.hpp"

namespace moris
{
    //--------------------------------------------------------------------------------------------------------------

    namespace hmr
    {
        class HMR;
        class Mesh;
        class Element;

        // User-defined refinement function
        // typedef sint  ( *Refinement_Function ) (
        //        hmr::Element           * aElement,
        //        const Matrix< DDRMat > & aElementLocalValues);
    }    // namespace hmr

    //--------------------------------------------------------------------------------------------------------------

    namespace mtk
    {
        class Field;
    }

    namespace wrk
    {
        class Performer;

        //--------------------------------------------------------------------------------------------------------------

        struct Refinement_Parameters
        {
            //! Field names
            Cell< std::string > mFieldNames;

            //! Refinement level
            Cell< Matrix< DDSMat > > mRefinementLevel;

            //! Refinement Pattern
            Cell< Matrix< DDSMat > > mRefinementPattern;

            Cell< Matrix< DDSMat > > mRefinemenCopytPatternToPattern;

            hmr::Refinement_Function mRefinementFunction = nullptr;

            moris::Cell< std::string > mRefinementFunctionName;

            moris::Cell< hmr::Refinement_Function_2 > mRefinementFunction_2;
        };

        //--------------------------------------------------------------------------------------------------------------

        class Refinement_Mini_Performer
        {
            //--------------------------------------------------------------------------------------------------------------

          private:
            Refinement_Parameters mParameters;

            std::shared_ptr< Library_IO > mLibrary = nullptr;

            //--------------------------------------------------------------------------------------------------------------

          public:
            //--------------------------------------------------------------------------------------------------------------

            Refinement_Mini_Performer(){};

            //--------------------------------------------------------------------------------------------------------------

            Refinement_Mini_Performer(
                    ParameterList                &aParameterlist,
                    std::shared_ptr< Library_IO > aLibrary = nullptr );

            //--------------------------------------------------------------------------------------------------------------

            ~Refinement_Mini_Performer(){};

            //--------------------------------------------------------------------------------------------------------------

            void perform_refinement(
                    Cell< std::shared_ptr< mtk::Field > > &aFields,
                    std::shared_ptr< hmr::HMR >            aHMR );

            //--------------------------------------------------------------------------------------------------------------

            uint perform_refinement_low_level_elements(
                    Cell< std::shared_ptr< mtk::Field > > &aFields,
                    std::shared_ptr< hmr::HMR >            aHMR );

            //--------------------------------------------------------------------------------------------------------------

            void perform_refinement_based_on_working_pattern(
                    Cell< std::shared_ptr< mtk::Field > > &aFields,
                    std::shared_ptr< hmr::HMR >            aHMR );

            void perform_refinement_2(
                    Cell< std::shared_ptr< mtk::Field > > &aFields,
                    std::shared_ptr< hmr::HMR >            aHMR );

            //--------------------------------------------------------------------------------------------------------------

            /**
             * This function reorganizes the input data. The input data is structured such that it allows for an intuitive use by the user.
             * However, the structure provided by this function allows for more straight forward use in the refinement function
             *
             * @param aPatternForRefinement Pattern which will be refined
             * @param aFieldsForRefinement  Fields which are used to refined each pattern
             * @param aRefinements          Refinements for each field per pattern
             */
            void prepare_input_for_refinement(
                    Cell< moris_index >                                      &aPatternForRefinement,
                    moris::Cell< moris::Cell< std::string > >                &aFieldsForRefinement,
                    moris::Cell< moris::Cell< uint > >                       &aRefinements,
                    moris::Cell< sint >                                      &aMaxRefinementPerPattern,
                    moris::Cell< moris::Cell< hmr::Refinement_Function_2 > > &aRefinementFunctions );

            //--------------------------------------------------------------------------------------------------------------
            // FIXME stuff below this line will be deleted soon
            //----------------------------------------------------------------------------------------------------------------

            /**
             * Performs refinement using HMR based on the information provided by the cell of performers
             *
             * @param aHMR Active instance of HMR.
             * @param aPerformers Performers which can provide refinement data.
             * @param aSimultaneous If true (default), refinement steps are generated using all performers simultaneously.
             * If false, one performer will have all refinement steps performed before moving to the next performer.
             */
            void perform_refinement_old(
                    std::shared_ptr< hmr::HMR >          aHMR,
                    Cell< std::shared_ptr< Performer > > aPerformers,
                    bool                                 aSimultaneous = true );

            //--------------------------------------------------------------------------------------------------------------

            /**
             * Queues a single refinement step. This is only intended to be called by \ref moris::wrk::perform_refinement()
             *
             * @param aHMR Active instance of HMR.
             * @param aMesh Mesh generated by HMR, can potentially have some refinement performed on it already.
             * @param aPerformer A performer which can provide refinement data.
             * @param aRefinementNumber The refinement number being queued.
             * @return
             */
            void queue_single_refinement( std::shared_ptr< hmr::HMR > aHMR,
                    std::shared_ptr< hmr::Mesh >                      aMesh,
                    std::shared_ptr< Performer >                      aPerformer,
                    uint                                              aRefinementNumber,
                    uint                                              aMeshIndex );

            //--------------------------------------------------------------------------------------------------------------

            uint queue_low_level_elements_for_refinement(
                    std::shared_ptr< hmr::HMR >  aHMR,
                    std::shared_ptr< hmr::Mesh > aMesh,
                    std::shared_ptr< Performer > aPerformer,
                    uint                         aMeshIndex );

            //--------------------------------------------------------------------------------------------------------------

            uint get_max_refinement_level( const Cell< std::shared_ptr< Performer > > &aPerformers );

            //--------------------------------------------------------------------------------------------------------------

            void get_all_refinement_mesh_indices(
                    const Cell< std::shared_ptr< Performer > > &aPerformers,
                    moris::Matrix< DDSMat >                    &aAllPatternMap,
                    moris::uint                                &aNumPattern );

            //--------------------------------------------------------------------------------------------------------------

        };    // end class: Refinement_Mini_Performer

        //--------------------------------------------------------------------------------------------------------------

    }    // namespace wrk
}    // namespace moris
