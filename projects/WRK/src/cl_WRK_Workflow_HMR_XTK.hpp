/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_WRK_Workflow_HMR_XTK.hpp
 *
 */

#pragma once

#include "cl_WRK_Workflow.hpp"
#include "cl_OPT_Criteria_Interface.hpp"
#include "moris_typedefs.hpp"    //MRS/COR/src
#include "cl_Vector.hpp"

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_Communication_Tools.hpp"

namespace xtk
{
    class Model;
}

namespace moris
{
    class Library_IO;
    //------------------------------------------------------------------------------
    namespace hmr
    {
        class HMR;
    }
    namespace mtk
    {
        class Mesh_Manager;
    }
    namespace ge
    {
        class GEN_Geometry_Engine;
    }
    namespace mdl
    {
        class Model;
    }
    namespace opt
    {
        class Manager;
    }

    namespace wrk
    {
        class Performer_Manager;
        //------------------------------------------------------------------------------
        // Naming convention here means the background mesh is constructed by HMR and the integration mesh is constructed by XTK
        class Workflow_HMR_XTK : public Workflow
        {
          private:

          public:
            //------------------------------------------------------------------------------
            /**
             * constructor
             */
            Workflow_HMR_XTK( wrk::Performer_Manager* aPerformerManager );

            //------------------------------------------------------------------------------
            /**
             * destructor
             */
            ~Workflow_HMR_XTK() override{};

            //------------------------------------------------------------------------------
            /**
             * Initializes the vectors of ADV values, lower bounds, and upper bounds
             */
            void initialize(
                    Vector< real >& aADVs,
                    Vector< real >& aLowerBounds,
                    Vector< real >& aUpperBounds,
                    Matrix< IdMat >& aIjklIDs ) override;

            //------------------------------------------------------------------------------
            /**
             * Gets the criteria values given a new set of ADVs
             *
             * @return vector of criteria
             */
            Vector< real > perform( Vector< real >& aNewADVs ) override;

            //------------------------------------------------------------------------------
            /**
             * Gets the derivative of the criteria with respect to the advs
             *
             * @return matrix d(criteria)_i/d(adv)_j
             */
            Matrix< DDRMat > compute_dcriteria_dadv() override;

            //------------------------------------------------------------------------------
            /**
             * @brief Performs T-Matrix and MPC outputs (used for project work) if requested
             * by user.
             *
             * @param aMTKPerformer pointer to MTK Mesh_Manager
             * @param aXTKPerformer pointer to XTK_Model
             */
            bool
            output_T_matrices(
                    const std::shared_ptr< mtk::Mesh_Manager >& aMTKPerformer,
                    xtk::Model* const &                         aXTKPerformer );
        };
        //------------------------------------------------------------------------------
    }    // namespace wrk
} /* namespace moris */
