/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_WRK_Workflow_STK_XTK.hpp
 *
 */

#pragma once

#include "cl_WRK_Workflow.hpp"
#include "moris_typedefs.hpp"                       //MRS/COR/src
#include "cl_Vector.hpp"

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

namespace moris
{
    namespace mtk
    {
        class Interpolation_Mesh;
        class Integration_Mesh;
    }
    namespace wrk
    {
        class Performer_Manager;
        //------------------------------------------------------------------------------
        // Naming convention here means the background mesh is constructed by STK and the integration mesh is constructed by XTK
        class Workflow_STK_XTK : public Workflow
        {
            private:
                std::shared_ptr<mtk::Interpolation_Mesh> mIpMesh;
                std::shared_ptr<mtk::Integration_Mesh>   mIgMesh;
            public:

                //------------------------------------------------------------------------------
                /**
                 * constructor
                 */
                Workflow_STK_XTK( wrk::Performer_Manager * aPerformerManager );

                //------------------------------------------------------------------------------
                /**
                 * destructor
                 */
                ~Workflow_STK_XTK() override{};

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

                void
                create_xtk();

                void
                create_stk( Vector< Vector< Parameter_List > > & aParameterLists);

        };
        //------------------------------------------------------------------------------
    } /* namespace mdl */
} /* namespace moris */
