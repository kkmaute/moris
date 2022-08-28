/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_WRK_Workflow_STK_XTK.hpp
 *
 */

#ifndef PROJECTS_FEM_MDL_SRC_CL_WRK_WORKFLOW_STK_XTK_HPP_
#define PROJECTS_FEM_MDL_SRC_CL_WRK_WORKFLOW_STK_XTK_HPP_

#include "cl_WRK_Workflow.hpp"
#include "typedefs.hpp"                       //MRS/COR/src
#include "cl_Cell.hpp"

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
                ~Workflow_STK_XTK(){};

                //------------------------------------------------------------------------------
                /**
                 * Initializes the vectors of ADV values, lower bounds, and upper bounds
                 */
                void initialize(
                        Matrix<DDRMat>& aADVs,
                        Matrix<DDRMat>& aLowerBounds,
                        Matrix<DDRMat>& aUpperBounds,
                        Matrix< IdMat  >& aIjklIDs);

                //------------------------------------------------------------------------------
                /**
                 * Gets the criteria values given a new set of ADVs
                 *
                 * @return vector of criteria
                 */
                Matrix< DDRMat > perform( Matrix< DDRMat >& aNewADVs );

                //------------------------------------------------------------------------------
                /**
                 * Gets the derivative of the criteria with respect to the advs
                 *
                 * @return matrix d(criteria)_i/d(adv)_j
                 */
                Matrix<DDRMat> compute_dcriteria_dadv();

                void
                create_xtk();

                void
                create_stk(Cell< Cell<ParameterList> > & aParameterLists);

        };
        //------------------------------------------------------------------------------
    } /* namespace mdl */
} /* namespace moris */

#endif /* PROJECTS_FEM_MDL_SRC_CL_WRK_WORKFLOW_STK_XTK_HPP_ */

