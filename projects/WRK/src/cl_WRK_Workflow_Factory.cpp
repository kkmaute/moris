/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_WRK_Workflow_Factory.cpp
 *
 */

#include "cl_WRK_Workflow_Factory.hpp"
#include "cl_WRK_Workflow_HMR_XTK.hpp"
#include "cl_WRK_Workflow_STK_XTK.hpp"
#include "cl_WRK_Workflow_STK_FEM.hpp"
namespace moris
{
    namespace wrk
    {
        std::shared_ptr<wrk::Workflow>
        create_workflow(std::string const      & aWRKFlowType,
                        wrk::Performer_Manager * aPerformerManager)
        {
            if(aWRKFlowType.compare("HMR_XTK") == 0)
            {
                return std::make_shared<Workflow_HMR_XTK>( aPerformerManager );
            }
            else if(aWRKFlowType.compare("STK_XTK") == 0)
            {
                return std::make_shared<Workflow_STK_XTK>( aPerformerManager );
            }
            else if(aWRKFlowType.compare("STK_FEM") == 0)
            {
                return std::make_shared<Workflow_STK_FEM>( aPerformerManager );
            }
            else
            {
                MORIS_ERROR(0,"Invalid Workflow Type");
                return nullptr;
            }
        }
    }
}

