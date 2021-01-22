/*
 * cl_WRK_Workflow_Factory.hpp
 *
 *  Created on: Jan 22, 2021
 *      Author: doble
 */
#ifndef PROJECTS_FEM_MDL_SRC_CL_WRK_WORKFLOW_FACTORY_HPP_
#define PROJECTS_FEM_MDL_SRC_CL_WRK_WORKFLOW_FACTORY_HPP_

#include <memory>
#include "cl_WRK_Workflow.hpp"
#include "cl_WRK_Performer_Manager.hpp"
namespace moris
{
    namespace wrk
    {   
        /*!
        *  Construct the desired workflow
        *  @param[in] aWRKFlowType  String which specifies the desired workflow
        * 
        */
        std::shared_ptr<wrk::Workflow>
        create_workflow(std::string const & aWRKFlowType,
                        wrk::Performer_Manager * aPerformerManager);
    }
}
#endif /*PROJECTS_FEM_MDL_SRC_CL_WRK_WORKFLOW_FACTORY_HPP_*/
