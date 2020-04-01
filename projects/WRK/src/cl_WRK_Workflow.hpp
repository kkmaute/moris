/*
 * cl_WRK_Workflow.hpp
 *
 *  Created on: Feb 19, 2020
 *      Author: schmidt
 */

#ifndef PROJECTS_FEM_MDL_SRC_CL_WRK_WORKFLOW_HPP_
#define PROJECTS_FEM_MDL_SRC_CL_WRK_WORKFLOW_HPP_

#include "typedefs.hpp"                       //MRS/COR/src
#include "cl_Cell.hpp"                        //MRS/CON/src

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
//------------------------------------------------------------------------------

        class Performer_Manager
        {
            std::shared_ptr< Library_IO > mLibrary;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------
            /**
             * constructor
             * @param[ in ] aMesh          mesh for this problem
             */
            Performer_Manager( std::shared_ptr< Library_IO > aLibrary );

//------------------------------------------------------------------------------
            /**
             * destructor
             */
            ~Performer_Manager(){};

            void initialize();

            void perform_refinement();

            void perform();

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace mdl */
} /* namespace moris */


#endif /* PROJECTS_FEM_MDL_SRC_CL_WRK_PERFORMER_MANAGER_HPP_ */
