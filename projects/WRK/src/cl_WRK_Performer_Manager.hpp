/*
 * cl_WRK_Performer_Manager.hpp
 *
 *  Created on: Feb 19, 2020
 *      Author: schmidt
 */

#ifndef PROJECTS_FEM_MDL_SRC_CL_WRK_PERFORMER_MANAGER_HPP_
#define PROJECTS_FEM_MDL_SRC_CL_WRK_PERFORMER_MANAGER_HPP_

#include "typedefs.hpp"                       //MRS/COR/src
#include "cl_Cell.hpp"                        //MRS/CON/src

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_MTK_Enums.hpp"

namespace moris
{
class Library_IO;
//------------------------------------------------------------------------------
    namespace mtk
    {
       class Mesh_Manager;
    }

    namespace vis
    {
       class Output_Manager;
    }

    namespace fem
    {
        class IWG;
        class Node_Base;
        class Cell;
        class Set;
        class Field_Interpolator;
        class Set_User_Info;
    }

    namespace dla
    {
        class Linear_Solver;
        class Linear_Solver_Algorithm;
    }

    namespace MSI
    {
        class Model_Solver_Interface;
        class MSI_Solver_Interface;
        class Equation_Set;
        class Equation_Object;
        class Equation_Model;
        class Design_Variable_Interface;
        enum class Dof_Type;
    }
    namespace tsa
    {
        class Time_Solver;
        class Time_Solver_Algorithm;
    }
    namespace sol
    {
        class SOL_Warehouse;
    }
    namespace wrk
    {
//------------------------------------------------------------------------------

        class Performer_Manager
        {
            std::string mInputFile;

            std::shared_ptr< Library_IO > mLibrary;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------
            /**
             * constructor
             * @param[ in ] aMesh          mesh for this problem
             */
            Performer_Manager( const std::string & aInputFilePath );

//------------------------------------------------------------------------------
            /**
             * destructor
             */
            ~Performer_Manager();

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace mdl */
} /* namespace moris */


#endif /* PROJECTS_FEM_MDL_SRC_CL_WRK_PERFORMER_MANAGER_HPP_ */
