/*
 * cl_FEM_Model.hpp
 *
 *  Created on: Aug 22, 2018
 *      Author: schmidt
 */

#ifndef PROJECTS_FEM_MDL_SRC_CL_FEM_MODEL_HPP_
#define PROJECTS_FEM_MDL_SRC_CL_FEM_MODEL_HPP_

#include "typedefs.hpp"                       //MRS/COR/src
#include "cl_Cell.hpp"                        //MRS/CON/src

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_MTK_Enums.hpp"

#include "cl_MSI_Equation_Model.hpp"

namespace moris
{

//------------------------------------------------------------------------------
    namespace mtk
    {
       class Mesh_Manager;
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

    namespace MSI
    {
        class Model_Solver_Interface;
        class MSI_Solver_Interface;
        class Equation_Set;
        class Equation_Object;
        class Design_Variable_Interface;
        enum class Dof_Type;
    }

    namespace fem
    {
//------------------------------------------------------------------------------

        class FEM_Model : public  MSI::Equation_Model
        {
            // pointer to reference mesh
            mtk::Mesh_Manager* mMeshManager = nullptr;
            moris_index        mMeshPairIndex;

            // list of IP node pointers
            moris::Cell< fem::Node_Base* > mIPNodes;

            // list of IP cell pointers
            moris::Cell< fem::Cell* > mIPCells;

            moris::Cell< moris::real >       mQi;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------
            /**
             * constructor
             * @param[ in ] aMesh          mesh for this problem
             * @param[ in ] aMeshPairIndex mesh pair index
             * @param[ in ] aSetInfo       cell of set user info
             */
            FEM_Model(       mtk::Mesh_Manager                 * aMeshManager,
                       const moris_index                       & aMeshPairIndex,
                             moris::Cell< fem::Set_User_Info > & aSetInfo );

//------------------------------------------------------------------------------
            /**
             * destructor
             */
            ~FEM_Model();

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace mdl */
} /* namespace moris */


#endif /* PROJECTS_FEM_MDL_SRC_CL_MDL_MODEL_HPP_ */
