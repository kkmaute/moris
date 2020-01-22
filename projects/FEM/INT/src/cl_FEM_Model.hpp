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

        class FEM_Model
        {
            // pointer to reference mesh
            mtk::Mesh_Manager* mMeshManager = nullptr;
            moris_index        mMeshPairIndex;

            // list of IP node pointers
            moris::Cell< fem::Node_Base* > mIPNodes;

            // list of IG node pointers
             moris::Cell< fem::Node_Base* > mIGNodes;

            // list of IP cell pointers
            moris::Cell< fem::Cell* > mIPCells;

            // list of IG cell pointers
            moris::Cell< fem::Cell* > mIGCells;

            // list of FEM sets
            moris::Cell< MSI::Equation_Set * > mFemSets;

            // list of FEM clusters
            moris::Cell< MSI::Equation_Object* > mFemClusters;

            map< moris_index, moris_index >   mMeshSetToFemSetMap;

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
            /**
             * get equation sets for test
             */
            moris::Cell< MSI::Equation_Set * > & get_equation_sets( )
            {
                return mFemSets;
            };

//------------------------------------------------------------------------------
            /**
             * get equation objects
             */
            moris::Cell< MSI::Equation_Object * > & get_equation_objects( )
            {
                return mFemClusters;
            };

//------------------------------------------------------------------------------
            /**
             * MTK set to fem set index map
             */
            map< moris_index, moris_index > & get_mesh_set_to_fem_set_index_map( )
            {
                return mMeshSetToFemSetMap;
            };

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace mdl */
} /* namespace moris */


#endif /* PROJECTS_FEM_MDL_SRC_CL_MDL_MODEL_HPP_ */
