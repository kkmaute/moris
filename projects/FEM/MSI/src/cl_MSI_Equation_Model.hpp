/*
 * cl_MSI_Model.hpp
 *
 *  Created on: Aug 22, 2018
 *      Author: schmidt
 */

#ifndef PROJECTS_FEM_MDL_SRC_CL_MSI_MODEL_HPP_
#define PROJECTS_FEM_MDL_SRC_CL_MSIMODEL_HPP_

#include "typedefs.hpp"                       //MRS/COR/src
#include "cl_Cell.hpp"                        //MRS/CON/src

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

namespace moris
{

//------------------------------------------------------------------------------

    namespace MSI
    {
        class Model_Solver_Interface;
        class MSI_Solver_Interface;
        class Equation_Set;
        class Equation_Object;
        class Design_Variable_Interface;
        enum class Dof_Type;
//------------------------------------------------------------------------------

        class Equation_Model
        {
        protected:
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
             */
            Equation_Model(){};

//------------------------------------------------------------------------------
            /**
             * destructor
             */
            virtual ~Equation_Model(){};

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
    } /* namespace MSI */
} /* namespace moris */


#endif /* PROJECTS_FEM_MDL_SRC_CL_MSI_MODEL_HPP_ */
