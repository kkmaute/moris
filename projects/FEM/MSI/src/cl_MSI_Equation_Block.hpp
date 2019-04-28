/*
 * cl_MSI_Equation_Block.hpp
 *
 *  Created on: Apr 10, 2019
 *      Author: schmidt
 */

#ifndef SRC_FEM_CL_MSI_EQUATION_BLOCK_HPP_
#define SRC_FEM_CL_MSI_EQUATION_BLOCK_HPP_

#include "assert.h"
#include "cl_Communication_Tools.hpp"               //FEM/INT/src

namespace moris
{
namespace mtk
{
   class Cell;
}
    namespace MSI
    {
    class Model_Solver_Interface;
    class Equation_Object;
    enum class Dof_Type;
//------------------------------------------------------------------------------
    /**
     * \brief element block class that communicates with the mesh interface
     */
    class Equation_Block
    {
    protected:
        Cell< MSI::Equation_Object* > mElements;

        // map of the element active dof types
        moris::Cell< enum MSI::Dof_Type > mEqnObjDofTypeList; // List of dof types of this equation obj

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        /**
         * trivial constructor
         */
        Equation_Block( ){};

//------------------------------------------------------------------------------

        /**
         * trivial destructor
         */
        virtual ~Equation_Block(){};

//------------------------------------------------------------------------------

//        void delete_pointers();

//------------------------------------------------------------------------------

        virtual void finalize( MSI::Model_Solver_Interface * aModelSolverInterface )
        {
            MORIS_ERROR(false,"Equation_Block::finalize(), not implemented");
        };

//------------------------------------------------------------------------------

        Cell< MSI::Equation_Object * > & get_equation_object_list()
        {
            return mElements;
        };

//------------------------------------------------------------------------------

        moris::Cell< enum MSI::Dof_Type > & get_unique_dof_type_list()
        {
            return mEqnObjDofTypeList;
        }
    };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_MSI_EQUATION_BLOCK_HPP_ */
