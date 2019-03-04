/*
 * cl_FEM_Element_HJ.hpp
 *
 *  Created on: Feb 18, 2019
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_ELEMENT_HJ_HPP_
#define SRC_FEM_CL_FEM_ELEMENT_HJ_HPP_

#include "assert.h"

#include "typedefs.hpp"           //MRS/COR/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_FEM_Element.hpp" //FEM/INT/src

namespace moris
{

//------------------------------------------------------------------------------

    namespace fem
    {
//------------------------------------------------------------------------------

    /**
     * \brief Element_HamiltonJacobi class that inherits from the Element class
     */
    class Element_HJ : public Element
    {

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------
        /**
         * constructor
         */
        Element_HJ( Cell< IWG * >         aListOfIWGs,
                    Cell< MSI::Dof_Type > aListOfNodePDOFsTypes,
                    Cell< MSI::Dof_Type > aListOfElementPDOFsTypes );

//------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~Element_HJ(){};

//------------------------------------------------------------------------------
        /**
         * compute the residual
         */
        void compute_residual();

//------------------------------------------------------------------------------
        /**
         * compute the jacobian
         */
        void compute_jacobian();

//------------------------------------------------------------------------------
        /**
         * compute the jacobian and the residual
         */
        void compute_jacobian_and_residual();

//------------------------------------------------------------------------------
    };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */




#endif /* SRC_FEM_CL_FEM_ELEMENT_HJ_HPP_ */
