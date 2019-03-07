/*
 * cl_FEM_Element_Factory.hpp
 *
 *  Created on: Feb 18, 2019
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_ELEMENT_FACTORY_HPP_
#define SRC_FEM_CL_FEM_ELEMENT_FACTORY_HPP_

#include "assert.h"

#include "typedefs.hpp"               //MRS/COR/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_FEM_Enums.hpp"           //FEM/INT/src
#include "cl_FEM_IWG.hpp"             //FEM/INT/src
#include "cl_FEM_Node.hpp"            //FEM/INT/src
#include "cl_MSI_Equation_Object.hpp" //FEM/MSI/src
#include "cl_MTK_Cell.hpp" //MTK
#include "cl_FEM_Element.hpp"         //FEM/INT/src


namespace moris
{
//------------------------------------------------------------------------------
    namespace fem
    {
//------------------------------------------------------------------------------

    /**
     * \brief element factory
     */
    class Element_Factory
    {

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        /**
         * constructor
         */
        Element_Factory();

//------------------------------------------------------------------------------

        /**
         * trivial destructor
         */
        ~Element_Factory();

//------------------------------------------------------------------------------
        /**
         * create element
         */
        MSI::Equation_Object * create_element( Element_Type         aElementType,
                                               mtk::Cell          * aCell,
                                               Cell< IWG* >       & aIWGs,
                                               Cell< Node_Base* > & aNodes  );

    };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_ELEMENT_FACTORY_HPP_ */
