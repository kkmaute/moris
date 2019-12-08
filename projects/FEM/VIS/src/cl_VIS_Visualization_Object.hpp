/*
 * cl_VIS_Visualization_Object.hpp
 *
 *  Created on: Dez 02, 2019
 *      Author: schmidt
 */
#ifndef SRC_FEM_CL_VIS_VISUALIZATION_OBJECT_HPP_
#define SRC_FEM_CL_VIS_VISUALIZATION_OBJECT_HPP_

#include "cl_Cell.hpp"
#include "cl_Communication_Tools.hpp"
#include "cl_Communication_Manager.hpp"

#include "cl_MSI_Equation_Object.hpp"


namespace moris
{
    namespace vis
    {
        class Visualization_Object
        {
        protected:

            MSI::Equation_Object * mEquationObject =nullptr;

        public:
            Visualization_Object( MSI::Equation_Object * aEquationObject ) : mEquationObject( aEquationObject )
            {

            };

//-----------------------------------------------------------------------------------------------------------

            ~Visualization_Object()
            {

            };

            virtual Matrix< DDRMat > get_element_values() = 0;

            virtual Matrix< IndexMat > get_element_inds() = 0;

//-----------------------------------------------------------------------------------------------------------

//            void create_visualization_nodes()

//-----------------------------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------------------------

        };
    } /* namespace VIS */
} /* namespace moris */

#endif /* SRC_FEM_CL_VIS_VISUALIZATION_OBJECT_HPP_ */
