/*
 * cl_VIS_Visualization_Set.hpp
 *
 *  Created on: Dez 02, 2019
 *      Author: schmidt
 */
#ifndef SRC_FEM_CL_VIS_VISUALIZATION_SET_HPP_
#define SRC_FEM_CL_VIS_VISUALIZATION_SET_HPP_

#include "cl_Cell.hpp"
#include "cl_Communication_Tools.hpp"
#include "cl_Communication_Manager.hpp"

#include "cl_MSI_Equation_Set.hpp"
#include "cl_VIS_Visualization_Object.hpp"

namespace moris
{
    namespace vis
    {
        class Visualization_Set
        {
        protected:
            MSI::Equation_Set * mEquationSet = nullptr;

            Cell< vis::Visualization_Object * > mVisualizationObjects;

        public:
            Visualization_Set(  MSI::Equation_Set * aEquationSet ) : mEquationSet( aEquationSet )
            {
//                uint get_num_equation_objects()
//                Cell< MSI::Equation_Object * > & get_equation_object_list()
            };

//-----------------------------------------------------------------------------------------------------------

            ~Visualization_Set()
            {

            };

//-----------------------------------------------------------------------------------------------------------

            virtual Matrix< DDRMat > calculate_elemental_values() = 0;

//-----------------------------------------------------------------------------------------------------------

            virtual void create_visualization_objects() = 0;

//-----------------------------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------------------------

        };
    } /* namespace VIS */
} /* namespace moris */

#endif /* SRC_FEM_CL_VIS_VISUALIZATION_SET_HPP_ */
