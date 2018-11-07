/*
 * MSI_Adof_Order_Hack.hpp
 *
 *  Created on: Oct 26, 2018
 *      Author: messe
 */

#ifndef PROJECTS_FEM_MSI_SRC_MSI_ADOF_ORDER_HACK_HPP_
#define PROJECTS_FEM_MSI_SRC_MSI_ADOF_ORDER_HACK_HPP_

#include "typedefs.hpp"

namespace moris
{
    namespace MSI
    {
        /**
         * fixme: We have not fully decided yet how to call the
         *        B-Spline orders from within the model.
         *        Since Christian does not want to mess around with
         *        Mathias' API, this global variable is created
         *        as a temporary hack.
         */
        uint gAdofOrderHack = 0;
    }
}

#endif /* PROJECTS_FEM_MSI_SRC_MSI_ADOF_ORDER_HACK_HPP_ */
