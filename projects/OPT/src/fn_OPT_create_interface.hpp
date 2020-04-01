//
// Created by christopherson on 3/4/20.
//

#ifndef MORIS_FN_OPT_CREATE_INTERFACE_HPP
#define MORIS_FN_OPT_CREATE_INTERFACE_HPP

#include "cl_OPT_Interface.hpp"
#include "cl_Param_List.hpp"

namespace moris
{
    namespace opt
    {
        /**
         * Creates an instance of the specified Problem class and returns a shared pointer to it
         *
         * @param aParameterList moris::ParameterList
         * @return std::shared_ptr<Problem>
         */
        std::shared_ptr<Interface> create_interface(ParameterList aParameterList);
    }
}

#endif //MORIS_FN_OPT_CREATE_INTERFACE_HPP
