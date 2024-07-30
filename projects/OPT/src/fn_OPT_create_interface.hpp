/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_OPT_create_interface.hpp
 *
 */

#ifndef MORIS_FN_OPT_CREATE_INTERFACE_HPP
#define MORIS_FN_OPT_CREATE_INTERFACE_HPP

#include "cl_OPT_Criteria_Interface.hpp"
#include "cl_Parameter_List.hpp"

namespace moris
{
    namespace opt
    {
        /**
         * Creates an instance of an Interface class or Interface_Manager and returns a shared pointer to it.
         *
         * @param aParameterLists Parameter lists for individual interfaces
         * @return Interface class (can be manager)
         */
        std::shared_ptr<Criteria_Interface> create_interface(
                Vector< Parameter_List > aParameterLists,
                Vector<std::shared_ptr<Criteria_Interface>> aInterfaces = Vector<std::shared_ptr<Criteria_Interface>>(0));

        /**
         * Creates an instance of the specified Interface class and returns a shared pointer to it
         *
         * @param aParameterList A single interface parameter list
         * @return Interface class
         */
        std::shared_ptr<Criteria_Interface> create_interface( Parameter_List aParameterList);

    }
}

#endif //MORIS_FN_OPT_CREATE_INTERFACE_HPP

