//
// Created by christopherson on 5/19/20.
//

#ifndef MORIS_FN_GEN_CREATE_PROPERTY_HPP
#define MORIS_FN_GEN_CREATE_PROPERTY_HPP

#include "cl_GEN_Property.hpp"
#include "cl_Param_List.hpp"
#include "fn_Exec_load_user_library.hpp"

namespace moris
{
    namespace ge
    {
        /**
         * Creates an instance of the specified Property class and returns a shared pointer to it
         *
         * @param aPropertyParameterList Parameter list for creating a Property class
         * @param aADVs Reference to the initial adv vector
         * @param aLibrary pointer to library for loading user-defined functions
         * @return Pointer to specific Property class
         */
        std::shared_ptr<Property> create_property( ParameterList aPropertyParameterList,
                                                   Matrix<DDRMat>& aADVs,
                                                   Cell<std::shared_ptr<Property>> aPropertyDependencies,
                                                   std::shared_ptr<moris::Library_IO> aLibrary = nullptr );
    }
}

#endif //MORIS_FN_GEN_CREATE_PROPERTY_HPP
