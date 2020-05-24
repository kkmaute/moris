//
// Created by christopherson on 5/19/20.
//

#ifndef MORIS_FN_GEN_CREATE_PROPERTIES_HPP
#define MORIS_FN_GEN_CREATE_PROPERTIES_HPP

#include "cl_GEN_Property.hpp"
#include "cl_Param_List.hpp"
#include "fn_Exec_load_user_library.hpp"

namespace moris
{
    namespace ge
    {
        /**
         * Higher-level call for creating properties, which ensures that all property dependencies are resolved correctly
         *
         * @param aPropertyParameterLists Parameter lists for creating Property classes
         * @param aADVs Reference to the initial adv vector
         * @param aLibrary pointer to library for loading user-defined functions
         * @return Pointer to specific Property class
         */
        Cell<std::shared_ptr<Property>> create_properties( Cell<ParameterList> aPropertyParameterList,
                                                           Matrix<DDRMat>& aADVs,
                                                           std::shared_ptr<moris::Library_IO> aLibrary = nullptr );
    }
}

#endif //MORIS_FN_GEN_CREATE_PROPERTIES_HPP
