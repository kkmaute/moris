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
         * Higher-level call for creating a cell of properties, which ensures that all property dependencies are
         * resolved correctly.
         *
         * @param aPropertyParameterLists Parameter lists for creating Property classes
         * @param aADVs Reference to the initial adv vector
         * @param aLibrary pointer to library for loading user-defined functions
         * @return Pointer to specific Property class
         */
        Cell<std::shared_ptr<Property>> create_properties(
                Cell<ParameterList>         aPropertyParameterLists,
                Matrix<DDRMat>&             aADVs,
                std::shared_ptr<Library_IO> aLibrary = nullptr);

        /**
         * Higher-level call for creating a cell of properties, which ensures that all property dependencies are
         * resolved correctly.
         *
         * @param aPropertyParameterLists Parameter lists for creating Property classes
         * @param aOwnedADVs Distributed owned ADVs
         * @param aLibrary pointer to library for loading user-defined functions
         * @return Pointer to specific Property class
         */
        Cell<std::shared_ptr<Property>> create_properties(
                Cell<ParameterList>         aPropertyParameterLists,
                sol::Dist_Vector*           aOwnedADVs,
                std::shared_ptr<Library_IO> aLibrary = nullptr);

        /**
         * Creates an instance of the specified Property class and returns a shared pointer to it.
         *
         * @param aPropertyParameterList Parameter list for creating a Property class
         * @param aADVs Reference to the initial adv vector
         * @param aLibrary pointer to library for loading user-defined functions
         * @return Pointer to specific Property class
         */
        std::shared_ptr<Property> create_property(
                ParameterList                   aPropertyParameterList,
                Matrix<DDRMat>&                 aADVs,
                Cell<std::shared_ptr<Property>> aPropertyDependencies,
                std::shared_ptr<Library_IO>     aLibrary = nullptr);

        /**
         * Creates an instance of the specified Property class and returns a shared pointer to it.
         *
         * @param aPropertyParameterList Parameter list for creating a Property class
         * @param aOwnedADVs Distributed owned ADVs
         * @param aLibrary pointer to library for loading user-defined functions
         * @return Pointer to specific Property class
         */
        std::shared_ptr<Property> create_property(
                ParameterList                   aPropertyParameterList,
                sol::Dist_Vector*               aOwnedADVs,
                Cell<std::shared_ptr<Property>> aPropertyDependencies,
                std::shared_ptr<Library_IO>     aLibrary = nullptr);

        /**
         * Sets the property variables which depend on ADVs. Used by create_property().
         *
         * @param aPropertyParameterList Parameter list for creating a property class
         * @param aNumADVs The number of total ADVs
         * @param aPropertyVariableIndices Indices of property variables to fill with ADVs
         * @param aADVIndices Indices of ADVs for filling the property variables
         */
        void set_property_variable_inputs(
                ParameterList  aPropertyParameterList,
                uint           aNumADVs,
                Matrix<DDUMat>& aPropertyVariableIndices,
                Matrix<DDUMat>& aADVIndices);
    }
}

#endif //MORIS_FN_GEN_CREATE_PROPERTIES_HPP
