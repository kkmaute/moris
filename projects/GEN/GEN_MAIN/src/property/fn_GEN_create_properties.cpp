#include "fn_GEN_create_properties.hpp"
#include "fn_GEN_create_properties.hpp"
#include "fn_Parsing_Tools.hpp"

#include "cl_GEN_Scaled_Field.hpp"
#include "cl_GEN_Discrete_Property.hpp"
#include "cl_GEN_User_Defined_Property.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Cell<std::shared_ptr<Property>> create_properties(
                Cell<ParameterList>             aPropertyParameterLists,
                Matrix<DDRMat>&                 aADVs,
                Cell<std::shared_ptr<Geometry>> aGeometries,
                std::shared_ptr<Library_IO>     aLibrary)
        {
            // Initialize
            uint tNumProperties = aPropertyParameterLists.size();
            Cell<std::shared_ptr<Property>> tProperties(tNumProperties);
            Cell<std::string> tPropertyNames(tNumProperties);
            Cell<Cell<std::string>> tNeededFieldNames(tNumProperties);
            Cell<Cell<std::shared_ptr<Field>>> tNeededFields(tNumProperties);
            
            // Fill names, dependencies
            for (uint tPropertyIndex = 0; tPropertyIndex < tNumProperties; tPropertyIndex++)
            {
                tPropertyNames(tPropertyIndex) = aPropertyParameterLists(tPropertyIndex).get<std::string>("name");
                tNeededFieldNames(tPropertyIndex) = 
                        string_to_cell<std::string>(aPropertyParameterLists(tPropertyIndex).get<std::string>("dependencies"));
                tNeededFields(tPropertyIndex).resize(tNeededFieldNames(tPropertyIndex).size());
            }
            
            // Build based on dependencies (this is not optimally efficient, but doesn't need to be)
            bool tBuild;
            uint tNumPropertiesLeft = tNumProperties;
            uint tLoopCount = 0;
            while (tNumPropertiesLeft > 0)
            {
                for (uint tBuildPropertyIndex = 0; tBuildPropertyIndex < tNumProperties; tBuildPropertyIndex++)
                {
                    // Check if property needs to be built
                    if (tProperties(tBuildPropertyIndex) == nullptr)
                    {
                        tBuild = true;
                        
                        // Check if property dependencies are built
                        if (tNeededFieldNames(tBuildPropertyIndex).size() > 0)
                        {
                            // Loop over dependencies
                            for (uint tDependencyIndex = 0; tDependencyIndex < tNeededFieldNames(tBuildPropertyIndex).size(); tDependencyIndex++)
                            {
                                // Checking each property by name
                                for (uint tCheckPropertyIndex = 0; tCheckPropertyIndex < tNumProperties; tCheckPropertyIndex++)
                                {
                                    // Name match found
                                    if (tNeededFieldNames(tBuildPropertyIndex)(tDependencyIndex) == tPropertyNames(tCheckPropertyIndex))
                                    {
                                        // Dependency is built already
                                        if (tProperties(tCheckPropertyIndex))
                                        {
                                            tNeededFields(tBuildPropertyIndex)(tDependencyIndex) = tProperties(tCheckPropertyIndex);
                                        }

                                        // Dependency is not built, cannot build current property
                                        else
                                        {
                                            tBuild = false;
                                        }
                                    }
                                }
                            }
                        }

                        // Build
                        if (tBuild)
                        {
                            // Loop over dependencies, this time to add generic field dependencies
                            for (uint tDependencyIndex = 0; tDependencyIndex < tNeededFieldNames(tBuildPropertyIndex).size(); tDependencyIndex++)
                            {
                                // Checking each field by name
                                for (uint tCheckFieldIndex = 0; tCheckFieldIndex < aGeometries.size(); tCheckFieldIndex++)
                                {
                                    // Name match found
                                    if (tNeededFieldNames(tBuildPropertyIndex)(tDependencyIndex) == aGeometries(tCheckFieldIndex)->get_name())
                                    {
                                        tNeededFields(tBuildPropertyIndex)(tDependencyIndex) = aGeometries(tCheckFieldIndex);
                                    }
                                }
                            }

                            // Build property and decrement remaining properties to build
                            tProperties(tBuildPropertyIndex) = create_property(
                                    aPropertyParameterLists(tBuildPropertyIndex), 
                                    aADVs,
                                    tNeededFields(tBuildPropertyIndex),
                                    aLibrary);
                            tNumPropertiesLeft--;
                        }
                    }
                }
                tLoopCount++;
                MORIS_ERROR(tLoopCount <= tNumProperties, "In fn_GEN_create_properties, a circular property dependency was detected. Exiting.");
            }
            
            return tProperties;
        }

        //--------------------------------------------------------------------------------------------------------------

        Cell<std::shared_ptr<Property>> create_properties(
                Cell<ParameterList>             aPropertyParameterLists,
                sol::Dist_Vector*               aOwnedADVs,
                Cell<std::shared_ptr<Geometry>> aGeometries,
                std::shared_ptr<Library_IO>     aLibrary)
        {
            // Initialize
            uint tNumProperties = aPropertyParameterLists.size();
            Cell<std::shared_ptr<Property>> tProperties(tNumProperties);
            Cell<std::string> tPropertyNames(tNumProperties);
            Cell<Cell<std::string>> tNeededFieldNames(tNumProperties);
            Cell<Cell<std::shared_ptr<Field>>> tNeededFields(tNumProperties);

            // Fill names, dependencies
            for (uint tPropertyIndex = 0; tPropertyIndex < tNumProperties; tPropertyIndex++)
            {
                tPropertyNames(tPropertyIndex) = aPropertyParameterLists(tPropertyIndex).get<std::string>("name");
                tNeededFieldNames(tPropertyIndex) =
                        string_to_cell<std::string>(aPropertyParameterLists(tPropertyIndex).get<std::string>("dependencies"));
                tNeededFields(tPropertyIndex).resize(tNeededFieldNames(tPropertyIndex).size());
            }

            // Build based on dependencies (this is not optimally efficient, but doesn't need to be)
            bool tBuild;
            uint tNumPropertiesLeft = tNumProperties;
            uint tLoopCount = 0;
            while (tNumPropertiesLeft > 0)
            {
                for (uint tBuildPropertyIndex = 0; tBuildPropertyIndex < tNumProperties; tBuildPropertyIndex++)
                {
                    // Check if property needs to be built
                    if (tProperties(tBuildPropertyIndex) == nullptr)
                    {
                        tBuild = true;

                        // Check if property dependencies are built
                        if (tNeededFieldNames(tBuildPropertyIndex).size() > 0)
                        {
                            // Loop over dependencies
                            for (uint tDependencyIndex = 0; tDependencyIndex < tNeededFieldNames(tBuildPropertyIndex).size(); tDependencyIndex++)
                            {
                                // Checking each property by name
                                for (uint tCheckPropertyIndex = 0; tCheckPropertyIndex < tNumProperties; tCheckPropertyIndex++)
                                {
                                    // Name match found
                                    if (tNeededFieldNames(tBuildPropertyIndex)(tDependencyIndex) == tPropertyNames(tCheckPropertyIndex))
                                    {
                                        // Dependency is built already
                                        if (tProperties(tCheckPropertyIndex))
                                        {
                                            tNeededFields(tBuildPropertyIndex)(tDependencyIndex) = tProperties(tCheckPropertyIndex);
                                        }

                                        // Dependency is not built, cannot build current property
                                        else
                                        {
                                            tBuild = false;
                                        }
                                    }
                                }
                            }
                        }

                        // Build
                        if (tBuild)
                        {
                            // Loop over dependencies, this time to add generic field dependencies
                            for (uint tDependencyIndex = 0; tDependencyIndex < tNeededFieldNames(tBuildPropertyIndex).size(); tDependencyIndex++)
                            {
                                // Checking each field by name
                                for (uint tCheckFieldIndex = 0; tCheckFieldIndex < aGeometries.size(); tCheckFieldIndex++)
                                {
                                    // Name match found
                                    if (tNeededFieldNames(tBuildPropertyIndex)(tDependencyIndex) == aGeometries(tCheckFieldIndex)->get_name())
                                    {
                                        tNeededFields(tBuildPropertyIndex)(tDependencyIndex) = aGeometries(tCheckFieldIndex);
                                    }
                                }
                            }

                            // Build property and decrement remaining properties to build
                            tProperties(tBuildPropertyIndex) = create_property(
                                    aPropertyParameterLists(tBuildPropertyIndex),
                                    aOwnedADVs,
                                    tNeededFields(tBuildPropertyIndex),
                                    aLibrary);
                            tNumPropertiesLeft--;
                        }
                    }
                }
                tLoopCount++;
                MORIS_ERROR(tLoopCount <= tNumProperties, "In fn_GEN_create_properties, a circular property dependency was detected. Exiting.");
            }

            return tProperties;
        }

        //--------------------------------------------------------------------------------------------------------------

        std::shared_ptr<Property> create_property(
                ParameterList                aPropertyParameterList,
                Matrix<DDRMat>&              aADVs,
                Cell<std::shared_ptr<Field>> aFieldDependencies,
                std::shared_ptr<Library_IO>  aLibrary)
        {
            // Property type/name
            std::string tPropertyType = aPropertyParameterList.get<std::string>("type");
            std::string tPropertyName = aPropertyParameterList.get<std::string>("name");

            // Property inputs
            Matrix<DDUMat> tPropertyVariableIndices(0, 0);
            Matrix<DDUMat> tADVIndices(0, 0);
            set_property_variable_inputs(aPropertyParameterList, aADVs.length(), tPropertyVariableIndices, tADVIndices);
            Matrix<DDRMat> tConstantParameters(0, 0);

            // Get constant parameters
            string_to_mat(aPropertyParameterList.get<std::string>("constant_parameters"), tConstantParameters);

            // Build Property
            if (tPropertyType == "scaled_field")
            {
                return std::make_shared<Scaled_Field>(
                        aADVs,
                        tPropertyVariableIndices,
                        tADVIndices,
                        tConstantParameters,
                        aFieldDependencies,
                        tPropertyName);
            }
            else if (tPropertyType == "discrete")
            {
                return std::make_shared<Discrete_Property>(
                        aADVs,
                        tPropertyVariableIndices,
                        tADVIndices,
                        tConstantParameters,
                        tPropertyName);
            }
            else if (tPropertyType == "user_defined")
            {
                // Check if library is given
                MORIS_ERROR(aLibrary != nullptr, "Library must be given in order to create a user-defined geometry.");

                // Get sensitivity function if needed
                std::string tSensitivityFunctionName = aPropertyParameterList.get<std::string>("sensitivity_function_name");
                MORIS_GEN_SENSITIVITY_FUNCTION tSensitivityFunction =
                        (tSensitivityFunctionName == "" ? nullptr : aLibrary->load_gen_sensitivity_function(tSensitivityFunctionName));

                return std::make_shared<User_Defined_Property>(
                        aADVs,
                        tPropertyVariableIndices,
                        tADVIndices,
                        tConstantParameters,
                        aFieldDependencies,
                        aLibrary->load_gen_field_function(aPropertyParameterList.get<std::string>("field_function_name")),
                        tSensitivityFunction,
                        tPropertyName);
            }
            else
            {
                MORIS_ERROR(false, tPropertyType.append(" is not recognized as a valid Property type in fn_GEN_create_property.").c_str());
                return nullptr;
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        std::shared_ptr<Property> create_property(
                ParameterList                aPropertyParameterList,
                sol::Dist_Vector*            aOwnedADVs,
                Cell<std::shared_ptr<Field>> aFieldDependencies,
                std::shared_ptr<Library_IO>  aLibrary)
        {
            // Property type/name
            std::string tPropertyType = aPropertyParameterList.get<std::string>("type");
            std::string tPropertyName = aPropertyParameterList.get<std::string>("name");

            // Property inputs
            Matrix<DDUMat> tPropertyVariableIndices(0, 0);
            Matrix<DDUMat> tADVIndices(0, 0);
            Matrix<DDRMat> tConstantParameters(0, 0);

            // Get constant/variable parameters
            set_property_variable_inputs(aPropertyParameterList, aOwnedADVs->vec_local_length(), tPropertyVariableIndices, tADVIndices);
            string_to_mat(aPropertyParameterList.get<std::string>("constant_parameters"), tConstantParameters);

            // Build Property
            if (tPropertyType == "scaled_field")
            {
                return std::make_shared<Scaled_Field>(
                        aOwnedADVs,
                        tPropertyVariableIndices,
                        tADVIndices,
                        tConstantParameters,
                        aFieldDependencies,
                        tPropertyName);
            }
            else if (tPropertyType == "discrete")
            {
                return std::make_shared<Discrete_Property>(
                        aOwnedADVs,
                        tPropertyVariableIndices,
                        tADVIndices,
                        tConstantParameters,
                        tPropertyName);
            }
            else if (tPropertyType == "user_defined")
            {
                // Check if library is given
                MORIS_ERROR(aLibrary != nullptr, "Library must be given in order to create a user-defined geometry.");

                // Get sensitivity function if needed
                std::string tSensitivityFunctionName = aPropertyParameterList.get<std::string>("sensitivity_function_name");
                MORIS_GEN_SENSITIVITY_FUNCTION tSensitivityFunction =
                        (tSensitivityFunctionName == "" ? nullptr : aLibrary->load_gen_sensitivity_function(tSensitivityFunctionName));

                return std::make_shared<User_Defined_Property>(
                        aOwnedADVs,
                        tPropertyVariableIndices,
                        tADVIndices,
                        tConstantParameters,
                        aFieldDependencies,
                        aLibrary->load_gen_field_function(aPropertyParameterList.get<std::string>("field_function_name")),
                        tSensitivityFunction,
                        tPropertyName);
            }
            else
            {
                MORIS_ERROR(false, tPropertyType.append(" is not recognized as a valid Property type in fn_GEN_create_property.").c_str());
                return nullptr;
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void set_property_variable_inputs(
                ParameterList  aPropertyParameterList,
                uint           aNumADVs,
                Matrix<DDUMat>& aPropertyVariableIndices,
                Matrix<DDUMat>& aADVIndices)
        {
            bool tFillVariables = false;
            bool tFillADVs = false;

            // Determine if variable or ADV indices need to be filled (specified by "all")
            if (aPropertyParameterList.get<std::string>("property_variable_indices") == "all")
            {
                tFillVariables = true;
            }
            else
            {
                string_to_mat(aPropertyParameterList.get<std::string>("property_variable_indices"), aPropertyVariableIndices);
            }
            if (aPropertyParameterList.get<std::string>("adv_indices") == "all")
            {
                tFillADVs = true;
            }
            else
            {
                string_to_mat(aPropertyParameterList.get<std::string>("adv_indices"), aADVIndices);
            }

            // Perform fill
            if (tFillVariables and tFillADVs)
            {
                aPropertyVariableIndices.resize(aNumADVs, 1);
                aADVIndices.resize(aNumADVs, 1);
                for (uint tIndex = 0; tIndex < aNumADVs; tIndex++)
                {
                    aPropertyVariableIndices(tIndex) = tIndex;
                    aADVIndices(tIndex) = tIndex;
                }
            }
            else if (tFillVariables)
            {
                aPropertyVariableIndices.resize(aADVIndices.length(), 1);
                for (uint tIndex = 0; tIndex < aADVIndices.length(); tIndex++)
                {
                    aPropertyVariableIndices(tIndex) = tIndex;
                }
            }
            else if (tFillADVs)
            {
                aADVIndices.resize(aPropertyVariableIndices.length(), 1);
                for (uint tIndex = 0; tIndex < aPropertyVariableIndices.length(); tIndex++)
                {
                    aADVIndices(tIndex) = tIndex;
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
