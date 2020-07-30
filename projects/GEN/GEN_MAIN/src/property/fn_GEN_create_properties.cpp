#include "fn_GEN_create_properties.hpp"
#include "fn_GEN_create_properties.hpp"
#include "fn_Parsing_Tools.hpp"

#include "cl_GEN_Discrete_Property.hpp"
#include "cl_GEN_User_Defined_Property.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Cell<std::shared_ptr<Property>> create_properties(
                Cell<ParameterList>                aPropertyParameterLists,
                Matrix<DDRMat>&                    aADVs,
                std::shared_ptr<moris::Library_IO> aLibrary)
        {
            // Initialize
            uint tNumProperties = aPropertyParameterLists.size();
            Cell<std::shared_ptr<Property>> tProperties(tNumProperties);
            Cell<std::string> tPropertyNames(tNumProperties);
            Cell<Cell<std::string>> tPropertyDependencies(tNumProperties);
            Cell<Cell<std::shared_ptr<Property>>> tNeededProperties(tNumProperties);
            
            // Fill names, dependencies
            for (uint tPropertyIndex = 0; tPropertyIndex < tNumProperties; tPropertyIndex++)
            {
                tPropertyNames(tPropertyIndex) = aPropertyParameterLists(tPropertyIndex).get<std::string>("name");
                tPropertyDependencies(tPropertyIndex) = 
                        string_to_cell<std::string>(aPropertyParameterLists(tPropertyIndex).get<std::string>("dependencies"));
                tNeededProperties(tPropertyIndex).resize(tPropertyDependencies(tPropertyIndex).size());
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
                        
                        // Check if dependencies are built
                        if (tPropertyDependencies(tBuildPropertyIndex).size() > 0)
                        {
                            for (uint tDependencyIndex = 0; tDependencyIndex < tPropertyDependencies(tBuildPropertyIndex).size(); tDependencyIndex++)
                            {
                                for (uint tCheckPropertyIndex = 0; tCheckPropertyIndex < tNumProperties; tCheckPropertyIndex++)
                                {
                                    if (tPropertyDependencies(tBuildPropertyIndex)(tDependencyIndex) == tPropertyNames(tCheckPropertyIndex)
                                    && tProperties(tCheckPropertyIndex) == nullptr)
                                    {
                                        tBuild = false;
                                    }
                                }
                            }
                        }

                        // Build
                        if (tBuild)
                        {
                            tProperties(tBuildPropertyIndex) = create_property(aPropertyParameterLists(tBuildPropertyIndex), aADVs, tNeededProperties(tBuildPropertyIndex), aLibrary);
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
                ParameterList                      aPropertyParameterList,
                Matrix<DDRMat>&                    aADVs,
                Cell<std::shared_ptr<Property>>    aPropertyDependencies,
                std::shared_ptr<moris::Library_IO> aLibrary)
        {
            // Property type
            std::string tPropertyType = aPropertyParameterList.get<std::string>("type");

            // Property inputs
            Matrix<DDUMat> tPropertyVariableIndices(0, 0);
            Matrix<DDUMat> tADVIndices(0, 0);
            Matrix<DDRMat> tConstantParameters(0, 0);
            bool tFillVariables = false;
            bool tFillADVs = false;

            // Determine if variable or ADV indices need to be filled (specified by "all")
            if (aPropertyParameterList.get<std::string>("property_variable_indices") == "all")
            {
                tFillVariables = true;
            }
            else
            {
                string_to_mat(aPropertyParameterList.get<std::string>("property_variable_indices"), tPropertyVariableIndices);
            }
            if (aPropertyParameterList.get<std::string>("adv_indices") == "all")
            {
                tFillADVs = true;
            }
            else
            {
                string_to_mat(aPropertyParameterList.get<std::string>("adv_indices"), tADVIndices);
            }

            // Perform fill
            if (tFillVariables and tFillADVs)
            {
                tPropertyVariableIndices.resize(aADVs.length(), 1);
                tADVIndices.resize(aADVs.length(), 1);
                for (uint tIndex = 0; tIndex < aADVs.length(); tIndex++)
                {
                    tPropertyVariableIndices(tIndex) = tIndex;
                    tADVIndices(tIndex) = tIndex;
                }
            }
            else if (tFillVariables)
            {
                tPropertyVariableIndices.resize(tADVIndices.length(), 1);
                for (uint tIndex = 0; tIndex < tADVIndices.length(); tIndex++)
                {
                    tPropertyVariableIndices(tIndex) = tIndex;
                }
            }
            else if (tFillADVs)
            {
                tADVIndices.resize(tPropertyVariableIndices.length(), 1);
                for (uint tIndex = 0; tIndex < tPropertyVariableIndices.length(); tIndex++)
                {
                    tADVIndices(tIndex) = tIndex;
                }
            }

            // Get constant parameters
            string_to_mat(aPropertyParameterList.get<std::string>("constant_parameters"), tConstantParameters);

            // Build Property
            if (tPropertyType == "discrete")
            {
                return std::make_shared<Discrete_Property>(aADVs, tPropertyVariableIndices, tADVIndices, tConstantParameters);
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
                        aADVs, tPropertyVariableIndices, tADVIndices, tConstantParameters, aPropertyDependencies,
                        aLibrary->load_gen_field_function(aPropertyParameterList.get<std::string>("field_function_name")),
                        tSensitivityFunction);
            }
            else
            {
                MORIS_ERROR(false, tPropertyType.append(" is not recognized as a valid Property type in fn_GEN_create_property.").c_str());
                return nullptr;
            }
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
