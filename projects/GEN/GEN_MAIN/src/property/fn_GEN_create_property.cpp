//
// Created by christopherson on 5/19/20.
//

#include "fn_GEN_create_property.hpp"
#include "fn_Parsing_Tools.hpp"

#include "cl_GEN_Discrete_Property.hpp"
#include "cl_GEN_User_Defined_Property.hpp"

namespace moris
{
    namespace ge
    {
        std::shared_ptr<Property> create_property(ParameterList aPropertyParameterList,
                                                  Matrix<DDRMat>& aADVs,
                                                  Cell<std::shared_ptr<Property>> aPropertyDependencies,
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
            print(tPropertyVariableIndices, "prop ind");
            print(tADVIndices, "adv ind");
            if (tPropertyType == "discrete")
            {
                return std::make_shared<Discrete_Property>(aADVs, tPropertyVariableIndices, tADVIndices, tConstantParameters);
            }
            else if (tPropertyType == "user_defined")
            {
                return std::make_shared<User_Defined_Property>(
                        aADVs, tPropertyVariableIndices, tADVIndices, tConstantParameters, aPropertyDependencies,
                        aLibrary->load_gen_field_function(aPropertyParameterList.get<std::string>("field_function_name")),
                        aLibrary->load_gen_sensitivity_function(aPropertyParameterList.get<std::string>("sensitivity_function_name")));
            }
            else
            {
                MORIS_ERROR(false, tPropertyType.append(" is not recognized as a valid Property type in fn_GEN_create_property.").c_str());
                return nullptr;
            }
        }
    }
}
