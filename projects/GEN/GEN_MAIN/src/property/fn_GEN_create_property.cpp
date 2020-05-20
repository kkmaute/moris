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
            Matrix<DDUMat> tPropertyVariableIndices(1, 1);
            Matrix<DDUMat> tADVIndices(1, 1);
            Matrix<DDRMat> tConstantParameters(1, 1);

            // Get from parameter list
            string_to_mat(aPropertyParameterList.get<std::string>("property_variable_indices"), tPropertyVariableIndices);
            string_to_mat(aPropertyParameterList.get<std::string>("adv_indices"), tADVIndices);
            string_to_mat(aPropertyParameterList.get<std::string>("constant_parameters"), tConstantParameters);

            // Build Property
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
