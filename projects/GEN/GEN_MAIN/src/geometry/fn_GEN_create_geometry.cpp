//
// Created by christopherson on 4/13/20.
//

#include "fn_GEN_create_geometry.hpp"
#include "fn_Parsing_Tools.hpp"

#include "cl_GEN_Circle.hpp"

namespace moris
{
    namespace ge
    {
        std::shared_ptr<Geometry> create_geometry(ParameterList aGeometryParameterList, Matrix<DDRMat>& aADVs)
        {
            // Geometry type
            std::string tGeometryType = aGeometryParameterList.get<std::string>("type");

            // Geometry inputs
            Matrix<DDUMat> tGeometryVariableIndices(1, 1);
            Matrix<DDUMat> tADVIndices(1, 1);
            Matrix<DDRMat> tConstantParameters(1, 1);

            // Get from parameter list
            string_to_mat(aGeometryParameterList.get<std::string>("geometry_variable_indices"), tGeometryVariableIndices);
            string_to_mat(aGeometryParameterList.get<std::string>("adv_indices"), tADVIndices);
            string_to_mat(aGeometryParameterList.get<std::string>("constant_parameters"), tConstantParameters);

            // Build Geometry
            if (tGeometryType == "circle")
            {
                return std::make_shared<Circle>(aADVs, tGeometryVariableIndices, tADVIndices, tConstantParameters);
            }
            else
            {
                MORIS_ERROR(false, tGeometryType.append(" is not recognized as a valid Geometry type in fn_GEN_create_geometry.").c_str());
                return nullptr;
            }
        }
    }
}
