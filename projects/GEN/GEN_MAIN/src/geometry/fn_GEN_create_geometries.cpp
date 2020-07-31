#include "fn_GEN_create_geometries.hpp"
#include "fn_Parsing_Tools.hpp"

#include "cl_GEN_Circle.hpp"
#include "cl_GEN_Sphere.hpp"
#include "cl_GEN_Plane.hpp"
#include "cl_GEN_User_Defined_Geometry.hpp"
#include "cl_GEN_Level_Set.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Cell<std::shared_ptr<Geometry>> create_geometries(
                Cell<ParameterList>         aGeometryParameterLists,
                Matrix<DDRMat>&             aADVs,
                std::shared_ptr<Library_IO> aLibrary)
        {
            // Create geometry cell
            Cell<std::shared_ptr<Geometry>> tGeometries(aGeometryParameterLists.size());

            // Create individual geometries
            for (uint tGeometryIndex = 0; tGeometryIndex < aGeometryParameterLists.size(); tGeometryIndex++)
            {
                // Use parameter list to create geometry
                tGeometries(tGeometryIndex) = create_geometry(aGeometryParameterLists(tGeometryIndex), aADVs, aLibrary);
            }

            return tGeometries;
        }

        //--------------------------------------------------------------------------------------------------------------

        std::shared_ptr<Geometry> create_geometry(
                ParameterList               aGeometryParameterList,
                Matrix<DDRMat>&             aADVs,
                std::shared_ptr<Library_IO> aLibrary)
        {
            // Geometry type
            std::string tGeometryType = aGeometryParameterList.get<std::string>("type");

            // Geometry inputs
            Matrix<DDUMat> tGeometryVariableIndices(0, 0);
            Matrix<DDUMat> tADVIndices(0, 0);
            set_geometry_variable_inputs(aGeometryParameterList, aADVs.length(), tGeometryVariableIndices, tADVIndices);

            // Get constant parameters
            Matrix<DDRMat> tConstantParameters =
                    string_to_mat<DDRMat>(aGeometryParameterList.get<std::string>("constant_parameters"));

            // Get refinement info
            sint tNumRefinements = aGeometryParameterList.get<sint>("number_of_refinements");
            sint tRefinementFunctionIndex = aGeometryParameterList.get<sint>("refinement_function_index");

            // Get level set info
            sint tBSplineMeshIndex = aGeometryParameterList.get<sint>("bspline_mesh_index");
            real tLevelSetLowerBound = aGeometryParameterList.get<real>("level_set_lower_bound");
            real tLevelSetUpperBound = aGeometryParameterList.get<real>("level_set_upper_bound");

            // Build Geometry
            if (tGeometryType == "circle")
            {
                return std::make_shared<Circle>(aADVs,
                                                tGeometryVariableIndices,
                                                tADVIndices,
                                                tConstantParameters,
                                                tNumRefinements,
                                                tRefinementFunctionIndex,
                                                tBSplineMeshIndex,
                                                tLevelSetLowerBound,
                                                tLevelSetUpperBound);
            }
            else if (tGeometryType == "sphere")
            {
                return std::make_shared<Sphere>(aADVs,
                                                tGeometryVariableIndices,
                                                tADVIndices,
                                                tConstantParameters,
                                                tNumRefinements,
                                                tRefinementFunctionIndex,
                                                tBSplineMeshIndex,
                                                tLevelSetLowerBound,
                                                tLevelSetUpperBound);
            }
            else if (tGeometryType == "plane")
            {
                return std::make_shared<Plane>(aADVs,
                                               tGeometryVariableIndices,
                                               tADVIndices,
                                               tConstantParameters,
                                               tNumRefinements,
                                               tRefinementFunctionIndex,
                                               tBSplineMeshIndex,
                                               tLevelSetLowerBound,
                                               tLevelSetUpperBound);
            }
            else if (tGeometryType == "user_defined")
            {
                // Check if library is given
                MORIS_ERROR(aLibrary != nullptr, "Library must be given in order to create a user-defined geometry.");

                // Get sensitivity function if needed
                std::string tSensitivityFunctionName = aGeometryParameterList.get<std::string>("sensitivity_function_name");
                MORIS_GEN_SENSITIVITY_FUNCTION tSensitivityFunction =
                        (tSensitivityFunctionName == "" ? nullptr : aLibrary->load_gen_sensitivity_function(tSensitivityFunctionName));

                // Create user-defined geometry
                return std::make_shared<User_Defined_Geometry>(aADVs,
                                                               tGeometryVariableIndices,
                                                               tADVIndices,
                                                               tConstantParameters,
                                                               aLibrary->load_gen_field_function(aGeometryParameterList.get<std::string>("field_function_name")),
                                                               tSensitivityFunction,
                                                               tNumRefinements,
                                                               tRefinementFunctionIndex,
                                                               tBSplineMeshIndex,
                                                               tLevelSetLowerBound,
                                                               tLevelSetUpperBound);
            }
            else
            {
                MORIS_ERROR(false, tGeometryType.append(" is not recognized as a valid Geometry type in fn_GEN_create_geometry.").c_str());
                return nullptr;
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void set_geometry_variable_inputs(
                ParameterList  aGeometryParameterList,
                uint           aNumADVs,
                Matrix<DDUMat>& aGeometryVariableIndices,
                Matrix<DDUMat>& aADVIndices)
        {
            bool tFillVariables = false;
            bool tFillADVs = false;

            // Determine if variable or ADV indices need to be filled (specified by "all")
            if (aGeometryParameterList.get<std::string>("geometry_variable_indices") == "all")
            {
                tFillVariables = true;
            }
            else
            {
                string_to_mat(aGeometryParameterList.get<std::string>("geometry_variable_indices"), aGeometryVariableIndices);
            }
            if (aGeometryParameterList.get<std::string>("adv_indices") == "all")
            {
                tFillADVs = true;
            }
            else
            {
                string_to_mat(aGeometryParameterList.get<std::string>("adv_indices"), aADVIndices);
            }

            // Perform fill
            if (tFillVariables and tFillADVs)
            {
                aGeometryVariableIndices.resize(aNumADVs, 1);
                aADVIndices.resize(aNumADVs, 1);
                for (uint tIndex = 0; tIndex < aNumADVs; tIndex++)
                {
                    aGeometryVariableIndices(tIndex) = tIndex;
                    aADVIndices(tIndex) = tIndex;
                }
            }
            else if (tFillVariables)
            {
                aGeometryVariableIndices.resize(aADVIndices.length(), 1);
                for (uint tIndex = 0; tIndex < aADVIndices.length(); tIndex++)
                {
                    aGeometryVariableIndices(tIndex) = tIndex;
                }
            }
            else if (tFillADVs)
            {
                aADVIndices.resize(aGeometryVariableIndices.length(), 1);
                for (uint tIndex = 0; tIndex < aGeometryVariableIndices.length(); tIndex++)
                {
                    aADVIndices(tIndex) = tIndex;
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
