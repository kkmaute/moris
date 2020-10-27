#include "fn_GEN_create_geometries.hpp"
#include "fn_Parsing_Tools.hpp"

#include "cl_GEN_Circle.hpp"
#include "cl_GEN_Superellipse.hpp"
#include "cl_GEN_Sphere.hpp"
#include "cl_GEN_Superellipsoid.hpp"
#include "cl_GEN_Plane.hpp"
#include "cl_GEN_User_Defined_Geometry.hpp"
#include "cl_GEN_Level_Set.hpp"
#include "cl_GEN_Multigeometry.hpp"
#include "cl_GEN_Swiss_Cheese_Slice.hpp"

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
            Cell<std::shared_ptr<Geometry>> tGeometries(0);
            Cell<std::shared_ptr<Multigeometry>> tMultigeometries(0);

            // Create individual geometries
            for (uint tGeometryIndex = 0; tGeometryIndex < aGeometryParameterLists.size(); tGeometryIndex++)
            {
                // Create geometry
                std::shared_ptr<Geometry> tGeometry = create_geometry(aGeometryParameterLists(tGeometryIndex), aADVs, aLibrary);

                // Determine if to add to multigeometry
                bool tMultigeometryFound = false;
                std::string tGeometryName = tGeometry->get_name();
                if (tGeometryName != "")
                {
                    // Loop to see if this multigeometry ID exists already
                    for (uint tMultigeometryIndex = 0; tMultigeometryIndex < tMultigeometries.size(); tMultigeometryIndex++)
                    {
                        if (tMultigeometries(tMultigeometryIndex)->get_name() == tGeometryName)
                        {
                            tMultigeometryFound = true;
                            tMultigeometries(tMultigeometryIndex)->add_geometry(tGeometry);
                            break;
                        }
                    }

                    // Check for creating new multigeometry
                    if (not tMultigeometryFound)
                    {
                        for (uint tCreatedGeometryIndex = 0; tCreatedGeometryIndex < tGeometries.size(); tCreatedGeometryIndex++)
                        {
                            if (tGeometries(tCreatedGeometryIndex)->get_name() == tGeometryName)
                            {
                                tMultigeometryFound = true;
                                tMultigeometries.push_back(std::make_shared<Multigeometry>(
                                        Cell<std::shared_ptr<Geometry>>({tGeometries(tCreatedGeometryIndex), tGeometry}),
                                        tGeometryName));
                                tGeometries.erase(tCreatedGeometryIndex);
                                break;
                            }
                        }
                    }
                }

                // If no multigeometry, add as regular geometry
                if (not tMultigeometryFound)
                {
                    tGeometries.push_back(tGeometry);
                }
            }

            // Add multigeometries at the end
            for (uint tMultigeometryIndex = 0; tMultigeometryIndex < tMultigeometries.size(); tMultigeometryIndex++)
            {
                tGeometries.push_back(tMultigeometries(tMultigeometryIndex));
            }

            return tGeometries;
        }

        //--------------------------------------------------------------------------------------------------------------

        Cell<std::shared_ptr<Geometry>> create_geometries(
                Cell<ParameterList>         aGeometryParameterLists,
                sol::Dist_Vector*           aOwnedADVs,
                std::shared_ptr<Library_IO> aLibrary)
        {
            // Create geometry cell
            Cell<std::shared_ptr<Geometry>> tGeometries(0);
            Cell<std::shared_ptr<Multigeometry>> tMultigeometries(0);

            // Create individual geometries
            for (uint tGeometryIndex = 0; tGeometryIndex < aGeometryParameterLists.size(); tGeometryIndex++)
            {
                // Create geometry
                std::shared_ptr<Geometry> tGeometry = create_geometry(aGeometryParameterLists(tGeometryIndex), aOwnedADVs, aLibrary);

                // Determine if to add to multigeometry
                bool tMultigeometryFound = false;
                std::string tGeometryName = tGeometry->get_name();
                if (tGeometryName != "")
                {
                    // Loop to see if this multigeometry ID exists already
                    for (uint tMultigeometryIndex = 0; tMultigeometryIndex < tMultigeometries.size(); tMultigeometryIndex++)
                    {
                        if (tMultigeometries(tMultigeometryIndex)->get_name() == tGeometryName)
                        {
                            tMultigeometryFound = true;
                            tMultigeometries(tMultigeometryIndex)->add_geometry(tGeometry);
                            break;
                        }
                    }

                    // Check for creating new multigeometry
                    if (not tMultigeometryFound)
                    {
                        for (uint tCreatedGeometryIndex = 0; tCreatedGeometryIndex < tGeometries.size(); tCreatedGeometryIndex++)
                        {
                            if (tGeometries(tCreatedGeometryIndex)->get_name() == tGeometryName)
                            {
                                tMultigeometryFound = true;
                                tMultigeometries.push_back(std::make_shared<Multigeometry>(
                                        Cell<std::shared_ptr<Geometry>>({tGeometries(tCreatedGeometryIndex), tGeometry}),
                                        tGeometryName));
                                tGeometries.erase(tCreatedGeometryIndex);
                                break;
                            }
                        }
                    }
                }

                // If no multigeometry, add as regular geometry
                if (not tMultigeometryFound)
                {
                    tGeometries.push_back(tGeometry);
                }
            }

            // Add multigeometries at the end
            for (uint tMultigeometryIndex = 0; tMultigeometryIndex < tMultigeometries.size(); tMultigeometryIndex++)
            {
                tGeometries.push_back(tMultigeometries(tMultigeometryIndex));
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
            Matrix<DDRMat> tConstantParameters(0, 0);

            // If not a swiss cheese, get ADV inputs
            if (tGeometryType.compare(0, 12, "swiss_cheese"))
            {
                // ADV parameters
                set_geometry_variable_inputs(aGeometryParameterList, aADVs.length(), tGeometryVariableIndices, tADVIndices);

                // Constant parameters
                tConstantParameters = string_to_mat<DDRMat>(aGeometryParameterList.get<std::string>("constant_parameters"));
            }

            // Get name
            std::string tGeometryName = aGeometryParameterList.get<std::string>("name");

            // Get refinement info
            sint tNumRefinements = aGeometryParameterList.get<sint>("number_of_refinements");
            sint tRefinementFunctionIndex = aGeometryParameterList.get<sint>("refinement_function_index");

            // Get level set info
            sint tBSplineMeshIndex = aGeometryParameterList.get<sint>("bspline_mesh_index");
            real tBSplineLowerBound = aGeometryParameterList.get<real>("bspline_lower_bound");
            real tBSplineUpperBound = aGeometryParameterList.get<real>("bspline_upper_bound");

            // Build Geometry
            if (tGeometryType == "circle")
            {
                return std::make_shared<Circle>(
                        aADVs,
                        tGeometryVariableIndices,
                        tADVIndices,
                        tConstantParameters,
                        tGeometryName,
                        tNumRefinements,
                        tRefinementFunctionIndex,
                        tBSplineMeshIndex,
                        tBSplineLowerBound,
                        tBSplineUpperBound);
            }
            else if (tGeometryType == "superellipse")
            {
                return std::make_shared<Superellipse>(
                        aADVs,
                        tGeometryVariableIndices,
                        tADVIndices,
                        tConstantParameters,
                        tGeometryName,
                        tNumRefinements,
                        tRefinementFunctionIndex,
                        tBSplineMeshIndex,
                        tBSplineLowerBound,
                        tBSplineUpperBound);
            }
            else if (tGeometryType == "sphere")
            {
                return std::make_shared<Sphere>(
                        aADVs,
                        tGeometryVariableIndices,
                        tADVIndices,
                        tConstantParameters,
                        tGeometryName,
                        tNumRefinements,
                        tRefinementFunctionIndex,
                        tBSplineMeshIndex,
                        tBSplineLowerBound,
                        tBSplineUpperBound);
            }
            else if (tGeometryType == "superellipsoid")
            {
                return std::make_shared<Superellipsoid>(
                        aADVs,
                        tGeometryVariableIndices,
                        tADVIndices,
                        tConstantParameters,
                        tGeometryName,
                        tNumRefinements,
                        tRefinementFunctionIndex,
                        tBSplineMeshIndex,
                        tBSplineLowerBound,
                        tBSplineUpperBound);
            }
            else if (tGeometryType == "plane")
            {
                return std::make_shared<Plane>(
                        aADVs,
                        tGeometryVariableIndices,
                        tADVIndices,
                        tConstantParameters,
                        tGeometryName,
                        tNumRefinements,
                        tRefinementFunctionIndex,
                        tBSplineMeshIndex,
                        tBSplineLowerBound,
                        tBSplineUpperBound);
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
                return std::make_shared<User_Defined_Geometry>(
                        aADVs,
                        tGeometryVariableIndices,
                        tADVIndices,
                        tConstantParameters,
                        aLibrary->load_gen_field_function(aGeometryParameterList.get<std::string>("field_function_name")),
                        tSensitivityFunction,
                        tGeometryName,
                        tNumRefinements,
                        tRefinementFunctionIndex,
                        tBSplineMeshIndex,
                        tBSplineLowerBound,
                        tBSplineUpperBound);
            }
            else if (tGeometryType == "swiss_cheese_slice")
            {
                // Check for definition
                uint tNumXHoles = (uint)aGeometryParameterList.get<sint>("number_of_x_holes");
                uint tNumYHoles = (uint)aGeometryParameterList.get<sint>("number_of_y_holes");
                real tTargetXSpacing = aGeometryParameterList.get<real>("target_x_spacing");
                real tTargetYSpacing = aGeometryParameterList.get<real>("target_y_spacing");

                MORIS_ERROR((tNumXHoles > 1 and tNumYHoles > 1) or (tTargetXSpacing and tTargetYSpacing),
                            "In a swiss cheese parameter list, you must specify either a number of holes > 1 or a target "
                            "spacing in each direction.");

                if (tNumXHoles)
                {
                    return std::make_shared<Swiss_Cheese_Slice>(
                            aGeometryParameterList.get<real>("left_bound"),
                            aGeometryParameterList.get<real>("right_bound"),
                            aGeometryParameterList.get<real>("bottom_bound"),
                            aGeometryParameterList.get<real>("top_bound"),
                            tNumXHoles,
                            tNumYHoles,
                            aGeometryParameterList.get<real>("hole_x_semidiameter"),
                            aGeometryParameterList.get<real>("hole_y_semidiameter"),
                            aGeometryParameterList.get<real>("superellipse_exponent"),
                            aGeometryParameterList.get<real>("superellipse_scaling"),
                            aGeometryParameterList.get<real>("superellipse_regularization"),
                            aGeometryParameterList.get<real>("superellipse_shift"),
                            aGeometryParameterList.get<real>("row_offset"),
                            tGeometryName,
                            tNumRefinements,
                            tRefinementFunctionIndex,
                            tBSplineMeshIndex,
                            tBSplineLowerBound,
                            tBSplineUpperBound);
                }
                else
                {
                    return std::make_shared<Swiss_Cheese_Slice>(
                            aGeometryParameterList.get<real>("left_bound"),
                            aGeometryParameterList.get<real>("right_bound"),
                            aGeometryParameterList.get<real>("bottom_bound"),
                            aGeometryParameterList.get<real>("top_bound"),
                            tTargetXSpacing,
                            tTargetYSpacing,
                            aGeometryParameterList.get<real>("hole_x_semidiameter"),
                            aGeometryParameterList.get<real>("hole_y_semidiameter"),
                            aGeometryParameterList.get<real>("superellipse_exponent"),
                            aGeometryParameterList.get<real>("superellipse_scaling"),
                            aGeometryParameterList.get<real>("superellipse_regularization"),
                            aGeometryParameterList.get<real>("superellipse_shift"),
                            aGeometryParameterList.get<real>("row_offset"),
                            aGeometryParameterList.get<bool>("allow_less_than_target_spacing"),
                            tGeometryName,
                            tNumRefinements,
                            tRefinementFunctionIndex,
                            tBSplineMeshIndex,
                            tBSplineLowerBound,
                            tBSplineUpperBound);
                }

            }
            else
            {
                MORIS_ERROR(false, tGeometryType.append(" is not recognized as a valid Geometry type in fn_GEN_create_geometry.").c_str());
                return nullptr;
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        std::shared_ptr<Geometry> create_geometry(
                ParameterList               aGeometryParameterList,
                sol::Dist_Vector*           aOwnedADVs,
                std::shared_ptr<Library_IO> aLibrary)
        {
            // Geometry type
            std::string tGeometryType = aGeometryParameterList.get<std::string>("type");

            // Geometry inputs
            Matrix<DDUMat> tGeometryVariableIndices(0, 0);
            Matrix<DDUMat> tADVIndices(0, 0);
            Matrix<DDRMat> tConstantParameters(0, 0);

            // If not a swiss cheese, get ADV inputs
            if (tGeometryType.compare(0, 12, "swiss_cheese"))
            {
                // ADV parameters
                set_geometry_variable_inputs(aGeometryParameterList, aOwnedADVs->vec_local_length(), tGeometryVariableIndices, tADVIndices);

                // Constant parameters
                tConstantParameters = string_to_mat<DDRMat>(aGeometryParameterList.get<std::string>("constant_parameters"));
            }

            // Get name
            std::string tGeometryName = aGeometryParameterList.get<std::string>("name");

            // Get refinement info
            sint tNumRefinements = aGeometryParameterList.get<sint>("number_of_refinements");
            sint tRefinementFunctionIndex = aGeometryParameterList.get<sint>("refinement_function_index");

            // Get level set info
            sint tBSplineMeshIndex = aGeometryParameterList.get<sint>("bspline_mesh_index");
            real tBSplineLowerBound = aGeometryParameterList.get<real>("bspline_lower_bound");
            real tBSplineUpperBound = aGeometryParameterList.get<real>("bspline_upper_bound");

            // Build Geometry
            if (tGeometryType == "circle")
            {
                return std::make_shared<Circle>(
                        aOwnedADVs,
                        tGeometryVariableIndices,
                        tADVIndices,
                        tConstantParameters,
                        tGeometryName,
                        tNumRefinements,
                        tRefinementFunctionIndex,
                        tBSplineMeshIndex,
                        tBSplineLowerBound,
                        tBSplineUpperBound);
            }
            else if (tGeometryType == "superellipse")
            {
                return std::make_shared<Superellipse>(
                        aOwnedADVs,
                        tGeometryVariableIndices,
                        tADVIndices,
                        tConstantParameters,
                        tGeometryName,
                        tNumRefinements,
                        tRefinementFunctionIndex,
                        tBSplineMeshIndex,
                        tBSplineLowerBound,
                        tBSplineUpperBound);
            }
            else if (tGeometryType == "sphere")
            {
                return std::make_shared<Sphere>(
                        aOwnedADVs,
                        tGeometryVariableIndices,
                        tADVIndices,
                        tConstantParameters,
                        tGeometryName,
                        tNumRefinements,
                        tRefinementFunctionIndex,
                        tBSplineMeshIndex,
                        tBSplineLowerBound,
                        tBSplineUpperBound);
            }
            else if (tGeometryType == "superellipsoid")
            {
                return std::make_shared<Superellipsoid>(
                        aOwnedADVs,
                        tGeometryVariableIndices,
                        tADVIndices,
                        tConstantParameters,
                        tGeometryName,
                        tNumRefinements,
                        tRefinementFunctionIndex,
                        tBSplineMeshIndex,
                        tBSplineLowerBound,
                        tBSplineUpperBound);
            }
            else if (tGeometryType == "plane")
            {
                return std::make_shared<Plane>(
                        aOwnedADVs,
                        tGeometryVariableIndices,
                        tADVIndices,
                        tConstantParameters,
                        tGeometryName,
                        tNumRefinements,
                        tRefinementFunctionIndex,
                        tBSplineMeshIndex,
                        tBSplineLowerBound,
                        tBSplineUpperBound);
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
                return std::make_shared<User_Defined_Geometry>(
                        aOwnedADVs,
                        tGeometryVariableIndices,
                        tADVIndices,
                        tConstantParameters,
                        aLibrary->load_gen_field_function(aGeometryParameterList.get<std::string>("field_function_name")),
                        tSensitivityFunction,
                        tGeometryName,
                        tNumRefinements,
                        tRefinementFunctionIndex,
                        tBSplineMeshIndex,
                        tBSplineLowerBound,
                        tBSplineUpperBound);
            }
            else if (tGeometryType == "swiss_cheese_slice")
            {
                // Check for definition
                uint tNumXHoles = (uint)aGeometryParameterList.get<sint>("number_of_x_holes");
                uint tNumYHoles = (uint)aGeometryParameterList.get<sint>("number_of_y_holes");
                real tTargetXSpacing = aGeometryParameterList.get<real>("target_x_spacing");
                real tTargetYSpacing = aGeometryParameterList.get<real>("target_y_spacing");

                MORIS_ERROR((tNumXHoles > 1 and tNumYHoles > 1) or (tTargetXSpacing and tTargetYSpacing),
                            "In a swiss cheese parameter list, you must specify either a number of holes > 1 or a target "
                            "spacing in each direction.");

                if (tNumXHoles)
                {
                    return std::make_shared<Swiss_Cheese_Slice>(
                            aGeometryParameterList.get<real>("left_bound"),
                            aGeometryParameterList.get<real>("right_bound"),
                            aGeometryParameterList.get<real>("bottom_bound"),
                            aGeometryParameterList.get<real>("top_bound"),
                            tNumXHoles,
                            tNumYHoles,
                            aGeometryParameterList.get<real>("hole_x_semidiameter"),
                            aGeometryParameterList.get<real>("hole_y_semidiameter"),
                            aGeometryParameterList.get<real>("superellipse_exponent"),
                            aGeometryParameterList.get<real>("superellipse_scaling"),
                            aGeometryParameterList.get<real>("superellipse_regularization"),
                            aGeometryParameterList.get<real>("superellipse_shift"),
                            aGeometryParameterList.get<real>("row_offset"),
                            tGeometryName,
                            tNumRefinements,
                            tRefinementFunctionIndex,
                            tBSplineMeshIndex,
                            tBSplineLowerBound,
                            tBSplineUpperBound);
                }
                else
                {
                    return std::make_shared<Swiss_Cheese_Slice>(
                            aGeometryParameterList.get<real>("left_bound"),
                            aGeometryParameterList.get<real>("right_bound"),
                            aGeometryParameterList.get<real>("bottom_bound"),
                            aGeometryParameterList.get<real>("top_bound"),
                            tTargetXSpacing,
                            tTargetYSpacing,
                            aGeometryParameterList.get<real>("hole_x_semidiameter"),
                            aGeometryParameterList.get<real>("hole_y_semidiameter"),
                            aGeometryParameterList.get<real>("superellipse_exponent"),
                            aGeometryParameterList.get<real>("superellipse_scaling"),
                            aGeometryParameterList.get<real>("superellipse_regularization"),
                            aGeometryParameterList.get<real>("superellipse_shift"),
                            aGeometryParameterList.get<real>("row_offset"),
                            aGeometryParameterList.get<bool>("allow_less_than_target_spacing"),
                            tGeometryName,
                            tNumRefinements,
                            tRefinementFunctionIndex,
                            tBSplineMeshIndex,
                            tBSplineLowerBound,
                            tBSplineUpperBound);
                }

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
