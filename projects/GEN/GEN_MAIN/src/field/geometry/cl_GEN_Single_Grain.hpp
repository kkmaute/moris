#ifndef MORIS_CL_GEN_GRAIN_HPP
#define MORIS_CL_GEN_GRAIN_HPP

#include "cl_GEN_Geometry.hpp"
#include "cl_GEN_Field_Analytic.hpp"
#include "fn_Exec_load_user_library.hpp"

namespace moris
{
    namespace ge
    {
        class Single_Grain : public Geometry, public Field_Analytic
        {

        private:

                std::shared_ptr<Geometry>      mVoxelGeometrie = nullptr;
                uint                           mIndex = MORIS_UINT_MAX;


        public:

            /**
             * Constructor, sets the pointers to advs and constant parameters for evaluations.
             *
             * @param aADVs Reference to the full advs
             * @param aGeometryVariableIndices Indices of geometry variables to be filled by the ADVs
             * @param aADVIndices The indices of the ADV vector to fill in the geometry variables
             * @param aConstantParameters The constant parameters not filled by ADVs
             * @param aName Name of this field for identification
             * @param aNumRefinements The number of refinement steps to use for this field
             * @param aRefinementMeshIndices Indices of meshes to perform refinement on
             * @param aRefinementFunctionIndex The index of a user-defined refinement function (-1 = default refinement)
             */
                Single_Grain(
                    Matrix<DDRMat>&                aADVs,
                    Matrix<DDUMat>                 aGeometryVariableIndices,
                    Matrix<DDUMat>                 aADVIndices,
                    Matrix<DDRMat>                 aConstantParameters,
                    std::shared_ptr<Geometry>      aVoxelGeometrie,
                    uint                           aIndex,
                    std::string                    aName = "",
                    Matrix<DDSMat>                 aNumRefinements = {{}},
                    Matrix<DDSMat>                 aRefinementMeshIndices = {{}},
                    sint                           aRefinementFunctionIndex = -1,
                    sint                           aBSplineMeshIndex = -2);

            /**
             * Constructor, sets the field variable pointers to ADVs and constant parameters for evaluations.
             *
             * @param aOwnedADVs Pointer to the owned distributed ADVs
             * @param aFieldVariableIndices Indices of geometry variables to be filled by the ADVs
             * @param aADVIndices The indices of the ADV vector to fill in the geometry variables
             * @param aConstantParameters The constant parameters not filled by ADVs
             * @param aName Name of this field for identification
             * @param aNumRefinements The number of refinement steps to use for this field
             * @param aRefinementMeshIndices Indices of meshes to perform refinement on
             * @param aRefinementFunctionIndex The index of a user-defined refinement function (-1 = default refinement)
             * @param aBSplineMeshIndex Index of a B-spline mesh for discretization (-2 = none, -1 = store nodal values)
             */
                Single_Grain(
                    sol::Dist_Vector* aOwnedADVs,
                    Matrix<DDUMat>                 aGeometryVariableIndices,
                    Matrix<DDUMat>                 aADVIndices,
                    Matrix<DDRMat>                 aConstantParameters,
                    std::shared_ptr<Geometry>      aVoxelGeometrie,
                    uint                           aIndex,
                    std::string                    aName = "",
                    Matrix<DDSMat>                 aNumRefinements = {{}},
                    Matrix<DDSMat>                 aRefinementMeshIndices = {{}},
                    sint                           aRefinementFunctionIndex = -1,
                    sint                           aBSplineMeshIndex = -2);

            /**
             * Constructor with only constant parameters
             *
             * @param aConstantParameters The constant parameters not filled by ADVs
             * @param aName Name of this field for identification
             * @param aNumRefinements The number of refinement steps to use for this field
             * @param aRefinementMeshIndices Indices of meshes to perform refinement on
             * @param aRefinementFunctionIndex The index of a user-defined refinement function (-1 = default refinement)
             * @param aBSplineMeshIndex Index of a B-spline mesh for discretization (-2 = none, -1 = store nodal values)
             */
                Single_Grain(
                    Matrix<DDRMat>           aConstantParameters,
                    std::shared_ptr<Geometry>      aVoxelGeometrie,
                    uint                           aIndex,
                    std::string              aName = "",
                    Matrix<DDSMat>           aNumRefinements = {{}},
                    Matrix<DDSMat>           aRefinementMeshIndices = {{}},
                    sint                     aRefinementFunctionIndex = -1,
                    sint                     aBSplineMeshIndex = -2);

            /**
             * Given a node coordinate, returns the field value.
             *
             * @param aCoordinates Coordinate values
             * @return Distance to this geometry
             */
            real get_field_value(const Matrix<DDRMat>& aCoordinates);

            /**
             * Given a node coordinate, evaluates the sensitivity of the geometry field with respect to all of the
             * geometry variables.
             *
             * @param aCoordinates Coordinate values
             * @return Vector of sensitivities
             */
            const Matrix<DDRMat>& get_field_sensitivities(const Matrix<DDRMat>& aCoordinates);

        private:


        };
    }
}

#endif //MORIS_CL_GEN_VOXEL_HPP
