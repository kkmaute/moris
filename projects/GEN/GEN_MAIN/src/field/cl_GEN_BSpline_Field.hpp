#ifndef MORIS_CL_GEN_BSPLINE_FIELD_HPP
#define MORIS_CL_GEN_BSPLINE_FIELD_HPP

#include "cl_GEN_Field_Discrete_Integration.hpp"

namespace moris
{
    namespace ge
    {
        class BSpline_Field : public Field_Discrete_Integration
        {

        private:
            Matrix<DDSMat> mSharedADVIds;
            mtk::Interpolation_Mesh* mMesh;
            sol::Dist_Vector* mOwnedNodalValues = nullptr;
            sol::Dist_Vector* mSharedNodalValues = nullptr;

        public:
            /**
             * Constructor where ADVs are added based on an input field and a B-spline mesh.
             *
             * @param aOwnedADVs Pointer to the owned distributed ADVs
             * @param aCoefficientIndices Coefficient indices to be mapped to
             * @param aSharedADVIds All owned and shared ADV IDs for this B-spline field
             * @param aADVOffsetID Offset in the owned ADV IDs for pulling ADV IDs
             * @param aMesh The mesh pointer where the B-spline information can be obtained
             * @param aField Field for initializing the B-spline level set discretization
             */
            BSpline_Field(
                    sol::Dist_Vector*        aOwnedADVs,
                    const Matrix<DDUMat>&    aCoefficientIndices,
                    const Matrix<DDSMat>&    aSharedADVIds,
                    uint                     aADVOffsetID,
                    mtk::Interpolation_Mesh* aMesh,
                    std::shared_ptr<Field>   aField);

            /**
             * Destructor
             */
            ~BSpline_Field();

            /**
             * Given a node index, returns the field value.
             *
             * @param aNodeIndex Node index
             * @return Distance to this geometry
             */
            real get_field_value(uint aNodeIndex);

            /**
             * Given a node index, evaluates the sensitivity of the geometry field with respect to all of the
             * geometry variables.
             *
             * @param aNodeIndex Node index
             * @return Vector of sensitivities
             */
            const Matrix<DDRMat>& get_field_sensitivities(uint aNodeIndex);

            /**
             * Gets the IDs of ADVs which this field depends on for evaluations.
             *
             * @param aNodeIndex Node index
             * @param aCoordinates Node coordinates
             * @return Determining ADV IDs at this node
             */
            virtual Matrix<DDSMat> get_determining_adv_ids(uint aNodeIndex);

            /**
             * Imports the local ADVs required from the full owned ADV distributed vector, and recomputes nodal values.
             *
             * @param aOwnedADVs Full owned distributed ADV vector
             */
            void import_advs(sol::Dist_Vector* aOwnedADVs);

        protected:

            /**
             * Gets the mesh that this field depends on.
             *
             * @return Mesh, default nullptr
             */
            mtk::Interpolation_Mesh* get_mesh();

        private:

            /**
             * Maps the level set field from nodes to B-splines for the given geometry.
             *
             * @return Target field
             */
            Matrix<DDRMat> map_to_bsplines(std::shared_ptr<Field> aField);

        };
    }
}

#endif //MORIS_CL_GEN_BSPLINE_FIELD_HPP
