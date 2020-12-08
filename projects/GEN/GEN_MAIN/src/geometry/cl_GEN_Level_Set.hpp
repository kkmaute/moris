#ifndef MORIS_CL_GEN_LEVEL_SET_HPP
#define MORIS_CL_GEN_LEVEL_SET_HPP

#include "cl_GEN_Geometry.hpp"
#include "cl_GEN_Field_Discrete_Integration.hpp"
#include "cl_MTK_Mesh_Core.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"

namespace moris
{
    namespace ge
    {

        class Level_Set : public Geometry, public Field_Discrete_Integration
        {

        private:
            mtk::Interpolation_Mesh* mMesh;
            sol::Dist_Vector* mOwnedNodalValues = nullptr;
            sol::Dist_Vector* mSharedNodalValues = nullptr;

        public:

            /**
             * Constructor where ADVs are added based on an input field and a B-spline mesh.
             *
             * @param aOwnedADVs Pointer to the owned distributed ADVs
             * @param aOwnedADVIds All owned ADV IDs on this processor
             * @param aSharedADVIds All owned and shared ADV IDs for this B-spline field
             * @param aOwnedADVIdsOffset Offset in the owned ADV IDs for pulling ADV IDs
             * @param aMesh The mesh pointer where the B-spline information can be obtained
             * @param aGeometry Geometry for initializing the B-spline level set discretization
             */
            Level_Set(sol::Dist_Vector*         aOwnedADVs,
                      const Matrix<DDSMat>&     aOwnedADVIds,
                      const Matrix<DDSMat>&     aSharedADVIds,
                      uint                      aOwnedADVIdsOffset,
                      mtk::Interpolation_Mesh*  aMesh,
                      std::shared_ptr<Geometry> aGeometry);

            /**
             * Destructor
             */
            ~Level_Set();

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

            /**
             * Function for determining if this geometry is to be used for seeding a B-spline level set field.
             *
             * @return false
             */
            bool conversion_to_bsplines();

        private:

            /**
             * Maps the level set field from nodes to B-splines for the given geometry.
             *
             * @return Target field
             */
            Matrix<DDRMat> map_to_bsplines(std::shared_ptr<Geometry> aGeometry);

        };
    }
}

#endif //MORIS_CL_GEN_LEVEL_SET_HPP
